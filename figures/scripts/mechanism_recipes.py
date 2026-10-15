#!/usr/bin/env python3
"""从显式源表绘制机制证据；不拟合模型、不补齐缺失测量。"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".cache/matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
from figure_style import PROFILE, COLORS, DIVERGING, apply_style, save_figure


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def require(frame, columns, unique=()):
    missing = set(columns) - set(frame.columns)
    if missing:
        raise ValueError(f"缺少源表字段: {sorted(missing)}")
    if unique and frame.duplicated(list(unique)).any():
        raise ValueError(f"源表标识重复: {unique}")
    if frame.empty:
        raise ValueError("没有可绘制的源数据；不能生成空白成功图")


def cmap():
    result = LinearSegmentedColormap.from_list("evidence_direction", DIVERGING)
    result.set_bad(PROFILE["scales"]["missing_color"])
    return result


def numeric(frame, columns):
    for name in columns:
        before = frame[name]
        value = pd.to_numeric(before, errors="coerce")
        if (before.notna() & value.isna()).any() or np.isinf(value).any():
            raise ValueError(f"{name} 含非数值或无限值")
        frame[name] = value


def text_order(frame, field):
    if frame[field].isna().any():
        raise ValueError(f"{field} 标识缺失")
    return list(dict.fromkeys(frame[field].astype(str)))


def drivers(data, secondary):
    require(data, ["pathway", "context", "gene", "log2fc", "leading_edge", "NES"], ["pathway", "context", "gene"])
    require(secondary, ["context_a", "context_b", "jaccard", "shared_n", "universe_n"])
    numeric(data, ["log2fc", "leading_edge", "NES"])
    numeric(secondary, ["jaccard", "shared_n", "universe_n"])
    if data.pathway.nunique() != 1:
        raise ValueError("一图一个程序；调用者须显式声明显示程序")
    if not data.leading_edge.dropna().isin([0, 1]).all():
        raise ValueError("leading_edge 只能是0/1/缺测")
    if ((secondary.jaccard < 0) | (secondary.jaccard > 1)).any():
        raise ValueError("Jaccard超出0–1")
    genes, contexts = text_order(data, "gene"), text_order(data, "context")
    if len(genes) > 24:
        raise ValueError("超过24基因，须由调用者记录明确显示子集")
    expression = data.pivot(index="gene", columns="context", values="log2fc").reindex(index=genes, columns=contexts)
    leading = data.pivot(index="gene", columns="context", values="leading_edge").reindex(index=genes, columns=contexts)
    nes_check = data.groupby("context").NES.nunique()
    if (nes_check > 1).any():
        raise ValueError("同一情境出现不同NES")
    nes = data.groupby("context", sort=False).NES.first().reindex(contexts)
    fig = plt.figure(figsize=(183 / 25.4, min(230, 80 + 5 * len(genes)) / 25.4))
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 1.05], height_ratios=[1, max(5, len(genes))], wspace=.7, hspace=.18)
    top = fig.add_subplot(grid[0, 0])
    left = fig.add_subplot(grid[1, 0])
    right = fig.add_subplot(grid[:, 1])
    fig.subplots_adjust(left=.22, right=.94, top=.88, bottom=.15)
    top.imshow(nes.to_numpy()[None, :], cmap=cmap(), norm=TwoSlopeNorm(0, -3, 3), aspect="auto")
    top.set(yticks=[0], yticklabels=["NES"], xticks=[])
    for j, value in enumerate(nes):
        if np.isfinite(value):
            top.text(j, 0, f"{value:+.2f}", ha="center", va="center", fontsize=7)
    im = left.imshow(expression, cmap=cmap(), norm=TwoSlopeNorm(0, -2, 2), aspect="auto")
    y, x = np.where(leading.to_numpy() == 1)
    left.scatter(x, y, s=11, facecolors="none", edgecolors="#374151", linewidths=.6)
    left.set(yticks=range(len(genes)), yticklabels=genes, xticks=range(len(contexts)), xticklabels=contexts)
    left.tick_params(axis="x", rotation=40)
    cb = fig.colorbar(im, ax=left, orientation="horizontal", fraction=.055, pad=.13, ticks=[-2, 0, 2])
    cb.set_label("RNA log2FC; circles: leading edge")
    cb.ax.set_xticklabels(["≤−2", "0", "≥2"])
    pairs = secondary.context_a.astype(str) + " / " + secondary.context_b.astype(str)
    right.barh(range(len(pairs)), secondary.jaccard, color=PROFILE["palettes"]["accent"]["teal"], height=.58)
    for i, row in enumerate(secondary.itertuples()):
        right.text(float(row.jaccard) + .025, i, f"{int(row.shared_n)} shared", va="center", fontsize=7)
    right.set(yticks=range(len(pairs)), yticklabels=pairs, xlim=(0, 1), xlabel="Leading-edge Jaccard", title="Shared drivers across contexts")
    right.invert_yaxis()
    for label,ax in zip('ABC',[top,left,right]):
        ax.text(-.08,1.03,label,transform=ax.transAxes,fontweight='bold',va='bottom',ha='right')
    title = str(data.pathway.iloc[0]).removeprefix("HALLMARK_").replace("_", " ")
    if 'comparison' in data: title += '\n' + str(data.comparison.iloc[0])
    fig.suptitle(title, y=.97, fontsize=8)
    fig.text(.02, .02, "Independent contexts; overlap background = jointly tested members of this pathway.", fontsize=7)
    return fig


def regulatory(data, secondary=None):
    columns = ["tf_activity_delta", "target_rna_log2fc", "motif_deviation_delta", "methylation_delta"]
    require(data, ["edge_id", "tf", "target", "motif_label", "methylation_qualified"] + columns, ["edge_id"])
    numeric(data, columns)
    if len(data) > 24:
        raise ValueError("超过24条边，调用者必须声明显示子集")
    qualified = data.methylation_qualified.astype(str).str.lower().isin(["true", "1"])
    if data.loc[~qualified, "methylation_delta"].notna().any():
        raise ValueError("未通过直接覆盖的甲基化差值不能显示为观测值")
    fig, axes = plt.subplots(1, 4, figsize=(183 / 25.4, min(230, 60 + 5 * len(data)) / 25.4), sharey=True)
    fig.subplots_adjust(left=.31, right=.97, top=.82, bottom=.2, wspace=.45)
    titles = ["RNA-derived\nTF score Δ", "Target RNA\nlog2FC", "Motif accessibility\ndeviation Δ", "Direct CpG\nfraction Δ"]
    limits = [max(1, np.nanmax(np.abs(data[columns[0]]))) if data[columns[0]].notna().any() else 1, 2,
              max(1, np.nanmax(np.abs(data[columns[2]]))) if data[columns[2]].notna().any() else 1, 1]
    labels = [f"{r.tf} → {r.target}" for r in data.itertuples()]
    for i, (ax, col, title, limit) in enumerate(zip(axes, columns, titles, limits)):
        im = ax.imshow(data[col].to_numpy()[:, None], cmap=cmap(), norm=TwoSlopeNorm(0, -limit, limit), aspect="auto")
        ax.set(xticks=[], yticks=range(len(data)), title=title)
        if i == 0:
            ax.set_yticklabels(labels)
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=.07, pad=.04, ticks=[-limit, 0, limit])
        cb.ax.tick_params(labelsize=7)
        ax.text(-.10,1.13,'ABCD'[i],transform=ax.transAxes,fontweight='bold',va='bottom',ha='right')
    title="Regulatory evidence availability: RNA / ATAC / direct methylation"
    if 'comparison' in data: title += '\n' + str(data.comparison.iloc[0])
    fig.suptitle(title, y=.98, fontsize=8)
    fig.text(.02, .045, "Arrows denote supplied prior edges, not validated binding. Grey = missing or ineligible.\nTF scores and target RNA are one modality; motif deviation is not TF occupancy.", fontsize=7)
    return fig


def cpg(data, secondary=None):
    require(data, ["locus_id", "context", "sample_id", "condition", "chr", "position", "methylated_reads", "total_reads"],
            ["locus_id", "context", "sample_id", "chr", "position"])
    numeric(data, ["position", "methylated_reads", "total_reads"])
    if data.locus_id.nunique() != 1 or data.chr.nunique() != 1:
        raise ValueError("一图一个同染色体位点")
    m, n = data.methylated_reads, data.total_reads
    if ((m < 0) | (n < 0) | (m > n) | (m.notna() & (m % 1 != 0)) | (n.notna() & (n % 1 != 0))).any():
        raise ValueError("甲基化read数/总覆盖无效")
    if (m.notna() != n.notna()).any():
        raise ValueError("read数和覆盖必须同时存在或同时缺失")
    data = data.copy()
    data["fraction"] = (m / n).where(n > 0)
    data["row"] = data.context.astype(str) + " | " + data.sample_id.astype(str)
    rows = text_order(data, "row")
    positions = sorted(data.position.unique())
    if len(positions) > 500:
        raise ValueError("超过500 CpG；须先指定更小的可追溯位点窗口")
    pct = data.pivot(index="row", columns="position", values="fraction").reindex(index=rows, columns=positions)
    cov = data.pivot(index="row", columns="position", values="total_reads").reindex(index=rows, columns=positions)
    fig, axes = plt.subplots(2, 1, figsize=(183 / 25.4, min(235, 75 + len(rows) * 7) / 25.4), sharex=True)
    fig.subplots_adjust(left=.26, right=.89, top=.87, bottom=.24, hspace=.25)
    cm = LinearSegmentedColormap.from_list("methylation", ["#F7F7F7", PROFILE["palettes"]["accent"]["purple"]])
    cm.set_bad(PROFILE["scales"]["missing_color"])
    im = axes[0].imshow(pct, aspect="auto", cmap=cm, vmin=0, vmax=1)
    fig.colorbar(im, ax=axes[0], fraction=.025, pad=.025, label="Methylated fraction", ticks=[0, .5, 1])
    depth_cmap = plt.get_cmap("Blues").copy()
    depth_cmap.set_bad(PROFILE["scales"]["missing_color"])
    im2 = axes[1].imshow(np.log2(cov + 1), aspect="auto", cmap=depth_cmap, vmin=0)
    fig.colorbar(im2, ax=axes[1], fraction=.025, pad=.025, label="log2(depth + 1)")
    for ax in axes:
        ax.set(yticks=range(len(rows)), yticklabels=rows)
    yy, xx = np.where((cov.to_numpy() > 0) & (cov.to_numpy() < 10))
    axes[0].scatter(xx, yy, marker="x", s=10, color="#374151", linewidths=.5)
    every = max(1, int(np.ceil(len(positions) / 12)))
    ticks = list(range(0, len(positions), every))
    chrom=str(data.chr.iloc[0]);chrom=chrom if chrom.startswith('chr') else 'chr'+chrom
    assembly=str(data.genome_build.iloc[0])+', ' if 'genome_build' in data else ''
    axes[1].set(xticks=ticks, xticklabels=[str(int(positions[i])) for i in ticks], xlabel=f"{chrom} CpG coordinate ({assembly}1-based)")
    axes[1].tick_params(axis="x", rotation=50)
    axes[0].set_title("Direct per-sample methylation")
    axes[1].set_title("Sequencing coverage")
    for label,ax in zip('AB',axes):
        ax.text(-.08,1.03,label,transform=ax.transAxes,fontweight='bold',va='bottom',ha='right')
    title=str(data.locus_id.iloc[0])
    if 'candidate_symbol' in data:title=str(data.candidate_symbol.iloc[0])+' | '+title
    fig.suptitle(title, y=.97, fontsize=8)
    depth_label = "<10"
    if "source_mincov" in data and data.source_mincov.nunique() == 1:
        minimum = int(data.source_mincov.iloc[0])
        if 0 < minimum < 10:
            depth_label = f"{minimum}–9"
    fig.text(.02, .025, f"Grey = unavailable (including censored source records); × = observed depth {depth_label}.\nCpG columns are categorical; stages are different individuals. Bisulfite does not separate 5mC/5hmC.", fontsize=7)
    return fig


def reactions(data, secondary):
    require(data, ["node_id", "label", "x", "y", "measurement_state", "effect"], ["node_id"])
    require(secondary, ["source", "target", "label", "edge_type", "source_url"])
    numeric(data, ["x", "y", "effect"])
    if data[['x','y']].isna().any().any():
        raise ValueError("节点坐标必须完整")
    ids = set(data.node_id)
    if not set(secondary.source).union(secondary.target) <= ids:
        raise ValueError("反应边端点缺失")
    if secondary.source_url.isna().any() or secondary.source_url.astype(str).str.strip().eq("").any():
        raise ValueError("反应边必须有独立可核验的来源")
    allowed = {"measured_level1", "candidate_level2", "unmeasured", "biochemical_assay"}
    if not set(data.measurement_state) <= allowed:
        raise ValueError("未知化合物测量状态")
    if data.loc[data.measurement_state.eq("unmeasured"), "effect"].notna().any():
        raise ValueError("未测节点不能带效应")
    fig, ax = plt.subplots(figsize=(183 / 25.4, 170 / 25.4))
    fig.subplots_adjust(left=.05, right=.95, top=.88, bottom=.18)
    lookup = data.set_index("node_id")
    ax.set(xlim=(data.x.min() - .7, data.x.max() + .7), ylim=(data.y.min() - .6, data.y.max() + .6))
    node_artists={}
    for row in data.itertuples():
        observed = row.measurement_state in {"measured_level1", "biochemical_assay"}
        color = cmap()(Normalize(-2, 2)(row.effect)) if np.isfinite(row.effect) else PROFILE["palettes"]["neutral"]["light"]
        node_artists[row.node_id]=ax.text(row.x, row.y, row.label, ha="center", va="center", fontsize=8,zorder=3,
                bbox={"boxstyle": "round,pad=.4", "facecolor": color,
                      "edgecolor": "#374151" if observed else "#B7BEC6",
                      "linestyle": "--" if row.measurement_state == "candidate_level2" else "-", "linewidth": .8})
    # Render node boxes before clipped arrows so each PDF/PNG renderer updates their bounds.
    fig.canvas.draw()
    for index, edge in enumerate(secondary.itertuples()):
        a, b = lookup.loc[edge.source], lookup.loc[edge.target]
        curve = getattr(edge, "curvature", 0)
        curve = 0 if pd.isna(curve) else float(curve)
        reversible = str(getattr(edge, "reversible", False)).lower() in {"true", "1"}
        arrow = FancyArrowPatch((a.x, a.y), (b.x, b.y), arrowstyle="<->" if reversible else "-|>", mutation_scale=9,
                                linewidth=PROFILE["typography"]["line_mm"] * 72 / 25.4,
                                color="#6B7280", shrinkA=3, shrinkB=3,
                                patchA=node_artists[edge.source].get_bbox_patch(),
                                patchB=node_artists[edge.target].get_bbox_patch(),
                                connectionstyle=f"arc3,rad={curve}",zorder=4)
        ax.add_patch(arrow)
        dx = getattr(edge, "label_dx", 0)
        dy = getattr(edge, "label_dy", .13)
        dx = 0 if pd.isna(dx) else float(dx)
        dy = .13 if pd.isna(dy) else float(dy)
        ax.text((a.x + b.x) / 2 + dx, (a.y + b.y) / 2 + dy, str(edge.label), fontsize=7,
                ha="center", va="center", zorder=5, bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    ax.axis("off")
    title="Carbon allocation: sourced biochemical connections"
    if 'comparison' in data:title+='\nRNA effects: '+str(data.comparison.iloc[0])
    fig.suptitle(title, y=.97, fontsize=8)
    legend=str(data.effect_legend.iloc[0]) if 'effect_legend' in data else 'Edge labels follow the frozen source; transcripts are not enzyme activity.'
    fig.text(.04, .04, "Arrows = sourced reactions/transport, not flux; compartment pools were not measured.\nWhite nodes = no validated concentration; dashed border = putative level 2 identity.\n"+legend, fontsize=7)
    return fig


RENDERERS = {"drivers": drivers, "regulatory": regulatory, "cpg": cpg, "reactions": reactions}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", choices=list(RENDERERS), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--secondary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--data-kind", required=True, choices=["real_PFOS_source", "real_project_source", "public_source", "synthetic_example"])
    args = parser.parse_args()
    if args.kind in {"drivers", "reactions"} and args.secondary is None:
        parser.error("此图型需要--secondary来源表")
    apply_style()
    data = pd.read_csv(args.input, sep="\t")
    other = pd.read_csv(args.secondary, sep="\t") if args.secondary else None
    figure = RENDERERS[args.kind](data, other)
    args.output.mkdir(parents=True, exist_ok=True)
    outputs = save_figure(figure, args.output / "figure.pdf", height_mm=figure.get_figheight() * 25.4)
    plt.close(figure)
    inputs = {}
    for source, name in [(args.input, "source.tsv"), (args.secondary, "secondary.tsv")]:
        if source is not None:
            destination = args.output / name
            if source.resolve() != destination.resolve():
                shutil.copy2(source, destination)
            inputs[name] = sha(destination)
    profile_path = Path(__file__).resolve().parent / getattr(__import__("figure_style"), "PROFILE_RELATIVE")
    receipt = {"status": "COMPLETED", "kind": args.kind, "data_kind": args.data_kind,
               "intended_use": "internal_exploratory_preview", "publication_figure_gate": "not_run",
               "scientific_coordinate_count": {'drivers':3,'regulatory':4,'cpg':2,'reactions':1}[args.kind],
               "composition_method": "source_data_shared_canvas", "input_media_files": [],
               "atomic_panel_geometry_comparison": "not_run",
               "statistics_recomputed": False, "input_sha256": inputs,
               "script_sha256": sha(__file__), "profile_sha256": sha(profile_path),
               "style_sha256": sha(Path(__file__).with_name("figure_style.py")),
               "outputs": {key: sha(value) for key, value in outputs.items()}}
    (args.output / "figure.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "COMPLETED", "kind": args.kind, "rows": len(data)}))


if __name__ == "__main__":
    main()

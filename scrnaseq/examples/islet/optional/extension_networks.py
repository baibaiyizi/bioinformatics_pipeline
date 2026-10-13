from __future__ import annotations

import hashlib
import json
import math
import os
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial.distance import jensenshannon

from .plotting import DIVERGING_CMAP, CELL_TYPE_COLORS, save_figure
from .figure_style import COLORS
from .runtime import atomic_json, atomic_text, atomic_tsv, sha256


MIN_CELLS = 20
EXPR_PROP = 0.10
LIANA_RESOURCE = "mouseconsensus"
LIANA_PERMUTATIONS = 1000
LIANA_SEED = 20260905
PRIORITY_MAGNITUDE_RANK = 0.05
SCENIC_SEEDS = (1907, 2909, 3911)
SCENIC_MAX_CELLS_PER_STRATUM = 500
SCENIC_EDGE_TOP_PER_TF = 100
SCENIC_EDGE_MIN_SEEDS = 2
SCENIC_WORKERS = 16


def _utc_now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def _file_record(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {"path": str(resolved), "sha256": sha256(resolved), "size_bytes": resolved.stat().st_size}


def _write_execution_inputs(path: Path, files: Iterable[Path]) -> dict[str, Any]:
    payload = {"files": [_file_record(item) for item in files], "started_at": _utc_now()}
    atomic_json(path, payload)
    return payload


def _inputs_unchanged(payload: dict[str, Any]) -> bool:
    for item in payload["files"]:
        path = Path(item["path"])
        if not path.is_file() or path.stat().st_size != item["size_bytes"] or sha256(path) != item["sha256"]:
            return False
    return True


def _relative(path: Path, root: Path) -> str:
    return str(path.resolve().relative_to(root.resolve()))


def _subunits(entity: str) -> tuple[str, ...]:
    return tuple(part for part in str(entity).split("_") if part)


def _validate_adata(data: ad.AnnData) -> None:
    required_obs = {"library_id", "condition", "cell_type"}
    missing = sorted(required_obs - set(data.obs.columns))
    if missing:
        raise ValueError(f"网络分析缺少 obs 字段: {missing}")
    if "counts" not in data.layers:
        raise ValueError("网络分析要求整数 counts layer")
    counts = data.layers["counts"]
    if not sp.issparse(counts):
        raise ValueError("counts layer 必须为稀疏矩阵")
    if counts.data.size and (np.any(counts.data < 0) or not np.allclose(counts.data, np.rint(counts.data))):
        raise ValueError("counts layer 必须为非负整数 UMI")
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("barcode 与 gene index 必须唯一")


def _write_liana_resource_snapshot(path: Path) -> pd.DataFrame:
    import liana as li

    resource = li.resource.select_resource(LIANA_RESOURCE).copy()
    resource = resource[["ligand", "receptor"]].drop_duplicates().sort_values(["ligand", "receptor"])
    resource["ligand_subunits"] = resource["ligand"].map(lambda x: ";".join(_subunits(x)))
    resource["receptor_subunits"] = resource["receptor"].map(lambda x: ";".join(_subunits(x)))
    resource["species"] = "Mus musculus"
    resource["liana_resource"] = LIANA_RESOURCE
    atomic_tsv(path, resource)
    return resource


def _write_cellchat_pathway_snapshot(path: Path, rscript: Path) -> pd.DataFrame:
    """Export only the mouse CellChat pathway labels; LIANA remains the scoring method."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    code = r'''
suppressPackageStartupMessages(library(CellChat))
args <- commandArgs(trailingOnly=TRUE)
db <- CellChatDB.mouse
expand_entity <- function(x) {
  if (x %in% rownames(db$complex)) {
    vals <- unlist(db$complex[x, , drop=FALSE], use.names=FALSE)
    vals <- vals[!is.na(vals) & nzchar(vals)]
    return(paste(vals, collapse="_"))
  }
  x
}
out <- db$interaction[, c("interaction_name", "pathway_name", "ligand", "receptor", "annotation", "evidence")]
out$ligand_complex <- vapply(out$ligand, expand_entity, character(1))
out$receptor_complex <- vapply(out$receptor, expand_entity, character(1))
out <- unique(out[, c("ligand_complex", "receptor_complex", "pathway_name", "annotation", "evidence", "interaction_name")])
out <- out[order(out$ligand_complex, out$receptor_complex, out$pathway_name), ]
write.table(out, file=args[[1]], sep="\t", row.names=FALSE, quote=FALSE, na="")
'''
    try:
        subprocess.run([str(rscript), "-e", code, str(temporary)], check=True, text=True)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _detection_by_group(data: ad.AnnData) -> dict[tuple[str, str], float]:
    counts = data.layers["counts"].tocsr()
    result: dict[tuple[str, str], float] = {}
    labels = data.obs["cell_type"].astype(str).to_numpy()
    for label in sorted(set(labels)):
        idx = np.flatnonzero(labels == label)
        props = np.asarray((counts[idx] > 0).mean(axis=0)).ravel()
        result.update(((label, str(gene)), float(value)) for gene, value in zip(data.var_names, props, strict=True))
    return result


def _complex_gate(row: pd.Series, detection: dict[tuple[str, str], float], expr_prop: float = EXPR_PROP) -> tuple[bool, str, str]:
    ligand_props = [detection.get((str(row["source"]), gene), 0.0) for gene in _subunits(row["ligand_complex"])]
    receptor_props = [detection.get((str(row["target"]), gene), 0.0) for gene in _subunits(row["receptor_complex"])]
    ligand_min = min(ligand_props, default=0.0)
    receptor_min = min(receptor_props, default=0.0)
    passed = ligand_min >= expr_prop and receptor_min >= expr_prop
    return passed, ";".join(f"{x:.6g}" for x in ligand_props), ";".join(f"{x:.6g}" for x in receptor_props)


def _module_reusable(output_dir: Path, manifest_name: str, execution_name: str) -> dict[str, Any] | None:
    manifest_path, execution_path = output_dir / manifest_name, output_dir / execution_name
    if not manifest_path.is_file() or not execution_path.is_file():
        return None
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    execution = json.loads(execution_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "completed" or not _inputs_unchanged(execution):
        return None
    expected = [output_dir / path for path in manifest.get("tables", [])]
    for figure in manifest.get("figures", []):
        expected.extend(output_dir / figure[key] for key in ("png", "pdf", "source"))
    if not expected or any(not path.is_file() or path.stat().st_size == 0 for path in expected):
        return None
    return manifest


def _liana_library(
    lib: ad.AnnData,
    resource: pd.DataFrame,
    *,
    expr_prop: float,
    n_perms: int,
    seed: int,
) -> pd.DataFrame:
    import liana as li

    counts_by_type = lib.obs["cell_type"].astype(str).value_counts()
    eligible = counts_by_type[counts_by_type >= MIN_CELLS].index
    lib = lib[lib.obs["cell_type"].astype(str).isin(eligible)].copy()
    detection = _detection_by_group(lib)
    result = li.mt.rank_aggregate(
        lib,
        groupby="cell_type",
        resource=resource[["ligand", "receptor"]],
        expr_prop=expr_prop,
        min_cells=MIN_CELLS,
        use_raw=False,
        return_all_lrs=False,
        n_perms=n_perms,
        seed=seed,
        n_jobs=min(SCENIC_WORKERS, os.cpu_count() or 1),
        inplace=False,
        verbose=True,
    )
    if result is None or result.empty:
        raise RuntimeError("LIANA 未产生候选")
    gates = result.apply(lambda row: _complex_gate(row, detection, expr_prop), axis=1, result_type="expand")
    gates.columns = ["all_subunits_pass", "ligand_subunit_expr_props", "receptor_subunit_expr_props"]
    result = pd.concat([result.reset_index(drop=True), gates], axis=1)
    return result[result["all_subunits_pass"]].copy()


def _balanced_communication_data(data: ad.AnnData) -> tuple[dict[str, ad.AnnData], pd.DataFrame]:
    libraries = sorted(data.obs["library_id"].astype(str).unique())
    known = data.obs["cell_type"].astype(str) != "unknown"
    sizes = data.obs.loc[known].groupby(["library_id", "cell_type"], observed=True).size()
    common_types = sorted(set.intersection(*[set(sizes.loc[library].index.astype(str)) for library in libraries]))
    rng = np.random.default_rng(20260905)
    output: dict[str, ad.AnnData] = {}
    sampling_rows: list[dict[str, Any]] = []
    obs_library = data.obs["library_id"].astype(str).to_numpy()
    obs_type = data.obs["cell_type"].astype(str).to_numpy()
    for library in libraries:
        selected: list[int] = []
        for cell_type in common_types:
            target = min(500, min(int(sizes.loc[(item, cell_type)]) for item in libraries))
            candidates = np.flatnonzero((obs_library == library) & (obs_type == cell_type))
            chosen = np.sort(rng.choice(candidates, target, replace=False))
            selected.extend(chosen.tolist())
            sampling_rows.append({"library_id": library, "cell_type": cell_type, "available_cells": len(candidates), "selected_cells": target, "seed": 20260905})
        output[library] = data[np.asarray(selected, dtype=int)].copy()
    return output, pd.DataFrame(sampling_rows)


def _annotate_pathways(frame: pd.DataFrame, pathway_db: pd.DataFrame) -> pd.DataFrame:
    lookup: dict[tuple[str, str], list[str]] = {}
    for row in pathway_db.itertuples(index=False):
        lookup.setdefault((str(row.ligand_complex), str(row.receptor_complex)), []).append(str(row.pathway_name))
    output = frame.copy()
    output["pathway"] = [
        ";".join(sorted(set(lookup.get((str(ligand), str(receptor)), [])))) or "unannotated_mouseconsensus"
        for ligand, receptor in zip(output["ligand_complex"], output["receptor_complex"], strict=True)
    ]
    return output


def _plot_communication_bubbles(full: pd.DataFrame, figure_dir: Path, source_dir: Path) -> list[dict[str, str]]:
    figures: list[dict[str, str]] = []
    ranked = full.sort_values(["library_id", "magnitude_rank"]).groupby("library_id", observed=True).head(20).copy()
    ranked["interaction"] = ranked["ligand_complex"].astype(str) + " → " + ranked["receptor_complex"].astype(str)
    ranked["cell_pair"] = ranked["source"].astype(str) + " → " + ranked["target"].astype(str)
    ranked["plot_strength"] = -np.log10(ranked["magnitude_rank"].clip(lower=1e-12))
    stem = "communication_lr_bubble"
    atomic_tsv(source_dir / f"{stem}.tsv", ranked)
    libraries = list(dict.fromkeys(ranked["library_id"].astype(str)))
    fig, axes = plt.subplots(1, len(libraries), figsize=(6.8 * len(libraries), 7.0), squeeze=False)
    for ax, library in zip(axes[0], libraries, strict=True):
        sub = ranked[ranked["library_id"].astype(str) == library].copy().reset_index(drop=True)
        pairs = list(dict.fromkeys(sub["cell_pair"]))
        interactions = list(reversed(list(dict.fromkeys(sub["interaction"]))))
        xmap, ymap = {x: i for i, x in enumerate(pairs)}, {x: i for i, x in enumerate(interactions)}
        points = ax.scatter(
            sub["cell_pair"].map(xmap), sub["interaction"].map(ymap),
            s=25 + 32 * sub["plot_strength"].clip(upper=8),
            c=sub["magnitude_rank"], cmap="viridis_r", vmin=0, vmax=max(0.05, float(sub["magnitude_rank"].max())),
            edgecolor="#333333", linewidth=0.3,
        )
        ax.set_xticks(range(len(pairs)), pairs, rotation=55, ha="right", fontsize=7)
        ax.set_yticks(range(len(interactions)), interactions, fontsize=7)
        ax.set_title(library)
        ax.set_xlabel("Sender → receiver")
        fig.colorbar(points, ax=ax, label="LIANA consensus magnitude rank")
    fig.suptitle("Per-library ligand–receptor candidates (descriptive)")
    save_figure(fig, figure_dir / stem, 300)
    figures.append({"title": "各文库配体—受体候选气泡图", "png": f"figures/{stem}.png", "pdf": f"figures/{stem}.pdf", "source": f"source_data/{stem}.tsv", "interpretation": "LIANA 共识 magnitude rank 越低表示文库内候选优先级越高；该图不是 PFOS 与 WT 的条件检验。"})

    annotated = full[full["pathway"] != "unannotated_mouseconsensus"].copy()
    annotated["pathway"] = annotated["pathway"].str.split(";")
    annotated = annotated.explode("pathway")
    pathway = (
        annotated.groupby(["library_id", "source", "target", "pathway"], observed=True)
        .agg(candidate_interactions=("ligand_complex", "size"), best_magnitude_rank=("magnitude_rank", "min"))
        .reset_index()
    )
    pathway["cell_pair"] = pathway["source"].astype(str) + " → " + pathway["target"].astype(str)
    pathway["plot_strength"] = -np.log10(pathway["best_magnitude_rank"].clip(lower=1e-12))
    pathway = pathway.sort_values(["library_id", "best_magnitude_rank"]).groupby("library_id", observed=True).head(25)
    stem = "communication_pathway_bubble"
    atomic_tsv(source_dir / f"{stem}.tsv", pathway)
    fig, axes = plt.subplots(1, len(libraries), figsize=(6.8 * len(libraries), 7.0), squeeze=False)
    for ax, library in zip(axes[0], libraries, strict=True):
        sub = pathway[pathway["library_id"].astype(str) == library].copy()
        pairs = list(dict.fromkeys(sub["cell_pair"]))
        pathways = list(reversed(list(dict.fromkeys(sub["pathway"]))))
        xmap, ymap = {x: i for i, x in enumerate(pairs)}, {x: i for i, x in enumerate(pathways)}
        points = ax.scatter(sub["cell_pair"].map(xmap), sub["pathway"].map(ymap), s=30 + 18 * sub["candidate_interactions"], c=sub["best_magnitude_rank"], cmap="viridis_r", vmin=0, vmax=0.05, edgecolor="#333333", linewidth=0.3)
        ax.set_xticks(range(len(pairs)), pairs, rotation=55, ha="right", fontsize=7)
        ax.set_yticks(range(len(pathways)), pathways, fontsize=7)
        ax.set_title(library)
        ax.set_xlabel("Sender → receiver")
        fig.colorbar(points, ax=ax, label="Best magnitude rank")
    fig.suptitle("CellChatDB.mouse pathway labels for LIANA candidates")
    save_figure(fig, figure_dir / stem, 300)
    figures.append({"title": "各文库通路标注的通讯候选", "png": f"figures/{stem}.png", "pdf": f"figures/{stem}.pdf", "source": f"source_data/{stem}.tsv", "interpretation": "通路名称来自 CellChatDB.mouse，候选分数来自 LIANA；结果仅用于提出可验证假说。"})
    return figures


def _plot_communication_network(summary: pd.DataFrame, figure_dir: Path, source_dir: Path) -> dict[str, str]:
    stem = "communication_network"
    atomic_tsv(source_dir / f"{stem}.tsv", summary)
    libraries = list(dict.fromkeys(summary["library_id"].astype(str)))
    cell_types = sorted(set(summary["source"].astype(str)) | set(summary["target"].astype(str)))
    theta = np.linspace(0, 2 * np.pi, len(cell_types), endpoint=False)
    pos = {cell: np.array([np.cos(angle), np.sin(angle)]) for cell, angle in zip(cell_types, theta, strict=True)}
    fig, axes = plt.subplots(1, len(libraries), figsize=(6.5 * len(libraries), 6.0), squeeze=False)
    max_weight = max(float(summary["priority_candidates"].max()), 1.0)
    for ax, library in zip(axes[0], libraries, strict=True):
        sub = summary[summary["library_id"].astype(str) == library]
        for row in sub.itertuples(index=False):
            if row.priority_candidates <= 0:
                continue
            start, end = pos[str(row.source)], pos[str(row.target)]
            ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "-|>", "color": "#555555", "alpha": 0.25 + 0.65 * row.priority_candidates / max_weight, "lw": 0.5 + 4 * row.priority_candidates / max_weight, "shrinkA": 17, "shrinkB": 17})
        for cell in cell_types:
            xy = pos[cell]
            ax.scatter(*xy, s=520, color=CELL_TYPE_COLORS.get(cell, "#999999"), edgecolor="white", linewidth=1.2, zorder=3)
            ax.text(*xy, cell, ha="center", va="center", fontsize=8, zorder=4)
        ax.set_title(library)
        ax.set_xlim(-1.35, 1.35)
        ax.set_ylim(-1.35, 1.35)
        ax.axis("off")
    fig.suptitle("Priority LR candidates between sender and receiver cell types")
    save_figure(fig, figure_dir / stem, 300)
    return {"title": "各文库通讯候选网络", "png": f"figures/{stem}.png", "pdf": f"figures/{stem}.pdf", "source": f"source_data/{stem}.tsv", "interpretation": "边宽表示 LIANA magnitude rank 前 5% 的候选数量，不代表互作概率或差异信号。"}


def _draw_communication_sensitivity(sensitivity: pd.DataFrame, target_stem: Path, *, library_conditions: dict[str, str]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
    scenarios = list(dict.fromkeys(sensitivity["scenario"].astype(str)))
    libraries = list(dict.fromkeys(sensitivity["library_id"].astype(str)))
    x = np.arange(len(scenarios)); width = 0.36
    for index, library in enumerate(libraries):
        subset = sensitivity[sensitivity["library_id"].astype(str) == library].set_index("scenario").reindex(scenarios)
        color = COLORS[library_conditions[library]]
        axes[0].bar(x + (index - 0.5) * width, subset["priority_candidates"], width, label=library, color=color)
        axes[1].bar(x + (index - 0.5) * width, subset["priority_overlap_with_primary"], width, label=library, color=color)
    for ax in axes:
        ax.set_xticks(x, scenarios, rotation=35, ha="right"); ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("Priority LR candidates"); axes[1].set_ylabel("Jaccard overlap with primary")
    axes[1].set_ylim(0, 1); fig.suptitle("Technical sensitivity of communication candidates (descriptive)")
    save_figure(fig, target_stem, 300)


def run_communication(data: ad.AnnData, output_dir: Path, input_h5ad: Path, *, rscript: Path) -> dict[str, Any]:
    import liana as li

    reusable = _module_reusable(output_dir, "communication_manifest.json", "communication_execution_inputs.json")
    if reusable is not None:
        return reusable
    table_dir, figure_dir, source_dir = (output_dir / name for name in ("tables", "figures", "source_data"))
    for path in (table_dir, figure_dir, source_dir):
        path.mkdir(parents=True, exist_ok=True)
    resource_path = source_dir / "liana_mouseconsensus_resource.tsv"
    pathway_path = source_dir / "cellchatdb_mouse_pathway_annotations.tsv"
    resource = _write_liana_resource_snapshot(resource_path)
    pathway_db = _write_cellchat_pathway_snapshot(pathway_path, rscript)
    module_source = Path(__file__)
    execution_path = output_dir / "communication_execution_inputs.json"
    inputs = _write_execution_inputs(execution_path, [input_h5ad, module_source, resource_path, pathway_path])

    excluded = (
        data.obs.groupby(["library_id", "condition"], observed=True)["cell_type"]
        .agg(total_cells="size", unknown_cells=lambda x: int((x.astype(str) == "unknown").sum()))
        .reset_index()
    )
    excluded["included_known_cells"] = excluded["total_cells"] - excluded["unknown_cells"]
    excluded["exclusion_rule"] = "cell_type == unknown excluded from primary communication analysis"
    atomic_tsv(table_dir / "communication_excluded_cells.tsv", excluded)

    frames: list[pd.DataFrame] = []
    for library in sorted(data.obs["library_id"].astype(str).unique()):
        mask = (data.obs["library_id"].astype(str) == library) & (data.obs["cell_type"].astype(str) != "unknown")
        lib = data[mask].copy()
        result = _liana_library(lib, resource, expr_prop=EXPR_PROP, n_perms=LIANA_PERMUTATIONS, seed=LIANA_SEED)
        result.insert(0, "library_id", library)
        result.insert(1, "condition", str(lib.obs["condition"].astype(str).iloc[0]))
        result["priority_candidate"] = result["magnitude_rank"] <= PRIORITY_MAGNITUDE_RANK
        result["candidate_rule"] = f"magnitude_rank <= {PRIORITY_MAGNITUDE_RANK}; no between-condition test"
        result["inferential_scope"] = "within-library co-expression candidate; no PFOS-versus-WT P/FDR"
        frames.append(result)
    full = _annotate_pathways(pd.concat(frames, ignore_index=True), pathway_db)
    full_path = table_dir / "communication_lr_full.tsv"
    atomic_tsv(full_path, full)

    summary = (
        full.groupby(["library_id", "condition", "source", "target"], observed=True)
        .agg(eligible_lr=("ligand_complex", "size"), priority_candidates=("priority_candidate", "sum"), best_magnitude_rank=("magnitude_rank", "min"), median_magnitude_rank=("magnitude_rank", "median"))
        .reset_index()
    )
    summary["scope"] = "descriptive within-library summary"
    atomic_tsv(table_dir / "communication_sender_receiver_summary.tsv", summary)

    candidates = full[full["priority_candidate"]].copy()
    key = ["source", "target", "ligand_complex", "receptor_complex"]
    catalog = (
        candidates.groupby(key, observed=True)
        .agg(libraries=("library_id", lambda x: ";".join(sorted(set(map(str, x))))), conditions=("condition", lambda x: ";".join(sorted(set(map(str, x))))), n_libraries=("library_id", "nunique"), best_magnitude_rank=("magnitude_rank", "min"), pathways=("pathway", lambda x: ";".join(sorted(set(map(str, x))))))
        .reset_index()
    )
    catalog["catalog_status"] = np.where(catalog["n_libraries"] == data.obs["library_id"].nunique(), "shared", "library-specific")
    catalog["scope"] = "top-5% magnitude-rank catalog; presence/absence only; no differential test"
    atomic_tsv(table_dir / "communication_shared_unique_candidates.tsv", catalog)
    pathway_summary = (
        full.assign(pathway=full["pathway"].str.split(";")).explode("pathway")
        .groupby(["library_id", "condition", "source", "target", "pathway"], observed=True)
        .agg(eligible_lr=("ligand_complex", "size"), priority_candidates=("priority_candidate", "sum"), best_magnitude_rank=("magnitude_rank", "min"))
        .reset_index()
    )
    pathway_summary["scope"] = "CellChatDB.mouse label over LIANA candidates; descriptive"
    atomic_tsv(table_dir / "communication_pathway_summary.tsv", pathway_summary)

    sensitivity_frames: list[pd.DataFrame] = []
    primary_keys: dict[str, set[tuple[str, str, str, str]]] = {}
    key_columns = ["source", "target", "ligand_complex", "receptor_complex"]
    for library, group in full.groupby("library_id", observed=True):
        priority = group[group["priority_candidate"]]
        primary_keys[str(library)] = set(priority[key_columns].itertuples(index=False, name=None))
        sensitivity_frames.append(pd.DataFrame([{
            "library_id": library, "scenario": "expr0.10_allcells", "expr_prop": EXPR_PROP,
            "cells": int((data.obs["library_id"].astype(str) == str(library)).sum()), "eligible_lr": len(group),
            "priority_candidates": int(group["priority_candidate"].sum()), "priority_overlap_with_primary": 1.0,
            "permutations": LIANA_PERMUTATIONS, "scope": "technical sensitivity; no condition test",
        }]))
    sensitivity_lr: list[pd.DataFrame] = []
    for threshold in (0.05, 0.20):
        for library in sorted(data.obs["library_id"].astype(str).unique()):
            mask = (data.obs["library_id"].astype(str) == library) & (data.obs["cell_type"].astype(str) != "unknown")
            result = _liana_library(data[mask].copy(), resource, expr_prop=threshold, n_perms=100, seed=LIANA_SEED)
            result["priority_candidate"] = result["magnitude_rank"] <= PRIORITY_MAGNITUDE_RANK
            current = set(result.loc[result["priority_candidate"], key_columns].itertuples(index=False, name=None))
            union = current | primary_keys[library]
            sensitivity_frames.append(pd.DataFrame([{
                "library_id": library, "scenario": f"expr{threshold:.2f}_allcells", "expr_prop": threshold,
                "cells": int(mask.sum()), "eligible_lr": len(result), "priority_candidates": int(result["priority_candidate"].sum()),
                "priority_overlap_with_primary": len(current & primary_keys[library]) / len(union) if union else 1.0,
                "permutations": 100, "scope": "technical sensitivity; no condition test",
            }]))
            result.insert(0, "library_id", library); result.insert(1, "scenario", f"expr{threshold:.2f}_allcells")
            sensitivity_lr.append(result)
    balanced, sampling = _balanced_communication_data(data)
    atomic_tsv(table_dir / "communication_sensitivity_sampling.tsv", sampling)
    for library, lib in balanced.items():
        result = _liana_library(lib, resource, expr_prop=EXPR_PROP, n_perms=100, seed=LIANA_SEED)
        result["priority_candidate"] = result["magnitude_rank"] <= PRIORITY_MAGNITUDE_RANK
        current = set(result.loc[result["priority_candidate"], key_columns].itertuples(index=False, name=None))
        union = current | primary_keys[library]
        sensitivity_frames.append(pd.DataFrame([{
            "library_id": library, "scenario": "expr0.10_equalcells", "expr_prop": EXPR_PROP,
            "cells": lib.n_obs, "eligible_lr": len(result), "priority_candidates": int(result["priority_candidate"].sum()),
            "priority_overlap_with_primary": len(current & primary_keys[library]) / len(union) if union else 1.0,
            "permutations": 100, "scope": "equal-cell technical sensitivity; no condition test",
        }]))
        result.insert(0, "library_id", library); result.insert(1, "scenario", "expr0.10_equalcells")
        sensitivity_lr.append(result)
    sensitivity = pd.concat(sensitivity_frames, ignore_index=True)
    atomic_tsv(table_dir / "communication_sensitivity.tsv", sensitivity)
    pd.concat(sensitivity_lr, ignore_index=True).to_csv(table_dir / "communication_sensitivity_lr.tsv.gz", sep="\t", index=False, compression="gzip")
    stem = "communication_sensitivity"
    atomic_tsv(source_dir / f"{stem}.tsv", sensitivity)
    _draw_communication_sensitivity(sensitivity, figure_dir / stem, library_conditions=data.obs[["library_id", "condition"]].drop_duplicates().set_index("library_id")["condition"].astype(str).to_dict())

    figures = _plot_communication_bubbles(full, figure_dir, source_dir)
    figures.append(_plot_communication_network(summary, figure_dir, source_dir))
    figures.append({"title": "通讯候选表达阈值与等细胞数敏感性", "png": f"figures/{stem}.png", "pdf": f"figures/{stem}.pdf", "source": f"source_data/{stem}.tsv", "interpretation": "比较 expr_prop=0.05/0.10/0.20 及按文库×细胞类型等细胞数后的候选数量与主分析重叠；属于技术敏感性，不是条件检验。"})
    if not _inputs_unchanged(inputs):
        raise RuntimeError("communication 输入或资源在执行期间发生变化")
    metrics = {
        "libraries": int(full["library_id"].nunique()),
        "resource_pairs": int(len(resource)),
        "eligible_lr_rows": int(len(full)),
        "priority_lr_rows": int(full["priority_candidate"].sum()),
        "shared_priority_candidates": int((catalog["catalog_status"] == "shared").sum()),
        "unknown_cells_excluded": int(excluded["unknown_cells"].sum()),
        "sensitivity_scenarios": int(sensitivity["scenario"].nunique()),
    }
    manifest = {
        "status": "completed",
        "method": f"LIANA+ {li.__version__} rank_aggregate per library; {LIANA_RESOURCE}; all complex subunits expr_prop>={EXPR_PROP}; min_cells>={MIN_CELLS}; CellChatDB.mouse pathway labels only",
        "metrics": metrics,
        "figures": figures,
        "tables": [
            "tables/communication_lr_full.tsv", "tables/communication_sender_receiver_summary.tsv",
            "tables/communication_pathway_summary.tsv", "tables/communication_shared_unique_candidates.tsv",
            "tables/communication_excluded_cells.tsv", "tables/communication_sensitivity.tsv",
            "tables/communication_sensitivity_sampling.tsv", "tables/communication_sensitivity_lr.tsv.gz",
        ],
        "resources": {"liana_mouseconsensus_sha256": sha256(resource_path), "cellchatdb_mouse_pathways_sha256": sha256(pathway_path)},
        "completed_at": _utc_now(),
    }
    atomic_json(output_dir / "communication_manifest.json", manifest)
    return manifest


def _prepare_scenic_inputs(data: ad.AnnData, source_dir: Path) -> tuple[Path, Path, pd.DataFrame, list[str]]:
    import loompy

    known = data.obs["cell_type"].astype(str) != "unknown"
    counts = data.obs.loc[known].groupby(["library_id", "cell_type"], observed=True).size()
    eligible = counts[counts >= MIN_CELLS]
    libraries = sorted(data.obs["library_id"].astype(str).unique())
    cell_types = sorted(set.intersection(*[set(eligible.loc[library].index.astype(str)) for library in libraries]))
    strata = [(library, cell_type) for library in libraries for cell_type in cell_types]
    per_stratum = min(SCENIC_MAX_CELLS_PER_STRATUM, min(int(counts.loc[stratum]) for stratum in strata))
    rng = np.random.default_rng(20260905)
    selected: list[int] = []
    rows: list[dict[str, Any]] = []
    obs_library = data.obs["library_id"].astype(str).to_numpy()
    obs_type = data.obs["cell_type"].astype(str).to_numpy()
    for library, cell_type in strata:
        candidates = np.flatnonzero((obs_library == library) & (obs_type == cell_type))
        chosen = np.sort(rng.choice(candidates, size=per_stratum, replace=False))
        selected.extend(chosen.tolist())
        rows.extend({"barcode": str(data.obs_names[index]), "library_id": library, "cell_type": cell_type, "selected": True, "balance_cells_per_stratum": per_stratum} for index in chosen)
    selected_array = np.asarray(selected, dtype=int)
    train_counts = data.layers["counts"][selected_array].tocsr()
    min_samples = max(1, math.ceil(0.01 * len(selected_array)))
    min_counts = max(1, math.ceil(3 * 0.01 * len(selected_array)))
    detected = np.asarray((train_counts > 0).sum(axis=0)).ravel()
    totals = np.asarray(train_counts.sum(axis=0)).ravel()
    gene_mask = (detected >= min_samples) & (totals >= min_counts)
    genes = data.var_names[gene_mask].astype(str).tolist()
    if len(genes) < 1000:
        raise RuntimeError(f"SCENIC 过滤后仅 {len(genes)} genes")
    train_path = source_dir / "scenic_training_expression.loom"
    all_path = source_dir / "scenic_all_cells_expression.loom"
    gene_idx = np.flatnonzero(gene_mask)
    loompy.create(str(train_path), train_counts[:, gene_idx].T.tocsc(), {"Gene": np.asarray(genes)}, {"CellID": data.obs_names[selected_array].astype(str).to_numpy()})
    loompy.create(str(all_path), data.layers["counts"][:, gene_idx].T.tocsc(), {"Gene": np.asarray(genes)}, {"CellID": data.obs_names.astype(str).to_numpy()})
    barcodes = pd.DataFrame(rows).sort_values(["library_id", "cell_type", "barcode"])
    atomic_tsv(source_dir / "scenic_training_barcodes.tsv", barcodes)
    filter_table = pd.DataFrame({"gene": data.var_names.astype(str), "training_detected_cells": detected, "training_total_umi": totals, "retained": gene_mask, "min_detected_cells": min_samples, "min_total_umi": min_counts})
    atomic_tsv(source_dir / "scenic_gene_filter.tsv", filter_table)
    return train_path, all_path, barcodes, genes


def _regulon_specificity(auc: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for library in sorted(metadata["library_id"].astype(str).unique()):
        lib_index = metadata.index[metadata["library_id"].astype(str) == library]
        values = auc.loc[lib_index].to_numpy(dtype=float)
        for cell_type in sorted(metadata.loc[lib_index, "cell_type"].astype(str).unique()):
            label = (metadata.loc[lib_index, "cell_type"].astype(str).to_numpy() == cell_type).astype(float)
            label /= label.sum()
            for col, regulon in enumerate(auc.columns):
                profile = np.maximum(values[:, col], 0)
                profile = profile / profile.sum() if profile.sum() > 0 else np.full(len(profile), 1 / len(profile))
                rss = 1.0 - float(jensenshannon(profile, label, base=2.0))
                group_values = values[label > 0, col]
                rows.append({"library_id": library, "condition": str(metadata.loc[lib_index, "condition"].astype(str).iloc[0]), "cell_type": cell_type, "regulon": str(regulon), "cells": int((label > 0).sum()), "auc_mean": float(np.mean(group_values)), "auc_median": float(np.median(group_values)), "auc_sd": float(np.std(group_values)), "rss": rss, "scope": "within-library cell-type specificity; no condition test"})
    return pd.DataFrame(rows)


def _draw_regulon_activity(heat: pd.DataFrame) -> plt.Figure:
    top = heat.groupby("regulon", observed=True)["rss"].max().nlargest(20).index
    pivot = heat.pivot(index="regulon", columns="stratum", values="auc_zscore").reindex(index=top)
    fig, ax = plt.subplots(figsize=(max(9.0, 0.55 * pivot.shape[1]), max(6.0, 0.30 * pivot.shape[0])))
    image = ax.imshow(pivot.to_numpy(), aspect="auto", cmap=DIVERGING_CMAP, vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(pivot.shape[1]), pivot.columns, rotation=55, ha="right", fontsize=7)
    ax.set_yticks(range(pivot.shape[0]), pivot.index, fontsize=7)
    ax.set_title("Pooled-consensus regulon activity")
    fig.colorbar(image, ax=ax, label="Within-regulon AUC z-score")
    return fig


def _plot_regulons(summary: pd.DataFrame, edges: pd.DataFrame, figure_dir: Path, source_dir: Path) -> list[dict[str, str]]:
    figures: list[dict[str, str]] = []
    top = summary.groupby("regulon", observed=True)["rss"].max().nlargest(20).index
    heat = summary[summary["regulon"].isin(top)].copy()
    heat["stratum"] = heat["library_id"].astype(str) + " | " + heat["cell_type"].astype(str)
    heat["auc_zscore"] = heat.groupby("regulon", observed=True)["auc_mean"].transform(lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0.0)
    stem = "regulon_activity_heatmap"
    atomic_tsv(source_dir / f"{stem}.tsv", heat)
    fig = _draw_regulon_activity(heat)
    save_figure(fig, figure_dir / stem, 300)
    figures.append({"title": "Pooled consensus regulon 活性热图", "png": f"figures/{stem}.png", "pdf": f"figures/{stem}.pdf", "source": f"source_data/{stem}.tsv", "interpretation": "同一个 pooled consensus 网络投射到全部细胞；各列为文库×细胞类型的描述性汇总。"})

    stability = edges.groupby("TF", observed=True).agg(stable_edges=("target", "size"), mean_seed_support=("seed_support", "mean"), median_importance=("importance_median", "median")).reset_index().sort_values(["stable_edges", "mean_seed_support"], ascending=False).head(30)
    stem = "regulon_tf_stability"
    atomic_tsv(source_dir / f"{stem}.tsv", stability)
    fig, ax = plt.subplots(figsize=(8.5, max(5.0, 0.22 * len(stability))))
    ordered = stability.sort_values("stable_edges")
    bars = ax.barh(ordered["TF"], ordered["stable_edges"], color="#4C78A8")
    for bar, support in zip(bars, ordered["mean_seed_support"], strict=True):
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2, f" {support:.2f}/3", va="center", fontsize=7)
    ax.set_xlabel("Edges retained in at least 2 of 3 seed top-100 lists")
    ax.set_ylabel("Reference-list candidate regulator")
    ax.set_title("GRNBoost2 random-seed stability")
    save_figure(fig, figure_dir / stem, 300)
    figures.append({"title": "GRNBoost2 转录因子稳定性", "png": f"figures/{stem}.png", "pdf": f"figures/{stem}.pdf", "source": f"source_data/{stem}.tsv", "interpretation": "稳定性指每个 TF 的前 100 条边在三个随机种子之间的重复出现，不是生物学重复。"})

    top_tf = summary.groupby("regulon", observed=True)["rss"].max().nlargest(8).index.astype(str)
    network = edges[edges["TF"].astype(str).isin([x.replace("(+)", "").replace("(-)", "") for x in top_tf])].sort_values(["TF", "importance_median"], ascending=[True, False]).groupby("TF", observed=True).head(5)
    stem = "regulon_top_network"
    atomic_tsv(source_dir / f"{stem}.tsv", network)
    fig, ax = plt.subplots(figsize=(10.0, 7.0))
    tfs = sorted(network["TF"].astype(str).unique())
    targets = sorted(network["target"].astype(str).unique())
    y_tf = {tf: y for tf, y in zip(tfs, np.linspace(0.08, 0.92, max(len(tfs), 1)), strict=True)}
    y_target = {gene: y for gene, y in zip(targets, np.linspace(0.03, 0.90, max(len(targets), 1)), strict=True)}
    for row in network.itertuples(index=False):
        ax.plot([0.2, 0.8], [y_tf[str(row.TF)], y_target[str(row.target)]], color="#888888", alpha=0.25 + 0.25 * row.seed_support, lw=0.5 + 2 * row.importance_median / max(float(network["importance_median"].max()), 1e-12))
    for tf, y in y_tf.items():
        ax.scatter(0.2, y, s=260, color="#E69F00", zorder=3); ax.text(0.2, y, tf, ha="center", va="center", fontsize=7)
    for gene, y in y_target.items():
        ax.scatter(0.8, y, s=170, color="#56B4E9", zorder=3); ax.text(0.8, y, gene, ha="center", va="center", fontsize=7)
    ax.text(0.2, 0.985, "Candidate regulator", ha="center", va="top", weight="bold")
    ax.text(0.8, 0.985, "Motif-supported target", ha="center", va="top", weight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    fig.suptitle("High-priority pooled-consensus regulon edges", y=0.995)
    fig.subplots_adjust(top=0.91)
    save_figure(fig, figure_dir / stem, 300)
    figures.append({"title": "高优先级 pooled-consensus 调控网络", "png": f"figures/{stem}.png", "pdf": f"figures/{stem}.pdf", "source": f"source_data/{stem}.tsv", "interpretation": "图中边在多个 GRNBoost2 种子中重复并通过 motif 剪枝，但仍只是计算调控候选。"})
    return figures


def run_regulon(data: ad.AnnData, output_dir: Path, resources_dir: Path, input_h5ad: Path) -> dict[str, Any]:
    reusable = _module_reusable(output_dir, "regulon_manifest.json", "regulon_execution_inputs.json")
    if reusable is not None:
        return reusable
    table_dir, figure_dir, source_dir = (output_dir / name for name in ("tables", "figures", "source_data"))
    for path in (table_dir, figure_dir, source_dir):
        path.mkdir(parents=True, exist_ok=True)
    project_root = Path(__file__).resolve().parents[2]
    runner = project_root / "scripts" / "extension_scenic_runner.py"
    lock = project_root / "scripts" / "extension_scenic_requirements.lock"
    scenic_python = project_root / ".venv-pyscenic" / "bin" / "python"
    if not scenic_python.is_file():
        raise FileNotFoundError(f"独立 pySCENIC 环境不存在: {scenic_python}")
    rankings = sorted(resources_dir.glob("mm10_*genes_vs_motifs.rankings.feather"))
    motif = resources_dir / "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
    tf_list = resources_dir / "allTFs_mm.txt"
    if len(rankings) != 2 or not motif.is_file() or not tf_list.is_file():
        raise FileNotFoundError("pySCENIC mouse v10 资源不完整")
    train_loom, all_loom, barcodes, genes = _prepare_scenic_inputs(data, source_dir)
    execution_path = output_dir / "regulon_execution_inputs.json"
    inputs = _write_execution_inputs(execution_path, [input_h5ad, Path(__file__), runner, lock, tf_list, motif, *rankings, train_loom, all_loom])
    work = output_dir / "source_data" / "scenic_work"
    work.mkdir(exist_ok=True)
    command = [
        str(scenic_python), str(runner), "--training-loom", str(train_loom), "--all-cells-loom", str(all_loom),
        "--tf-list", str(tf_list), "--motif-annotations", str(motif), "--rankings", *map(str, rankings),
        "--output-dir", str(work), "--workers", str(SCENIC_WORKERS), "--top-per-tf", str(SCENIC_EDGE_TOP_PER_TF),
        "--min-seeds", str(SCENIC_EDGE_MIN_SEEDS), "--seeds", *map(str, SCENIC_SEEDS),
    ]
    log_path = output_dir / "scenic_execution.log"
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=project_root, stdout=log, stderr=subprocess.STDOUT, text=True, check=True)
    stable_edges = pd.read_csv(work / "stable_edges.tsv", sep="\t")
    regulons = pd.read_csv(work / "motif_supported_regulons.tsv", sep="\t")
    auc = pd.read_csv(work / "regulon_auc_all_cells.tsv.gz", sep="\t", index_col=0)
    if set(data.obs_names) != set(auc.index.astype(str)):
        raise RuntimeError("AUCell barcode 与正式 H5AD 不一致")
    auc = auc.reindex(data.obs_names.astype(str))
    metadata = data.obs[["library_id", "condition", "cell_type"]].copy()
    metadata.index = metadata.index.astype(str)
    summary = _regulon_specificity(auc, metadata)
    atomic_tsv(table_dir / "regulon_seed_stable_edges.tsv", stable_edges)
    atomic_tsv(table_dir / "regulon_motif_supported_targets.tsv", regulons)
    atomic_tsv(table_dir / "regulon_library_celltype_auc_rss.tsv", summary)
    auc.to_csv(table_dir / "regulon_auc_all_cells.tsv.gz", sep="\t", compression="gzip")
    intersection = json.loads((work / "resource_intersection.json").read_text(encoding="utf-8"))
    atomic_json(table_dir / "regulon_resource_intersection.json", intersection)
    figures = _plot_regulons(summary, stable_edges.merge(regulons[["TF", "target"]].drop_duplicates(), on=["TF", "target"], how="inner"), figure_dir, source_dir)
    if not _inputs_unchanged(inputs):
        raise RuntimeError("regulon 输入或资源在执行期间发生变化")
    metrics = {
        "training_cells": int(len(barcodes)), "balanced_cells_per_library_celltype": int(barcodes["balance_cells_per_stratum"].iloc[0]),
        "training_genes": int(len(genes)), "grn_seeds": len(SCENIC_SEEDS), "stable_edges": int(len(stable_edges)),
        "motif_supported_edges": int(len(regulons)), "regulons": int(auc.shape[1]), "scored_cells": int(auc.shape[0]),
        **intersection,
    }
    manifest = {
        "status": "completed",
        "method": "pooled library×cell-type balanced pySCENIC consensus: three seeded GRNBoost2 top-100-per-TF edge lists, edges recurring in >=2 seeds, mouse-v10 cisTarget motif pruning, AUCell on all cells",
        "metrics": metrics,
        "figures": figures,
        "tables": ["tables/regulon_seed_stable_edges.tsv", "tables/regulon_motif_supported_targets.tsv", "tables/regulon_library_celltype_auc_rss.tsv", "tables/regulon_auc_all_cells.tsv.gz", "tables/regulon_resource_intersection.json"],
        "resources": {path.name: sha256(path) for path in [tf_list, motif, *rankings]},
        "completed_at": _utc_now(),
    }
    atomic_json(output_dir / "regulon_manifest.json", manifest)
    return manifest


def _validate_figure_contract(output_dir: Path, figures: list[dict[str, str]]) -> None:
    for figure in figures:
        for key in ("png", "pdf", "source"):
            path = output_dir / figure[key]
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(f"图形契约缺失: {path}")
        if not figure.get("title") or not figure.get("interpretation"):
            raise RuntimeError("图形 manifest 缺少 title/interpretation")


def run_networks(input_h5ad: Path, output_dir: Path, resources_dir: Path, *, rscript: Path) -> dict[str, Any]:
    input_h5ad, output_dir, resources_dir = Path(input_h5ad), Path(output_dir), Path(resources_dir)
    if not input_h5ad.is_file():
        raise FileNotFoundError(input_h5ad)
    output_dir.mkdir(parents=True, exist_ok=True)
    data = ad.read_h5ad(input_h5ad)
    _validate_adata(data)
    communication = run_communication(data, output_dir, input_h5ad, rscript=rscript)
    regulon = run_regulon(data, output_dir, resources_dir, input_h5ad)
    figures = communication["figures"] + regulon["figures"]
    _validate_figure_contract(output_dir, figures)
    manifest = {
        "status": "completed",
        "method": {"communication": communication["method"], "regulon": regulon["method"]},
        "metrics": {"communication": communication["metrics"], "regulon": regulon["metrics"]},
        "figures": figures,
        "tables": communication["tables"] + regulon["tables"],
        "claim_boundary": "Exploratory candidate generation from one library per condition; no condition P/FDR, differential communication, differential regulon, causality, or functional interaction claim.",
        "completed_at": _utc_now(),
    }
    stage_inputs = []
    for name in ("communication_execution_inputs.json", "regulon_execution_inputs.json"):
        stage_payload = json.loads((output_dir / name).read_text(encoding="utf-8"))
        stage_inputs.extend(stage_payload["files"])
    unique_inputs = {item["path"]: item for item in stage_inputs}
    combined_inputs = {"files": list(unique_inputs.values()), "started_at": min(
        json.loads((output_dir / name).read_text(encoding="utf-8"))["started_at"]
        for name in ("communication_execution_inputs.json", "regulon_execution_inputs.json")
    ), "completed_at": _utc_now(), "stage_code_hashes_verified": True}
    if not _inputs_unchanged(combined_inputs):
        raise RuntimeError("总网络模块输入/代码/资源哈希复核失败")
    atomic_json(output_dir / "execution_inputs.json", combined_inputs)
    atomic_json(output_dir / "report_manifest.json", manifest)
    project_root = Path(__file__).resolve().parents[2]
    finalizer = project_root / "scripts" / "extension_scenic_figure_finalize.py"
    subprocess.run([sys.executable, str(finalizer), str(output_dir)], cwd=project_root, check=True)
    final_inputs = json.loads((output_dir / "execution_inputs.json").read_text(encoding="utf-8"))
    if final_inputs.get("verification_status") != "PASS":
        raise RuntimeError("网络模块最终输入/代码/资源哈希未通过")
    return json.loads((output_dir / "report_manifest.json").read_text(encoding="utf-8"))

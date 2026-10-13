"""从已签名数值源重绘全部现役图；历史 runs 始终只读。"""
from __future__ import annotations

from datetime import datetime, timezone
import html
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from .plot_lifecycle import PlotRun

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import WorkflowConfig
from .plotting import categorical_embedding, save_figure, stacked_composition
from .runtime import atomic_json, atomic_text, sha256


def _stage(stem: str) -> str:
    if stem.startswith(("PFOS_Lib", "WT_Lib")) or stem in {"qc_retention_flow", "soupx_hormone_audit", "mt_threshold_sensitivity"}:
        return "01_qc"
    if stem in {"umap_library", "umap_condition", "umap_leiden"}:
        return "02_atlas"
    if stem in {"umap_cell_type", "marker_evidence_dotplot", "marker_evidence_heatmap", "annotation_evidence_map"}:
        return "03_annotation"
    if stem == "evidence_gate_matrix":
        return "00_audit"
    return "04_descriptive"


class DisplayContext:
    """Route reads to the sealed producer and writes to the unpublished result."""
    def __init__(self, config: WorkflowConfig, source: Path, output: Path):
        self.config, self.source, self.output = config, source, output

    def path(self, *parts: str) -> Path:
        if parts[0] in {"07_source_data", "08_figures"}:
            target = self.output / _stage(Path(parts[-1]).stem) / ("tables" if parts[0] == "07_source_data" else "plots")
            target.mkdir(parents=True, exist_ok=True)
            return target.joinpath(*parts[1:])
        return self.source.joinpath(*parts)


def _base_entries(source: Path) -> list[dict]:
    report = source / "09_reports/PFOS_islet_scrna_report_zh.html"
    entries = []
    for stem, body in re.findall(r'<figure id="fig-([^"]+)">(.*?)</figure>', report.read_text(), re.S):
        source_ref = re.search(r'href="([^"]+\.tsv)"', body).group(1)
        title = re.search(r'<h3>(.*?)</h3>', body).group(1)
        interpretation = re.sub(r'<[^>]+>', '', re.search(r'<figcaption>(.*?)<span', body, re.S).group(1)).strip()
        entries.append({"stem": stem, "title": html.unescape(title), "interpretation": html.unescape(interpretation),
                        "source_path": (report.parent / source_ref).resolve(), "stage": _stage(stem)})
    if not entries:
        raise RuntimeError("基础报告没有可追溯图件")
    return entries


def _verify_producer(path: Path, hash_file: Path, required: list[Path]) -> dict[str, str]:
    hashes = pd.read_csv(hash_file, sep="\t", dtype={"path": str, "sha256": str}).set_index("path")
    verified = {}
    for item in required:
        relative = str(item.relative_to(path))
        if relative not in hashes.index or not item.is_file():
            raise FileNotFoundError(f"绘图来源不在原运行签名中：{relative}")
        actual = sha256(item)
        if actual != hashes.at[relative, "sha256"]:
            raise RuntimeError(f"绘图来源与原运行签名不一致：{relative}")
        verified[str(item)] = actual
    return verified




def _render_base(config: WorkflowConfig, source: Path, output: Path, entries: list[dict]) -> None:
    from . import report_assets as assets
    from .qc import draw_qc_source
    context = DisplayContext(config, source, output)
    dpi = 300
    data = ad.read_h5ad(source / "04_annotation/annotated.h5ad")
    if set(data.obs["library_id"].astype(str)) != set(config.samples()["library_id"].astype(str)):
        raise RuntimeError("当前样本登记与冻结 H5AD 的文库身份不一致")
    for threshold_path in (source / "02_qc").glob("*/qc_thresholds.tsv"):
        threshold = pd.read_csv(threshold_path, sep="\t")["max_pct_mt"].iloc[0]
        if float(threshold) != float(config.data["qc"]["max_pct_mt"]):
            raise RuntimeError("当前 QC 标注阈值与冻结来源不一致")
    assets._plot_marker_evidence(context, data, dpi)
    for function in (assets._plot_annotation_evidence, assets._plot_focal_genes,
                     assets._plot_modules, assets._plot_balanced, assets._plot_evidence_gate,
                     assets._plot_qc_retention, assets._plot_soupx, assets._plot_mt_sensitivity):
        function(context, dpi)
    for entry in entries:
        stem = entry["stem"]
        target = output / entry["stage"] / "plots" / stem
        frame = pd.read_csv(entry["source_path"], sep="\t")
        if stem.startswith(("PFOS_Lib", "WT_Lib")):
            library = "PFOS_Lib01" if stem.startswith("PFOS") else "WT_Lib01"
            called = len(pd.read_csv(source / "07_source_data" / f"{library}_called_cell_qc.tsv", sep="\t"))
            draw_qc_source(frame, stem, config, called, target)
        elif stem.startswith("umap_"):
            column = "cell_type" if stem == "umap_cell_type" else "label"
            titles = {"umap_library": "Library", "umap_condition": "Condition (descriptive)",
                      "umap_leiden": "Leiden clusters", "umap_cell_type": "Marker-first cell type"}
            categorical_embedding(frame[["UMAP1", "UMAP2"]].to_numpy(), frame[column], titles[stem], target, dpi)
        elif stem == "captured_composition":
            stacked_composition(frame, target, dpi)
        # Preserve the exact, signed plot source rather than rewriting its serialization.
        table = output / entry["stage"] / "tables" / f"{stem}.tsv"
        table.parent.mkdir(parents=True, exist_ok=True)
        if table.is_file():
            pd.testing.assert_frame_equal(pd.read_csv(table, sep="\t"), frame,
                                          check_dtype=False, rtol=1e-5, atol=1e-6)
        shutil.copy2(entry["source_path"], table)


def _render_beta(source: Path, output: Path, entries: list[dict]) -> None:
    from . import extension_beta as beta
    thresholds = json.loads((source / "beta_identity_review_thresholds.json").read_text())
    for entry in entries:
        stem = Path(entry["png"]).stem
        frame = pd.read_csv(source / entry["source"], sep="\t", low_memory=False)
        if stem == "beta_umap_library":
            fig = beta._draw_umap_categorical(frame, "library_id", "Beta cells by library", {"PFOS_Lib01": beta.CONDITION_COLORS["PFOS"], "WT_Lib01": beta.CONDITION_COLORS["WT"]})
        elif stem == "beta_umap_clusters":
            fig = beta._draw_umap_categorical(frame, beta.PREFERRED_CLUSTER_KEY, "Beta-state clusters (Leiden 0.5)", beta.BETA_CLUSTER_COLORS)
        elif stem == "beta_cell_cycle_phase_umap":
            fig = beta._draw_umap_categorical(frame, "beta_cell_cycle_phase", "Score-derived cell-cycle phase", beta.CELL_CYCLE_COLORS)
        elif stem in {"beta_umap_state_scores", "beta_target_feature_umap", "beta_cell_cycle_umap", "beta_relative_state_distance"}:
            cols = [f"beta_score_{name}" for name in beta.STATE_GENE_SETS] if stem == "beta_umap_state_scores" else ["Mpc1", "Mpc2", "Atg7"] if stem == "beta_target_feature_umap" else ["beta_S_score", "beta_G2M_score"] if stem == "beta_cell_cycle_umap" else ["beta_relative_state_distance"]
            fig = beta._draw_continuous_grid(frame, cols, entry["title"], ncols=1 if stem == "beta_relative_state_distance" else 3)
        elif stem == "beta_identity_hormone_qc_audit":
            parts = [frame[frame["row_kind"].eq(kind)].dropna(axis=1, how="all") for kind in ("cell_qc", "gene_layer_summary", "cluster_review_summary")]
            fig = beta._draw_identity_audit(*parts, thresholds)
        elif stem == "beta_cell_cycle_composition":
            fig = beta._draw_composition(frame, "beta_cell_cycle_phase", "Captured score-derived cell-cycle phase composition", "Cell-cycle phase", beta.CELL_CYCLE_COLORS)
        elif stem == "beta_cluster_composition":
            fig = beta._draw_composition(frame, colors=beta.BETA_CLUSTER_COLORS)
        else:
            functions = {"beta_marker_dotplot": beta._draw_marker_dot, "beta_marker_heatmap": beta._draw_marker_heatmap,
                         "beta_target_violin": beta._draw_target_violin, "beta_qc_score_associations": beta._draw_qc_heatmap,
                         "beta_core_exclude_cluster4_sensitivity": beta._draw_cluster4_sensitivity,
                         "beta_paga_undirected": beta._draw_paga, "beta_trajectory_sensitivity": beta._draw_sensitivity,
                         "beta_balanced_subsample_sensitivity": beta._draw_balanced}
            fig = functions[stem](frame)
        save_figure(fig, output / "plots" / stem)


def _render_enrichment(source: Path, output: Path) -> None:
    from . import extension_enrichment as enrichment
    ora = pd.read_csv(source / "tables/ora_full.tsv", sep="\t", dtype={"group": str}, low_memory=False)
    status = pd.read_csv(source / "tables/ora_status.tsv", sep="\t", dtype={"group": str})
    resources = source.parents[2] / "reference/extension_gene_sets"
    for level in ("cluster", "cell_type", "beta_cluster"):
        groups = sorted(status.loc[status["analysis_level"].eq(level), "group"].astype(str).unique(), key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
        for group in groups:
            enrichment._plot_group_ora(ora, level, group, output)
        enrichment._overview_heatmap(ora, level, output)
        enrichment._term_network(ora, level, output, resources)
    enrichment._render_descriptive_figures(*[pd.read_csv(source / "tables" / name, sep="\t", low_memory=False) for name in (
        "descriptive_ranked_pathway_es.tsv", "descriptive_pathway_es_detection05_stability.tsv", "descriptive_pathway_es_equal_cell_stability.tsv")], output)


def _render_networks(source: Path, output: Path) -> None:
    from . import extension_networks as networks
    read = lambda name: pd.read_csv(source / "tables" / name, sep="\t", low_memory=False)
    plots, tables = output / "plots", output / "tables"
    networks._plot_communication_bubbles(read("communication_lr_full.tsv"), plots, tables)
    networks._plot_communication_network(read("communication_sender_receiver_summary.tsv"), plots, tables)
    networks._draw_communication_sensitivity(read("communication_sensitivity.tsv"), plots / "communication_sensitivity",
        library_conditions=read("communication_sender_receiver_summary.tsv")[["library_id", "condition"]].drop_duplicates().set_index("library_id")["condition"].astype(str).to_dict())
    edges = read("regulon_seed_stable_edges.tsv").merge(
        read("regulon_motif_supported_targets.tsv")[["TF", "target"]].drop_duplicates(), on=["TF", "target"], how="inner")
    networks._plot_regulons(read("regulon_library_celltype_auc_rss.tsv"), edges, plots, tables)


def replot(config: WorkflowConfig, source_run: Path, extension_run: Path, check: bool = False) -> dict:
    with PlotRun(config.root) as run:
        return _replot(config, source_run, extension_run, check, run)


def _replot(config: WorkflowConfig, source_run: Path, extension_run: Path, check: bool, run: PlotRun) -> dict:
    source = config.resolve(source_run)
    extension = config.resolve(extension_run)
    if json.loads((source / "00_audit/audit_status.json").read_text()).get("status") != "PASS":
        raise RuntimeError("基础 run 尚未审计通过")
    base = _base_entries(source)
    source_files = [item["source_path"] for item in base]
    source_files += [source / "04_annotation/annotated.h5ad"]
    source_files += [p for directory in ("00_audit", "02_qc", "04_annotation", "05_descriptive", "06_optional") for p in (source / directory).rglob("*.tsv") if p.name != "output_hashes.tsv"]
    source_files += [source / "01_cellranger/metrics_summary.tsv"]
    source_files += list((source / "02_qc").glob("*/*status.json"))
    before = _verify_producer(source, source / "00_audit/output_hashes.tsv", source_files)
    bundles = []
    for original, stage in (("04_beta_states", "05_beta_state"), ("05_enrichment", "06_enrichment"), ("06_networks", "07_networks")):
        directory = extension / original
        manifest = json.loads((directory / "report_manifest.json").read_text())
        if manifest.get("status") != "completed":
            raise RuntimeError(f"现役扩展模块未完成：{original}")
        required = [directory / entry["source"] for entry in manifest["figures"]]
        required += [p for p in (directory / "tables").glob("*") if p.is_file() and p.suffix in {".tsv", ".csv", ".gz", ".json"}]
        required += [directory / "report_manifest.json"]
        if original == "04_beta_states":
            required += [directory / "beta_identity_review_thresholds.json"]
        before.update(_verify_producer(directory, directory / "module_output_hashes.tsv", required))
        bundles.append((directory, stage, manifest["figures"]))
    resource_provenance = json.loads((extension / "05_enrichment/provenance.json").read_text())["resources"]
    for name, expected in resource_provenance["files"].items():
        path = config.root / "reference/extension_gene_sets" / name
        if not path.is_file() or sha256(path) != expected:
            raise RuntimeError(f"富集绘图资源与冻结签名不同：{name}")
        before[str(path)] = expected
    from .figure_style import PROFILE_RELATIVE
    profile_path = (Path(__file__).parent / PROFILE_RELATIVE).resolve(strict=True)
    source_config = [config.source, config.resolve(config.data["samples"]), config.resolve(config.data["markers"]), profile_path]
    before.update({str(path): sha256(path) for path in source_config})
    code_before = {p.name: sha256(p) for p in Path(__file__).parent.glob("*.py")}
    counts = {"base": len(base), **{stage: len(entries) for _, stage, entries in bundles}}
    if check:
        return {"status": "PASS", "mode": "plots-only", "figure_counts": counts, "sources_verified": len(before)}
    result = config.root / "result"
    work = run.stage(".result.plots-")
    output = work / "result"
    output.mkdir()
    receipt = {"status": "STARTED", "mode": "plots-only", "started_at": datetime.now(timezone.utc).isoformat(),
               "source_run": str(source.relative_to(config.root)), "extension_run": str(extension.relative_to(config.root)),
               "profile_path": str(__import__("os").path.relpath(profile_path, config.root)), "profile_sha256": sha256(profile_path),
               "source_signature": before, "code_signature": code_before, "statistics_recomputed": False,
               "scientific_status": {"base_audit": "PASS", "extension_terminal_audit": "MISSING", "velocity_count": "FAILED", "condition_inference": "NO-GO_n1_per_condition"}}
    try:
        for stage in ("00_audit", "01_qc", "02_atlas", "03_annotation", "04_descriptive", "05_beta_state", "06_enrichment", "07_networks"):
            (output / stage / "plots").mkdir(parents=True, exist_ok=True)
            (output / stage / "tables").mkdir(exist_ok=True)
        atomic_json(output / "00_audit/replot.json", receipt)
        _render_base(config, source, output, base)
        all_entries = []
        for item in base:
            all_entries.append({"figure_id": f"{source.name}/{item['stem']}", "title": item["title"], "interpretation": item["interpretation"],
                                "png": f"{item['stage']}/plots/{item['stem']}.png", "pdf": f"{item['stage']}/plots/{item['stem']}.pdf", "source": f"{item['stage']}/tables/{item['stem']}.tsv"})
        for directory, stage, entries in bundles:
            target = output / stage
            if stage == "05_beta_state":
                _render_beta(directory, target, entries)
            elif stage == "06_enrichment":
                _render_enrichment(directory, target)
            else:
                _render_networks(directory, target)
            for entry in entries:
                entry = dict(entry)
                if stage == "06_enrichment" and "_ora_dotplot" in entry["png"]:
                    entry["interpretation"] = entry["interpretation"].replace("气泡横轴为富集倍数、大小为命中marker数、颜色为−log10(nominal P)，边框颜色区分数据库。", "横轴为命中数/query_size；面积正比命中数；颜色为−log10(nominal P)，固定0–10，≥10封顶；条目标签标明数据库。")
                if stage == "06_enrichment" and "_ora_barplot" in entry["png"]:
                    entry["interpretation"] = entry["interpretation"].replace("条长为marker基因比例，颜色区分数据库。", "条长及颜色为−log10(nominal P)，色阶固定0–10、≥10封顶；条目标签标明数据库。")
                rel = Path(entry["png"]).relative_to("figures").with_suffix("")
                # Enrichment draws through its existing figures/source_data layout, then publishes the canonical layout.
                if stage == "06_enrichment":
                    sections = {"cluster": "01_cluster", "cell_type": "02_cell_type", "beta_cluster": "03_beta_cluster", "descriptive": "04_descriptive"}
                    section = Path(sections.get(str(rel.parent), str(rel.parent)))
                    dest = target / section
                    (dest / "plots").mkdir(parents=True, exist_ok=True)
                    (dest / "tables").mkdir(exist_ok=True)
                    for suffix in (".png", ".pdf"):
                        (target / "figures" / rel.with_suffix(suffix)).rename(dest / "plots" / rel.with_suffix(suffix).name)
                    source_target = dest / "tables" / Path(entry["source"]).name
                    media_base = Path(stage) / section / "plots" / rel.name
                else:
                    source_target = target / "tables" / Path(entry["source"]).name
                    media_base = Path(stage) / "plots" / rel.name
                shutil.copy2(directory / entry["source"], source_target)
                all_entries.append({**entry, "figure_id": f"{extension.name}/{directory.name}/{rel}",
                                    "png": str(media_base.with_suffix(".png")), "pdf": str(media_base.with_suffix(".pdf")),
                                    "source": str(source_target.relative_to(output))})
            if stage == "06_enrichment":
                shutil.rmtree(target / "figures")
                shutil.rmtree(target / "source_data")
            # Keep complete numerical result tables with exact bytes, once per owner stage.
            for table in (directory / "tables").glob("*"):
                if table.is_file() and table.suffix in {".tsv", ".csv", ".json", ".gz"}:
                    destination = target / "tables" / table.name
                    if destination.exists() and sha256(destination) != sha256(table):
                        destination = destination.with_name("complete_" + table.name)
                    shutil.copy2(table, destination)
        changed = [name for name, expected in before.items() if sha256(Path(name)) != expected]
        if changed:
            raise RuntimeError(f"原运行来源被改变，拒绝发布：{changed}")
        if {p.name: sha256(p) for p in Path(__file__).parent.glob("*.py")} != code_before:
            raise RuntimeError("绘图源码在执行期间发生变化，拒绝发布")
        actual = {str(path.relative_to(output)) for path in output.rglob("*.png")}
        expected = {entry["png"] for entry in all_entries}
        if actual != expected:
            raise RuntimeError(f"图件覆盖不一致：missing={expected-actual}, extra={actual-expected}")
        from PIL import Image
        for entry in all_entries:
            if not all((output / entry[key]).is_file() for key in ("png", "pdf", "source")):
                raise RuntimeError(f"图件未配对或来源缺失：{entry['figure_id']}")
            with Image.open(output / entry["png"]) as im:
                if im.mode != "RGB" or abs(im.width - 2161) > 1:
                    raise RuntimeError(f"图件规格不合格：{entry['figure_id']} {im.mode} {im.size}")
        for entry in all_entries:
            stem = Path(entry["png"]).stem
            if stem in BETA_DESCRIPTIONS:
                entry["title"], entry["interpretation"] = BETA_DESCRIPTIONS[stem]
        frame = pd.DataFrame(all_entries)
        frame.to_csv(output / "00_audit/tables/figure_manifest.tsv", sep="\t", index=False)
        _write_gallery(output, all_entries)
        receipt.update(status="COMPLETED", finished_at=datetime.now(timezone.utc).isoformat(),
                       figure_pairs=len(all_entries), figure_counts=counts, source_unchanged=True,
                       code_signature=code_before,
                       output_signature={str(p.relative_to(output)): sha256(p) for p in output.rglob("*") if p.is_file() and p != output / "00_audit/replot.json"})
        atomic_json(output / "00_audit/replot.json", receipt)
        run.publish(output, result)
        return receipt
    except BaseException as error:
        receipt.update(status="FAILED", error=str(error))
        raise


BETA_DESCRIPTIONS = {
    "beta_umap_library": ("β细胞按文库的 UMAP", "未做批次整合；文库与条件完全混杂。"),
    "beta_umap_clusters": ("β细胞状态聚类", "Leiden 0.5 为表达状态分区，仅用于描述与 marker 注释。"),
    "beta_umap_state_scores": ("β细胞连续程序评分", "程序评分使用表达匹配的对照基因，保持连续量，不强制指定最高分状态。"),
    "beta_marker_dotplot": ("β细胞状态 marker 气泡图", "点大小为检出比例，颜色为群间标准化平均表达。"),
    "beta_marker_heatmap": ("β细胞状态 marker 热图", "展示原定的合格优先 marker 和冻结程序锚点；完整 marker 表保留。"),
    "beta_target_feature_umap": ("焦点基因在 β细胞中的定位", "表达为 normalize-total log1p，只支持细胞定位。"),
    "beta_target_violin": ("焦点基因的文库分布", "无轮廓填色表示捕获细胞表达密度；点为中位数、细线为Q1–Q3。每条件一个文库，仅作描述，不作动物级推断。"),
    "beta_cell_cycle_umap": ("β细胞周期程序评分", "周期评分不证明细胞增殖或停滞。"),
    "beta_cell_cycle_phase_umap": ("β细胞评分定义的周期阶段", "G1/S/G2M 为评分定义，不是独立验证的增殖表型。"),
    "beta_cell_cycle_composition": ("捕获的评分定义周期组成", "各文库内比例没有条件检验 P 值。"),
    "beta_cluster_composition": ("捕获的 β状态组成", "比例描述捕获细胞，不能估计动物层面的组成效应。"),
    "beta_qc_score_associations": ("β程序与 QC 的相关性诊断", "Spearman 系数为文库内诊断，没有条件 P 值。"),
    "beta_identity_hormone_qc_audit": ("β身份与非同源激素 QC 审计", "保留原 β标签；复杂度、残留激素、校正敏感性及双细胞指标均为中性复核证据。"),
    "beta_core_exclude_cluster4_sensitivity": ("全部 β细胞与排除状态群 4 的敏感性", "全部原标注 β细胞为主要结果；排除群 4 仅为技术敏感性，不是修改 QC 规则。"),
    "beta_paga_undirected": ("无方向的 β状态图", "PAGA 边权表示群连接，不包含时间或因果方向。"),
    "beta_relative_state_distance": ("假设生成的相对状态距离", "方向性拟时未通过门槛；距离只描述冻结表达图上的相对几何位置。"),
    "beta_trajectory_sensitivity": ("β状态几何敏感性", "比较邻居、主成分、随机种子及起点扰动；不解释为成年胰岛发育时间。"),
    "beta_balanced_subsample_sensitivity": ("β状态等细胞数抽样敏感性", "等细胞数只评估技术稳定性，不能创造生物学重复。"),
}


def _write_gallery(output: Path, entries: list[dict]) -> None:
    lines = ["<!doctype html><html lang='zh'><meta charset='utf-8'><title>PFOS 胰岛单细胞完整图件</title>",
             "<style>body{font:16px sans-serif;max-width:1100px;margin:auto;padding:30px}img{max-width:100%}figure{margin:36px 0}</style>",
             "<h1>PFOS 胰岛单细胞完整图件</h1><p>每条件一个文库，文库与条件完全混杂。所有条件比较保持描述性；marker ORA 用于注释。velocity 计数失败，扩展整体终审仍为 MISSING。</p>"]
    for entry in entries:
        lines.append(f"<figure id='{html.escape(entry['figure_id'])}'><h2>{html.escape(entry['title'])}</h2><img src='{entry['png']}'><figcaption>生物学意义：描述胰岛细胞类型、状态或计算候选。数据分布特征：{html.escape(entry['interpretation'])} 最终结论：沿用原运行的证据边界，本次重绘不增加条件推断。 <a href='{entry['pdf']}'>PDF</a> · <a href='{entry['source']}'>源数据</a></figcaption></figure>")
    atomic_text(output / "index.html", "\n".join(lines) + "</html>\n")

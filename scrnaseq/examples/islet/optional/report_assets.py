from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np
import pandas as pd

from .plotting import DIVERGING_CMAP, CELL_TYPE_COLORS, save_figure
from .figure_style import COLORS
from .runtime import RunContext, atomic_text, atomic_tsv


ALWAYS_REPORT_FIGURES = (
    "marker_evidence_dotplot",
    "marker_evidence_heatmap",
    "annotation_evidence_map",
    "focal_gene_localization",
    "module_score_profile",
    "balanced_subsampling_stability",
    "evidence_gate_matrix",
)

COMPANY_REPORT_DECISIONS = (
    ("项目与分析流程", "项目说明、实验与生信流程", "adopt", "输入冻结、阶段状态、软件与输出哈希", "补充可审计的数据谱系"),
    ("QC 与细胞过滤", "过滤前后 QC 图", "adopt_expand", "barcode-rank、阈值、mt 敏感性、性别 marker 与保留率", "按文库报告并冻结阈值"),
    ("多细胞过滤", "doublet 分数与过滤", "adopt_expand", "按 chemistry/recovered-cell 表校准 Scrublet", "报告 expected rate、threshold 和预测数"),
    ("ambient RNA", "报告结构中未见专门审计", "current_only", "SoupX 自动估计及胰岛激素校正前后审计", "胰岛数据的关键污染控制"),
    ("聚类与降维", "UMAP 与 cluster", "adopt_adapt", "未整合表达的 PCA/邻接图/Leiden/UMAP", "sample 与 condition 完全混杂，拒绝批次整合"),
    ("cluster marker", "dotplot、heatmap、富集", "adopt_adapt", "marker dotplot/heatmap；|log2FC|>=0.58 且 P<0.05 仅用于注释", "不能把 cluster marker 检验挪作条件推断"),
    ("SingleR", "cluster/cell 自动注释", "adapt", "cluster 聚合、通用参考、弱辅助且不覆盖 marker", "避免泛参考自动覆盖胰岛 marker"),
    ("参考图谱注释", "未见独立胰岛 atlas 冲突门", "current_only", "MIA 小鼠胰岛 atlas 辅助、冲突 veto、低置信保留 unknown", "同时保存 marker/MIA/SingleR 证据"),
    ("细胞类型组成", "频数与比例", "adapt", "captured composition", "只描述两份文库，不称显著组成变化"),
    ("条件差异表达", "按细胞类型的 PFOS-vs-WT DEG/火山图", "reject_replace", "pseudobulk counts、描述性 log2 ratio、检测率差与等细胞数敏感性", "每组一库，禁止条件 P/FDR/DEG"),
    ("差异基因富集", "GO/Pathway 条件富集", "reject_replace", "预定义生物学模块连续评分", "无合格条件 DEG，不能生成条件富集显著性"),
    ("平均表达", "细胞群平均表达", "adopt_expand", "library×cell type pseudobulk 与 Mpc1/Mpc2/Atg7 定位", "保存整数 counts、均值与检测率"),
    ("CellChat", "通讯网络与通路比较", "NO-GO", "不运行条件通讯比较", "动物层面 n=1/组；细胞不能充当重复"),
    ("SCENIC", "条件/细胞群 regulon", "NO-GO", "不运行条件特异 regulon", "条件与文库完全混杂；仅可另行审计 pooled 候选"),
    ("细胞周期", "细胞周期预测", "not_in_scope", "当前未运行", "不是预先定义的成年胰岛终点；不得事后解释为 PFOS 效应"),
    ("Monocle3", "拟时与轨迹基因", "NO-GO", "不运行成年胰岛发育轨迹", "无验证的起点、方向或真实时间轴"),
    ("inferCNV", "拷贝数推断", "NO-GO", "不运行", "非肿瘤胰岛、无 CNV 假设与合格参考；不属于当前问题"),
    ("scVelo/PAGA", "velocity 与关键基因", "NO-GO", "不运行条件 velocity", "无独立审计的 spliced/unspliced 契约且 n=1/组"),
    ("复现与主张审计", "报告结构中未见逐图源数据和 claim map", "current_only", "同名 TSV/PNG/PDF、H5AD 契约、eligibility 与 claim-result map", "把可复现性和结论边界作为正式结果"),
)


def _natural(values: pd.Series | list[str] | np.ndarray) -> list[str]:
    unique = {str(value) for value in values}
    return sorted(unique, key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value))


def _read_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def report_figure_stems(context: RunContext) -> list[str]:
    stems = list(ALWAYS_REPORT_FIGURES)
    if context.path("01_cellranger", "metrics_summary.tsv").is_file() and context.path("02_qc", "qc_summary.tsv").is_file():
        stems.append("qc_retention_flow")
    if list(context.path("02_qc").glob("*/soupx_hormone_audit.tsv")):
        stems.append("soupx_hormone_audit")
    if list(context.path("02_qc").glob("*/mt_sensitivity.tsv")):
        stems.append("mt_threshold_sensitivity")
    return stems


def report_asset_outputs(context: RunContext) -> list[Path]:
    outputs: list[Path] = []
    for stem in report_figure_stems(context):
        outputs.extend([
            context.path("07_source_data", f"{stem}.tsv"),
            context.path("08_figures", f"{stem}.png"),
            context.path("08_figures", f"{stem}.pdf"),
        ])
    outputs.extend([
        context.path("00_audit", "company_report_structure_comparison.tsv"),
        context.path("09_reports", "company_report_structure_gap_analysis.md"),
        context.path("09_reports", "eligible_cluster_markers_top.tsv"),
    ])
    return outputs


def _company_structure_assets(context: RunContext) -> None:
    columns = ["company_module", "presentation_reviewed", "decision", "current_delivery", "rationale"]
    frame = pd.DataFrame(COMPANY_REPORT_DECISIONS, columns=columns)
    frame.insert(0, "company_data_used", False)
    frame.insert(0, "source_scope", "report.html headings and image references only")
    atomic_tsv(context.path("00_audit", "company_report_structure_comparison.tsv"), frame)
    lines = [
        "# 公司报告结构学习与差距处理",
        "",
        "- 审阅范围：仅 `report.html` 的标题、章节和图件引用；未读取或导入公司表达矩阵、结果表、RDS、文本数据或数值结果。",
        "- 公司数据用于本分析：`false`。",
        "- `.archive/vendor_2023` 权限：只读参考，正式结果全部来自当前项目冻结 FASTQ 和本流程产物。",
        "- 决策原则：展示形式可以采用；不满足生物学单位和设计资格的统计分析不复制。",
        "",
        "| 公司模块 | 决策 | 当前交付 | 理由 |",
        "|---|---|---|---|",
    ]
    for row in frame.itertuples(index=False):
        lines.append(f"| {row.company_module} | {row.decision} | {row.current_delivery} | {row.rationale} |")
    atomic_text(context.path("09_reports", "company_report_structure_gap_analysis.md"), "\n".join(lines) + "\n")


def _top_cluster_markers(context: RunContext) -> None:
    path = context.path("03_atlas", "cluster_markers.tsv")
    columns = ["cluster", "names", "logfoldchanges", "pvals", "pvals_adj", "pct_nz_group", "pct_nz_reference", "eligible_for_annotation", "test_scope"]
    if not path.is_file():
        atomic_tsv(context.path("09_reports", "eligible_cluster_markers_top.tsv"), pd.DataFrame(columns=columns))
        return
    markers = pd.read_csv(path, sep="\t")
    eligible = markers.loc[markers["eligible_for_annotation"].astype(bool)].copy()
    eligible = eligible.loc[~eligible["names"].astype(str).str.match(r"^(mt-|Rp[sl]|Gm\d+|Rik$)", case=False, na=False)]
    eligible = eligible.sort_values(["cluster", "scores"], ascending=[True, False]).groupby("cluster", observed=True, sort=False).head(5)
    atomic_tsv(context.path("09_reports", "eligible_cluster_markers_top.tsv"), eligible.reindex(columns=columns))


def _plot_marker_evidence(context: RunContext, data: ad.AnnData, dpi: int) -> None:
    preferred = float(context.config.data["preprocessing"]["preferred_resolution"])
    cluster_key = f"leiden_{preferred:g}"
    clusters = _natural(data.obs[cluster_key].astype(str))
    marker_map: dict[str, list[str]] = {}
    for level, sets in context.config.markers().items():
        for label, spec in sets.items():
            for gene in [item for item in spec.get("positive", []) if item in data.var_names][:3]:
                marker_map.setdefault(gene, []).append(f"{level}:{label}")
    genes = list(marker_map)
    if not genes:
        raise RuntimeError("报告 marker 图无可用配置基因")
    expression = data[:, genes].X
    counts = data.layers["counts"][:, [data.var_names.get_loc(gene) for gene in genes]]
    cluster_values = data.obs[cluster_key].astype(str).to_numpy()
    rows = []
    means = []
    detections = []
    for cluster in clusters:
        mask = cluster_values == cluster
        means.append(np.asarray(expression[mask].mean(axis=0)).ravel())
        detections.append(np.asarray((counts[mask] > 0).mean(axis=0)).ravel())
    mean_matrix = np.vstack(means)
    detection_matrix = np.vstack(detections)
    gene_mean = mean_matrix.mean(axis=0)
    gene_sd = mean_matrix.std(axis=0)
    z_matrix = np.divide(mean_matrix - gene_mean, gene_sd, out=np.zeros_like(mean_matrix), where=gene_sd > 0)
    for row_index, cluster in enumerate(clusters):
        for gene_index, gene in enumerate(genes):
            rows.append({
                "cluster": cluster,
                "gene": gene,
                "configured_positive_for": ";".join(marker_map[gene]),
                "mean_log1p": float(mean_matrix[row_index, gene_index]),
                "detection_rate": float(detection_matrix[row_index, gene_index]),
                "cluster_zscore": float(z_matrix[row_index, gene_index]),
                "scope": "configured marker visualization; cluster annotation only",
            })
    source = pd.DataFrame(rows)
    for stem in ("marker_evidence_dotplot", "marker_evidence_heatmap"):
        atomic_tsv(context.path("07_source_data", f"{stem}.tsv"), source)

    width = max(11.0, len(genes) * 0.34)
    fig = plt.figure(figsize=(width, max(5.2, len(clusters) * 0.42) + 0.65))
    layout = fig.add_gridspec(2, 2, width_ratios=[1, 0.045], height_ratios=[1, 0.14])
    ax = fig.add_subplot(layout[0, 0])
    colorbar_ax = fig.add_subplot(layout[0, 1])
    legend_ax = fig.add_subplot(layout[1, :])
    legend_ax.set_axis_off()
    x = np.tile(np.arange(len(genes)), len(clusters))
    y = np.repeat(np.arange(len(clusters)), len(genes))
    scatter = ax.scatter(x, y, s=120 * detection_matrix.ravel(), c=z_matrix.ravel(), cmap=DIVERGING_CMAP, vmin=-2.5, vmax=2.5, linewidths=0.15, edgecolors="#444444")
    ax.set_xticks(np.arange(len(genes)), genes, rotation=70, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(clusters)), clusters)
    ax.set_xlabel("Configured positive marker")
    ax.set_ylabel(f"Leiden {preferred:g} cluster")
    ax.set_title("Marker evidence: color = cluster z-score; area = detection rate")
    legend_ax.legend(handles=[ax.scatter([], [], s=120*f, color="#B7BEC6", label=f"{100*f:g}%") for f in (.25, .5, 1)], title="Detected cells (linear point area)", frameon=False, loc="center", ncol=3)
    ax.invert_yaxis()
    colorbar = fig.colorbar(scatter, cax=colorbar_ax, label="Within-gene cluster z-score", ticks=[-2.5, 0, 2.5])
    colorbar.ax.set_yticklabels(["≤−2.5", "0", "≥2.5"])
    save_figure(fig, context.path("08_figures", "marker_evidence_dotplot"), dpi)

    fig, ax = plt.subplots(figsize=(max(7.0, len(clusters) * 0.55), max(8.0, len(genes) * 0.28)))
    image = ax.imshow(z_matrix.T, cmap=DIVERGING_CMAP, vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(np.arange(len(clusters)), clusters)
    ax.set_yticks(np.arange(len(genes)), genes, fontsize=8)
    ax.set_xlabel(f"Leiden {preferred:g} cluster")
    ax.set_ylabel("Configured positive marker")
    ax.set_title("Marker evidence heatmap (within-gene standardized)")
    colorbar = fig.colorbar(image, ax=ax, label="Within-gene cluster z-score", shrink=0.7, ticks=[-2.5, 0, 2.5])
    colorbar.ax.set_yticklabels(["≤−2.5", "0", "≥2.5"])
    save_figure(fig, context.path("08_figures", "marker_evidence_heatmap"), dpi)


def _label_color(label: str) -> str:
    value = str(label).lower().replace(" ", "_")
    for cell_type, color in CELL_TYPE_COLORS.items():
        if cell_type in value:
            return color
    if any(token in value for token in ("macroph", "monocyte", "lymph", "t_cell", "b_cell")):
        return CELL_TYPE_COLORS["immune"]
    if any(token in value for token in ("fibro", "stellate", "pericyte")):
        return CELL_TYPE_COLORS["mesenchymal"]
    return "#E6E6E6" if value not in {"unknown", "not_available", "nan"} else CELL_TYPE_COLORS["unknown"]


def _plot_annotation_evidence(context: RunContext, dpi: int) -> None:
    decisions = pd.read_csv(context.path("04_annotation", "annotation_decisions.tsv"), sep="\t", dtype=str)
    singler_path = context.path("04_annotation", "singler_support.tsv")
    if singler_path.is_file():
        singler = pd.read_csv(singler_path, sep="\t", dtype=str)[["cluster", "singler_pruned_label"]]
        decisions = decisions.merge(singler, on="cluster", how="left")
    else:
        decisions["singler_pruned_label"] = "not_available"
    decisions["singler_pruned_label"] = decisions["singler_pruned_label"].fillna("not_available")
    order = {cluster: index for index, cluster in enumerate(_natural(decisions["cluster"]))}
    decisions = decisions.sort_values("cluster", key=lambda column: column.astype(str).map(order))
    source = decisions[["cluster", "marker_label", "mia_label", "singler_pruned_label", "final_label", "annotation_review_status", "decision_reason"]].copy()
    source["singler_role"] = "generic weak auxiliary; does not overwrite marker/MIA decision"
    atomic_tsv(context.path("07_source_data", "annotation_evidence_map.tsv"), source)
    columns = ["marker_label", "mia_label", "singler_pruned_label", "final_label"]
    labels = ["Marker hierarchy", "MIA islet atlas", "SingleR generic", "Final retained"]
    fig, ax = plt.subplots(figsize=(11.5, max(5.0, len(source) * 0.52)))
    ax.set_xlim(0, len(columns))
    ax.set_ylim(0, len(source))
    for y, row in enumerate(source.itertuples(index=False)):
        for x, column in enumerate(columns):
            value = str(getattr(row, column))
            color = _label_color(value)
            ax.add_patch(plt.Rectangle((x, y), 1, 1, facecolor=color, edgecolor="white", linewidth=1.5))
            shown = value if len(value) <= 18 else value[:16] + "…"
            red, green, blue = to_rgb(color)
            text_color = "white" if 0.2126 * red + 0.7152 * green + 0.0722 * blue < 0.48 else "black"
            ax.text(x + 0.5, y + 0.5, shown, ha="center", va="center", fontsize=7, color=text_color)
    ax.set_xticks(np.arange(len(columns)) + 0.5, labels)
    ax.set_yticks(np.arange(len(source)) + 0.5, source["cluster"].astype(str))
    ax.set_ylabel("Leiden cluster")
    ax.set_title("Annotation evidence audit (SingleR is non-overwriting weak support)")
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)
    save_figure(fig, context.path("08_figures", "annotation_evidence_map"), dpi)


def _plot_focal_genes(context: RunContext, dpi: int) -> None:
    frame = pd.read_csv(context.path("05_descriptive", "focal_gene_localization.tsv"), sep="\t")
    atomic_tsv(context.path("07_source_data", "focal_gene_localization.tsv"), frame)
    libraries = list(dict.fromkeys(frame["library_id"].astype(str)))
    cell_types = [item for item in CELL_TYPE_COLORS if item in set(frame["cell_type"].astype(str))]
    cell_types.extend(sorted(set(frame["cell_type"].astype(str)) - set(cell_types)))
    genes = list(dict.fromkeys(frame["gene"].astype(str)))
    mean_max = max(float(frame["mean_log1p"].max()), 1e-6)
    fig, axes = plt.subplots(1, len(libraries), figsize=(max(8.5, len(libraries) * 5.2), 4.3), squeeze=False, sharey=True)
    plotted = None
    for axis, library in zip(axes[0], libraries, strict=True):
        subset = frame[frame["library_id"].astype(str) == library]
        for y, gene in enumerate(genes):
            for x, cell_type in enumerate(cell_types):
                row = subset[(subset["gene"].astype(str) == gene) & (subset["cell_type"].astype(str) == cell_type)]
                if row.empty:
                    continue
                plotted = axis.scatter(x, y, s=200 * float(row.iloc[0]["detection_rate"]), c=[float(row.iloc[0]["mean_log1p"])], cmap="viridis", vmin=0, vmax=mean_max, edgecolors="#333333", linewidths=0.3)
        axis.set_xticks(np.arange(len(cell_types)), cell_types, rotation=40, ha="right")
        axis.set_yticks(np.arange(len(genes)), genes)
        axis.set_title(library)
        axis.set_xlabel("Cell type")
        axis.invert_yaxis()
    axes[0, 0].set_ylabel("Focal gene")
    if plotted is not None:
        fig.colorbar(plotted, ax=axes.ravel().tolist(), label="Mean normalize-total log1p expression", shrink=0.78)
    fig.suptitle("Mpc1/Mpc2/Atg7 localization: size = detection rate")
    save_figure(fig, context.path("08_figures", "focal_gene_localization"), dpi)


def _plot_modules(context: RunContext, dpi: int) -> None:
    wide = pd.read_csv(context.path("05_descriptive", "module_score_summary.tsv"), sep="\t")
    mean_columns = [column for column in wide if column.startswith("module_") and column.endswith("_mean")]
    rows = []
    for row in wide.itertuples(index=False):
        values = row._asdict()
        for column in mean_columns:
            rows.append({
                "library_id": values["library_id"],
                "condition": values["condition"],
                "cell_type": values["cell_type"],
                "module": column.removeprefix("module_").removesuffix("_mean"),
                "mean_score": float(values[column]),
            })
    frame = pd.DataFrame(rows)
    frame["visual_zscore"] = frame.groupby("module", observed=True)["mean_score"].transform(
        lambda values: (values - values.mean()) / values.std(ddof=0) if values.std(ddof=0) > 0 else 0.0
    )
    frame["scope"] = "module-wise z-score for visualization; raw mean retained; no condition inference"
    atomic_tsv(context.path("07_source_data", "module_score_profile.tsv"), frame)
    libraries = list(dict.fromkeys(frame["library_id"].astype(str)))
    modules = sorted(frame["module"].astype(str).unique())
    cell_types = [item for item in CELL_TYPE_COLORS if item in set(frame["cell_type"].astype(str))]
    cell_types.extend(sorted(set(frame["cell_type"].astype(str)) - set(cell_types)))
    fig, axes = plt.subplots(1, len(libraries), figsize=(max(10.0, len(libraries) * 5.6), max(5.3, len(modules) * 0.5)), squeeze=False, sharey=True)
    image = None
    for axis, library in zip(axes[0], libraries, strict=True):
        pivot = frame[frame["library_id"].astype(str) == library].pivot(index="module", columns="cell_type", values="visual_zscore").reindex(index=modules, columns=cell_types)
        image = axis.imshow(pivot.to_numpy(dtype=float), cmap=DIVERGING_CMAP, vmin=-2.5, vmax=2.5, aspect="auto")
        axis.set_xticks(np.arange(len(cell_types)), cell_types, rotation=40, ha="right")
        axis.set_yticks(np.arange(len(modules)), [item.replace("_", " ") for item in modules], fontsize=8)
        axis.set_title(library)
        axis.set_xlabel("Cell type")
    axes[0, 0].set_ylabel("Predefined module")
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), label="Within-module z-score", shrink=0.76)
    fig.suptitle("Continuous module-score profiles (descriptive)")
    save_figure(fig, context.path("08_figures", "module_score_profile"), dpi)


def _plot_balanced(context: RunContext, dpi: int) -> None:
    frame = pd.read_csv(context.path("05_descriptive", "balanced_subsampling_stability.tsv"), sep="\t")
    atomic_tsv(context.path("07_source_data", "balanced_subsampling_stability.tsv"), frame)
    targets = list(dict.fromkeys(frame["target"].astype(str)))
    cell_types = [item for item in CELL_TYPE_COLORS if item in set(frame["cell_type"].astype(str))]
    cell_types.extend(sorted(set(frame["cell_type"].astype(str)) - set(cell_types)))
    median = frame.pivot(index="target", columns="cell_type", values="median_case_minus_control").reindex(index=targets, columns=cell_types)
    stability = frame.pivot(index="target", columns="cell_type", values="direction_stability").reindex(index=targets, columns=cell_types)
    values = median.to_numpy(dtype=float)
    finite = np.abs(values[np.isfinite(values)])
    limit = float(np.quantile(finite, 0.95)) if finite.size else 1.0
    limit = max(limit, 1e-6)
    fig, axes = plt.subplots(1, 2, figsize=(max(12.0, len(cell_types) * 1.35), max(7.0, len(targets) * 0.45)), sharey=True)
    first = axes[0].imshow(values, cmap=DIVERGING_CMAP, vmin=-limit, vmax=limit, aspect="auto")
    second = axes[1].imshow(stability.to_numpy(dtype=float), cmap="viridis", vmin=0.5, vmax=1.0, aspect="auto")
    for axis in axes:
        axis.set_xticks(np.arange(len(cell_types)), cell_types, rotation=40, ha="right")
        axis.set_yticks(np.arange(len(targets)), [item.removeprefix("module_").replace("_", " ") for item in targets], fontsize=8)
        axis.set_xlabel("Cell type")
    axes[0].set_title("Median PFOS - WT after equal-cell sampling")
    axes[0].set_ylabel("Target")
    axes[1].set_title("Directional consistency across resamples")
    fig.colorbar(first, ax=axes[0], label="Median descriptive difference", shrink=0.75)
    fig.colorbar(second, ax=axes[1], label="Direction consistency", shrink=0.75)
    fig.suptitle("Technical subsampling sensitivity; not biological replication")
    save_figure(fig, context.path("08_figures", "balanced_subsampling_stability"), dpi)


def _plot_evidence_gate(context: RunContext, dpi: int) -> None:
    eligibility = pd.read_csv(context.path("00_audit", "analysis_eligibility.tsv"), sep="\t", dtype=str)
    optional = pd.read_csv(context.path("06_optional", "optional_module_status.tsv"), sep="\t", dtype=str)
    first = eligibility.rename(columns={"analysis": "analysis_item", "reason": "rationale"})[["analysis_item", "status", "statistical_unit", "rationale"]]
    second = optional.rename(columns={"module": "analysis_item", "reason": "rationale", "scope": "statistical_unit"})[["analysis_item", "status", "statistical_unit", "rationale"]]
    frame = pd.concat([first, second], ignore_index=True).drop_duplicates("analysis_item", keep="first")
    frame["status_class"] = frame["status"].map(lambda value: "GO" if value == "GO" else ("NO-GO" if value == "NO-GO" else "SUPPORT"))
    frame["scope"] = "eligibility decision; not a biological result"
    atomic_tsv(context.path("07_source_data", "evidence_gate_matrix.tsv"), frame)
    colors = {"GO": "#009E73", "NO-GO": "#D55E00", "SUPPORT": "#56B4E9"}
    fig, ax = plt.subplots(figsize=(12.5, max(6.0, len(frame) * 0.43)))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(frame))
    for y, row in enumerate(frame.itertuples(index=False)):
        ax.add_patch(plt.Rectangle((0, y), 1, 1, facecolor=colors[row.status_class], alpha=0.85, edgecolor="white"))
        ax.text(0.02, y + 0.5, str(row.analysis_item), ha="left", va="center", fontsize=8)
        ax.text(0.98, y + 0.5, str(row.status), ha="right", va="center", fontsize=8, fontweight="bold")
    ax.set_yticks([])
    ax.set_xticks([])
    ax.invert_yaxis()
    ax.set_title("Analysis eligibility: biological unit and design gate")
    for spine in ax.spines.values():
        spine.set_visible(False)
    save_figure(fig, context.path("08_figures", "evidence_gate_matrix"), dpi)


def _plot_qc_retention(context: RunContext, dpi: int) -> None:
    metrics = pd.read_csv(context.path("01_cellranger", "metrics_summary.tsv"), sep="\t", dtype=str)
    qc = pd.read_csv(context.path("02_qc", "qc_summary.tsv"), sep="\t")
    rows = []
    for _, item in metrics.iterrows():
        library = str(item["library_id"])
        called = int(str(item["Estimated Number of Cells"]).replace(",", ""))
        final = int(qc.set_index("library_id").loc[library, "cells"])
        scrublet = _read_json(context.path("02_qc", library, "scrublet_status.json"))
        rows.append({
            "library_id": library,
            "cellranger_called_cells": called,
            "scrublet_predicted_doublets": int(scrublet.get("predicted", 0)),
            "final_qc_singlets": final,
            "final_retention_fraction": final / called,
            "removed_after_cell_calling": called - final,
            "scope": "technical retention; not condition composition inference",
        })
    frame = pd.DataFrame(rows)
    atomic_tsv(context.path("07_source_data", "qc_retention_flow.tsv"), frame)
    x = np.arange(len(frame))
    width = 0.34
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    first = ax.bar(x - width / 2, frame["cellranger_called_cells"], width, color="#56B4E9", label="Cell Ranger called")
    second = ax.bar(x + width / 2, frame["final_qc_singlets"], width, color="#009E73", label="Final QC singlets")
    ax.bar_label(first, fmt="%d", padding=3, fontsize=8)
    ax.bar_label(second, fmt="%d", padding=3, fontsize=8)
    ax.set_xticks(x, frame["library_id"])
    ax.set_ylabel("Cells")
    ax.set_title("Cell-calling to final-QC retention")
    ax.legend(frameon=False)
    for index, row in frame.iterrows():
        ax.text(index, max(row["cellranger_called_cells"], row["final_qc_singlets"]) * 0.88, f"retained {row['final_retention_fraction']:.1%}\nScrublet {row['scrublet_predicted_doublets']:,}", ha="center", va="top", fontsize=8)
    save_figure(fig, context.path("08_figures", "qc_retention_flow"), dpi)


def _plot_soupx(context: RunContext, dpi: int) -> None:
    frames = [pd.read_csv(path, sep="\t") for path in sorted(context.path("02_qc").glob("*/soupx_hormone_audit.tsv"))]
    frame = pd.concat(frames, ignore_index=True)
    status = {path.parent.name: _read_json(path) for path in context.path("02_qc").glob("*/soupx_status.json")}
    frame["soupx_status"] = frame["library_id"].map(lambda value: status.get(str(value), {}).get("status", "unknown"))
    frame["rho_median"] = frame["library_id"].map(lambda value: status.get(str(value), {}).get("rho_median", np.nan))
    frame["fallback_used"] = frame["library_id"].map(lambda value: status.get(str(value), {}).get("fallback", False))
    atomic_tsv(context.path("07_source_data", "soupx_hormone_audit.tsv"), frame)
    genes = list(dict.fromkeys(frame["gene"].astype(str)))
    libraries = list(dict.fromkeys(frame["library_id"].astype(str)))
    x = np.arange(len(genes))
    width = 0.8 / len(libraries)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
    for index, library in enumerate(libraries):
        subset = frame[frame["library_id"].astype(str) == library].set_index("gene").reindex(genes)
        offset = (index - (len(libraries) - 1) / 2) * width
        color = COLORS[context.config.samples().set_index("library_id").loc[library, "condition"]]
        axes[0].bar(x + offset, subset["removed_fraction"], width, label=library, color=color)
        axes[1].plot(x, subset["raw_detection"], marker="o", label=f"{library} raw", color=color)
        axes[1].plot(x, subset["corrected_detection"], marker="s", linestyle="--", label=f"{library} corrected", color=color)
    axes[0].set_xticks(x, genes)
    axes[0].set_ylabel("Removed fraction of UMI counts")
    axes[0].set_title("SoupX hormone-count audit")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set_xticks(x, genes)
    axes[1].set_ylim(0, 1.03)
    axes[1].set_ylabel("Detection fraction across called cells")
    axes[1].set_title("Detection before and after correction")
    axes[1].legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle("Ambient RNA correction audit; auto-estimated rho, no fixed fallback")
    save_figure(fig, context.path("08_figures", "soupx_hormone_audit"), dpi)


def _plot_mt_sensitivity(context: RunContext, dpi: int) -> None:
    frames = [pd.read_csv(path, sep="\t") for path in sorted(context.path("02_qc").glob("*/mt_sensitivity.tsv"))]
    frame = pd.concat(frames, ignore_index=True)
    frame["scope"] = "mitochondrial threshold sensitivity only; all other QC criteria held outside this count"
    atomic_tsv(context.path("07_source_data", "mt_threshold_sensitivity.tsv"), frame)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
    for library, subset in frame.groupby("library_id", sort=False):
        subset = subset.sort_values("max_pct_mt")
        color = COLORS[context.config.samples().set_index("library_id").loc[library, "condition"]]
        axes[0].plot(subset["max_pct_mt"], subset["cells_passing_mt_only"], marker="o", label=library, color=color)
        axes[1].plot(subset["max_pct_mt"], subset["fraction_passing_mt_only"], marker="o", label=library, color=color)
    axes[0].set_ylabel("Called cells passing mt-only gate")
    axes[1].set_ylabel("Fraction passing mt-only gate")
    axes[1].set_ylim(0, 1)
    for axis in axes:
        axis.set_xlabel("Maximum mitochondrial counts (%)")
        axis.set_xticks(sorted(frame["max_pct_mt"].unique()))
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle("Mitochondrial threshold sensitivity (15/20/25%)")
    save_figure(fig, context.path("08_figures", "mt_threshold_sensitivity"), dpi)


def generate_report_assets(context: RunContext) -> list[str]:
    dpi = int(context.config.data["report"]["dpi"])
    stems = report_figure_stems(context)
    _company_structure_assets(context)
    _top_cluster_markers(context)
    data = ad.read_h5ad(context.path("04_annotation", "annotated.h5ad"))
    _plot_marker_evidence(context, data, dpi)
    _plot_annotation_evidence(context, dpi)
    _plot_focal_genes(context, dpi)
    _plot_modules(context, dpi)
    _plot_balanced(context, dpi)
    _plot_evidence_gate(context, dpi)
    if "qc_retention_flow" in stems:
        _plot_qc_retention(context, dpi)
    if "soupx_hormone_audit" in stems:
        _plot_soupx(context, dpi)
    if "mt_threshold_sensitivity" in stems:
        _plot_mt_sensitivity(context, dpi)
    return stems

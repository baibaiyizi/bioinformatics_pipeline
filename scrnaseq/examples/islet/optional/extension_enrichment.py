from __future__ import annotations

# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) marker ORA dot/bar → cross-type inherit → param inherit
# (b) ORA overview heatmap → heatmap visual parameters → param inherit
# (c) term-overlap network → cross-type inherit → param inherit
# (d) descriptive ES lollipop/sensitivity heatmap → cross-type inherit → param inherit
# RULE: No production asset has the same enrichment statistics and data structure.

import hashlib
import json
import math
import os
import re
import subprocess
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc

from .runtime import atomic_json, atomic_text, atomic_tsv, sha256


# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED = "#B2182B"
GREY = "#999999"
BLACK = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})


def save_cns_figure(fig: plt.Figure, filename: str) -> None:
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)


DATABASE_COLORS = {
    "GO_BP": "#0072B2",
    "GO_CC": "#009E73",
    "GO_MF": "#CC79A7",
    "KEGG": "#D55E00",
}
DATABASE_ORDER = ["GO_BP", "GO_CC", "GO_MF", "KEGG"]
LEVEL_CN = {"cluster": "全图谱 cluster", "cell_type": "细胞类型", "beta_cluster": "β细胞子群"}
MARKER_LOG2FC = 0.58
MARKER_PVALUE = 0.05
MIN_GENE_SET_SIZE = 10
MAX_GENE_SET_SIZE = 500
FIGURE_DPI = 300
SEED = 42


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _slug(value: object) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_.")
    return cleaned or "unnamed"


def _write_log(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"{_utc_now()}\t{message}\n")


def _save_figure(fig: plt.Figure, figure_stem: Path) -> tuple[Path, Path]:
    figure_stem.parent.mkdir(parents=True, exist_ok=True)
    png = figure_stem.with_suffix(".png")
    pdf = figure_stem.with_suffix(".pdf")
    save_cns_figure(fig, str(figure_stem))
    plt.close(fig)
    return png, pdf


def _marker_contract(frame: pd.DataFrame, analysis_level: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required = {
        "group",
        "names",
        "logfoldchanges",
        "pvals",
        "pvals_adj",
        "pct_nz_group",
        "pct_nz_reference",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"marker table lacks required columns: {missing}")
    work = frame.copy()
    work["group"] = work["group"].astype(str)
    work["names"] = work["names"].astype(str)
    for column in ("logfoldchanges", "pvals", "pvals_adj", "pct_nz_group", "pct_nz_reference"):
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work["tested_and_detectable"] = (
        np.isfinite(work["pvals"])
        & np.isfinite(work["logfoldchanges"])
        & ((work["pct_nz_group"] > 0) | (work["pct_nz_reference"] > 0))
    )
    work["selected_positive_marker"] = (
        work["tested_and_detectable"]
        & (work["logfoldchanges"] >= MARKER_LOG2FC)
        & (work["pvals"] < MARKER_PVALUE)
    )
    work["analysis_level"] = analysis_level
    work["selection_scope"] = (
        "positive marker: log2FC>=0.58 and nominal P<0.05; annotation only; "
        "not a PFOS-WT condition test"
    )
    universe = (
        work.loc[work["tested_and_detectable"], ["analysis_level", "group", "names"]]
        .drop_duplicates()
        .rename(columns={"names": "gene_symbol"})
        .sort_values(["analysis_level", "group", "gene_symbol"])
    )
    query = (
        work.loc[work["selected_positive_marker"], ["analysis_level", "group", "names"]]
        .drop_duplicates()
        .rename(columns={"names": "gene_symbol"})
        .sort_values(["analysis_level", "group", "gene_symbol"])
    )
    summary = (
        work.groupby(["analysis_level", "group"], observed=True)
        .agg(
            genes_tested=("names", "nunique"),
            genes_detectable=("tested_and_detectable", "sum"),
            positive_markers=("selected_positive_marker", "sum"),
        )
        .reset_index()
    )
    summary["marker_log2fc_threshold"] = MARKER_LOG2FC
    summary["marker_nominal_pvalue_threshold"] = MARKER_PVALUE
    summary["p_adjusted_used_for_marker_selection"] = False
    summary["scope"] = "marker annotation only; cells are not biological replicates"
    return work, query, universe, summary


def _cluster_marker_table(source_run: Path) -> pd.DataFrame:
    path = source_run / "03_atlas" / "cluster_markers.tsv"
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    if "cluster" not in frame.columns:
        raise ValueError("cluster marker table lacks cluster")
    return frame.rename(columns={"cluster": "group"})


def _completed_beta_cluster_markers(extended_run: Path) -> tuple[pd.DataFrame | None, list[Path]]:
    module_root = extended_run / "04_beta_states"
    marker_path = module_root / "tables" / "beta_cluster_markers_all.tsv"
    if not marker_path.exists():
        return None, []
    manifest_path = module_root / "report_manifest.json"
    execution_path = module_root / "execution_inputs.json"
    if not manifest_path.is_file() or not execution_path.is_file():
        raise RuntimeError("beta marker exists without completed beta module manifests")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    execution = json.loads(execution_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "completed" or execution.get("verification_status") != "PASS":
        raise RuntimeError("beta module is not in a completed, input-verified state")
    frame = pd.read_csv(marker_path, sep="\t")
    if "cluster" not in frame:
        raise ValueError("beta cluster marker table lacks cluster")
    return frame.rename(columns={"cluster": "group"}), [marker_path, manifest_path, execution_path]


def _compute_cell_type_markers(input_h5ad: Path, output_path: Path, log_path: Path) -> pd.DataFrame:
    if output_path.is_file():
        cached = pd.read_csv(output_path, sep="\t")
        required = {"group", "names", "logfoldchanges", "pvals", "pvals_adj", "pct_nz_group", "pct_nz_reference"}
        if required <= set(cached.columns):
            _write_log(log_path, f"reused cell-type marker table: {output_path}")
            return cached
    _write_log(log_path, f"loading H5AD for cell-type marker computation: {input_h5ad}")
    data = ad.read_h5ad(input_h5ad)
    if "cell_type" not in data.obs:
        raise ValueError("H5AD lacks obs['cell_type']")
    groups = sorted(set(data.obs["cell_type"].astype(str)) - {"unknown"})
    if len(groups) < 2:
        raise ValueError("fewer than two known cell types")
    cell_counts = data.obs["cell_type"].astype(str).value_counts()
    groups = [group for group in groups if int(cell_counts[group]) >= 20]
    if len(groups) < 2:
        raise ValueError("fewer than two known cell types with at least 20 cells")
    _write_log(log_path, f"rank_genes_groups start: groups={','.join(groups)}; method=wilcoxon; reference=rest")
    sc.tl.rank_genes_groups(
        data,
        groupby="cell_type",
        groups=groups,
        reference="rest",
        method="wilcoxon",
        use_raw=True,
        pts=True,
        n_genes=data.raw.n_vars if data.raw is not None else data.n_vars,
        key_added="extension_cell_type_markers",
    )
    tables: list[pd.DataFrame] = []
    for group in groups:
        frame = sc.get.rank_genes_groups_df(data, group=group, key="extension_cell_type_markers")
        frame.insert(0, "group", group)
        frame["cells_in_group"] = int(cell_counts[group])
        frame["cells_in_reference"] = int(data.n_obs - cell_counts[group])
        frame["test_scope"] = "cell_type marker annotation only; unknown cells remain in rest reference"
        tables.append(frame)
    result = pd.concat(tables, ignore_index=True)
    atomic_tsv(output_path, result)
    _write_log(log_path, f"rank_genes_groups complete: {result.shape[0]} rows")
    del data
    return result


def _resource_files(resources_dir: Path) -> list[Path]:
    return [
        resources_dir / "symbol_to_entrez.tsv",
        resources_dir / "go_term2gene.tsv",
        resources_dir / "go_term2name.tsv",
        resources_dir / "kegg_term2gene.tsv",
        resources_dir / "kegg_term2name.tsv",
        resources_dir / "resource_metadata.tsv",
    ]


def _go_names_are_complete(path: Path) -> bool:
    try:
        names = pd.read_csv(path, sep="\t", dtype=str)
    except (OSError, pd.errors.ParserError):
        return False
    required = {"term_id", "term_name", "ontology"}
    if not required <= set(names.columns) or names.empty:
        return False
    return bool(
        names["term_id"].notna().all()
        and names["term_name"].notna().all()
        and names["term_name"].str.len().gt(0).all()
        and names["ontology"].isin({"BP", "CC", "MF"}).all()
    )


def _run_r(command: list[str], log_path: Path, cwd: Path) -> None:
    _write_log(log_path, "$ " + " ".join(command))
    env = os.environ.copy()
    env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    with log_path.open("a", encoding="utf-8") as handle:
        subprocess.run(command, cwd=cwd, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT, text=True)


def _freeze_execution_inputs(paths: Iterable[Path], output_path: Path) -> dict[str, Any]:
    resolved = sorted({Path(path).resolve() for path in paths})
    missing = [str(path) for path in resolved if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"execution inputs missing: {missing}")
    snapshot = {
        "files": [{"path": str(path), "sha256": sha256(path)} for path in resolved],
        "started_at": _utc_now(),
        "verification_status": "pending",
    }
    atomic_json(output_path, snapshot)
    return snapshot


def _verify_execution_inputs(snapshot: dict[str, Any], output_path: Path) -> None:
    changed = []
    for item in snapshot["files"]:
        path = Path(item["path"])
        observed = sha256(path) if path.is_file() else "MISSING"
        if observed != item["sha256"]:
            changed.append({"path": str(path), "expected": item["sha256"], "observed": observed})
    snapshot["verified_at"] = _utc_now()
    snapshot["verification_status"] = "PASS" if not changed else "FAIL"
    snapshot["changed"] = changed
    atomic_json(output_path, snapshot)
    if changed:
        raise RuntimeError(f"execution inputs changed during enrichment: {changed}")


def _ensure_resources(
    resources_dir: Path,
    gene_symbols: Iterable[str],
    r_script: Path,
    log_path: Path,
    project_root: Path,
) -> dict[str, Any]:
    resources_dir.mkdir(parents=True, exist_ok=True)
    genes = pd.DataFrame({"gene_symbol": sorted(set(map(str, gene_symbols)))})
    genes_path = resources_dir / "resource_gene_symbols.tsv"
    atomic_tsv(genes_path, genes)
    genes_hash = sha256(genes_path)
    manifest_path = resources_dir / "resource_manifest.json"
    current_files = _resource_files(resources_dir)
    can_reuse = False
    if manifest_path.is_file() and all(path.is_file() for path in current_files):
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        can_reuse = old.get("source_gene_symbols_sha256") == genes_hash and all(
            old.get("files", {}).get(path.name) == sha256(path) for path in current_files
        ) and _go_names_are_complete(resources_dir / "go_term2name.tsv")
    if not can_reuse:
        _run_r(
            [
                "Rscript",
                str(r_script),
                "build-resources",
                "--resources-dir",
                str(resources_dir),
                "--genes",
                str(genes_path),
            ],
            log_path,
            project_root,
        )
    file_hashes = {path.name: sha256(path) for path in current_files}
    metadata = pd.read_csv(resources_dir / "resource_metadata.tsv", sep="\t")
    manifest = {
        "status": "frozen",
        "created_or_verified_at_utc": _utc_now(),
        "source_gene_symbols_sha256": genes_hash,
        "source_gene_count": int(len(genes)),
        "resources_dir": str(resources_dir.resolve()),
        "files": file_hashes,
        "metadata": dict(zip(metadata["key"], metadata["value"], strict=True)),
        "resource_scope": "Mus musculus genes detectable in the formal PFOS islet marker comparisons",
        "sources": {
            "GO": "org.Mm.eg.db GOALL/ONTOLOGYALL plus GO.db names",
            "KEGG": "KEGGREST live mouse pathway link/list endpoints, frozen locally",
        },
    }
    atomic_json(manifest_path, manifest)
    _write_log(log_path, f"resource status={'reused' if can_reuse else 'rebuilt'}; gene_hash={genes_hash}")
    return manifest


def _ora_non_name_fingerprint(path: Path) -> str:
    """Hash the serialized ORA table after removing only the term_name field."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        header = handle.readline().rstrip(b"\r\n").split(b"\t")
        if b"term_name" not in header:
            raise ValueError("ORA table lacks term_name")
        name_index = header.index(b"term_name")
        keep = [index for index in range(len(header)) if index != name_index]
        digest.update(b"\t".join(header[index] for index in keep) + b"\n")
        for line in handle:
            fields = line.rstrip(b"\r\n").split(b"\t")
            if len(fields) != len(header):
                raise ValueError("malformed ORA TSV while computing repair fingerprint")
            digest.update(b"\t".join(fields[index] for index in keep) + b"\n")
    return digest.hexdigest()


def _go_name_repair_baseline(output_dir: Path, resources_dir: Path) -> dict[str, Any] | None:
    ora_path = output_dir / "tables" / "ora_full.tsv"
    names_path = resources_dir / "go_term2name.tsv"
    required = [
        output_dir / "tables" / "marker_query.tsv",
        output_dir / "tables" / "marker_detectable_universe.tsv",
        resources_dir / "go_term2gene.tsv",
    ]
    if not ora_path.is_file() or not names_path.is_file() or not all(path.is_file() for path in required):
        return None
    ora_names = pd.read_csv(ora_path, sep="\t", usecols=["database", "term_name"])
    missing_go_names = int(
        (ora_names["database"].str.startswith("GO_") & ora_names["term_name"].isna()).sum()
    )
    if missing_go_names == 0:
        return None
    return {
        "detected_issue": "GO term names missing from a prior completed enrichment output",
        "prior_missing_go_term_name_rows": missing_go_names,
        "prior_non_name_ora_sha256": _ora_non_name_fingerprint(ora_path),
        "prior_input_sha256": {str(path.resolve()): sha256(path) for path in required},
    }


def _audit_go_name_repair(
    baseline: dict[str, Any] | None,
    output_dir: Path,
    resources_dir: Path,
) -> None:
    if baseline is None:
        return
    ora_path = output_dir / "tables" / "ora_full.tsv"
    required = [
        output_dir / "tables" / "marker_query.tsv",
        output_dir / "tables" / "marker_detectable_universe.tsv",
        resources_dir / "go_term2gene.tsv",
    ]
    ora_names = pd.read_csv(ora_path, sep="\t", usecols=["database", "term_name"])
    final_missing = int(
        (ora_names["database"].str.startswith("GO_") & ora_names["term_name"].isna()).sum()
    )
    final_non_name = _ora_non_name_fingerprint(ora_path)
    final_inputs = {str(path.resolve()): sha256(path) for path in required}
    audit = {
        **baseline,
        "final_missing_go_term_name_rows": final_missing,
        "final_non_name_ora_sha256": final_non_name,
        "final_input_sha256": final_inputs,
        "non_name_ora_unchanged": final_non_name == baseline["prior_non_name_ora_sha256"],
        "query_universe_term2gene_unchanged": final_inputs == baseline["prior_input_sha256"],
    }
    audit["status"] = "PASS" if (
        final_missing == 0
        and audit["non_name_ora_unchanged"]
        and audit["query_universe_term2gene_unchanged"]
    ) else "FAIL"
    atomic_json(output_dir / "tables" / "go_name_repair_audit.json", audit)
    if audit["status"] != "PASS":
        raise RuntimeError(f"GO name repair changed statistical content: {audit}")


def _term_label(row: pd.Series, width: int = 54) -> str:
    name = str(row["term_name"])
    if len(name) > width:
        name = name[: width - 1] + "…"
    return f"{row['database']}: {name}"


def _top_theme_text(
    frame: pd.DataFrame, count: int = 3, p_column: str = "pvalue", ascending: bool = True
) -> str:
    if frame.empty:
        return "无"
    ordered = frame.sort_values([p_column, "term_name"], ascending=[ascending, True])
    return "、".join(ordered["term_name"].astype(str).drop_duplicates().head(count))


def _display_terms(group_result: pd.DataFrame, per_database: int = 6) -> pd.DataFrame:
    eligible = group_result.loc[
        (group_result["pvalue"] < 0.05)
        & (group_result["overlap_count"] >= 2)
        & (group_result["fold_enrichment"] > 0)
    ].copy()
    if eligible.empty:
        return eligible
    eligible["database_rank"] = eligible.groupby("database", observed=True)["pvalue"].rank(
        method="first"
    )
    selected = eligible.loc[eligible["database_rank"] <= per_database].copy()
    selected["minus_log10_nominal_p"] = -np.log10(selected["pvalue"].clip(lower=np.finfo(float).tiny))
    selected["gene_ratio"] = selected["overlap_count"] / selected["query_size"]
    selected["display_label"] = selected.apply(_term_label, axis=1)
    selected["display_label"] = np.where(
        selected["passes_bh_0_05"].astype(bool), "[BH] " + selected["display_label"], selected["display_label"]
    )
    return selected.sort_values(["database", "pvalue", "p_adjust_bh", "term_id"])


def _empty_figure(title: str, message: str, stem: Path) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    ax.axis("off")
    ax.text(0.5, 0.58, title, ha="center", va="center", fontsize=10, fontweight="bold")
    ax.text(0.5, 0.38, message, ha="center", va="center", fontsize=8, color="#555555")
    return _save_figure(fig, stem)


def _plot_group_ora(
    result: pd.DataFrame,
    analysis_level: str,
    group: str,
    output_dir: Path,
) -> list[dict[str, str]]:
    group_result = result.loc[
        (result["analysis_level"] == analysis_level) & (result["group"].astype(str) == str(group))
    ].copy()
    shown = _display_terms(group_result)
    figure_dir = output_dir / "figures" / analysis_level
    source_dir = output_dir / "source_data" / analysis_level
    slug = _slug(group)
    entries: list[dict[str, str]] = []
    for kind in ("dotplot", "barplot"):
        stem_name = f"{analysis_level}_{slug}_ora_{kind}"
        figure_stem = figure_dir / stem_name
        source_path = source_dir / f"{stem_name}.tsv"
        atomic_tsv(source_path, shown)
        if shown.empty:
            png, pdf = _empty_figure(
                f"{analysis_level}: {group}",
                "No terms passed nominal P < 0.05 with at least two overlapping marker genes.",
                figure_stem,
            )
        else:
            plot = shown.sort_values(["database", "pvalue"], ascending=[True, False]).reset_index(drop=True)
            height = min(9.0, max(3.8, 1.6 + 0.30 * len(plot)))
            fig, ax = plt.subplots(figsize=(7.2, height))
            y = np.arange(len(plot))
            colors = [DATABASE_COLORS.get(db, "#777777") for db in plot["database"]]
            if kind == "dotplot":
                sizes = 22 + 16 * np.sqrt(plot["overlap_count"].to_numpy(dtype=float))
                scatter = ax.scatter(
                    plot["fold_enrichment"],
                    y,
                    s=sizes,
                    c=plot["minus_log10_nominal_p"],
                    cmap="viridis",
                    edgecolors=colors,
                    linewidths=1.0,
                )
                ax.set_xlabel("Fold enrichment")
                colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
                colorbar.set_label("−log10(nominal P)", fontsize=8)
            else:
                ax.barh(y, plot["gene_ratio"], color=colors, alpha=0.88)
                ax.set_xlabel("Marker gene ratio")
            ax.set_yticks(y, plot["display_label"])
            ax.set_ylabel("")
            ax.set_title(f"Exploratory marker ORA — {analysis_level}: {group}\nnominal display; annotation only")
            ax.grid(axis="x", color="#E0E0E0", linewidth=0.3, alpha=0.35)
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(axis="both", labelsize=7)
            fig.tight_layout()
            png, pdf = _save_figure(fig, figure_stem)
        entries.append(
            {
                "title": f"{LEVEL_CN.get(analysis_level, analysis_level)} {group} marker ORA {'气泡图' if kind == 'dotplot' else '条形图'}",
                "png": str(png.relative_to(output_dir)),
                "pdf": str(pdf.relative_to(output_dir)),
                "source": str(source_path.relative_to(output_dir)),
                "interpretation": (
                    f"实际展示 {len(shown)} 个 nominal P<0.05 且命中≥2基因的条目，其中 "
                    f"{int(shown['passes_bh_0_05'].astype(bool).sum()) if not shown.empty else 0} 个同时通过族内BH；"
                    f"最突出条目为 {_top_theme_text(shown)}。"
                    + (
                        "气泡横轴为富集倍数、大小为命中marker数、颜色为−log10(nominal P)，边框颜色区分数据库。"
                        if kind == "dotplot"
                        else "条长为marker基因比例，颜色区分数据库。"
                    )
                    + "这些P值仅用于细胞群功能注释，不是PFOS–WT条件推断。"
                ),
            }
        )
    return entries


def _overview_heatmap(result: pd.DataFrame, analysis_level: str, output_dir: Path) -> dict[str, str]:
    subset = result.loc[
        (result["analysis_level"] == analysis_level)
        & (result["pvalue"] < 0.05)
        & (result["overlap_count"] >= 2)
    ].copy()
    stem_name = f"{analysis_level}_ora_overview_heatmap"
    figure_stem = output_dir / "figures" / stem_name
    source_path = output_dir / "source_data" / f"{stem_name}.tsv"
    if subset.empty:
        atomic_tsv(source_path, subset)
        png, pdf = _empty_figure(
            f"{analysis_level} ORA overview", "No terms passed the frozen display rule.", figure_stem
        )
    else:
        term_key = subset["database"].astype(str) + "::" + subset["term_id"].astype(str)
        subset["term_key"] = term_key
        term_summary = (
            subset.groupby("term_key", observed=True)
            .agg(min_nominal_p=("pvalue", "min"), min_bh=("p_adjust_bh", "min"), group_count=("group", "nunique"), term_name=("term_name", "first"), database=("database", "first"))
            .reset_index()
            .sort_values(["group_count", "min_nominal_p"], ascending=[False, True])
            .head(30)
        )
        shown = subset.loc[subset["term_key"].isin(term_summary["term_key"])].copy()
        shown["minus_log10_nominal_p"] = -np.log10(shown["pvalue"].clip(lower=np.finfo(float).tiny))
        atomic_tsv(source_path, shown)
        groups = sorted(shown["group"].astype(str).unique(), key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
        terms = term_summary["term_key"].tolist()[::-1]
        matrix = shown.pivot_table(index="term_key", columns="group", values="minus_log10_nominal_p", aggfunc="max").reindex(index=terms, columns=groups).fillna(0)
        labels = []
        label_map = term_summary.set_index("term_key")
        for key in matrix.index:
            row = label_map.loc[key]
            label = f"{row['database']}: {row['term_name']}"
            labels.append(label if len(label) <= 58 else label[:57] + "…")
        height = min(9.0, max(4.5, 1.8 + 0.28 * len(matrix)))
        fig, ax = plt.subplots(figsize=(7.2, height))
        image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="magma", vmin=0)
        ax.set_xticks(np.arange(len(groups)), groups, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(labels)), labels)
        ax.set_xlabel("Cluster" if analysis_level == "cluster" else "Cell type")
        ax.set_title(f"Exploratory {analysis_level} marker ORA overview\nnominal display; annotation only")
        ax.tick_params(axis="both", labelsize=7)
        cbar = fig.colorbar(image, ax=ax, pad=0.02)
        cbar.set_label("−log10(nominal P); 0 = not displayed", fontsize=8)
        fig.tight_layout()
        png, pdf = _save_figure(fig, figure_stem)
    return {
        "title": f"{LEVEL_CN.get(analysis_level, analysis_level)} marker ORA总览热图",
        "png": str(png.relative_to(output_dir)),
        "pdf": str(pdf.relative_to(output_dir)),
        "source": str(source_path.relative_to(output_dir)),
        "interpretation": (
            f"实际显示 {shown['term_key'].nunique() if not subset.empty else 0} 个条目，覆盖 "
            f"{shown['group'].nunique() if not subset.empty else 0} 个群；颜色为−log10(nominal P)。"
            f"最突出条目为 {_top_theme_text(shown if not subset.empty else subset)}。"
            "条目先按覆盖群数、再按最小nominal P选取，完整BH保留在全表中，热图仅作注释总览。"
        ),
    }


def _resource_term_members(resources_dir: Path) -> dict[str, set[str]]:
    go = pd.read_csv(resources_dir / "go_term2gene.tsv", sep="\t", dtype=str)
    kegg = pd.read_csv(resources_dir / "kegg_term2gene.tsv", sep="\t", dtype=str)
    members: dict[str, set[str]] = {}
    for ontology in ("BP", "CC", "MF"):
        db = f"GO_{ontology}"
        for term_id, frame in go.loc[go["ontology"] == ontology].groupby("term_id", observed=True):
            members[f"{db}::{term_id}"] = set(frame["entrez_id"])
    for term_id, frame in kegg.groupby("term_id", observed=True):
        members[f"KEGG::{term_id}"] = set(frame["entrez_id"])
    return members


def _term_network(
    result: pd.DataFrame, analysis_level: str, output_dir: Path, resources_dir: Path
) -> tuple[dict[str, str], pd.DataFrame]:
    significant = result.loc[
        (result["analysis_level"] == analysis_level)
        & (result["pvalue"] < 0.05)
        & (result["overlap_count"] >= 2)
    ].copy()
    stem_name = f"{analysis_level}_ora_term_overlap_network"
    figure_stem = output_dir / "figures" / stem_name
    source_path = output_dir / "source_data" / f"{stem_name}.tsv"
    representatives = pd.DataFrame()
    if significant.empty:
        atomic_tsv(source_path, significant)
        png, pdf = _empty_figure(
            f"{analysis_level} term-overlap network", "No nominal-P-supported terms available.", figure_stem
        )
    else:
        significant["term_key"] = significant["database"].astype(str) + "::" + significant["term_id"].astype(str)
        nodes = (
            significant.groupby("term_key", observed=True)
            .agg(
                database=("database", "first"),
                term_id=("term_id", "first"),
                term_name=("term_name", "first"),
                min_nominal_p=("pvalue", "min"),
                min_bh=("p_adjust_bh", "min"),
                group_count=("group", "nunique"),
                groups=("group", lambda x: ";".join(sorted(set(map(str, x))))),
            )
            .reset_index()
            .sort_values(["group_count", "min_nominal_p"], ascending=[False, True])
            .head(40)
        )
        members = _resource_term_members(resources_dir)
        graph = nx.Graph()
        for row in nodes.itertuples(index=False):
            graph.add_node(row.term_key, **row._asdict())
        keys = nodes["term_key"].tolist()
        edge_rows: list[dict[str, Any]] = []
        for i, left in enumerate(keys):
            left_members = members.get(left, set())
            for right in keys[i + 1 :]:
                right_members = members.get(right, set())
                shared = left_members & right_members
                union = left_members | right_members
                jaccard = len(shared) / len(union) if union else 0.0
                if len(shared) >= 3 and jaccard >= 0.20:
                    graph.add_edge(left, right, jaccard=jaccard, shared_genes=len(shared))
                    edge_rows.append(
                        {
                            "record_type": "edge",
                            "source": left,
                            "target": right,
                            "jaccard": jaccard,
                            "shared_gene_count": len(shared),
                        }
                    )
        components = list(nx.connected_components(graph))
        component_map = {node: index + 1 for index, component in enumerate(components) for node in component}
        node_rows = nodes.copy()
        node_rows.insert(0, "record_type", "node")
        node_rows["component"] = node_rows["term_key"].map(component_map)
        node_rows["representative"] = False
        rep_indices = node_rows.groupby("component", observed=True)["min_nominal_p"].idxmin()
        node_rows.loc[rep_indices, "representative"] = True
        representatives = node_rows.loc[node_rows["representative"]].copy()
        combined = pd.concat([node_rows, pd.DataFrame(edge_rows)], ignore_index=True, sort=False)
        atomic_tsv(source_path, combined)

        fig, ax = plt.subplots(figsize=(7.2, 5.4))
        ax.axis("off")
        ordered_components = sorted(components, key=len, reverse=True)
        n_columns = max(1, math.ceil(math.sqrt(len(ordered_components))))
        n_rows = max(1, math.ceil(len(ordered_components) / n_columns))
        positions: dict[str, np.ndarray] = {}
        largest_component = max(map(len, ordered_components))
        for component_index, component in enumerate(ordered_components):
            subgraph = graph.subgraph(component)
            relative_scale = max(0.42, math.sqrt(len(component) / largest_component))
            local = nx.spring_layout(
                subgraph,
                seed=SEED + component_index,
                k=1.15 / math.sqrt(max(1, len(component))),
                iterations=200,
                scale=0.68 * relative_scale,
            )
            column = component_index % n_columns
            row = component_index // n_columns
            center = np.array(
                [1.75 * (column - (n_columns - 1) / 2), -1.75 * (row - (n_rows - 1) / 2)]
            )
            positions.update({node: np.asarray(point) + center for node, point in local.items()})
        edge_widths = [0.5 + 4 * graph.edges[edge]["jaccard"] for edge in graph.edges]
        nx.draw_networkx_edges(graph, positions, ax=ax, width=edge_widths, edge_color="#B8BEC4", alpha=0.65)
        node_colors = [DATABASE_COLORS.get(graph.nodes[node]["database"], "#777777") for node in graph.nodes]
        node_sizes = [70 + 35 * min(12, -math.log10(max(graph.nodes[node]["min_nominal_p"], 1e-300))) for node in graph.nodes]
        nx.draw_networkx_nodes(graph, positions, ax=ax, node_color=node_colors, node_size=node_sizes, edgecolors="white", linewidths=0.5)
        labels = {}
        for node in graph.nodes:
            if bool(node_rows.set_index("term_key").loc[node, "representative"]):
                name = str(graph.nodes[node]["term_name"])
                labels[node] = name if len(name) <= 32 else name[:31] + "…"
        label_positions = {
            node: np.asarray(positions[node]) + np.array([0.0, 0.22]) for node in labels
        }
        nx.draw_networkx_labels(
            graph,
            label_positions,
            labels=labels,
            ax=ax,
            font_size=6,
            verticalalignment="bottom",
        )
        ax.set_title(f"Exploratory {analysis_level} marker ORA term-overlap network\nnominal display; labels mark component representatives")
        fig.tight_layout()
        png, pdf = _save_figure(fig, figure_stem)
    entry = {
        "title": f"{LEVEL_CN.get(analysis_level, analysis_level)} marker ORA去冗余网络",
        "png": str(png.relative_to(output_dir)),
        "pdf": str(pdf.relative_to(output_dir)),
        "source": str(source_path.relative_to(output_dir)),
        "interpretation": (
            f"网络含 {len(nodes) if not significant.empty else 0} 个nominal-P条目、"
            f"{len(edge_rows) if not significant.empty else 0} 条重叠边和"
            f"{len(components) if not significant.empty else 0} 个连通分量；"
            f"最突出功能为 {_top_theme_text(significant)}。节点颜色区分数据库、大小表示−log10(nominal P)，"
            "边要求共享≥3基因且Jaccard≥0.20；完整BH保留在源表，网络只用于识别冗余功能主题。"
        ),
    }
    return entry, representatives


def _running_sum_es_from_positions(scores: np.ndarray, hit_positions: np.ndarray) -> tuple[float, np.ndarray]:
    n_total = len(scores)
    hit_positions = np.sort(np.asarray(hit_positions, dtype=int))
    n_hits = len(hit_positions)
    if n_hits == 0 or n_hits == n_total:
        return float("nan"), np.array([], dtype=int)
    weights = np.abs(scores[hit_positions])
    if float(weights.sum()) == 0:
        weights = np.ones_like(weights)
    normalized = weights / weights.sum()
    misses_before = hit_positions - np.arange(n_hits)
    running_after_hit = np.cumsum(normalized) - misses_before / (n_total - n_hits)
    running_before_hit = running_after_hit - normalized
    high_index = int(np.argmax(running_after_hit))
    low_index = int(np.argmin(running_before_hit))
    high = float(running_after_hit[high_index])
    low = float(running_before_hit[low_index])
    if abs(high) >= abs(low):
        return high, hit_positions[: high_index + 1]
    return low, hit_positions[low_index:]


def _running_sum_es(scores: np.ndarray, hits: np.ndarray) -> tuple[float, np.ndarray]:
    return _running_sum_es_from_positions(scores, np.flatnonzero(np.asarray(hits, dtype=bool)))


def _descriptive_resources(resources_dir: Path) -> tuple[pd.DataFrame, dict[str, tuple[pd.DataFrame, pd.DataFrame]]]:
    symbol_map = pd.read_csv(resources_dir / "symbol_to_entrez.tsv", sep="\t", dtype=str)
    go = pd.read_csv(resources_dir / "go_term2gene.tsv", sep="\t", dtype=str)
    go_names = pd.read_csv(resources_dir / "go_term2name.tsv", sep="\t", dtype=str)
    kegg = pd.read_csv(resources_dir / "kegg_term2gene.tsv", sep="\t", dtype=str)
    kegg_names = pd.read_csv(resources_dir / "kegg_term2name.tsv", sep="\t", dtype=str)
    resources = {
        "GO_BP": (go.loc[go["ontology"] == "BP", ["term_id", "entrez_id"]], go_names.loc[go_names["ontology"] == "BP"]),
        "GO_CC": (go.loc[go["ontology"] == "CC", ["term_id", "entrez_id"]], go_names.loc[go_names["ontology"] == "CC"]),
        "GO_MF": (go.loc[go["ontology"] == "MF", ["term_id", "entrez_id"]], go_names.loc[go_names["ontology"] == "MF"]),
        "KEGG": (kegg, kegg_names),
    }
    return symbol_map, resources


def _compute_descriptive_es(
    effects: pd.DataFrame,
    symbol_map: pd.DataFrame,
    resources: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    detection_threshold: float,
    variant: str,
) -> pd.DataFrame:
    selected = effects.loc[
        np.isfinite(effects["log2_ratio_case_minus_control"])
        & (effects[["case_detection_rate", "control_detection_rate"]].max(axis=1) >= detection_threshold)
        & (effects["cell_type"].astype(str) != "unknown")
    ].copy()
    selected = selected.merge(
        symbol_map[["gene_symbol", "entrez_id"]], left_on="gene", right_on="gene_symbol", how="inner"
    )
    rows: list[dict[str, Any]] = []
    for cell_type, group in selected.groupby("cell_type", observed=True):
        ranked = group.sort_values(
            ["log2_ratio_case_minus_control", "gene"], ascending=[False, True]
        ).drop_duplicates("entrez_id", keep="first")
        ranked_ids = ranked["entrez_id"].to_numpy(dtype=str)
        scores = ranked["log2_ratio_case_minus_control"].to_numpy(dtype=float)
        id_to_symbol = dict(zip(ranked["entrez_id"], ranked["gene"], strict=True))
        rank_position = {item: index for index, item in enumerate(ranked_ids)}
        ranked_set = set(ranked_ids)
        for database, (term2gene, term2name) in resources.items():
            filtered = term2gene.loc[term2gene["entrez_id"].isin(ranked_set)].drop_duplicates()
            name_map = dict(zip(term2name["term_id"], term2name["term_name"], strict=False))
            for term_id, members_frame in filtered.groupby("term_id", observed=True):
                members = set(members_frame["entrez_id"])
                if not MIN_GENE_SET_SIZE <= len(members) <= MAX_GENE_SET_SIZE:
                    continue
                hit_positions = np.array(sorted(rank_position[item] for item in members), dtype=int)
                es, leading_positions = _running_sum_es_from_positions(scores, hit_positions)
                leading_ids = ranked_ids[leading_positions]
                rows.append(
                    {
                        "cell_type": str(cell_type),
                        "database": database,
                        "term_id": term_id,
                        "term_name": name_map.get(term_id, term_id),
                        "ranked_gene_count": len(ranked_ids),
                        "gene_set_size_in_rank": len(members),
                        "descriptive_es": es,
                        "leading_edge_count": len(leading_ids),
                        "leading_edge_entrez_ids": ";".join(leading_ids),
                        "leading_edge_gene_symbols": ";".join(id_to_symbol[item] for item in leading_ids),
                        "detection_threshold": detection_threshold,
                        "analysis_variant": variant,
                        "ranking_metric": "PFOS-minus-WT two-library descriptive log2 CPM ratio",
                        "scope": "unpermuted, unnormalized descriptive ES; no condition P/FDR; exploratory localization only",
                    }
                )
    return pd.DataFrame(rows)


def _set_from_semicolon(value: object) -> set[str]:
    if pd.isna(value) or not str(value):
        return set()
    return set(str(value).split(";"))


def _detection_sensitivity(primary: pd.DataFrame, sensitivity: pd.DataFrame) -> pd.DataFrame:
    keys = ["cell_type", "database", "term_id", "term_name"]
    left = primary[keys + ["ranked_gene_count", "gene_set_size_in_rank", "descriptive_es", "leading_edge_entrez_ids"]].rename(
        columns={
            "ranked_gene_count": "primary_ranked_gene_count",
            "gene_set_size_in_rank": "primary_gene_set_size",
            "descriptive_es": "primary_descriptive_es",
            "leading_edge_entrez_ids": "primary_leading_edge_entrez_ids",
        }
    )
    right = sensitivity[keys + ["ranked_gene_count", "gene_set_size_in_rank", "descriptive_es", "leading_edge_entrez_ids"]].rename(
        columns={
            "ranked_gene_count": "detection05_ranked_gene_count",
            "gene_set_size_in_rank": "detection05_gene_set_size",
            "descriptive_es": "detection05_descriptive_es",
            "leading_edge_entrez_ids": "detection05_leading_edge_entrez_ids",
        }
    )
    merged = left.merge(right, on=keys, how="left", validate="one_to_one")
    valid = np.isfinite(merged["detection05_descriptive_es"])
    merged["direction_concordant"] = np.where(
        valid,
        np.sign(merged["primary_descriptive_es"]) == np.sign(merged["detection05_descriptive_es"]),
        np.nan,
    )
    merged["absolute_es_difference"] = (
        merged["detection05_descriptive_es"] - merged["primary_descriptive_es"]
    ).abs()
    merged["leading_edge_jaccard"] = [
        (len(a & b) / len(a | b) if a | b else np.nan)
        for a, b in zip(
            merged["primary_leading_edge_entrez_ids"].map(_set_from_semicolon),
            merged["detection05_leading_edge_entrez_ids"].map(_set_from_semicolon),
            strict=True,
        )
    ]
    merged["scope"] = "detection>=0.05 pathway-ES sensitivity; no P/FDR/NES or biological replication"
    return merged


def _balanced_cell_pathway_sensitivity(
    input_h5ad: Path,
    primary: pd.DataFrame,
    symbol_map: pd.DataFrame,
    resources: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
    repeats: int = 50,
    max_cells: int = 1000,
) -> pd.DataFrame:
    data = ad.read_h5ad(input_h5ad)
    if "counts" not in data.layers:
        raise ValueError("H5AD lacks SoupX-audited counts layer for balanced pathway sensitivity")
    if "condition" not in data.obs or "cell_type" not in data.obs:
        raise ValueError("H5AD lacks condition/cell_type for balanced pathway sensitivity")
    counts = data.layers["counts"].tocsr()
    var_index = pd.Series(np.arange(data.n_vars), index=data.var_names.astype(str))
    mapped = symbol_map.loc[symbol_map["gene_symbol"].isin(var_index.index), ["gene_symbol", "entrez_id"]].copy()
    mapped["var_index"] = mapped["gene_symbol"].map(var_index).astype(int)
    mapped = mapped.sort_values(["entrez_id", "gene_symbol"]).drop_duplicates("entrez_id")
    mapped_entrez = mapped["entrez_id"].to_numpy(dtype=str)
    mapped_symbols = mapped["gene_symbol"].to_numpy(dtype=str)
    mapped_var_index = mapped["var_index"].to_numpy(dtype=int)
    mapped_position = {entrez: index for index, entrez in enumerate(mapped_entrez)}
    term_member_positions: dict[str, np.ndarray] = {}
    for database, (term2gene, _) in resources.items():
        for term_id, frame in term2gene.groupby("term_id", observed=True):
            positions = sorted(
                {mapped_position[item] for item in frame["entrez_id"] if item in mapped_position}
            )
            term_member_positions[f"{database}::{term_id}"] = np.asarray(positions, dtype=np.int32)
    condition = data.obs["condition"].astype(str).to_numpy()
    cell_types = data.obs["cell_type"].astype(str).to_numpy()
    output_rows: list[dict[str, Any]] = []
    for cell_type, primary_terms in primary.groupby("cell_type", observed=True):
        case_indices = np.flatnonzero((cell_types == str(cell_type)) & (condition == "PFOS"))
        control_indices = np.flatnonzero((cell_types == str(cell_type)) & (condition == "WT"))
        balanced_n = min(max_cells, len(case_indices), len(control_indices))
        if balanced_n < 20:
            continue
        primary_terms = primary_terms.reset_index(drop=True)
        primary_leading = primary_terms["leading_edge_entrez_ids"].map(_set_from_semicolon).tolist()
        primary_es = primary_terms["descriptive_es"].to_numpy(dtype=float)
        es_values = np.full((len(primary_terms), repeats), np.nan, dtype=np.float32)
        leading_jaccard = np.full((len(primary_terms), repeats), np.nan, dtype=np.float32)
        rng = np.random.default_rng(SEED + sum(map(ord, str(cell_type))))
        for repeat in range(repeats):
            case_pick = rng.choice(case_indices, size=balanced_n, replace=False)
            control_pick = rng.choice(control_indices, size=balanced_n, replace=False)
            case_sum = np.asarray(counts[case_pick].sum(axis=0)).ravel()
            control_sum = np.asarray(counts[control_pick].sum(axis=0)).ravel()
            case_cpm = case_sum / max(float(case_sum.sum()), 1.0) * 1e6
            control_cpm = control_sum / max(float(control_sum.sum()), 1.0) * 1e6
            score = np.log2((case_cpm + 1.0) / (control_cpm + 1.0))
            detected = (case_sum[mapped_var_index] > 0) | (control_sum[mapped_var_index] > 0)
            ids = mapped_entrez[detected]
            symbols = mapped_symbols[detected]
            mapped_scores = score[mapped_var_index[detected]]
            order = np.lexsort((symbols, -mapped_scores))
            detected_base_positions = np.flatnonzero(detected)
            ranked_base_positions = detected_base_positions[order]
            ranked_ids = ids[order]
            ranked_scores = mapped_scores[order]
            base_to_rank = np.full(len(mapped_entrez), -1, dtype=np.int32)
            base_to_rank[ranked_base_positions] = np.arange(len(ranked_ids), dtype=np.int32)
            for row_index, term in enumerate(primary_terms.itertuples(index=False)):
                member_positions = term_member_positions.get(
                    f"{term.database}::{term.term_id}", np.empty(0, dtype=np.int32)
                )
                positions = base_to_rank[member_positions]
                positions = np.sort(positions[positions >= 0])
                if not MIN_GENE_SET_SIZE <= len(positions) <= MAX_GENE_SET_SIZE:
                    continue
                es, leading_positions = _running_sum_es_from_positions(ranked_scores, positions)
                es_values[row_index, repeat] = es
                repeat_leading = set(ranked_ids[leading_positions])
                union = primary_leading[row_index] | repeat_leading
                leading_jaccard[row_index, repeat] = len(primary_leading[row_index] & repeat_leading) / len(union) if union else np.nan
        for row_index, term in enumerate(primary_terms.itertuples(index=False)):
            values = es_values[row_index]
            jaccard = leading_jaccard[row_index]
            valid = np.isfinite(values)
            valid_jaccard = np.isfinite(jaccard)
            output_rows.append(
                {
                    "cell_type": str(cell_type),
                    "database": term.database,
                    "term_id": term.term_id,
                    "term_name": term.term_name,
                    "primary_descriptive_es": primary_es[row_index],
                    "balanced_cells_per_condition": balanced_n,
                    "requested_repeats": repeats,
                    "valid_repeats": int(valid.sum()),
                    "median_balanced_es": float(np.nanmedian(values)) if valid.any() else np.nan,
                    "q10_balanced_es": float(np.nanquantile(values, 0.10)) if valid.any() else np.nan,
                    "q90_balanced_es": float(np.nanquantile(values, 0.90)) if valid.any() else np.nan,
                    "direction_stability": float(np.mean(np.sign(values[valid]) == np.sign(primary_es[row_index]))) if valid.any() else np.nan,
                    "median_leading_edge_jaccard": float(np.nanmedian(jaccard)) if valid_jaccard.any() else np.nan,
                    "q10_leading_edge_jaccard": float(np.nanquantile(jaccard, 0.10)) if valid_jaccard.any() else np.nan,
                    "q90_leading_edge_jaccard": float(np.nanquantile(jaccard, 0.90)) if valid_jaccard.any() else np.nan,
                    "seed": SEED + sum(map(ord, str(cell_type))),
                    "scope": "equal-cell technical pathway-ES sensitivity; no P/FDR/NES or biological replication",
                }
            )
    del data
    return pd.DataFrame(output_rows)


def _render_descriptive_figures(
    result: pd.DataFrame,
    detection_comparison: pd.DataFrame,
    balanced: pd.DataFrame,
    output_dir: Path,
) -> list[dict[str, str]]:
    figures: list[dict[str, str]] = []
    for cell_type, frame in result.groupby("cell_type", observed=True):
        selected_rows = []
        for _, db_frame in frame.groupby("database", observed=True):
            selected_rows.append(db_frame.nlargest(3, "descriptive_es"))
            selected_rows.append(db_frame.nsmallest(3, "descriptive_es"))
        shown = pd.concat(selected_rows, ignore_index=True).drop_duplicates(["database", "term_id"])
        shown["display_label"] = shown.apply(_term_label, axis=1)
        shown = shown.sort_values("descriptive_es")
        stem_name = f"cell_type_{_slug(cell_type)}_descriptive_pathway_es"
        source_path = output_dir / "source_data" / "descriptive" / f"{stem_name}.tsv"
        atomic_tsv(source_path, shown)
        height = min(10.5, max(4.0, 1.6 + 0.29 * len(shown)))
        fig, ax = plt.subplots(figsize=(7.2, height))
        y = np.arange(len(shown))
        colors = [DATABASE_COLORS.get(item, "#777777") for item in shown["database"]]
        ax.axvline(0, color="#555555", linewidth=0.6)
        ax.hlines(y, 0, shown["descriptive_es"], colors=colors, linewidth=1.2)
        ax.scatter(shown["descriptive_es"], y, c=colors, s=26, zorder=3)
        ax.set_yticks(y, shown["display_label"])
        ax.set_xlabel("Unnormalized descriptive enrichment score")
        ax.set_title(f"Two-library descriptive pathway ranking — {cell_type}")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="both", labelsize=7)
        fig.tight_layout()
        png, pdf = _save_figure(fig, output_dir / "figures" / "descriptive" / stem_name)
        figures.append(
            {
                "title": f"{cell_type}细胞两文库PFOS−WT描述性通路排序",
                "png": str(png.relative_to(output_dir)),
                "pdf": str(pdf.relative_to(output_dir)),
                "source": str(source_path.relative_to(output_dir)),
                "interpretation": (
                    f"实际显示 {len(shown)} 个条目；颜色区分数据库，横轴及线段方向表示未标准化描述性ES的正负。"
                    f"绝对ES较突出的功能包括 {_top_theme_text(shown.assign(_abs_es=shown['descriptive_es'].abs()), p_column='_abs_es', ascending=False)}。"
                    "每个数据库在该细胞类型的PFOS−WT排名中选最高和最低3项；未做置换或NES标准化，"
                    "不产生P/FDR，也不支持动物层面的条件推断。"
                ),
            }
        )
    sensitivity_summary = (
        detection_comparison.groupby("cell_type", observed=True)
        .agg(
            detection05_terms_compared=("detection05_descriptive_es", lambda x: int(np.isfinite(x).sum())),
            detection05_direction_concordance=("direction_concordant", "mean"),
            detection05_median_leading_jaccard=("leading_edge_jaccard", "median"),
        )
        .reset_index()
        .merge(
            balanced.groupby("cell_type", observed=True)
            .agg(
                equal_cell_terms_compared=("valid_repeats", lambda x: int((x > 0).sum())),
                equal_cell_terms_valid_repeats_0=("valid_repeats", lambda x: int((x == 0).sum())),
                equal_cell_terms_valid_repeats_1_to_49=("valid_repeats", lambda x: int(x.between(1, 49).sum())),
                equal_cell_terms_valid_repeats_50=("valid_repeats", lambda x: int((x == 50).sum())),
                equal_cell_median_direction_stability=("direction_stability", "median"),
                equal_cell_median_leading_jaccard=("median_leading_edge_jaccard", "median"),
            )
            .reset_index(),
            on="cell_type",
            how="outer",
        )
    )
    sensitivity_summary["scope"] = "pathway-ES technical sensitivity only; no P/FDR/NES"
    sensitivity_stem = "descriptive_pathway_es_sensitivity_summary"
    sensitivity_source = output_dir / "source_data" / "descriptive" / f"{sensitivity_stem}.tsv"
    atomic_tsv(sensitivity_source, sensitivity_summary)
    metric_columns = [
        "detection05_direction_concordance",
        "detection05_median_leading_jaccard",
        "equal_cell_median_direction_stability",
        "equal_cell_median_leading_jaccard",
    ]
    matrix = sensitivity_summary.set_index("cell_type")[metric_columns]
    fig, ax = plt.subplots(figsize=(7.2, max(3.8, 1.8 + 0.45 * len(matrix))))
    image = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(
        np.arange(4),
        ["Detection≥5%\ndirection", "Detection≥5%\nleading edge", "Equal-cell\ndirection", "Equal-cell\nleading edge"],
    )
    ax.set_yticks(np.arange(len(matrix)), matrix.index.astype(str))
    ax.set_title("Two-library descriptive pathway ES sensitivity")
    ax.tick_params(axis="both", labelsize=7)
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("Concordance / stability (0–1)", fontsize=8)
    fig.tight_layout()
    png, pdf = _save_figure(fig, output_dir / "figures" / "descriptive" / sensitivity_stem)
    figures.append(
        {
            "title": "两文库描述性通路ES敏感性总览",
            "png": str(png.relative_to(output_dir)),
            "pdf": str(pdf.relative_to(output_dir)),
            "source": str(sensitivity_source.relative_to(output_dir)),
            "interpretation": (
                f"实际比较 {len(sensitivity_summary)} 个细胞类型。颜色依次表示检测率≥5%重算后的ES方向一致率和leading-edge Jaccard，"
                "以及50次等细胞数重采样中有效重复的ES方向稳定率和leading-edge Jaccard；"
                f"全部50次有效 {int((balanced['valid_repeats'] == 50).sum())} 项、1–49次有效 "
                f"{int(balanced['valid_repeats'].between(1, 49).sum())} 项、0次有效 {int((balanced['valid_repeats'] == 0).sum())} 项。"
                "分位数和稳定率以各term的valid_repeats为分母；全部为技术敏感性，不产生P/FDR/NES。"
            ),
        }
    )
    return figures


def _descriptive_ranked_pathways(
    input_h5ad: Path, source_run: Path, resources_dir: Path, output_dir: Path
) -> tuple[pd.DataFrame, list[dict[str, str]], list[str]]:
    effect_path = source_run / "05_descriptive" / "condition_descriptive_gene_effects.tsv"
    if not effect_path.is_file():
        return pd.DataFrame(), [], []
    effects = pd.read_csv(effect_path, sep="\t")
    required = {
        "cell_type",
        "gene",
        "log2_ratio_case_minus_control",
        "case_detection_rate",
        "control_detection_rate",
    }
    if not required <= set(effects.columns):
        raise ValueError(f"descriptive effect table lacks {sorted(required - set(effects.columns))}")
    symbol_map, resources = _descriptive_resources(resources_dir)
    result = _compute_descriptive_es(effects, symbol_map, resources, 0.0, "primary_any_detection")
    result_path = output_dir / "tables" / "descriptive_ranked_pathway_es.tsv"
    atomic_tsv(result_path, result)
    detection05 = _compute_descriptive_es(effects, symbol_map, resources, 0.05, "detection_rate_at_least_0.05")
    detection05_path = output_dir / "tables" / "descriptive_ranked_pathway_es_detection05.tsv"
    atomic_tsv(detection05_path, detection05)
    detection_comparison = _detection_sensitivity(result, detection05)
    detection_comparison_path = output_dir / "tables" / "descriptive_pathway_es_detection05_stability.tsv"
    atomic_tsv(detection_comparison_path, detection_comparison)
    balanced = _balanced_cell_pathway_sensitivity(input_h5ad, result, symbol_map, resources)
    balanced_path = output_dir / "tables" / "descriptive_pathway_es_equal_cell_stability.tsv"
    atomic_tsv(balanced_path, balanced)
    figures = _render_descriptive_figures(result, detection_comparison, balanced, output_dir)
    tables = [
        str(detection05_path.relative_to(output_dir)),
        str(detection_comparison_path.relative_to(output_dir)),
        str(balanced_path.relative_to(output_dir)),
    ]
    return result, figures, tables


def _validate_outputs(output_dir: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    missing: list[str] = []
    for figure in manifest["figures"]:
        for key in ("png", "pdf", "source"):
            if not (output_dir / figure[key]).is_file():
                missing.append(figure[key])
    for table in manifest["tables"]:
        if not (output_dir / table).is_file():
            missing.append(table)
    ora = pd.read_csv(
        output_dir / "tables" / "ora_full.tsv",
        sep="\t",
        dtype={"group": str},
        low_memory=False,
    )
    missing_term_names = int(ora["term_name"].isna().sum())
    if missing_term_names:
        raise RuntimeError(f"ORA contains {missing_term_names} missing term names")
    forbidden = [column for column in ora.columns if re.search(r"condition.*(pvalue|fdr)|pfos.*(pvalue|fdr)", column, re.I)]
    if forbidden:
        raise RuntimeError(f"condition-inference columns found in ORA: {forbidden}")
    status = pd.read_csv(
        output_dir / "tables" / "ora_status.tsv",
        sep="\t",
        dtype={"group": str},
        low_memory=False,
    )
    ora_family_count = int(
        ora[["analysis_level", "group", "database"]].drop_duplicates().shape[0]
    )
    status_family_count = int(len(status))
    if ora_family_count != status_family_count:
        raise RuntimeError(
            f"ORA/status family mismatch: {ora_family_count} versus {status_family_count}"
        )
    descriptive_tables = [
        "descriptive_ranked_pathway_es.tsv",
        "descriptive_ranked_pathway_es_detection05.tsv",
        "descriptive_pathway_es_detection05_stability.tsv",
        "descriptive_pathway_es_equal_cell_stability.tsv",
    ]
    for name in descriptive_tables:
        frame = pd.read_csv(output_dir / "tables" / name, sep="\t", low_memory=False)
        if frame["term_name"].isna().any():
            raise RuntimeError(f"{name} contains missing term names")
        forbidden_descriptive = [
            column
            for column in frame.columns
            if re.search(r"(^|_)(p|pvalue|padj|fdr|nes)($|_)", column, re.I)
        ]
        if forbidden_descriptive:
            raise RuntimeError(f"{name} contains forbidden inference columns: {forbidden_descriptive}")
    validation = {
        "status": "PASS" if not missing else "FAIL",
        "missing_paths": missing,
        "ora_rows": int(len(ora)),
        "ora_families": ora_family_count,
        "ora_status_families": status_family_count,
        "missing_ora_term_names": missing_term_names,
        "descriptive_tables_checked": descriptive_tables,
        "figure_count": len(manifest["figures"]),
        "paired_figure_sources": int(sum((output_dir / item["source"]).is_file() for item in manifest["figures"])),
        "checked_at_utc": _utc_now(),
    }
    if missing:
        raise RuntimeError(f"enrichment output validation missing {missing}")
    atomic_json(output_dir / "validation.json", validation)
    return validation


def _method_plan_text(resources_dir: Path) -> str:
    return f"""# PFOS islet marker enrichment method

- Input organism/reference: Mus musculus, GRCm39/Ensembl 115 expression reference; annotation resources are frozen separately.
- Marker query: positive marker log2FC >= {MARKER_LOG2FC} and nominal P < {MARKER_PVALUE}; adjusted P is retained but is not the marker selection rule.
- Background: genes with finite marker statistics and nonzero detection in the target or reference cells, recomputed for every cluster/cell-type comparison; each database further restricts this to mapped genes.
- ORA: R `stats::phyper` upper-tail hypergeometric probability; BH over every term in one group/database with {MIN_GENE_SET_SIZE}-{MAX_GENE_SET_SIZE} background genes, including zero-overlap terms. Raw P values are cross-checked against `clusterProfiler::enricher` for overlapping terms.
- GO: BP/CC/MF from `org.Mm.eg.db` GOALL/ONTOLOGYALL and `GO.db`; KEGG: mouse pathway links/names from KEGG REST, frozen in `{resources_dir}`.
- Inference boundary: marker ORA P/BH values support atlas cluster, cell-type, or beta-subcluster annotation only. They do not test PFOS versus WT and cells are not biological replicates.
- Descriptive PFOS-WT pathway ranking: full detectable cell-type gene ranking by two-library log2 CPM ratio, standard weighted running-sum ES without permutations or normalization. It has no P/FDR and cannot support a population condition effect.
- Descriptive sensitivity: repeat the full pathway ES after requiring maximum PFOS/WT detection rate >=0.05, then attempt pathway ES for 50 seed-fixed equal-cell subsamples per cell type (up to 1,000 cells per condition). Save the actual valid-repeat denominator, ES direction stability and leading-edge Jaccard; a repeat is not evaluable when its detected term size falls outside 10-{MAX_GENE_SET_SIZE}. Neither sensitivity creates P/FDR/NES.
- Display: the user-frozen project rule is nominal ORA P < 0.05 with >=2 marker genes. Plot labels say nominal/exploratory/annotation-only, and `[BH]` marks the subset also passing within-family BH. Complete test families and a nominal-versus-BH retention summary remain in tables. Network edges require >=3 shared genes and Jaccard >=0.20.

## Knowledge reuse

| knowledge_id | candidate_idea | decision | reason | adaptation | validation |
|---|---|---|---|---|---|
| 2026-08-28-network-pharmacology-core-api-loop | explicit species background for ORA | adopt | background universe is required for valid hypergeometric sampling | use each marker comparison's detectable/tested genes | mapping and family sizes saved |
| 2026-09-03-mfuzz-short-omics-time-series | full tables separate from display and adaptive figure height | adapt | prevents top-term figures from redefining the tested family | fixed nominal-P/overlap display rule and variable height | source TSV/PDF/PNG parity |
| nine small module scores | substitute module scores for enrichment | reject | modules do not represent GO/KEGG over-representation | none | no module score input read |
"""


def finalize_enrichment_display(output_dir: Path) -> dict[str, Any]:
    """Re-render figures from sealed numerical tables without recomputing statistics."""
    output_dir = Path(output_dir).resolve()
    snapshot_dir = output_dir / "provenance_snapshots"
    numerical_source = snapshot_dir / "extension_enrichment_numerical_as_used.py"
    numerical_execution = snapshot_dir / "execution_inputs_numerical_as_used.json"
    numerical_manifest = snapshot_dir / "report_manifest_numerical_as_used.json"
    numerical_validation = snapshot_dir / "validation_numerical_as_used.json"
    required_snapshots = [
        numerical_source,
        numerical_execution,
        numerical_manifest,
        numerical_validation,
    ]
    missing_snapshots = [str(path) for path in required_snapshots if not path.is_file()]
    if missing_snapshots:
        raise FileNotFoundError(f"display finalizer lacks numerical snapshots: {missing_snapshots}")
    numerical_inputs = json.loads(numerical_execution.read_text(encoding="utf-8"))
    numerical_source_entry = next(
        (
            item
            for item in numerical_inputs["files"]
            if item["path"].endswith("src/pfos_scrna/extension_enrichment.py")
        ),
        None,
    )
    if numerical_source_entry is None or sha256(numerical_source) != numerical_source_entry["sha256"]:
        raise RuntimeError("numerical source snapshot does not match the executed source hash")
    if numerical_inputs.get("verification_status") != "PASS" or numerical_inputs.get("changed"):
        raise RuntimeError("numerical execution snapshot was not input-verified")

    core_names = [
        "cluster_marker_selection.tsv",
        "cell_type_markers.tsv",
        "cell_type_marker_selection.tsv",
        "beta_cluster_marker_selection.tsv",
        "marker_query.tsv",
        "marker_detectable_universe.tsv",
        "marker_input_summary.tsv",
        "ora_full.tsv",
        "ora_status.tsv",
        "clusterprofiler_validation.tsv",
        "gene_id_mapping.tsv",
        "ora_nominal_bh_support_summary.tsv",
        "descriptive_ranked_pathway_es.tsv",
        "descriptive_ranked_pathway_es_detection05.tsv",
        "descriptive_pathway_es_detection05_stability.tsv",
        "descriptive_pathway_es_equal_cell_stability.tsv",
        "go_name_repair_audit.json",
    ]
    core_paths = [output_dir / "tables" / name for name in core_names]
    missing_core = [str(path) for path in core_paths if not path.is_file()]
    if missing_core:
        raise FileNotFoundError(f"display finalizer lacks numerical tables: {missing_core}")
    core_before = {str(path): sha256(path) for path in core_paths}

    resources_dir = Path(
        json.loads((output_dir / "provenance.json").read_text(encoding="utf-8"))["resources"]["resources_dir"]
    )
    r_script_entry = next(
        item for item in numerical_inputs["files"] if item["path"].endswith("scripts/extension_enrichment.R")
    )
    if sha256(Path(r_script_entry["path"])) != r_script_entry["sha256"]:
        raise RuntimeError("R numerical implementation changed before display finalization")
    atomic_text(output_dir / "method_plan.md", _method_plan_text(resources_dir))

    finalizer_input_path = output_dir / "execution_inputs.json"
    finalizer_snapshot = _freeze_execution_inputs(
        [
            *required_snapshots,
            *core_paths,
            output_dir / "provenance.json",
            Path(__file__),
            Path(r_script_entry["path"]),
            resources_dir / "resource_manifest.json",
            *_resource_files(resources_dir),
        ],
        finalizer_input_path,
    )

    ora = pd.read_csv(
        output_dir / "tables" / "ora_full.tsv",
        sep="\t",
        dtype={"group": str},
        low_memory=False,
    )
    status = pd.read_csv(
        output_dir / "tables" / "ora_status.tsv",
        sep="\t",
        dtype={"group": str},
        low_memory=False,
    )
    figures: list[dict[str, str]] = []
    analysis_levels = [
        level for level in ("cluster", "cell_type", "beta_cluster")
        if level in set(status["analysis_level"])
    ]
    for level in analysis_levels:
        groups = sorted(
            status.loc[status["analysis_level"] == level, "group"].astype(str).unique(),
            key=lambda value: (
                not value.isdigit(),
                int(value) if value.isdigit() else value,
            ),
        )
        for group in groups:
            figures.extend(_plot_group_ora(ora, level, group, output_dir))
        figures.append(_overview_heatmap(ora, level, output_dir))
        network_entry, representatives = _term_network(ora, level, output_dir, resources_dir)
        figures.append(network_entry)
        atomic_tsv(
            output_dir / "tables" / f"{level}_ora_representative_terms.tsv",
            representatives,
        )

    descriptive = pd.read_csv(
        output_dir / "tables" / "descriptive_ranked_pathway_es.tsv",
        sep="\t",
        low_memory=False,
    )
    detection = pd.read_csv(
        output_dir / "tables" / "descriptive_pathway_es_detection05_stability.tsv",
        sep="\t",
        low_memory=False,
    )
    balanced = pd.read_csv(
        output_dir / "tables" / "descriptive_pathway_es_equal_cell_stability.tsv",
        sep="\t",
        low_memory=False,
    )
    figures.extend(_render_descriptive_figures(descriptive, detection, balanced, output_dir))

    core_after = {str(path): sha256(path) for path in core_paths}
    if core_after != core_before:
        raise RuntimeError("display finalization changed sealed numerical tables")
    base_manifest = json.loads(numerical_manifest.read_text(encoding="utf-8"))
    metrics = dict(base_manifest["metrics"])
    metrics.update(
        {
            "ora_families": int(len(status)),
            "ora_nominal_display_terms": int(
                ((ora["pvalue"] < 0.05) & (ora["overlap_count"] >= 2)).sum()
            ),
            "descriptive_equal_cell_valid_repeats_0": int((balanced["valid_repeats"] == 0).sum()),
            "descriptive_equal_cell_valid_repeats_1_to_49": int(
                balanced["valid_repeats"].between(1, 49).sum()
            ),
            "descriptive_equal_cell_valid_repeats_50": int((balanced["valid_repeats"] == 50).sum()),
        }
    )
    provenance_paths = [
        "display_finalization.json",
        "provenance_snapshots/extension_enrichment_numerical_as_used.py",
        "provenance_snapshots/execution_inputs_numerical_as_used.json",
        "provenance_snapshots/report_manifest_numerical_as_used.json",
        "provenance_snapshots/validation_numerical_as_used.json",
    ]
    tables = list(dict.fromkeys([*base_manifest["tables"], *provenance_paths]))
    finalization = {
        "status": "pending",
        "scope": "display-only finalization; no marker, ORA, or descriptive ES recomputation",
        "started_at_utc": _utc_now(),
        "numerical_source_snapshot": str(numerical_source.relative_to(output_dir)),
        "numerical_source_sha256": sha256(numerical_source),
        "numerical_execution_snapshot": str(numerical_execution.relative_to(output_dir)),
        "numerical_execution_snapshot_sha256": sha256(numerical_execution),
        "finalizer_source": str(Path(__file__).resolve()),
        "finalizer_source_sha256": sha256(Path(__file__)),
        "sealed_numerical_tables_before": core_before,
        "sealed_numerical_tables_after": core_after,
        "sealed_numerical_tables_unchanged": True,
    }
    atomic_json(output_dir / "display_finalization.json", finalization)
    manifest = {
        **base_manifest,
        "status": "completed",
        "method": base_manifest["method"] + "; display-only finalizer from sealed numerical tables",
        "metrics": metrics,
        "figures": figures,
        "tables": tables,
        "display_finalization": {
            "status": "PASS",
            "record": "display_finalization.json",
            "numerical_tables_unchanged": True,
        },
    }
    validation = _validate_outputs(output_dir, manifest)
    manifest["validation"] = validation
    _verify_execution_inputs(finalizer_snapshot, finalizer_input_path)
    finalization.update(
        {
            "status": "PASS",
            "completed_at_utc": _utc_now(),
            "execution_inputs": "execution_inputs.json",
            "execution_inputs_sha256": sha256(finalizer_input_path),
            "validation": validation,
        }
    )
    atomic_json(output_dir / "display_finalization.json", finalization)
    atomic_json(output_dir / "report_manifest.json", manifest)
    return manifest


def run_enrichment(
    input_h5ad: Path,
    source_run: Path,
    output_dir: Path,
    resources_dir: Path,
) -> dict[str, Any]:
    """Run marker ORA and descriptive full-rank pathway localization.

    `source_run` is read-only. All new results are written under `output_dir`, and
    database resources are frozen under `resources_dir`.
    """
    input_h5ad = Path(input_h5ad).resolve()
    source_run = Path(source_run).resolve()
    output_dir = Path(output_dir).resolve()
    resources_dir = Path(resources_dir).resolve()
    if not input_h5ad.is_file():
        raise FileNotFoundError(input_h5ad)
    project_root = Path(__file__).resolve().parents[2]
    r_script = project_root / "scripts" / "extension_enrichment.R"
    if not r_script.is_file():
        raise FileNotFoundError(r_script)
    for folder in ("tables", "figures", "source_data", "logs"):
        (output_dir / folder).mkdir(parents=True, exist_ok=True)
    go_name_repair_baseline = _go_name_repair_baseline(output_dir, resources_dir)
    log_path = output_dir / "logs" / "enrichment_run.log"
    _write_log(log_path, "enrichment extension started")

    cluster_raw = _cluster_marker_table(source_run)
    cluster_table, cluster_query, cluster_universe, cluster_summary = _marker_contract(cluster_raw, "cluster")
    cluster_path = output_dir / "tables" / "cluster_marker_selection.tsv"
    atomic_tsv(cluster_path, cluster_table)

    cell_marker_path = output_dir / "tables" / "cell_type_markers.tsv"
    cell_raw = _compute_cell_type_markers(input_h5ad, cell_marker_path, log_path)
    cell_table, cell_query, cell_universe, cell_summary = _marker_contract(cell_raw, "cell_type")
    cell_selection_path = output_dir / "tables" / "cell_type_marker_selection.tsv"
    atomic_tsv(cell_selection_path, cell_table)

    marker_tables = [cluster_query, cell_query]
    universe_tables = [cluster_universe, cell_universe]
    summary_tables = [cluster_summary, cell_summary]
    beta_input_paths: list[Path] = []
    beta_raw, beta_input_paths = _completed_beta_cluster_markers(output_dir.parent)
    beta_selection_path = output_dir / "tables" / "beta_cluster_marker_selection.tsv"
    if beta_raw is not None:
        beta_table, beta_query, beta_universe, beta_summary = _marker_contract(beta_raw, "beta_cluster")
        atomic_tsv(beta_selection_path, beta_table)
        marker_tables.append(beta_query)
        universe_tables.append(beta_universe)
        summary_tables.append(beta_summary)
        _write_log(log_path, f"included completed beta-cluster marker layer: {len(beta_table)} rows")
    query = pd.concat(marker_tables, ignore_index=True)
    universe = pd.concat(universe_tables, ignore_index=True)
    summary = pd.concat(summary_tables, ignore_index=True)
    query_path = output_dir / "tables" / "marker_query.tsv"
    universe_path = output_dir / "tables" / "marker_detectable_universe.tsv"
    summary_path = output_dir / "tables" / "marker_input_summary.tsv"
    atomic_tsv(query_path, query)
    atomic_tsv(universe_path, universe)
    atomic_tsv(summary_path, summary)

    resource_manifest = _ensure_resources(
        resources_dir,
        universe["gene_symbol"],
        r_script,
        log_path,
        project_root,
    )
    execution_inputs_path = output_dir / "execution_inputs.json"
    execution_snapshot = _freeze_execution_inputs(
        [
            input_h5ad,
            source_run / "03_atlas" / "cluster_markers.tsv",
            source_run / "05_descriptive" / "condition_descriptive_gene_effects.tsv",
            cell_marker_path,
            query_path,
            universe_path,
            Path(__file__),
            r_script,
            project_root / "pyproject.toml",
            project_root / "uv.lock",
            resources_dir / "resource_manifest.json",
            *_resource_files(resources_dir),
            *beta_input_paths,
        ],
        execution_inputs_path,
    )
    _run_r(
        [
            "Rscript",
            str(r_script),
            "run-ora",
            "--resources-dir",
            str(resources_dir),
            "--query",
            str(query_path),
            "--universe",
            str(universe_path),
            "--output-dir",
            str(output_dir / "tables"),
            "--min-size",
            str(MIN_GENE_SET_SIZE),
            "--max-size",
            str(MAX_GENE_SET_SIZE),
        ],
        log_path,
        project_root,
    )
    ora = pd.read_csv(output_dir / "tables" / "ora_full.tsv", sep="\t", dtype={"group": str})
    _audit_go_name_repair(go_name_repair_baseline, output_dir, resources_dir)
    status = pd.read_csv(output_dir / "tables" / "ora_status.tsv", sep="\t", dtype={"group": str})
    support_summary = (
        ora.groupby(["analysis_level", "group", "database"], observed=True)
        .agg(
            tested_terms=("term_id", "nunique"),
            terms_with_overlap=("overlap_count", lambda x: int((x > 0).sum())),
            nominal_p_lt_0_05=("passes_nominal_p_0_05", "sum"),
            bh_lt_0_05=("passes_bh_0_05", "sum"),
        )
        .reset_index()
    )
    support_summary["bh_retention_among_nominal"] = np.where(
        support_summary["nominal_p_lt_0_05"] > 0,
        support_summary["bh_lt_0_05"] / support_summary["nominal_p_lt_0_05"],
        np.nan,
    )
    support_summary["display_rule"] = "user-frozen nominal ORA P<0.05 and overlap>=2"
    support_summary["bh_role"] = "multiplicity sensitivity summary; never relabeled as nominal significance"
    atomic_tsv(output_dir / "tables" / "ora_nominal_bh_support_summary.tsv", support_summary)
    figures: list[dict[str, str]] = []
    analysis_levels = ["cluster", "cell_type"] + (["beta_cluster"] if beta_raw is not None else [])
    for level in analysis_levels:
        for group in sorted(status.loc[status["analysis_level"] == level, "group"].astype(str).unique(), key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x)):
            figures.extend(_plot_group_ora(ora, level, group, output_dir))
        figures.append(_overview_heatmap(ora, level, output_dir))
        network_entry, representatives = _term_network(ora, level, output_dir, resources_dir)
        figures.append(network_entry)
        atomic_tsv(output_dir / "tables" / f"{level}_ora_representative_terms.tsv", representatives)

    descriptive, descriptive_figures, descriptive_sensitivity_tables = _descriptive_ranked_pathways(
        input_h5ad, source_run, resources_dir, output_dir
    )
    figures.extend(descriptive_figures)
    table_paths = [
        "tables/cluster_marker_selection.tsv",
        "tables/cell_type_markers.tsv",
        "tables/cell_type_marker_selection.tsv",
        "tables/marker_query.tsv",
        "tables/marker_detectable_universe.tsv",
        "tables/marker_input_summary.tsv",
        "tables/ora_full.tsv",
        "tables/ora_status.tsv",
        "tables/clusterprofiler_validation.tsv",
        "tables/gene_id_mapping.tsv",
        "tables/ora_nominal_bh_support_summary.tsv",
        "tables/cluster_ora_representative_terms.tsv",
        "tables/cell_type_ora_representative_terms.tsv",
    ]
    if beta_raw is not None:
        table_paths.extend(
            ["tables/beta_cluster_marker_selection.tsv", "tables/beta_cluster_ora_representative_terms.tsv"]
        )
    if (output_dir / "tables" / "go_name_repair_audit.json").is_file():
        table_paths.append("tables/go_name_repair_audit.json")
    if not descriptive.empty:
        table_paths.append("tables/descriptive_ranked_pathway_es.tsv")
        table_paths.extend(descriptive_sensitivity_tables)

    atomic_text(output_dir / "method_plan.md", _method_plan_text(resources_dir))

    provenance = {
        "status": "completed",
        "started_or_recorded_at_utc": _utc_now(),
        "input_h5ad": str(input_h5ad),
        "input_h5ad_sha256": sha256(input_h5ad),
        "cluster_markers": str(source_run / "03_atlas" / "cluster_markers.tsv"),
        "cluster_markers_sha256": sha256(source_run / "03_atlas" / "cluster_markers.tsv"),
        "source_run": str(source_run),
        "output_dir": str(output_dir),
        "resources": resource_manifest,
        "parameters": {
            "marker_min_log2fc": MARKER_LOG2FC,
            "marker_nominal_pvalue": MARKER_PVALUE,
            "min_gene_set_size": MIN_GENE_SET_SIZE,
            "max_gene_set_size": MAX_GENE_SET_SIZE,
            "ora_bh_family": "analysis_level x group x database, all terms in size range",
            "figure_dpi": FIGURE_DPI,
            "seed": SEED,
        },
        "software": {
            "python": os.sys.version.split()[0],
            "anndata": version("anndata"),
            "scanpy": version("scanpy"),
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "networkx": nx.__version__,
        },
        "inference_boundary": "n=1 library per condition; no PFOS-WT condition P/FDR",
    }
    atomic_json(output_dir / "provenance.json", provenance)
    table_paths.extend(["method_plan.md", "provenance.json", "execution_inputs.json"])

    equal_cell = (
        pd.read_csv(
            output_dir / "tables" / "descriptive_pathway_es_equal_cell_stability.tsv",
            sep="\t",
            low_memory=False,
        )
        if not descriptive.empty
        else pd.DataFrame()
    )
    metrics = {
        "cluster_groups": int(summary.loc[summary["analysis_level"] == "cluster", "group"].nunique()),
        "cell_types": int(summary.loc[summary["analysis_level"] == "cell_type", "group"].nunique()),
        "beta_clusters": int(summary.loc[summary["analysis_level"] == "beta_cluster", "group"].nunique()),
        "marker_queries": int(len(query)),
        "ora_rows_full_family": int(len(ora)),
        "ora_families": int(len(status)),
        "ora_fdr_terms": int((ora["p_adjust_bh"] < 0.05).sum()),
        "ora_nominal_p_terms": int((ora["pvalue"] < 0.05).sum()),
        "descriptive_es_rows": int(len(descriptive)),
        "descriptive_detection05_stability_rows": int(
            len(pd.read_csv(output_dir / "tables" / "descriptive_pathway_es_detection05_stability.tsv", sep="\t"))
        ) if not descriptive.empty else 0,
        "descriptive_equal_cell_stability_rows": int(
            len(equal_cell)
        ) if not descriptive.empty else 0,
        "descriptive_equal_cell_valid_repeats_0": int((equal_cell["valid_repeats"] == 0).sum()) if not equal_cell.empty else 0,
        "descriptive_equal_cell_valid_repeats_1_to_49": int(equal_cell["valid_repeats"].between(1, 49).sum()) if not equal_cell.empty else 0,
        "descriptive_equal_cell_valid_repeats_50": int((equal_cell["valid_repeats"] == 50).sum()) if not equal_cell.empty else 0,
        "condition_inference_p_or_fdr_columns": 0,
    }
    manifest = {
        "status": "completed",
        "method": (
            "comparison-specific detectable-background marker ORA using frozen mouse GO/KEGG; "
            "R hypergeometric/BH with clusterProfiler cross-check; optional unpermuted descriptive pathway ES"
        ),
        "metrics": metrics,
        "figures": figures,
        "tables": table_paths,
        "scope": "cluster/cell-type functional annotation and two-library descriptive pathway localization",
        "inference_boundary": "n=1 library per condition; marker P/BH is annotation-only; descriptive ES has no P/FDR",
    }
    validation = _validate_outputs(output_dir, manifest)
    manifest["validation"] = validation
    _verify_execution_inputs(execution_snapshot, execution_inputs_path)
    atomic_json(output_dir / "report_manifest.json", manifest)
    _write_log(log_path, f"enrichment extension completed: {metrics}")
    return manifest

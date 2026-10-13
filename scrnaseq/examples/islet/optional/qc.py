from __future__ import annotations

import logging
import math
import subprocess
from pathlib import Path
import re

import anndata as ad
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import scrublet as scr

from .cellranger import cellranger_outs
from .runtime import RunContext, atomic_json, atomic_tsv, run_command
from .plotting import save_figure
from .figure_style import COLORS, PROFILE


LOGGER = logging.getLogger("pfos_scrna")
HORMONES = ("Ins1", "Ins2", "Gcg", "Sst", "Ppy", "Iapp")


def _mad(values: np.ndarray) -> float:
    median = np.nanmedian(values)
    return float(np.nanmedian(np.abs(values - median)))


def expected_doublet_rate(recovered_cells: int) -> float:
    recovered = np.array([1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 15000, 20000], dtype=float)
    rates = np.array([0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.061, 0.076, 0.113, 0.151], dtype=float)
    return float(np.clip(np.interp(recovered_cells, recovered, rates), 0.004, 0.15))


def _soupx_counts(context: RunContext, library_id: str, data: ad.AnnData) -> tuple[sp.csr_matrix, dict[str, object]]:
    outs = cellranger_outs(context, library_id)
    exchange = context.path("02_qc", library_id, "soupx")
    exchange.mkdir(parents=True, exist_ok=True)
    command = [context.config.tool("rscript"), context.config.root / "scripts" / "run_soupx.R", outs / "raw_feature_bc_matrix", outs / "filtered_feature_bc_matrix", exchange]
    try:
        run_command(command, log_file=context.path("logs", f"soupx_{library_id}.log"))
        matrix = scipy.io.mmread(exchange / "corrected.mtx").tocsr().T
        features = pd.read_csv(exchange / "features.tsv", sep="\t", header=None, dtype=str)
        barcodes = pd.read_csv(exchange / "barcodes.tsv", header=None, dtype=str)[0].tolist()
        ids = features[0].tolist()
        data_ids = data.var["gene_ids"].astype(str).tolist()
        if ids != data_ids or barcodes != data.obs_names.tolist() or matrix.shape != data.shape:
            raise RuntimeError("SoupX 输出与 filtered AnnData 索引不一致")
        if matrix.data.size and (matrix.data.min() < 0 or not np.allclose(matrix.data, np.rint(matrix.data))):
            raise RuntimeError("SoupX corrected counts 不是非负整数")
        increase = matrix - data.layers["counts_raw"]
        if (increase.data > 0).any() or matrix.sum() > data.layers["counts_raw"].sum():
            raise RuntimeError("SoupX corrected counts 出现大于 raw 的单元或总计数")
        summary = pd.read_csv(exchange / "summary.tsv", sep="\t").iloc[0].to_dict()
        for key in ("rho_median", "rho_mean"):
            value = float(summary[key])
            if not np.isfinite(value) or not 0 <= value <= 1:
                raise RuntimeError(f"SoupX {key} 无效: {value}")
        summary.update({"status": "completed", "fallback": False})
        return matrix.astype(np.int64), summary
    except Exception as exc:
        LOGGER.warning("SoupX %s 自动估计失败，明确保留 raw counts: %s", library_id, exc)
        return data.layers["counts_raw"].copy(), {"status": "failed_raw_retained", "fallback": False, "reason": str(exc)}


def _scrublet(data: ad.AnnData, seed: int, chemistry: str) -> dict[str, object]:
    try:
        if not chemistry.startswith("SC3Pv3"):
            raise RuntimeError(f"没有为检测 chemistry 冻结 doublet-rate 表: {chemistry}")
        rate = expected_doublet_rate(data.n_obs)
        runner = scr.Scrublet(data.layers["counts_raw"], expected_doublet_rate=rate, random_state=seed)
        max_pcs = max(2, min(30, data.n_obs - 2, data.n_vars - 2))
        score, prediction = runner.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=max_pcs, verbose=False)
        data.obs["doublet_score"] = score
        data.obs["predicted_doublet"] = prediction.astype(bool)
        return {"status": "completed", "chemistry": chemistry, "rate_source": "10x_SC3Pv3_recovered_cell_table", "expected_rate": rate, "threshold": float(runner.threshold_), "predicted": int(prediction.sum())}
    except Exception as exc:
        data.obs["doublet_score"] = np.nan
        data.obs["predicted_doublet"] = False
        return {"status": "failed_no_doublets_removed", "chemistry": chemistry, "reason": str(exc)}


def _hormone_audit(data: ad.AnnData, library_id: str) -> pd.DataFrame:
    rows = []
    for gene in HORMONES:
        if gene not in data.var_names:
            rows.append({"library_id": library_id, "gene": gene, "present": False})
            continue
        index = data.var_names.get_loc(gene)
        raw = np.asarray(data.layers["counts_raw"][:, index].toarray()).ravel()
        corrected = np.asarray(data.layers["counts"][:, index].toarray()).ravel()
        rows.append({"library_id": library_id, "gene": gene, "present": True, "raw_total": int(raw.sum()), "corrected_total": int(corrected.sum()), "raw_detection": float((raw > 0).mean()), "corrected_detection": float((corrected > 0).mean()), "removed_fraction": float(1 - corrected.sum() / raw.sum()) if raw.sum() else 0.0})
    return pd.DataFrame(rows)


def _barcode_rank(raw_h5: Path) -> np.ndarray:
    with h5py.File(raw_h5, "r") as handle:
        indptr = handle["matrix/indptr"][:]
        dataset = handle["matrix/data"]
        totals = np.zeros(len(indptr) - 1, dtype=np.int64)
        block = 100_000
        for first in range(0, len(totals), block):
            last = min(first + block, len(totals))
            offsets = indptr[first:last + 1]
            values = dataset[int(offsets[0]):int(offsets[-1])]
            cumulative = np.concatenate([np.array([0], dtype=np.int64), np.cumsum(values, dtype=np.int64)])
            local = offsets - offsets[0]
            totals[first:last] = cumulative[local[1:]] - cumulative[local[:-1]]
    return np.sort(totals[totals > 0])[::-1]


def draw_qc_source(frame: pd.DataFrame, stem: str, config: WorkflowConfig, called: int, target: Path) -> None:
    library = re.sub(r"_(barcode_rank|called_cell_qc|doublet_audit)$", "", stem)
    if stem.endswith("barcode_rank"):
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        ax.loglog(frame["barcode_rank"], frame["umi_counts"], color=COLORS[stem.split("_", 1)[0]], linewidth=1)
        # 使用同字号的普通数值，避免默认数学上标缩小成4.9pt。
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
        ax.axvline(called, color=PROFILE["palettes"]["neutral"]["dark"], linestyle="--", linewidth=1, label=f"Cell Ranger called={called:,}")
        ax.set(xlabel="Barcode rank", ylabel="UMI counts", title=f"{library} barcode rank")
        ax.legend(frameon=False)
    elif stem.endswith("called_cell_qc"):
        fig, axes = plt.subplots(2, 2, figsize=(9, 7))
        axes[0, 0].hist(np.log10(frame["total_counts"] + 1), bins=60, color=COLORS[stem.split("_", 1)[0]])
        axes[0, 0].set_xlabel("log10(UMI+1)")
        axes[0, 1].hist(frame["n_genes_by_counts"], bins=60, color=COLORS[stem.split("_", 1)[0]])
        axes[0, 1].set_xlabel("Detected genes")
        axes[1, 0].hist(frame["pct_counts_mt"], bins=60, color=COLORS[stem.split("_", 1)[0]])
        axes[1, 0].axvline(float(config.data["qc"]["max_pct_mt"]), color="black", linestyle="--")
        axes[1, 0].set_xlabel("Mitochondrial counts (%)")
        axes[1, 1].scatter(np.log10(frame["total_counts"] + 1), frame["n_genes_by_counts"], s=3, alpha=.35, linewidths=0, rasterized=True, color=COLORS[stem.split("_", 1)[0]])
        axes[1, 1].set(xlabel="log10(UMI+1)", ylabel="Detected genes")
        fig.suptitle(f"{library} called-cell QC")
    else:
        fig, axes = plt.subplots(1, 2, figsize=(9, 4))
        axes[0].hist(frame["doublet_score"].dropna(), bins=60, color=COLORS[stem.split("_", 1)[0]])
        axes[0].set(xlabel="Scrublet score", ylabel="Called cells")
        for predicted, color, marker in [(False, COLORS[stem.split("_", 1)[0]], "o"), (True, PROFILE["palettes"]["accent"]["warm"], "^")]:
            subset = frame.loc[frame["predicted_doublet"].astype(bool).eq(predicted)]
            axes[1].scatter(np.log10(subset["total_counts"] + 1), subset["doublet_score"], color=color, marker=marker, label="Predicted doublet" if predicted else "Other called cells", s=4, alpha=.5, linewidths=0, rasterized=True)
        axes[1].legend(frameon=False)
        axes[1].set(xlabel="log10(UMI+1)", ylabel="Scrublet score")
        fig.suptitle(f"{library} doublet audit")
    save_figure(fig, target)



def _qc_figures(context: RunContext, library_id: str, data: ad.AnnData, raw_h5: Path) -> int:
    ranks = _barcode_rank(raw_h5)
    sources = {
        "barcode_rank": pd.DataFrame({"barcode_rank": np.arange(1, len(ranks) + 1), "umi_counts": ranks}),
        "called_cell_qc": pd.DataFrame({"barcode": data.obs_names.astype(str), **{column: data.obs[column].to_numpy() for column in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb")}}),
        "doublet_audit": pd.DataFrame({"barcode": data.obs_names.astype(str), **{column: data.obs[column].to_numpy() for column in ("total_counts", "doublet_score", "predicted_doublet")}}),
    }
    for kind, frame in sources.items():
        stem = f"{library_id}_{kind}"
        draw_qc_source(frame, stem, context.config, data.n_obs, context.path("08_figures", stem))
        atomic_tsv(context.path("07_source_data", f"{stem}.tsv"), frame)
    return len(ranks)


def _sex_marker_audit(data: ad.AnnData, library_id: str) -> pd.DataFrame:
    rows = []
    for gene in ("Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Xist"):
        if gene not in data.var_names:
            rows.append({"library_id": library_id, "gene": gene, "present": False})
            continue
        values = np.asarray(data.layers["counts"][:, data.var_names.get_loc(gene)].toarray()).ravel()
        rows.append({"library_id": library_id, "gene": gene, "present": True, "total_count": int(values.sum()), "detection_rate": float((values > 0).mean())})
    return pd.DataFrame(rows)


def _qc_one(context: RunContext, row: object) -> ad.AnnData:
    library_id = row.library_id
    outs = cellranger_outs(context, library_id)
    data = sc.read_10x_h5(outs / "filtered_feature_bc_matrix.h5", gex_only=True)
    data.var_names_make_unique()
    data.X = data.X.tocsr().astype(np.int64)
    data.layers["counts_raw"] = data.X.copy()
    data.layers["counts"], soupx = _soupx_counts(context, library_id, data)
    atomic_json(context.path("02_qc", library_id, "soupx_status.json"), soupx)
    atomic_tsv(context.path("02_qc", library_id, "soupx_hormone_audit.tsv"), _hormone_audit(data, library_id))
    symbols = data.var_names.astype(str)
    data.var["gene_symbol"] = symbols
    data.var["mt"] = symbols.str.startswith("mt-")
    data.var["ribo"] = symbols.str.match(r"^Rp[sl]")
    data.var["hb"] = symbols.str.match(r"^Hb[ab]")
    data.X = data.layers["counts"]
    sc.pp.calculate_qc_metrics(data, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None, log1p=False)
    housekeeping = [gene for gene in context.config.data["qc"]["housekeeping"] if gene in data.var_names]
    if housekeeping:
        data.obs["housekeeping_detected"] = np.asarray((data[:, housekeeping].X > 0).sum(axis=1)).ravel()
    else:
        data.obs["housekeeping_detected"] = 0
    chemistry_table = pd.read_csv(context.path("01_cellranger", "chemistry_detection.tsv"), sep="\t", dtype=str)
    detected_chemistry = chemistry_table.set_index("library_id").loc[library_id, "detected_chemistry"]
    scrublet = _scrublet(data, int(context.config.data["project"]["seed"]), detected_chemistry)
    atomic_json(context.path("02_qc", library_id, "scrublet_status.json"), scrublet)
    raw_nonzero_barcodes = _qc_figures(context, library_id, data, outs / "raw_feature_bc_matrix.h5")
    atomic_tsv(context.path("02_qc", library_id, "sex_marker_audit.tsv"), _sex_marker_audit(data, library_id))
    atomic_tsv(context.path("02_qc", library_id, "cell_calling_audit.tsv"), pd.DataFrame([{"library_id": library_id, "raw_nonzero_barcodes": int(raw_nonzero_barcodes), "cellranger_called_cells": int(data.n_obs), "called_fraction_of_nonzero_barcodes": float(data.n_obs / max(raw_nonzero_barcodes, 1)), "scope": "cell-calling audit; Cell Ranger call retained before downstream QC"}]))

    qc = context.config.data["qc"]
    log_counts = np.log1p(data.obs["total_counts"].to_numpy())
    log_genes = np.log1p(data.obs["n_genes_by_counts"].to_numpy())
    nmad = float(qc["mad_n"])
    count_low = max(float(qc["min_counts"]), math.expm1(np.nanmedian(log_counts) - nmad * _mad(log_counts)))
    count_high = math.expm1(np.nanmedian(log_counts) + nmad * _mad(log_counts))
    gene_low = max(float(qc["min_genes"]), math.expm1(np.nanmedian(log_genes) - nmad * _mad(log_genes)))
    gene_high = math.expm1(np.nanmedian(log_genes) + nmad * _mad(log_genes))
    data.obs["high_count_outlier"] = data.obs["total_counts"] > count_high
    data.obs["qc_pass"] = (
        data.obs["total_counts"].between(count_low, count_high)
        & data.obs["n_genes_by_counts"].between(gene_low, gene_high)
        & (data.obs["pct_counts_mt"] <= float(qc["max_pct_mt"]))
        & (data.obs["pct_counts_hb"] <= float(qc["max_pct_hb"]))
        & (data.obs["housekeeping_detected"] >= int(qc["min_housekeeping"]))
        & ~data.obs["predicted_doublet"]
    )
    thresholds = pd.DataFrame([{"library_id": library_id, "count_low": count_low, "count_high": count_high, "gene_low": gene_low, "gene_high": gene_high, "max_pct_mt": qc["max_pct_mt"], "max_pct_hb": qc["max_pct_hb"], "min_housekeeping": qc["min_housekeeping"], "housekeeping_role": "audit_only" if int(qc["min_housekeeping"]) == 0 else "filter", "mad_n": nmad, "thresholds_frozen_before_condition_view": True}])
    atomic_tsv(context.path("02_qc", library_id, "qc_thresholds.tsv"), thresholds)
    sensitivity = pd.DataFrame([{"library_id": library_id, "max_pct_mt": threshold, "cells_passing_mt_only": int((data.obs["pct_counts_mt"] <= threshold).sum()), "fraction_passing_mt_only": float((data.obs["pct_counts_mt"] <= threshold).mean())} for threshold in qc["mt_sensitivity"]])
    atomic_tsv(context.path("02_qc", library_id, "mt_sensitivity.tsv"), sensitivity)
    for column, value in row._asdict().items():
        data.obs[column] = detected_chemistry if column == "chemistry" else str(value)
    data.obs_names = pd.Index([f"{library_id}|{barcode}" for barcode in data.obs_names])
    data.obs_names.name = "cell_id"
    data.X = data.layers["counts"]
    data.write_h5ad(context.path("02_qc", library_id, "called_cells_qc.h5ad"), compression="gzip")
    return data[data.obs["qc_pass"]].copy()


def run_qc(context: RunContext) -> Path:
    output = context.path("02_qc", "qc_combined.h5ad")
    figure_stems = [f"{row.library_id}_{name}" for row in context.config.samples().itertuples(index=False) for name in ("barcode_rank", "called_cell_qc", "doublet_audit")]
    figure_outputs = [context.path("08_figures", f"{stem}.{extension}") for stem in figure_stems for extension in ("png", "pdf")]
    source_outputs = [context.path("07_source_data", f"{stem}.tsv") for stem in figure_stems]
    with context.stage("qc", [output, *figure_outputs, *source_outputs]) as execute:
        if not execute:
            return output
        filtered = [_qc_one(context, row) for row in context.config.samples().itertuples(index=False)]
        combined = ad.concat(filtered, join="inner", merge="same", uns_merge="unique", index_unique=None)
        if not combined.obs_names.is_unique:
            raise RuntimeError("拼接后 cell barcode 不唯一")
        combined.X = combined.layers["counts"].copy()
        combined.write_h5ad(output, compression="gzip")
        summary = combined.obs.groupby("library_id", observed=True).agg(cells=("qc_pass", "size"), median_counts=("total_counts", "median"), median_genes=("n_genes_by_counts", "median"), median_pct_mt=("pct_counts_mt", "median")).reset_index()
        atomic_tsv(context.path("02_qc", "qc_summary.tsv"), summary)
    return output

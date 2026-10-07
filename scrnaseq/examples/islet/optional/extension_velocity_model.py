"""Per-library beta velocity diagnostics; no condition-level dynamics test."""
from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path

import anndata as ad
from anndata.io import read_loom
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scvelo as scv
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

from .extension_velocity import _fingerprints, _verify_fingerprints, loom_barcode
from .plotting import save_figure
from .runtime import atomic_json, atomic_text, atomic_tsv


# Operational screening criteria, fixed before inspecting the splice matrices.
# Passing means "eligible for limited exploratory projection", not validated fate.
HORMONE_GENES = ("Ins1", "Ins2", "Gcg", "Sst", "Ppy", "Iapp")
CRITERIA = {"min_shared_counts": 20, "min_cells_each_layer": 50, "min_coverage_genes": 100, "min_r2": 0.01,
            "min_likelihood": 0.1, "min_velocity_genes": 30, "min_fit_genes": 30,
            "min_median_confidence": 0.5, "min_neighbor_cosine": 0.5,
            "max_abs_qc_rho": 0.5, "min_cross_model_cosine": 0.3,
            "excluded_hormone_genes": list(HORMONE_GENES)}


def model_gene_mask(index: pd.Index) -> np.ndarray:
    return np.asarray(~index.str.startswith(("mt-", "Rpl", "Rps")) & ~index.isin(HORMONE_GENES))


def row_cosine(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    num = np.einsum("ij,ij->i", first, second)
    den = np.linalg.norm(first, axis=1) * np.linalg.norm(second, axis=1)
    return np.divide(num, den, out=np.full(len(num), np.nan), where=den > 0)


def projection_eligible(summary: dict, sensitivity: list[dict]) -> bool:
    return bool(summary["dynamic_eligible_genes"] >= CRITERIA["min_fit_genes"]
                and summary["median_confidence"] >= CRITERIA["min_median_confidence"]
                and summary["median_dynamic_confidence"] >= CRITERIA["min_median_confidence"]
                and summary["median_cross_model_cosine"] >= CRITERIA["min_cross_model_cosine"]
                and len(sensitivity) == 2
                and all(np.isfinite(row[key]) and row[key] >= CRITERIA["min_neighbor_cosine"] for row in sensitivity
                        for key in ("median_cosine", "dynamic_median_cosine"))
                and all(row["dynamic_median_confidence"] >= CRITERIA["min_median_confidence"]
                        and row["cross_model_median_cosine"] >= CRITERIA["min_cross_model_cosine"]
                        and all(np.isfinite(row[key]) and row[key] < CRITERIA["max_abs_qc_rho"] for key in
                                ("abs_dynamic_length_total_counts_rho", "abs_dynamic_length_pct_counts_mt_rho"))
                        for row in sensitivity)
                and all(np.isfinite(summary[key]) and summary[key] < CRITERIA["max_abs_qc_rho"]
                        for key in ("abs_length_total_counts_rho", "abs_length_pct_counts_mt_rho",
                                    "abs_dynamic_length_total_counts_rho", "abs_dynamic_length_pct_counts_mt_rho")))


def strict_velocity_genes(data: ad.AnnData, key: str = "velocity") -> np.ndarray:
    """Undo scVelo's automatic r2-percentile relaxation on weak signal."""
    mask = (data.var[f"{key}_r2"].to_numpy() > CRITERIA["min_r2"])
    mask &= np.isfinite(data.var[f"{key}_gamma"]) & (data.var[f"{key}_gamma"] > 0)
    mask &= model_gene_mask(data.var_names)
    data.var[f"{key}_genes"] = mask
    return np.asarray(mask)


def _load_splice(loom: Path, beta: ad.AnnData, library: str) -> ad.AnnData:
    data = read_loom(loom, X_name="spliced", sparse=True, cleanup=False)
    if any(str(value).split(":", 1)[0] != library for value in data.obs_names):
        raise ValueError("Loom sample prefix does not match the requested library")
    data.obs_names = [loom_barcode(str(x), library) for x in data.obs_names]
    if not data.obs_names.is_unique:
        raise ValueError("Non-unique loom barcodes")
    target_obs = beta.obs_names[beta.obs.library_id.eq(library)]
    if not set(target_obs) <= set(data.obs_names):
        raise ValueError("Beta barcodes missing in splice matrix")
    if "Accession" not in data.var:
        raise ValueError("Loom missing Ensembl accession")
    data.var_names = data.var.Accession.astype(str)
    if not data.var_names.is_unique or not beta.var.gene_ids.is_unique:
        raise ValueError("Non-unique Ensembl mapping")
    gene_map = beta.var.reset_index(names="symbol").set_index("gene_ids")
    common = data.var_names.intersection(gene_map.index, sort=False)
    data = data[target_obs, common].copy()
    data.obs = beta.obs.loc[target_obs].copy()
    data.var["gene_ids"] = data.var_names
    data.var_names = gene_map.loc[common, "symbol"].to_numpy()
    data.obsm["X_umap"] = beta[target_obs].obsm["X_umap"].copy()
    for layer in ("spliced", "unspliced", "ambiguous"):
        matrix = sp.csr_matrix(data.layers[layer])
        if not np.all(np.isfinite(matrix.data)) or np.any(matrix.data < 0) or not np.all(matrix.data == np.floor(matrix.data)):
            raise ValueError(f"Invalid integer splice counts: {layer}")
        data.layers[layer] = matrix.astype(np.uint32)
    return data


def _save_plot(folder: Path, stem: str, figure, table: pd.DataFrame, title: str, interpretation: str) -> dict:
    save_figure(figure, folder / "figures" / stem, dpi=300)
    atomic_tsv(folder / "source_data" / f"{stem}.tsv", table)
    return {"title": title, "png": f"figures/{stem}.png", "pdf": f"figures/{stem}.pdf",
            "source": f"source_data/{stem}.tsv", "interpretation": interpretation}


def _prepare(data: ad.AnnData, n_neighbors: int) -> ad.AnnData:
    data = data[:, model_gene_mask(data.var_names)].copy()
    # scVelo 0.3.4 delegates HVG/log transforms to Scanpy.
    scv.pp.filter_and_normalize(data, min_shared_counts=CRITERIA["min_shared_counts"])
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=min(2000, data.n_vars), flavor="seurat", subset=True)
    sc.pp.pca(data, n_comps=min(40, data.n_vars - 1), svd_solver="arpack", random_state=42)
    sc.pp.neighbors(data, n_neighbors=n_neighbors, n_pcs=min(40, data.obsm["X_pca"].shape[1]), random_state=42)
    scv.pp.moments(data, n_pcs=None, n_neighbors=None)
    scv.tl.velocity(data, mode="stochastic", min_r2=CRITERIA["min_r2"])
    strict_velocity_genes(data)
    return data


@threadpool_limits.wrap(limits=1)
def model_velocity(beta_h5ad: Path, output_dir: Path) -> dict:
    # Dynamical fitting starts at most four workers, each with one BLAS thread.
    os.environ.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"})
    output_dir = output_dir.resolve()
    root = Path(__file__).resolve().parents[2]
    for name in ("tables", "figures", "source_data"):
        (output_dir / name).mkdir(parents=True, exist_ok=True)
    criteria_path = output_dir / "prespecified_screening_criteria.json"
    if criteria_path.exists() and json.loads(criteria_path.read_text()) != CRITERIA:
        raise ValueError("Prespecified velocity criteria changed; use a new run-id")
    atomic_json(criteria_path, CRITERIA)
    beta = ad.read_h5ad(beta_h5ad)
    libraries = sorted(beta.obs.library_id.astype(str).unique())
    looms = {lib: output_dir / "counts" / lib / f"{lib}.loom" for lib in libraries}
    frozen = _fingerprints([beta_h5ad, Path(__file__), criteria_path, root / "src/pfos_scrna/extension_velocity.py",
                            root / "uv.lock", root / "environments/velocity/uv.lock", *looms.values()])
    frozen["started_at"] = datetime.now().isoformat()
    atomic_json(output_dir / "execution_inputs.json", frozen)
    figures, summaries, tables, objects = [], [], [], []
    for library in libraries:
        data = _load_splice(looms[library], beta, library)
        cell_qc = pd.DataFrame(index=data.obs_names)
        cell_qc.index.name = "barcode"
        for layer in ("spliced", "unspliced", "ambiguous"):
            cell_qc[layer] = np.asarray(data.layers[layer].sum(axis=1)).ravel()
        cell_qc["unspliced_fraction"] = cell_qc.unspliced / cell_qc.sum(axis=1).clip(lower=1)
        cell_qc[["UMAP1", "UMAP2"]] = data.obsm["X_umap"]
        review_columns = [column for column in ("beta_leiden_0.5", "beta_low_identity_candidate",
                          "beta_low_complexity_candidate", "beta_residual_noncognate_candidate") if column in data.obs]
        cell_qc[review_columns] = data.obs[review_columns]
        gene_qc = pd.DataFrame(index=data.var_names)
        gene_qc.index.name = "gene"
        for layer in ("spliced", "unspliced", "ambiguous"):
            gene_qc[f"{layer}_counts"] = np.asarray(data.layers[layer].sum(axis=0)).ravel()
            gene_qc[f"{layer}_cells"] = np.asarray((data.layers[layer] > 0).sum(axis=0)).ravel()
        eligible = (gene_qc.spliced_cells >= CRITERIA["min_cells_each_layer"]) & (gene_qc.unspliced_cells >= CRITERIA["min_cells_each_layer"])
        gene_qc["excluded_from_model"] = ~model_gene_mask(gene_qc.index)
        gene_qc["hormone_audit_only"] = gene_qc.index.isin(HORMONE_GENES)
        eligible &= ~gene_qc.excluded_from_model
        gene_qc["coverage_eligible"] = eligible
        for stem, frame in ((f"{library}_cell_splice_qc", cell_qc), (f"{library}_gene_splice_qc", gene_qc)):
            atomic_tsv(output_dir / "tables" / f"{stem}.tsv", frame.reset_index())
            tables.append(f"tables/{stem}.tsv")
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        axes[0].scatter(np.log1p(cell_qc.spliced), np.log1p(cell_qc.unspliced), s=2, rasterized=True)
        axes[0].set(xlabel="log1p spliced UMIs", ylabel="log1p unspliced UMIs", title=library)
        axes[1].hist(cell_qc.unspliced_fraction, bins=40, color="#56B4E9")
        axes[1].set(xlabel="Unspliced fraction", ylabel="Captured beta cells")
        figures.append(_save_plot(output_dir, f"{library}_splice_coverage", fig, cell_qc.reset_index(),
                                 f"{library} β候选细胞剪接覆盖",
                                 f"左图每点为一个上游β标签细胞，共{data.n_obs:,}个；右图为未剪接比例分布，中位数{cell_qc.unspliced_fraction.median():.3f}。"
                                 f"{int(eligible.sum()):,}个基因满足两层各≥50细胞检出及非mt/ribosomal/预设六激素基因门。计数来自自有BAM，不等同SoupX校正矩阵；覆盖不是处理效应，低身份支持细胞仍保留审查标记。"))
        hormone_qc = gene_qc.loc[gene_qc.hormone_audit_only].reset_index()
        if not hormone_qc.empty:
            hormone_path = f"tables/{library}_excluded_hormone_splice_qc.tsv"
            atomic_tsv(output_dir / hormone_path, hormone_qc)
            tables.append(hormone_path)
            fig, ax = plt.subplots(figsize=(8, 4))
            positions = np.arange(len(hormone_qc))
            for j, layer in enumerate(("spliced", "unspliced", "ambiguous")):
                ax.bar(positions + (j-1)*.24, hormone_qc[f"{layer}_counts"], .24, label=layer)
            ax.set_xticks(positions, hormone_qc.gene)
            ax.set_yscale("symlog", linthresh=1)
            ax.set(ylabel="Summed BAM-derived UMIs (symlog)", xlabel="Genes retained for QC only")
            ax.legend(frameon=False)
            figures.append(_save_plot(output_dir, f"{library}_excluded_hormone_qc", fig, hormone_qc,
                                     f"{library} 六个激素基因的原始剪接计数审查",
                                     "三种颜色为spliced/unspliced/ambiguous层，纵轴为symlog总UMI。原始矩阵完整保留六基因，但因未做剪接层ambient校正，在查看真实矩阵前已将其排除主速度模型；不能把这些计数的大小解读为分泌速率或方向。"))
        data.write_h5ad(output_dir / f"{library}_splice_counts.h5ad", compression="gzip")
        objects.append(f"{library}_splice_counts.h5ad")
        summary = {"library_id": library, "beta_cells": data.n_obs, "mapped_genes": data.n_vars,
                   "coverage_genes": int(eligible.sum()), "status": "NO-GO_insufficient_splice_coverage",
                   "temporal_or_condition_inference": False}
        if eligible.sum() < CRITERIA["min_coverage_genes"]:
            summaries.append(summary)
            continue
        restricted = data[:, eligible.to_numpy()].copy()
        fit = _prepare(restricted, 30)
        mask = strict_velocity_genes(fit)
        summary["stochastic_velocity_genes"] = int(mask.sum())
        if mask.sum() < CRITERIA["min_velocity_genes"]:
            summary["status"] = "NO-GO_insufficient_unrelaxed_velocity_genes"
            atomic_tsv(output_dir / "tables" / f"{library}_stochastic_gene_fits.tsv", fit.var.reset_index(names="gene"))
            tables.append(f"tables/{library}_stochastic_gene_fits.tsv")
            summaries.append(summary)
            continue
        scv.tl.velocity_graph(fit, n_jobs=4, show_progress_bar=False)
        scv.tl.velocity_confidence(fit)
        sensitivity, neighbor_models = [], {}
        for neighbors in (15, 50):
            other = _prepare(restricted, neighbors)
            neighbor_models[neighbors] = other
            common = fit.var_names[mask].intersection(other.var_names[other.var.velocity_genes])
            if len(common) >= 10:
                cosine = row_cosine(fit[:, common].layers["velocity"], other[:, common].layers["velocity"])
                sensitivity.append({"neighbors": neighbors, "shared_velocity_genes": len(common),
                                    "median_cosine": float(np.nanmedian(cosine))})
            else:
                sensitivity.append({"neighbors": neighbors, "shared_velocity_genes": len(common), "median_cosine": np.nan})
        # Fit strongest supported stochastic genes; selection is blind to condition labels.
        dynamic_genes = fit.var.loc[mask].sort_values("velocity_r2", ascending=False).head(300).index.tolist()
        scv.tl.recover_dynamics(fit, var_names=dynamic_genes, n_jobs=4, max_iter=20, show_progress_bar=False)
        scv.tl.velocity(fit, mode="dynamical", vkey="velocity_dynamic", min_r2=CRITERIA["min_r2"], min_likelihood=CRITERIA["min_likelihood"])
        valid_fit = np.isfinite(fit.var.fit_likelihood) & (fit.var.fit_likelihood >= CRITERIA["min_likelihood"]) & (fit.var.fit_r2 > CRITERIA["min_r2"])
        valid_fit &= model_gene_mask(fit.var_names)
        fit.var["velocity_dynamic_genes"] = valid_fit
        summary["median_dynamic_confidence"] = np.nan
        summary["median_cross_model_cosine"] = np.nan
        for covariate in ("total_counts", "pct_counts_mt"):
            summary[f"abs_dynamic_length_{covariate}_rho"] = np.nan
        if valid_fit.sum() >= CRITERIA["min_fit_genes"]:
            scv.tl.velocity_graph(fit, vkey="velocity_dynamic", n_jobs=4, show_progress_bar=False)
            scv.tl.velocity_confidence(fit, vkey="velocity_dynamic")
            summary["median_dynamic_confidence"] = float(np.nanmedian(fit.obs.velocity_dynamic_confidence))
            for covariate in ("total_counts", "pct_counts_mt"):
                summary[f"abs_dynamic_length_{covariate}_rho"] = abs(float(spearmanr(fit.obs.velocity_dynamic_length, fit.obs[covariate]).statistic))
            common_mask = np.asarray(valid_fit) & mask
            if common_mask.sum() >= CRITERIA["min_fit_genes"]:
                summary["median_cross_model_cosine"] = float(np.nanmedian(row_cosine(
                    fit.layers["velocity"][:, common_mask], fit.layers["velocity_dynamic"][:, common_mask])))
        # Dynamic stability must be measured on refitted dynamic vectors, not borrowed
        # from the stochastic screening graph or its confidence values.
        for row in sensitivity:
            other = neighbor_models[row["neighbors"]]
            genes = [gene for gene in dynamic_genes if gene in other.var_names]
            scv.tl.recover_dynamics(other, var_names=genes, n_jobs=4, max_iter=20, show_progress_bar=False)
            scv.tl.velocity(other, mode="dynamical", vkey="velocity_dynamic", min_r2=CRITERIA["min_r2"], min_likelihood=CRITERIA["min_likelihood"])
            other_valid = (other.var.fit_likelihood >= CRITERIA["min_likelihood"]) & (other.var.fit_r2 > CRITERIA["min_r2"])
            other_valid &= model_gene_mask(other.var_names)
            other.var["velocity_dynamic_genes"] = other_valid
            row.update(dynamic_median_confidence=np.nan, cross_model_median_cosine=np.nan,
                       abs_dynamic_length_total_counts_rho=np.nan, abs_dynamic_length_pct_counts_mt_rho=np.nan)
            if other_valid.sum() >= CRITERIA["min_fit_genes"]:
                scv.tl.velocity_graph(other, vkey="velocity_dynamic", n_jobs=4, show_progress_bar=False)
                scv.tl.velocity_confidence(other, vkey="velocity_dynamic")
                row["dynamic_median_confidence"] = float(np.nanmedian(other.obs.velocity_dynamic_confidence))
                for covariate in ("total_counts", "pct_counts_mt"):
                    row[f"abs_dynamic_length_{covariate}_rho"] = abs(float(spearmanr(other.obs.velocity_dynamic_length, other.obs[covariate]).statistic))
                within_mask = np.asarray(other_valid & other.var.velocity_genes)
                if within_mask.sum() >= CRITERIA["min_fit_genes"]:
                    row["cross_model_median_cosine"] = float(np.nanmedian(row_cosine(
                        other.layers["velocity"][:, within_mask], other.layers["velocity_dynamic"][:, within_mask])))
            common = fit.var_names[valid_fit].intersection(other.var_names[other_valid])
            row["dynamic_shared_genes"] = len(common)
            row["dynamic_median_cosine"] = (float(np.nanmedian(row_cosine(
                fit[:, common].layers["velocity_dynamic"], other[:, common].layers["velocity_dynamic"])))
                if len(common) >= CRITERIA["min_fit_genes"] else np.nan)
            alt_path = f"tables/{library}_dynamic_fits_neighbors{row['neighbors']}.tsv"
            atomic_tsv(output_dir / alt_path, other.var.reset_index(names="gene"))
            tables.append(alt_path)
        atomic_tsv(output_dir / "tables" / f"{library}_neighbor_sensitivity.tsv", pd.DataFrame(sensitivity))
        tables.append(f"tables/{library}_neighbor_sensitivity.tsv")
        del neighbor_models
        summary.update(dynamic_fitted_genes=len(dynamic_genes), dynamic_eligible_genes=int(valid_fit.sum()),
                       median_confidence=float(np.nanmedian(fit.obs.velocity_confidence)),
                       min_neighbor_cosine=float(pd.DataFrame(sensitivity).median_cosine.min()),
                       min_dynamic_neighbor_cosine=float(pd.DataFrame(sensitivity).dynamic_median_cosine.min()))
        for covariate in ("total_counts", "pct_counts_mt"):
            summary[f"abs_length_{covariate}_rho"] = abs(float(spearmanr(fit.obs.velocity_length, fit.obs[covariate]).statistic))
        pass_gate = projection_eligible(summary, sensitivity)
        summary["status"] = "PASS_exploratory_local_projection" if pass_gate else "NO-GO_directional_projection_screen"
        fit_table = fit.var.reset_index(names="gene")
        atomic_tsv(output_dir / "tables" / f"{library}_gene_dynamics.tsv", fit_table)
        tables.append(f"tables/{library}_gene_dynamics.tsv")
        source = fit.obs[["velocity_length", "velocity_confidence", "total_counts", "pct_counts_mt"]].copy()
        source[review_columns] = fit.obs[review_columns]
        for column in ("velocity_dynamic_length", "velocity_dynamic_confidence"):
            if column in fit.obs:
                source[column] = fit.obs[column]
        source[["UMAP1", "UMAP2"]] = fit.obsm["X_umap"]
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        image = axes[0].scatter(source.UMAP1, source.UMAP2, c=source.velocity_confidence, s=2, cmap="viridis", rasterized=True)
        fig.colorbar(image, ax=axes[0], label="Neighbor velocity confidence")
        axes[0].set(xlabel="Beta atlas UMAP1", ylabel="Beta atlas UMAP2")
        axes[1].scatter(fit.var.velocity_r2, fit.var.fit_likelihood, s=9, color="#0072B2")
        axes[1].set(xlabel="Stochastic R²", ylabel="Dynamical likelihood")
        # One figure's source includes both cell and gene panels, explicit row_kind.
        source.index.name = "barcode"
        combined = pd.concat([source.reset_index().assign(row_kind="cell"), fit_table.assign(row_kind="gene")], ignore_index=True)
        figures.append(_save_plot(output_dir, f"{library}_dynamics_diagnostics", fig, combined,
                                 f"{library} 动力学拟合与置信度",
                                 f"左图颜色表示stochastic邻域速度一致性，中位数{summary['median_confidence']:.3f}；右图逐基因展示稳态R²与动力学似然。"
                                 f"实际拟合{len(dynamic_genes)}个基因，{int(valid_fit.sum())}个通过固定拟合门；dynamic自身置信度中位数{summary['median_dynamic_confidence']:.3f}。"
                                 f"预设诊断状态：{summary['status']}。置信度高不能单独证明真实命运方向。"))
        # Show QC associations and identity flags without post-hoc removal or new subtypes.
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        for ax, covariate, xlabel in zip(axes, ("total_counts", "pct_counts_mt"),
                                        ("log1p parent total counts", "Parent mt counts (%)")):
            x = np.log1p(source[covariate]) if covariate == "total_counts" else source[covariate]
            color = np.where(source["beta_low_identity_candidate"], "#E69F00", "#0072B2") if "beta_low_identity_candidate" in source else "#0072B2"
            ax.scatter(x, source.velocity_length, c=color, s=2, alpha=.4, rasterized=True)
            ax.set(xlabel=xlabel, ylabel="Stochastic velocity length")
        figures.append(_save_plot(output_dir, f"{library}_velocity_qc_associations", fig, source.reset_index(),
                                 f"{library} 速度与测序质量关联",
                                 f"横轴分别为上游总UMI及mt比例，纵轴为stochastic速度长度；|Spearman ρ|分别为"
                                 f"{summary['abs_length_total_counts_rho']:.3f}、{summary['abs_length_pct_counts_mt_rho']:.3f}。"
                                 "如存在身份复核标记，橙色为低β身份支持候选、蓝色为其余细胞；这是复核标签，不是新细胞亚型。dynamic模型的对应关联另列于资格表，不能借用本图通过门控。"))
        if "beta_leiden_0.5" in source:
            grouped = source.groupby("beta_leiden_0.5", observed=True)
            measurements = [column for column in source if column.startswith("velocity_") or column in ("total_counts", "pct_counts_mt")]
            group_qc = grouped[measurements].median().add_prefix("median_")
            group_qc.insert(0, "cells", grouped.size())
            for column in review_columns:
                if column.endswith("_candidate"):
                    group_qc[f"{column}_fraction"] = grouped[column].mean()
            group_qc.index.name = "beta_cluster"
            group_path = f"tables/{library}_cluster_velocity_qc.tsv"
            atomic_tsv(output_dir / group_path, group_qc.reset_index())
            tables.append(group_path)
        phase_genes = fit.var.loc[fit.var.fit_likelihood.notna()].sort_values("fit_likelihood", ascending=False).head(6).index
        fig, axes = plt.subplots(2, 3, figsize=(12, 7))
        phase_rows = []
        for ax, gene in zip(axes.flat, phase_genes):
            i = fit.var_names.get_loc(gene)
            s = np.asarray(fit.layers["Ms"][:, i]).ravel()
            u = np.asarray(fit.layers["Mu"][:, i]).ravel()
            ax.scatter(s, u, s=2, rasterized=True, alpha=0.4)
            ax.set(title=gene, xlabel="Spliced moment", ylabel="Unspliced moment")
            phase_rows.append(pd.DataFrame({"barcode": fit.obs_names, "gene": gene, "Ms": s, "Mu": u}))
        if phase_rows:
            for ax in list(axes.flat)[len(phase_rows):]:
                ax.set_visible(False)
            figures.append(_save_plot(output_dir, f"{library}_phase_portraits", fig, pd.concat(phase_rows, ignore_index=True),
                                     f"{library} 剪接相图", "各小图展示一个基因的剪接/未剪接邻域矩；按拟合似然选取作模型诊断，不能将此选取当独立验证或PFOS作用证据。"))
        else:
            plt.close(fig)
        if pass_gate:
            scv.tl.velocity_embedding(fit, basis="umap", vkey="velocity_dynamic")
            coords = fit.obsm["X_umap"]
            vectors = fit.obsm["velocity_dynamic_umap"]
            projected = pd.DataFrame({"barcode": fit.obs_names, "UMAP1": coords[:, 0], "UMAP2": coords[:, 1],
                                      "v1": vectors[:, 0], "v2": vectors[:, 1]})
            fig, ax = plt.subplots(figsize=(7, 6))
            ax.scatter(coords[:, 0], coords[:, 1], s=2, c="#bbbbbb", rasterized=True)
            selected = np.random.default_rng(42).choice(len(coords), min(800, len(coords)), replace=False)
            ax.quiver(coords[selected, 0], coords[selected, 1], vectors[selected, 0], vectors[selected, 1], alpha=0.5, width=0.002)
            ax.set(xlabel="Beta atlas UMAP1", ylabel="Beta atlas UMAP2", title=f"{library}: exploratory local projection")
            figures.append(_save_plot(output_dir, f"{library}_local_velocity", fig, projected,
                                     f"{library} β细胞内部速度投影", "仅在预设覆盖、拟合与敏感性门通过后绘制；箭头为表达空间估计的低维投影，不是成年发育轨迹或条件差异，不能据此证明真实时间/命运。"))
        fit.uns["velocity_interpretation_gate"] = summary
        fit.write_h5ad(output_dir / f"{library}_beta_velocity.h5ad", compression="gzip")
        objects.append(f"{library}_beta_velocity.h5ad")
        summaries.append(summary)
    atomic_tsv(output_dir / "tables/velocity_eligibility.tsv", pd.DataFrame(summaries))
    tables.append("tables/velocity_eligibility.tsv")
    atomic_text(output_dir / "methods.md", """# RNA velocity方法与边界

使用自有Cell Ranger BAM、同版GTF、GRCm39 repeat mask和uint32计数。按文库独立分析上游β标签候选，保留低复杂度/低身份支持等复核字段，不假设全部候选都是纯β细胞。velocyto计数环境独立uv锁定；scVelo在主环境运行。中断且未SHA封存的中间BAM不得在统一入口自动复用。

基因按≥50细胞同时检出spliced/unspliced预筛，剔除mt/ribosomal及Ins1/Ins2/Gcg/Sst/Ppy/Iapp六激素基因。六基因原始三层计数保留并另列QC；排除只针对主模型，不等同判定每条激素转录本为污染，也不能据其计数解释分泌方向。该规则在查看真实剪接矩阵前冻结，依据本项目激素校正审查和人胰岛原始研究提示的测量风险，属于项目适配而非通用velocity标准。

主邻居30，敏感性15/50；stochastic R²固定>0.01，禁用库自动降低阈值的结果。最多300个支持基因执行动力学拟合。动态模型在15/30/50邻居下分别重拟合，每种参数独立计算其confidence、总计数/mt关联及与对应stochastic向量的一致性，并比较不同邻居下dynamic向量的稳定性。动态投影不能借用stochastic诊断通过门控。筛查阈值详见prespecified_screening_criteria.json，是预设项目操作门，不是经过验证的通用判别标准；不输出条件P/FDR。

[scVelo官方文档](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html)；[velocyto官方CLI](https://velocyto.org/velocyto.py/tutorial/cli.html)；[人胰岛ambient RNA原始研究](https://link.springer.com/article/10.1038/s44318-026-00744-w)。
""")
    _verify_fingerprints(frozen)
    frozen["verification_status"] = "PASS"
    atomic_json(output_dir / "execution_inputs.json", frozen)
    tables.extend(["methods.md", "prespecified_screening_criteria.json", "execution_inputs.json"])
    if (output_dir / "counts/counting_status.tsv").exists():
        tables.append("counts/counting_status.tsv")
    manifest = {"status": "completed", "method": "自有BAM剪接计数与逐库β细胞scVelo覆盖/拟合/邻居敏感性诊断",
                "metrics": {"libraries": summaries}, "figures": figures, "tables": tables, "objects": objects,
                "condition_inference": "NO-GO_n1_per_condition"}
    atomic_json(output_dir / "report_manifest.json", manifest)
    return manifest

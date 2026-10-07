#!/usr/bin/env python3
"""Tissue-agnostic spatial transcriptomics pipeline.

The module keeps imports lazy so validation and help work before scientific
packages are installed. Heavy optional methods fail with explicit messages.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable


STEPS = [
    "01_validate_inputs",
    "02_build_objects",
    "03_qc_filter",
    "04_normalize_features",
    "05_integrate_cluster",
    "06_spatial_neighbors_stats",
    "07_spatial_domains",
    "08_cell_type_annotation",
    "09_deconvolution",
    "10_segment_objects",
    "11_programs_descriptive",
    "12_neighborhood_communication",
    "13_multi_sample_summary",
    "14_report_index",
]


DEFAULT_CONFIG: dict[str, Any] = {
    "project": {
        "sample_key": "sample",
        "group_key": "group",
        "batch_key": "batch",
        "spatial_key": "spatial",
        "random_seed": 20260629,
        "organism": "mouse",
    },
    "paths": {
        "processed_dir": "02processed_h5ad",
        "integrated_dir": "03integrated_h5ad",
        "deconvolution_dir": "04deconvolution",
        "result_dir": "result",
        "cache_dir": ".cache/spatial",
        "marker_gene_file": "config/marker_genes.tsv",
    },
    "input": {
        "visium_count_file": "filtered_feature_bc_matrix.h5",
        "load_images": True,
        "matrix_orientation": "genes_by_spots",
        "make_var_names_unique": True,
    },
    "qc": {
        "min_counts": 500,
        "min_genes": 200,
        "max_mito_percent": 25,
        "min_cells_per_gene": 3,
        "spatial_local_outlier_k": 8,
        "spatial_low_count_ratio": 0.25,
        "mito_prefix": {"mouse": "mt-", "human": "MT-"},
        "ribo_prefix": {"mouse": ["Rpl", "Rps"], "human": ["RPL", "RPS"]},
    },
    "normalization": {
        "target_sum": 10000,
        "hvg_flavor": "cell_ranger",
        "n_top_genes": 3000,
        "regress_out": ["total_counts", "pct_counts_mt"],
        "scale_max_value": 10,
        "n_pcs": 50,
    },
    "integration": {
        "method": "none",
        "n_neighbors": 15,
        "leiden_resolution": 0.8,
        "umap_min_dist": 0.3,
    },
    "spatial": {
        "coord_type": "generic",
        "n_neighs": 6,
        "radius": None,
        "moran_top_genes": 1000,
        "moran_permutations": 100,
        "nhood_permutations": 1000,
    },
    "domains": {
        "method": "expression_spatial_leiden",
        "resolution": 0.8,
        "spatial_weight": 0.25,
        "min_segment_spots": 5,
    },
    "annotation": {
        "method": "marker_score",
        "min_marker_genes": 2,
        "unknown_margin": 0.05,
    },
    "deconvolution": {
        "method": "marker_score",
        "reference_h5ad": "",
        "reference_cell_type_key": "cell_type",
        "reference_batch_key": "batch",
        "cell2location": {
            "n_cells_per_location": 8,
            "detection_alpha": 20,
            "max_epochs_reference": 250,
            "max_epochs_spatial": 30000,
            "use_gpu": "auto",
        },
    },
    "segments": {
        "source_key": "spatial_domain",
        "layer_key": "cell_type_pred",
        "distance_threshold": "auto",
        "min_spots": 5,
        "compute_scaled_distance": True,
    },
    "programs": {
        "method": "nmf",
        "n_components": 12,
        "top_genes_per_program": 50,
    },
    "communication": {
        "cluster_key": "cell_type_pred",
        "min_spots_per_group": 10,
        "run_squidpy_ligrec": False,
    },
    "plotting": {
        "palette": ["#A1C9F4", "#FFB482", "#8DE5A1", "#FF9F9B", "#B39FDB", "#FDFFB6"],
        "dpi": 300,
        "save_pdf": True,
        "save_png": True,
    },
}


@dataclass
class Context:
    project_root: Path
    config: dict[str, Any]
    manifest: list[dict[str, str]]

    def path(self, *parts: str | Path) -> Path:
        return self.project_root.joinpath(*map(str, parts))

    @property
    def result_dir(self) -> Path:
        return self.path(self.config["paths"]["result_dir"])

    def step_dir(self, step: str) -> Path:
        return self.result_dir / step

    def processed_dir(self) -> Path:
        return self.path(self.config["paths"]["processed_dir"])

    def integrated_dir(self) -> Path:
        return self.path(self.config["paths"]["integrated_dir"])

    def deconv_dir(self) -> Path:
        return self.path(self.config["paths"]["deconvolution_dir"])


def deep_update(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    out = json.loads(json.dumps(base))
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_update(out[key], value)
        else:
            out[key] = value
    return out


def load_config(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"MISSING config: {path}")
    try:
        import yaml
    except ImportError as exc:
        raise SystemExit("PyYAML is required to read spatial_config.yaml") from exc
    with path.open("r", encoding="utf-8") as handle:
        user_config = yaml.safe_load(handle)
    if not isinstance(user_config, dict) or not user_config:
        raise ValueError("配置必须是非空映射")
    config = deep_update(DEFAULT_CONFIG, user_config)
    for section, methods in {"annotation": {"marker_score"}, "programs": {"nmf"}, "domains": {"expression_spatial_leiden"}, "integration": {"none", "harmony"}, "deconvolution": {"marker_score", "cell2location"}}.items():
        if config[section]["method"] not in methods:
            raise ValueError(f"{section}.method 未实现: {config[section]['method']}")
    return config


def load_manifest(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"Manifest not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [{k: (v or "").strip() for k, v in row.items()} for row in reader]
    if not rows:
        raise SystemExit(f"Manifest is empty: {path}")
    return rows


def active_rows(ctx: Context) -> list[dict[str, str]]:
    rows = []
    for row in ctx.manifest:
        if row.get("exclude", "0").strip().lower() in {"1", "true", "yes"}:
            continue
        rows.append(row)
    if not rows:
        raise SystemExit("No active samples in manifest")
    return rows


def require_pkg(name: str, pip_name: str | None = None):
    try:
        return __import__(name)
    except ImportError as exc:
        pkg = pip_name or name
        raise SystemExit(f"Required package '{pkg}' is not installed for this step") from exc


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, obj: Any) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    ensure_dir(path.parent)
    if fieldnames is None:
        keys: list[str] = []
        for row in rows:
            for key in row:
                if key not in keys:
                    keys.append(key)
        fieldnames = keys
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def savefig_versions(fig, plot_base: Path, ctx: Context) -> None:
    ensure_dir(plot_base.parent)
    plotting = ctx.config["plotting"]
    if plotting.get("save_pdf", True):
        fig.savefig(plot_base.with_suffix(".pdf"), bbox_inches="tight")
    if plotting.get("save_png", True):
        fig.savefig(plot_base.with_suffix(".png"), dpi=plotting.get("dpi", 300), bbox_inches="tight")


def sparse_to_dense(x):
    import scipy.sparse as sp

    return x.toarray() if sp.issparse(x) else x


def read_table(path: Path):
    import pandas as pd

    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, index_col=0)


def resolve_input_path(ctx: Context, row: dict[str, str], key: str = "input_path") -> Path | None:
    value = row.get(key, "")
    if not value:
        return None
    path = Path(value)
    return path if path.is_absolute() else ctx.project_root / path


def latest_existing(paths: list[Path]) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise SystemExit("None of these expected files exist: " + ", ".join(map(str, paths)))


def read_adata(path: Path):
    ad = require_pkg("anndata")
    if not path.exists():
        raise SystemExit(f"AnnData file not found: {path}")
    return ad.read_h5ad(path)


def write_adata(adata, path: Path) -> None:
    ensure_dir(path.parent)
    adata.write_h5ad(path, compression="gzip")


def sample_key(ctx: Context) -> str:
    return ctx.config["project"]["sample_key"]


def group_key(ctx: Context) -> str:
    return ctx.config["project"]["group_key"]


def batch_key(ctx: Context) -> str:
    return ctx.config["project"]["batch_key"]


def spatial_key(ctx: Context) -> str:
    return ctx.config["project"]["spatial_key"]


def step_01_validate_inputs(ctx: Context, step: str, *, write: bool = True) -> None:
    rows_out: list[dict[str, Any]] = []
    seen: set[str] = set()
    valid_types = {"visium_dir", "h5ad", "matrix_coords"}
    for row in active_rows(ctx):
        sample = row.get("sample", "")
        input_type = row.get("input_type", "")
        status = "ok"
        messages: list[str] = []
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
            status = "error"
            messages.append("missing sample")
        if sample in seen:
            status = "error"
            messages.append("duplicate sample")
        seen.add(sample)
        if not row.get(group_key(ctx), "").strip():
            status = "error"
            messages.append("missing explicit group")
        if row.get("organism") and row["organism"] != ctx.config["project"]["organism"]:
            status = "error"
            messages.append("manifest organism differs from project")
        if input_type not in valid_types:
            status = "error"
            messages.append(f"input_type must be one of {sorted(valid_types)}")
        if input_type in {"visium_dir", "h5ad"}:
            input_path = resolve_input_path(ctx, row)
            if input_path is None or not input_path.exists():
                status = "error"
                messages.append(f"input_path not found: {input_path}")
        if input_type == "h5ad" and input_path is not None and not input_path.is_file():
            status = "error"
            messages.append("h5ad input_path must be a file")
        if input_type == "visium_dir" and input_path is not None:
            count = input_path / (row.get("count_file") or ctx.config["input"]["visium_count_file"])
            if not input_path.is_dir() or not count.is_file():
                status = "error"
                messages.append("missing Visium directory/count file")
        if input_type == "matrix_coords":
            for col in ["count_file", "spatial_file"]:
                p = resolve_input_path(ctx, row, col)
                if p is None or not p.is_file():
                    status = "error"
                    messages.append(f"{col} not found: {p}")
        row_out = dict(row)
        row_out["resolved_input_path"] = str(resolve_input_path(ctx, row) or "")
        row_out["status"] = status
        row_out["message"] = "; ".join(messages)
        rows_out.append(row_out)

    outdir = ctx.step_dir(step)
    if not rows_out:
        raise ValueError("No active samples")
    if write:
        write_tsv(outdir / "tables" / "spatial_input_status.tsv", rows_out)
        write_json(outdir / "tables" / "pipeline_config_resolved.json", ctx.config)
    if any(row["status"] == "error" for row in rows_out):
        raise SystemExit("Input validation failed: " + json.dumps(rows_out, ensure_ascii=False))


def _add_manifest_metadata(adata, row: dict[str, str]) -> None:
    for key, value in row.items():
        if key in {"input_path", "count_file", "spatial_file", "image_file", "reference_h5ad", "notes"}:
            continue
        if key and value:
            adata.obs[key] = value


def _read_matrix_coords(ctx: Context, row: dict[str, str]):
    import anndata as ad
    import pandas as pd

    count_file = resolve_input_path(ctx, row, "count_file")
    spatial_file = resolve_input_path(ctx, row, "spatial_file")
    counts = read_table(count_file)
    orientation = ctx.config["input"]["matrix_orientation"]
    if orientation == "spots_by_genes":
        obs = counts.index.astype(str)
        var = counts.columns.astype(str)
        matrix = counts.values
    else:
        obs = counts.columns.astype(str)
        var = counts.index.astype(str)
        matrix = counts.T.values
    adata = ad.AnnData(matrix)
    adata.obs_names = obs
    adata.var_names = var
    spatial_df = pd.read_csv(spatial_file, sep=None, engine="python")
    barcode_col = "barcode" if "barcode" in spatial_df.columns else spatial_df.columns[0]
    x_col = "x" if "x" in spatial_df.columns else spatial_df.columns[1]
    y_col = "y" if "y" in spatial_df.columns else spatial_df.columns[2]
    spatial_df = spatial_df.set_index(barcode_col).loc[adata.obs_names]
    adata.obsm[spatial_key(ctx)] = spatial_df[[x_col, y_col]].to_numpy()
    return adata


def step_02_build_objects(ctx: Context, step: str) -> None:
    ad = require_pkg("anndata")
    import numpy as np

    ensure_dir(ctx.processed_dir())
    built_rows: list[dict[str, Any]] = []
    adatas = []
    for row in active_rows(ctx):
        sample = row["sample"]
        input_type = row["input_type"]
        if input_type == "visium_dir":
            sc = require_pkg("scanpy")
            kwargs = {
                "path": str(resolve_input_path(ctx, row) or ""),
                "count_file": row.get("count_file") or ctx.config["input"]["visium_count_file"],
                "library_id": row.get("library_id") or sample,
                "load_images": bool(ctx.config["input"].get("load_images", True)),
            }
            adata = sc.read_visium(**kwargs)
        elif input_type == "h5ad":
            adata = ad.read_h5ad(resolve_input_path(ctx, row))
        elif input_type == "matrix_coords":
            adata = _read_matrix_coords(ctx, row)
        else:
            raise SystemExit(f"Unsupported input_type: {input_type}")

        matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
        values = matrix.data if hasattr(matrix, "tocsr") else np.asarray(matrix)
        if not np.all(np.isfinite(values)) or np.any(values < 0) or not np.allclose(values, np.rint(values)):
            raise ValueError(f"{sample}: raw counts must be nonnegative finite integers")
        if not adata.obs_names.is_unique:
            raise ValueError(f"{sample}: duplicate spot identifiers")
        coords = adata.obsm.get(spatial_key(ctx))
        if coords is None or np.asarray(coords).shape != (adata.n_obs,2) or not np.all(np.isfinite(coords)):
            raise ValueError(f"{sample}: finite n_spots × 2 spatial coordinates required")
        adata.X = matrix.copy()
        adata.layers["counts"] = matrix.copy()
        if ctx.config["input"].get("make_var_names_unique", True):
            adata.var_names_make_unique()
        _add_manifest_metadata(adata, row)
        adata.obs[sample_key(ctx)] = sample
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()
        adata.obs_names = [f"{sample}:{x}" for x in adata.obs_names]
        out_file = ctx.processed_dir() / f"{sample}.raw.h5ad"
        write_adata(adata, out_file)
        adatas.append(adata)
        built_rows.append(
            {
                "sample": sample,
                "n_spots": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "has_spatial": spatial_key(ctx) in adata.obsm,
                "output_h5ad": str(out_file),
            }
        )

    combined = ad.concat(adatas, join="outer", merge="same", index_unique=None, fill_value=0)
    combined.var_names_make_unique()
    if "counts" not in combined.layers:
        combined.layers["counts"] = combined.X.copy()
    write_adata(combined, ctx.processed_dir() / "all_samples.raw.h5ad")
    write_tsv(ctx.step_dir(step) / "tables" / "built_objects.tsv", built_rows)


def step_03_qc_filter(ctx: Context, step: str) -> None:
    sc = require_pkg("scanpy")
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from sklearn.neighbors import NearestNeighbors

    adata = read_adata(ctx.processed_dir() / "all_samples.raw.h5ad")
    organism = ctx.config["project"].get("organism", "mouse")
    qc_conf = ctx.config["qc"]
    mito_prefix = qc_conf["mito_prefix"].get(organism, "MT-")
    ribo_prefix = tuple(qc_conf["ribo_prefix"].get(organism, ["RPL", "RPS"]))
    gene_names = adata.var_names.astype(str)
    adata.var["mt"] = [g.startswith(mito_prefix) for g in gene_names]
    adata.var["ribo"] = [g.startswith(ribo_prefix) for g in gene_names]
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=None)

    keep = (
        (adata.obs["total_counts"] >= qc_conf["min_counts"])
        & (adata.obs["n_genes_by_counts"] >= qc_conf["min_genes"])
        & (adata.obs["pct_counts_mt"] <= qc_conf["max_mito_percent"])
    )
    adata.obs["qc_pass_basic"] = keep
    adata.obs["spatial_low_count_outlier"] = False
    if spatial_key(ctx) in adata.obsm:
        coords = adata.obsm[spatial_key(ctx)]
        totals = adata.obs["total_counts"].to_numpy()
        for sample in adata.obs[sample_key(ctx)].astype(str).unique():
            idx = np.where(adata.obs[sample_key(ctx)].astype(str).to_numpy() == sample)[0]
            if len(idx) <= qc_conf["spatial_local_outlier_k"]:
                continue
            nn = NearestNeighbors(n_neighbors=qc_conf["spatial_local_outlier_k"] + 1).fit(coords[idx])
            neigh = nn.kneighbors(coords[idx], return_distance=False)[:, 1:]
            local_median = np.median(totals[idx][neigh], axis=1)
            flags = totals[idx] < (local_median * qc_conf["spatial_low_count_ratio"])
            adata.obs.iloc[idx, adata.obs.columns.get_loc("spatial_low_count_outlier")] = flags
    adata.obs["qc_pass"] = adata.obs["qc_pass_basic"] & ~adata.obs["spatial_low_count_outlier"]
    adata = adata[adata.obs["qc_pass"].to_numpy(), :].copy()
    sc.pp.filter_genes(adata, min_cells=qc_conf["min_cells_per_gene"])
    write_adata(adata, ctx.processed_dir() / "all_samples.qc.h5ad")

    summary = (
        adata.obs.groupby(sample_key(ctx), observed=True)
        .agg(
            spots=("total_counts", "size"),
            median_counts=("total_counts", "median"),
            median_genes=("n_genes_by_counts", "median"),
            median_mito_pct=("pct_counts_mt", "median"),
        )
        .reset_index()
    )
    summary.to_csv(ctx.step_dir(step) / "tables" / "sample_qc_summary.tsv", sep="\t", index=False)
    obs_cols = [sample_key(ctx), group_key(ctx), "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "spatial_low_count_outlier"]
    obs_cols = [c for c in obs_cols if c in adata.obs.columns]
    adata.obs[obs_cols].to_csv(ctx.step_dir(step) / "tables" / "spot_qc_metrics.tsv", sep="\t")

    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, axes = plt.subplots(1, 3, figsize=(12, 3))
    for ax, col in zip(axes, ["total_counts", "n_genes_by_counts", "pct_counts_mt"]):
        sns.boxplot(data=adata.obs, x=sample_key(ctx), y=col, ax=ax, color="#A1C9F4")
        ax.tick_params(axis="x", rotation=45)
    savefig_versions(fig, ctx.step_dir(step) / "plots" / "sample_qc_boxplots", ctx)
    plt.close(fig)


def step_04_normalize_features(ctx: Context, step: str) -> None:
    sc = require_pkg("scanpy")

    adata = read_adata(ctx.processed_dir() / "all_samples.qc.h5ad")
    norm = ctx.config["normalization"]
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=norm["target_sum"])
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor=norm["hvg_flavor"], n_top_genes=norm["n_top_genes"],
                                batch_key=sample_key(ctx) if sample_key(ctx) in adata.obs else None)
    regressors = [c for c in norm.get("regress_out", []) if c in adata.obs.columns]
    if regressors:
        sc.pp.regress_out(adata, regressors)
    sc.pp.scale(adata, max_value=norm["scale_max_value"])
    pca_features = int(adata.var["highly_variable"].sum())
    n_comps = min(norm["n_pcs"], adata.n_obs - 1, pca_features - 1)
    if n_comps < 1:
        raise ValueError("PCA requires at least 2 spots and 2 highly variable genes")
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
    write_adata(adata, ctx.processed_dir() / "all_samples.normalized.h5ad")
    hvg = adata.var.reset_index().rename(columns={"index": "gene"})
    hvg.to_csv(ctx.step_dir(step) / "tables" / "highly_variable_genes.tsv", sep="\t", index=False)


def step_05_integrate_cluster(ctx: Context, step: str) -> None:
    sc = require_pkg("scanpy")
    import pandas as pd

    ensure_dir(ctx.integrated_dir())
    adata = read_adata(ctx.processed_dir() / "all_samples.normalized.h5ad")
    integ = ctx.config["integration"]
    use_rep = "X_pca"
    if integ["method"] == "harmony":
        if batch_key(ctx) not in adata.obs:
            raise ValueError("Harmony requires explicit batch metadata")
        import scanpy.external as sce
        sce.pp.harmony_integrate(adata, key=batch_key(ctx), basis="X_pca", adjusted_basis="X_pca_harmony")
        use_rep = "X_pca_harmony"
    elif integ["method"] != "none":
        raise ValueError("integration.method must be none or harmony")
    sc.pp.neighbors(adata, n_neighbors=integ["n_neighbors"], use_rep=use_rep)
    sc.tl.umap(adata, min_dist=integ["umap_min_dist"])
    sc.tl.leiden(adata, resolution=integ["leiden_resolution"], key_added="leiden")
    write_adata(adata, ctx.integrated_dir() / "all_samples.clustered.h5ad")
    counts = adata.obs.groupby([sample_key(ctx), "leiden"], observed=True).size().reset_index(name="n_spots")
    counts.to_csv(ctx.step_dir(step) / "tables" / "leiden_counts_by_sample.tsv", sep="\t", index=False)

    import matplotlib.pyplot as plt

    for color in ["leiden", group_key(ctx), sample_key(ctx)]:
        if color not in adata.obs:
            continue
        sc.pl.umap(adata, color=color, show=False)
        fig = plt.gcf()
        savefig_versions(fig, ctx.step_dir(step) / "plots" / f"umap_{color}", ctx)
        plt.close(fig)


def step_06_spatial_neighbors_stats(ctx: Context, step: str) -> None:
    sc = require_pkg("scanpy")
    sq = require_pkg("squidpy")
    import numpy as np
    import pandas as pd

    adata = read_adata(ctx.integrated_dir() / "all_samples.clustered.h5ad")
    sp_conf = ctx.config["spatial"]
    kwargs: dict[str, Any] = {"coord_type": sp_conf.get("coord_type", "generic")}
    if sp_conf.get("radius"):
        kwargs["radius"] = sp_conf["radius"]
    else:
        kwargs["n_neighs"] = sp_conf["n_neighs"]
    sq.gr.spatial_neighbors(adata, spatial_key=spatial_key(ctx), library_key=sample_key(ctx), **kwargs)
    cluster_key = "leiden"
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, n_perms=sp_conf["nhood_permutations"])
    if "highly_variable" in adata.var:
        hv_mask = adata.var["highly_variable"].fillna(False).to_numpy()
    else:
        hv_mask = np.ones(adata.n_vars, dtype=bool)
    genes = adata.var_names[hv_mask][: sp_conf["moran_top_genes"]]
    if len(genes) == 0:
        genes = adata.var_names[: min(sp_conf["moran_top_genes"], adata.n_vars)]
    sq.gr.spatial_autocorr(adata, mode="moran", genes=list(genes), n_perms=sp_conf["moran_permutations"])
    if "moranI" in adata.uns:
        adata.uns["moranI"].to_csv(ctx.step_dir(step) / "tables" / "spatial_autocorr_moranI.tsv", sep="\t")
    key = f"{cluster_key}_nhood_enrichment"
    if key in adata.uns:
        z = pd.DataFrame(
            adata.uns[key]["zscore"],
            index=adata.obs[cluster_key].cat.categories,
            columns=adata.obs[cluster_key].cat.categories,
        )
        z.to_csv(ctx.step_dir(step) / "tables" / "neighborhood_enrichment_zscore.tsv", sep="\t")
    write_adata(adata, ctx.integrated_dir() / "all_samples.spatial_graph.h5ad")


def step_07_spatial_domains(ctx: Context, step: str) -> None:
    sc = require_pkg("scanpy")
    import numpy as np
    import pandas as pd
    from sklearn.preprocessing import StandardScaler

    adata = read_adata(latest_existing([ctx.integrated_dir() / "all_samples.spatial_graph.h5ad", ctx.integrated_dir() / "all_samples.clustered.h5ad"]))
    if "X_pca_harmony" in adata.obsm:
        expr = adata.obsm["X_pca_harmony"]
    else:
        expr = adata.obsm["X_pca"]
    coords = adata.obsm.get(spatial_key(ctx))
    if coords is None:
        raise SystemExit(f"Missing obsm['{spatial_key(ctx)}'] for spatial domain analysis")
    x_expr = StandardScaler().fit_transform(expr[:, : min(30, expr.shape[1])])
    x_spatial = StandardScaler().fit_transform(coords) * float(ctx.config["domains"]["spatial_weight"])
    adata.obsm["X_expr_spatial"] = np.hstack([x_expr, x_spatial])
    sc.pp.neighbors(adata, use_rep="X_expr_spatial", n_neighbors=ctx.config["integration"]["n_neighbors"], key_added="domain_neighbors")
    sc.tl.leiden(adata, resolution=ctx.config["domains"]["resolution"], neighbors_key="domain_neighbors", key_added="spatial_domain")
    counts = adata.obs.groupby([sample_key(ctx), "spatial_domain"], observed=True).size().reset_index(name="n_spots")
    counts.to_csv(ctx.step_dir(step) / "tables" / "spatial_domain_counts.tsv", sep="\t", index=False)
    write_adata(adata, ctx.integrated_dir() / "all_samples.domains.h5ad")

    import matplotlib.pyplot as plt

    sc.pl.umap(adata, color="spatial_domain", show=False)
    savefig_versions(plt.gcf(), ctx.step_dir(step) / "plots" / "umap_spatial_domain", ctx)
    plt.close()


def load_marker_table(ctx: Context):
    import pandas as pd

    path = ctx.path(ctx.config["paths"]["marker_gene_file"])
    if not path.exists():
        raise SystemExit(f"Marker gene file not found: {path}")
    return pd.read_csv(path, sep="\t")


def step_08_cell_type_annotation(ctx: Context, step: str) -> None:
    sc = require_pkg("scanpy")
    import numpy as np
    import pandas as pd

    adata = read_adata(ctx.integrated_dir() / "all_samples.domains.h5ad")
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=ctx.config["normalization"]["target_sum"])
    adata.uns.pop("log1p", None)
    sc.pp.log1p(adata)
    markers = load_marker_table(ctx)
    tissue = ctx.config["project"].get("default_tissue", "generic")
    markers = markers[(markers["tissue"].isin(["generic", tissue]))]
    score_cols: list[str] = []
    for cell_type, sub in markers.groupby("cell_type"):
        genes = [g for g in sub["gene"].astype(str).unique() if g in adata.var_names]
        if len(genes) < ctx.config["annotation"]["min_marker_genes"]:
            continue
        col = f"score_{cell_type}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=col, use_raw=False)
        score_cols.append(col)
    if not score_cols:
        adata.obs["cell_type_pred"] = "unknown"
    else:
        scores = adata.obs[score_cols]
        best = scores.idxmax(axis=1).str.replace("^score_", "", regex=True)
        margin = scores.max(axis=1) - scores.apply(lambda row: row.nlargest(2).iloc[-1] if len(row) > 1 else 0, axis=1)
        best = best.where(margin >= ctx.config["annotation"]["unknown_margin"], "unknown")
        adata.obs["cell_type_pred"] = pd.Categorical(best)
        scores.to_csv(ctx.step_dir(step) / "tables" / "marker_scores_by_spot.tsv", sep="\t")
    adata.obs.groupby([sample_key(ctx), "cell_type_pred"], observed=True).size().reset_index(name="n_spots").to_csv(
        ctx.step_dir(step) / "tables" / "cell_type_counts.tsv", sep="\t", index=False
    )
    write_adata(adata, ctx.integrated_dir() / "all_samples.annotated.h5ad")


def _export_marker_scores(ctx: Context, adata, outdir: Path) -> None:
    import pandas as pd

    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if not score_cols:
        pd.DataFrame(index=adata.obs_names).to_csv(outdir / "tables" / "marker_scores_by_spot.tsv", sep="\t")
        return
    abundance = adata.obs[[sample_key(ctx)] + score_cols].copy()
    abundance.columns = [sample_key(ctx)] + [c.replace("score_", "") for c in score_cols]
    abundance.to_csv(outdir / "tables" / "marker_scores_by_spot.tsv", sep="\t")


def _run_cell2location(ctx: Context, adata, outdir: Path) -> None:
    sc = require_pkg("scanpy")
    c2l = require_pkg("cell2location")
    from cell2location.models import Cell2location, RegressionModel

    conf = ctx.config["deconvolution"]
    c2l_conf = conf["cell2location"]
    ref_path = conf.get("reference_h5ad") or ""
    if not ref_path:
        raise SystemExit("deconvolution.reference_h5ad is required for Cell2location")
    ref_path = Path(ref_path)
    if not ref_path.is_absolute():
        ref_path = ctx.project_root / ref_path
    if not ref_path.exists():
        raise SystemExit(f"Cell2location reference_h5ad not found: {ref_path}")

    ref = sc.read_h5ad(ref_path)
    cell_type_key = conf["reference_cell_type_key"]
    ref_batch_key = conf.get("reference_batch_key", "")
    if cell_type_key not in ref.obs:
        raise SystemExit(f"Reference cell type key not found in obs: {cell_type_key}")
    if "counts" in ref.layers:
        ref.X = ref.layers["counts"].copy()
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    shared = ref.var_names.intersection(adata.var_names)
    if len(shared) < 500:
        raise SystemExit(f"Too few shared genes for Cell2location: {len(shared)}")
    ref = ref[:, shared].copy()
    adata = adata[:, shared].copy()

    setup_kwargs: dict[str, Any] = {"labels_key": cell_type_key}
    if ref_batch_key and ref_batch_key in ref.obs:
        setup_kwargs["batch_key"] = ref_batch_key
    RegressionModel.setup_anndata(ref, **setup_kwargs)
    ref_model = RegressionModel(ref)
    use_gpu = c2l_conf.get("use_gpu", "auto")
    ref_model.train(max_epochs=int(c2l_conf["max_epochs_reference"]), use_gpu=use_gpu)
    ref = ref_model.export_posterior(ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500, "use_gpu": use_gpu})

    factor_names = ref.uns["mod"]["factor_names"]
    cols = [f"means_per_cluster_mu_fg_{name}" for name in factor_names]
    inf_aver = ref.varm["means_per_cluster_mu_fg"][cols].copy()
    inf_aver.columns = factor_names

    spatial_batch = sample_key(ctx) if sample_key(ctx) in adata.obs else None
    if spatial_batch:
        Cell2location.setup_anndata(adata=adata, batch_key=spatial_batch)
    else:
        Cell2location.setup_anndata(adata=adata)
    model = Cell2location(
        adata,
        cell_state_df=inf_aver,
        N_cells_per_location=int(c2l_conf["n_cells_per_location"]),
        detection_alpha=float(c2l_conf["detection_alpha"]),
    )
    model.train(max_epochs=int(c2l_conf["max_epochs_spatial"]), batch_size=None, train_size=1, use_gpu=use_gpu)
    adata = model.export_posterior(
        adata,
        sample_kwargs={"num_samples": 1000, "batch_size": adata.n_obs, "use_gpu": use_gpu},
    )
    abundance_key = "q05_cell_abundance_w_sf"
    if abundance_key in adata.obsm:
        abundance = adata.obsm[abundance_key]
        abundance.to_csv(outdir / "tables" / "cell2location_q05_cell_abundance.tsv", sep="\t")
        for col in abundance.columns:
            adata.obs[f"c2l_{col}"] = abundance[col].to_numpy()
    ensure_dir(ctx.deconv_dir() / "cell2location_model")
    model.save(str(ctx.deconv_dir() / "cell2location_model"), overwrite=True)
    write_adata(adata, ctx.deconv_dir() / "all_samples.cell2location.h5ad")


def step_09_deconvolution(ctx: Context, step: str) -> None:
    ensure_dir(ctx.deconv_dir())
    adata = read_adata(ctx.integrated_dir() / "all_samples.annotated.h5ad")
    method = ctx.config["deconvolution"].get("method", "marker_score")
    if method == "marker_score":
        _export_marker_scores(ctx, adata, ctx.step_dir(step))
    elif method == "cell2location":
        _run_cell2location(ctx, adata, ctx.step_dir(step))
        adata = read_adata(ctx.deconv_dir() / "all_samples.cell2location.h5ad")
    else:
        raise SystemExit(f"Unsupported deconvolution.method: {method}")
    write_json(ctx.step_dir(step) / "tables" / "method_status.json", {"method": method, "output_kind": "marker_score_not_abundance" if method == "marker_score" else "estimated_cell_abundance", "cross_sample_inference": "not_implemented"})
    write_adata(adata, ctx.deconv_dir() / "all_samples.deconvolution_ready.h5ad")


def _auto_distance(coords) -> float:
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    if len(coords) < 3:
        return 0.0
    nn = NearestNeighbors(n_neighbors=2).fit(coords)
    d = nn.kneighbors(coords, return_distance=True)[0][:, 1]
    return float(np.median(d) * 1.5)


def step_10_segment_objects(ctx: Context, step: str) -> None:
    require_pkg("anndata")
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    from sklearn.neighbors import radius_neighbors_graph

    adata = read_adata(latest_existing([ctx.deconv_dir() / "all_samples.deconvolution_ready.h5ad", ctx.integrated_dir() / "all_samples.annotated.h5ad"]))
    seg_conf = ctx.config["segments"]
    source = seg_conf["source_key"]
    if source not in adata.obs:
        raise ValueError(f"segments.source_key missing: {source}")
    coords = adata.obsm.get(spatial_key(ctx))
    if coords is None:
        raise SystemExit(f"Missing obsm['{spatial_key(ctx)}'] for segmentation")
    segment_ids = np.array(["other"] * adata.n_obs, dtype=object)
    scaled_dist = np.full(adata.n_obs, np.nan)
    rows: list[dict[str, Any]] = []

    for sample in adata.obs[sample_key(ctx)].astype(str).unique():
        sample_idx = np.where(adata.obs[sample_key(ctx)].astype(str).to_numpy() == sample)[0]
        if len(sample_idx) < seg_conf["min_spots"]:
            continue
        distance_threshold = seg_conf["distance_threshold"]
        radius = _auto_distance(coords[sample_idx]) if distance_threshold == "auto" else float(distance_threshold)
        if radius <= 0:
            continue
        for label in adata.obs.iloc[sample_idx][source].astype(str).unique():
            idx = sample_idx[adata.obs.iloc[sample_idx][source].astype(str).to_numpy() == label]
            if len(idx) < seg_conf["min_spots"]:
                continue
            graph = radius_neighbors_graph(coords[idx], radius=radius, mode="connectivity", include_self=False)
            n_comp, comp = connected_components(graph, directed=False)
            for comp_id in range(n_comp):
                comp_idx = idx[comp == comp_id]
                if len(comp_idx) < seg_conf["min_spots"]:
                    continue
                seg_id = f"{sample}_{label}_{comp_id}"
                segment_ids[comp_idx] = seg_id
                cent = coords[comp_idx].mean(axis=0)
                dist = np.sqrt(((coords[comp_idx] - cent) ** 2).sum(axis=1))
                denom = dist.max() if dist.max() > 0 else 1.0
                scaled_dist[comp_idx] = dist / denom
                rows.append(
                    {
                        "segment_id": seg_id,
                        "sample": sample,
                        "source_label": label,
                        "n_spots": len(comp_idx),
                        "centroid_x": cent[0],
                        "centroid_y": cent[1],
                        "radius_scaled_from": denom,
                    }
                )
    adata.obs["segment_id"] = pd.Categorical(segment_ids)
    adata.obs["scaled_dist"] = scaled_dist
    write_tsv(ctx.step_dir(step) / "tables" / "segment_summary.tsv", rows)
    if rows:
        layer = seg_conf.get("layer_key", "cell_type_pred")
        comp = adata.obs.groupby(["segment_id", layer], observed=True).size().reset_index(name="n_spots")
        comp.to_csv(ctx.step_dir(step) / "tables" / "segment_layer_composition.tsv", sep="\t", index=False)
    write_adata(adata, ctx.integrated_dir() / "all_samples.segmented.h5ad")


def step_11_programs_descriptive(ctx: Context, step: str) -> None:
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import NMF

    adata = read_adata(latest_existing([ctx.integrated_dir() / "all_samples.segmented.h5ad", ctx.integrated_dir() / "all_samples.annotated.h5ad"]))
    prog = ctx.config["programs"]
    if "highly_variable" in adata.var:
        genes_mask = adata.var["highly_variable"].fillna(False).to_numpy()
    else:
        genes_mask = pd.Series(False, index=adata.var_names).to_numpy()
    if genes_mask.sum() < 10:
        genes_mask = np.ones(adata.n_vars, dtype=bool)
    if "counts" not in adata.layers:
        raise ValueError("程序及描述汇总需要原始 counts 层，不能使用 scale 后的 X")
    counts = adata.layers["counts"]
    values = counts.data if hasattr(counts, "tocsr") else np.asarray(counts)
    if not np.all(np.isfinite(values)) or np.any(values < 0) or not np.allclose(values, np.rint(values)):
        raise ValueError("counts 必须是非负有限整数")
    totals = np.asarray(counts.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise ValueError("存在零 counts 的 spot")
    normalized = sparse_to_dense(counts) * (ctx.config["normalization"]["target_sum"] / totals[:, None])
    log_expr = np.log1p(normalized)
    x = log_expr[:, genes_mask]
    if min(x.shape) < 3:
        raise ValueError("NMF requires at least 3 spots and 3 genes")
    n_comp = min(int(prog["n_components"]), max(2, min(x.shape) - 1))
    nmf = NMF(n_components=n_comp, init="nndsvda", random_state=ctx.config["project"]["random_seed"], max_iter=500)
    w = nmf.fit_transform(x)
    write_json(ctx.step_dir(step) / "tables" / "nmf_fit_status.json",
               {"n_iter": int(nmf.n_iter_), "max_iter": int(nmf.max_iter),
                "status": "iteration_limit_reached" if nmf.n_iter_ >= nmf.max_iter else "converged",
                "reconstruction_error": float(nmf.reconstruction_err_)})
    h = nmf.components_
    prog_cols = [f"program_{i+1:02d}" for i in range(n_comp)]
    for i, col in enumerate(prog_cols):
        adata.obs[col] = w[:, i]
    genes = adata.var_names[genes_mask].to_numpy()
    top_rows: list[dict[str, Any]] = []
    for i, col in enumerate(prog_cols):
        order = np.argsort(h[i])[::-1][: prog["top_genes_per_program"]]
        for rank, j in enumerate(order, 1):
            top_rows.append({"program": col, "rank": rank, "gene": genes[j], "weight": float(h[i, j])})
    write_tsv(ctx.step_dir(step) / "tables" / "nmf_program_top_genes.tsv", top_rows)
    pd.DataFrame(w, index=adata.obs_names, columns=prog_cols).to_csv(ctx.step_dir(step) / "tables" / "nmf_program_scores.tsv", sep="\t")

    # 每个样本仅输出描述汇总；spot 不作为条件比较的独立重复。
    rows = []
    program_rows = []
    for sample in adata.obs[sample_key(ctx)].astype(str).unique():
        mask = adata.obs[sample_key(ctx)].astype(str).to_numpy() == sample
        groups = adata.obs.loc[mask, group_key(ctx)].dropna().astype(str).unique()
        if len(groups) != 1:
            raise ValueError(f"sample {sample} must map to exactly one group")
        for j, gene in enumerate(adata.var_names):
            rows.append({"sample": sample, "group": groups[0], "gene": gene,
                         "n_spots": int(mask.sum()), "mean_log1p_normalized": float(log_expr[mask,j].mean()),
                         "sum_raw_counts": float(np.asarray(counts[mask,j].sum()).item())})
        program_rows.append({"sample": sample, "group": groups[0], "n_spots": int(mask.sum()),
                             **{col: float(w[mask,i].mean()) for i,col in enumerate(prog_cols)}})
    write_tsv(ctx.step_dir(step) / "tables" / "sample_expression_descriptive.tsv", rows)
    write_tsv(ctx.step_dir(step) / "tables" / "sample_program_scores.tsv", program_rows)
    write_json(ctx.step_dir(step) / "tables" / "inference_status.json",
               {"cross_sample_inference": "not_implemented", "unit": "sample summary", "pvalues_produced": False,
                "reason": "spots/sections are not independent biological replicates; no condition test implemented"})
    write_adata(adata, ctx.integrated_dir() / "all_samples.programs.h5ad")


def step_12_neighborhood_communication(ctx: Context, step: str) -> None:
    import numpy as np
    import pandas as pd
    from sklearn.neighbors import NearestNeighbors

    adata = read_adata(latest_existing([ctx.integrated_dir() / "all_samples.programs.h5ad", ctx.integrated_dir() / "all_samples.segmented.h5ad"]))
    key = ctx.config["communication"]["cluster_key"]
    if key not in adata.obs:
        raise ValueError(f"communication.cluster_key missing: {key}")
    coords = adata.obsm.get(spatial_key(ctx))
    if coords is None:
        raise SystemExit(f"Missing obsm['{spatial_key(ctx)}'] for communication step")
    labels = adata.obs[key].astype(str).to_numpy()
    pairs: dict[tuple[str, str], int] = {}
    for sample in adata.obs[sample_key(ctx)].astype(str).unique():
        idx = np.where(adata.obs[sample_key(ctx)].astype(str).to_numpy() == sample)[0]
        if len(idx) < 3:
            continue
        nn = NearestNeighbors(n_neighbors=min(ctx.config["spatial"]["n_neighs"] + 1, len(idx))).fit(coords[idx])
        neigh = nn.kneighbors(coords[idx], return_distance=False)[:, 1:]
        for i_local, nbrs in enumerate(neigh):
            a = labels[idx[i_local]]
            for j_local in nbrs:
                b = labels[idx[j_local]]
                pairs[(a, b)] = pairs.get((a, b), 0) + 1
    rows = [{"source": a, "target": b, "n_neighbor_edges": n} for (a, b), n in sorted(pairs.items())]
    write_tsv(ctx.step_dir(step) / "tables" / "spatial_contact_edges.tsv", rows)

    if ctx.config["communication"].get("run_squidpy_ligrec", False):
        sq = require_pkg("squidpy")
        sq.gr.ligrec(adata, cluster_key=key, n_perms=100)
        if "ligrec" in adata.uns:
            adata.uns["ligrec"]["pvalues"].to_csv(ctx.step_dir(step) / "tables" / "squidpy_ligrec_pvalues.tsv", sep="\t")


def step_13_multi_sample_summary(ctx: Context, step: str) -> None:
    import pandas as pd

    adata = read_adata(latest_existing([ctx.integrated_dir() / "all_samples.programs.h5ad", ctx.integrated_dir() / "all_samples.segmented.h5ad", ctx.integrated_dir() / "all_samples.annotated.h5ad"]))
    summaries: list[dict[str, Any]] = []
    for col in ["leiden", "spatial_domain", "cell_type_pred", "segment_id"]:
        if col not in adata.obs:
            continue
        tab = adata.obs.groupby([sample_key(ctx), col], observed=True).size().reset_index(name="n_spots")
        tab.to_csv(ctx.step_dir(step) / "tables" / f"{col}_by_sample.tsv", sep="\t", index=False)
        summaries.append({"summary": col, "levels": int(tab[col].nunique()), "rows": int(tab.shape[0])})
    write_tsv(ctx.step_dir(step) / "tables" / "summary_index.tsv", summaries)


def step_14_report_index(ctx: Context, step: str) -> None:
    status_file = ctx.result_dir / "step_status.json"
    states = json.loads(status_file.read_text()) if status_file.exists() else {}
    lines = ["# 空间转录组结果索引", "", "依据当前输入签名及执行状态生成；文件存在本身不表示步骤完成。", "", "|模块|执行状态|表格|图文件|", "|---|---|---:|---:|"]
    for module in STEPS:
        mdir = ctx.result_dir / module
        status = states.get(module, {}).get("status", "not_run")
        if module == step:
            status = "completed"
        tables = sorted((mdir / "tables").glob("*")) if (mdir / "tables").exists() else []
        plots = sorted((mdir / "plots").glob("*")) if (mdir / "plots").exists() else []
        lines.append(f"|{module}|{status}|{len(tables)}|{len(plots)}|")
    lines.extend(["", "## 解释边界", "",
        "- 本流程没有实现正式的条件间推断；第 11 步只输出样本描述，不产生条件 P 值或 logFC。",
        "- spot、segment、切片均不能替代独立生物学样本。",
        "- marker score 是表达评分，不能解释为细胞丰度。",
        "- 空间域和 segment 为计算标签，需要组织学及 marker 审阅。",
        "- Cell2location 是未在本轮验证的可选分支；RCTD、SPOTlight 等方法未接入本入口。", ""])
    ensure_dir(ctx.result_dir)
    (ctx.result_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


STEP_FUNCTIONS: dict[str, Callable[[Context, str], None]] = {
    "01_validate_inputs": step_01_validate_inputs,
    "02_build_objects": step_02_build_objects,
    "03_qc_filter": step_03_qc_filter,
    "04_normalize_features": step_04_normalize_features,
    "05_integrate_cluster": step_05_integrate_cluster,
    "06_spatial_neighbors_stats": step_06_spatial_neighbors_stats,
    "07_spatial_domains": step_07_spatial_domains,
    "08_cell_type_annotation": step_08_cell_type_annotation,
    "09_deconvolution": step_09_deconvolution,
    "10_segment_objects": step_10_segment_objects,
    "11_programs_descriptive": step_11_programs_descriptive,
    "12_neighborhood_communication": step_12_neighborhood_communication,
    "13_multi_sample_summary": step_13_multi_sample_summary,
    "14_report_index": step_14_report_index,
}


def build_context(args: argparse.Namespace) -> Context:
    root = Path(args.project_root).resolve()
    config_path = Path(args.config)
    config = load_config(config_path if config_path.is_absolute() else root / config_path)
    manifest_path = Path(args.manifest or config["paths"]["manifest"])
    manifest = load_manifest(manifest_path if manifest_path.is_absolute() else root / manifest_path)
    return Context(project_root=root, config=config, manifest=manifest)


def bind_run(ctx: Context, *, write: bool) -> None:
    """A configuration/input change needs a new output tree, not leftover H5ADs."""
    step_01_validate_inputs(ctx, "01_validate_inputs", write=False)
    paths = {Path(__file__).resolve()}
    marker = ctx.path(ctx.config["paths"]["marker_gene_file"])
    if marker.is_file():
        paths.add(marker)
    ref = ctx.config["deconvolution"].get("reference_h5ad")
    if ref:
        paths.add(ctx.path(ref))
    for row in active_rows(ctx):
        if row["input_type"] == "visium_dir":
            base = resolve_input_path(ctx, row)
            paths.add(base / (row.get("count_file") or ctx.config["input"]["visium_count_file"]))
            if (base / "spatial").is_dir():
                paths.update(p for p in (base / "spatial").rglob("*") if p.is_file())
        else:
            keys = ("input_path",) if row["input_type"] == "h5ad" else ("count_file", "spatial_file")
            paths.update(resolve_input_path(ctx, row, key) for key in keys)
    hashes = {}
    for path in sorted(paths):
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(8*1024*1024), b""):
                h.update(block)
        hashes[str(path.resolve())] = h.hexdigest()
    payload = {"config": ctx.config, "manifest": ctx.manifest, "inputs": hashes}
    signature = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    target = ctx.result_dir / "run_signature.json"
    if target.is_file():
        if json.loads(target.read_text())["signature"] != signature:
            raise ValueError("STALE: 配置、源码或输入已变；使用新的输出目录/项目副本")
    else:
        outputs = [ctx.processed_dir(),ctx.integrated_dir(),ctx.deconv_dir()]
        existing = [p for folder in outputs if folder.exists() for p in folder.rglob("*.h5ad")]
        existing += [p for p in ctx.result_dir.glob("*/tables/*") if p.is_file()]
        if existing:
            raise ValueError("UNVERIFIED: 已有结果没有输入签名，不能作为本次输入")
        if write:
            write_json(target, {"signature": signature, **payload})


def execute_step(ctx: Context, step: str) -> None:
    status_file = ctx.result_dir / "step_status.json"
    states = json.loads(status_file.read_text()) if status_file.exists() else {}
    position = STEPS.index(step)
    # Each alternative set names the producers of files this step may consume.
    required = {"03_qc_filter": ["02_build_objects"], "04_normalize_features": ["03_qc_filter"],
                "05_integrate_cluster": ["04_normalize_features"], "06_spatial_neighbors_stats": ["05_integrate_cluster"],
                "07_spatial_domains": ["05_integrate_cluster"], "08_cell_type_annotation": ["07_spatial_domains"],
                "09_deconvolution": ["08_cell_type_annotation"], "10_segment_objects": ["08_cell_type_annotation"],
                "11_programs_descriptive": ["08_cell_type_annotation"],
                "12_neighborhood_communication": ["11_programs_descriptive", "10_segment_objects"],
                "13_multi_sample_summary": ["11_programs_descriptive", "10_segment_objects", "08_cell_type_annotation"]}
    if step != "14_report_index":
        invalid = [name for name in STEPS[:position] if states.get(name, {}).get("status") in {"failed", "running"}]
        if invalid:
            raise ValueError("上游失败或尚未完成，拒绝使用遗留对象: " + ",".join(invalid))
        producers = required.get(step, [])
        if producers and not any(states.get(name, {}).get("status") == "completed" for name in producers):
            raise ValueError("输入生产步骤未完成: " + ",".join(producers))
    downstream = [name for name in STEPS[position+1:] if states.get(name,{}).get("status") == "completed"]
    if downstream:
        raise ValueError("下游已经完成；重跑上游需新输出树，防止复用旧下游: " + ",".join(downstream))
    ensure_dir(ctx.step_dir(step) / "tables")
    ensure_dir(ctx.step_dir(step) / "plots")
    states[step] = {"status":"running"}; write_json(status_file,states)
    try:
        STEP_FUNCTIONS[step](ctx, step)
    except BaseException as exc:
        states[step] = {"status":"failed", "reason":str(exc)}; write_json(status_file,states)
        raise
    states[step] = {"status":"completed"}; write_json(status_file,states)


def cmd_check(args: argparse.Namespace) -> None:
    ctx = build_context(args)
    bind_run(ctx, write=False)
    print("PASS inputs; no outputs created")


def cmd_list(_: argparse.Namespace) -> None:
    for i, step in enumerate(STEPS, 1):
        print(f"{i:2d}  {step}")


def cmd_run_step(args: argparse.Namespace) -> None:
    ctx = build_context(args)
    step = args.step
    if step.isdigit():
        step = STEPS[int(step) - 1]
    if step not in STEP_FUNCTIONS:
        raise SystemExit(f"Unknown step: {args.step}")
    bind_run(ctx, write=True)
    execute_step(ctx, step)


def cmd_run_all(args: argparse.Namespace) -> None:
    ctx = build_context(args)
    bind_run(ctx, write=True)
    for step in STEPS:
        execute_step(ctx, step)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p_list = sub.add_parser("list", help="List pipeline steps")
    p_list.set_defaults(func=cmd_list)
    for name in ["run-step", "run-all", "check"]:
        p = sub.add_parser(name)
        p.add_argument("--project-root", default=".")
        p.add_argument("--config", default="config/spatial_config.yaml")
        p.add_argument("--manifest", default=None)
        if name == "run-step":
            p.add_argument("--step", required=True)
            p.set_defaults(func=cmd_run_step)
        elif name == "check":
            p.set_defaults(func=cmd_check)
        else:
            p.set_defaults(func=cmd_run_all)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv or sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()

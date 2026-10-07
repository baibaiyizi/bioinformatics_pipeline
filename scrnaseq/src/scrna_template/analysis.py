from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

from .plotting import categorical_embedding
from .runtime import RunContext, atomic_tsv





def run_atlas(context: RunContext, qc_path: str | None = None) -> str:
    input_path = context.path("02_qc", "qc_combined.h5ad") if qc_path is None else context.config.resolve(qc_path)
    output = context.path("03_atlas", "atlas.h5ad")
    marker_output = context.path("03_atlas", "cluster_markers.tsv")
    preferred = float(context.config.data["preprocessing"]["preferred_resolution"])
    preferred_key = f"leiden_{preferred:g}"
    figure_outputs = [context.path("08_figures", f"umap_{name}.{extension}") for name in ("library", "condition", "leiden") for extension in ("png", "pdf")]
    if not context.config.data["report"].get("figures", True):
        figure_outputs = []
    source_outputs = [context.path("07_source_data", f"umap_{name}.tsv") for name in ("library", "condition", "leiden")]
    with context.stage("atlas", [output, marker_output, *source_outputs, *figure_outputs]) as execute:
        if not execute:
            return str(output)
        data = sc.read_h5ad(input_path)
        data.X = data.layers["counts"].copy().tocsr()
        sc.pp.filter_genes(data, min_cells=int(context.config.data["qc"]["min_cells_per_gene"]))
        sc.pp.normalize_total(data, target_sum=float(context.config.data["preprocessing"]["target_sum"]))
        sc.pp.log1p(data)
        data.raw = data
        n_top = min(int(context.config.data["preprocessing"]["n_hvg"]), max(50, data.n_vars - 1))
        sc.pp.highly_variable_genes(data, n_top_genes=n_top, flavor="seurat", subset=False)
        excluded = data.var["mt"] | data.var_names.isin(context.config.data["preprocessing"].get("exclude_hvg_genes", []))
        data.var.loc[excluded, "highly_variable"] = False
        if int(data.var["highly_variable"].sum()) < 20:
            raise RuntimeError("排除 mt/显式基因后 HVG 少于 20")
        n_pcs = min(int(context.config.data["preprocessing"]["n_pcs"]), data.n_obs - 2, int(data.var["highly_variable"].sum()) - 1)
        if n_pcs < 2:
            raise RuntimeError("细胞数不足以构建 PCA")
        sc.pp.pca(data, n_comps=n_pcs, mask_var="highly_variable", random_state=int(context.config.data["project"]["seed"]))
        neighbors = min(int(context.config.data["preprocessing"]["n_neighbors"]), max(2, data.n_obs - 1))
        sc.pp.neighbors(data, n_neighbors=neighbors, n_pcs=n_pcs, random_state=int(context.config.data["project"]["seed"]))
        sc.tl.umap(data, random_state=int(context.config.data["project"]["seed"]))
        for resolution in context.config.data["preprocessing"]["leiden_resolutions"]:
            sc.tl.leiden(data, resolution=float(resolution), key_added=f"leiden_{float(resolution):g}", random_state=int(context.config.data["project"]["seed"]), flavor="igraph", n_iterations=2)
        if preferred_key not in data.obs:
            raise RuntimeError(f"preferred Leiden 不存在: {preferred_key}")
        sc.tl.rank_genes_groups(data, groupby=preferred_key, method="wilcoxon", use_raw=True, pts=True)
        tables = []
        for cluster in data.obs[preferred_key].cat.categories:
            frame = sc.get.rank_genes_groups_df(data, group=cluster)
            frame.insert(0, "cluster", str(cluster))
            frame["eligible_for_annotation"] = (frame["logfoldchanges"].abs() >= float(context.config.data["statistics"]["marker_min_log2fc"])) & (frame["pvals"] < float(context.config.data["statistics"]["marker_pvalue"]))
            frame["test_scope"] = "cluster_marker_annotation_only"
            tables.append(frame)
        markers = pd.concat(tables, ignore_index=True)
        atomic_tsv(marker_output, markers)
        dpi = int(context.config.data["report"]["dpi"])
        if context.config.data["report"].get("figures", True):
            categorical_embedding(data.obsm["X_umap"], data.obs["library_id"], "Library", context.path("08_figures", "umap_library"), dpi)
            categorical_embedding(data.obsm["X_umap"], data.obs["condition"], "Condition (descriptive)", context.path("08_figures", "umap_condition"), dpi)
            categorical_embedding(data.obsm["X_umap"], data.obs[preferred_key], f"Leiden {preferred:g}", context.path("08_figures", "umap_leiden"), dpi)
        base_source = pd.DataFrame({"cell_id": data.obs_names.astype(str), "UMAP1": data.obsm["X_umap"][:, 0], "UMAP2": data.obsm["X_umap"][:, 1]})
        for name, values in (("library", data.obs["library_id"]), ("condition", data.obs["condition"]), ("leiden", data.obs[preferred_key])):
            source = base_source.copy()
            source["label"] = values.astype(str).to_numpy()
            atomic_tsv(context.path("07_source_data", f"umap_{name}.tsv"), source)
        data.uns["analysis_contract"] = {"integration": "not_used", "condition_inference": "see_inference_status", "marker_min_log2fc": float(context.config.data["statistics"]["marker_min_log2fc"]), "marker_pvalue": float(context.config.data["statistics"]["marker_pvalue"])}
        data.write_h5ad(output, compression="gzip")
    return str(output)

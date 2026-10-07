from __future__ import annotations

import logging
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from .plotting import categorical_embedding
from .runtime import RunContext, atomic_json, atomic_text, atomic_tsv, run_command


LOGGER = logging.getLogger("pfos_scrna")


def _mean_score(data: ad.AnnData, genes: list[str]) -> np.ndarray:
    present = [gene for gene in genes if gene in data.var_names]
    if not present:
        return np.full(data.n_obs, np.nan)
    return np.asarray(data[:, present].X.mean(axis=1)).ravel()


def _marker_labels(data: ad.AnnData, marker_config: dict[str, object], cluster_key: str, min_fraction: float, min_margin: float) -> pd.DataFrame:
    rows = []
    counts = data.layers["counts"]
    cluster_values = data.obs[cluster_key].astype(str)
    clusters = sorted(cluster_values.unique(), key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value))
    for level, marker_sets in marker_config.items():
        marker_genes = list(dict.fromkeys(
            gene
            for spec in marker_sets.values()
            for gene in [*spec.get("positive", []), *spec.get("negative", [])]
            if gene in data.var_names
        ))
        if not marker_genes:
            raise RuntimeError(f"{level} 没有可用 marker")
        centroids = np.vstack([
            np.asarray(data[cluster_values.to_numpy() == cluster, marker_genes].X.mean(axis=0)).ravel()
            for cluster in clusters
        ])
        means = centroids.mean(axis=0)
        scales = centroids.std(axis=0)
        standardized = np.divide(centroids - means, scales, out=np.zeros_like(centroids, dtype=float), where=scales > 0)
        gene_index = {gene: index for index, gene in enumerate(marker_genes)}
        score_columns: list[tuple[str, str, list[str], np.ndarray]] = []
        for label, spec in marker_sets.items():
            positive = [gene for gene in spec.get("positive", []) if gene in gene_index]
            negative = [gene for gene in spec.get("negative", []) if gene in gene_index]
            if not positive:
                score = np.full(len(clusters), -np.inf)
            else:
                score = standardized[:, [gene_index[gene] for gene in positive]].mean(axis=1)
                if negative:
                    score -= standardized[:, [gene_index[gene] for gene in negative]].mean(axis=1)
            column = f"marker_score_{level}_{label}"
            data.obs[column] = cluster_values.map(dict(zip(clusters, score, strict=True))).astype(float)
            score_columns.append((label, column, positive, score))
        for cluster_index, cluster in enumerate(clusters):
            mask = cluster_values.to_numpy() == cluster
            ranked = sorted(((label, float(score[cluster_index])) for label, _, _, score in score_columns), key=lambda item: item[1], reverse=True)
            label, top = ranked[0]
            second = ranked[1][1] if len(ranked) > 1 else -np.inf
            positive = next(genes for candidate, _, genes, _ in score_columns if candidate == label)
            support = 0.0
            if positive:
                subset = counts[mask][:, [data.var_names.get_loc(gene) for gene in positive]]
                detection = np.asarray((subset > 0).mean(axis=0)).ravel()
                support = float((detection >= 0.1).mean())
            confidence = float(top - second)
            accepted = support >= min_fraction and confidence >= min_margin
            rows.append({"level": level, "cluster": cluster, "marker_label": label if accepted else "unknown", "top_candidate": label, "marker_score": top, "score_margin": confidence, "positive_marker_support": support, "accepted": accepted, "score_method": "cluster_mean_log1p_gene_zscore_positive_minus_negative"})
    return pd.DataFrame(rows)


def _mia_broad_label(label: str) -> str:
    value = str(label).strip().lower().replace(" ", "_")
    if any(token in value for token in ("beta", "alpha", "delta", "gamma", "epsilon", "endocrine")):
        return "endocrine"
    if "duct" in value:
        return "ductal"
    if any(token in value for token in ("acinar", "exocrine")):
        return "acinar"
    if any(token in value for token in ("endothelial", "vascular")):
        return "endothelial"
    if any(token in value for token in ("stellate", "mesench", "pericyte", "fibroblast")):
        return "mesenchymal"
    if any(token in value for token in ("immune", "macroph", "monocyte", "lymph", "t_cell", "b_cell", "mast")):
        return "immune"
    if any(token in value for token in ("schwann", "neural", "neuron", "glial")):
        return "neural"
    return "unknown"


def _mia_subtype_compatible(marker_label: str, mia_label: str) -> bool:
    marker = str(marker_label).lower()
    mia = str(mia_label).lower()
    if mia in {"endocrine", "unknown", "not_available"}:
        return True
    if "beta_delta" in mia or ("beta" in mia and "delta" in mia):
        return marker in {"beta", "delta"}
    return marker in mia


def _annotation_decisions(evidence: pd.DataFrame, mia: pd.DataFrame) -> pd.DataFrame:
    broad = evidence[evidence["level"] == "broad"].set_index("cluster")
    endocrine = evidence[evidence["level"] == "endocrine"].set_index("cluster")
    mia_indexed = mia.assign(cluster=mia["cluster"].astype(str)).set_index("cluster") if not mia.empty else pd.DataFrame()
    rows = []
    for cluster in sorted(broad.index.astype(str), key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value)):
        broad_label = str(broad.loc[cluster, "marker_label"])
        subtype_label = str(endocrine.loc[cluster, "marker_label"]) if broad_label == "endocrine" and cluster in endocrine.index else "not_applicable"
        marker_label = subtype_label if broad_label == "endocrine" else broad_label
        final_label = marker_label
        status = "marker_only_reference_unavailable"
        reason = "marker hierarchy accepted; MIA unavailable"
        mia_label = "not_available"
        mia_broad = "unknown"
        if marker_label == "unknown":
            status = "low_confidence_unknown"
            reason = "marker support or margin below frozen threshold"
        elif not mia.empty and cluster in mia_indexed.index:
            mia_label = str(mia_indexed.loc[cluster, "mia_label"])
            mia_broad = _mia_broad_label(mia_label)
            if mia_broad != "unknown" and mia_broad != broad_label:
                final_label = "unknown"
                status = "broad_marker_mia_conflict_unknown"
                reason = f"marker broad={broad_label}; MIA broad={mia_broad}; MIA vetoes but does not replace"
            elif broad_label == "endocrine" and not _mia_subtype_compatible(marker_label, mia_label):
                final_label = "unknown"
                status = "endocrine_marker_mia_conflict_unknown"
                reason = f"marker subtype={marker_label}; MIA={mia_label}; MIA vetoes but does not replace"
            else:
                status = "marker_mia_concordant"
                reason = "marker hierarchy accepted and MIA compatible"
        rows.append({"cluster": cluster, "broad_marker_label": broad_label, "endocrine_marker_label": subtype_label, "marker_label": marker_label, "mia_label": mia_label, "mia_broad_label": mia_broad, "final_label": final_label, "annotation_review_status": status, "decision_reason": reason})
    return pd.DataFrame(rows)


def _mia_support(context: RunContext, data: ad.AnnData, cluster_key: str) -> tuple[pd.DataFrame, dict[str, object]]:
    atlas_path = Path(context.config.data["reference"]["mia_atlas"])
    label_key = context.config.data["annotation"]["mia_label_key"]
    per_label = int(context.config.data["annotation"]["mia_cells_per_label"])
    try:
        reference = ad.read_h5ad(atlas_path, backed="r")
        if label_key not in reference.obs:
            raise KeyError(label_key)
        symbol_key = "gene_symbol" if "gene_symbol" in reference.var else "gene_symbol_FINAL"
        symbols = reference.var[symbol_key].astype(str).to_numpy()
        query_symbols = set(data.var_names.astype(str))
        unique = ~pd.Index(symbols).duplicated(keep="first")
        indices = np.flatnonzero(unique & np.isin(symbols, list(query_symbols)))
        if len(indices) < 200:
            raise RuntimeError(f"MIA/query common genes too few: {len(indices)}")
        labels = reference.obs[label_key].astype(str)
        rng = np.random.default_rng(int(context.config.data["project"]["seed"]))
        selected = []
        for label in sorted(labels.unique()):
            members = np.flatnonzero(labels.to_numpy() == label)
            selected.extend(rng.choice(members, size=min(per_label, len(members)), replace=False).tolist())
        reference_memory = reference[selected, indices].to_memory()
        selected_symbols = symbols[indices]
        reference_memory.var_names = pd.Index(selected_symbols)
        reference_labels = labels.iloc[selected].to_numpy()
        query_indices = [data.var_names.get_loc(gene) for gene in selected_symbols]
        query_matrix = data.X[:, query_indices]
        ref_matrix = reference_memory.X
        rows = []
        ref_centroids = {label: np.asarray(ref_matrix[reference_labels == label].mean(axis=0)).ravel() for label in sorted(set(reference_labels))}
        for cluster in sorted(data.obs[cluster_key].astype(str).unique()):
            mask = data.obs[cluster_key].astype(str).to_numpy() == cluster
            query_centroid = np.asarray(query_matrix[mask].mean(axis=0)).ravel()
            scores = []
            for label, centroid in ref_centroids.items():
                correlation = float(np.corrcoef(query_centroid, centroid)[0, 1]) if np.std(query_centroid) and np.std(centroid) else -1.0
                scores.append((label, correlation))
            scores.sort(key=lambda item: item[1], reverse=True)
            rows.append({"cluster": cluster, "mia_label": scores[0][0], "mia_correlation": scores[0][1], "mia_margin": scores[0][1] - scores[1][1] if len(scores) > 1 else np.nan})
        reference.file.close()
        return pd.DataFrame(rows), {"status": "completed", "method": "MIA atlas sampled centroid correlation", "atlas": str(atlas_path), "cells_per_label": per_label, "common_genes": len(indices), "overwrites_marker": False}
    except Exception as exc:
        return pd.DataFrame(columns=["cluster", "mia_label", "mia_correlation", "mia_margin"]), {"status": "failed_marker_retained", "reason": str(exc), "overwrites_marker": False}


def _singler_support(context: RunContext, data: ad.AnnData, cluster_key: str) -> tuple[pd.DataFrame, dict[str, object]]:
    exchange = context.path("04_annotation", "singler_exchange")
    exchange.mkdir(parents=True, exist_ok=True)
    matrix_path = exchange / "matrix.mtx"
    feature_path = exchange / "features.txt"
    barcode_path = exchange / "barcodes.txt"
    cluster_path = exchange / "clusters.tsv"
    output = context.path("04_annotation", "singler_support.tsv")
    try:
        cluster_values = data.obs[cluster_key].astype(str).to_numpy()
        clusters = sorted(set(cluster_values), key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value))
        cluster_index = {cluster: index for index, cluster in enumerate(clusters)}
        membership = sp.csc_matrix(
            (
                np.ones(data.n_obs, dtype=np.float32),
                (np.arange(data.n_obs), np.asarray([cluster_index[value] for value in cluster_values])),
            ),
            shape=(data.n_obs, len(clusters)),
        )
        expression = data.X if sp.issparse(data.X) else sp.csr_matrix(data.X)
        cluster_sums = (expression.T @ membership).tocoo()
        scipy.io.mmwrite(matrix_path, cluster_sums)
        atomic_text(feature_path, "\n".join(data.var_names.astype(str)) + "\n")
        atomic_text(barcode_path, "\n".join(clusters) + "\n")
        atomic_tsv(cluster_path, pd.DataFrame({"barcode": clusters, "cluster": clusters}))
        command_env = os.environ.copy()
        removed_proxy = []
        for key in ("HTTP_PROXY", "HTTPS_PROXY", "http_proxy", "https_proxy"):
            if command_env.get(key, "").rstrip("/") == "http://127.0.0.1:7890":
                removed_proxy.append(key)
                command_env.pop(key)
        command_env["PFOS_SINGLER_THREADS"] = str(int(context.config.data["resources"]["singler_threads"]))
        output.unlink(missing_ok=True)
        run_command([context.config.data["tools"]["rscript"], context.config.root / "scripts" / "run_singler.R", matrix_path, feature_path, barcode_path, cluster_path, output], log_file=context.path("logs", "singler.log"), env=command_env)
        result = pd.read_csv(output, sep="\t", dtype=str)
        required = {"cluster", "singler_label", "singler_pruned_label", "singler_delta_next"}
        if not required <= set(result.columns) or result["cluster"].duplicated().any() or set(result["cluster"]) != set(clusters):
            raise RuntimeError("SingleR cluster output contract failed")
        return result, {"status": "completed", "reference": "celldex::MouseRNAseqData", "method": "SingleR on cluster-summed normalize-total log1p expression", "input_cells": data.n_obs, "input_clusters": len(clusters), "role": "weak_auxiliary", "overwrites_marker": False, "threads": int(command_env["PFOS_SINGLER_THREADS"]), "removed_inactive_loopback_proxy": removed_proxy}
    except Exception as exc:
        return pd.DataFrame(columns=["cluster", "singler_label", "singler_pruned_label", "singler_delta_next"]), {"status": "failed_marker_retained", "reason": str(exc), "overwrites_marker": False}


def _module_scores(context: RunContext, data: ad.AnnData) -> pd.DataFrame:
    rows = []
    for set_id, frame in context.config.gene_sets().groupby("set_id", sort=True):
        genes = frame["gene_symbol"].tolist()
        present = [gene for gene in genes if gene in data.var_names]
        column = f"module_{set_id}"
        data.obs[column] = _mean_score(data, present)
        rows.append({"set_id": set_id, "configured_genes": len(genes), "present_genes": len(present), "genes_present": ";".join(present), "score_definition": "mean normalize-total log1p expression; continuous, no forced class"})
    return pd.DataFrame(rows)


def run_annotation(context: RunContext, atlas_path: str | None = None) -> str:
    input_path = context.path("03_atlas", "atlas.h5ad") if atlas_path is None else Path(atlas_path)
    output = context.path("04_annotation", "annotated.h5ad")
    evidence_path = context.path("04_annotation", "annotation_evidence.tsv")
    decision_path = context.path("04_annotation", "annotation_decisions.tsv")
    module_path = context.path("04_annotation", "module_score_contract.tsv")
    figure_outputs = [context.path("08_figures", f"umap_cell_type.{extension}") for extension in ("png", "pdf")]
    source_output = context.path("07_source_data", "umap_cell_type.tsv")
    with context.stage("annotation", [output, evidence_path, decision_path, module_path, source_output, *figure_outputs]) as execute:
        if not execute:
            return str(output)
        data = ad.read_h5ad(input_path)
        preferred = float(context.config.data["preprocessing"]["preferred_resolution"])
        cluster_key = f"leiden_{preferred:g}"
        evidence = _marker_labels(data, context.config.markers(), cluster_key, float(context.config.data["annotation"]["min_marker_fraction"]), float(context.config.data["annotation"]["min_score_margin"]))
        mia, mia_status = _mia_support(context, data, cluster_key) if context.config.data["annotation"]["mia_support"] else (pd.DataFrame(), {"status": "skipped"})
        singler, singler_status = _singler_support(context, data, cluster_key) if context.config.data["annotation"]["singleR"] else (pd.DataFrame(), {"status": "skipped"})
        atomic_json(context.path("04_annotation", "reference_annotation_status.json"), {"mia": mia_status, "singler": singler_status})
        if not mia.empty:
            evidence = evidence.merge(mia, on="cluster", how="left")
        if not singler.empty:
            evidence = evidence.merge(singler, on="cluster", how="left")
        evidence["final_rule"] = "marker hierarchy is primary; MIA may veto incompatible labels but never supplies a replacement; low confidence or conflict remains unknown"
        atomic_tsv(evidence_path, evidence)
        decisions = _annotation_decisions(evidence, mia)
        atomic_tsv(decision_path, decisions)
        decision_index = decisions.set_index("cluster")
        cluster_values = data.obs[cluster_key].astype(str)
        data.obs["marker_label"] = cluster_values.map(decision_index["marker_label"]).fillna("unknown").astype("category")
        data.obs["cell_type"] = cluster_values.map(decision_index["final_label"]).fillna("unknown").astype("category")
        if not mia.empty:
            mia_by_cluster = mia.set_index("cluster")
            data.obs["mia_label"] = data.obs[cluster_key].astype(str).map(mia_by_cluster["mia_label"]).fillna("not_available")
            data.obs["mia_support"] = data.obs[cluster_key].astype(str).map(mia_by_cluster["mia_correlation"]).astype(float)
        else:
            data.obs["mia_label"] = "not_available"
            data.obs["mia_support"] = np.nan
        if not singler.empty:
            singler_by_cluster = singler.set_index("cluster")
            data.obs["singler_label"] = data.obs[cluster_key].astype(str).map(singler_by_cluster["singler_pruned_label"]).fillna("not_available")
        else:
            data.obs["singler_label"] = "not_available"
        broad_evidence = evidence[evidence["level"] == "broad"].set_index("cluster")
        endocrine_evidence = evidence[evidence["level"] == "endocrine"].set_index("cluster")
        confidence = cluster_values.map(broad_evidence["score_margin"]).astype(float)
        endocrine_mask = cluster_values.map(broad_evidence["marker_label"]).astype(str) == "endocrine"
        confidence.loc[endocrine_mask] = cluster_values.loc[endocrine_mask].map(endocrine_evidence["score_margin"]).astype(float)
        data.obs["annotation_confidence"] = confidence
        data.obs["annotation_source"] = "cluster_zscored_marker_hierarchy_with_non_overwriting_mia_veto"
        data.obs["annotation_review_status"] = cluster_values.map(decision_index["annotation_review_status"]).fillna("low_confidence_unknown")
        data.obs["condition_inference_eligible"] = False
        hormones = [gene for gene in ("Ins1", "Ins2", "Gcg", "Sst", "Ppy") if gene in data.var_names]
        if hormones:
            detection = np.column_stack([np.asarray(data.layers["counts"][:, data.var_names.get_loc(gene)].toarray()).ravel() > 0 for gene in hormones])
            data.obs["lineage_conflict_review"] = detection.sum(axis=1) >= 2
        else:
            data.obs["lineage_conflict_review"] = False
        module_contract = _module_scores(context, data)
        atomic_tsv(module_path, module_contract)
        dpi = int(context.config.data["report"]["dpi"])
        categorical_embedding(data.obsm["X_umap"], data.obs["cell_type"], "Marker-first cell type", context.path("08_figures", "umap_cell_type"), dpi)
        atomic_tsv(source_output, pd.DataFrame({"cell_id": data.obs_names.astype(str), "UMAP1": data.obsm["X_umap"][:, 0], "UMAP2": data.obsm["X_umap"][:, 1], "cell_type": data.obs["cell_type"].astype(str), "marker_label": data.obs["marker_label"].astype(str), "mia_label": data.obs["mia_label"].astype(str), "annotation_review_status": data.obs["annotation_review_status"].astype(str)}))
        data.write_h5ad(output, compression="gzip")
    return str(output)

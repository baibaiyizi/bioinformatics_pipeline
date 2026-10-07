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


LOGGER = logging.getLogger("scrna_template")


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



def run_annotation(context, atlas_path):
    output = context.path('04_annotation', 'annotated.h5ad')
    evidence_path = context.path('04_annotation', 'annotation_evidence.tsv')
    with context.stage('annotation', [output, evidence_path]):
        data = ad.read_h5ad(atlas_path)
        key = f"leiden_{float(context.config.data['preprocessing']['preferred_resolution']):g}"
        marker_sets = context.config.markers()['cell_types']
        if marker_sets:
            spec = context.config.data['annotation']
            evidence = _marker_labels(data, {'cell_types': marker_sets}, key, float(spec['min_fraction']), float(spec['min_margin']))
            labels = evidence.set_index('cluster')['marker_label']
            data.obs['cell_type'] = data.obs[key].astype(str).map(labels).fillna('unknown')
            evidence['scope'] = 'marker-based annotation support; requires biological review'
        else:
            data.obs['cell_type'] = 'unknown'
            evidence = pd.DataFrame([{'cluster': c, 'marker_label': 'unknown', 'accepted': False, 'scope': 'markers not configured'} for c in data.obs[key].astype(str).unique()])
        atomic_tsv(evidence_path, evidence)
        data.write_h5ad(output, compression='gzip')
    return output

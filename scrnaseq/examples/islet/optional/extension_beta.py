from __future__ import annotations

import json
import logging
import platform
import sys
from dataclasses import dataclass
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from time import monotonic
from typing import Any, Callable

import anndata as ad
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from scipy.sparse.csgraph import connected_components
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

from .plotting import DIVERGING_CMAP, save_figure, category_color
from .figure_style import PROFILE, COLORS
from .runtime import atomic_json, atomic_text, atomic_tsv, sha256


LOGGER = logging.getLogger("pfos_scrna.extension_beta")
SEED = 42
MAIN_NEIGHBORS = 30
MAIN_PCS = 40
LEIDEN_RESOLUTIONS = (0.3, 0.5, 0.8)
PREFERRED_CLUSTER_KEY = "beta_leiden_0.5"
CONDITION_COLORS = COLORS
BETA_CLUSTER_COLORS = PROFILE["palettes"]["beta_clusters"]
CELL_CYCLE_COLORS = PROFILE["palettes"]["cell_cycle"]

# Frozen before looking at condition differences. These programs are scored as
# continuous covariates with expression-matched control genes. They are not
# mutually exclusive classes and no cell is assigned to its maximum score.
STATE_GENE_SETS: dict[str, tuple[str, ...]] = {
    "maturity": (
        "Mafa", "Ucn3", "Slc2a2", "Pdx1", "Nkx6-1", "Nkx2-2", "Isl1",
        "Pcsk1", "Ero1b", "Iapp", "Cpe", "Abcc8", "Kcnj11",
    ),
    "gsis": (
        "Gck", "Slc2a2", "Abcc8", "Kcnj11", "Cacna1d", "Cacna1c",
        "Snap25", "Stx1a", "Syt7", "Pcsk1", "Pcsk2", "Cpe",
    ),
    "er_stress": (
        "Hspa5", "Ddit3", "Xbp1", "Atf4", "Atf6", "Ern1", "Eif2ak3",
        "Herpud1", "Hyou1", "Pdia4", "Dnajb9", "Manf",
    ),
    "oxidative_stress": (
        "Nfe2l2", "Hmox1", "Nqo1", "Gpx1", "Gpx2", "Sod1", "Sod2",
        "Txnrd1", "Prdx1", "Gclc", "Gclm", "Srxn1",
    ),
    "dedifferentiation": (
        "Aldh1a3", "Sox9", "Neurog3", "Gast", "Krt8", "Krt18", "Vim",
        "S100a10", "Cd81", "Cldn4", "Igfbp7", "Ldha",
    ),
}

# Root selection and directional diagnostics use disjoint genes. This prevents
# the state-distance gate from validating a direction with the same expression
# variables that defined the root.
ROOT_ANCHOR_GENE_SETS: dict[str, tuple[str, ...]] = {
    "root_maturity": ("Mafa", "Ucn3", "Slc2a2", "Pdx1", "Nkx6-1"),
    "root_stress": ("Hspa5", "Ddit3", "Xbp1", "Hmox1", "Nqo1", "Sod2", "Aldh1a3", "Sox9", "Gast"),
}

HELDOUT_DIRECTION_GENE_SETS: dict[str, tuple[str, ...]] = {
    "direction_maturity_heldout": ("Nkx2-2", "Isl1", "Pcsk1", "Ero1b", "Iapp", "Cpe", "Abcc8", "Kcnj11"),
    "direction_stress_heldout": (
        "Atf4", "Atf6", "Ern1", "Eif2ak3", "Herpud1", "Hyou1", "Pdia4", "Dnajb9", "Manf",
        "Gpx1", "Gpx2", "Sod1", "Txnrd1", "Prdx1", "Gclc", "Gclm", "Srxn1",
        "Neurog3", "Krt8", "Krt18", "Vim", "S100a10", "Cd81", "Cldn4", "Igfbp7", "Ldha",
    ),
}

S_GENES: tuple[str, ...] = (
    "Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung",
    "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1", "Cenpu",
    "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76",
    "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51",
    "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm",
    "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", "E2f8",
)

G2M_GENES: tuple[str, ...] = (
    "Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a",
    "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf",
    "Tacc3", "Pimreg", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb",
    "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp",
    "Cdca3", "Jpt1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1",
    "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmr",
    "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf",
    "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa",
)

HVG_EXCLUDE = {
    "Ins1", "Ins2", "Iapp", "Gcg", "Sst", "Ppy", "Pyy", "Ghrl",
}

# These panels are an annotation/QC audit, not a replacement cell-typing rule.
# The identity transcription factors are evaluated separately from abundant
# beta-cell products because insulin-family transcripts can also be ambient.
BETA_IDENTITY_GENES: tuple[str, ...] = (
    "Pdx1", "Nkx6-1", "Mafa", "Isl1", "Nkx2-2", "Ucn3", "Slc2a2",
)
BETA_PRODUCT_GENES: tuple[str, ...] = ("Ins1", "Ins2", "Iapp")
NONCOGNATE_HORMONES: tuple[str, ...] = ("Gcg", "Sst", "Ppy", "Pyy", "Ghrl")
ALTERNATIVE_LINEAGE_TFS: tuple[str, ...] = ("Arx", "Hhex")
IDENTITY_REVIEW_CLUSTER = "4"


@dataclass(frozen=True)
class TrajectoryRun:
    name: str
    neighbor_count: int
    n_pcs: int
    seed: int
    root_index: int


def _ensure_sparse_csr(matrix: Any, name: str) -> sp.csr_matrix:
    if not sp.issparse(matrix):
        raise TypeError(f"{name} must remain sparse")
    result = matrix.tocsr()
    if result.data.size and (not np.isfinite(result.data).all() or (result.data < 0).any()):
        raise ValueError(f"{name} contains non-finite or negative values")
    return result


def _load_beta(input_h5ad: Path) -> ad.AnnData:
    source = ad.read_h5ad(input_h5ad, backed="r")
    try:
        required_obs = {"cell_type", "library_id", "condition", "total_counts", "n_genes_by_counts", "pct_counts_mt"}
        missing = required_obs.difference(source.obs.columns)
        if missing:
            raise ValueError(f"input missing required obs columns: {sorted(missing)}")
        if not {"counts", "counts_raw"}.issubset(source.layers.keys()):
            raise ValueError("input must contain counts and counts_raw layers")
        mask = source.obs["cell_type"].astype(str).eq("beta").to_numpy()
        if int(mask.sum()) < 100:
            raise ValueError(f"too few beta cells: {int(mask.sum())}")
        beta = source[mask, :].to_memory()
    finally:
        source.file.close()

    if not beta.obs_names.is_unique or not beta.var_names.is_unique:
        raise ValueError("beta output requires globally unique barcodes and genes")
    beta.X = _ensure_sparse_csr(beta.X, "X")
    beta.layers["counts"] = _ensure_sparse_csr(beta.layers["counts"], "counts")
    beta.layers["counts_raw"] = _ensure_sparse_csr(beta.layers["counts_raw"], "counts_raw")
    if beta.layers["counts"].data.size and not np.allclose(beta.layers["counts"].data, np.rint(beta.layers["counts"].data), atol=1e-8):
        raise ValueError("counts layer must contain integer-valued UMI counts")
    if beta.layers["counts_raw"].data.size and not np.allclose(beta.layers["counts_raw"].data, np.rint(beta.layers["counts_raw"].data), atol=1e-8):
        raise ValueError("counts_raw layer must contain integer-valued UMI counts")

    # Remove parent-atlas representations whose coordinates and graphs are no
    # longer valid after lineage subsetting. Expression and count layers remain.
    beta.obsm.clear()
    beta.varm.clear()
    beta.obsp.clear()
    beta.uns.clear()
    beta.obs["beta_analysis_scope"] = "all_QC-passing_beta_cells"
    beta.obs["condition_inference_eligible"] = False
    return beta


def _present_genes(data: ad.AnnData, genes: tuple[str, ...], minimum: int, label: str) -> list[str]:
    present = [gene for gene in genes if gene in data.var_names]
    if len(present) < minimum:
        raise ValueError(f"{label} has only {len(present)} present genes; requires at least {minimum}")
    return present


def _score_programs(data: ad.AnnData) -> pd.DataFrame:
    all_signal_genes = set().union(
        *STATE_GENE_SETS.values(),
        *ROOT_ANCHOR_GENE_SETS.values(),
        *HELDOUT_DIRECTION_GENE_SETS.values(),
        S_GENES,
        G2M_GENES,
    )
    pool = [
        gene for gene in data.var_names
        if not str(gene).startswith("mt-")
        and not str(gene).startswith("Rpl")
        and not str(gene).startswith("Rps")
        and (gene not in HVG_EXCLUDE or gene in all_signal_genes)
    ]
    rows: list[dict[str, Any]] = []
    for set_id, configured in STATE_GENE_SETS.items():
        present = _present_genes(data, configured, minimum=4, label=set_id)
        score_name = f"beta_score_{set_id}"
        sc.tl.score_genes(
            data,
            present,
            ctrl_size=max(50, len(present)),
            gene_pool=pool,
            n_bins=25,
            score_name=score_name,
            random_state=SEED,
            ctrl_as_ref=False,
            use_raw=False,
        )
        rows.append({
            "set_id": set_id,
            "configured_genes": len(configured),
            "present_genes": len(present),
            "genes_present": ";".join(present),
            "genes_missing": ";".join(gene for gene in configured if gene not in present),
            "score_definition": "Scanpy 1.12.4 score_genes on normalize-total log1p X; expression-bin matched controls; ctrl_as_ref=False; random_state=42",
            "interpretation": "continuous exploratory program score; no maximum-score cell classification",
        })

    root_genes = set().union(*ROOT_ANCHOR_GENE_SETS.values())
    heldout_genes = set().union(*HELDOUT_DIRECTION_GENE_SETS.values())
    if root_genes & heldout_genes:
        raise ValueError("root anchors and held-out direction genes must be disjoint")
    for set_id, configured in {**ROOT_ANCHOR_GENE_SETS, **HELDOUT_DIRECTION_GENE_SETS}.items():
        present = _present_genes(data, configured, minimum=3, label=set_id)
        score_name = f"beta_{set_id}_score"
        sc.tl.score_genes(
            data,
            present,
            ctrl_size=max(50, len(present)),
            gene_pool=pool,
            n_bins=25,
            score_name=score_name,
            random_state=SEED,
            ctrl_as_ref=False,
            use_raw=False,
        )
        rows.append({
            "set_id": set_id,
            "configured_genes": len(configured),
            "present_genes": len(present),
            "genes_present": ";".join(present),
            "genes_missing": ";".join(gene for gene in configured if gene not in present),
            "score_definition": "Scanpy 1.12.4 score_genes on normalize-total log1p X; expression-bin matched controls; ctrl_as_ref=False; random_state=42",
            "interpretation": "root-only anchor" if set_id.startswith("root_") else "held-out direction diagnostic; genes disjoint from root anchors",
        })

    s_present = _present_genes(data, S_GENES, minimum=10, label="cell_cycle_S")
    g2m_present = _present_genes(data, G2M_GENES, minimum=10, label="cell_cycle_G2M")
    sc.tl.score_genes_cell_cycle(
        data,
        s_genes=s_present,
        g2m_genes=g2m_present,
        random_state=SEED,
        ctrl_as_ref=False,
        use_raw=False,
    )
    data.obs.rename(columns={"S_score": "beta_S_score", "G2M_score": "beta_G2M_score", "phase": "beta_cell_cycle_phase"}, inplace=True)
    rows.extend([
        {
            "set_id": "cell_cycle_S",
            "configured_genes": len(S_GENES),
            "present_genes": len(s_present),
            "genes_present": ";".join(s_present),
            "genes_missing": ";".join(gene for gene in S_GENES if gene not in s_present),
            "score_definition": "Scanpy 1.12.4 score_genes_cell_cycle on normalize-total log1p X; random_state=42",
            "interpretation": "descriptive cell-cycle score; phase does not establish proliferation",
        },
        {
            "set_id": "cell_cycle_G2M",
            "configured_genes": len(G2M_GENES),
            "present_genes": len(g2m_present),
            "genes_present": ";".join(g2m_present),
            "genes_missing": ";".join(gene for gene in G2M_GENES if gene not in g2m_present),
            "score_definition": "Scanpy 1.12.4 score_genes_cell_cycle on normalize-total log1p X; random_state=42",
            "interpretation": "descriptive cell-cycle score; phase does not establish proliferation",
        },
    ])
    return pd.DataFrame(rows)


def _build_atlas(data: ad.AnnData) -> None:
    n_top = min(3000, data.n_vars - 1)
    sc.pp.highly_variable_genes(data, n_top_genes=n_top, flavor="seurat", subset=False)
    excluded = data.var_names.isin(HVG_EXCLUDE)
    for column in ("mt", "ribo", "hb"):
        if column in data.var:
            excluded = excluded | data.var[column].astype(bool).to_numpy()
    data.var.loc[excluded, "highly_variable"] = False
    if int(data.var["highly_variable"].sum()) < 100:
        raise ValueError("fewer than 100 eligible beta-lineage HVGs")
    n_pcs = min(MAIN_PCS, data.n_obs - 2, int(data.var["highly_variable"].sum()) - 1)
    sc.pp.pca(data, n_comps=n_pcs, mask_var="highly_variable", random_state=SEED)
    sc.pp.neighbors(data, n_neighbors=MAIN_NEIGHBORS, n_pcs=n_pcs, random_state=SEED)
    sc.tl.umap(data, random_state=SEED)
    for resolution in LEIDEN_RESOLUTIONS:
        sc.tl.leiden(
            data,
            resolution=resolution,
            key_added=f"beta_leiden_{resolution:g}",
            random_state=SEED,
            flavor="igraph",
            n_iterations=2,
        )


def _cluster_markers(data: ad.AnnData) -> pd.DataFrame:
    sc.tl.rank_genes_groups(data, groupby=PREFERRED_CLUSTER_KEY, method="wilcoxon", use_raw=False, pts=True)
    frames: list[pd.DataFrame] = []
    for cluster in data.obs[PREFERRED_CLUSTER_KEY].cat.categories:
        frame = sc.get.rank_genes_groups_df(data, group=cluster)
        frame.insert(0, "cluster", str(cluster))
        frame["eligible_for_annotation"] = frame["logfoldchanges"].abs().ge(0.58) & frame["pvals"].lt(0.05)
        frame["test_scope"] = "beta_cluster_marker_annotation_only"
        frame["condition_comparison"] = False
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _expression_vector(data: ad.AnnData, gene: str) -> np.ndarray:
    matrix = data[:, gene].X
    return np.asarray(matrix.toarray()).ravel() if sp.issparse(matrix) else np.asarray(matrix).ravel()


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    finite = np.isfinite(x) & np.isfinite(y)
    if int(finite.sum()) < 3 or np.nanstd(x[finite]) == 0 or np.nanstd(y[finite]) == 0:
        return float("nan")
    return float(spearmanr(x[finite], y[finite]).statistic)


def _root_candidates(data: ad.AnnData, count: int = 3) -> list[int]:
    maturity = data.obs["beta_root_maturity_score"].rank(pct=True).to_numpy(float)
    stress = data.obs["beta_root_stress_score"].rank(pct=True).to_numpy(float)
    root_score = maturity - stress
    eligible = np.flatnonzero(maturity >= 0.75)
    if eligible.size < count:
        eligible = np.arange(data.n_obs)
    order = eligible[np.argsort(root_score[eligible])[::-1]]
    chosen: list[int] = []
    pca = np.asarray(data.obsm["X_pca"][:, : min(10, data.obsm["X_pca"].shape[1])])
    for candidate in order:
        if not chosen:
            chosen.append(int(candidate))
        else:
            distances = [np.linalg.norm(pca[candidate] - pca[item]) for item in chosen]
            if min(distances) >= np.nanmedian(np.linalg.norm(pca - pca[chosen[0]], axis=1)) * 0.05:
                chosen.append(int(candidate))
        if len(chosen) == count:
            break
    for candidate in order:
        if len(chosen) == count:
            break
        if int(candidate) not in chosen:
            chosen.append(int(candidate))
    return chosen


def _graph_metrics(connectivity: sp.spmatrix, libraries: np.ndarray) -> dict[str, float | int]:
    graph = connectivity.tocsr()
    n_components, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels)
    upper = sp.triu(graph, k=1).tocoo()
    cross = float(np.mean(libraries[upper.row] != libraries[upper.col])) if upper.nnz else float("nan")
    return {
        "connected_components": int(n_components),
        "largest_component_fraction": float(sizes.max() / graph.shape[0]),
        "cross_library_edge_fraction": cross,
    }


def _trajectory_variant(data: ad.AnnData, spec: TrajectoryRun) -> tuple[np.ndarray, pd.Series, dict[str, float | int]]:
    n_pcs = min(spec.n_pcs, data.obsm["X_pca"].shape[1])
    temp = ad.AnnData(X=sp.csr_matrix((data.n_obs, 1)), obs=data.obs[["library_id"]].copy())
    temp.obsm["X_pca"] = np.asarray(data.obsm["X_pca"][:, :n_pcs], dtype=np.float32)
    sc.pp.neighbors(temp, n_neighbors=spec.neighbor_count, n_pcs=n_pcs, use_rep="X_pca", random_state=spec.seed)
    sc.tl.leiden(temp, resolution=0.5, key_added="variant_cluster", random_state=spec.seed, flavor="igraph", n_iterations=2)
    sc.tl.diffmap(temp, n_comps=15)
    temp.uns["iroot"] = int(spec.root_index)
    sc.tl.dpt(temp, n_dcs=10)
    pseudotime = temp.obs["dpt_pseudotime"].to_numpy(float)
    metrics = _graph_metrics(temp.obsp["connectivities"], temp.obs["library_id"].astype(str).to_numpy())
    return pseudotime, temp.obs["variant_cluster"].copy(), metrics


def _trajectory_analysis(data: ad.AnnData) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    roots = _root_candidates(data)
    specs = [
        TrajectoryRun("baseline", 30, 40, 42, roots[0]),
        TrajectoryRun("neighbors_15", 15, 40, 42, roots[0]),
        TrajectoryRun("neighbors_50", 50, 40, 42, roots[0]),
        TrajectoryRun("pcs_30", 30, 30, 42, roots[0]),
        TrajectoryRun("seed_137", 30, 40, 137, roots[0]),
        TrajectoryRun("alternate_root_1", 30, 40, 42, roots[1]),
        TrajectoryRun("alternate_root_2", 30, 40, 42, roots[2]),
    ]
    results: dict[str, np.ndarray] = {}
    clusterings: dict[str, pd.Series] = {}
    graph_stats: dict[str, dict[str, float | int]] = {}
    for spec in specs:
        LOGGER.info("trajectory sensitivity: %s", spec.name)
        pt, clustering, metrics = _trajectory_variant(data, spec)
        results[spec.name] = pt
        clusterings[spec.name] = clustering
        graph_stats[spec.name] = metrics

    baseline = results["baseline"]
    baseline_clusters = clusterings["baseline"].astype(str).to_numpy()
    rows: list[dict[str, Any]] = []
    for spec in specs:
        pt = results[spec.name]
        rows.append({
            "variant": spec.name,
            "neighbors": spec.neighbor_count,
            "n_pcs": spec.n_pcs,
            "seed": spec.seed,
            "root_cell_id": str(data.obs_names[spec.root_index]),
            "root_library": str(data.obs.iloc[spec.root_index]["library_id"]),
            "finite_fraction": float(np.isfinite(pt).mean()),
            "spearman_vs_baseline": 1.0 if spec.name == "baseline" else _spearman(baseline, pt),
            "cluster_ari_vs_baseline": 1.0 if spec.name == "baseline" else float(adjusted_rand_score(baseline_clusters, clusterings[spec.name].astype(str).to_numpy())),
            **graph_stats[spec.name],
            "scope": "technical trajectory sensitivity; no condition inference",
        })
    sensitivity = pd.DataFrame(rows)

    data.obs["beta_relative_state_distance"] = baseline
    data.obs["beta_trajectory_root"] = False
    data.obs.iloc[roots[0], data.obs.columns.get_loc("beta_trajectory_root")] = True
    data.uns["iroot"] = roots[0]
    sc.tl.paga(data, groups=PREFERRED_CLUSTER_KEY, model="v1.2")
    data.uns["paga"]["model_used"] = "v1.2"

    finite = np.isfinite(baseline)
    finite_fraction = float(finite.mean())
    maturity_rho = _spearman(baseline, data.obs["beta_direction_maturity_heldout_score"].to_numpy(float))
    stress_rho = _spearman(baseline, data.obs["beta_direction_stress_heldout_score"].to_numpy(float))
    s_rho = _spearman(baseline, data.obs["beta_S_score"].to_numpy(float))
    g2m_rho = _spearman(baseline, data.obs["beta_G2M_score"].to_numpy(float))
    variant_correlations = sensitivity.loc[sensitivity["variant"] != "baseline", "spearman_vs_baseline"].to_numpy(float)

    qbin = pd.qcut(pd.Series(baseline[finite]), q=5, labels=False, duplicates="drop")
    bin_frame = pd.DataFrame({
        "bin": qbin.to_numpy(),
        "library_id": data.obs.loc[finite, "library_id"].astype(str).to_numpy(),
    })
    coverage = bin_frame.groupby(["bin", "library_id"], observed=True).size().rename("n").reset_index()
    coverage["fraction_within_bin"] = coverage["n"] / coverage.groupby("bin")["n"].transform("sum")
    library_count_by_bin = coverage.groupby("bin")["library_id"].nunique()
    all_libraries_each_bin = bool((library_count_by_bin == data.obs["library_id"].nunique()).all())
    minimum_library_fraction = float(coverage["fraction_within_bin"].min()) if all_libraries_each_bin else 0.0

    baseline_graph = graph_stats["baseline"]
    criteria = {
        "largest_component_fraction_ge_0.95": float(baseline_graph["largest_component_fraction"]) >= 0.95,
        "finite_pseudotime_fraction_ge_0.95": finite_fraction >= 0.95,
        "all_sensitivity_correlations_ge_0.60": bool(np.isfinite(variant_correlations).all() and (variant_correlations >= 0.60).all()),
        "median_sensitivity_correlation_ge_0.75": bool(np.isfinite(variant_correlations).all() and np.median(variant_correlations) >= 0.75),
        "both_libraries_in_each_pseudotime_quintile": all_libraries_each_bin,
        "minimum_library_fraction_each_quintile_ge_0.05": minimum_library_fraction >= 0.05,
        "heldout_maturity_decreases_along_distance_rho_le_minus_0.15": bool(np.isfinite(maturity_rho) and maturity_rho <= -0.15),
        "heldout_stress_increases_along_distance_rho_ge_0.15": bool(np.isfinite(stress_rho) and stress_rho >= 0.15),
        "cell_cycle_not_dominant_abs_rho_lt_0.50": bool(max(abs(s_rho), abs(g2m_rho)) < 0.50),
        "paga_contains_intercluster_edge": bool(data.uns["paga"]["connectivities"].nnz > 0),
    }
    failures = [name for name, passed in criteria.items() if not passed]
    passed = not failures
    gate = {
        "status": "PASS_exploratory_relative_state_distance" if passed else "NO-GO_directional_pseudotime",
        "directional_trajectory_eligible": passed,
        "scope": "beta-lineage cross-sectional state geometry only",
        "root_definition": "highest root-maturity anchor rank minus root-stress anchor rank among root-maturity top quartile; condition was not used; root genes are disjoint from direction-diagnostic genes",
        "root_cell_id": str(data.obs_names[roots[0]]),
        "root_library_recorded_for_audit_only": str(data.obs.iloc[roots[0]]["library_id"]),
        "criteria": criteria,
        "failed_criteria": failures,
        "metrics": {
            **baseline_graph,
            "finite_pseudotime_fraction": finite_fraction,
            "minimum_library_fraction_each_quintile": minimum_library_fraction,
            "heldout_maturity_spearman_rho": maturity_rho,
            "heldout_stress_spearman_rho": stress_rho,
            "S_score_spearman_rho": s_rho,
            "G2M_score_spearman_rho": g2m_rho,
            "minimum_sensitivity_spearman": float(np.nanmin(variant_correlations)),
            "median_sensitivity_spearman": float(np.nanmedian(variant_correlations)),
        },
        "interpretation": (
            "DPT is retained only as a hypothesis-generating relative state distance within adult beta cells; it is not developmental time or a PFOS-induced transition."
            if passed else
            "Directional DPT interpretation is rejected. Retain the undirected PAGA state graph and the computed distances only as diagnostics; do not claim development or PFOS-driven transition."
        ),
        "condition_inference": "NO-GO: one library per condition",
        "temporal_or_developmental_trajectory_eligible": False,
        "directional_state_axis_eligible": passed,
    }
    coverage["scope"] = "descriptive library coverage along relative state distance; no condition test"
    return gate, sensitivity, coverage


def _cluster_and_library_tables(data: ad.AnnData) -> dict[str, pd.DataFrame]:
    cluster_counts = data.obs.groupby(["library_id", PREFERRED_CLUSTER_KEY], observed=True).size().rename("n_cells").reset_index()
    cluster_counts["captured_fraction_within_library"] = cluster_counts["n_cells"] / cluster_counts.groupby("library_id", observed=True)["n_cells"].transform("sum")
    cluster_counts["scope"] = "captured beta-cell composition; descriptive only; library is not an animal replicate"

    score_cols = [f"beta_score_{name}" for name in STATE_GENE_SETS] + ["beta_S_score", "beta_G2M_score"]
    summary_rows: list[dict[str, Any]] = []
    for keys, group in data.obs.groupby(["library_id", PREFERRED_CLUSTER_KEY], observed=True):
        for score in score_cols:
            values = group[score].to_numpy(float)
            summary_rows.append({
                "library_id": str(keys[0]), "beta_cluster": str(keys[1]), "score": score,
                "n_cells": len(group), "mean": float(np.mean(values)), "median": float(np.median(values)),
                "q25": float(np.quantile(values, 0.25)), "q75": float(np.quantile(values, 0.75)),
                "scope": "descriptive continuous score summary; no condition test",
            })
    score_summary = pd.DataFrame(summary_rows)

    qc_rows: list[dict[str, Any]] = []
    qc_cols = [column for column in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "doublet_score") if column in data.obs]
    for library, group in data.obs.groupby("library_id", observed=True):
        for qc in qc_cols:
            for score in score_cols:
                qc_rows.append({
                    "library_id": str(library), "qc_metric": qc, "score": score,
                    "n_cells": len(group), "spearman_rho": _spearman(group[qc].to_numpy(float), group[score].to_numpy(float)),
                    "scope": "within-library QC association diagnostic; no P-value and no condition inference",
                })
    qc_associations = pd.DataFrame(qc_rows)

    phase = data.obs.groupby(["library_id", "beta_cell_cycle_phase"], observed=True).size().rename("n_cells").reset_index()
    phase["captured_fraction_within_library"] = phase["n_cells"] / phase.groupby("library_id", observed=True)["n_cells"].transform("sum")
    phase["scope"] = "descriptive score-derived phase; does not establish proliferation"
    return {"cluster_composition": cluster_counts, "score_summary": score_summary, "qc_associations": qc_associations, "cell_cycle": phase}


def _layer_gene_matrix(data: ad.AnnData, layer: str, genes: list[str]) -> sp.csr_matrix:
    positions = data.var_names.get_indexer(genes)
    if np.any(positions < 0):
        raise ValueError(f"missing genes in {layer} identity audit: {[genes[i] for i in np.flatnonzero(positions < 0)]}")
    return sp.csr_matrix(data.layers[layer][:, positions])


def _identity_uncertainty_audit(data: ad.AnnData) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Add neutral beta-identity review flags without changing parent labels.

    Thresholds are fixed relative to all parent-labelled beta cells. Candidate
    flags are screening aids only; none of them reclassifies a cell.
    """
    beta_identity = _present_genes(data, BETA_IDENTITY_GENES, 3, "beta_identity_audit")
    beta_products = _present_genes(data, BETA_PRODUCT_GENES, 2, "beta_product_audit")
    noncognate = _present_genes(data, NONCOGNATE_HORMONES, 3, "noncognate_hormone_audit")
    alternative_tfs = [gene for gene in ALTERNATIVE_LINEAGE_TFS if gene in data.var_names]

    corrected_identity = _layer_gene_matrix(data, "counts", beta_identity)
    corrected_products = _layer_gene_matrix(data, "counts", beta_products)
    corrected_noncognate = _layer_gene_matrix(data, "counts", noncognate)
    raw_noncognate = _layer_gene_matrix(data, "counts_raw", noncognate)
    corrected_alternative = (
        _layer_gene_matrix(data, "counts", alternative_tfs)
        if alternative_tfs else sp.csr_matrix((data.n_obs, 0), dtype=np.float32)
    )
    corrected_total = np.asarray(data.layers["counts"].sum(axis=1)).ravel().astype(float)
    raw_total = np.asarray(data.layers["counts_raw"].sum(axis=1)).ravel().astype(float)
    corrected_noncognate_umi = np.asarray(corrected_noncognate.sum(axis=1)).ravel().astype(float)
    raw_noncognate_umi = np.asarray(raw_noncognate.sum(axis=1)).ravel().astype(float)

    thresholds = {
        "low_complexity_quantile": 0.05,
        "n_genes_by_counts_le": float(np.quantile(data.obs["n_genes_by_counts"].to_numpy(float), 0.05)),
        "total_counts_le": float(np.quantile(data.obs["total_counts"].to_numpy(float), 0.05)),
        "pct_counts_mt_ge": float(np.quantile(data.obs["pct_counts_mt"].to_numpy(float), 0.95)),
        "doublet_score_ge": float(np.quantile(data.obs["doublet_score"].to_numpy(float), 0.95)) if "doublet_score" in data.obs else None,
        "low_beta_identity_detected_genes_le": 1,
        "residual_noncognate_detected_genes_ge": 2,
        "residual_noncognate_corrected_umi_fraction_ge": 0.002,
        "ambient_sensitive_removed_fraction_ge": 0.50,
        "mixed_lineage_requires_alternative_tf_detected_ge": 1,
        "review_cluster": IDENTITY_REVIEW_CLUSTER,
        "scope": "candidate annotation/QC review only; no parent-label replacement or condition inference",
    }
    beta_identity_detected = np.asarray((corrected_identity > 0).sum(axis=1)).ravel().astype(int)
    beta_products_detected = np.asarray((corrected_products > 0).sum(axis=1)).ravel().astype(int)
    noncognate_detected = np.asarray((corrected_noncognate > 0).sum(axis=1)).ravel().astype(int)
    alternative_detected = np.asarray((corrected_alternative > 0).sum(axis=1)).ravel().astype(int)
    corrected_fraction = np.divide(
        corrected_noncognate_umi,
        corrected_total,
        out=np.zeros(data.n_obs, dtype=float),
        where=corrected_total > 0,
    )
    raw_fraction = np.divide(raw_noncognate_umi, raw_total, out=np.zeros(data.n_obs, dtype=float), where=raw_total > 0)
    removed_fraction = np.divide(
        raw_noncognate_umi - corrected_noncognate_umi,
        raw_noncognate_umi,
        out=np.zeros(data.n_obs, dtype=float),
        where=raw_noncognate_umi > 0,
    )
    removed_fraction = np.clip(removed_fraction, 0, 1)

    low_complexity = (
        data.obs["n_genes_by_counts"].to_numpy(float) <= thresholds["n_genes_by_counts_le"]
    ) | (data.obs["total_counts"].to_numpy(float) <= thresholds["total_counts_le"])
    high_mt = data.obs["pct_counts_mt"].to_numpy(float) >= thresholds["pct_counts_mt_ge"]
    high_doublet_score = (
        data.obs["doublet_score"].to_numpy(float) >= thresholds["doublet_score_ge"]
        if thresholds["doublet_score_ge"] is not None else np.zeros(data.n_obs, dtype=bool)
    )
    predicted_doublet = (
        data.obs["predicted_doublet"].astype(bool).to_numpy()
        if "predicted_doublet" in data.obs else np.zeros(data.n_obs, dtype=bool)
    )
    low_identity = beta_identity_detected <= thresholds["low_beta_identity_detected_genes_le"]
    residual_noncognate = (
        noncognate_detected >= thresholds["residual_noncognate_detected_genes_ge"]
    ) & (corrected_fraction >= thresholds["residual_noncognate_corrected_umi_fraction_ge"])
    ambient_sensitive = (raw_noncognate_umi > 0) & (
        removed_fraction >= thresholds["ambient_sensitive_removed_fraction_ge"]
    )
    mixed_lineage = residual_noncognate & (alternative_detected >= 1) & (beta_identity_detected >= 2)
    annotation_uncertain = low_complexity | low_identity | residual_noncognate | mixed_lineage | predicted_doublet

    review_status = np.full(data.n_obs, "not_flagged", dtype=object)
    review_status[high_doublet_score & annotation_uncertain] = "high_doublet_score_review"
    review_status[ambient_sensitive & residual_noncognate] = "ambient_sensitive_residual_hormone_review"
    review_status[residual_noncognate] = "residual_noncognate_hormone_review"
    review_status[low_identity] = "low_beta_identity_support_review"
    review_status[low_complexity & low_identity] = "low_complexity_low_beta_identity_review"
    review_status[low_complexity & residual_noncognate] = "low_complexity_residual_hormone_uncertain"
    review_status[mixed_lineage] = "mixed_lineage_candidate"

    audit_values: dict[str, Any] = {
        "beta_identity_genes_detected": beta_identity_detected,
        "beta_product_genes_detected": beta_products_detected,
        "beta_noncognate_hormones_detected": noncognate_detected,
        "beta_alternative_lineage_tfs_detected": alternative_detected,
        "beta_noncognate_umi_raw": raw_noncognate_umi,
        "beta_noncognate_umi_corrected": corrected_noncognate_umi,
        "beta_noncognate_fraction_raw": raw_fraction,
        "beta_noncognate_fraction_corrected": corrected_fraction,
        "beta_noncognate_removed_fraction": removed_fraction,
        "beta_low_complexity_candidate": low_complexity,
        "beta_high_mt_candidate": high_mt,
        "beta_high_doublet_score_candidate": high_doublet_score,
        "beta_low_identity_candidate": low_identity,
        "beta_residual_noncognate_candidate": residual_noncognate,
        "beta_ambient_sensitive_candidate": ambient_sensitive,
        "beta_mixed_lineage_candidate": mixed_lineage,
        "beta_annotation_uncertainty_candidate": annotation_uncertain,
        "beta_annotation_review_candidate": review_status,
    }
    for column, values in audit_values.items():
        data.obs[column] = values

    selected_obs = [
        "library_id", "condition", PREFERRED_CLUSTER_KEY, "cell_type", "total_counts",
        "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "doublet_score",
        "predicted_doublet", "annotation_confidence", "annotation_review_status",
        "mia_label", "mia_support", "singler_label", "lineage_conflict_review",
        *audit_values,
    ]
    selected_obs = list(dict.fromkeys(column for column in selected_obs if column in data.obs))
    cell_audit = data.obs[selected_obs].copy().reset_index(names="cell_id")
    cell_audit.rename(columns={PREFERRED_CLUSTER_KEY: "beta_cluster", "cell_type": "parent_cell_type"}, inplace=True)
    cell_audit["scope"] = "neutral candidate flag; parent beta label retained; no condition inference"

    gene_rows: list[dict[str, Any]] = []
    gene_groups = {
        "beta_identity_tf_or_maturity": beta_identity,
        "beta_product": beta_products,
        "noncognate_hormone": noncognate,
        "alternative_lineage_tf": alternative_tfs,
    }
    clusters = data.obs[PREFERRED_CLUSTER_KEY].astype(str).to_numpy()
    for layer in ("counts_raw", "counts"):
        for evidence_group, genes in gene_groups.items():
            if not genes:
                continue
            matrix = _layer_gene_matrix(data, layer, genes)
            for cluster in sorted(np.unique(clusters)):
                mask = clusters == cluster
                subset = matrix[mask]
                means = np.asarray(subset.mean(axis=0)).ravel()
                detections = np.asarray((subset > 0).mean(axis=0)).ravel()
                for gene, mean_umi, detection in zip(genes, means, detections, strict=True):
                    gene_rows.append({
                        "beta_cluster": cluster, "layer": layer, "evidence_group": evidence_group,
                        "gene": gene, "n_cells": int(mask.sum()), "mean_umi": float(mean_umi),
                        "detection_fraction": float(detection),
                        "scope": "cluster annotation/QC evidence; no condition comparison",
                    })
    gene_summary = pd.DataFrame(gene_rows)

    cluster_rows: list[dict[str, Any]] = []
    for cluster, group in cell_audit.groupby("beta_cluster", observed=True):
        residual_rate = float(group["beta_residual_noncognate_candidate"].mean())
        low_complexity_rate = float(group["beta_low_complexity_candidate"].mean())
        low_identity_rate = float(group["beta_low_identity_candidate"].mean())
        mixed_rate = float(group["beta_mixed_lineage_candidate"].mean())
        ambient_rate = float(group["beta_ambient_sensitive_candidate"].mean())
        if low_complexity_rate >= 0.50 and residual_rate >= 0.50:
            status = "low_complexity_with_residual_hormones_beta_identity_uncertain"
        elif mixed_rate >= 0.10:
            status = "mixed_lineage_candidate"
        elif residual_rate >= 0.50:
            status = "residual_noncognate_hormones_need_review"
        elif low_identity_rate >= 0.50:
            status = "low_beta_identity_support_need_review"
        else:
            status = "no_cluster_level_identity_flag"
        cluster_rows.append({
            "beta_cluster": str(cluster), "n_cells": len(group),
            "median_total_counts": float(group["total_counts"].median()),
            "median_n_genes_by_counts": float(group["n_genes_by_counts"].median()),
            "median_pct_counts_mt": float(group["pct_counts_mt"].median()),
            "median_doublet_score": float(group["doublet_score"].median()) if "doublet_score" in group else np.nan,
            "predicted_doublet_fraction": float(group["predicted_doublet"].astype(bool).mean()) if "predicted_doublet" in group else np.nan,
            "low_complexity_candidate_fraction": low_complexity_rate,
            "high_mt_candidate_fraction": float(group["beta_high_mt_candidate"].mean()),
            "high_doublet_score_candidate_fraction": float(group["beta_high_doublet_score_candidate"].mean()),
            "low_beta_identity_candidate_fraction": low_identity_rate,
            "residual_noncognate_candidate_fraction": residual_rate,
            "ambient_sensitive_candidate_fraction": ambient_rate,
            "mixed_lineage_candidate_fraction": mixed_rate,
            "median_beta_identity_genes_detected": float(group["beta_identity_genes_detected"].median()),
            "median_raw_noncognate_umi": float(group["beta_noncognate_umi_raw"].median()),
            "median_corrected_noncognate_umi": float(group["beta_noncognate_umi_corrected"].median()),
            "median_noncognate_removed_fraction": float(group["beta_noncognate_removed_fraction"].median()),
            "cluster_review_status": status,
            "parent_label_action": "retained_without_overwrite",
            "scope": "cluster-level annotation/QC review; no condition inference",
        })
    cluster_summary = pd.DataFrame(cluster_rows)
    return cell_audit, gene_summary, cluster_summary, thresholds


def _cluster4_exclusion_sensitivity(data: ad.AnnData) -> pd.DataFrame:
    """Describe how summaries change when review cluster 4 is omitted.

    The full 14,304-cell analysis remains primary; this is not a second QC rule.
    """
    cluster = data.obs[PREFERRED_CLUSTER_KEY].astype(str).to_numpy()
    score_columns = [f"beta_score_{name}" for name in STATE_GENE_SETS] + ["beta_S_score", "beta_G2M_score"]
    rows: list[dict[str, Any]] = []
    library_values = [str(value) for value in data.obs["library_id"].astype(str).unique()] + ["all_libraries_pooled"]
    target_expression = {gene: _expression_vector(data, gene) for gene in ("Mpc1", "Mpc2", "Atg7") if gene in data.var_names}
    for library in library_values:
        library_mask = np.ones(data.n_obs, dtype=bool) if library == "all_libraries_pooled" else data.obs["library_id"].astype(str).to_numpy() == library
        core_mask = library_mask & (cluster != IDENTITY_REVIEW_CLUSTER)
        for score in score_columns:
            values = data.obs[score].to_numpy(float)
            all_value = float(values[library_mask].mean())
            core_value = float(values[core_mask].mean()) if core_mask.any() else np.nan
            rows.append({
                "library_id": library, "metric": f"mean_{score}", "all_beta_value": all_value,
                "exclude_cluster4_value": core_value, "difference_exclude_minus_all": core_value - all_value,
                "all_beta_cells": int(library_mask.sum()), "exclude_cluster4_cells": int(core_mask.sum()),
            })
        for gene, values in target_expression.items():
            for statistic, function in (("mean_log1p", np.mean), ("detection_fraction", lambda x: np.mean(x > 0))):
                all_value = float(function(values[library_mask]))
                core_value = float(function(values[core_mask])) if core_mask.any() else np.nan
                rows.append({
                    "library_id": library, "metric": f"{gene}_{statistic}", "all_beta_value": all_value,
                    "exclude_cluster4_value": core_value, "difference_exclude_minus_all": core_value - all_value,
                    "all_beta_cells": int(library_mask.sum()), "exclude_cluster4_cells": int(core_mask.sum()),
                })
    result = pd.DataFrame(rows)
    result["review_cluster"] = IDENTITY_REVIEW_CLUSTER
    result["primary_analysis"] = "all parent-labelled beta cells"
    result["scope"] = "descriptive technical sensitivity only; exclusion is not a revised QC rule or condition test"
    return result


def _balanced_subsampling(data: ad.AnnData, repeats: int = 20, max_cells: int = 3000) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    score_cols = [f"beta_score_{name}" for name in STATE_GENE_SETS] + ["beta_S_score", "beta_G2M_score"]
    library_indices = {str(lib): np.flatnonzero(data.obs["library_id"].astype(str).to_numpy() == str(lib)) for lib in data.obs["library_id"].unique()}
    draw_n = min(max_cells, *(len(index) for index in library_indices.values()))
    rows: list[dict[str, Any]] = []
    for repeat in range(repeats):
        for library, indices in library_indices.items():
            chosen = rng.choice(indices, size=draw_n, replace=False)
            for score in score_cols:
                full = data.obs.iloc[indices][score].to_numpy(float)
                sampled = data.obs.iloc[chosen][score].to_numpy(float)
                rows.append({
                    "repeat": repeat + 1, "library_id": library, "score": score, "cells_sampled": draw_n,
                    "full_mean": float(full.mean()), "subsample_mean": float(sampled.mean()),
                    "mean_difference": float(sampled.mean() - full.mean()),
                    "full_median": float(np.median(full)), "subsample_median": float(np.median(sampled)),
                    "scope": "equal-cell technical sensitivity only; not a biological replicate analysis",
                })
    return pd.DataFrame(rows)


def _theme(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7, width=0.5)
    ax.title.set_fontsize(8)
    ax.xaxis.label.set_size(8)
    ax.yaxis.label.set_size(8)


def _save_with_source(output_dir: Path, stem: str, source: pd.DataFrame, draw: Callable[[], plt.Figure]) -> dict[str, str]:
    source_path = output_dir / "source_data" / f"{stem}.tsv"
    atomic_tsv(source_path, source)
    figure = draw()
    save_figure(figure, output_dir / "figures" / stem, dpi=300)
    return {
        "png": f"figures/{stem}.png", "pdf": f"figures/{stem}.pdf", "source": f"source_data/{stem}.tsv",
    }


def _umap_source(data: ad.AnnData, columns: list[str]) -> pd.DataFrame:
    frame = pd.DataFrame({"cell_id": data.obs_names.astype(str), "UMAP1": data.obsm["X_umap"][:, 0], "UMAP2": data.obsm["X_umap"][:, 1]})
    for column in columns:
        frame[column] = data.obs[column].astype(str).to_numpy() if not pd.api.types.is_numeric_dtype(data.obs[column]) else data.obs[column].to_numpy()
    return frame


def _draw_umap_categorical(source: pd.DataFrame, column: str, title: str, colors: dict[str, str] | None = None) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    categories = [key for key in (colors or {}).keys() if key in set(source[column].astype(str))] if colors else sorted(source[column].astype(str).unique())
    for index, category in enumerate(categories):
        sub = source[source[column].astype(str).eq(category)]
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=3, alpha=0.72, linewidths=0, rasterized=True, label=category, color=colors[category] if colors is not None else category_color(category))
    ax.set(title=title, xlabel="UMAP1", ylabel="UMAP2")
    ax.set_xticks([]); ax.set_yticks([])
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, markerscale=3)
    _theme(ax)
    return fig


def _draw_continuous_grid(
    source: pd.DataFrame,
    columns: list[str],
    title: str,
    ncols: int = 3,
) -> plt.Figure:
    if ncols < 1:
        raise ValueError("ncols must be positive")
    nrows = int(np.ceil(len(columns) / ncols))
    width = 14.0 if ncols == 3 else 4.7 * ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(width, 4.2 * nrows), squeeze=False)
    for ax, column in zip(axes.flat, columns, strict=False):
        values = source[column].to_numpy(float)
        order = np.argsort(values)
        scatter = ax.scatter(source["UMAP1"].to_numpy()[order], source["UMAP2"].to_numpy()[order], c=values[order], s=3, cmap="viridis", linewidths=0, rasterized=True)
        ax.set_title(column.replace("beta_score_", "").replace("beta_", ""))
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(scatter, ax=ax, shrink=0.72)
    for ax in axes.flat[len(columns):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    return fig


def _marker_summary(data: ad.AnnData, markers: pd.DataFrame) -> pd.DataFrame:
    top = markers.loc[markers["eligible_for_annotation"] & ~markers["names"].astype(str).str.match(r"^(mt-|Rp[sl]|Ins[12]$)")].copy()
    genes: list[str] = []
    for cluster in sorted(top["cluster"].astype(str).unique()):
        genes.extend(top.loc[top["cluster"].astype(str).eq(cluster)].sort_values("scores", ascending=False)["names"].astype(str).head(3))
    anchors = [gene for values in STATE_GENE_SETS.values() for gene in values if gene in data.var_names]
    genes = list(dict.fromkeys(genes + anchors))[:45]
    rows: list[dict[str, Any]] = []
    clusters = data.obs[PREFERRED_CLUSTER_KEY].astype(str)
    for gene in genes:
        expression = _expression_vector(data, gene)
        for cluster in sorted(clusters.unique()):
            values = expression[clusters.to_numpy() == cluster]
            rows.append({"cluster": cluster, "gene": gene, "mean_log1p": float(values.mean()), "detection_fraction": float((values > 0).mean()), "n_cells": len(values)})
    frame = pd.DataFrame(rows)
    frame["zscore_cluster_mean"] = frame.groupby("gene")["mean_log1p"].transform(lambda value: (value - value.mean()) / (value.std(ddof=0) if value.std(ddof=0) else 1.0))
    return frame


def _draw_marker_dot(source: pd.DataFrame) -> plt.Figure:
    genes = list(dict.fromkeys(source["gene"]))
    clusters = sorted(source["cluster"].astype(str).unique())
    fig, ax = plt.subplots(figsize=(max(10.0, len(genes) * 0.30), max(4.5, len(clusters) * 0.55)))
    xmap = {gene: index for index, gene in enumerate(genes)}
    ymap = {cluster: index for index, cluster in enumerate(clusters)}
    points = ax.scatter(source["gene"].map(xmap), source["cluster"].astype(str).map(ymap), s=200 * source["detection_fraction"], c=source["zscore_cluster_mean"], cmap=DIVERGING_CMAP, vmin=-2.5, vmax=2.5, linewidths=0.2, edgecolors="#444444")
    ax.set_xticks(range(len(genes)), genes, rotation=90)
    ax.set_yticks(range(len(clusters)), clusters)
    ax.set(xlabel="Gene", ylabel="Beta cluster", title="Beta-cluster marker and program evidence")
    fig.colorbar(points, ax=ax, label="Cluster mean z-score", shrink=0.8)
    _theme(ax); fig.tight_layout()
    return fig


def _draw_marker_heatmap(source: pd.DataFrame) -> plt.Figure:
    pivot = source.pivot(index="gene", columns="cluster", values="zscore_cluster_mean")
    fig, ax = plt.subplots(figsize=(max(6.5, pivot.shape[1] * 0.65), max(8.0, pivot.shape[0] * 0.22)))
    sns.heatmap(pivot, cmap=DIVERGING_CMAP, center=0, vmin=-2.5, vmax=2.5, ax=ax, cbar_kws={"label": "Cluster mean z-score"})
    ax.set(title="Beta-cluster marker heatmap", xlabel="Beta cluster", ylabel="Gene")
    fig.tight_layout()
    return fig


def _draw_target_violin(source: pd.DataFrame) -> plt.Figure:
    melted = source.melt(id_vars=["cell_id", "library_id", PREFERRED_CLUSTER_KEY], value_vars=["Mpc1", "Mpc2", "Atg7"], var_name="gene", value_name="log1p_expression")
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.5), sharey=False)
    for gene, ax in zip(("Mpc1", "Mpc2", "Atg7"), axes, strict=True):
        subset = melted[melted["gene"].eq(gene)]
        sns.violinplot(data=subset, x="library_id", y="log1p_expression", hue="library_id", palette={"PFOS_Lib01": CONDITION_COLORS["PFOS"], "WT_Lib01": CONDITION_COLORS["WT"]}, cut=0, inner=None, linewidth=0, alpha=0.30, legend=False, ax=ax)
        for position, library in enumerate(subset["library_id"].drop_duplicates()):
            values = subset.loc[subset["library_id"].eq(library), "log1p_expression"].to_numpy()
            values = values[np.isfinite(values)]
            if len(values):
                q1, median, q3 = np.quantile(values, [.25, .5, .75], method="linear")
                ax.vlines(position, q1, q3, color=PROFILE["palettes"]["neutral"]["dark"], linewidth=.25 * 72 / 25.4, zorder=4)
                ax.scatter(position, median, s=(1.1 * 72 / 25.4)**2, color=PROFILE["palettes"]["neutral"]["dark"], linewidths=0, zorder=5)
        ax.set(title=gene, xlabel="Library", ylabel="Normalize-total log1p expression")
        ax.tick_params(axis="x", rotation=20)
        _theme(ax)
    fig.suptitle("Focal-gene expression in beta cells (descriptive)", fontsize=10)
    fig.tight_layout()
    return fig


def _draw_composition(
    source: pd.DataFrame,
    category_column: str = PREFERRED_CLUSTER_KEY,
    title: str = "Captured beta-cluster composition",
    legend_title: str = "Beta cluster",
    colors: dict[str, str] | None = None,
) -> plt.Figure:
    pivot = source.pivot(index="library_id", columns=category_column, values="captured_fraction_within_library").fillna(0)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    plot_colors = [colors[str(category)] for category in pivot.columns] if colors is not None else None
    pivot.plot.bar(stacked=True, ax=ax, color=plot_colors, width=0.72)
    ax.set(title=title, xlabel="Library", ylabel="Captured fraction", ylim=(0, 1))
    ax.legend(title=legend_title, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    _theme(ax); fig.tight_layout()
    return fig


def _draw_identity_audit(
    cell_audit: pd.DataFrame,
    gene_summary: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    thresholds: dict[str, Any],
) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 10.5))
    corrected = gene_summary.loc[
        gene_summary["layer"].eq("counts")
        & gene_summary["evidence_group"].isin(["beta_identity_tf_or_maturity", "noncognate_hormone"])
    ].copy()
    genes = [gene for gene in (*BETA_IDENTITY_GENES, *NONCOGNATE_HORMONES) if gene in set(corrected["gene"])]
    clusters = sorted(corrected["beta_cluster"].astype(str).unique())
    xmap = {gene: index for index, gene in enumerate(genes)}
    ymap = {cluster: index for index, cluster in enumerate(clusters)}
    points = axes[0, 0].scatter(
        corrected["gene"].map(xmap), corrected["beta_cluster"].astype(str).map(ymap),
        s=200 * corrected["detection_fraction"], c=np.log1p(corrected["mean_umi"]),
        cmap="viridis", linewidths=0.2, edgecolors="#444444",
    )
    axes[0, 0].set_xticks(range(len(genes)), genes, rotation=90)
    axes[0, 0].set_yticks(range(len(clusters)), clusters)
    axes[0, 0].set(title="Corrected beta-identity and noncognate-hormone evidence", xlabel="Gene", ylabel="Beta cluster")
    fig.colorbar(points, ax=axes[0, 0], label="log1p mean corrected UMI", shrink=0.78)

    hormone = gene_summary.loc[gene_summary["evidence_group"].eq("noncognate_hormone")].groupby(
        ["beta_cluster", "layer"], observed=True
    )["mean_umi"].sum().rename("mean_summed_noncognate_umi").reset_index()
    hormone_pivot = hormone.pivot(index="beta_cluster", columns="layer", values="mean_summed_noncognate_umi").fillna(0)
    hormone_pivot.plot.bar(ax=axes[0, 1], color={"counts_raw": PROFILE["palettes"]["neutral"]["mid"], "counts": PROFILE["palettes"]["accent"]["warm"]})
    axes[0, 1].set(title="Raw versus corrected noncognate hormone burden", xlabel="Beta cluster", ylabel="Mean summed UMI")
    axes[0, 1].legend(title="Layer", frameon=False)

    review_cluster = cell_audit["beta_cluster"].astype(str).eq(IDENTITY_REVIEW_CLUSTER)
    axes[1, 0].scatter(
        np.log1p(cell_audit.loc[~review_cluster, "total_counts"]),
        cell_audit.loc[~review_cluster, "n_genes_by_counts"], s=2, alpha=0.18,
        c="#999999", rasterized=True, label="Other beta clusters",
    )
    axes[1, 0].scatter(
        np.log1p(cell_audit.loc[review_cluster, "total_counts"]),
        cell_audit.loc[review_cluster, "n_genes_by_counts"], s=3, alpha=0.42,
        c=BETA_CLUSTER_COLORS[IDENTITY_REVIEW_CLUSTER], rasterized=True, label=f"Cluster {IDENTITY_REVIEW_CLUSTER}",
    )
    axes[1, 0].axvline(np.log1p(thresholds["total_counts_le"]), color="#222222", linestyle="--", linewidth=0.8)
    axes[1, 0].axhline(thresholds["n_genes_by_counts_le"], color="#222222", linestyle="--", linewidth=0.8)
    axes[1, 0].set(title="Complexity audit", xlabel="log1p total counts", ylabel="Genes detected")
    axes[1, 0].legend(frameon=False, markerscale=3)

    flag_columns = {
        "low complexity": "low_complexity_candidate_fraction",
        "low beta identity": "low_beta_identity_candidate_fraction",
        "residual noncognate": "residual_noncognate_candidate_fraction",
        "ambient sensitive": "ambient_sensitive_candidate_fraction",
        "mixed lineage": "mixed_lineage_candidate_fraction",
        "high doublet score": "high_doublet_score_candidate_fraction",
    }
    flag_long = cluster_summary.melt(
        id_vars="beta_cluster", value_vars=list(flag_columns.values()),
        var_name="flag", value_name="fraction",
    )
    reverse_names = {value: key for key, value in flag_columns.items()}
    flag_long["flag"] = flag_long["flag"].map(reverse_names)
    sns.barplot(data=flag_long, x="beta_cluster", y="fraction", hue="flag", hue_order=list(flag_columns), palette=dict(zip(flag_columns, PROFILE["palettes"]["accent"].values())), ax=axes[1, 1])
    axes[1, 1].set(title="Neutral annotation-review candidate rates", xlabel="Beta cluster", ylabel="Fraction", ylim=(0, 1))
    axes[1, 1].legend(title="Candidate flag", frameon=False, fontsize=7, ncol=2)
    for ax in axes.flat:
        _theme(ax)
    fig.tight_layout()
    return fig


def _draw_cluster4_sensitivity(source: pd.DataFrame) -> plt.Figure:
    selected_metrics = [
        "mean_beta_score_maturity", "mean_beta_score_gsis", "mean_beta_score_er_stress",
        "mean_beta_score_oxidative_stress", "mean_beta_score_dedifferentiation",
        "Mpc1_mean_log1p", "Mpc2_mean_log1p", "Atg7_mean_log1p",
    ]
    subset = source.loc[
        source["metric"].isin(selected_metrics) & source["library_id"].ne("all_libraries_pooled")
    ].copy()
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    sns.barplot(
        data=subset, y="metric", x="difference_exclude_minus_all", hue="library_id",
        palette={"PFOS_Lib01": CONDITION_COLORS["PFOS"], "WT_Lib01": CONDITION_COLORS["WT"]}, ax=ax,
    )
    ax.axvline(0, color="#444444", linewidth=0.8)
    ax.set(
        title=f"Descriptive sensitivity after excluding review cluster {IDENTITY_REVIEW_CLUSTER}",
        xlabel="Exclude-cluster-4 value minus all-beta value", ylabel="Metric",
    )
    ax.legend(title="Library", frameon=False)
    _theme(ax); fig.tight_layout()
    return fig


def _draw_qc_heatmap(source: pd.DataFrame) -> plt.Figure:
    pivot = source.pivot(index=["library_id", "qc_metric"], columns="score", values="spearman_rho")
    fig, ax = plt.subplots(figsize=(10.5, max(4.5, pivot.shape[0] * 0.42)))
    sns.heatmap(pivot, cmap=DIVERGING_CMAP, center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", ax=ax, cbar_kws={"label": "Spearman rho"})
    ax.set(title="Within-library QC associations with beta-state scores", xlabel="Continuous score", ylabel="Library / QC metric")
    fig.tight_layout()
    return fig


def _paga_source(data: ad.AnnData) -> pd.DataFrame:
    categories = list(data.obs[PREFERRED_CLUSTER_KEY].cat.categories.astype(str))
    coords = pd.DataFrame(data.obsm["X_umap"], columns=["x", "y"]).assign(cluster=data.obs[PREFERRED_CLUSTER_KEY].astype(str).to_numpy()).groupby("cluster").mean()
    counts = data.obs[PREFERRED_CLUSTER_KEY].astype(str).value_counts()
    matrix = data.uns["paga"]["connectivities"].toarray()
    rows: list[dict[str, Any]] = []
    for cluster in categories:
        rows.append({"row_type": "node", "source_cluster": cluster, "target_cluster": "", "weight": np.nan, "x_source": coords.loc[cluster, "x"], "y_source": coords.loc[cluster, "y"], "x_target": np.nan, "y_target": np.nan, "n_cells": int(counts[cluster])})
    for i, source in enumerate(categories):
        for j in range(i + 1, len(categories)):
            weight = float(matrix[i, j])
            if weight >= 0.05:
                target = categories[j]
                rows.append({"row_type": "edge", "source_cluster": source, "target_cluster": target, "weight": weight, "x_source": coords.loc[source, "x"], "y_source": coords.loc[source, "y"], "x_target": coords.loc[target, "x"], "y_target": coords.loc[target, "y"], "n_cells": np.nan})
    frame = pd.DataFrame(rows)
    frame["interpretation"] = "undirected PAGA connectivity; edge threshold 0.05; no temporal direction"
    return frame


def _draw_paga(source: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    edges = source[source["row_type"].eq("edge")]
    for row in edges.itertuples():
        ax.plot([row.x_source, row.x_target], [row.y_source, row.y_target], color="#888888", linewidth=0.4 + 4 * row.weight, alpha=0.8, zorder=1)
    nodes = source[source["row_type"].eq("node")]
    size = 30 + 500 * np.sqrt(nodes["n_cells"] / nodes["n_cells"].max())
    node_colors = [BETA_CLUSTER_COLORS[str(cluster)] for cluster in nodes["source_cluster"]]
    ax.scatter(nodes["x_source"], nodes["y_source"], s=size, c=node_colors, edgecolors="white", linewidths=0.7, zorder=2)
    for row in nodes.itertuples():
        ax.text(row.x_source, row.y_source, row.source_cluster, ha="center", va="center", fontsize=8, zorder=3)
    ax.set_title("Undirected beta-state PAGA graph")
    ax.axis("off"); fig.tight_layout()
    return fig


def _draw_sensitivity(source: pd.DataFrame) -> plt.Figure:
    subset = source[source["variant"].ne("baseline")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    axes[0].barh(subset["variant"], subset["spearman_vs_baseline"], color="#56B4E9")
    axes[0].axvline(0.60, color="#D55E00", linestyle="--", linewidth=1)
    axes[0].set(xlabel="Spearman rho vs baseline distance", title="Distance stability", xlim=(-1, 1))
    axes[1].barh(subset["variant"], subset["cluster_ari_vs_baseline"], color="#009E73")
    axes[1].set(xlabel="Adjusted Rand index", title="Cluster stability", xlim=(0, 1))
    for ax in axes: _theme(ax)
    fig.tight_layout()
    return fig


def _draw_balanced(source: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    sns.boxplot(data=source, x="score", y="mean_difference", hue="library_id", palette={"PFOS_Lib01": CONDITION_COLORS["PFOS"], "WT_Lib01": CONDITION_COLORS["WT"]}, showfliers=False, ax=ax)
    ax.axhline(0, color="#555555", linewidth=0.8)
    ax.set(title="Equal-cell subsampling sensitivity", xlabel="Continuous score", ylabel="Subsample mean - full-library mean")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(title="Library", frameon=False)
    _theme(ax); fig.tight_layout()
    return fig


def _write_plots(
    data: ad.AnnData,
    output_dir: Path,
    markers: pd.DataFrame,
    tables: dict[str, pd.DataFrame],
    sensitivity: pd.DataFrame,
    gate: dict[str, Any],
    balanced: pd.DataFrame,
    identity_cells: pd.DataFrame,
    identity_genes: pd.DataFrame,
    identity_clusters: pd.DataFrame,
    identity_thresholds: dict[str, Any],
    cluster4_sensitivity: pd.DataFrame,
) -> list[dict[str, str]]:
    figures: list[dict[str, str]] = []

    def add(stem: str, title: str, interpretation: str, source: pd.DataFrame, draw: Callable[[], plt.Figure]) -> None:
        files = _save_with_source(output_dir, stem, source, draw)
        figures.append({"title": title, **files, "interpretation": interpretation})

    source = _umap_source(data, ["library_id"])
    add("beta_umap_library", "Beta-cell UMAP by library", "No integration was applied; library and condition are completely confounded.", source, lambda: _draw_umap_categorical(source, "library_id", "Beta cells by library", {"PFOS_Lib01": CONDITION_COLORS["PFOS"], "WT_Lib01": CONDITION_COLORS["WT"]}))
    source = _umap_source(data, [PREFERRED_CLUSTER_KEY])
    add("beta_umap_clusters", "Beta-cell UMAP by cluster", "Leiden 0.5 is an expression-state partition used for description and marker annotation.", source, lambda: _draw_umap_categorical(source, PREFERRED_CLUSTER_KEY, "Beta-state clusters (Leiden 0.5)", BETA_CLUSTER_COLORS))

    score_columns = [f"beta_score_{name}" for name in STATE_GENE_SETS]
    source = _umap_source(data, score_columns)
    add("beta_umap_state_scores", "Continuous beta-state scores", "Scores use expression-matched control genes; cells were not forced into a highest-score class.", source, lambda: _draw_continuous_grid(source, score_columns, "Continuous beta-state program scores"))

    marker_source = _marker_summary(data, markers)
    add("beta_marker_dotplot", "Beta-cluster marker dot plot", "Point size is detection fraction and color is cluster-standardized mean expression.", marker_source, lambda: _draw_marker_dot(marker_source))
    add("beta_marker_heatmap", "Beta-cluster marker heat map", "Displayed genes are annotation-eligible top markers plus frozen program anchors; the full marker table is retained.", marker_source, lambda: _draw_marker_heatmap(marker_source))

    target_source = _umap_source(data, ["library_id", PREFERRED_CLUSTER_KEY])
    for gene in ("Mpc1", "Mpc2", "Atg7"):
        target_source[gene] = _expression_vector(data, gene) if gene in data.var_names else np.nan
    add("beta_target_feature_umap", "Mpc1/Mpc2/Atg7 localization", "Expression is normalize-total log1p and supports beta-cell localization only.", target_source, lambda: _draw_continuous_grid(target_source, ["Mpc1", "Mpc2", "Atg7"], "Focal-gene localization in beta cells"))
    add("beta_target_violin", "Mpc1/Mpc2/Atg7 library distributions", "Library-level distributions are descriptive because each condition contains one library.", target_source, lambda: _draw_target_violin(target_source))

    cycle_source = _umap_source(data, ["beta_cell_cycle_phase", "beta_S_score", "beta_G2M_score"])
    add("beta_cell_cycle_umap", "Beta-cell cycle scores", "Score-derived phase is descriptive; low cycling evidence is not proof of proliferation or arrest.", cycle_source, lambda: _draw_continuous_grid(cycle_source, ["beta_S_score", "beta_G2M_score"], "Cell-cycle program scores"))
    add("beta_cell_cycle_phase_umap", "Beta-cell score-derived phase", "G1/S/G2M labels are score-derived descriptions and do not establish proliferation.", cycle_source, lambda: _draw_umap_categorical(cycle_source, "beta_cell_cycle_phase", "Score-derived cell-cycle phase", CELL_CYCLE_COLORS))
    phase_source = tables["cell_cycle"]
    add("beta_cell_cycle_composition", "Captured cell-cycle phase composition", "Fractions describe score-derived phases within each library and have no condition P-value.", phase_source, lambda: _draw_composition(phase_source, "beta_cell_cycle_phase", "Captured score-derived cell-cycle phase composition", "Cell-cycle phase", CELL_CYCLE_COLORS))

    composition = tables["cluster_composition"]
    add("beta_cluster_composition", "Captured beta-cluster composition", "Fractions are captured-cell descriptions and do not estimate animal-level composition effects.", composition, lambda: _draw_composition(composition, colors=BETA_CLUSTER_COLORS))
    qc_source = tables["qc_associations"]
    add("beta_qc_score_associations", "QC-score association audit", "Spearman coefficients are within-library diagnostics without P-values.", qc_source, lambda: _draw_qc_heatmap(qc_source))

    audit_source = pd.concat([
        identity_cells.assign(row_kind="cell_qc"),
        identity_genes.assign(row_kind="gene_layer_summary"),
        identity_clusters.assign(row_kind="cluster_review_summary"),
    ], ignore_index=True, sort=False)
    add(
        "beta_identity_hormone_qc_audit", "Beta-identity and noncognate-hormone QC audit",
        "Parent beta labels remain unchanged. Low complexity, residual noncognate hormones, ambient-correction sensitivity, alternative-lineage TFs and doublet metrics are neutral review evidence.",
        audit_source,
        lambda: _draw_identity_audit(identity_cells, identity_genes, identity_clusters, identity_thresholds),
    )
    add(
        "beta_core_exclude_cluster4_sensitivity", "All-beta versus exclude-cluster-4 sensitivity",
        "All parent-labelled beta cells remain primary; omitting review cluster 4 is a descriptive technical sensitivity and does not define a revised QC rule.",
        cluster4_sensitivity,
        lambda: _draw_cluster4_sensitivity(cluster4_sensitivity),
    )

    paga_source = _paga_source(data)
    add("beta_paga_undirected", "Undirected beta-state graph", "PAGA edge weights summarize cluster connectivity and carry no temporal or causal direction.", paga_source, lambda: _draw_paga(paga_source))
    distance_source = _umap_source(data, ["beta_relative_state_distance", "library_id", PREFERRED_CLUSTER_KEY])
    distance_note = gate["interpretation"]
    add("beta_relative_state_distance", "DPT relative state distance", distance_note, distance_source, lambda: _draw_continuous_grid(distance_source, ["beta_relative_state_distance"], "Hypothesis-generating relative state distance", ncols=1))
    add("beta_trajectory_sensitivity", "Trajectory sensitivity", "Neighbor, PC, seed and root perturbations are compared with the frozen baseline; thresholds are in trajectory_gate.json.", sensitivity, lambda: _draw_sensitivity(sensitivity))
    add("beta_balanced_subsample_sensitivity", "Equal-cell subsampling sensitivity", "Equal cell counts assess technical stability only and do not create biological replication.", balanced, lambda: _draw_balanced(balanced))
    return figures


def _write_documents(
    input_h5ad: Path,
    output_dir: Path,
    data: ad.AnnData,
    gate: dict[str, Any],
    figures: list[dict[str, str]],
    identity_clusters: pd.DataFrame,
    elapsed: float,
) -> None:
    method_plan = """# Beta-lineage extension method plan

## Frozen design

- Input: existing annotated H5AD, read only; beta cells selected by `cell_type == beta`.
- Unit: one PFOS library and one WT library; all condition inference remains ineligible.
- Representation: all QC-passing beta cells, sparse normalize-total log1p X, 3,000 eligible HVGs, PCA, unintegrated kNN graph and UMAP.
- Clustering: Leiden resolutions 0.3, 0.5 and 0.8; 0.5 is the descriptive display partition.
- Marker threshold: `abs(log2FC) >= 0.58` and raw cluster-marker `P < 0.05`, for annotation only.
- State scores: frozen maturity, GSIS, ER-stress, oxidative-stress and dedifferentiation programs with expression-bin matched control genes. Scores stay continuous.
- Cell cycle: Scanpy S/G2M scores and score-derived phase are descriptive.
- Identity/QC review: parent beta labels remain primary and unchanged; raw/corrected noncognate hormones, beta-identity TFs, complexity, mitochondrial fraction and doublet evidence generate neutral candidate flags. Cluster 4 exclusion is sensitivity only.
- State geometry: undirected PAGA is always reportable. DPT direction requires all predeclared connectivity, library-coverage, sensitivity, program-direction and cell-cycle gates.

## Knowledge reuse record

knowledge_id | candidate_idea | decision | reason | adaptation | validation
--- | --- | --- | --- | --- | ---
pmid-37697055 | mouse islet beta-state subclustering and multi-evidence annotation | adapt | species and tissue match; PFOS is absent from the reference | preserve abstention and use within-lineage programs without automatic state transfer | marker/program concordance and sensitivity
2026-09-01-pfos-beta-cell-response-boundaries | maturity, secretion and stress axes | background_only | direct studies constrain candidate axes but do not measure this F1 endpoint | score transcript programs only and avoid functional claims | QC correlations and library coverage
2026-08-30-replicate-aware-single-cell-and-computational-perturbation | sample-level inference gate | adopt | one library per condition cannot support population inference | prohibit condition P/FDR and label all composition/state contrasts descriptive | output schema audit
Scanpy official DPT/PAGA documentation | graph-based relative distance and coarse connectivity | adapt | snapshot geometry needs a root and does not encode real time | mature-high/stress-low root independent of condition; neighbor/PC/seed/root sensitivity | trajectory_gate.json
"""
    atomic_text(output_dir / "method_plan.md", method_plan)
    provenance = f"""# Provenance

- Started/finished: {datetime.now().isoformat(timespec='seconds')}
- Input: `{input_h5ad}`
- Input SHA-256: `{sha256(input_h5ad)}`
- Output: `{output_dir}`
- Cells/genes: {data.n_obs} / {data.n_vars}
- Python: {sys.version.split()[0]} ({platform.platform()})
- scanpy: {version('scanpy')}; anndata: {version('anndata')}; numpy: {version('numpy')}; scipy: {version('scipy')}
- Random seed: {SEED}
- Main graph: {MAIN_NEIGHBORS} neighbors, {MAIN_PCS} PCs, no sample integration
- Elapsed seconds: {elapsed:.3f}
- Exit status: completed
"""
    atomic_text(output_dir / "provenance.md", provenance)
    claim_map = pd.DataFrame([
        {"claim_id": "BETA-01", "claim": "Beta cells form descriptive expression-state partitions", "input": str(input_h5ad), "script": "src/pfos_scrna/extension_beta.py", "table": "tables/beta_cluster_composition.tsv", "figure": "figures/beta_umap_clusters.png", "test": "tests/test_extension_beta.py", "scope": "captured cells", "status": "supported_descriptive"},
        {"claim_id": "BETA-02", "claim": "Mpc1/Mpc2/Atg7 are localized within beta-state geometry", "input": str(input_h5ad), "script": "src/pfos_scrna/extension_beta.py", "table": "source_data/beta_target_feature_umap.tsv", "figure": "figures/beta_target_feature_umap.png", "test": "tests/test_extension_beta.py", "scope": "expression localization", "status": "supported_descriptive"},
        {"claim_id": "BETA-03", "claim": "Continuous maturity, GSIS and stress programs can be mapped", "input": str(input_h5ad), "script": "src/pfos_scrna/extension_beta.py", "table": "tables/beta_score_summary.tsv", "figure": "figures/beta_umap_state_scores.png", "test": "tests/test_extension_beta.py", "scope": "program scores", "status": "supported_descriptive"},
        {"claim_id": "BETA-04", "claim": "Directional trajectory interpretation", "input": str(input_h5ad), "script": "src/pfos_scrna/extension_beta.py", "table": "trajectory_gate.json", "figure": "figures/beta_paga_undirected.png", "test": "tests/test_extension_beta.py", "scope": "adult beta-state geometry", "status": gate["status"]},
        {"claim_id": "BETA-05", "claim": "PFOS causes beta-state, composition or trajectory change", "input": str(input_h5ad), "script": "src/pfos_scrna/extension_beta.py", "table": "trajectory_gate.json", "figure": "", "test": "tests/test_extension_beta.py", "scope": "animal-level condition inference", "status": "NO-GO_n1_per_condition"},
        {"claim_id": "BETA-06", "claim": "Parent beta labels include a cluster requiring identity/QC review", "input": str(input_h5ad), "script": "src/pfos_scrna/extension_beta.py", "table": "tables/beta_identity_cluster_review.tsv", "figure": "figures/beta_identity_hormone_qc_audit.png", "test": "tests/test_extension_beta.py", "scope": "neutral annotation review; no automatic relabeling", "status": "supported_review_flag"},
    ])
    atomic_tsv(output_dir / "claim_result_map.tsv", claim_map)
    atomic_text(output_dir / "unresolved.md", "# Unresolved\n\n- Independent animal replication is required for PFOS-WT inference.\n- Cross-sectional DPT requires orthogonal temporal or perturbational validation before temporal interpretation.\n- Cluster 4 remains a parent-labelled beta candidate with low-complexity and residual noncognate-hormone uncertainty; orthogonal annotation evidence is required before relabeling or exclusion.\n")
    review_row = identity_clusters.loc[identity_clusters["beta_cluster"].astype(str).eq(IDENTITY_REVIEW_CLUSTER)]
    review_metrics = json.loads(review_row.to_json(orient="records"))[0] if len(review_row) else {"cluster_review_status": "review_cluster_absent"}
    manifest = {
        "module": "beta_states",
        "status": "completed",
        "method": "unintegrated beta-lineage Scanpy atlas with continuous program scores, cell-cycle audit, PAGA and gated DPT sensitivity",
        "metrics": {"cells": data.n_obs, "genes": data.n_vars, "libraries": int(data.obs["library_id"].nunique()), "trajectory_status": gate["status"], "identity_review_cluster": review_metrics, "elapsed_seconds": elapsed},
        "figures": figures,
        "tables": sorted(str(path.relative_to(output_dir)) for path in (output_dir / "tables").glob("*.tsv")),
        "objects": ["beta_states.h5ad", "trajectory_gate.json", "sensitivity.tsv", "beta_identity_review_thresholds.json"],
        "condition_inference": "NO-GO_n1_per_condition",
    }
    atomic_json(output_dir / "report_manifest.json", manifest)


def run_beta(input_h5ad: Path, output_dir: Path) -> dict[str, Any]:
    """Run the complete, descriptive beta-lineage extension.

    The input is read only. Every formal result is written below ``output_dir``.
    No condition-level hypothesis test is performed.
    """
    input_h5ad = input_h5ad.resolve()
    output_dir = output_dir.resolve()
    if not input_h5ad.is_file():
        raise FileNotFoundError(input_h5ad)
    for directory in (output_dir, output_dir / "tables", output_dir / "figures", output_dir / "source_data"):
        directory.mkdir(parents=True, exist_ok=True)
    project_root = Path(__file__).resolve().parents[2]
    protected_paths = [input_h5ad, Path(__file__).resolve(), project_root / "pyproject.toml", project_root / "uv.lock"]
    protected_paths = [path for path in protected_paths if path.is_file()]
    execution_inputs = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "files": [
            {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}
            for path in protected_paths
        ],
    }
    atomic_json(output_dir / "execution_inputs.json", execution_inputs)
    start = monotonic()
    log_lines = [f"{datetime.now().isoformat(timespec='seconds')} START input={input_h5ad}"]

    data = _load_beta(input_h5ad)
    log_lines.append(f"loaded beta cells={data.n_obs} genes={data.n_vars}")
    gene_contract = _score_programs(data)
    _build_atlas(data)
    markers = _cluster_markers(data)
    gate, sensitivity, trajectory_coverage = _trajectory_analysis(data)
    tables = _cluster_and_library_tables(data)
    balanced = _balanced_subsampling(data)
    identity_cells, identity_genes, identity_clusters, identity_thresholds = _identity_uncertainty_audit(data)
    cluster4_sensitivity = _cluster4_exclusion_sensitivity(data)

    atomic_tsv(output_dir / "tables" / "beta_gene_set_contract.tsv", gene_contract)
    atomic_tsv(output_dir / "tables" / "beta_cluster_markers_all.tsv", markers)
    atomic_tsv(output_dir / "tables" / "beta_cluster_composition.tsv", tables["cluster_composition"])
    atomic_tsv(output_dir / "tables" / "beta_score_summary.tsv", tables["score_summary"])
    atomic_tsv(output_dir / "tables" / "beta_qc_score_associations.tsv", tables["qc_associations"])
    atomic_tsv(output_dir / "tables" / "beta_cell_cycle_composition.tsv", tables["cell_cycle"])
    atomic_tsv(output_dir / "tables" / "beta_trajectory_library_coverage.tsv", trajectory_coverage)
    atomic_tsv(output_dir / "tables" / "beta_balanced_subsample_sensitivity.tsv", balanced)
    atomic_tsv(output_dir / "tables" / "beta_identity_cell_candidates.tsv", identity_cells)
    atomic_tsv(output_dir / "tables" / "beta_identity_gene_layer_summary.tsv", identity_genes)
    atomic_tsv(output_dir / "tables" / "beta_identity_cluster_review.tsv", identity_clusters)
    atomic_tsv(output_dir / "tables" / "beta_core_exclude_cluster4_sensitivity.tsv", cluster4_sensitivity)
    atomic_tsv(output_dir / "sensitivity.tsv", sensitivity)
    atomic_json(output_dir / "trajectory_gate.json", gate)
    atomic_json(output_dir / "beta_identity_review_thresholds.json", identity_thresholds)

    data.uns["beta_analysis_contract"] = {
        "scope": "adult beta-cell descriptive state geometry",
        "integration": "not_used_condition_sample_confounding",
        "condition_inference": "NO-GO_n1_per_condition",
        "marker_min_abs_log2fc": 0.58,
        "marker_raw_pvalue": 0.05,
        "marker_scope": "cluster_annotation_only",
        "state_scores": "continuous_no_forced_class",
        "trajectory_status": gate["status"],
        "identity_review": "neutral_candidate_flags_parent_labels_retained",
        "primary_beta_scope": "all_parent-labelled_beta_cells",
        "cluster4_exclusion": "technical_descriptive_sensitivity_only",
    }
    figures = _write_plots(
        data, output_dir, markers, tables, sensitivity, gate, balanced,
        identity_cells, identity_genes, identity_clusters, identity_thresholds, cluster4_sensitivity,
    )
    data.strings_to_categoricals()
    data.write_h5ad(output_dir / "beta_states.h5ad", compression="gzip")
    final_hashes = {str(path): sha256(path) for path in protected_paths}
    changed = [row["path"] for row in execution_inputs["files"] if final_hashes[row["path"]] != row["sha256"]]
    if changed:
        raise RuntimeError(f"protected execution inputs changed during beta run: {changed}")
    execution_inputs["verified_unchanged_at_completion"] = datetime.now().isoformat(timespec="seconds")
    execution_inputs["verification_status"] = "PASS"
    atomic_json(output_dir / "execution_inputs.json", execution_inputs)
    elapsed = monotonic() - start
    _write_documents(input_h5ad, output_dir, data, gate, figures, identity_clusters, elapsed)
    log_lines.extend([
        f"trajectory_status={gate['status']}",
        f"figures={len(figures)} paired PNG/PDF/source TSV",
        f"elapsed_seconds={elapsed:.3f}",
        f"{datetime.now().isoformat(timespec='seconds')} COMPLETED",
    ])
    atomic_text(output_dir / "run_log.md", "# Run log\n\n```text\n" + "\n".join(log_lines) + "\n```\n")
    return {
        "status": "completed",
        "cells": data.n_obs,
        "genes": data.n_vars,
        "trajectory_status": gate["status"],
        "output_dir": str(output_dir),
        "figures": len(figures),
        "elapsed_seconds": elapsed,
    }

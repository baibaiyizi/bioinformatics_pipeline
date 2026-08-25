#!/usr/bin/env python3
"""Validate paired RNA/ATAC AnnData, perform modality QC and compute ATAC LSI."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize


class ContractError(ValueError):
    pass


def resolve(base: Path, raw: str) -> Path:
    path = Path(raw).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    required = {"rna_h5ad", "atac_h5ad", "paired", "sample_column", "output_dir"}
    missing = sorted(required - set(config))
    if missing:
        raise ContractError(f"配置缺少字段: {', '.join(missing)}")
    return config


def matrix_metrics(matrix: Any) -> tuple[np.ndarray, np.ndarray]:
    if sp.issparse(matrix):
        counts = np.asarray(matrix.sum(axis=1)).ravel()
        features = np.asarray((matrix > 0).sum(axis=1)).ravel()
    else:
        arr = np.asarray(matrix)
        counts = arr.sum(axis=1)
        features = (arr > 0).sum(axis=1)
    return counts.astype(float), features.astype(int)


def validate_counts(adata: ad.AnnData, label: str, sample_column: str) -> None:
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ContractError(f"{label} 的 obs_names/var_names 必须唯一")
    if sample_column not in adata.obs:
        raise ContractError(f"{label}.obs 缺少 sample_column: {sample_column}")
    if adata.obs[sample_column].isna().any():
        raise ContractError(f"{label}.obs[{sample_column}] 存在缺失")
    matrix = adata.X
    values = matrix.data if sp.issparse(matrix) else np.asarray(matrix).ravel()
    if not np.isfinite(values).all() or (values < 0).any():
        raise ContractError(f"{label}.X 必须是非负有限 counts")
    if not np.allclose(values, np.round(values)):
        raise ContractError(f"{label}.X 不是整数型原始 counts")


def filter_modality(adata: ad.AnnData, min_counts: int, min_features: int, prefix: str) -> tuple[ad.AnnData, pd.DataFrame]:
    counts, features = matrix_metrics(adata.X)
    keep = (counts >= min_counts) & (features >= min_features)
    metrics = pd.DataFrame({
        "cell_id": adata.obs_names.astype(str),
        f"{prefix}_total_counts": counts,
        f"{prefix}_n_features": features,
        f"{prefix}_pass_qc": keep,
    })
    result = adata[keep].copy()
    result.obs[f"{prefix}_total_counts"] = counts[keep]
    result.obs[f"{prefix}_n_features"] = features[keep]
    result.layers["counts"] = result.X.copy()
    return result, metrics


def tfidf_lsi(matrix: Any, n_components: int, seed: int) -> np.ndarray:
    X = sp.csr_matrix(matrix, dtype=float)
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    if (row_sums <= 0).any():
        raise ContractError("ATAC QC 后仍存在总 counts 为 0 的细胞")
    tf = sp.diags(1.0 / row_sums) @ X
    detected = np.asarray((X > 0).sum(axis=0)).ravel()
    idf = np.log1p(X.shape[0] / np.maximum(detected, 1))
    tfidf = normalize(tf @ sp.diags(idf), norm="l2", axis=1)
    max_components = min(tfidf.shape[0] - 1, tfidf.shape[1] - 1)
    if max_components < 2:
        raise ContractError("ATAC 矩阵过小，无法计算至少两个 LSI 分量")
    n_components = min(n_components, max_components)
    return TruncatedSVD(n_components=n_components, random_state=seed).fit_transform(tfidf)


def run(config_path: Path) -> Path:
    config_path = config_path.resolve()
    config = load_config(config_path)
    base = config_path.parent
    rna_path = resolve(base, config["rna_h5ad"])
    atac_path = resolve(base, config["atac_h5ad"])
    out = resolve(base, config["output_dir"])
    if not rna_path.is_file() or not atac_path.is_file():
        raise ContractError("RNA 或 ATAC h5ad 不存在")
    rna = ad.read_h5ad(rna_path)
    atac = ad.read_h5ad(atac_path)
    sample_column = str(config["sample_column"])
    validate_counts(rna, "RNA", sample_column)
    validate_counts(atac, "ATAC", sample_column)
    if bool(config["paired"]):
        if set(rna.obs_names) != set(atac.obs_names):
            raise ContractError("paired=true 时 RNA 与 ATAC 细胞集合必须完全一致；禁止静默取交集")
        atac = atac[rna.obs_names].copy()
        if not rna.obs[sample_column].astype(str).equals(atac.obs[sample_column].astype(str)):
            raise ContractError("同一 barcode 在 RNA 与 ATAC 中的 sample_id 不一致")
    out.mkdir(parents=True, exist_ok=True)
    rna_qc, rna_metrics = filter_modality(rna, int(config.get("rna_min_counts", 1)), int(config.get("rna_min_features", 1)), "rna")
    atac_qc, atac_metrics = filter_modality(atac, int(config.get("atac_min_counts", 1)), int(config.get("atac_min_features", 1)), "atac")
    if bool(config["paired"]):
        shared = rna_qc.obs_names.intersection(atac_qc.obs_names, sort=False)
        if len(shared) == 0:
            raise ContractError("两模态 QC 后没有共享细胞")
        rna_qc = rna_qc[shared].copy()
        atac_qc = atac_qc[shared].copy()
    lsi = tfidf_lsi(atac_qc.X, int(config.get("lsi_components", 30)), int(config.get("seed", 20260824)))
    atac_qc.obsm["X_lsi"] = lsi
    rna_qc.write_h5ad(out / "rna_qc.h5ad", compression="gzip")
    atac_qc.write_h5ad(out / "atac_qc_lsi.h5ad", compression="gzip")
    rna_metrics.merge(atac_metrics, on="cell_id", how="outer").to_csv(out / "cell_qc.csv", index=False)
    pd.DataFrame({"cell_id": atac_qc.obs_names.astype(str), **{f"LSI{i+1}": lsi[:, i] for i in range(lsi.shape[1])}}).to_csv(out / "atac_lsi.csv", index=False)
    h5mu_status = "not_requested"
    if bool(config.get("write_h5mu", True)):
        if importlib.util.find_spec("mudata") is None:
            h5mu_status = "dependency_missing"
        else:
            import mudata as md
            md.MuData({"rna": rna_qc, "atac": atac_qc}).write_h5mu(out / "paired_multiome.h5mu")
            h5mu_status = "completed"
    sample_counts = rna_qc.obs[sample_column].astype(str).value_counts().rename_axis("sample_id").reset_index(name="n_cells_after_joint_qc")
    sample_counts.to_csv(out / "sample_summary.csv", index=False)
    audit = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "rna_input": str(rna_path), "rna_sha256": sha256(rna_path),
        "atac_input": str(atac_path), "atac_sha256": sha256(atac_path),
        "paired": bool(config["paired"]), "sample_column": sample_column,
        "n_rna_input": rna.n_obs, "n_atac_input": atac.n_obs,
        "n_joint_after_qc": int(len(rna_qc)) if bool(config["paired"]) else None,
        "lsi_components": int(lsi.shape[1]), "seed": int(config.get("seed", 20260824)),
        "h5mu_status": h5mu_status,
        "multivi_status": "not_run; requires explicit design, scvi-tools and resource audit",
        "inference_boundary": "cell-level embeddings do not create biological replicates; DE/DA defaults to sample pseudobulk",
    }
    (out / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "effective_config.yml").write_text(yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    try:
        print(run(args.config))
    except (ContractError, OSError, yaml.YAMLError) as exc:
        print(f"[CONTRACT_ERROR] {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

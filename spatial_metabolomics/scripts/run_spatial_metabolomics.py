#!/usr/bin/env python3
"""Auditable spatial metabolomics preprocessing and within-sample spatial statistics."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from sklearn.neighbors import NearestNeighbors


class ContractError(ValueError):
    pass


def resolve(base: Path, raw: str | None) -> Path | None:
    if raw in (None, ""):
        return None
    path = Path(str(raw)).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    required = {"mode", "input", "output_dir"}
    missing = sorted(required - set(config))
    if missing:
        raise ContractError(f"配置缺少字段: {', '.join(missing)}")
    if config["mode"] not in {"long_table", "imzml"}:
        raise ContractError("mode 仅支持 long_table 或 imzml")
    return config


def load_long_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"sample_id", "pixel_id", "x", "y", "mz", "intensity"}
    missing = sorted(required - set(frame))
    if missing:
        raise ContractError(f"长表缺少列: {', '.join(missing)}")
    if "roi" not in frame:
        frame["roi"] = "unassigned"
    for col in ["x", "y", "mz", "intensity"]:
        frame[col] = pd.to_numeric(frame[col], errors="raise")
    if frame[["sample_id", "pixel_id", "x", "y", "mz", "intensity"]].isna().any().any():
        raise ContractError("光谱长表关键字段存在缺失")
    if (frame["mz"] <= 0).any() or (frame["intensity"] < 0).any():
        raise ContractError("m/z 必须为正，intensity 必须非负")
    coords_per_pixel = frame.groupby(["sample_id", "pixel_id"])[["x", "y"]].nunique()
    if (coords_per_pixel > 1).any().any():
        raise ContractError("同一 sample_id/pixel_id 对应多个坐标")
    return frame


def load_imzml_manifest(path: Path) -> pd.DataFrame:
    if importlib.util.find_spec("pyimzml") is None:
        raise ContractError("mode=imzml 需要安装 pyimzML；当前环境缺失")
    from pyimzml.ImzMLParser import ImzMLParser
    manifest = pd.read_csv(path)
    required = {"sample_id", "imzml", "ibd"}
    missing = sorted(required - set(manifest))
    if missing:
        raise ContractError(f"imzML manifest 缺少列: {', '.join(missing)}")
    rows: list[pd.DataFrame] = []
    for record in manifest.to_dict("records"):
        imzml = resolve(path.parent, record["imzml"])
        ibd = resolve(path.parent, record["ibd"])
        if imzml is None or ibd is None or not imzml.is_file() or not ibd.is_file():
            raise ContractError(f"样本 {record['sample_id']} 的 imzML/ibd 不完整")
        expected_ibd = imzml.with_suffix(".ibd")
        if ibd.resolve() != expected_ibd.resolve():
            raise ContractError("pyimzML 按同 basename 读取 .ibd；manifest 中 imzML/ibd basename 不一致")
        parser = ImzMLParser(str(imzml))
        for index, coord in enumerate(parser.coordinates):
            mz, intensity = parser.getspectrum(index)
            rows.append(pd.DataFrame({
                "sample_id": str(record["sample_id"]), "pixel_id": f"{record['sample_id']}:{index}",
                "x": coord[0], "y": coord[1], "roi": str(record.get("roi", "unassigned")),
                "mz": mz, "intensity": intensity,
            }))
    if not rows:
        raise ContractError("imzML manifest 未读取到任何光谱")
    return pd.concat(rows, ignore_index=True)


def build_matrices(frame: pd.DataFrame, decimals: int, min_features: int, tic_quantile: float):
    data = frame.copy()
    data["mz_bin"] = data["mz"].round(decimals)
    grouped = data.groupby(["sample_id", "pixel_id", "x", "y", "roi", "mz_bin"], as_index=False)["intensity"].sum()
    matrix = grouped.pivot_table(index=["sample_id", "pixel_id", "x", "y", "roi"], columns="mz_bin", values="intensity", fill_value=0.0, aggfunc="sum")
    values = matrix.to_numpy(dtype=float)
    tic = values.sum(axis=1)
    detected = (values > 0).sum(axis=1)
    sample_ids = matrix.index.get_level_values("sample_id")
    thresholds = pd.Series(tic, index=matrix.index).groupby(level="sample_id").quantile(tic_quantile)
    tic_cutoff = np.asarray([thresholds.loc[s] for s in sample_ids])
    keep = (tic > 0) & (tic >= tic_cutoff) & (detected >= min_features)
    qc = matrix.index.to_frame(index=False)
    qc["tic"] = tic
    qc["detected_features"] = detected
    qc["tic_cutoff_within_sample"] = tic_cutoff
    qc["pass_qc"] = keep
    filtered = matrix.loc[keep].copy()
    fvalues = filtered.to_numpy(dtype=float)
    scale = np.divide(1_000_000.0, fvalues.sum(axis=1), out=np.zeros(len(fvalues)), where=fvalues.sum(axis=1) > 0)
    normalized = pd.DataFrame(np.log1p(fvalues * scale[:, None]), index=filtered.index, columns=filtered.columns)
    return filtered, normalized, qc


def morans_i(values: np.ndarray, coords: np.ndarray, k: int) -> float:
    n = len(values)
    if n < 4 or np.allclose(values, values[0]):
        return float("nan")
    k = min(k, n - 1)
    indices = NearestNeighbors(n_neighbors=k + 1).fit(coords).kneighbors(return_distance=False)[:, 1:]
    weights = np.zeros((n, n), dtype=float)
    for i, neighbors in enumerate(indices):
        weights[i, neighbors] = 1.0
    weights = np.maximum(weights, weights.T)
    z = values - values.mean()
    denominator = float(np.sum(z * z))
    s0 = float(weights.sum())
    return float(n / s0 * np.sum(weights * np.outer(z, z)) / denominator) if denominator > 0 and s0 > 0 else float("nan")


def spatial_statistics(normalized: pd.DataFrame, k: int, top_features: int, permutations: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for sample_id in normalized.index.get_level_values("sample_id").unique():
        sub = normalized.xs(sample_id, level="sample_id", drop_level=False)
        coords = sub.index.to_frame(index=False)[["x", "y"]].to_numpy(dtype=float)
        variances = sub.var(axis=0).sort_values(ascending=False)
        for feature in variances.head(top_features).index:
            values = sub[feature].to_numpy(dtype=float)
            observed = morans_i(values, coords, k)
            if not np.isfinite(observed) or permutations <= 0:
                p_value = float("nan")
            else:
                null = np.asarray([morans_i(rng.permutation(values), coords, k) for _ in range(permutations)])
                p_value = float((1 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1))
            rows.append({"sample_id": sample_id, "mz_bin": feature, "moran_i": observed, "permutation_p": p_value, "n_pixels": len(sub)})
    return pd.DataFrame(rows)


def validate_annotations(path: Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    frame = pd.read_csv(path)
    required = {"mz_bin", "annotation", "identification_level", "standard_confirmed", "source_version"}
    missing = sorted(required - set(frame))
    if missing:
        raise ContractError(f"注释表缺少列: {', '.join(missing)}")
    levels = pd.to_numeric(frame["identification_level"], errors="raise")
    if not levels.isin([1, 2, 3, 4]).all():
        raise ContractError("identification_level 仅允许 1-4")
    confirmed = frame["standard_confirmed"].astype(str).str.lower().isin(["true", "1", "yes"])
    if ((levels == 1) & ~confirmed).any():
        raise ContractError("level 1 注释必须有 standard_confirmed=true")
    return frame


def run(config_path: Path) -> Path:
    config_path = config_path.resolve()
    config = load_config(config_path)
    input_path = resolve(config_path.parent, config["input"])
    out = resolve(config_path.parent, config["output_dir"])
    annotation_path = resolve(config_path.parent, config.get("annotation_table"))
    assert input_path is not None and out is not None
    if not input_path.is_file():
        raise ContractError(f"输入不存在: {input_path}")
    frame = load_long_table(input_path) if config["mode"] == "long_table" else load_imzml_manifest(input_path)
    annotations = validate_annotations(annotation_path)
    matrix, normalized, qc = build_matrices(frame, int(config.get("mz_decimals", 3)), int(config.get("min_detected_features", 3)), float(config.get("min_tic_quantile", 0.01)))
    if matrix.empty:
        raise ContractError("QC 后无像素保留")
    out.mkdir(parents=True, exist_ok=True)
    qc.to_csv(out / "pixel_qc.csv", index=False)
    matrix.reset_index().to_csv(out / "feature_matrix.csv", index=False)
    normalized.reset_index().to_csv(out / "normalized_log_tic.csv", index=False)
    roi_summary = normalized.groupby(level=["sample_id", "roi"]).mean().reset_index()
    roi_summary.to_csv(out / "roi_sample_summary.csv", index=False)
    spatial = spatial_statistics(normalized, int(config.get("spatial_k", 4)), int(config.get("spatial_top_features", 50)), int(config.get("permutations", 99)), int(config.get("seed", 20260824)))
    spatial.to_csv(out / "within_sample_spatial_stats.csv", index=False)
    if annotations is not None:
        annotations.to_csv(out / "annotation_evidence.csv", index=False)
    audit = {
        "created_at": datetime.now(timezone.utc).isoformat(), "mode": config["mode"],
        "input": str(input_path), "input_sha256": sha256(input_path),
        "n_samples": int(matrix.index.get_level_values("sample_id").nunique()),
        "n_pixels_input": int(qc.shape[0]), "n_pixels_pass": int(qc["pass_qc"].sum()),
        "n_features": int(matrix.shape[1]), "mz_decimals": int(config.get("mz_decimals", 3)),
        "seed": int(config.get("seed", 20260824)), "annotation_status": "provided" if annotations is not None else "not_provided",
        "spatial_inference": "within-sample exploratory only; pixels/ROIs are not biological replicates",
        "cross_sample_inference": "not_performed",
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
    except (ContractError, OSError, yaml.YAMLError, pd.errors.ParserError) as exc:
        print(f"[CONTRACT_ERROR] {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Leakage-resistant nested validation for frozen, sample-level feature tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import sklearn
import yaml
from sklearn.base import clone
from sklearn.dummy import DummyClassifier, DummyRegressor
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.feature_selection import SelectKBest, f_classif, f_regression
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    brier_score_loss,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import (
    GridSearchCV,
    GroupKFold,
    KFold,
    StratifiedGroupKFold,
    StratifiedKFold,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


class ContractError(ValueError):
    """Raised when the frozen feature table violates the analysis contract."""


@dataclass(frozen=True)
class Dataset:
    frame: pd.DataFrame
    features: list[str]
    target: pd.Series
    groups: pd.Series | None
    ids: pd.Series


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_path(config_path: Path, raw: str | None) -> Path | None:
    if raw in (None, ""):
        return None
    path = Path(raw).expanduser()
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    required = {"input_csv", "target_column", "task", "feature_columns", "output_dir"}
    missing = sorted(required - set(config))
    if missing:
        raise ContractError(f"配置缺少字段: {', '.join(missing)}")
    if config["task"] not in {"classification", "regression"}:
        raise ContractError("task 仅支持 classification 或 regression")
    if config.get("model", "logistic" if config["task"] == "classification" else "ridge") not in {
        "logistic", "ridge", "random_forest"
    }:
        raise ContractError("model 仅支持 logistic、ridge 或 random_forest")
    return config


def load_dataset(path: Path, config: dict[str, Any], *, external: bool = False) -> Dataset:
    if not path.is_file():
        raise ContractError(f"输入文件不存在: {path}")
    frame = pd.read_csv(path)
    target_col = config["target_column"]
    group_col = config.get("group_column")
    id_col = config.get("id_column")
    features = list(config["feature_columns"])
    if not features or len(features) != len(set(features)):
        raise ContractError("feature_columns 必须是非空且无重复的显式列表")
    required = set(features) | {target_col}
    if group_col:
        required.add(group_col)
    if id_col:
        required.add(id_col)
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ContractError(f"输入表缺少列: {', '.join(missing)}")
    if frame[target_col].isna().any():
        raise ContractError("目标列存在缺失值")
    non_numeric = [c for c in features if not pd.api.types.is_numeric_dtype(frame[c])]
    if non_numeric:
        raise ContractError(f"特征必须为数值列: {', '.join(non_numeric)}")
    values = frame[features].to_numpy(dtype=float)
    if np.isinf(values).any():
        raise ContractError("特征包含正负无穷值")
    if id_col:
        if frame[id_col].isna().any():
            raise ContractError(f"id_column {id_col} 存在缺失值")
        ids = frame[id_col].astype(str)
        if ids.duplicated().any():
            raise ContractError(f"id_column {id_col} 不是唯一键")
    else:
        ids = pd.Series([f"row_{i}" for i in range(len(frame))], index=frame.index)
    if group_col and frame[group_col].isna().any():
        raise ContractError("group_column 存在缺失值")
    groups = frame[group_col].astype(str) if group_col else None
    target = frame[target_col]
    if config["task"] == "classification":
        labels = list(pd.unique(target))
        if len(labels) != 2:
            raise ContractError(f"最小闭环仅支持二分类，检测到 {len(labels)} 类")
        positive = config.get("positive_label")
        if positive is None:
            positive = sorted(labels, key=str)[-1]
        matches = [x for x in labels if str(x) == str(positive)]
        if len(matches) != 1:
            raise ContractError(f"positive_label {positive!r} 不在目标列中")
        target = (target == matches[0]).astype(int)
    elif not pd.api.types.is_numeric_dtype(target):
        raise ContractError("回归目标必须是数值列")
    if len(frame) < 8 and not external:
        raise ContractError("开发集少于 8 行，不能进行可靠的嵌套交叉验证 smoke loop")
    return Dataset(frame, features, target, groups, ids)


def make_splitter(task: str, splits: int, seed: int, grouped: bool):
    if splits < 2:
        raise ContractError("交叉验证折数必须至少为 2")
    if grouped:
        # StratifiedGroupKFold 的 shuffle 只会打乱近似同分布的组；当每组是纯类别时，
        # shuffle 反而可能生成单类别 fold。固定顺序优先保证可计算性和可复现性。
        return StratifiedGroupKFold(splits, shuffle=False) if task == "classification" else GroupKFold(splits)
    return StratifiedKFold(splits, shuffle=True, random_state=seed) if task == "classification" else KFold(splits, shuffle=True, random_state=seed)


def validate_splits(data: Dataset, config: dict[str, Any]) -> None:
    outer = int(config.get("outer_splits", 5))
    inner = int(config.get("inner_splits", 3))
    if data.groups is not None:
        n_groups = data.groups.nunique()
        if n_groups < max(outer, inner):
            raise ContractError(f"独立组数 {n_groups} 少于所需折数 {max(outer, inner)}")
    elif config["task"] == "classification":
        least = int(data.target.value_counts().min())
        if least < max(outer, inner):
            raise ContractError(f"较少类别仅 {least} 个样本，少于所需折数 {max(outer, inner)}")


def materialize_splits(splitter: Any, X: pd.DataFrame, y: np.ndarray, groups: np.ndarray | None, task: str, label: str):
    splits = list(splitter.split(X, y, groups))
    if task == "classification":
        for fold, (train_idx, test_idx) in enumerate(splits, start=1):
            if len(np.unique(y[train_idx])) < 2 or len(np.unique(y[test_idx])) < 2:
                raise ContractError(f"{label} 第 {fold} 折出现单类别；请减少折数、增加独立组或重新审视分组设计")
    return splits


def build_pipeline(config: dict[str, Any]) -> tuple[Pipeline, dict[str, list[Any]], str]:
    task = config["task"]
    seed = int(config.get("seed", 20260824))
    k = config.get("feature_k", "all")
    if k != "all":
        k = int(k)
        if k < 1:
            raise ContractError("feature_k 必须为 all 或正整数")
    selector = SelectKBest(f_classif if task == "classification" else f_regression, k=k)
    model_name = config.get("model", "logistic" if task == "classification" else "ridge")
    if task == "classification" and model_name == "logistic":
        estimator = LogisticRegression(max_iter=5000, class_weight="balanced", random_state=seed)
        grid = {"model__C": [0.1, 1.0, 10.0]}
        scoring = "roc_auc"
    elif task == "regression" and model_name == "ridge":
        estimator = Ridge()
        grid = {"model__alpha": [0.1, 1.0, 10.0]}
        scoring = "neg_root_mean_squared_error"
    elif task == "classification" and model_name == "random_forest":
        estimator = RandomForestClassifier(class_weight="balanced", random_state=seed, n_jobs=1)
        grid = {"model__n_estimators": [100, 300], "model__max_depth": [None, 5]}
        scoring = "roc_auc"
    elif task == "regression" and model_name == "random_forest":
        estimator = RandomForestRegressor(random_state=seed, n_jobs=1)
        grid = {"model__n_estimators": [100, 300], "model__max_depth": [None, 5]}
        scoring = "neg_root_mean_squared_error"
    else:
        raise ContractError(f"模型 {model_name} 与任务 {task} 不兼容")
    pipeline = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("selector", selector),
        ("model", estimator),
    ])
    return pipeline, grid, scoring


def classification_metrics(y: np.ndarray, score: np.ndarray) -> dict[str, float]:
    pred = (score >= 0.5).astype(int)
    return {
        "roc_auc": float(roc_auc_score(y, score)),
        "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "brier": float(brier_score_loss(y, score)),
    }


def regression_metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "rmse": float(math.sqrt(mean_squared_error(y, pred))),
        "mae": float(mean_absolute_error(y, pred)),
        "r2": float(r2_score(y, pred)),
    }


def predict_score(model: Any, X: pd.DataFrame, task: str) -> np.ndarray:
    if task == "classification":
        return np.asarray(model.predict_proba(X)[:, 1], dtype=float)
    return np.asarray(model.predict(X), dtype=float)


def decision_curve(y: np.ndarray, score: np.ndarray) -> pd.DataFrame:
    rows = []
    n = len(y)
    prevalence = float(np.mean(y))
    for threshold in np.linspace(0.05, 0.95, 19):
        pred = score >= threshold
        tp = int(np.sum(pred & (y == 1)))
        fp = int(np.sum(pred & (y == 0)))
        weight = threshold / (1 - threshold)
        rows.append({
            "threshold": threshold,
            "net_benefit_model": tp / n - fp / n * weight,
            "net_benefit_all": prevalence - (1 - prevalence) * weight,
            "net_benefit_none": 0.0,
        })
    return pd.DataFrame(rows)


def calibration_table(y: np.ndarray, score: np.ndarray) -> pd.DataFrame:
    bins = pd.qcut(pd.Series(score), q=min(10, len(np.unique(score))), duplicates="drop")
    table = pd.DataFrame({"observed": y, "predicted": score, "bin": bins.astype(str)})
    return table.groupby("bin", observed=True).agg(n=("observed", "size"), mean_predicted=("predicted", "mean"), observed_rate=("observed", "mean")).reset_index()


def nested_validate(data: Dataset, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, Pipeline, dict[str, Any]]:
    validate_splits(data, config)
    seed = int(config.get("seed", 20260824))
    outer = make_splitter(config["task"], int(config.get("outer_splits", 5)), seed, data.groups is not None)
    pipeline, grid, scoring = build_pipeline(config)
    X = data.frame[data.features]
    y = data.target.to_numpy()
    groups = None if data.groups is None else data.groups.to_numpy()
    predictions = np.full(len(data.frame), np.nan)
    fold_assignment = np.full(len(data.frame), -1, dtype=int)
    fold_rows: list[dict[str, Any]] = []
    split_iter = materialize_splits(outer, X, y, groups, config["task"], "外层交叉验证")
    for fold, (train_idx, test_idx) in enumerate(split_iter, start=1):
        inner = make_splitter(config["task"], int(config.get("inner_splits", 3)), seed + fold, groups is not None)
        inner_splits = materialize_splits(
            inner,
            X.iloc[train_idx],
            y[train_idx],
            None if groups is None else groups[train_idx],
            config["task"],
            f"外层第 {fold} 折的内层交叉验证",
        )
        search = GridSearchCV(clone(pipeline), grid, scoring=scoring, cv=inner_splits, n_jobs=1, refit=True, error_score="raise")
        search.fit(X.iloc[train_idx], y[train_idx])
        score = predict_score(search.best_estimator_, X.iloc[test_idx], config["task"])
        predictions[test_idx] = score
        fold_assignment[test_idx] = fold
        metrics = classification_metrics(y[test_idx], score) if config["task"] == "classification" else regression_metrics(y[test_idx], score)
        fold_rows.append({"fold": fold, "n_train": len(train_idx), "n_test": len(test_idx), **metrics, "best_params": json.dumps(search.best_params_, sort_keys=True)})
    if np.isnan(predictions).any():
        raise RuntimeError("外层交叉验证未覆盖全部开发样本")
    inner_all = make_splitter(config["task"], int(config.get("inner_splits", 3)), seed + 1000, groups is not None)
    final_inner_splits = materialize_splits(inner_all, X, y, groups, config["task"], "最终调参")
    final_search = GridSearchCV(clone(pipeline), grid, scoring=scoring, cv=final_inner_splits, n_jobs=1, refit=True, error_score="raise")
    final_search.fit(X, y)
    pred_frame = pd.DataFrame({"id": data.ids.astype(str), "outer_fold": fold_assignment, "observed": y, "prediction": predictions})
    if data.groups is not None:
        pred_frame["group"] = data.groups.astype(str).to_numpy()
    overall = classification_metrics(y, predictions) if config["task"] == "classification" else regression_metrics(y, predictions)
    return pred_frame, pd.DataFrame(fold_rows), final_search.best_estimator_, {"overall": overall, "final_best_params": final_search.best_params_}


def permutation_table(model: Pipeline, data: Dataset, config: dict[str, Any]) -> pd.DataFrame:
    scoring = "roc_auc" if config["task"] == "classification" else "neg_root_mean_squared_error"
    result = permutation_importance(model, data.frame[data.features], data.target, scoring=scoring, n_repeats=10, random_state=int(config.get("seed", 20260824)), n_jobs=1)
    return pd.DataFrame({"feature": data.features, "importance_mean": result.importances_mean, "importance_sd": result.importances_std}).sort_values("importance_mean", ascending=False)


def maybe_shap(model: Pipeline, data: Dataset, output_dir: Path, enabled: bool) -> str:
    if not enabled:
        return "not_requested"
    try:
        import shap  # type: ignore
    except ImportError:
        return "dependency_missing"
    transformed = model[:-1].transform(data.frame[data.features])
    names = np.asarray(data.features)[model.named_steps["selector"].get_support()].tolist()
    try:
        explainer = shap.Explainer(model.named_steps["model"], transformed)
        values = explainer(transformed)
        np.savez_compressed(output_dir / "shap_values.npz", values=np.asarray(values.values), base_values=np.asarray(values.base_values), feature_names=np.asarray(names))
    except Exception as exc:  # explicit failure artifact, never silently replace SHAP
        (output_dir / "shap_error.txt").write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")
        return "failed"
    return "completed"


def write_model_card(path: Path, config: dict[str, Any], summary: dict[str, Any], shap_status: str, external_status: str) -> None:
    metric_lines = "\n".join(f"- {k}: {v:.6g}" for k, v in summary["overall"].items())
    path.write_text(
        "# 模型卡\n\n"
        f"- 任务：{config['task']}\n"
        f"- 模型：{config.get('model', 'default')}\n"
        f"- 验证：嵌套交叉验证；外层 {config.get('outer_splits', 5)} 折，内层 {config.get('inner_splits', 3)} 折\n"
        f"- 分组保护：{'启用' if config.get('group_column') else '未提供分组列'}\n"
        f"- 外部评估：{external_status}\n"
        f"- SHAP：{shap_status}\n"
        "- 推断边界：当前性能为预测关联，不证明因果或临床效用。\n\n"
        "## 开发集外层交叉验证性能\n\n"
        f"{metric_lines}\n",
        encoding="utf-8",
    )


def run(config_path: Path) -> Path:
    config_path = config_path.resolve()
    config = load_config(config_path)
    input_path = resolve_path(config_path, config["input_csv"])
    external_path = resolve_path(config_path, config.get("external_csv"))
    output_dir = resolve_path(config_path, config["output_dir"])
    assert input_path is not None and output_dir is not None
    output_dir.mkdir(parents=True, exist_ok=True)
    data = load_dataset(input_path, config)
    predictions, folds, model, summary = nested_validate(data, config)
    predictions.to_csv(output_dir / "oof_predictions.csv", index=False)
    folds.to_csv(output_dir / "fold_metrics.csv", index=False)
    if config["task"] == "classification":
        y = predictions["observed"].to_numpy(dtype=int)
        score = predictions["prediction"].to_numpy(dtype=float)
        calibration_table(y, score).to_csv(output_dir / "calibration.csv", index=False)
        decision_curve(y, score).to_csv(output_dir / "decision_curve.csv", index=False)
        dummy = DummyClassifier(strategy="prior").fit(np.zeros((len(y), 1)), y)
        summary["dummy_baseline"] = classification_metrics(y, dummy.predict_proba(np.zeros((len(y), 1)))[:, 1])
    else:
        y = predictions["observed"].to_numpy(dtype=float)
        dummy = DummyRegressor(strategy="mean").fit(np.zeros((len(y), 1)), y)
        summary["dummy_baseline"] = regression_metrics(y, dummy.predict(np.zeros((len(y), 1))))
    permutation_table(model, data, config).to_csv(output_dir / "permutation_importance.csv", index=False)
    joblib.dump(model, output_dir / "final_model.joblib")
    external_status = "not_provided"
    if external_path is not None:
        external = load_dataset(external_path, config, external=True)
        if external.features != data.features:
            raise ContractError("外部数据特征顺序与开发数据不一致")
        score = predict_score(model, external.frame[external.features], config["task"])
        external_predictions = pd.DataFrame({"id": external.ids.astype(str), "observed": external.target, "prediction": score})
        external_predictions.to_csv(output_dir / "external_predictions.csv", index=False)
        summary["external"] = classification_metrics(external.target.to_numpy(), score) if config["task"] == "classification" else regression_metrics(external.target.to_numpy(), score)
        external_status = "completed"
    shap_status = maybe_shap(model, data, output_dir, bool(config.get("run_shap", False)))
    audit = {
        "created_at": utc_now(),
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "external_input": None if external_path is None else str(external_path),
        "n_rows": len(data.frame),
        "n_features": len(data.features),
        "feature_columns": data.features,
        "target_column": config["target_column"],
        "group_column": config.get("group_column"),
        "id_column": config.get("id_column"),
        "seed": int(config.get("seed", 20260824)),
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {"numpy": np.__version__, "pandas": pd.__version__, "scikit_learn": sklearn.__version__},
        "leakage_guards": ["imputation_in_pipeline", "scaling_in_pipeline", "selection_in_pipeline", "nested_cv", "group_aware_if_configured"],
        "shap_status": shap_status,
    }
    (output_dir / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (output_dir / "effective_config.yml").write_text(yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8")
    write_model_card(output_dir / "model_card.md", config, summary, shap_status, external_status)
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    try:
        output = run(args.config)
    except (ContractError, FileNotFoundError, pd.errors.ParserError, yaml.YAMLError) as exc:
        print(f"[CONTRACT_ERROR] {exc}", file=sys.stderr)
        return 2
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

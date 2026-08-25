from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


SCRIPT = Path(__file__).parents[1] / "scripts" / "run_ml.py"
SPEC = importlib.util.spec_from_file_location("run_ml", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class MLPipelineTests(unittest.TestCase):
    def make_data(self, root: Path, *, leaked_id: bool = False) -> tuple[Path, Path]:
        rng = np.random.default_rng(7)
        rows = []
        for subject in range(20):
            outcome = subject % 2
            for repeat in range(2):
                rows.append({
                    "sample_id": "duplicate" if leaked_id else f"s{subject}_{repeat}",
                    "subject_id": f"subject_{subject}",
                    "outcome": "case" if outcome else "control",
                    "feature_1": outcome + rng.normal(scale=0.35),
                    "feature_2": rng.normal(),
                    "feature_3": rng.normal(),
                })
        data_path = root / "features.csv"
        pd.DataFrame(rows).to_csv(data_path, index=False)
        config = {
            "input_csv": str(data_path),
            "target_column": "outcome",
            "id_column": "sample_id",
            "group_column": "subject_id",
            "task": "classification",
            "positive_label": "case",
            "feature_columns": ["feature_1", "feature_2", "feature_3"],
            "model": "logistic",
            "outer_splits": 4,
            "inner_splits": 3,
            "feature_k": 2,
            "seed": 9,
            "output_dir": str(root / "result"),
            "run_shap": False,
        }
        config_path = root / "config.yml"
        config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
        return data_path, config_path

    def test_grouped_nested_cv_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, config = self.make_data(root)
            output = MODULE.run(config)
            expected = {
                "audit.json", "summary.json", "fold_metrics.csv", "oof_predictions.csv",
                "calibration.csv", "decision_curve.csv", "permutation_importance.csv",
                "final_model.joblib", "model_card.md", "effective_config.yml",
            }
            self.assertTrue(expected.issubset({p.name for p in output.iterdir()}))
            audit = json.loads((output / "audit.json").read_text())
            self.assertEqual(audit["n_rows"], 40)
            self.assertIn("nested_cv", audit["leakage_guards"])
            predictions = pd.read_csv(output / "oof_predictions.csv")
            self.assertEqual(len(predictions), 40)
            self.assertTrue((predictions.groupby("group")["outer_fold"].nunique() == 1).all())

    def test_duplicate_ids_fail_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, config_path = self.make_data(root, leaked_id=True)
            config = MODULE.load_config(config_path)
            with self.assertRaisesRegex(MODULE.ContractError, "不是唯一键"):
                MODULE.load_dataset(Path(config["input_csv"]), config)

    def test_regression_with_external_evaluation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = np.random.default_rng(11)
            train = pd.DataFrame({"sample_id": [f"d{i}" for i in range(36)], "x1": rng.normal(size=36), "x2": rng.normal(size=36)})
            train["outcome"] = 2 * train["x1"] - train["x2"] + rng.normal(scale=0.2, size=36)
            external = pd.DataFrame({"sample_id": [f"e{i}" for i in range(12)], "x1": rng.normal(size=12), "x2": rng.normal(size=12)})
            external["outcome"] = 2 * external["x1"] - external["x2"] + rng.normal(scale=0.2, size=12)
            train.to_csv(root / "train.csv", index=False)
            external.to_csv(root / "external.csv", index=False)
            config = {
                "input_csv": str(root / "train.csv"), "external_csv": str(root / "external.csv"),
                "target_column": "outcome", "id_column": "sample_id", "task": "regression",
                "feature_columns": ["x1", "x2"], "model": "ridge", "outer_splits": 3,
                "inner_splits": 3, "feature_k": "all", "seed": 3, "output_dir": str(root / "out"),
            }
            config_path = root / "config.yml"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            output = MODULE.run(config_path)
            self.assertEqual(len(pd.read_csv(output / "external_predictions.csv")), 12)
            summary = json.loads((output / "summary.json").read_text())
            self.assertIn("external", summary)


if __name__ == "__main__":
    unittest.main()

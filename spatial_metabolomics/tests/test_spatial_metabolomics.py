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


SCRIPT = Path(__file__).parents[1] / "scripts" / "run_spatial_metabolomics.py"
SPEC = importlib.util.spec_from_file_location("run_spatial_metabolomics", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class SpatialMetabolomicsTests(unittest.TestCase):
    def fixture(self, root: Path) -> Path:
        rng = np.random.default_rng(5)
        rows = []
        for sample in ["s1", "s2"]:
            for x in range(4):
                for y in range(4):
                    pixel = f"{sample}_{x}_{y}"
                    for mz in [100.001, 150.002, 200.003, 250.004]:
                        rows.append({"sample_id": sample, "pixel_id": pixel, "x": x, "y": y, "roi": "left" if x < 2 else "right", "mz": mz, "intensity": max(0, 20 + (30 if mz < 120 and x < 2 else 0) + rng.normal())})
        data = root / "spectra.csv"
        pd.DataFrame(rows).to_csv(data, index=False)
        config = {"mode": "long_table", "input": str(data), "output_dir": str(root / "out"), "mz_decimals": 3, "min_detected_features": 2, "min_tic_quantile": 0, "spatial_k": 3, "spatial_top_features": 4, "permutations": 9, "seed": 8}
        path = root / "config.yml"
        path.write_text(yaml.safe_dump(config), encoding="utf-8")
        return path

    def test_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = MODULE.run(self.fixture(Path(tmp)))
            self.assertTrue((out / "within_sample_spatial_stats.csv").is_file())
            audit = json.loads((out / "audit.json").read_text())
            self.assertEqual(audit["n_samples"], 2)
            self.assertEqual(audit["cross_sample_inference"], "not_performed")

    def test_false_level_one_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config_path = self.fixture(root)
            annotation = root / "annotations.csv"
            pd.DataFrame([{"mz_bin": 100.001, "annotation": "candidate", "identification_level": 1, "standard_confirmed": False, "source_version": "test"}]).to_csv(annotation, index=False)
            config = yaml.safe_load(config_path.read_text())
            config["annotation_table"] = str(annotation)
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
            with self.assertRaisesRegex(MODULE.ContractError, "standard_confirmed"):
                MODULE.run(config_path)


if __name__ == "__main__":
    unittest.main()

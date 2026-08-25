from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


SCRIPT = Path(__file__).parents[1] / "scripts" / "run_multiome.py"
SPEC = importlib.util.spec_from_file_location("run_multiome", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class MultiomeTests(unittest.TestCase):
    def fixture(self, root: Path, mismatch: bool = False) -> Path:
        rng = np.random.default_rng(4)
        cells = [f"cell_{i}" for i in range(24)]
        obs = pd.DataFrame({"sample_id": np.repeat(["s1", "s2", "s3", "s4"], 6), "group": np.repeat(["A", "A", "B", "B"], 6)}, index=cells)
        rna = ad.AnnData(sp.csr_matrix(rng.poisson(2, (24, 20))), obs=obs.copy(), var=pd.DataFrame(index=[f"g{i}" for i in range(20)]))
        atac_cells = cells[:-1] + (["other"] if mismatch else [cells[-1]])
        atac_obs = obs.copy()
        atac_obs.index = atac_cells
        atac = ad.AnnData(sp.csr_matrix(rng.poisson(1, (24, 30))), obs=atac_obs, var=pd.DataFrame(index=[f"p{i}" for i in range(30)]))
        rna.write_h5ad(root / "rna.h5ad")
        atac.write_h5ad(root / "atac.h5ad")
        config = {"rna_h5ad": str(root / "rna.h5ad"), "atac_h5ad": str(root / "atac.h5ad"), "paired": True, "sample_column": "sample_id", "group_column": "group", "output_dir": str(root / "out"), "lsi_components": 5, "write_h5mu": False}
        path = root / "config.yml"
        path.write_text(yaml.safe_dump(config), encoding="utf-8")
        return path

    def test_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = MODULE.run(self.fixture(Path(tmp)))
            self.assertTrue((out / "rna_qc.h5ad").is_file())
            self.assertTrue((out / "atac_qc_lsi.h5ad").is_file())
            audit = json.loads((out / "audit.json").read_text())
            self.assertEqual(audit["n_joint_after_qc"], 24)
            self.assertIn("do not create biological replicates", audit["inference_boundary"])

    def test_paired_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.fixture(Path(tmp), mismatch=True)
            with self.assertRaisesRegex(MODULE.ContractError, "禁止静默取交集"):
                MODULE.run(config)


if __name__ == "__main__":
    unittest.main()

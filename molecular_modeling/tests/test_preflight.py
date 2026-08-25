from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml


SCRIPT = Path(__file__).parents[1] / "scripts" / "preflight.py"
SPEC = importlib.util.spec_from_file_location("preflight", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


PDB = """HEADER    SYNTHETIC\nATOM      1  N   ALA A   1      11.104  13.207   9.347  1.00 20.00           N\nATOM      2  CA  ALA A   1      12.560  13.300   9.500  1.00 20.00           C\nEND\n"""


class PreflightTests(unittest.TestCase):
    def fixture(self, root: Path, replicas: int = 3) -> Path:
        (root / "receptor.pdb").write_text(PDB)
        ligand_rows = []
        for role in ["native", "positive_control", "decoy", "candidate"]:
            structure = root / f"{role}.sdf"
            structure.write_text(f"{role}\n  SCI\n\n  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 0 0 0 0 0\nM  V30 END CTAB\nM  END\n$$$$\n")
            ligand_rows.append({"ligand_id": role, "input_structure": structure.name, "protonation_state": "pH7.4", "formal_charge": "0", "charge_method": "AM1-BCC", "role": role})
        pd.DataFrame(ligand_rows).to_csv(root / "ligands.csv", index=False)
        config = {
            "receptor_pdb": str(root / "receptor.pdb"), "ligand_table": str(root / "ligands.csv"), "output_dir": str(root / "out"),
            "docking": {"engine": "vina", "box_center": [0, 0, 0], "box_size": [20, 20, 20], "seeds": [1, 2, 3], "redocking": True, "require_posebusters": True},
            "md": {"enabled": True, "replicas": replicas, "seeds": list(range(10, 10 + replicas)), "forcefield": "amber14", "water_model": "tip3pfb", "ligand_parameterization": "openff", "temperature_k": 300, "pressure_bar": 1, "production_ns": 10},
        }
        path = root / "config.yml"
        path.write_text(yaml.safe_dump(config), encoding="utf-8")
        return path

    def test_preflight_reports_blocked_without_dependencies(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = MODULE.run(self.fixture(Path(tmp)))
            audit = json.loads((out / "preflight.json").read_text())
            self.assertFalse(audit["execution_performed"])
            self.assertIn("docking score is not binding affinity", audit["claim_boundaries"])
            self.assertTrue((out / "execution_plan.md").is_file())

    def test_single_md_replica_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.fixture(Path(tmp), replicas=1)
            with self.assertRaisesRegex(MODULE.ContractError, "至少需要三个独立重复"):
                MODULE.run(config)


if __name__ == "__main__":
    unittest.main()

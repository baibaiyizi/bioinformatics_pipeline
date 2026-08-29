from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))

from sci_modeling.contract import ContractError, load_config, load_ligands
from sci_modeling.runtime import has_cuda_support, mdrun_used_gpu
from sci_modeling.workflow import _pose_rmsds, _write_amber_protein, _write_selected_pose, audit_receptor, md_production, preflight


PDB = """HEADER    SYNTHETIC HOLO
ATOM      1  N   ALA A   1      11.104  13.207   9.347  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.560  13.300   9.500  1.00 20.00           C
HETATM    3  C1  STI A 201      13.000  14.000  10.000  1.00 20.00           C
END
"""

SDF = """ligand
  SCI

  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 0 0 0 0 0
M  V30 END CTAB
M  END
$$$$
"""


class ModelingTests(unittest.TestCase):
    def test_gromacs_2026_cuda_version_field_is_detected(self):
        self.assertTrue(has_cuda_support("GROMACS version: 2026.3\nGPU support:         CUDA\n"))
        self.assertFalse(has_cuda_support("GROMACS version: 2026.3\nGPU support:         disabled\n"))
        self.assertTrue(mdrun_used_gpu("1 GPU selected for this run.\nMapping of GPU IDs"))
        self.assertFalse(mdrun_used_gpu("GPU support: CUDA\nNo compatible GPUs were detected"))

    def fixture(self, root: Path, *, docking_seeds=(1, 2, 3), md_replicas=3, protonation="pH7.4", roles=("native", "positive_control", "decoy", "candidate"), forcefield="ff14SB", production=None, metal=False, missing_residue=False) -> Path:
        pdb = PDB
        if metal:
            pdb = pdb.replace("END\n", "HETATM    4 ZN    ZN A 301      15.000  15.000  15.000  1.00 20.00          ZN\nEND\n")
        if missing_residue:
            pdb = pdb.replace("HEADER    SYNTHETIC HOLO\n", "HEADER    SYNTHETIC HOLO\nREMARK 465 MISSING RESIDUES\n")
        (root / "receptor.pdb").write_text(pdb)
        (root / "native.sdf").write_text(SDF)
        ligand_rows = []
        for index, role in enumerate(roles):
            structure = root / f"{role}_{index}.sdf"
            structure.write_text(SDF)
            ligand_rows.append({"ligand_id": f"{role}_{index}", "input_structure": structure.name, "role": role, "protonation_state": protonation, "formal_charge": "0", "charge_method": "AM1-BCC"})
        pd.DataFrame(ligand_rows).to_csv(root / "ligands.tsv", sep="\t", index=False)
        selected = ligand_rows[0]["ligand_id"] if ligand_rows else "none"
        config = {
            "project_id": "fixture",
            "receptor": {"input_pdb": str(root / "receptor.pdb"), "source": {"database": "RCSB_PDB", "id": "1IEP"}, "protein_chains": ["A"], "covalent_ligand": False, "protonation_method": "reduce_hbond_network"},
            "pocket": {"basis": "co_crystal_ligand", "native_ligand_resname": "STI", "native_reference_sdf": str(root / "native.sdf"), "center": [1, 2, 3], "size": [20, 20, 20]},
            "ligand_table": str(root / "ligands.tsv"), "output_dir": str(root / "out"),
            "docking": {"engine": "vina", "seeds": list(docking_seeds), "exhaustiveness": 8, "num_modes": 5, "require_posebusters": True, "redocking": {"enabled": True, "rmsd_threshold_angstrom": 2.0, "minimum_passing_seeds": 2}},
            "md": {"engine": "gromacs", "ligand_id": selected, "replicas": md_replicas, "seeds": list(range(10, 10 + md_replicas)), "parameterization": {"protein_forcefield": forcefield, "ligand_forcefield": "GAFF2", "charge_method": "AM1-BCC", "water_model": "TIP3P"}, "smoke": {"steps": 10}, "production": production if production is not None else {"enabled": False}},
            "analysis": {"ligand_selection": "resname LIG", "contact_cutoff_angstrom": 4.5, "key_distances": []},
        }
        path = root / "config.yml"
        path.write_text(yaml.safe_dump(config))
        return path

    def test_preflight_records_resource_limited_not_ready(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.fixture(Path(tmp))
            out = preflight(path)
            payload = json.loads((out / "preflight.json").read_text())
            production = next(row for row in payload["stages"] if row["stage"] == "md-production")
            self.assertEqual("resource-limited", production["status"])
            self.assertIn("short trajectory", payload["claim_boundaries"][1])

    def test_missing_native_or_controls_rejected(self):
        for roles in [("positive_control", "decoy", "candidate"), ("native", "candidate")]:
            with self.subTest(roles=roles), tempfile.TemporaryDirectory() as tmp:
                path = self.fixture(Path(tmp), roles=roles)
                config = load_config(path)
                with self.assertRaises(ContractError):
                    load_ligands(path, config)

    def test_three_docking_seeds_and_replicas_are_hard_gates(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ContractError, "三个不同 seed"):
                load_config(self.fixture(Path(tmp), docking_seeds=(1, 2)))
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ContractError, "三个独立 replica"):
                load_config(self.fixture(Path(tmp), md_replicas=2))

    def test_unclear_protonation_and_forcefield_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.fixture(Path(tmp), protonation="unknown")
            with self.assertRaisesRegex(ContractError, "质子化状态不明确"):
                load_ligands(path, load_config(path))
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ContractError, "ff14SB"):
                load_config(self.fixture(Path(tmp), forcefield="amber99sb"))

    def test_metal_blocks_automatic_parameterization(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.fixture(Path(tmp), metal=True)
            config = load_config(path)
            audit = audit_receptor(Path(config["receptor"]["input_pdb"]), "STI", ["A"])
            self.assertIn("metal:ZN", audit["automatic_parameterization_blockers"])

    def test_unresolved_missing_residue_review_blocks_dock(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = preflight(self.fixture(Path(tmp), missing_residue=True))
            stages = pd.read_csv(out / "stage_status.tsv", sep="\t").set_index("stage")
            self.assertEqual("review_required", stages.loc["preflight", "status"])
            self.assertEqual("blocked", stages.loc["dock", "status"])

    def test_production_requires_explicit_duration_and_pose(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ContractError, "duration_ns"):
                load_config(self.fixture(Path(tmp), production={"enabled": True, "approved_pose": "pose.sdf"}))
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ContractError, "approved_pose"):
                load_config(self.fixture(Path(tmp), production={"enabled": True, "duration_ns": 10}))

    def test_gpu_unavailable_forbids_production(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pose = root / "pose.sdf"
            pose.write_text(SDF)
            path = self.fixture(root, production={"enabled": True, "duration_ns": 1, "approved_pose": str(pose)})
            with patch("sci_modeling.workflow.dependencies", return_value={"gpu_ready": False}):
                with self.assertRaisesRegex(ContractError, "禁止生产 MD"):
                    md_production(path)

    def test_pose_selection_writes_one_explicit_conformer(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
            AllChem.EmbedMolecule(molecule, randomSeed=7)
            reference = root / "reference.sdf"
            poses = root / "poses.sdf"
            writer = Chem.SDWriter(str(reference))
            writer.write(molecule)
            writer.close()
            writer = Chem.SDWriter(str(poses))
            writer.write(molecule)
            shifted = Chem.Mol(molecule)
            shifted.GetConformer().SetAtomPosition(0, (5.0, 5.0, 5.0))
            writer.write(shifted)
            writer.close()

            values = _pose_rmsds(poses, reference)
            self.assertEqual(2, len(values))
            self.assertAlmostEqual(0.0, values[0], places=6)
            selected = root / "selected.sdf"
            _write_selected_pose(poses, 1, selected)
            selected_molecules = [mol for mol in Chem.SDMolSupplier(str(selected), removeHs=False) if mol]
            self.assertEqual(1, len(selected_molecules))

    def test_reduce_histidine_hydrogens_become_explicit_amber_states(self):
        reduced = """ATOM      1  ND1 HIS A  10      11.104  13.207   9.347  1.00 20.00           N
ATOM      2  HD1 HIS A  10      11.204  13.307   9.447  1.00 20.00           H
ATOM      3  NE2 HIS A  10      12.560  13.300   9.500  1.00 20.00           N
END
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, output = root / "reduced.pdb", root / "amber.pdb"
            source.write_text(reduced)
            assignments = _write_amber_protein(source, output, ["A"])
            self.assertEqual({"A:10": "HID"}, assignments)
            self.assertIn("HID A  10", output.read_text())


if __name__ == "__main__":
    unittest.main()

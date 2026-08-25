from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml


SCRIPT = Path(__file__).parents[1] / "scripts" / "run_network.py"
SPEC = importlib.util.spec_from_file_location("run_network", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class NetworkTests(unittest.TestCase):
    def fixture(self, root: Path, duplicate_evidence: bool = False) -> Path:
        nodes = pd.DataFrame([
            {"node_id": "drug:A", "node_type": "compound", "canonical_id": "CHEM:A", "species": ""},
            {"node_id": "gene:X", "node_type": "gene", "canonical_id": "ENSMUSG:X", "species": "10090"},
            {"node_id": "gene:Y", "node_type": "gene", "canonical_id": "ENSMUSG:Y", "species": "10090"},
            {"node_id": "disease:D", "node_type": "disease", "canonical_id": "MONDO:D", "species": ""},
        ])
        evidence_ids = ["e1", "e2", "e3"]
        edge_rows = [
            {"source": "drug:A", "target": "gene:X", "relation": "targets", "evidence_source": "curated_db", "evidence_id": evidence_ids[0], "source_version": "2026-01", "evidence_tier": "curated", "confidence": 0.9},
            {"source": "gene:X", "target": "gene:Y", "relation": "interacts", "evidence_source": "STRING", "evidence_id": evidence_ids[1], "source_version": "12.0", "evidence_tier": "computational", "confidence": 0.8},
            {"source": "gene:Y", "target": "disease:D", "relation": "associated", "evidence_source": "OpenTargets", "evidence_id": evidence_ids[2], "source_version": "2026-Q2", "evidence_tier": "genetic", "confidence": 0.95},
        ]
        if duplicate_evidence:
            edge_rows.append(dict(edge_rows[0]))
        edges = pd.DataFrame(edge_rows)
        nodes.to_csv(root / "nodes.csv", index=False)
        edges.to_csv(root / "edges.csv", index=False)
        config = {"nodes": str(root / "nodes.csv"), "edges": str(root / "edges.csv"), "output_dir": str(root / "out"), "robustness_repeats": 5, "edge_dropout": 0.1, "seed": 2}
        path = root / "config.yml"
        path.write_text(yaml.safe_dump(config), encoding="utf-8")
        return path

    def test_smoke_and_claim_boundary(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = MODULE.run(self.fixture(Path(tmp)))
            self.assertTrue((out / "evidence_graph.graphml").is_file())
            audit = json.loads((out / "audit.json").read_text())
            self.assertIn("not efficacy", audit["claim_boundary"])
            cards = (out / "candidate_evidence_cards.md").read_text()
            self.assertIn("不等于疗效", cards)

    def test_duplicate_evidence_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self.fixture(Path(tmp), duplicate_evidence=True)
            with self.assertRaisesRegex(MODULE.ContractError, "禁止重复计权"):
                MODULE.run(config)


if __name__ == "__main__":
    unittest.main()

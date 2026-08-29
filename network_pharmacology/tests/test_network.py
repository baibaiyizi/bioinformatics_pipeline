from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import httpx
import pandas as pd
import yaml

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))

from sci_network.adapters import ChEMBL, OpenTargets, PubChem
from sci_network.analysis import build_multidigraph, rank_candidates, reactome_ora
from sci_network.contract import ContractError, load_config
from sci_network.http import AuditClient, RetrievalError
from sci_network.pipeline import _one_mouse_ortholog, run


class NetworkTests(unittest.TestCase):
    def test_opentargets_pagination_and_count_reconciliation(self):
        calls = []

        def transport(**kwargs):
            page = kwargs["json"]["variables"]["page"]["index"]
            calls.append(page)
            rows = [{"target": {"id": f"ENSG{page}", "approvedSymbol": f"G{page}"}, "score": 0.9 - page / 10, "datasourceScores": []}]
            return 200, {"x-api-version": "v4"}, {"data": {"disease": {"associatedTargets": {"count": 2, "rows": rows}}}}

        with tempfile.TemporaryDirectory() as tmp:
            client = AuditClient(Path(tmp), transport=transport)
            rows = OpenTargets(client).disease_targets("MONDO_1", page_size=1, max_targets=10)
            self.assertEqual([0, 1], calls)
            self.assertEqual(2, len(rows))
            client.write_manifest()
            manifest = [json.loads(line) for line in (Path(tmp) / "provenance/requests.jsonl").read_text().splitlines()]
            self.assertEqual({0, 1}, {row["page"]["index"] for row in manifest})

    def test_opentargets_drug_targets_use_current_mechanism_schema(self):
        def transport(**kwargs):
            self.assertNotIn("linkedTargets", kwargs["json"]["query"])
            rows = [
                {"mechanismOfAction": "ABL inhibitor", "actionType": "INHIBITOR", "targets": [{"id": "ENSG1", "approvedSymbol": "ABL1"}]},
                {"mechanismOfAction": "BCR-ABL inhibitor", "actionType": "INHIBITOR", "targets": [{"id": "ENSG1", "approvedSymbol": "ABL1"}]},
            ]
            return 200, {}, {"data": {"drug": {"id": "CHEMBL1", "name": "X", "mechanismsOfAction": {"rows": rows}}}}

        with tempfile.TemporaryDirectory() as tmp:
            rows = OpenTargets(AuditClient(Path(tmp), transport=transport)).drug_targets("CHEMBL1")
            self.assertEqual(1, len(rows))
            self.assertEqual(["ABL inhibitor", "BCR-ABL inhibitor"], rows[0]["mechanisms"])

    def test_pagination_mismatch_and_api_failure_are_visible(self):
        def short_page(**kwargs):
            return 200, {}, {"data": {"disease": {"associatedTargets": {"count": 2, "rows": []}}}}

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RetrievalError, "分页不完整"):
                OpenTargets(AuditClient(Path(tmp), transport=short_page)).disease_targets("MONDO_1", 10, 10)

        def failed(**kwargs):
            return 500, {}, {"error": "service unavailable"}

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RetrievalError, "HTTP 500"):
                AuditClient(Path(tmp), transport=failed).request_json("source", "id", "GET", "https://example.test")

    def test_transient_503_is_retried_without_hiding_persistent_failure(self):
        calls = 0

        def transient(**kwargs):
            nonlocal calls
            calls += 1
            if calls < 3:
                return 503, {"retry-after": "0"}, {"error": "temporarily unavailable"}
            return 200, {}, {"ok": True}

        with tempfile.TemporaryDirectory() as tmp:
            payload, _ = AuditClient(Path(tmp), transport=transient).request_json("source", "id", "GET", "https://example.test")
            self.assertEqual({"ok": True}, payload)
            self.assertEqual(3, calls)

        def persistent(**kwargs):
            return 503, {"retry-after": "0"}, {"error": "still unavailable"}

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RetrievalError, "HTTP 503"):
                AuditClient(Path(tmp), transport=persistent).request_json("source", "id", "GET", "https://example.test")

    def test_transient_read_timeout_is_retried(self):
        calls = 0

        def transient(**kwargs):
            nonlocal calls
            calls += 1
            if calls == 1:
                raise httpx.ReadTimeout("temporary timeout")
            return 200, {}, {"ok": True}

        with tempfile.TemporaryDirectory() as tmp:
            payload, _ = AuditClient(Path(tmp), transport=transient).request_json("source", "id", "GET", "https://example.test")
            self.assertEqual({"ok": True}, payload)
            self.assertEqual(2, calls)

    def test_no_result_is_explicit_and_not_invented(self):
        def no_result(**kwargs):
            return 200, {}, {"data": {"disease": {"associatedTargets": {"count": 0, "rows": []}}}}

        with tempfile.TemporaryDirectory() as tmp:
            rows = OpenTargets(AuditClient(Path(tmp), transport=no_result)).disease_targets("MONDO_1", 10, 10)
            self.assertEqual([], rows)

    def test_http_204_json_endpoint_is_audited_as_empty_result(self):
        def no_content(**kwargs):
            return 204, {}, {}

        with tempfile.TemporaryDirectory() as tmp:
            client = AuditClient(Path(tmp), transport=no_content)
            payload, _ = client.request_json("rcsb", "no_hit", "POST", "https://example.test")
            self.assertEqual({}, payload)
            self.assertEqual(204, client.records[0]["status"])

    def test_pubchem_one_to_many_identity_rejected(self):
        def ambiguous(**kwargs):
            return 200, {}, {"PropertyTable": {"Properties": [{"CID": 1}, {"CID": 2}]}}

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RetrievalError, "非一对一"):
                PubChem(AuditClient(Path(tmp), transport=ambiguous)).identity("inchikey", "AAAA")

    def test_chembl_does_not_silently_truncate(self):
        def too_many(**kwargs):
            return 200, {}, {"page_meta": {"total_count": 1001}, "activities": []}

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RetrievalError, "禁止静默截断"):
                ChEMBL(AuditClient(Path(tmp), transport=too_many)).paged(
                    "activity", "activities", {"molecule_chembl_id": "CHEMBL1"}, max_records=1000
                )

    def test_chembl_uniprot_mapping_requires_one_direct_species_single_protein(self):
        def one_target(**kwargs):
            self.assertEqual("P00519", kwargs["params"]["target_components__accession"])
            self.assertEqual("SINGLE PROTEIN", kwargs["params"]["target_type"])
            row = {"target_chembl_id": "CHEMBL1862", "target_type": "SINGLE PROTEIN", "tax_id": 9606, "target_components": [{"accession": "P00519"}]}
            return 200, {}, {"page_meta": {"total_count": 1}, "targets": [row]}

        with tempfile.TemporaryDirectory() as tmp:
            target = ChEMBL(AuditClient(Path(tmp), transport=one_target)).single_protein_target("P00519", 9606)
            self.assertEqual("CHEMBL1862", target["target_chembl_id"])

        def ambiguous(**kwargs):
            rows = [
                {"target_chembl_id": f"CHEMBL{index}", "target_type": "SINGLE PROTEIN", "tax_id": 9606, "target_components": [{"accession": "P00519"}]}
                for index in (1, 2)
            ]
            return 200, {}, {"page_meta": {"total_count": 2}, "targets": rows}

        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(RetrievalError, "一对多"):
                ChEMBL(AuditClient(Path(tmp), transport=ambiguous)).single_protein_target("P00519", 9606)

    def test_duplicate_evidence_and_species_boundaries(self):
        nodes = pd.DataFrame([
            {"node_id": "d:D", "node_type": "disease", "canonical_id": "D", "species_taxid": "", "label": "D"},
            {"node_id": "t:T", "node_type": "target", "canonical_id": "T", "species_taxid": 9606, "label": "T"},
        ])
        edge = {"source": "d:D", "target": "t:T", "relation": "associated", "database": "Open Targets", "database_version": "v4", "evidence_id": "e1", "evidence_layer": "direct_species", "species_taxid": 9606, "score": 0.8, "lane": "disease_association", "network_type": "", "source_symbol": "", "target_symbol": "T"}
        with self.assertRaisesRegex(ContractError, "重复"):
            build_multidigraph(nodes, pd.DataFrame([edge, edge]))

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.yml"
            config = {"project_id": "x", "species_taxid": 10090, "disease_ids": [], "compounds": [], "target_seeds": [{"ensembl_gene_id": "ENSG0001", "evidence_source": "user", "evidence_id": "e"}], "output_dir": "out"}
            path.write_text(yaml.safe_dump(config))
            with self.assertRaisesRegex(ContractError, "不兼容"):
                load_config(path)

        class OneToMany:
            @staticmethod
            def mouse_orthologs(human_id):
                return [{"type": "ortholog_one2many", "target": {"id": "ENSMUSG1"}}]

        with self.assertRaisesRegex(ContractError, "不是一对一"):
            _one_mouse_ortholog(OneToMany(), "ENSG1")

    def test_reactome_requires_explicit_species_background(self):
        with self.assertRaisesRegex(ContractError, "背景宇宙"):
            reactome_ora({"P1"}, {}, 9606)
        result = reactome_ora({"P1"}, {"P1": [{"stId": "R-HSA-1", "displayName": "human"}, {"stId": "R-MMU-2", "displayName": "mouse"}], "P2": [{"stId": "R-HSA-1", "displayName": "human"}]}, 9606)
        self.assertEqual(["R-HSA-1"], result["pathway_id"].tolist())
        self.assertTrue(result["fdr_bh"].between(0, 1).all())

    def test_less_than_two_independent_lanes_is_not_ranked(self):
        targets = pd.DataFrame([{"target_id": "T1", "symbol": "G1", "evidence_layer": "direct_species"}])
        lanes = pd.DataFrame([{"target_id": "T1", "lane": "disease_association", "raw_score": 0.9}])
        ranked, _ = rank_candidates(targets, lanes)
        self.assertFalse(bool(ranked.loc[0, "eligible"]))
        self.assertTrue(pd.isna(ranked.loc[0, "priority_score"]))

        targets = pd.DataFrame([
            {"target_id": "T1", "symbol": "G1", "evidence_layer": "direct_species"},
            {"target_id": "T2", "symbol": "G2", "evidence_layer": "direct_species"},
        ])
        lanes = pd.DataFrame([
            {"target_id": "T1", "lane": "disease_association", "raw_score": 0.9},
            {"target_id": "T1", "lane": "structure_availability", "raw_score": 1.0},
            {"target_id": "T2", "lane": "disease_association", "raw_score": 0.8},
        ])
        ranked, _ = rank_candidates(targets, lanes)
        eligibility = ranked.set_index("target_id")["eligible"].to_dict()
        self.assertEqual({"T1": True, "T2": False}, eligibility)

    def test_end_to_end_fixture_keeps_string_network_types_separate(self):
        lookup = {"ENSG00000000001": "G1", "ENSG00000000002": "G2"}

        def transport(**kwargs):
            url = kwargs["url"]
            if "biomart/martservice" in url:
                if kwargs["params"].get("type") == "registry":
                    return 200, {}, '<MartRegistry><MartURLLocation database="ensembl_mart_116" /></MartRegistry>'
                query = kwargs["params"]["query"]
                target = next(target for target in lookup if target in query)
                if 'name="external_gene_name"' in query:
                    return 200, {}, f"{target}\t{lookup[target]}\n"
                if 'name="uniprotswissprot"' in query:
                    accession = "P00001" if target.endswith("1") else "P00002"
                    return 200, {}, f"{target}\t{accession}\n"
                raise AssertionError(query)
            if url.endswith("/api/json/version"):
                return 200, {}, [{"string_version": "12.0"}]
            if "rest.uniprot.org" in url:
                accession = url.rsplit("/", 1)[-1]
                return 200, {"x-uniprot-release": "fixture"}, {"primaryAccession": accession, "organism": {"taxonId": 9606}}
            if "string-db.org" in url:
                network_type = kwargs["params"]["network_type"]
                return 200, {}, [{"preferredName_A": "G1", "preferredName_B": "G2", "stringId_A": "9606.A", "stringId_B": "9606.B", "score": 0.9 if network_type == "physical" else 0.8}]
            if "rcsbsearch" in url:
                return 200, {}, {"total_count": 1, "result_set": [{"identifier": "1ABC"}]}
            raise AssertionError(url)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = {"project_id": "fixture", "species_taxid": 9606, "disease_ids": [], "compounds": [], "target_seeds": [
                {"ensembl_gene_id": "ENSG00000000001", "evidence_source": "user", "evidence_id": "s1"},
                {"ensembl_gene_id": "ENSG00000000002", "evidence_source": "user", "evidence_id": "s2"},
            ], "output_dir": "out", "retrieval": {"max_targets": 10, "max_structure_checks": 2, "page_size": 10, "timeout_seconds": 1, "string_required_score": 400}, "analysis": {"string_thresholds": [400, 700, 900], "reactome_background_uniprot": []}}
            path = root / "config.yml"
            path.write_text(yaml.safe_dump(config))
            out = run(path, transport=transport)
            edges = pd.read_csv(out / "evidence_edges.tsv", sep="\t")
            string_edges = edges[edges["database"] == "STRING"]
            self.assertEqual({"physical", "functional"}, set(string_edges["network_type"]))
            ranked = pd.read_csv(out / "candidate_priority.tsv", sep="\t")
            self.assertTrue(ranked["eligible"].all())
            self.assertTrue(pd.read_csv(out / "modeling_handoff.tsv", sep="\t").empty)
            self.assertTrue(pd.read_csv(out / "chembl_activity.tsv", sep="\t").empty)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import os
import pathlib
import tomllib
import unittest


AGENT_DIR = pathlib.Path(os.environ.get("SCI_AGENT_DIR", "/home/h1028/.codex/agents"))
WORKSPACE = pathlib.Path(__file__).parents[2]

ANALYSIS_ROLES = {
    "sci_analyst",
    "sci_transcriptomics",
    "sci_singlecell",
    "sci_spatialomics",
    "sci_epigenomics",
    "sci_metabolomics",
    "sci_multiomics",
    "sci_ml",
    "sci_network_pharmacology",
    "sci_molecular_modeling",
}


def configs():
    parsed = {}
    for path in sorted(AGENT_DIR.glob("sci*.toml")):
        data = tomllib.loads(path.read_text(encoding="utf-8"))
        parsed[data["name"]] = (path, data)
    return parsed


class SciAgentConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.agents = configs()

    def test_names_unique_and_complete(self):
        expected = ANALYSIS_ROLES | {"sci", "sci_evidence", "sci_writer", "sci_reviewer"}
        self.assertEqual(set(self.agents), expected)
        self.assertEqual(len(self.agents), 14)

    def test_analysis_role_runtime_contract(self):
        for name in ANALYSIS_ROLES:
            _, config = self.agents[name]
            self.assertEqual(config["model"], "gpt-5.6-sol", name)
            self.assertEqual(config["model_reasoning_effort"], "xhigh", name)
            self.assertEqual(config["sandbox_mode"], "workspace-write", name)
            instructions = config["developer_instructions"]
            self.assertIn("sci_analysis_contract.md", instructions, name)
            self.assertNotIn("spawn_agent", instructions, name)

    def test_reviewer_is_read_only(self):
        self.assertEqual(self.agents["sci_reviewer"][1]["sandbox_mode"], "read-only")

    def test_router_names_all_specialists_and_limits_concurrency(self):
        instructions = self.agents["sci"][1]["developer_instructions"]
        for role in ANALYSIS_ROLES:
            self.assertIn(role, instructions)
        self.assertIn("最多同时运行三个子 agent", instructions)
        self.assertIn("只有 `sci` 可以创建子 agent", instructions)

    def test_route_scenarios(self):
        router = self.agents["sci"][1]["developer_instructions"]
        cases = {
            "bulk/total/small RNA": ["sci_transcriptomics", "sci_analyst"],
            "scRNA": ["sci_singlecell", "pseudobulk"],
            "空间代谢组": ["sci_metabolomics", "sci_spatialomics"],
            "bulk ATAC": ["sci_epigenomics"],
            "跨组学整合": ["sci_multiomics"],
            "分类、回归": ["sci_ml"],
            "网络药理": ["sci_evidence", "sci_network_pharmacology"],
            "对接、MD": ["sci_molecular_modeling"],
        }
        for scenario, tokens in cases.items():
            self.assertIn(scenario, router)
            for token in tokens:
                self.assertIn(token, router)

    def test_negative_scientific_boundaries_are_explicit(self):
        text = "\n".join(data["developer_instructions"] for _, data in self.agents.values())
        shared = (AGENT_DIR / "sci_analysis_contract.md").read_text(encoding="utf-8")
        combined = text + "\n" + shared
        for forbidden_shortcut in [
            "细胞、spot、pixel",
            "全数据筛选特征",
            "网络中心性",
            "docking score",
            "单条短轨迹",
            "最高等级确定鉴定",
        ]:
            self.assertIn(forbidden_shortcut, combined)

    def test_workspace_documentation_names_new_roles(self):
        docs = (WORKSPACE / "skills" / "skill.md").read_text(encoding="utf-8")
        for role in ANALYSIS_ROLES:
            self.assertIn(role, docs)


if __name__ == "__main__":
    unittest.main()

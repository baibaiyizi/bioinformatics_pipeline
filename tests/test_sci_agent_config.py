from __future__ import annotations

import os
import pathlib
import tomllib
import unittest


AGENT_DIR = pathlib.Path(os.environ.get("SCI_AGENT_DIR", "/home/h1028/.codex/agents"))
WORKSPACE = pathlib.Path(__file__).parents[2]
LEARNING_ROOT = WORKSPACE / "sci_learning"

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

    def test_learning_loop_is_wired_into_router_and_evidence_agent(self):
        router = self.agents["sci"][1]["developer_instructions"]
        evidence = self.agents["sci_evidence"][1]["developer_instructions"]
        shared = (AGENT_DIR / "sci_analysis_contract.md").read_text(encoding="utf-8")

        for token in [
            "/home/h1028/workspace/sci_learning/",
            "knowledge_refs",
            "adopt",
            "adapt",
            "background_only",
            "reject",
            "knowledge_feedback",
        ]:
            self.assertIn(token, router)
            self.assertIn(token, shared)

        for token in [
            "new_to_library",
            "known",
            "uncertain",
            "abstract_only",
            "搜索片段",
            "学习卡不能作为正式引文",
        ]:
            self.assertIn(token, evidence)

    def test_learning_workspace_contract(self):
        required_files = [
            LEARNING_ROOT / "README.md",
            LEARNING_ROOT / "index.md",
            LEARNING_ROOT / "methods_index.md",
            LEARNING_ROOT / "cards" / "_template.md",
            LEARNING_ROOT / "repositories" / "_template.md",
            LEARNING_ROOT / "syntheses" / "_template.md",
        ]
        for path in required_files:
            self.assertTrue(path.is_file(), path)

        card = (LEARNING_ROOT / "cards" / "_template.md").read_text(encoding="utf-8")
        for token in [
            "learning_id",
            "full_text | abstract_only | web_page",
            "## 论证链",
            "## 值得学习的方法",
            "## 可迁移思路",
            "## 局限、替代解释和边界",
            "## 应用历史与实践反馈",
        ]:
            self.assertIn(token, card)

        methods = (LEARNING_ROOT / "methods_index.md").read_text(encoding="utf-8")
        for token in [
            "new_to_library",
            "known",
            "uncertain",
            "supported",
            "unsupported",
            "inconclusive",
            "not_run",
        ]:
            self.assertIn(token, methods)

        repository = (LEARNING_ROOT / "repositories" / "_template.md").read_text(
            encoding="utf-8"
        )
        for token in [
            "code_repository",
            "forge: \"github | gitlab\"",
            "canonical_url",
            "repository_id",
            "checked_commit",
            "metadata_only | readme_docs | source_audit | executed",
            "## 架构、入口与数据流",
            "## 依赖、测试与复现",
            "## 许可证、维护与风险",
            "## 关联论文和学习记录",
            "## 应用历史与实践反馈",
        ]:
            self.assertIn(token, repository)

    def test_code_repository_learning_boundaries(self):
        router = self.agents["sci"][1]["developer_instructions"]
        evidence = self.agents["sci_evidence"][1]["developer_instructions"]
        shared = (AGENT_DIR / "sci_analysis_contract.md").read_text(encoding="utf-8")
        combined = "\n".join([router, evidence, shared])
        for token in [
            "github-owner-repository",
            "gitlab-namespace-repository",
            "checked_commit",
            "metadata_only",
            "readme_docs",
            "source_audit",
            "executed",
            "许可证",
            "代码仓库本身不自动构成科学有效性证据",
        ]:
            self.assertIn(token, combined)

    def test_network_pharmacology_closed_loop_contract(self):
        instructions = self.agents["sci_network_pharmacology"][1][
            "developer_instructions"
        ]
        for token in [
            "Open Targets",
            "ChEMBL",
            "PubChem",
            "Reactome",
            "Ensembl",
            "MultiDiGraph",
            "human_ortholog",
            "direct_species",
            "STRING physical",
            "显式背景宇宙",
            "至少需要两条独立 lane",
            "留一证据层",
        ]:
            self.assertIn(token, instructions)
        self.assertNotIn("`primekg`", instructions)

    def test_vina_gromacs_closed_loop_contract(self):
        instructions = self.agents["sci_molecular_modeling"][1][
            "developer_instructions"
        ]
        for token in [
            "AutoDock Vina",
            "Meeko",
            "PoseBusters",
            "GROMACS 2026.3",
            "ff14SB",
            "GAFF2/AM1-BCC",
            "TIP3P",
            "ParmEd",
            "至少三个随机 seed",
            "至少两个 seed",
            "≤2 Å",
            "三个独立 seed",
            "resource-limited",
        ]:
            self.assertIn(token, instructions)
        self.assertIn("不安装或调用 DiffDock、OpenMM", instructions)


if __name__ == "__main__":
    unittest.main()

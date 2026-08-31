from __future__ import annotations

import os
import pathlib
import re
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

    def test_review_writing_and_visualization_learning_contracts(self):
        router = self.agents["sci"][1]["developer_instructions"]
        evidence = self.agents["sci_evidence"][1]["developer_instructions"]
        writer = self.agents["sci_writer"][1]["developer_instructions"]
        reviewer = self.agents["sci_reviewer"][1]["developer_instructions"]

        writing_contract = AGENT_DIR / "sci_review_writing_contract.md"
        figure_contract = AGENT_DIR / "sci_review_figure_contract.md"
        self.assertTrue(writing_contract.is_file())
        self.assertTrue(figure_contract.is_file())

        for token in ["writing", "visualization", "阶段 1–5"]:
            self.assertIn(token, router)
        for instructions in (writer, reviewer):
            self.assertIn("knowledge_refs", instructions)
            self.assertIn("knowledge_feedback", instructions)
            self.assertIn("sci_review_writing_contract.md", instructions)
            self.assertIn("sci_review_figure_contract.md", instructions)

        for token in [
            "nearest-review",
            "frozen baseline + dated delta + semantic propagation",
            "逐篇罗列",
            "proposal",
        ]:
            self.assertIn(token, writing_contract.read_text(encoding="utf-8"))

        for token in [
            "一图一主张",
            "MediaBox/CropBox",
            "FontFile",
            "deuteranopia",
            "程序化矢量",
            "target-review deconstruction",
            "orienting_overview",
            "review-specific boundary",
            "arrow paraphrase",
            "非定量综述示意图",
            "`basis` 只允许 `official | observed | inference`",
            "narrative_salience_encoding",
            "conceptual_schematic",
            "quantitative_figure",
            "实时目标期刊 > 实时出版社/平台",
            "semantic_asset_class",
            "qa_class",
            "evidence_origin",
            "claim_status",
            "relation_type",
            "figure_rule_ledger.tsv",
        ]:
            self.assertIn(token, figure_contract.read_text(encoding="utf-8"))

        self.assertIn("不得直接修改全局卡片或索引", evidence)
        self.assertIn("不得更新全局 `index.md` 或 `methods_index.md`", evidence)

    def test_learning_templates_accept_writing_and_visualization(self):
        readme = (LEARNING_ROOT / "README.md").read_text(encoding="utf-8")
        synthesis = (LEARNING_ROOT / "syntheses" / "_template.md").read_text(
            encoding="utf-8"
        )
        for token in ["writing", "visualization"]:
            self.assertIn(token, readme)
            self.assertIn(token, synthesis)

    def test_review_method_learning_ids_are_unique_and_syntheses_resolve(self):
        index = (LEARNING_ROOT / "index.md").read_text(encoding="utf-8")
        index_ids = re.findall(r"^\| ([^| ]+) \|", index, flags=re.M)
        self.assertEqual(len(index_ids), len(set(index_ids)))

        card_ids = []
        for path in sorted((LEARNING_ROOT / "cards").glob("*.md")):
            if path.name == "_template.md":
                continue
            match = re.search(
                r'^learning_id:\s*["\']?([^"\'\n]+)',
                path.read_text(encoding="utf-8"),
                flags=re.M,
            )
            self.assertIsNotNone(match, path)
            card_ids.append(match.group(1).strip())
        self.assertEqual(len(card_ids), len(set(card_ids)))
        self.assertTrue(set(card_ids).issubset(set(index_ids)))

        for name in [
            "2026-08-30-review-positioning-evidence-writing.md",
            "2026-08-30-review-figure-semantics-production-qa.md",
        ]:
            path = LEARNING_ROOT / "syntheses" / name
            self.assertTrue(path.is_file(), path)
            text = path.read_text(encoding="utf-8")
            source_section = text.split("## 实际精读来源", 1)[1].split(
                "## 应用反馈", 1
            )[0]
            referenced_ids = re.findall(r"`([^`]+)`", source_section)
            self.assertTrue(referenced_ids, path)
            for learning_id in referenced_ids:
                self.assertTrue(
                    (LEARNING_ROOT / "cards" / f"{learning_id}.md").is_file(),
                    f"{path}: unresolved learning ID {learning_id}",
                )

        rule_ledger = LEARNING_ROOT / "figure_rule_ledger.tsv"
        self.assertTrue(rule_ledger.is_file())
        ledger_lines = rule_ledger.read_text(encoding="utf-8").splitlines()
        self.assertEqual(
            ledger_lines[0].split("\t"),
            [
                "decision_id",
                "figure_id",
                "decision",
                "basis",
                "source_id",
                "locator",
                "scope",
                "review_date",
            ],
        )
        for line in ledger_lines[1:]:
            fields = line.split("\t")
            self.assertEqual(len(fields), 8, line)
            self.assertIn(fields[3], {"official", "observed", "inference"})

        review_figure_cards = {
            "url-www-nature-com-documents-natrev-artworkguide-pdf": "official",
            "url-academic-oup-com-humupd-pages-general": "official",
            "url-researcher-resources-acs-org-publish-author-guidelines": "official",
            "url-crosstalk-cell-com-hubfs-files-ga-guide-pdf": "official",
            "url-www-elsevier-com-about-policies-and-standards-generative-ai-policies-for-journals": "official",
            "url-www-elsevier-com-researcher-author-tools-and-resources-graphical-abstract": "official",
            "pmid-38300895": "observed",
            "pmid-42277276": "observed",
            "pmid-40580487": "observed",
            "pmid-38548833": "observed",
        }
        for learning_id, evidence_class in review_figure_cards.items():
            card_text = (LEARNING_ROOT / "cards" / f"{learning_id}.md").read_text(
                encoding="utf-8"
            )
            self.assertRegex(card_text, rf"(?m)^evidence_class: {evidence_class}$")

    def test_parallel_review_agents_are_not_global_learning_writers(self):
        evidence = self.agents["sci_evidence"][1]["developer_instructions"]
        router = self.agents["sci"][1]["developer_instructions"]
        self.assertIn("不得直接修改全局卡片或索引", evidence)
        self.assertIn("不得更新全局 `index.md` 或 `methods_index.md`", evidence)
        self.assertIn("全局学习库采用单写者", router)

    def test_learning_workspace_contract(self):
        required_files = [
            LEARNING_ROOT / "README.md",
            LEARNING_ROOT / "index.md",
            LEARNING_ROOT / "methods_index.md",
            LEARNING_ROOT / "figure_rule_ledger.tsv",
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

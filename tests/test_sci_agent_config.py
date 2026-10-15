from __future__ import annotations

import os
import pathlib
import re
import tomllib
import unittest


AGENT_DIR = pathlib.Path(os.environ.get("SCI_AGENT_DIR", pathlib.Path.home() / ".codex/agents"))
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
        if data["name"] in parsed:
            raise ValueError(f"Duplicate SCI role: {data['name']} ({path})")
        parsed[data["name"]] = (path, data)
    return parsed


class SciAgentConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.agents = configs()

    def require_terms(self, text, terms, owner):
        """Keep failures readable; this checks declarations, not agent behavior."""
        for term in terms:
            self.assertTrue(term in text, f"{owner}: missing declaration {term!r}")

    def contract(self, name):
        return (AGENT_DIR / name).read_text(encoding="utf-8")

    def test_shared_contract_references_resolve(self):
        sources = [path for path, _ in self.agents.values()]
        sources.extend(AGENT_DIR.glob("sci*contract.md"))
        for path in sources:
            refs = re.findall(r"\b(sci[a-z_]*_contract\.md)\b", path.read_text(encoding="utf-8"))
            for ref in refs:
                self.assertTrue((AGENT_DIR / ref).is_file(), f"{path.name}: missing {ref}")

    def test_all_roles_use_shared_scientific_core(self):
        core = self.contract("sci_core_contract.md")
        self.require_terms(core, ["contract_version:", "假说", "原始", "敏感", "授权"], "scientific core")
        for name, (_, data) in self.agents.items():
            self.assertIn("sci_core_contract.md", data["developer_instructions"], name)

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
        self.require_terms(instructions, ["并发", "当前会话"], "runtime delegation limits")

    def test_route_scenarios(self):
        router = self.agents["sci"][1]["developer_instructions"]
        cases = {
            "bulk/total/small RNA": ["sci_transcriptomics"],
            "scRNA": ["sci_singlecell"],
            "空间转录": ["sci_spatialomics"],
            "LC/GC-MS": ["sci_metabolomics"],
            "bulk ATAC": ["sci_epigenomics"],
            "跨组学整合": ["sci_multiomics"],
            "分类/回归": ["sci_ml"],
            "网络药理": ["sci_evidence", "sci_network_pharmacology"],
            "对接/MD": ["sci_molecular_modeling"],
        }
        for scenario, tokens in cases.items():
            rows = [line for line in router.splitlines()
                    if line.startswith("|") and scenario in line.split("|")[1]]
            self.assertEqual(len(rows), 1, f"Ambiguous or missing route: {scenario}")
            for token in tokens:
                self.assertIn(token, rows[0])

    def test_negative_scientific_boundaries_are_explicit(self):
        text = "\n".join(data["developer_instructions"] for _, data in self.agents.values())
        shared = (AGENT_DIR / "sci_analysis_contract.md").read_text(encoding="utf-8")
        combined = text + "\n" + shared + "\n" + self.contract("sci_core_contract.md")
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

        learning = (LEARNING_ROOT / "README.md").read_text(encoding="utf-8")
        card_template = (LEARNING_ROOT / "cards/_template.md").read_text(encoding="utf-8")
        self.require_terms(router, ["sci_learning/", "knowledge_refs"], "SCI retrieval entry")
        self.require_terms(shared, ["knowledge_refs", "knowledge_feedback"], "analysis handoff")
        self.require_terms(learning, ["adopt", "adapt", "background_only", "reject",
                                     "supported", "unsupported", "inconclusive", "not_run"],
                           "learning decisions and outcomes")

        for token in [
            "new_to_library",
            "known",
            "uncertain",
            "abstract_only",
            "搜索片段",
        ]:
            self.assertTrue(token in evidence + learning + card_template, f"learning authority: {token}")

    def test_quantitative_literature_learning_and_figure_value_gates(self):
        router = self.agents["sci"][1]["developer_instructions"]
        evidence = self.agents["sci_evidence"][1]["developer_instructions"]
        shared = (AGENT_DIR / "sci_analysis_contract.md").read_text(encoding="utf-8")

        for token in [
            "固定篇数",
            "精确发表日期",
            "规范化标题",
            "独立生物学单位",
            "引用数必须带数据库与快照日期",
        ]:
            self.assertTrue(token in evidence, f"evidence role: {token}")
        for token in [
            "方法采用以问题为单位",
            "estimand",
            "同一结果换画法不增加独立证据",
            "各门输入数、通过数、失败原因",
            "竞争型或自身型零假设",
            "不能用强制联合批次校正制造共同效应",
            "失败运行和负结果必须保留",
        ]:
            self.assertTrue(token in shared, f"analysis authority: {token}")
        self.require_terms(shared, ["NO-GO", "供者/患者/动物/窝", "跨队列", "物理尺度"],
                           "metadata eligibility and cohort analysis")

    def test_review_writing_and_visualization_learning_contracts(self):
        router = self.agents["sci"][1]["developer_instructions"]
        evidence = self.agents["sci_evidence"][1]["developer_instructions"]
        writer = self.agents["sci_writer"][1]["developer_instructions"]
        reviewer = self.agents["sci_reviewer"][1]["developer_instructions"]

        writing_contract = AGENT_DIR / "sci_review_writing_contract.md"
        figure_contract = AGENT_DIR / "sci_review_figure_contract.md"
        self.assertTrue(writing_contract.is_file())
        self.assertTrue(figure_contract.is_file())

        self.require_terms(router, ["sci_writer", "sci_reviewer"], "writing/review delegation")
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
            "程序化矢量",
            "target-review deconstruction",
            "orienting_overview",
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
            self.assertTrue(token in figure_contract.read_text(encoding="utf-8"),
                            f"review semantics: {token}")

        self.require_terms(self.contract("sci_figure_contract.md"),
                           ["MediaBox/CropBox", "字体嵌入", "deuteranopia", "源数据", "画布"],
                           "shared figure production")
        self.assertIn("sci_figure_contract.md", figure_contract.read_text(encoding="utf-8"))
        self.assertNotIn("原生矢量面板按1.0比例拼版", figure_contract.read_text(encoding="utf-8"))

        self.require_terms(evidence, ["主SCI", "不得并发修改全局"], "evidence writing ownership")

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
        self.require_terms(evidence, ["不得并发修改全局", "主SCI", "串行"],
                           "evidence writing ownership")
        learning = (LEARNING_ROOT / "README.md").read_text(encoding="utf-8")
        self.require_terms(router + learning, ["单写者", "串行"], "global learning ownership")

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
        combined = "\n".join([router, evidence, shared,
                               (LEARNING_ROOT / "README.md").read_text(encoding="utf-8"),
                               (LEARNING_ROOT / "repositories/_template.md").read_text(encoding="utf-8")])
        for token in [
            "github-",
            "gitlab-",
            "checked_commit",
            "metadata_only",
            "readme_docs",
            "source_audit",
            "executed",
            "许可证",
        ]:
            self.assertTrue(token in combined, f"repository learning authority: {token}")

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

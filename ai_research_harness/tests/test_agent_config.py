from __future__ import annotations

import tomllib
from collections import Counter
from pathlib import Path

from ai_research_harness.gate import load_qualification_tasks
from ai_research_harness.runtime import ASSETS


AGENT_DIR = Path("/home/h1028/.codex/agents")


def test_agent_family_is_complete_and_unique() -> None:
    parsed = [tomllib.loads(path.read_text(encoding="utf-8")) for path in sorted(AGENT_DIR.glob("ai_research*.toml"))]
    expected = {
        "ai_research",
        "ai_research_evidence",
        "ai_research_methodologist",
        "ai_research_ml",
        "ai_research_engineer",
        "ai_research_systems",
        "ai_research_theory",
        "ai_research_security",
        "ai_research_hci",
        "ai_research_writer",
        "ai_research_reviewer",
        "ai_research_memory_evolver",
    }
    assert {item["name"] for item in parsed} == expected
    assert len(parsed) == 12
    assert all(item["model"] == "gpt-5.6-sol" for item in parsed)
    assert all("ai_research_contract.md" in item["developer_instructions"] for item in parsed)


def test_router_names_roles_and_bounds_delegation() -> None:
    router = tomllib.loads((AGENT_DIR / "ai_research.toml").read_text(encoding="utf-8"))["developer_instructions"]
    for role in [
        "ai_research_evidence",
        "ai_research_methodologist",
        "ai_research_ml",
        "ai_research_engineer",
        "ai_research_systems",
        "ai_research_theory",
        "ai_research_security",
        "ai_research_hci",
        "ai_research_writer",
        "ai_research_reviewer",
        "ai_research_memory_evolver",
    ]:
        assert role in router
    assert "最多同时运行三个" in router
    assert "只有 `ai_research` 可以创建子 agent" in router


def test_sandbox_boundaries() -> None:
    configs = {
        data["name"]: data
        for data in [tomllib.loads(path.read_text(encoding="utf-8")) for path in AGENT_DIR.glob("ai_research*.toml")]
    }
    for name in {
        "ai_research_evidence",
        "ai_research_methodologist",
        "ai_research_theory",
        "ai_research_security",
        "ai_research_hci",
        "ai_research_reviewer",
        "ai_research_memory_evolver",
    }:
        assert configs[name]["sandbox_mode"] == "read-only"
    assert "不得委派" in configs["ai_research_reviewer"]["developer_instructions"]
    assert "不得直接修改 champion" in configs["ai_research_memory_evolver"]["developer_instructions"]


def test_shared_contract_has_domain_and_claim_gates() -> None:
    contract = (AGENT_DIR / "ai_research_contract.md").read_text(encoding="utf-8")
    for token in [
        "independent_reproduction",
        "benchmark contamination",
        "额外重试",
        "p50/p95/p99",
        "证明义务",
        "威胁模型",
        "参与者",
        "E（Experiential memory）",
        "rho（Invocation policy）",
        "item/completed",
    ]:
        assert token in contract


def test_frozen_qualification_distribution() -> None:
    tasks = load_qualification_tasks(ASSETS / "evals" / "frozen_qualification.yaml")
    assert len(tasks) == 32
    assert Counter(task.domain for task in tasks) == {
        "evidence": 8,
        "empirical_ml": 8,
        "reproducibility": 4,
        "systems": 3,
        "theory": 3,
        "security": 3,
        "hci": 3,
    }
    assert len({task.task_id for task in tasks}) == 32
    assert all(task.gold_fragments for task in tasks)

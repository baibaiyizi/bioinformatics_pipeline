from __future__ import annotations

from pathlib import Path

import pytest

from ai_research_harness.memory import MemorySelector, find_gold_leakage, find_near_duplicates, load_bundle_cards
from ai_research_harness.models import Goal, SkillCard, WorkingMemory
from ai_research_harness.runtime import ASSETS


def state(tmp_path: Path, objective: str = "比较 LLM agent retry budget") -> WorkingMemory:
    return WorkingMemory(
        run_id="r",
        task_id="t",
        project_root=str(tmp_path),
        goals=[Goal(goal_id="g", objective=objective, acceptance_checks=["a"])],
    )


def test_recuris_bundle_loads_and_selector_respects_limits(tmp_path: Path) -> None:
    _, cards = load_bundle_cards(ASSETS / "memory" / "recuris-v1")
    selected = MemorySelector(max_cards=2, max_tokens=3500).select(state(tmp_path), cards)
    assert 1 <= len(selected) <= 2
    assert any(card.skill_id == "matched-agent-evaluation-budget" for card in selected)


def test_project_conflict_blocks(tmp_path: Path) -> None:
    _, cards = load_bundle_cards(ASSETS / "memory" / "recuris-v1")
    conflict = cards[0].model_copy(update={"scope": "project", "body": "different"})
    with pytest.raises(ValueError, match="conflicting"):
        MemorySelector().select(state(tmp_path, "paper evidence"), cards, [conflict])


def test_relevant_project_card_precedes_global_card(tmp_path: Path) -> None:
    _, cards = load_bundle_cards(ASSETS / "memory" / "recuris-v1")
    project = cards[0].model_copy(
        update={
            "skill_id": "project-agent-rule",
            "scope": "project",
            "domains": ["agent"],
            "triggers": ["agent"],
            "path": tmp_path / "project.md",
        }
    )
    selected = MemorySelector(max_cards=1).select(state(tmp_path), cards, [project])
    assert [card.skill_id for card in selected] == ["project-agent-rule"]


def test_duplicate_and_gold_detection(tmp_path: Path) -> None:
    card = SkillCard(
        skill_id="a",
        version=1,
        scope="global",
        domains=["agent"],
        triggers=["retry"],
        sources=["s"],
        evidence_requirements=["e"],
        path=tmp_path / "a.md",
        body="matched retry budget and exact frozen-secret-fragment",
        token_estimate=8,
    )
    clone = card.model_copy(update={"skill_id": "b", "path": tmp_path / "b.md"})
    assert find_near_duplicates([card, clone])
    assert find_gold_leakage([card], ["frozen-secret-fragment"])

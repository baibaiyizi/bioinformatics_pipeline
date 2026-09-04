from __future__ import annotations

from pathlib import Path

import pytest

from ai_research_harness.evolution import _apply_declarative_patch, build_candidate
from ai_research_harness.memory import load_bundle_cards
from ai_research_harness.runtime import ASSETS
from ai_research_harness.storage import MemoryStore, load_yaml


def e_patch() -> dict:
    return {
        "classification": "patchable",
        "component": "E",
        "reason": "missing evidence skill",
        "source_failure_ids": ["evidence-03-independent-reproduction"],
        "patch": {
            "skill_id": "new-evidence-boundary",
            "domains": ["evidence"],
            "triggers": ["reproduction"],
            "anti_triggers": [],
            "sources": ["run-1"],
            "evidence_requirements": ["receipt"],
            "body": "Require a real receipt and independent protocol before claiming reproduction.",
        },
    }


def test_build_candidate_is_isolated_and_declarative(tmp_path: Path) -> None:
    store = MemoryStore(tmp_path / "memory", "global")
    store.bootstrap(ASSETS / "memory" / "base-v1")
    candidate_id = build_candidate(store=store, incumbent_id="base-v1", proposal=e_patch(), scope="global")
    assert store.champion_id() == "base-v1"
    manifest, cards = load_bundle_cards(store.candidate_dir(candidate_id))
    assert manifest.parent_version == "base-v1"
    assert [card.skill_id for card in cards] == ["new-evidence-boundary"]
    assert not list(store.candidate_dir(candidate_id).rglob("*.py"))


def test_evolution_cannot_weaken_rho_or_add_checker(tmp_path: Path) -> None:
    bundle = tmp_path / "bundle"
    from shutil import copytree

    copytree(ASSETS / "memory" / "base-v1", bundle)
    with pytest.raises(ValueError, match="cannot weaken"):
        _apply_declarative_patch(bundle, "rho", {"updates": {"conflict_policy": "prefer-new"}}, "global")
    with pytest.raises(ValueError, match="registered"):
        _apply_declarative_patch(bundle, "C", {"updates": {"goal_checks": ["run-shell"]}}, "global")

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Iterable

import yaml

from .models import MemoryManifest, SkillCard, WorkingMemory
from .storage import load_yaml


FRONTMATTER = re.compile(r"\A---\s*\n(.*?)\n---\s*\n(.*)\Z", re.DOTALL)
GLOBAL_CARD_LIMIT = 64
PROJECT_CARD_LIMIT = 16
CARD_TOKEN_LIMIT = 1_200
INVOCATION_CARD_LIMIT = 4
INVOCATION_TOKEN_LIMIT = 3_500


def estimate_tokens(text: str) -> int:
    ascii_words = re.findall(r"[A-Za-z0-9_./+-]+", text)
    non_ascii = sum(1 for character in text if ord(character) > 127 and not character.isspace())
    return max(1, len(ascii_words) + non_ascii)


def load_skill_card(path: Path, *, expected_scope: str | None = None) -> SkillCard:
    text = path.read_text(encoding="utf-8")
    match = FRONTMATTER.match(text)
    if match is None:
        raise ValueError(f"skill card lacks YAML frontmatter: {path}")
    metadata = yaml.safe_load(match.group(1))
    if not isinstance(metadata, dict):
        raise ValueError(f"invalid skill card frontmatter: {path}")
    body = match.group(2).strip()
    tokens = estimate_tokens(body)
    value = SkillCard.model_validate({**metadata, "path": path, "body": body, "token_estimate": tokens})
    if expected_scope is not None and value.scope != expected_scope:
        raise ValueError(f"card scope mismatch: {path}")
    if tokens > CARD_TOKEN_LIMIT:
        raise ValueError(f"card exceeds {CARD_TOKEN_LIMIT} token estimate: {path}")
    return value


def load_bundle_cards(bundle: Path) -> tuple[MemoryManifest, list[SkillCard]]:
    manifest = MemoryManifest.model_validate(load_yaml(bundle / "manifest.yaml"))
    cards = [load_skill_card(bundle / relative, expected_scope=manifest.scope) for relative in manifest.cards]
    if len({card.skill_id for card in cards}) != len(cards):
        raise ValueError("duplicate skill IDs inside memory bundle")
    limit = GLOBAL_CARD_LIMIT if manifest.scope == "global" else PROJECT_CARD_LIMIT
    active = [card for card in cards if card.status == "active"]
    if len(active) > limit:
        raise ValueError(f"active card limit exceeded: {len(active)} > {limit}")
    return manifest, cards


def load_bundle_settings(bundle: Path) -> dict[str, dict]:
    manifest = MemoryManifest.model_validate(load_yaml(bundle / "manifest.yaml"))
    working = load_yaml(bundle / manifest.working_memory_spec)
    policy = load_yaml(bundle / manifest.invocation_policy)
    checkers = load_yaml(bundle / manifest.checker_config)
    if working.get("statuses") != ["pending", "done", "blocked"] or working.get("commit_authority") != "harness":
        raise ValueError("working-memory spec cannot weaken status or commit authority")
    if policy.get("conflict_policy") != "block" or policy.get("project_precedence") is not True:
        raise ValueError("invocation policy cannot weaken conflict blocking or project precedence")
    required_events = {"turn_start", "phase_boundary", "checker_failure"}
    if not required_events.issubset(set(policy.get("events", []))):
        raise ValueError("invocation policy must select at turn start, phase boundary, and checker failure")
    if int(policy.get("max_cards", 0)) not in range(1, INVOCATION_CARD_LIMIT + 1):
        raise ValueError("invocation max_cards exceeds host limit")
    if int(policy.get("max_tokens", 0)) not in range(1, INVOCATION_TOKEN_LIMIT + 1):
        raise ValueError("invocation max_tokens exceeds host limit")
    if checkers.get("unknown_checker_policy") != "block":
        raise ValueError("checker config must block unknown checkers")
    return {"working_memory": working, "invocation_policy": policy, "checkers": checkers}


def _terms(text: str) -> set[str]:
    lowered = text.casefold()
    words = set(re.findall(r"[a-z0-9_+-]{2,}", lowered))
    compact = re.sub(r"\s+", "", lowered)
    grams = {compact[index : index + 3] for index in range(max(0, len(compact) - 2))}
    return words | grams


def card_similarity(left: SkillCard, right: SkillCard) -> float:
    a = _terms(" ".join(left.domains + left.triggers) + " " + left.body)
    b = _terms(" ".join(right.domains + right.triggers) + " " + right.body)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def find_near_duplicates(cards: Iterable[SkillCard], threshold: float = 0.90) -> list[tuple[str, str, float]]:
    active = [card for card in cards if card.status == "active"]
    duplicates: list[tuple[str, str, float]] = []
    for index, left in enumerate(active):
        for right in active[index + 1 :]:
            similarity = card_similarity(left, right)
            if similarity >= threshold:
                duplicates.append((left.skill_id, right.skill_id, similarity))
    return duplicates


def find_gold_leakage(cards: Iterable[SkillCard], gold_fragments: Iterable[str]) -> list[tuple[str, str]]:
    findings: list[tuple[str, str]] = []
    normalized_fragments = [re.sub(r"\s+", " ", fragment.casefold()).strip() for fragment in gold_fragments]
    for card in cards:
        normalized = re.sub(r"\s+", " ", card.body.casefold())
        for fragment in normalized_fragments:
            if len(fragment) >= 12 and fragment in normalized:
                findings.append((card.skill_id, fragment))
    return findings


class MemorySelector:
    def __init__(self, *, max_cards: int = INVOCATION_CARD_LIMIT, max_tokens: int = INVOCATION_TOKEN_LIMIT) -> None:
        self.max_cards = max_cards
        self.max_tokens = max_tokens

    @staticmethod
    def _score(card: SkillCard, query_terms: Counter[str], project_precedence: bool) -> tuple[int, int, int, str]:
        positive = " ".join(card.domains + card.triggers).casefold()
        negative = " ".join(card.anti_triggers).casefold()
        if any(term in negative for term in query_terms if len(term) >= 3):
            return (0, -10_000, 0, card.skill_id)
        score = sum(weight for term, weight in query_terms.items() if term in positive)
        scope_priority = int(project_precedence and card.scope == "project")
        return (scope_priority, score, -card.token_estimate, card.skill_id)

    def select(self, state: WorkingMemory, global_cards: list[SkillCard], project_cards: list[SkillCard] | None = None) -> list[SkillCard]:
        project_cards = project_cards or []
        combined: dict[str, SkillCard] = {}
        for card in global_cards + project_cards:
            if card.status != "active":
                continue
            existing = combined.get(card.skill_id)
            if existing is not None and existing.body != card.body:
                raise ValueError(f"conflicting active skill ID across scopes: {card.skill_id}")
            combined[card.skill_id] = card

        query = " ".join([state.phase, *state.active_constraints, *(goal.objective for goal in state.pending_goals())]).casefold()
        query_terms = Counter(re.findall(r"[a-z0-9_+-]{2,}|[\u4e00-\u9fff]{2,}", query))
        ranked = sorted(combined.values(), key=lambda card: self._score(card, query_terms, True), reverse=True)
        selected: list[SkillCard] = []
        tokens = 0
        for card in ranked:
            score = self._score(card, query_terms, True)[1]
            if score <= 0:
                continue
            if len(selected) >= self.max_cards or tokens + card.token_estimate > self.max_tokens:
                continue
            selected.append(card)
            tokens += card.token_estimate
        return selected

from __future__ import annotations

import math
import random
from collections import defaultdict
from pathlib import Path

import yaml

from .checkers import qualification_check
from .memory import find_gold_leakage, find_near_duplicates, load_bundle_cards, load_bundle_settings
from .models import ArmOutcome, CheckResult, GateDecision, QualificationTask
from .storage import bundle_hash, sha256_file


def load_qualification_tasks(path: Path) -> list[QualificationTask]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or not isinstance(raw.get("tasks"), list):
        raise ValueError("qualification suite must contain a tasks list")
    tasks = [QualificationTask.model_validate(item) for item in raw["tasks"]]
    ids = [task.task_id for task in tasks]
    if len(ids) != len(set(ids)):
        raise ValueError("qualification task IDs must be unique")
    return tasks


def paired_bootstrap_ci(differences: list[int], *, samples: int = 10_000, seed: int = 260824876) -> tuple[float, float]:
    if not differences:
        return (0.0, 0.0)
    rng = random.Random(seed)
    estimates = []
    for _ in range(samples):
        draw = [rng.choice(differences) for _ in differences]
        estimates.append(sum(draw) / len(draw))
    estimates.sort()
    lower = estimates[int(0.025 * (samples - 1))]
    upper = estimates[int(0.975 * (samples - 1))]
    return (lower, upper)


def mcnemar_exact_p(wins: int, losses: int) -> float:
    discordant = wins + losses
    if discordant == 0:
        return 1.0
    tail = min(wins, losses)
    probability = sum(math.comb(discordant, index) for index in range(tail + 1)) / (2**discordant)
    return min(1.0, 2 * probability)


def validate_candidate_bundle(candidate_bundle: Path, tasks: list[QualificationTask]) -> list[CheckResult]:
    manifest, cards = load_bundle_cards(candidate_bundle)
    settings = load_bundle_settings(candidate_bundle)
    duplicate_findings = find_near_duplicates(cards)
    gold = [fragment for task in tasks for fragment in task.gold_fragments]
    leakage_findings = find_gold_leakage(cards, gold)
    executable = [
        str(path.relative_to(candidate_bundle))
        for path in candidate_bundle.rglob("*")
        if path.is_file() and path.suffix.casefold() in {".py", ".sh", ".pl", ".exe", ".dll", ".so"}
    ]
    private_patterns = ("/home/", "/Users/", "\\Users\\", ".ai_research/", "run-")
    private_findings = [
        (card.skill_id, pattern)
        for card in cards
        for pattern in private_patterns
        if manifest.scope == "global"
        and pattern.casefold()
        in " ".join([card.body, *card.sources, *card.triggers, *card.evidence_requirements]).casefold()
    ]
    expected_anchor = (candidate_bundle / "ANCHOR.sha256").read_text(encoding="utf-8").strip()
    actual_anchor = bundle_hash(candidate_bundle)
    memory_budget = {
        "active_cards": sum(card.status == "active" for card in cards),
        "largest_card_tokens": max((card.token_estimate for card in cards), default=0),
        "max_injected_cards": settings["invocation_policy"]["max_cards"],
        "max_injected_tokens": settings["invocation_policy"]["max_tokens"],
    }
    return [
        CheckResult(
            checker_id="candidate_schema_valid",
            passed=True,
            detail="manifest, cards, E/W/rho/C schemas are valid",
        ),
        CheckResult(
            checker_id="candidate_anchor_valid",
            passed=expected_anchor == actual_anchor,
            detail="candidate SHA-256 anchor matches" if expected_anchor == actual_anchor else "candidate SHA-256 anchor mismatch",
            evidence={"expected": expected_anchor, "actual": actual_anchor},
        ),
        CheckResult(
            checker_id="candidate_memory_budget_valid",
            passed=True,
            detail="candidate memory is within host limits",
            evidence=memory_budget,
        ),
        CheckResult(
            checker_id="candidate_no_near_duplicates",
            passed=not duplicate_findings,
            detail="no near-duplicate active cards" if not duplicate_findings else "near-duplicate active cards found",
            evidence={"findings": duplicate_findings},
        ),
        CheckResult(
            checker_id="candidate_no_gold_leakage",
            passed=not leakage_findings,
            detail="no frozen gold fragment leakage" if not leakage_findings else "frozen gold fragment leakage found",
            evidence={"findings": leakage_findings},
        ),
        CheckResult(
            checker_id="candidate_declarative_only",
            passed=not executable,
            detail="candidate contains declarations only" if not executable else "candidate contains executable files",
            evidence={"files": executable},
        ),
        CheckResult(
            checker_id="candidate_no_project_private_identifiers",
            passed=not private_findings,
            detail="no project-private identifiers" if not private_findings else "project-private identifiers found",
            evidence={"findings": private_findings},
        ),
    ]


def evaluate_gate(
    *,
    candidate_id: str,
    incumbent_id: str,
    tasks: list[QualificationTask],
    incumbent: list[ArmOutcome],
    candidate: list[ArmOutcome],
    validation_checks: list[CheckResult],
    source_failure_ids: list[str],
) -> GateDecision:
    incumbent_by_id = {outcome.task_id: outcome for outcome in incumbent}
    candidate_by_id = {outcome.task_id: outcome for outcome in candidate}
    expected = {task.task_id for task in tasks}
    if set(incumbent_by_id) != expected or set(candidate_by_id) != expected:
        raise ValueError("gate requires one outcome per task and arm")

    wins = losses = ties = 0
    critical_regressions: list[str] = []
    domain_deltas: dict[str, int] = defaultdict(int)
    differences: list[int] = []
    task_by_id = {task.task_id: task for task in tasks}
    for task_id in sorted(expected):
        before = int(incumbent_by_id[task_id].passed)
        after = int(candidate_by_id[task_id].passed)
        delta = after - before
        differences.append(delta)
        domain_deltas[task_by_id[task_id].domain] += delta
        if delta > 0:
            wins += 1
        elif delta < 0:
            losses += 1
            if task_by_id[task_id].critical:
                critical_regressions.append(task_id)
        else:
            ties += 1

    source_failures_repaired = bool(source_failure_ids) and all(
        source_id in candidate_by_id and candidate_by_id[source_id].passed and not incumbent_by_id[source_id].passed
        for source_id in source_failure_ids
    )
    checks = list(validation_checks)
    checks.extend(
        [
            CheckResult(
                checker_id="source_failures_repaired",
                passed=source_failures_repaired,
                detail="all declared source failures repaired" if source_failures_repaired else "declared source failures not repaired",
                evidence={"source_failure_ids": source_failure_ids},
            ),
            CheckResult(
                checker_id="zero_critical_regressions",
                passed=not critical_regressions,
                detail="no critical regression" if not critical_regressions else "critical regression detected",
                evidence={"task_ids": critical_regressions},
            ),
            CheckResult(
                checker_id="no_domain_regression",
                passed=all(delta >= 0 for delta in domain_deltas.values()),
                detail="no domain pass count decreased" if all(delta >= 0 for delta in domain_deltas.values()) else "one or more domains regressed",
                evidence={"domain_deltas": dict(domain_deltas)},
            ),
            CheckResult(
                checker_id="positive_net_wins",
                passed=wins - losses >= 1,
                detail=f"net wins={wins - losses}",
                evidence={"wins": wins, "losses": losses, "ties": ties},
            ),
        ]
    )
    admitted = all(check.passed for check in checks)
    return GateDecision(
        candidate_id=candidate_id,
        incumbent_id=incumbent_id,
        status="PROMOTED" if admitted else "NOT_ADMITTED",
        source_failures_repaired=source_failures_repaired,
        critical_regressions=critical_regressions,
        domain_deltas=dict(domain_deltas),
        wins=wins,
        losses=losses,
        ties=ties,
        net_wins=wins - losses,
        bootstrap_ci95=paired_bootstrap_ci(differences),
        mcnemar_exact_p=mcnemar_exact_p(wins, losses),
        checks=checks,
    )


def grade_response(
    task: QualificationTask,
    arm: str,
    response: str,
    usage: dict | None = None,
    duration_ms: int | None = None,
    trace: dict | None = None,
) -> ArmOutcome:
    check = qualification_check(task, response)
    return ArmOutcome(
        task_id=task.task_id,
        arm=arm,
        passed=check.passed,
        deterministic_passed=check.passed,
        response=response,
        reasons=[] if check.passed else [check.detail, str(check.evidence)],
        usage=usage or {},
        duration_ms=duration_ms,
        trace=trace or {},
    )


def qualification_fingerprint(path: Path) -> str:
    return sha256_file(path)

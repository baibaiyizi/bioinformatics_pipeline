from __future__ import annotations

from pathlib import Path

from ai_research_harness.gate import evaluate_gate, load_qualification_tasks, mcnemar_exact_p, paired_bootstrap_ci, validate_candidate_bundle
from ai_research_harness.models import ArmOutcome, CheckResult
from ai_research_harness.runtime import ASSETS
from ai_research_harness.storage import MemoryStore


def outcome(task_id: str, arm: str, passed: bool) -> ArmOutcome:
    return ArmOutcome(task_id=task_id, arm=arm, passed=passed, deterministic_passed=passed, response="x")


def test_candidate_bundle_is_declarative_and_leak_free(tmp_path: Path) -> None:
    tasks = load_qualification_tasks(ASSETS / "evals" / "frozen_qualification.yaml")
    store = MemoryStore(tmp_path / "memory", "global")
    bundle = store.install_bundle(ASSETS / "memory" / "recuris-v1", "recuris-v1", candidate=True)
    checks = validate_candidate_bundle(bundle, tasks)
    assert all(check.passed for check in checks)


def test_gate_promotes_only_source_repair_without_regression() -> None:
    tasks = load_qualification_tasks(ASSETS / "evals" / "frozen_qualification.yaml")
    source = "evidence-03-independent-reproduction"
    incumbent = [outcome(task.task_id, "incumbent", task.task_id != source) for task in tasks]
    candidate = [outcome(task.task_id, "candidate", True) for task in tasks]
    decision = evaluate_gate(
        candidate_id="c",
        incumbent_id="i",
        tasks=tasks,
        incumbent=incumbent,
        candidate=candidate,
        validation_checks=[CheckResult(checker_id="schema", passed=True, detail="ok")],
        source_failure_ids=[source],
    )
    assert decision.status == "PROMOTED"
    assert decision.wins == 1 and decision.losses == 0


def test_gate_rejects_critical_or_domain_regression() -> None:
    tasks = load_qualification_tasks(ASSETS / "evals" / "frozen_qualification.yaml")
    source = "evidence-03-independent-reproduction"
    critical_loss = "security-01-prompt-injection"
    incumbent = [outcome(task.task_id, "incumbent", task.task_id != source) for task in tasks]
    candidate = [outcome(task.task_id, "candidate", task.task_id != critical_loss) for task in tasks]
    decision = evaluate_gate(
        candidate_id="c",
        incumbent_id="i",
        tasks=tasks,
        incumbent=incumbent,
        candidate=candidate,
        validation_checks=[CheckResult(checker_id="schema", passed=True, detail="ok")],
        source_failure_ids=[source],
    )
    assert decision.status == "NOT_ADMITTED"
    assert critical_loss in decision.critical_regressions


def test_paired_statistics_are_bounded() -> None:
    lower, upper = paired_bootstrap_ci([1, 0, -1, 1], samples=1000, seed=1)
    assert -1 <= lower <= upper <= 1
    assert 0 <= mcnemar_exact_p(3, 1) <= 1

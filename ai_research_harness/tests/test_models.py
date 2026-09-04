from __future__ import annotations

from pathlib import Path

import pytest

from ai_research_harness.models import CheckResult, Goal, GoalStatus, Receipt, WorkingMemory


def completed_receipt(authoritative: bool = True) -> Receipt:
    return Receipt(
        receipt_id="r1",
        kind="turn",
        authoritative=authoritative,
        status="completed",
        item_id="turn-1",
    )


def test_goal_requires_authoritative_receipt_and_all_checks() -> None:
    goal = Goal(goal_id="g", objective="deliver", acceptance_checks=["a"])
    with pytest.raises(ValueError, match="authoritative"):
        goal.propose_done([completed_receipt(False)], [CheckResult(checker_id="a", passed=True, detail="ok")])
    with pytest.raises(ValueError, match="failed"):
        goal.propose_done([completed_receipt()], [CheckResult(checker_id="a", passed=False, detail="bad")])
    goal.propose_done([completed_receipt()], [CheckResult(checker_id="a", passed=True, detail="ok")])
    assert goal.status == GoalStatus.done


def test_artifact_receipt_requires_digest() -> None:
    with pytest.raises(ValueError, match="sha256"):
        Receipt(receipt_id="a", kind="artifact", authoritative=True, status="completed", artifact_path="/tmp/a")


def test_working_memory_rejects_relative_project_and_duplicate_goals(tmp_path: Path) -> None:
    goal = Goal(goal_id="g", objective="x", acceptance_checks=["a"])
    with pytest.raises(ValueError, match="absolute"):
        WorkingMemory(run_id="r", task_id="t", project_root="relative", goals=[goal])
    with pytest.raises(ValueError, match="unique"):
        WorkingMemory(run_id="r", task_id="t", project_root=str(tmp_path), goals=[goal, goal.model_copy()])


def test_only_pending_goal_can_block() -> None:
    goal = Goal(goal_id="g", objective="x", acceptance_checks=["a"])
    goal.block("missing input")
    assert goal.status == GoalStatus.blocked
    with pytest.raises(ValueError, match="pending"):
        goal.block("again")


def test_zero_exit_code_cannot_override_failed_artifact_checker() -> None:
    goal = Goal(goal_id="g", objective="produce valid artifact", acceptance_checks=["artifact_sha256"])
    command = Receipt(
        receipt_id="command",
        kind="command_execution",
        authoritative=True,
        status="completed",
        item_id="command-1",
        details={"exit_code": 0},
    )
    failed_check = CheckResult(checker_id="artifact_sha256", passed=False, detail="digest mismatch")
    with pytest.raises(ValueError, match="failed"):
        goal.propose_done([command], [failed_check])

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

from .models import CheckResult, QualificationTask, Receipt
from .storage import sha256_file


Checker = Callable[..., CheckResult]


class CheckerRegistry:
    def __init__(self) -> None:
        self._checkers: dict[str, Checker] = {}

    def register(self, checker_id: str, checker: Checker) -> None:
        if checker_id in self._checkers:
            raise ValueError(f"checker already registered: {checker_id}")
        self._checkers[checker_id] = checker

    def run(self, checker_id: str, **kwargs: Any) -> CheckResult:
        if checker_id not in self._checkers:
            raise KeyError(f"unknown checker: {checker_id}")
        return self._checkers[checker_id](**kwargs)

    @property
    def ids(self) -> set[str]:
        return set(self._checkers)


def check_agent_response_present(*, response: str | None, **_: Any) -> CheckResult:
    passed = bool(response and response.strip())
    return CheckResult(
        checker_id="agent_response_present",
        passed=passed,
        detail="final agent response is present" if passed else "final agent response is missing",
    )


def check_no_blocked_marker(*, response: str | None, **_: Any) -> CheckResult:
    text = (response or "").casefold()
    blocked_markers = ("status: blocked", '"status":"blocked"', "not_admitted")
    matches = [marker for marker in blocked_markers if marker in text]
    return CheckResult(
        checker_id="no_blocked_marker",
        passed=not matches,
        detail="no blocked marker" if not matches else f"blocked marker found: {matches}",
    )


def check_artifact(*, path: Path, expected_sha256: str | None = None, **_: Any) -> CheckResult:
    exists = path.is_file()
    actual = sha256_file(path) if exists else None
    passed = exists and (expected_sha256 is None or actual == expected_sha256)
    return CheckResult(
        checker_id="artifact_sha256",
        passed=passed,
        detail="artifact exists and digest matches" if passed else "artifact missing or digest mismatch",
        evidence={"path": str(path), "sha256": actual},
    )


def qualification_check(task: QualificationTask, response: str) -> CheckResult:
    text = response.casefold()
    missing_groups = [group for group in task.required_groups if not any(term.casefold() in text for term in group)]
    forbidden_hits = [term for term in task.forbidden if term.casefold() in text]
    passed = not missing_groups and not forbidden_hits
    return CheckResult(
        checker_id=f"qualification:{task.task_id}",
        passed=passed,
        detail="qualification rubric passed" if passed else "qualification rubric failed",
        evidence={"missing_groups": missing_groups, "forbidden_hits": forbidden_hits},
    )


def authoritative_receipts(receipts: list[Receipt]) -> list[Receipt]:
    return [receipt for receipt in receipts if receipt.authoritative and receipt.status == "completed"]


def default_registry() -> CheckerRegistry:
    registry = CheckerRegistry()
    registry.register("agent_response_present", check_agent_response_present)
    registry.register("no_blocked_marker", check_no_blocked_marker)
    registry.register("artifact_sha256", check_artifact)
    return registry

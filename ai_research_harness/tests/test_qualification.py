from __future__ import annotations

import asyncio
from pathlib import Path

from ai_research_harness.gate import load_qualification_tasks
from ai_research_harness.qualification import QualificationRunner, _budget_violations
from ai_research_harness.runtime import ASSETS
from ai_research_harness.sdk_adapter import SDKRunResult
from tests.fakes import FakeBackend


class ConcurrencyBackend(FakeBackend):
    def __init__(self) -> None:
        super().__init__()
        self.active = 0
        self.max_active = 0

    async def execute(self, **kwargs):
        self.active += 1
        self.max_active = max(self.max_active, self.active)
        try:
            await asyncio.sleep(0.01)
            return await super().execute(**kwargs)
        finally:
            self.active -= 1


def test_qualification_limits_total_sdk_turn_concurrency(tmp_path: Path) -> None:
    tasks = load_qualification_tasks(ASSETS / "evals" / "frozen_qualification.yaml")
    selected = [tasks[0], tasks[2]]
    backend = ConcurrencyBackend()
    runner = QualificationRunner(backend, concurrency=1)
    incumbent, candidate, trials = asyncio.run(
        runner.run(
            tasks=selected,
            incumbent_bundle=ASSETS / "memory" / "base-v1",
            candidate_bundle=ASSETS / "memory" / "recuris-v1",
            cwd=tmp_path,
        )
    )
    assert backend.max_active == 1
    assert len(incumbent) == len(candidate) == len(trials) == 2
    assert all(call["sandbox"] == "read_only" for call in backend.calls)


def test_qualification_rejects_tool_and_token_budget_overrun() -> None:
    result = SDKRunResult(
        thread_id="thread",
        turn_id="turn",
        status="completed",
        final_response="x",
        receipts=[],
        usage={"total": {"totalTokens": 100_001}},
        completed_items=[{"type": "commandExecution", "status": "completed", "exitCode": 0}],
    )
    violations = _budget_violations(result)
    assert any("tool action" in violation for violation in violations)
    assert any("token budget" in violation for violation in violations)

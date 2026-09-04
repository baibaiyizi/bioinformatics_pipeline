from __future__ import annotations

import asyncio
import hashlib
from pathlib import Path

import pytest

from ai_research_harness.models import GoalStatus, Receipt
from ai_research_harness.runtime import HarnessRuntime, initialize_project
from ai_research_harness.storage import RunStore
from tests.fakes import FakeBackend


class InterruptingBackend(FakeBackend):
    async def execute(self, **kwargs):
        callback = kwargs["trace_callback"]
        result = callback(
            "contextCompaction/completed",
            {"thread_id": "thread-recover", "turn_id": "turn-interrupted", "method": "contextCompaction/completed"},
        )
        if result is not None:
            await result
        raise asyncio.CancelledError()


class ReceiptBackend(FakeBackend):
    def __init__(self, receipt: Receipt | None) -> None:
        super().__init__(response="创建文件完成。")
        self.receipt = receipt

    async def execute(self, **kwargs):
        result = await super().execute(**kwargs)
        if self.receipt is not None:
            result.receipts.append(self.receipt)
        return result


def test_runtime_commits_done_from_receipt_and_checks(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AI_RESEARCH_MEMORY_ROOT", str(tmp_path / "global-memory"))
    project = tmp_path / "project"
    project.mkdir()
    backend = FakeBackend("已完成，并明确 evidence boundary。")
    state, result = asyncio.run(HarnessRuntime(backend).run(project, "评估 agent retry budget"))
    assert state.goals[0].status == GoalStatus.done
    assert result.thread_id == state.thread_id
    run_dir = RunStore(project).run_dir(state.run_id)
    assert (run_dir / "response.md").is_file()
    assert "state_commit" in (run_dir / "trace.jsonl").read_text(encoding="utf-8")
    assert backend.calls[0]["sandbox"] == "workspace_write"
    assert (RunStore(project).anchors / "qualification.yaml").is_file()
    assert (RunStore(project).anchors / "qualification.sha256").is_file()


def test_runtime_blocks_failed_turn(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AI_RESEARCH_MEMORY_ROOT", str(tmp_path / "global-memory"))
    project = tmp_path / "project"
    project.mkdir()
    backend = FakeBackend("failure", status="failed")
    state, _ = asyncio.run(HarnessRuntime(backend).run(project, "run benchmark"))
    assert state.goals[0].status == GoalStatus.blocked


def test_resume_reuses_thread_and_refuses_done_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AI_RESEARCH_MEMORY_ROOT", str(tmp_path / "global-memory"))
    project = tmp_path / "project"
    project.mkdir()
    backend = FakeBackend()
    runtime = HarnessRuntime(backend)
    state, _ = asyncio.run(runtime.run(project, "paper audit"))
    with pytest.raises(ValueError, match="no pending"):
        asyncio.run(runtime.resume(project, state.run_id))


def test_context_compaction_is_traced(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AI_RESEARCH_MEMORY_ROOT", str(tmp_path / "global-memory"))
    project = tmp_path / "project"
    project.mkdir()
    backend = FakeBackend(
        trace_events=[
            (
                "contextCompaction/completed",
                {"thread_id": "thread-1", "turn_id": "turn-1", "method": "contextCompaction/completed"},
            )
        ]
    )
    state, _ = asyncio.run(HarnessRuntime(backend).run(project, "long context"))
    trace = (RunStore(project).run_dir(state.run_id) / "trace.jsonl").read_text(encoding="utf-8")
    assert "contextCompaction/completed" in trace


def test_interrupted_process_recovers_thread_from_atomic_trace(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AI_RESEARCH_MEMORY_ROOT", str(tmp_path / "global-memory"))
    project = tmp_path / "project"
    project.mkdir()
    with pytest.raises(asyncio.CancelledError):
        asyncio.run(HarnessRuntime(InterruptingBackend()).run(project, "recover me"))
    run_store = RunStore(project)
    run_id = next(run_store.runs.iterdir()).name
    pending = run_store.load_state(run_id)
    assert pending.goals[0].status == GoalStatus.pending
    completed, _ = asyncio.run(HarnessRuntime(FakeBackend()).resume(project, run_id))
    assert completed.goals[0].status == GoalStatus.done
    assert completed.thread_id == "thread-recover"


def test_artifact_claim_requires_current_project_file_and_hash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AI_RESEARCH_MEMORY_ROOT", str(tmp_path / "global-memory"))
    project = tmp_path / "project"
    project.mkdir()
    ungrounded, _ = asyncio.run(HarnessRuntime(ReceiptBackend(None)).run(project, "create output"))
    assert ungrounded.goals[0].status == GoalStatus.blocked

    artifact = project / "result.txt"
    artifact.write_text("verified", encoding="utf-8")
    receipt = Receipt(
        receipt_id="artifact",
        kind="artifact",
        authoritative=True,
        status="completed",
        artifact_path=str(artifact),
        sha256=hashlib.sha256(b"verified").hexdigest(),
        details={"source": "item/completed:fileChange"},
    )
    grounded, _ = asyncio.run(HarnessRuntime(ReceiptBackend(receipt)).run(project, "create output"))
    assert grounded.goals[0].status == GoalStatus.done


def test_workspace_write_rejects_broad_project_root() -> None:
    with pytest.raises(ValueError, match="too broad"):
        initialize_project(Path("/"))

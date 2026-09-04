from __future__ import annotations

from pathlib import Path

import pytest

from ai_research_harness.runtime import ASSETS
from ai_research_harness.storage import MemoryStore, RunStore


def test_memory_store_bootstrap_promote_tamper_and_rollback(tmp_path: Path) -> None:
    store = MemoryStore(tmp_path / "memory", "global")
    assert store.bootstrap(ASSETS / "memory" / "base-v1") == "base-v1"
    store.install_bundle(ASSETS / "memory" / "recuris-v1", "recuris-v1", candidate=True)
    store.record_gate_decision("recuris-v1", "NOT_ADMITTED", "/tmp/report.json")
    assert '"status": "NOT_ADMITTED"' in store.ledger.read_text(encoding="utf-8")
    store.materialize_candidate_as_version("recuris-v1")
    store.promote("recuris-v1", reason="test")
    assert store.champion_id() == "recuris-v1"
    store.promote("base-v1", reason="explicit test rollback", action="rollback")
    assert store.champion_id() == "base-v1"
    (store.version_dir("base-v1") / "working_memory.yaml").write_text("tampered: true\n", encoding="utf-8")
    with pytest.raises(ValueError, match="hash mismatch"):
        store.champion_id()


def test_version_id_rejects_path_traversal(tmp_path: Path) -> None:
    store = MemoryStore(tmp_path / "memory", "global")
    with pytest.raises(ValueError, match="invalid"):
        store.version_dir("../../escape")


def test_last_thread_id_skips_truncated_trace_tail(tmp_path: Path) -> None:
    project = tmp_path / "project"
    project.mkdir()
    store = RunStore(project)
    store.initialize()
    run = store.runs / "run-one"
    run.mkdir()
    (run / "trace.jsonl").write_text(
        '{"thread_id":"thread-stable"}\n{"thread_id":',
        encoding="utf-8",
    )
    assert store.last_thread_id("run-one") == "thread-stable"

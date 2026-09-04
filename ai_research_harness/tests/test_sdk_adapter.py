import hashlib
import asyncio
from pathlib import Path

import pytest

from ai_research_harness.sdk_adapter import OpenAICodexBackend, _final_response, _receipts, sdk_codex_config


def test_final_response_prefers_final_phase() -> None:
    items = [
        {"type": "agentMessage", "text": "progress", "phase": "commentary"},
        {"type": "agentMessage", "text": "final", "phase": "final_answer"},
    ]
    assert _final_response(items) == "final"
    assert _final_response(items[:1]) is None


def test_completed_items_become_authoritative_receipts() -> None:
    receipts = _receipts(
        [
            {"type": "commandExecution", "id": "cmd", "status": "completed", "exitCode": 0},
            {"type": "fileChange", "id": "patch", "status": "declined"},
        ],
        "turn-1",
        "completed",
    )
    assert receipts[0].authoritative
    assert receipts[0].status == "completed"
    assert receipts[1].status == "declined"
    assert receipts[-1].kind == "turn"


def test_receipts_reject_in_progress_and_nonzero_command() -> None:
    receipts = _receipts(
        [
            {"type": "commandExecution", "id": "running", "status": "inProgress", "exitCode": None},
            {"type": "commandExecution", "id": "bad", "status": "completed", "exitCode": 7},
        ],
        "turn-1",
        "interrupted",
    )
    assert [receipt.status for receipt in receipts] == ["interrupted", "failed", "interrupted"]


def test_completed_file_change_gets_real_hash_and_ignores_missing_file(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact.txt"
    artifact.write_text("real bytes", encoding="utf-8")
    receipts = _receipts(
        [
            {
                "type": "fileChange",
                "id": "patch",
                "status": "completed",
                "changes": [{"path": "artifact.txt"}, {"path": "partial.tmp"}],
            }
        ],
        "turn-1",
        "completed",
        tmp_path,
    )
    artifacts = [receipt for receipt in receipts if receipt.kind == "artifact"]
    assert len(artifacts) == 1
    assert artifacts[0].artifact_path == str(artifact)
    assert artifacts[0].sha256 == hashlib.sha256(b"real bytes").hexdigest()


def test_sdk_config_disables_inherited_tools_and_experimental_api() -> None:
    config = sdk_codex_config()
    assert config.experimental_api is False
    assert "notify=[]" in config.config_overrides
    assert "mcp_servers={}" in config.config_overrides
    assert "plugins={}" in config.config_overrides
    assert 'web_search="disabled"' in config.config_overrides


def test_unknown_sandbox_cannot_escalate_to_workspace_write(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unsupported sandbox"):
        asyncio.run(
            OpenAICodexBackend().execute(
                prompt="x",
                cwd=tmp_path,
                developer_instructions="x",
                skill_paths=[],
                sandbox="danger_full",
            )
        )

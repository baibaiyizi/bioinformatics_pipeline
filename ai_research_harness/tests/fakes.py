from __future__ import annotations

from pathlib import Path
from typing import Any

from ai_research_harness.models import Receipt
from ai_research_harness.sdk_adapter import SDKRunResult


class FakeBackend:
    def __init__(
        self,
        response: str = "完成并保留证据边界。",
        status: str = "completed",
        trace_events: list[tuple[str, dict[str, Any]]] | None = None,
    ) -> None:
        self.response = response
        self.status = status
        self.trace_events = trace_events or []
        self.calls: list[dict[str, Any]] = []

    async def execute(
        self,
        *,
        prompt: str,
        cwd: Path,
        developer_instructions: str,
        skill_paths: list[tuple[str, Path]],
        thread_id: str | None = None,
        sandbox: str = "read_only",
        output_schema: dict[str, Any] | None = None,
        trace_callback=None,
    ) -> SDKRunResult:
        self.calls.append(
            {
                "prompt": prompt,
                "cwd": cwd,
                "developer_instructions": developer_instructions,
                "skill_paths": skill_paths,
                "thread_id": thread_id,
                "sandbox": sandbox,
                "output_schema": output_schema,
            }
        )
        if trace_callback is not None:
            events = [
                *self.trace_events,
                (
                    "item/completed",
                    {
                        "thread_id": thread_id or "thread-1",
                        "turn_id": f"turn-{len(self.calls)}",
                        "item_id": "msg-1",
                        "item_type": "agentMessage",
                        "item_status": "completed",
                    },
                ),
            ]
            for method, summary in events:
                result = trace_callback(method, summary)
                if result is not None:
                    await result
        turn_id = f"turn-{len(self.calls)}"
        return SDKRunResult(
            thread_id=thread_id or "thread-1",
            turn_id=turn_id,
            status=self.status,
            final_response=self.response,
            receipts=[
                Receipt(
                    receipt_id=f"{turn_id}:turn",
                    kind="turn",
                    authoritative=True,
                    status=self.status,
                    item_id=turn_id,
                    details={"source": "turn/completed"},
                )
            ],
            usage={"totalTokens": 12},
            duration_ms=10,
        )

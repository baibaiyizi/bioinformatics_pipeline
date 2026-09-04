from __future__ import annotations

import asyncio
import dataclasses
import hashlib
from collections.abc import Awaitable, Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol

from pydantic import BaseModel

from .models import Receipt


TraceCallback = Callable[[str, dict[str, Any]], Awaitable[None] | None]
SDK_STARTUP_TIMEOUT_SECONDS = 45


def sdk_codex_config() -> Any:
    """Build an isolated stable-SDK configuration without changing user config."""
    from openai_codex import CodexConfig

    return CodexConfig(
        config_overrides=(
            "notify=[]",
            "mcp_servers={}",
            "plugins={}",
            'web_search="disabled"',
            "features.memories=false",
        ),
        client_name="ai_research_harness",
        client_title="Codex AI Research Harness",
        experimental_api=False,
    )


@dataclass(slots=True)
class SDKRunResult:
    thread_id: str
    turn_id: str
    status: str
    final_response: str | None
    receipts: list[Receipt]
    usage: dict[str, Any] = field(default_factory=dict)
    duration_ms: int | None = None
    completed_items: list[dict[str, Any]] = field(default_factory=list)


class AgentBackend(Protocol):
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
        trace_callback: TraceCallback | None = None,
    ) -> SDKRunResult: ...


def _jsonable(value: Any) -> Any:
    if isinstance(value, BaseModel):
        return value.model_dump(mode="json", by_alias=True)
    if dataclasses.is_dataclass(value):
        return dataclasses.asdict(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def _unwrap_item(item: dict[str, Any]) -> dict[str, Any]:
    root = item.get("root")
    return root if isinstance(root, dict) else item


def _item_type(item: dict[str, Any]) -> str:
    return str(_unwrap_item(item).get("type", "unknown"))


def _item_status(item: dict[str, Any]) -> str:
    status = str(_unwrap_item(item).get("status", "completed"))
    mapping = {"inProgress": "interrupted", "completed": "completed", "failed": "failed", "declined": "declined"}
    return mapping.get(status, "failed")


def _final_response(items: list[dict[str, Any]]) -> str | None:
    for item in reversed(items):
        value = _unwrap_item(item)
        if value.get("type") != "agentMessage":
            continue
        text = value.get("text")
        if not isinstance(text, str):
            continue
        phase = value.get("phase")
        if phase in {"final_answer", "finalAnswer"}:
            return text
    return None


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _receipts(items: list[dict[str, Any]], turn_id: str, turn_status: str, cwd: Path | None = None) -> list[Receipt]:
    receipts: list[Receipt] = []
    kind_map = {
        "agentMessage": "agent_message",
        "commandExecution": "command_execution",
        "fileChange": "file_change",
        "mcpToolCall": "mcp_tool_call",
        "webSearch": "web_search",
    }
    for item in items:
        value = _unwrap_item(item)
        item_type = str(value.get("type", ""))
        if item_type not in kind_map:
            continue
        item_id = value.get("id")
        status = _item_status(value)
        exit_code = value.get("exitCode")
        if item_type == "commandExecution" and status == "completed":
            if not isinstance(exit_code, int) or exit_code != 0:
                status = "failed"
        receipts.append(
            Receipt(
                receipt_id=f"{turn_id}:{item_id or len(receipts)}",
                kind=kind_map[item_type],
                authoritative=True,
                status=status,
                item_id=str(item_id) if item_id is not None else None,
                details={"source": "item/completed", "item_type": item_type, "exit_code": exit_code},
            )
        )
        if item_type == "fileChange" and status == "completed" and cwd is not None:
            root = cwd.resolve()
            for index, change in enumerate(value.get("changes", [])):
                if not isinstance(change, dict) or not isinstance(change.get("path"), str):
                    continue
                path = Path(change["path"])
                target = (path if path.is_absolute() else root / path).resolve()
                if target != root and root not in target.parents:
                    continue
                if not target.is_file():
                    continue
                receipts.append(
                    Receipt(
                        receipt_id=f"{turn_id}:{item_id or 'file'}:{index}",
                        kind="artifact",
                        authoritative=True,
                        status="completed",
                        item_id=str(item_id) if item_id is not None else None,
                        artifact_path=str(target),
                        sha256=_file_digest(target),
                        details={"source": "item/completed:fileChange"},
                    )
                )
    normalized_turn_status = turn_status if turn_status in {"completed", "failed", "interrupted"} else "failed"
    receipts.append(
        Receipt(
            receipt_id=f"{turn_id}:turn",
            kind="turn",
            authoritative=True,
            status=normalized_turn_status,
            item_id=turn_id,
            details={"source": "turn/completed"},
        )
    )
    return receipts


class OpenAICodexBackend:
    """Thin adapter over the stable openai-codex Python SDK surface."""

    def __init__(self, *, model: str = "gpt-5.6-sol", effort: str = "xhigh") -> None:
        self.model = model
        self.effort = effort
        self._codex: Any = None

    async def __aenter__(self) -> "OpenAICodexBackend":
        from openai_codex import AsyncCodex

        self._codex = AsyncCodex(config=sdk_codex_config())
        try:
            async with asyncio.timeout(SDK_STARTUP_TIMEOUT_SECONDS):
                await self._codex.__aenter__()
        except BaseException:
            await self._codex.close()
            self._codex = None
            raise
        return self

    async def __aexit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        if self._codex is not None:
            await self._codex.__aexit__(exc_type, exc, traceback)
            self._codex = None

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
        trace_callback: TraceCallback | None = None,
    ) -> SDKRunResult:
        from openai_codex import ApprovalMode, Sandbox, SkillInput, TextInput

        if sandbox not in {"read_only", "workspace_write"}:
            raise ValueError(f"unsupported sandbox mode: {sandbox}")
        if self._codex is None:
            raise RuntimeError("OpenAICodexBackend must be used as an async context manager")
        sandbox_value = Sandbox.read_only if sandbox == "read_only" else Sandbox.workspace_write
        if thread_id is None:
            thread = await self._codex.thread_start(
                approval_mode=ApprovalMode.deny_all,
                cwd=str(cwd.resolve()),
                developer_instructions=developer_instructions,
                model=self.model,
                sandbox=sandbox_value,
                config={"model_reasoning_effort": self.effort},
                service_name="ai_research_harness",
            )
        else:
            thread = await self._codex.thread_resume(
                thread_id,
                approval_mode=ApprovalMode.deny_all,
                cwd=str(cwd.resolve()),
                developer_instructions=developer_instructions,
                model=self.model,
                sandbox=sandbox_value,
                config={"model_reasoning_effort": self.effort},
            )

        inputs: list[Any] = [TextInput(prompt)]
        inputs.extend(SkillInput(name=name, path=str(path.resolve())) for name, path in skill_paths)
        handle = await thread.turn(
            inputs,
            approval_mode=ApprovalMode.deny_all,
            model=self.model,
            sandbox=sandbox_value,
            output_schema=output_schema,
        )
        items: list[dict[str, Any]] = []
        usage: dict[str, Any] = {}
        turn_payload: dict[str, Any] = {}
        async for event in handle.stream():
            payload = _jsonable(event.payload)
            if not isinstance(payload, dict):
                payload = {"value": payload}
            if trace_callback is not None:
                summary = {"method": event.method, "thread_id": thread.id, "turn_id": handle.id}
                if event.method == "item/completed":
                    item = payload.get("item", {})
                    summary.update(
                        {
                            "item_id": item.get("id") if isinstance(item, dict) else None,
                            "item_type": _item_type(item) if isinstance(item, dict) else "unknown",
                            "item_status": _item_status(item) if isinstance(item, dict) else "unknown",
                        }
                    )
                    value = _unwrap_item(item) if isinstance(item, dict) else {}
                    if summary["item_type"] == "commandExecution":
                        summary["action"] = {
                            "command": value.get("command"),
                            "cwd": value.get("cwd"),
                            "exit_code": value.get("exitCode"),
                        }
                    elif summary["item_type"] == "fileChange":
                        summary["action"] = {
                            "paths": [
                                change.get("path")
                                for change in value.get("changes", [])
                                if isinstance(change, dict)
                            ]
                        }
                    elif summary["item_type"] == "webSearch":
                        summary["action"] = {"query": value.get("query")}
                    elif summary["item_type"] == "mcpToolCall":
                        summary["action"] = {
                            "server": value.get("server"),
                            "tool": value.get("tool"),
                        }
                elif event.method == "turn/completed":
                    summary["turn"] = payload.get("turn", {})
                elif "usage" in event.method.casefold():
                    summary["usage"] = payload
                callback_result = trace_callback(event.method, summary)
                if callback_result is not None:
                    await callback_result
            if event.method == "item/completed" and isinstance(payload.get("item"), dict):
                items.append(payload["item"])
            elif event.method == "thread/tokenUsage/updated":
                usage = payload.get("tokenUsage", payload)
            elif event.method == "turn/completed":
                turn_payload = payload.get("turn", payload)

        status = str(turn_payload.get("status", "failed"))
        duration = turn_payload.get("durationMs")
        final = _final_response(items)
        return SDKRunResult(
            thread_id=thread.id,
            turn_id=handle.id,
            status=status,
            final_response=final,
            receipts=_receipts(items, handle.id, status, cwd),
            usage=usage if isinstance(usage, dict) else {},
            duration_ms=int(duration) if isinstance(duration, (int, float)) else None,
            completed_items=items,
        )

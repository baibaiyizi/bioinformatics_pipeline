from __future__ import annotations

import asyncio
import hashlib
import json
import os
import platform
from collections import Counter
from pathlib import Path
from typing import Any

from .gate import evaluate_gate, grade_response, validate_candidate_bundle
from .memory import MemorySelector, load_bundle_cards
from .models import ArmOutcome, GateDecision, Goal, QualificationTask, WorkingMemory
from .sdk_adapter import AgentBackend, OpenAICodexBackend


EVALUATOR_INSTRUCTIONS = """
你在执行冻结的 AI/计算机科学研究资格题。直接回答题目，明确关键边界和理由。不得查找本机 qualification 文件、gold、rubric 或另一 arm 输出；不得修改文件。模型、工具、权限和预算在两 arm 间必须保持一致。
固定预算：单次回答只允许一个 turn、零工具动作、零任务内重试，最终回答不超过 6000 字符，总 token usage 不超过 100000。超过即判失败。
""".strip()

QUALIFICATION_TURN_TIMEOUT_SECONDS = 120
QUALIFICATION_TOKEN_LIMIT = 100_000
QUALIFICATION_RESPONSE_CHAR_LIMIT = 6_000

REVIEW_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "required": ["a_pass", "b_pass", "a_reason", "b_reason"],
    "properties": {
        "a_pass": {"type": "boolean"},
        "b_pass": {"type": "boolean"},
        "a_reason": {"type": "string"},
        "b_reason": {"type": "string"},
    },
}


class QualificationRunner:
    def __init__(self, backend: AgentBackend | None = None, *, concurrency: int = 3) -> None:
        self.backend = backend
        self.concurrency = concurrency
        self.selector = MemorySelector()
        self._turn_semaphore = asyncio.Semaphore(concurrency)

    @staticmethod
    def _state(task: QualificationTask) -> WorkingMemory:
        return WorkingMemory(
            run_id=f"qualification-{task.task_id}",
            task_id=task.task_id,
            project_root=str(Path.cwd().resolve()),
            phase="evaluation",
            goals=[Goal(goal_id="answer", objective=task.prompt, acceptance_checks=["agent_response_present"])],
            active_constraints=[task.domain, "matched budget", "frozen qualification"],
        )

    async def _arm(self, backend: AgentBackend, task: QualificationTask, arm: str, bundle: Path, cwd: Path) -> ArmOutcome:
        _, cards = load_bundle_cards(bundle)
        selected = self.selector.select(self._state(task), cards)
        try:
            async with self._turn_semaphore:
                async with asyncio.timeout(QUALIFICATION_TURN_TIMEOUT_SECONDS):
                    result = await backend.execute(
                        prompt=task.prompt,
                        cwd=cwd,
                        developer_instructions=EVALUATOR_INSTRUCTIONS,
                        skill_paths=[(card.skill_id, card.path) for card in selected],
                        sandbox="read_only",
                    )
        except TimeoutError:
            return ArmOutcome(
                task_id=task.task_id,
                arm=arm,
                passed=False,
                deterministic_passed=False,
                response="",
                reasons=[f"turn exceeded {QUALIFICATION_TURN_TIMEOUT_SECONDS}s budget"],
            )
        outcome = grade_response(
            task,
            arm,
            result.final_response or "",
            result.usage,
            result.duration_ms,
            trace={
                "thread_id": result.thread_id,
                "turn_id": result.turn_id,
                "status": result.status,
                "completed_items": result.completed_items,
                "input_sha256": hashlib.sha256(task.prompt.encode("utf-8")).hexdigest(),
            },
        )
        violations = _budget_violations(result)
        if result.status != "completed":
            violations.append(f"turn status={result.status}")
        if violations:
            outcome.passed = False
            outcome.deterministic_passed = False
            outcome.reasons.extend(violations)
        return outcome

    async def _semantic_review(
        self,
        backend: AgentBackend,
        task: QualificationTask,
        incumbent: ArmOutcome,
        candidate: ArmOutcome,
        cwd: Path,
        *,
        reverse: bool,
    ) -> tuple[bool, bool, list[str], dict[str, Any]]:
        first, second = (candidate, incumbent) if reverse else (incumbent, candidate)
        prompt = (
            "按冻结 rubric 独立判断两个匿名回答是否通过。不要猜测来源或身份。\n\n"
            f"题目：{task.prompt}\n\nRubric：{task.rubric}\n\n"
            f"回答 A：\n{first.response}\n\n回答 B：\n{second.response}"
        )
        try:
            async with self._turn_semaphore:
                async with asyncio.timeout(QUALIFICATION_TURN_TIMEOUT_SECONDS):
                    result = await backend.execute(
                        prompt=prompt,
                        cwd=cwd,
                        developer_instructions="你是只读 fresh reviewer，只输出 schema 要求的 JSON；不得使用工具。",
                        skill_paths=[],
                        sandbox="read_only",
                        output_schema=REVIEW_SCHEMA,
                    )
        except TimeoutError:
            return False, False, ["semantic reviewer exceeded fixed time budget"], {"status": "timeout"}
        reviewer_trace = {
            "thread_id": result.thread_id,
            "turn_id": result.turn_id,
            "status": result.status,
            "completed_items": result.completed_items,
            "response": result.final_response,
            "usage": result.usage,
            "duration_ms": result.duration_ms,
        }
        violations = _budget_violations(result)
        if result.status != "completed" or violations:
            return False, False, [f"semantic reviewer protocol failure: {violations or result.status}"], reviewer_trace
        try:
            payload = json.loads(result.final_response or "")
            a_pass = bool(payload["a_pass"])
            b_pass = bool(payload["b_pass"])
            reasons = [str(payload.get("a_reason", "")), str(payload.get("b_reason", ""))]
        except (json.JSONDecodeError, KeyError, TypeError) as exc:
            return False, False, [f"semantic reviewer output invalid: {exc}"], reviewer_trace
        if reverse:
            return b_pass, a_pass, reasons, reviewer_trace
        return a_pass, b_pass, reasons, reviewer_trace

    async def _pair(
        self,
        backend: AgentBackend,
        task: QualificationTask,
        incumbent_bundle: Path,
        candidate_bundle: Path,
        cwd: Path,
        *,
        repeat: int,
    ) -> tuple[ArmOutcome, ArmOutcome]:
        if repeat % 2:
            candidate, incumbent = await asyncio.gather(
                self._arm(backend, task, "candidate", candidate_bundle, cwd),
                self._arm(backend, task, "incumbent", incumbent_bundle, cwd),
            )
        else:
            incumbent, candidate = await asyncio.gather(
                self._arm(backend, task, "incumbent", incumbent_bundle, cwd),
                self._arm(backend, task, "candidate", candidate_bundle, cwd),
            )
        if task.semantic_review:
            reverse = int(hashlib.sha256(f"{task.task_id}:{repeat}".encode()).hexdigest(), 16) % 2 == 1
            inc_semantic, cand_semantic, reasons, reviewer_trace = await self._semantic_review(
                backend, task, incumbent, candidate, cwd, reverse=reverse
            )
            incumbent.semantic_passed = inc_semantic
            candidate.semantic_passed = cand_semantic
            incumbent.passed = incumbent.deterministic_passed and inc_semantic
            candidate.passed = candidate.deterministic_passed and cand_semantic
            incumbent.reasons.extend(reasons[:1])
            candidate.reasons.extend(reasons[1:])
            incumbent.trace["semantic_review"] = reviewer_trace
            candidate.trace["semantic_review"] = reviewer_trace
        return incumbent, candidate

    async def run(
        self,
        *,
        tasks: list[QualificationTask],
        incumbent_bundle: Path,
        candidate_bundle: Path,
        cwd: Path,
    ) -> tuple[list[ArmOutcome], list[ArmOutcome], dict[str, list[dict[str, Any]]]]:
        owns_backend = self.backend is None
        backend = self.backend or OpenAICodexBackend()
        trials: dict[str, list[dict[str, Any]]] = {}

        async def run_task(task: QualificationTask) -> tuple[ArmOutcome, ArmOutcome]:
            first = await self._pair(backend, task, incumbent_bundle, candidate_bundle, cwd, repeat=0)
            pairs = [first]
            if first[0].passed != first[1].passed:
                pairs.append(await self._pair(backend, task, incumbent_bundle, candidate_bundle, cwd, repeat=1))
                pairs.append(await self._pair(backend, task, incumbent_bundle, candidate_bundle, cwd, repeat=2))
            inc_majority = Counter(pair[0].passed for pair in pairs).most_common(1)[0][0]
            cand_majority = Counter(pair[1].passed for pair in pairs).most_common(1)[0][0]
            incumbent = pairs[0][0].model_copy(update={"passed": inc_majority})
            candidate = pairs[0][1].model_copy(update={"passed": cand_majority})
            trials[task.task_id] = [
                {"incumbent": pair[0].model_dump(mode="json"), "candidate": pair[1].model_dump(mode="json")}
                for pair in pairs
            ]
            return incumbent, candidate

        try:
            if owns_backend:
                await backend.__aenter__()  # type: ignore[attr-defined]
            pairs = await asyncio.gather(*(run_task(task) for task in tasks))
        finally:
            if owns_backend:
                await backend.__aexit__(None, None, None)  # type: ignore[attr-defined]
        return [pair[0] for pair in pairs], [pair[1] for pair in pairs], trials


async def run_live_gate(
    *,
    tasks: list[QualificationTask],
    incumbent_bundle: Path,
    candidate_bundle: Path,
    cwd: Path,
    source_failure_ids: list[str],
    concurrency: int = 3,
) -> tuple[GateDecision, dict[str, Any]]:
    runner = QualificationRunner(concurrency=concurrency)
    incumbent, candidate, trials = await runner.run(
        tasks=tasks,
        incumbent_bundle=incumbent_bundle,
        candidate_bundle=candidate_bundle,
        cwd=cwd,
    )
    decision = evaluate_gate(
        candidate_id=candidate_bundle.name,
        incumbent_id=incumbent_bundle.name,
        tasks=tasks,
        incumbent=incumbent,
        candidate=candidate,
        validation_checks=validate_candidate_bundle(candidate_bundle, tasks),
        source_failure_ids=source_failure_ids,
    )
    metrics = _gate_metrics(incumbent_bundle, candidate_bundle, incumbent, candidate)
    decision = decision.model_copy(update={"metrics": metrics})
    return decision, {
        "incumbent": [item.model_dump(mode="json") for item in incumbent],
        "candidate": [item.model_dump(mode="json") for item in candidate],
        "trials": trials,
        "metrics": metrics,
    }


def _total_tokens(usage: dict[str, Any]) -> int:
    direct = usage.get("totalTokens")
    if isinstance(direct, int):
        return direct
    total = usage.get("total")
    if isinstance(total, dict) and isinstance(total.get("totalTokens"), int):
        return int(total["totalTokens"])
    return 0


def _budget_violations(result) -> list[str]:
    violations: list[str] = []
    action_types = {"commandExecution", "fileChange", "mcpToolCall", "webSearch", "collabAgentToolCall"}
    actions = [
        item
        for item in result.completed_items
        if str(item.get("root", item).get("type", "")) in action_types
    ]
    if actions:
        violations.append(f"tool action budget exceeded: {len(actions)} > 0")
    total_tokens = _total_tokens(result.usage)
    if total_tokens > QUALIFICATION_TOKEN_LIMIT:
        violations.append(f"token budget exceeded: {total_tokens} > {QUALIFICATION_TOKEN_LIMIT}")
    if len(result.final_response or "") > QUALIFICATION_RESPONSE_CHAR_LIMIT:
        violations.append(
            f"response budget exceeded: {len(result.final_response or '')} > {QUALIFICATION_RESPONSE_CHAR_LIMIT} chars"
        )
    return violations


def _bundle_size(bundle: Path) -> int:
    return sum(path.stat().st_size for path in bundle.rglob("*") if path.is_file())


def _gate_metrics(
    incumbent_bundle: Path,
    candidate_bundle: Path,
    incumbent: list[ArmOutcome],
    candidate: list[ArmOutcome],
) -> dict[str, Any]:
    _, incumbent_cards = load_bundle_cards(incumbent_bundle)
    _, candidate_cards = load_bundle_cards(candidate_bundle)
    incumbent_bytes = _bundle_size(incumbent_bundle)
    candidate_bytes = _bundle_size(candidate_bundle)
    return {
        "fixed_budget": {
            "model": "gpt-5.6-sol",
            "reasoning_effort": "xhigh",
            "model_seed": "not_exposed_by_openai-codex-0.147.0",
            "bootstrap_seed": 260824876,
            "turns_per_attempt": 1,
            "task_retries": 0,
            "tool_actions": 0,
            "token_limit": QUALIFICATION_TOKEN_LIMIT,
            "response_char_limit": QUALIFICATION_RESPONSE_CHAR_LIMIT,
            "timeout_seconds": QUALIFICATION_TURN_TIMEOUT_SECONDS,
            "host": {
                "platform": platform.platform(),
                "machine": platform.machine(),
                "logical_cpu_count": os.cpu_count(),
                "python": platform.python_version(),
            },
        },
        "incumbent": {
            "tokens": sum(_total_tokens(item.usage) for item in incumbent),
            "duration_ms": sum(item.duration_ms or 0 for item in incumbent),
        },
        "candidate": {
            "tokens": sum(_total_tokens(item.usage) for item in candidate),
            "duration_ms": sum(item.duration_ms or 0 for item in candidate),
        },
        "memory_delta": {
            "cards": len(candidate_cards) - len(incumbent_cards),
            "bytes": candidate_bytes - incumbent_bytes,
        },
    }

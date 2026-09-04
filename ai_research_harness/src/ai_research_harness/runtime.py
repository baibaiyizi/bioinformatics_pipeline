from __future__ import annotations

import json
import os
import tomllib
import uuid
from pathlib import Path

from .checkers import default_registry
from .memory import MemorySelector, load_bundle_cards, load_bundle_settings
from .models import Goal, GoalStatus, TraceRecord, WorkingMemory
from .sdk_adapter import AgentBackend, OpenAICodexBackend, SDKRunResult
from .storage import MemoryStore, RunStore, atomic_write, model_sha256, sha256_file


PACKAGE_ROOT = Path(__file__).resolve().parent
ASSETS = PACKAGE_ROOT / "assets"
DEFAULT_AGENT_CONFIG = Path("/home/h1028/.codex/agents/ai_research.toml")
ARTIFACT_CLAIM_MARKERS = (
    "创建文件",
    "生成文件",
    "写入文件",
    "产物位于",
    "created file",
    "wrote file",
    "artifact:",
)


def global_memory_root() -> Path:
    configured = os.environ.get("AI_RESEARCH_MEMORY_ROOT")
    return Path(configured).expanduser() if configured else Path.home() / ".codex" / "ai_research_memory"


def load_orchestrator_instructions(path: Path = DEFAULT_AGENT_CONFIG) -> str:
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    return str(data["developer_instructions"])


def ensure_seed_memories(store: MemoryStore) -> tuple[str, str]:
    base = ASSETS / "memory" / "base-v1"
    recur = ASSETS / "memory" / "recuris-v1"
    incumbent = store.bootstrap(base)
    store.install_bundle(recur, "recuris-v1", candidate=True)
    return incumbent, "recuris-v1"


def initialize_project(project_root: Path) -> dict[str, str]:
    project_root = project_root.resolve()
    if not project_root.is_dir():
        raise FileNotFoundError(f"project directory does not exist: {project_root}")
    if project_root in {Path("/"), Path.home(), Path.home() / ".codex"}:
        raise ValueError(f"project root is too broad for workspace-write: {project_root}")
    runs = RunStore(project_root)
    runs.initialize()
    anchor_suite = runs.anchors / "qualification.yaml"
    anchor_digest = runs.anchors / "qualification.sha256"
    if not anchor_suite.exists():
        source_suite = ASSETS / "evals" / "frozen_qualification.yaml"
        atomic_write(anchor_suite, source_suite.read_text(encoding="utf-8"))
        atomic_write(anchor_digest, sha256_file(anchor_suite) + "\n")
    memory = MemoryStore(global_memory_root(), "global")
    incumbent, candidate = ensure_seed_memories(memory)
    project_memory = MemoryStore(runs.project_memory, "project")
    project_champion = project_memory.bootstrap(ASSETS / "memory" / "project-base-v1")
    return {
        "project": str(project_root),
        "global_champion": incumbent,
        "project_champion": project_champion,
        "project_qualification": str(anchor_suite),
        "bundled_candidate": candidate,
    }


class HarnessRuntime:
    def __init__(self, backend: AgentBackend | None = None) -> None:
        self.backend = backend
        self.registry = default_registry()
        self.selector = MemorySelector()

    def _load_cards(self, project: RunStore) -> tuple[list, list, dict[str, dict]]:
        global_store = MemoryStore(global_memory_root(), "global")
        ensure_seed_memories(global_store)
        global_bundle = global_store.version_dir(global_store.champion_id())
        _, global_cards = load_bundle_cards(global_bundle)
        settings = load_bundle_settings(global_bundle)
        project_store = MemoryStore(project.project_memory, "project")
        project_cards: list = []
        if project_store.pointer.exists():
            project_bundle = project_store.version_dir(project_store.champion_id())
            _, project_cards = load_bundle_cards(project_bundle)
            project_settings = load_bundle_settings(project_bundle)
            for component in ("working_memory", "invocation_policy", "checkers"):
                settings[component].update(project_settings[component])
        return global_cards, project_cards, settings

    @staticmethod
    def _new_state(project_root: Path, prompt: str) -> WorkingMemory:
        run_id = f"run-{uuid.uuid4().hex[:12]}"
        task_id = f"task-{uuid.uuid5(uuid.NAMESPACE_URL, prompt).hex[:12]}"
        return WorkingMemory(
            run_id=run_id,
            task_id=task_id,
            project_root=str(project_root.resolve()),
            goals=[
                Goal(
                    goal_id="deliverable",
                    objective=prompt,
                    acceptance_checks=["agent_response_present", "no_blocked_marker"],
                )
            ],
            active_constraints=[
                "evidence bounded AI and computer science research",
                "no unsupported completion claim",
                "matched evaluation budgets",
            ],
        )

    async def run(self, project_root: Path, prompt: str) -> tuple[WorkingMemory, SDKRunResult]:
        project_root = project_root.resolve()
        initialize_project(project_root)
        store = RunStore(project_root)
        state = self._new_state(project_root, prompt)
        store.create_run(state, prompt)
        return await self._execute(store, state, prompt, resume=False)

    async def resume(self, project_root: Path, run_id: str) -> tuple[WorkingMemory, SDKRunResult]:
        store = RunStore(project_root.resolve())
        state = store.load_state(run_id)
        if not state.pending_goals():
            raise ValueError("run has no pending goal to resume")
        if state.thread_id is None:
            state.thread_id = store.last_thread_id(run_id)
            store.save_state(state)
        prompt = (store.run_dir(run_id) / "prompt.md").read_text(encoding="utf-8")
        continuation = "继续当前任务。先核对 working memory 中尚未完成的目标和已有 receipt，不重复已完成动作。\n\n原始任务：\n" + prompt
        return await self._execute(store, state, continuation, resume=True)

    async def _execute(self, store: RunStore, state: WorkingMemory, prompt: str, *, resume: bool) -> tuple[WorkingMemory, SDKRunResult]:
        before = model_sha256(state)
        global_cards, project_cards, settings = self._load_cards(store)
        policy = settings["invocation_policy"]
        selector = MemorySelector(max_cards=int(policy["max_cards"]), max_tokens=int(policy["max_tokens"]))
        selected = selector.select(state, global_cards, project_cards)
        configured_checks = list(settings["checkers"].get("goal_checks", []))
        if not configured_checks or any(checker not in self.registry.ids for checker in configured_checks):
            raise ValueError("champion references an unknown or empty checker set")
        state.pending_goals()[0].acceptance_checks = configured_checks
        state.selected_skill_ids = [card.skill_id for card in selected]
        store.save_state(state)
        store.trace(
            TraceRecord(
                run_id=state.run_id,
                event="skills_selected",
                thread_id=state.thread_id,
                state_before_sha256=before,
                state_after_sha256=model_sha256(state),
                selected_skill_ids=state.selected_skill_ids,
            )
        )

        instructions = load_orchestrator_instructions() + "\n\n本任务通过完整 Harness 启动，已明确授权 ai_research 按共享契约使用必要专科；最多三个并行专科。"
        instructions += "\n\n当前 working memory（只读上下文，状态由 Harness 提交）：\n" + json.dumps(
            state.model_dump(mode="json"), ensure_ascii=False, sort_keys=True
        )
        instructions += "\n\n当前 E/W/rho/C 声明式规范：\n" + json.dumps(settings, ensure_ascii=False, sort_keys=True)
        pre_commit_sha = model_sha256(state)

        async def trace_event(method: str, summary: dict) -> None:
            store.trace(
                TraceRecord(
                    run_id=state.run_id,
                    event=method,
                    thread_id=summary.get("thread_id") or state.thread_id,
                    turn_id=summary.get("turn_id"),
                    item_id=summary.get("item_id"),
                    selected_skill_ids=state.selected_skill_ids,
                    details=summary,
                )
            )

        owns_backend = self.backend is None
        backend = self.backend or OpenAICodexBackend()
        try:
            if owns_backend:
                await backend.__aenter__()  # type: ignore[attr-defined]
            result = await backend.execute(
                prompt=prompt,
                cwd=Path(state.project_root),
                developer_instructions=instructions,
                skill_paths=[(card.skill_id, card.path) for card in selected],
                thread_id=state.thread_id if resume else None,
                sandbox="workspace_write",
                trace_callback=trace_event,
            )
        except Exception as exc:
            goal = state.pending_goals()[0]
            goal.block(f"Codex runtime failed: {type(exc).__name__}: {exc}")
            store.save_state(state)
            store.trace(
                TraceRecord(
                    run_id=state.run_id,
                    event="run_blocked",
                    thread_id=state.thread_id,
                    state_before_sha256=pre_commit_sha,
                    state_after_sha256=model_sha256(state),
                    selected_skill_ids=state.selected_skill_ids,
                    commit_decision="blocked",
                    details={"error_type": type(exc).__name__, "error": str(exc)},
                )
            )
            raise
        finally:
            if owns_backend:
                await backend.__aexit__(None, None, None)  # type: ignore[attr-defined]

        state.thread_id = result.thread_id
        state.turn_ids.append(result.turn_id)
        if "phase_boundary" in policy.get("events", []):
            phase_before = model_sha256(state)
            state.phase = "verification"
            selected = selector.select(state, global_cards, project_cards)
            state.selected_skill_ids = [card.skill_id for card in selected]
            store.trace(
                TraceRecord(
                    run_id=state.run_id,
                    event="rho_selected",
                    thread_id=result.thread_id,
                    turn_id=result.turn_id,
                    state_before_sha256=phase_before,
                    state_after_sha256=model_sha256(state),
                    selected_skill_ids=state.selected_skill_ids,
                    details={"trigger": "phase_boundary", "phase": state.phase},
                )
            )
        goal = state.pending_goals()[0]
        checks = [self.registry.run(checker, response=result.final_response) for checker in goal.acceptance_checks]
        failed_checks = [check.checker_id for check in checks if not check.passed]
        if failed_checks and "checker_failure" in policy.get("events", []):
            failure_before = model_sha256(state)
            state.active_constraints.extend(
                marker
                for marker in (f"checker_failure:{checker}" for checker in failed_checks)
                if marker not in state.active_constraints
            )
            selected = selector.select(state, global_cards, project_cards)
            state.selected_skill_ids = [card.skill_id for card in selected]
            store.trace(
                TraceRecord(
                    run_id=state.run_id,
                    event="rho_selected",
                    thread_id=result.thread_id,
                    turn_id=result.turn_id,
                    state_before_sha256=failure_before,
                    state_after_sha256=model_sha256(state),
                    selected_skill_ids=state.selected_skill_ids,
                    details={"trigger": "checker_failure", "failed_checkers": failed_checks},
                )
            )
        response_text = (result.final_response or "").casefold()
        artifact_claimed = any(marker.casefold() in response_text for marker in ARTIFACT_CLAIM_MARKERS) or any(
            receipt.kind == "file_change" and receipt.status == "completed" for receipt in result.receipts
        )
        project_root = Path(state.project_root).resolve()
        artifact_grounded = False
        for receipt in result.receipts:
            if receipt.kind != "artifact" or not receipt.authoritative or receipt.status != "completed":
                continue
            if not receipt.artifact_path or not receipt.sha256:
                continue
            artifact_path = Path(receipt.artifact_path).resolve()
            if artifact_path != project_root and project_root not in artifact_path.parents:
                continue
            if artifact_path.is_file() and sha256_file(artifact_path) == receipt.sha256:
                artifact_grounded = True
                break
        if result.status == "completed":
            try:
                if artifact_claimed and not artifact_grounded:
                    raise ValueError("artifact claim lacks a real file SHA-256 receipt")
                goal.propose_done(result.receipts, checks)
                decision = "done"
            except ValueError as exc:
                goal.block(f"completion evidence rejected: {exc}", result.receipts)
                decision = "blocked"
        else:
            goal.block(f"turn ended with status {result.status}", result.receipts)
            decision = "blocked"
        store.save_state(state)
        if result.final_response is not None:
            atomic_write(store.run_dir(state.run_id) / "response.md", result.final_response)
        atomic_write(store.run_dir(state.run_id) / "usage.json", json.dumps(result.usage, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
        store.trace(
            TraceRecord(
                run_id=state.run_id,
                event="state_commit",
                thread_id=result.thread_id,
                turn_id=result.turn_id,
                state_before_sha256=pre_commit_sha,
                state_after_sha256=model_sha256(state),
                selected_skill_ids=state.selected_skill_ids,
                receipt={
                    "turn_status": result.status,
                    "items": [receipt.model_dump(mode="json") for receipt in result.receipts],
                },
                checker_results=[check.model_dump(mode="json") for check in checks],
                commit_decision=decision,
            )
        )
        return state, result


def state_exit_code(state: WorkingMemory) -> int:
    if all(goal.status == GoalStatus.done for goal in state.goals):
        return 0
    return 2

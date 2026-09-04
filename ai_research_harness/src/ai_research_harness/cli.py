from __future__ import annotations

import argparse
import asyncio
import importlib.metadata
import json
import sys
import tempfile
import tomllib
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .evolution import propose_candidate
from .gate import load_qualification_tasks, qualification_fingerprint, validate_candidate_bundle
from .models import MemoryManifest
from .qualification import run_live_gate
from .runtime import ASSETS, HarnessRuntime, ensure_seed_memories, global_memory_root, initialize_project, state_exit_code
from .sdk_adapter import SDK_STARTUP_TIMEOUT_SECONDS, sdk_codex_config
from .storage import MemoryStore, RunStore, append_jsonl, atomic_write, load_yaml


EXIT_OK = 0
EXIT_BLOCKED = 2
EXIT_GATE_FAILED = 3
EXIT_ENVIRONMENT = 4
AGENT_DIR = Path("/home/h1028/.codex/agents")
AI_LEARNING_DIR = Path("/home/h1028/workspace/ai_learning")
QUALIFICATION = ASSETS / "evals" / "frozen_qualification.yaml"


def _print(value: Any) -> None:
    print(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True))


def _scope_store(scope: str, project: Path | None) -> MemoryStore:
    if scope == "global":
        return MemoryStore(global_memory_root(), "global")
    if project is None:
        raise ValueError("project scope requires --project")
    return MemoryStore(RunStore(project.resolve()).project_memory, "project")


def _bootstrap_scope(store: MemoryStore, scope: str) -> str:
    if scope == "global":
        incumbent, _ = ensure_seed_memories(store)
        return incumbent
    return store.bootstrap(ASSETS / "memory" / "project-base-v1")


async def _doctor() -> int:
    checks: list[dict[str, Any]] = []
    try:
        version = importlib.metadata.version("openai-codex")
        checks.append({"check": "openai-codex-version", "passed": version == "0.147.0", "observed": version})
        runtime_version = importlib.metadata.version("openai-codex-cli-bin")
        checks.append(
            {"check": "bundled-codex-runtime-version", "passed": runtime_version == "0.147.0", "observed": runtime_version}
        )
    except importlib.metadata.PackageNotFoundError:
        checks.append({"check": "openai-codex-version", "passed": False, "observed": "not installed"})
        _print({"status": "blocked", "checks": checks})
        return EXIT_ENVIRONMENT

    configs: dict[str, Path] = {}
    for path in sorted(AGENT_DIR.glob("ai_research*.toml")):
        data = tomllib.loads(path.read_text(encoding="utf-8"))
        configs[str(data["name"])] = path
    checks.append({"check": "agent-configs", "passed": len(configs) == 12, "count": len(configs)})
    checks.append({"check": "shared-contract", "passed": (AGENT_DIR / "ai_research_contract.md").is_file()})
    checks.append({"check": "ai-learning", "passed": (AI_LEARNING_DIR / "index.md").is_file()})

    tasks = load_qualification_tasks(QUALIFICATION)
    domains = Counter(task.domain for task in tasks)
    expected = {"evidence": 8, "empirical_ml": 8, "reproducibility": 4, "systems": 3, "theory": 3, "security": 3, "hci": 3}
    checks.append({"check": "frozen-qualification", "passed": len(tasks) == 32 and dict(domains) == expected, "count": len(tasks), "sha256": qualification_fingerprint(QUALIFICATION)})

    global_store = MemoryStore(global_memory_root(), "global")
    if global_store.pointer.exists():
        try:
            champion = global_store.champion_id()
            candidates = []
            if global_store.candidates.is_dir():
                for candidate in sorted(path for path in global_store.candidates.iterdir() if path.is_dir()):
                    global_store.verify_bundle(candidate)
                    candidates.append(candidate.name)
            checks.append(
                {"check": "global-memory-hash", "passed": True, "champion": champion, "candidates": candidates}
            )
        except Exception as exc:
            checks.append({"check": "global-memory-hash", "passed": False, "error": str(exc)})
    else:
        checks.append({"check": "global-memory-hash", "passed": True, "observed": "not initialized"})

    try:
        from openai_codex import AsyncCodex

        async with asyncio.timeout(SDK_STARTUP_TIMEOUT_SECONDS):
            async with AsyncCodex(config=sdk_codex_config()) as codex:
                account = await codex.account()
                models = await codex.models(include_hidden=True)
        account_json = account.model_dump(mode="json") if hasattr(account, "model_dump") else repr(account)
        models_json = models.model_dump(mode="json") if hasattr(models, "model_dump") else repr(models)
        available = "gpt-5.6-sol" in json.dumps(models_json, ensure_ascii=False)
        authenticated = "notLoggedIn" not in json.dumps(account_json, ensure_ascii=False)
        checks.append({"check": "codex-authentication", "passed": authenticated})
        checks.append({"check": "required-model", "passed": available, "model": "gpt-5.6-sol"})
    except Exception as exc:
        checks.append({"check": "codex-sdk-runtime", "passed": False, "error": f"{type(exc).__name__}: {exc}"})

    passed = all(check["passed"] for check in checks)
    _print({"status": "complete" if passed else "blocked", "checks": checks})
    return EXIT_OK if passed else EXIT_ENVIRONMENT


async def _run_gate(candidate_id: str, scope: str, project: Path | None, concurrency: int) -> int:
    store = _scope_store(scope, project)
    incumbent_id = _bootstrap_scope(store, scope)
    if scope == "global" and candidate_id == "recuris-v1":
        ensure_seed_memories(store)
    candidate_bundle = store.candidate_dir(candidate_id)
    if not candidate_bundle.is_dir():
        raise FileNotFoundError(f"candidate not found: {candidate_id}")

    if scope == "global":
        suite_path = QUALIFICATION
    else:
        assert project is not None
        suite_path = RunStore(project.resolve()).anchors / "qualification.yaml"
        if not suite_path.is_file():
            _print({"status": "NOT_ADMITTED", "reason": "project regression anchors are absent", "expected": str(suite_path)})
            return EXIT_GATE_FAILED
        suite_anchor = suite_path.with_suffix(".sha256")
        if not suite_anchor.is_file() or suite_anchor.read_text(encoding="utf-8").strip() != qualification_fingerprint(suite_path):
            _print({"status": "NOT_ADMITTED", "reason": "project regression anchor hash is missing or mismatched"})
            return EXIT_GATE_FAILED
    tasks = load_qualification_tasks(suite_path)
    try:
        manifest = MemoryManifest.model_validate(load_yaml(candidate_bundle / "manifest.yaml"))
        preflight_checks = validate_candidate_bundle(candidate_bundle, tasks)
    except (OSError, ValueError) as exc:
        _print({"status": "NOT_ADMITTED", "reason": f"candidate validation error: {type(exc).__name__}: {exc}"})
        return EXIT_GATE_FAILED
    failed_preflight = [check.model_dump(mode="json") for check in preflight_checks if not check.passed]
    if failed_preflight:
        _print({"status": "NOT_ADMITTED", "reason": "candidate preflight failed", "checks": failed_preflight})
        return EXIT_GATE_FAILED
    with tempfile.TemporaryDirectory(prefix="ai-research-gate-") as snapshot_dir:
        decision, raw = await run_live_gate(
            tasks=tasks,
            incumbent_bundle=store.version_dir(incumbent_id),
            candidate_bundle=candidate_bundle,
            cwd=Path(snapshot_dir),
            source_failure_ids=manifest.source_failure_ids,
            concurrency=concurrency,
        )
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    report_dir = store.root / "reports" / f"{timestamp}-{candidate_id}-vs-{incumbent_id}"
    report_dir.mkdir(parents=True, exist_ok=False)
    atomic_write(report_dir / "decision.json", decision.model_dump_json(indent=2) + "\n")
    atomic_write(report_dir / "outcomes.json", json.dumps(raw, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    atomic_write(report_dir / "suite.sha256", qualification_fingerprint(suite_path) + "\n")
    for task_id, trials in raw["trials"].items():
        for trial_index, trial in enumerate(trials):
            for arm in ("incumbent", "candidate"):
                append_jsonl(
                    report_dir / "trace.jsonl",
                    {
                        "task_id": task_id,
                        "trial": trial_index,
                        "arm": arm,
                        "outcome": trial[arm],
                    },
                )
    store.record_gate_decision(candidate_id, decision.status, str(report_dir / "decision.json"))
    if decision.status == "PROMOTED":
        store.materialize_candidate_as_version(candidate_id)
        store.promote(candidate_id, reason="frozen qualification gate passed", gate_report=str(report_dir / "decision.json"))
    _print({"status": decision.status, "report": str(report_dir), "decision": decision.model_dump(mode="json")})
    return EXIT_OK if decision.status == "PROMOTED" else EXIT_GATE_FAILED


async def _evolve(project: Path, run_id: str, scope: str, concurrency: int) -> int:
    project = project.resolve()
    initialize_project(project)
    store = _scope_store(scope, project)
    _bootstrap_scope(store, scope)
    candidate_id, proposal = await propose_candidate(project_root=project, run_id=run_id, scope=scope, store=store)
    if candidate_id is None:
        _print({"status": "blocked", "reason": "failure classified as non-memory-patchable", "proposal": proposal})
        return EXIT_BLOCKED
    return await _run_gate(candidate_id, scope, project, concurrency)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="run_ai_research.sh", description="Evidence-gated Codex AI research harness")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("doctor")

    init = subparsers.add_parser("init")
    init.add_argument("--project", type=Path, required=True)

    run = subparsers.add_parser("run")
    run.add_argument("--project", type=Path, required=True)
    prompt_group = run.add_mutually_exclusive_group(required=True)
    prompt_group.add_argument("--prompt")
    prompt_group.add_argument("--task-file", type=Path)

    resume = subparsers.add_parser("resume")
    resume.add_argument("--project", type=Path, required=True)
    resume.add_argument("--run-id", required=True)

    status = subparsers.add_parser("status")
    status.add_argument("--project", type=Path, required=True)
    status.add_argument("--run-id")

    evolve = subparsers.add_parser("evolve")
    evolve.add_argument("--project", type=Path, required=True)
    evolve.add_argument("--run-id", required=True)
    evolve.add_argument("--scope", choices=["project", "global"], default="project")
    evolve.add_argument("--concurrency", type=int, default=3, choices=range(1, 4))

    gate = subparsers.add_parser("gate")
    gate.add_argument("--candidate", required=True)
    gate.add_argument("--scope", choices=["project", "global"], default="global")
    gate.add_argument("--project", type=Path)
    gate.add_argument("--concurrency", type=int, default=3, choices=range(1, 4))

    rollback = subparsers.add_parser("rollback")
    rollback.add_argument("--version", required=True)
    rollback.add_argument("--scope", choices=["project", "global"], default="global")
    rollback.add_argument("--project", type=Path)
    return parser


async def async_main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "doctor":
        return await _doctor()
    if args.command == "init":
        _print({"status": "complete", **initialize_project(args.project)})
        return EXIT_OK
    if args.command == "run":
        prompt = args.prompt if args.prompt is not None else args.task_file.read_text(encoding="utf-8")
        state, result = await HarnessRuntime().run(args.project, prompt)
        _print({"status": state.goals[0].status, "run_id": state.run_id, "thread_id": result.thread_id, "response": result.final_response})
        return state_exit_code(state)
    if args.command == "resume":
        state, result = await HarnessRuntime().resume(args.project, args.run_id)
        _print({"status": state.goals[0].status, "run_id": state.run_id, "thread_id": result.thread_id, "response": result.final_response})
        return state_exit_code(state)
    if args.command == "status":
        store = RunStore(args.project.resolve())
        if args.run_id:
            state = store.load_state(args.run_id)
            _print(state.model_dump(mode="json"))
            return state_exit_code(state)
        else:
            runs = sorted(path.name for path in store.runs.iterdir() if path.is_dir()) if store.runs.is_dir() else []
            _print({"project": str(args.project.resolve()), "runs": runs})
        return EXIT_OK
    if args.command == "evolve":
        return await _evolve(args.project, args.run_id, args.scope, args.concurrency)
    if args.command == "gate":
        return await _run_gate(args.candidate, args.scope, args.project, args.concurrency)
    if args.command == "rollback":
        store = _scope_store(args.scope, args.project)
        store.promote(args.version, reason="explicit user rollback", action="rollback")
        _print({"status": "complete", "champion": store.champion_id(), "scope": args.scope})
        return EXIT_OK
    raise AssertionError(args.command)


def main(argv: list[str] | None = None) -> int:
    try:
        return asyncio.run(async_main(argv))
    except (FileNotFoundError, ValueError, KeyError) as exc:
        _print({"status": "blocked", "error": f"{type(exc).__name__}: {exc}"})
        return EXIT_BLOCKED
    except Exception as exc:
        _print({"status": "environment_error", "error": f"{type(exc).__name__}: {exc}"})
        return EXIT_ENVIRONMENT


if __name__ == "__main__":
    sys.exit(main())

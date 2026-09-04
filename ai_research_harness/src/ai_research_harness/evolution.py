from __future__ import annotations

import json
import shutil
import tempfile
import tomllib
import uuid
from pathlib import Path
from typing import Any

import yaml

from .memory import load_bundle_cards
from .models import GoalStatus, MemoryManifest, utc_now
from .sdk_adapter import AgentBackend, OpenAICodexBackend
from .storage import MemoryStore, RunStore, atomic_write, load_yaml, write_yaml


EVOLUTION_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "required": ["classification", "component", "reason", "source_failure_ids", "patch"],
    "properties": {
        "classification": {
            "type": "string",
            "enum": ["patchable", "model", "tool", "data", "environment", "evaluation"],
        },
        "component": {"type": ["string", "null"], "enum": ["E", "W", "rho", "C", None]},
        "reason": {"type": "string"},
        "source_failure_ids": {"type": "array", "items": {"type": "string"}},
        "patch": {"type": "object"},
    },
}

ALLOWED_W_KEYS = {"done_requires", "blocked_requires", "model_may_only_propose_state"}
ALLOWED_RHO_KEYS = {"events", "max_cards", "max_tokens", "project_precedence", "conflict_policy", "selection_basis"}
ALLOWED_C_KEYS = {"goal_checks", "unknown_checker_policy", "priority"}
REGISTERED_CHECKERS = {"agent_response_present", "no_blocked_marker"}


def _load_evolver_instructions() -> str:
    path = Path("/home/h1028/.codex/agents/ai_research_memory_evolver.toml")
    return str(tomllib.loads(path.read_text(encoding="utf-8"))["developer_instructions"])


def _safe_skill_id(value: str) -> str:
    allowed = set("abcdefghijklmnopqrstuvwxyz0123456789-_")
    normalized = value.casefold()
    if not normalized or any(character not in allowed for character in normalized):
        raise ValueError("skill_id must contain lowercase letters, digits, hyphen, or underscore")
    return normalized


def _apply_declarative_patch(bundle: Path, component: str, patch: dict[str, Any], scope: str) -> None:
    manifest_path = bundle / "manifest.yaml"
    manifest_raw = load_yaml(manifest_path)
    if component == "E":
        required = {
            "skill_id",
            "domains",
            "triggers",
            "anti_triggers",
            "sources",
            "evidence_requirements",
            "body",
        }
        if set(patch) != required:
            raise ValueError(f"E patch keys must be exactly {sorted(required)}")
        skill_id = _safe_skill_id(str(patch["skill_id"]))
        body = str(patch["body"]).strip()
        if not body:
            raise ValueError("skill body cannot be empty")
        metadata = {
            "skill_id": skill_id,
            "version": 1,
            "status": "active",
            "scope": scope,
            "component": "E",
            "domains": list(patch["domains"]),
            "triggers": list(patch["triggers"]),
            "anti_triggers": list(patch["anti_triggers"]),
            "sources": list(patch["sources"]),
            "evidence_requirements": list(patch["evidence_requirements"]),
        }
        relative = f"cards/{skill_id}.md"
        target = bundle / relative
        if target.exists():
            raise ValueError(f"automatic evolution cannot overwrite an existing card: {skill_id}")
        target.parent.mkdir(parents=True, exist_ok=True)
        atomic_write(target, "---\n" + yaml.safe_dump(metadata, allow_unicode=True, sort_keys=False) + "---\n\n" + body + "\n")
        manifest_raw["cards"] = [*manifest_raw.get("cards", []), relative]
    else:
        updates = patch.get("updates")
        if set(patch) != {"updates"} or not isinstance(updates, dict):
            raise ValueError(f"{component} patch must contain only an updates mapping")
        if component == "W":
            relative = str(manifest_raw["working_memory_spec"])
            allowed = ALLOWED_W_KEYS
        elif component == "rho":
            relative = str(manifest_raw["invocation_policy"])
            allowed = ALLOWED_RHO_KEYS
        elif component == "C":
            relative = str(manifest_raw["checker_config"])
            allowed = ALLOWED_C_KEYS
        else:
            raise ValueError(f"unsupported component: {component}")
        if not set(updates).issubset(allowed):
            raise ValueError(f"disallowed {component} keys: {sorted(set(updates) - allowed)}")
        config = load_yaml(bundle / relative)
        config.update(updates)
        if component == "rho":
            if int(config.get("max_cards", 0)) not in range(1, 5):
                raise ValueError("rho max_cards must be 1..4")
            if int(config.get("max_tokens", 0)) not in range(1, 3501):
                raise ValueError("rho max_tokens must be 1..3500")
            if config.get("conflict_policy") != "block" or config.get("project_precedence") is not True:
                raise ValueError("rho cannot weaken conflict blocking or project isolation")
            required_events = {"turn_start", "phase_boundary", "checker_failure"}
            if not required_events.issubset(set(config.get("events", []))):
                raise ValueError("rho cannot remove mandatory selection events")
        if component == "C":
            checks = set(config.get("goal_checks", []))
            if not checks or not checks.issubset(REGISTERED_CHECKERS):
                raise ValueError("C may only select registered deterministic checkers")
            if config.get("unknown_checker_policy") != "block":
                raise ValueError("C cannot weaken unknown checker policy")
        write_yaml(bundle / relative, config)
    write_yaml(manifest_path, manifest_raw)


def build_candidate(
    *,
    store: MemoryStore,
    incumbent_id: str,
    proposal: dict[str, Any],
    scope: str,
) -> str:
    if proposal.get("classification") != "patchable" or proposal.get("component") not in {"E", "W", "rho", "C"}:
        raise ValueError(f"failure is not memory-patchable: {proposal.get('classification')}")
    candidate_id = f"candidate-{uuid.uuid4().hex[:12]}"
    temporary = Path(tempfile.mkdtemp(prefix=f".{candidate_id}.", dir=store.root))
    try:
        shutil.copytree(store.version_dir(incumbent_id), temporary / candidate_id)
        bundle = temporary / candidate_id
        (bundle / "ANCHOR.sha256").unlink(missing_ok=True)
        manifest = load_yaml(bundle / "manifest.yaml")
        manifest.update(
            {
                "version_id": candidate_id,
                "scope": scope,
                "parent_version": incumbent_id,
                "source_failure_ids": list(proposal.get("source_failure_ids", [])),
                "created_at": utc_now(),
            }
        )
        write_yaml(bundle / "manifest.yaml", manifest)
        _apply_declarative_patch(bundle, str(proposal["component"]), dict(proposal["patch"]), scope)
        store.install_bundle(bundle, candidate_id, candidate=True)
        load_bundle_cards(store.candidate_dir(candidate_id))
    finally:
        shutil.rmtree(temporary, ignore_errors=True)
    return candidate_id


async def propose_candidate(
    *,
    project_root: Path,
    run_id: str,
    scope: str,
    store: MemoryStore,
    backend: AgentBackend | None = None,
) -> tuple[str | None, dict[str, Any]]:
    run_store = RunStore(project_root)
    state = run_store.load_state(run_id)
    if not any(goal.status == GoalStatus.blocked for goal in state.goals):
        return None, {
            "classification": "evaluation",
            "component": None,
            "reason": "run has no recorded failure to attribute",
            "source_failure_ids": [],
            "patch": {},
        }
    trace_path = run_store.run_dir(run_id) / "trace.jsonl"
    trace_lines = trace_path.read_text(encoding="utf-8").splitlines()[-40:]
    incumbent_id = store.champion_id()
    manifest = load_yaml(store.version_dir(incumbent_id) / "manifest.yaml")
    prompt = (
        "根据以下失败状态和 trace 进行一次最小 E/W/rho/C 定位。只有确属记忆问题时返回 patchable。"
        "patch 只能符合宿主 schema，不得改代码、模型、权限、预算、gate 或冻结题。\n\n"
        f"scope={scope}\nworking_state={json.dumps(state.model_dump(mode='json'), ensure_ascii=False)}\n"
        f"incumbent_manifest={json.dumps(manifest, ensure_ascii=False)}\n"
        f"trace={json.dumps(trace_lines, ensure_ascii=False)}"
    )
    owns_backend = backend is None
    active_backend = backend or OpenAICodexBackend()
    try:
        if owns_backend:
            await active_backend.__aenter__()  # type: ignore[attr-defined]
        result = await active_backend.execute(
            prompt=prompt,
            cwd=project_root,
            developer_instructions=_load_evolver_instructions(),
            skill_paths=[],
            sandbox="read_only",
            output_schema=EVOLUTION_SCHEMA,
        )
    finally:
        if owns_backend:
            await active_backend.__aexit__(None, None, None)  # type: ignore[attr-defined]
    try:
        proposal = json.loads(result.final_response or "")
    except json.JSONDecodeError as exc:
        raise ValueError(f"memory evolver returned invalid JSON: {exc}") from exc
    classification = proposal.get("classification")
    if classification != "patchable":
        return None, proposal
    candidate_id = build_candidate(store=store, incumbent_id=incumbent_id, proposal=proposal, scope=scope)
    return candidate_id, proposal

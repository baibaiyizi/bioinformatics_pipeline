from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Iterable

import yaml
from pydantic import BaseModel

from .models import MemoryManifest, TraceRecord, WorkingMemory, utc_now


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value: Any) -> bytes:
    if isinstance(value, BaseModel):
        value = value.model_dump(mode="json")
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def model_sha256(value: BaseModel) -> str:
    return sha256_bytes(canonical_json(value))


def atomic_write(path: Path, data: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def append_jsonl(path: Path, value: BaseModel | dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = value.model_dump(mode="json") if isinstance(value, BaseModel) else value
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected YAML mapping: {path}")
    return value


def write_yaml(path: Path, value: BaseModel | dict[str, Any]) -> None:
    payload = value.model_dump(mode="json") if isinstance(value, BaseModel) else value
    atomic_write(path, yaml.safe_dump(payload, allow_unicode=True, sort_keys=False))


def bundle_hash(root: Path, *, exclude: Iterable[str] = ("ANCHOR.sha256",)) -> str:
    excluded = set(exclude)
    digest = hashlib.sha256()
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative in excluded:
            continue
        digest.update(relative.encode("utf-8") + b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


class RunStore:
    def __init__(self, project_root: Path) -> None:
        self.project_root = project_root.resolve()
        self.root = self.project_root / ".ai_research"
        self.runs = self.root / "runs"
        self.project_memory = self.root / "memory"
        self.anchors = self.root / "regression_anchors"

    def initialize(self) -> None:
        for directory in (self.runs, self.project_memory / "versions", self.project_memory / "candidates", self.anchors):
            directory.mkdir(parents=True, exist_ok=True)

    def run_dir(self, run_id: str) -> Path:
        return self.runs / run_id

    def create_run(self, state: WorkingMemory, prompt: str) -> Path:
        target = self.run_dir(state.run_id)
        if target.exists():
            raise FileExistsError(f"run already exists: {state.run_id}")
        target.mkdir(parents=True)
        atomic_write(target / "prompt.md", prompt)
        self.save_state(state)
        append_jsonl(target / "trace.jsonl", TraceRecord(run_id=state.run_id, event="run_created"))
        return target

    def save_state(self, state: WorkingMemory) -> None:
        state.touch()
        write_yaml(self.run_dir(state.run_id) / "working_memory.yaml", state)

    def load_state(self, run_id: str) -> WorkingMemory:
        return WorkingMemory.model_validate(load_yaml(self.run_dir(run_id) / "working_memory.yaml"))

    def trace(self, record: TraceRecord) -> None:
        append_jsonl(self.run_dir(record.run_id) / "trace.jsonl", record)

    def last_thread_id(self, run_id: str) -> str | None:
        trace_path = self.run_dir(run_id) / "trace.jsonl"
        if not trace_path.is_file():
            return None
        for line in reversed(trace_path.read_text(encoding="utf-8").splitlines()):
            try:
                value = json.loads(line)
            except json.JSONDecodeError:
                continue
            thread_id = value.get("thread_id")
            if isinstance(thread_id, str) and thread_id:
                return thread_id
        return None


class MemoryStore:
    def __init__(self, root: Path, scope: str) -> None:
        self.root = root.resolve()
        self.scope = scope
        self.versions = self.root / "versions"
        self.candidates = self.root / "candidates"
        self.pointer = self.root / "CHAMPION"
        self.ledger = self.root / "ledger.jsonl"

    def initialize(self) -> None:
        self.versions.mkdir(parents=True, exist_ok=True)
        self.candidates.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _validate_id(version_id: str) -> str:
        allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.")
        if not version_id or any(character not in allowed for character in version_id):
            raise ValueError(f"invalid version ID: {version_id!r}")
        return version_id

    def version_dir(self, version_id: str) -> Path:
        return self.versions / self._validate_id(version_id)

    def candidate_dir(self, candidate_id: str) -> Path:
        return self.candidates / self._validate_id(candidate_id)

    def install_bundle(self, source: Path, version_id: str, *, candidate: bool = False) -> Path:
        self.initialize()
        target = self.candidate_dir(version_id) if candidate else self.version_dir(version_id)
        if target.exists():
            self.verify_bundle(target)
            return target
        parent = target.parent
        temporary = Path(tempfile.mkdtemp(prefix=f".{version_id}.", dir=parent))
        try:
            for item in source.iterdir():
                destination = temporary / item.name
                if item.is_dir():
                    shutil.copytree(item, destination)
                else:
                    shutil.copy2(item, destination)
            manifest_path = temporary / "manifest.yaml"
            manifest = MemoryManifest.model_validate(load_yaml(manifest_path))
            if manifest.version_id != version_id or manifest.scope != self.scope:
                raise ValueError("bundle manifest identity does not match target store")
            atomic_write(temporary / "ANCHOR.sha256", bundle_hash(temporary) + "\n")
            os.replace(temporary, target)
        finally:
            if temporary.exists():
                shutil.rmtree(temporary)
        return target

    def verify_bundle(self, bundle: Path) -> str:
        anchor = bundle / "ANCHOR.sha256"
        if not anchor.is_file():
            raise ValueError(f"missing bundle anchor: {bundle}")
        expected = anchor.read_text(encoding="utf-8").strip()
        actual = bundle_hash(bundle)
        if actual != expected:
            raise ValueError(f"bundle hash mismatch: {bundle.name}")
        MemoryManifest.model_validate(load_yaml(bundle / "manifest.yaml"))
        return actual

    def champion_id(self) -> str:
        if not self.pointer.is_file():
            raise FileNotFoundError("memory champion is not initialized")
        version_id = self.pointer.read_text(encoding="utf-8").strip()
        self.verify_bundle(self.version_dir(version_id))
        return version_id

    def promote(self, version_id: str, *, reason: str, gate_report: str | None = None, action: str = "promote") -> None:
        bundle = self.version_dir(version_id)
        anchor = self.verify_bundle(bundle)
        previous = self.pointer.read_text(encoding="utf-8").strip() if self.pointer.exists() else None
        atomic_write(self.pointer, version_id + "\n")
        append_jsonl(
            self.ledger,
            {
                "timestamp": utc_now(),
                "action": action,
                "scope": self.scope,
                "from": previous,
                "to": version_id,
                "anchor": anchor,
                "reason": reason,
                "gate_report": gate_report,
            },
        )

    def record_gate_decision(self, candidate_id: str, status: str, gate_report: str) -> None:
        if status not in {"PROMOTED", "NOT_ADMITTED"}:
            raise ValueError(f"invalid gate decision: {status}")
        self._validate_id(candidate_id)
        append_jsonl(
            self.ledger,
            {
                "timestamp": utc_now(),
                "action": "gate_decision",
                "scope": self.scope,
                "candidate": candidate_id,
                "status": status,
                "champion_at_decision": self.champion_id(),
                "gate_report": gate_report,
            },
        )

    def bootstrap(self, base_bundle: Path) -> str:
        manifest = MemoryManifest.model_validate(load_yaml(base_bundle / "manifest.yaml"))
        self.install_bundle(base_bundle, manifest.version_id)
        if not self.pointer.exists():
            self.promote(manifest.version_id, reason="initial immutable baseline", action="bootstrap")
        return self.champion_id()

    def materialize_candidate_as_version(self, candidate_id: str) -> Path:
        source = self.candidate_dir(candidate_id)
        self.verify_bundle(source)
        target = self.version_dir(candidate_id)
        if target.exists():
            self.verify_bundle(target)
            return target
        temporary = Path(tempfile.mkdtemp(prefix=f".{candidate_id}.", dir=self.versions))
        shutil.rmtree(temporary)
        shutil.copytree(source, temporary)
        os.replace(temporary, target)
        return target

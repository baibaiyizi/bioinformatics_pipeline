from __future__ import annotations

import hashlib
import json
import logging
import os
import platform
import shlex
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from time import monotonic
from typing import Any, Iterator, Sequence

import pandas as pd

from .config import WorkflowConfig


LOGGER = logging.getLogger("scrna_template")


def configure_logging(log_file: Path | None = None, verbose: bool = False) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stdout)]
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        handle.write(value)
        temporary = Path(handle.name)
    temporary.replace(path)


def atomic_json(path: Path, value: Any) -> None:
    atomic_text(path, json.dumps(value, indent=2, ensure_ascii=False, default=str) + "\n")


def atomic_tsv(path: Path, frame: pd.DataFrame, index: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        frame.to_csv(handle, sep="\t", index=index)
        temporary = Path(handle.name)
    temporary.replace(path)


def sha256(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def digest_files(paths: Sequence[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted({item.resolve() for item in paths}):
        digest.update(str(path).encode())
        digest.update(sha256(path).encode())
    return digest.hexdigest()


def reference_signature(reference: Path) -> str:
    relative = (
        "reference.json",
        "fasta/genome.fa.fai",
        "genes/genes.gtf.gz",
        "star/genomeParameters.txt",
    )
    paths = [reference / item for item in relative]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"reference 契约文件缺失: {missing}")
    return digest_files(paths)


def project_code_signature(root: Path) -> str:
    paths: list[Path] = []
    for relative in ("pyproject.toml", "uv.lock", "run_scrna.sh"):
        path = root / relative
        if path.is_file():
            paths.append(path)
    for folder in (root / "src", root / "scripts", root / "config"):
        if folder.exists():
            paths.extend(path for path in folder.rglob("*") if path.is_file())
    return digest_files(paths)


def run_command(
    command: Sequence[str | Path],
    *,
    cwd: Path | None = None,
    log_file: Path | None = None,
    env: dict[str, str] | None = None,
) -> None:
    values = [str(value) for value in command]
    LOGGER.info("Command: %s", shlex.join(values))
    if log_file is None:
        subprocess.run(values, cwd=cwd, env=env, check=True)
        return
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with log_file.open("a", encoding="utf-8") as handle:
        handle.write(f"$ {shlex.join(values)}\n")
        handle.flush()
        subprocess.run(values, cwd=cwd, env=env, check=True, stdout=handle, stderr=subprocess.STDOUT, text=True)


RUN_DIRS = (
    "00_audit",
    "01_cellranger",
    "02_qc",
    "03_atlas",
    "04_annotation",
    "05_descriptive",
    "06_optional",
    "07_source_data",
    "08_figures",
    "09_reports",
    "logs",
    "tmp",
)


@dataclass
class RunContext:
    config: WorkflowConfig
    run_id: str
    run_dir: Path
    base_signature: str
    statuses: list[dict[str, Any]] = field(default_factory=list)

    @classmethod
    def create(cls, config: WorkflowConfig, run_id: str, input_paths: list[Path]) -> "RunContext":
        import re
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", run_id):
            raise ValueError("run-id 必须是安全目录名")
        run_dir = config.run_root / run_id
        if run_dir.exists():
            raise RuntimeError("run-id 已存在；使用新 run-id，禁止把旧文件当当前结果")
        signature_paths = [config.source, config.samples_path, config.markers_path, *input_paths]
        code_paths = list(Path(__file__).parent.glob("*.py"))
        base_payload = {"inputs_config_sha256": digest_files(signature_paths), "code_sha256": digest_files(code_paths)}
        signature = hashlib.sha256(json.dumps(base_payload, sort_keys=True).encode()).hexdigest()
        run_dir.mkdir(parents=True)
        for name in RUN_DIRS:
            (run_dir / name).mkdir()
        context = cls(config, run_id, run_dir, signature)
        context.write_manifest(base_payload)
        atomic_tsv(context.path("00_audit", "input_sha256.tsv"), pd.DataFrame([{"path": str(path), "sha256": sha256(path)} for path in signature_paths]))
        return context

    def path(self, *parts: str) -> Path:
        return self.run_dir.joinpath(*parts)

    def stage_signature(self, stage: str) -> str:
        return hashlib.sha256(f"{self.base_signature}:{stage}".encode()).hexdigest()

    def state_path(self, stage: str) -> Path:
        return self.path("00_audit", "stages", f"{stage}.json")

    def reusable(self, stage: str, outputs: Sequence[Path]) -> bool:
        state_path = self.state_path(stage)
        if not state_path.is_file():
            return False
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("signature") != self.stage_signature(stage):
            raise RuntimeError(f"{stage} 已有输出与当前 input/config/reference/code signature 不匹配，请使用新 run-id")
        return state.get("status") == "completed" and all(path.exists() for path in outputs)

    @contextmanager
    def stage(self, stage: str, outputs: Sequence[Path] = ()) -> Iterator[bool]:
        if self.reusable(stage, outputs):
            self.record_status(stage, "skipped", "signature 与输出契约完全匹配，安全 resume")
            yield False
            return
        started = datetime.now().isoformat(timespec="seconds")
        start = monotonic()
        atomic_json(self.state_path(stage), {"stage": stage, "status": "running", "signature": self.stage_signature(stage), "started": started})
        try:
            yield True
            missing = [str(path) for path in outputs if not path.exists()]
            if missing:
                raise RuntimeError(f"{stage} 输出契约缺失: {missing}")
        except Exception as exc:
            elapsed = monotonic() - start
            atomic_json(self.state_path(stage), {"stage": stage, "status": "failed", "signature": self.stage_signature(stage), "message": str(exc), "elapsed_seconds": elapsed})
            self.record_status(stage, "failed", str(exc), elapsed)
            raise
        else:
            missing = [str(path) for path in outputs if not path.exists()]
            if missing:
                raise RuntimeError(f"{stage} 输出契约缺失: {missing}")
            elapsed = monotonic() - start
            atomic_json(self.state_path(stage), {"stage": stage, "status": "completed", "signature": self.stage_signature(stage), "elapsed_seconds": elapsed, "outputs": [str(path) for path in outputs]})
            self.record_status(stage, "completed", elapsed_seconds=elapsed)

    def record_status(self, stage: str, status: str, message: str = "", elapsed_seconds: float | None = None) -> None:
        self.statuses.append({"stage": stage, "status": status, "message": message, "elapsed_seconds": elapsed_seconds, "timestamp": datetime.now().isoformat(timespec="seconds")})
        atomic_tsv(self.path("00_audit", "stage_status.tsv"), pd.DataFrame(self.statuses))

    def record_no_go(self, stage: str, message: str) -> None:
        atomic_json(self.state_path(stage), {"stage": stage, "status": "NO-GO", "signature": self.stage_signature(stage), "message": message})
        self.record_status(stage, "NO-GO", message)

    def write_manifest(self, signatures: dict[str, str]) -> None:
        packages: dict[str, str] = {}
        for package in ("anndata", "scanpy", "numpy", "pandas", "scipy", "scrublet"):
            try:
                packages[package] = version(package)
            except PackageNotFoundError:
                packages[package] = "not-installed"
        atomic_json(self.path("00_audit", "run_manifest.json"), {
            "run_id": self.run_id,
            "created": datetime.now().isoformat(timespec="seconds"),
            "project_root": str(self.config.root),
            "python": sys.version,
            "platform": platform.platform(),
            "pid": os.getpid(),
            "packages": packages,
            "signatures": signatures,
            "base_signature": self.base_signature,
        })

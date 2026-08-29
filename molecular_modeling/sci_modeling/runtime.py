from __future__ import annotations

import importlib.util
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any


CUDA_GROMACS_MARKER = Path("/home/h1028/workspace/model/environment_manifests/verification/gromacs-cuda-smoke.json")


def has_cuda_support(gmx_version_text: str) -> bool:
    return re.search(r"^GPU support:\s+CUDA\s*$", gmx_version_text, flags=re.MULTILINE) is not None


def mdrun_used_gpu(mdrun_log_text: str) -> bool:
    return re.search(r"\b[1-9][0-9]* GPUs? selected for this run\.", mdrun_log_text) is not None


def dependencies() -> dict[str, Any]:
    gmx = shutil.which("gmx")
    gmx_text = ""
    if gmx:
        result = subprocess.run([gmx, "--version"], text=True, capture_output=True, check=False)
        gmx_text = result.stdout + result.stderr
    gpu_nodes = all(Path(path).exists() for path in ["/dev/nvidia0", "/dev/nvidiactl", "/dev/nvidia-uvm"])
    gpu_kernel_visible = any(Path("/proc/driver/nvidia/gpus").glob("*/information"))
    smi = subprocess.run(["nvidia-smi", "-L"], text=True, capture_output=True, check=False) if shutil.which("nvidia-smi") else None
    gromacs_cuda_build = has_cuda_support(gmx_text)
    gpu_platform_ready = bool(gpu_nodes and smi and smi.returncode == 0 and gromacs_cuda_build)
    cuda_marker_valid = False
    if gpu_platform_ready and CUDA_GROMACS_MARKER.is_file():
        try:
            marker = json.loads(CUDA_GROMACS_MARKER.read_text(encoding="utf-8"))
            cuda_marker_valid = marker.get("gmx") == str(Path(gmx).resolve()) and marker.get("gmx_version") == gmx_text.strip()
        except (OSError, json.JSONDecodeError):
            cuda_marker_valid = False
    return {
        "vina": shutil.which("vina"),
        "meeko_receptor": shutil.which("mk_prepare_receptor.py"),
        "meeko_ligand": shutil.which("mk_prepare_ligand.py"),
        "meeko_export": shutil.which("mk_export.py"),
        "reduce": shutil.which("reduce"),
        "posebusters": importlib.util.find_spec("posebusters") is not None,
        "rdkit": importlib.util.find_spec("rdkit") is not None,
        "antechamber": shutil.which("antechamber"),
        "parmchk2": shutil.which("parmchk2"),
        "tleap": shutil.which("tleap"),
        "parmed": importlib.util.find_spec("parmed") is not None,
        "mdanalysis": importlib.util.find_spec("MDAnalysis") is not None,
        "gmx": gmx,
        "gmx_version": gmx_text.strip(),
        "gromacs_cuda_build": gromacs_cuda_build,
        "gpu_kernel_visible": gpu_kernel_visible,
        "gpu_nodes": gpu_nodes,
        "nvidia_smi": bool(smi and smi.returncode == 0),
        "gpu_platform_ready": gpu_platform_ready,
        "gromacs_cuda_smoke": cuda_marker_valid,
        "gpu_ready": bool(gpu_platform_ready and cuda_marker_valid),
    }


class CommandError(RuntimeError):
    pass


class Runner:
    def run(self, args: list[str], cwd: Path, *, stdin: str | None = None, log: Path | None = None) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(args, cwd=cwd, input=stdin, text=True, capture_output=True, check=False)
        if log:
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(json.dumps({"command": args, "cwd": str(cwd), "returncode": result.returncode}, ensure_ascii=False) + "\n" + result.stdout + "\n" + result.stderr, encoding="utf-8")
        if result.returncode != 0:
            raise CommandError(f"命令失败 ({result.returncode}): {' '.join(args)}")
        return result

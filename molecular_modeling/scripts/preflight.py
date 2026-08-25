#!/usr/bin/env python3
"""Scientific preflight for docking and MD; does not execute simulations."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


class ContractError(ValueError):
    pass


def resolve(base: Path, raw: str) -> Path:
    path = Path(raw).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    required = {"receptor_pdb", "ligand_table", "output_dir", "docking", "md"}
    missing = sorted(required - set(config))
    if missing:
        raise ContractError(f"配置缺少字段: {', '.join(missing)}")
    docking = config["docking"]
    for field in ["engine", "box_center", "box_size", "seeds", "redocking"]:
        if field not in docking:
            raise ContractError(f"docking 缺少字段: {field}")
    if docking["engine"] not in {"vina", "diffdock"}:
        raise ContractError("docking.engine 仅支持 vina 或 diffdock")
    if len(docking["box_center"]) != 3 or len(docking["box_size"]) != 3 or any(float(x) <= 0 for x in docking["box_size"]):
        raise ContractError("对接盒中心/尺寸必须是三个数，尺寸必须为正")
    if len(set(docking["seeds"])) < 3:
        raise ContractError("对接至少需要三个不同随机种子")
    md = config["md"]
    if bool(md.get("enabled", False)):
        required_md = {"replicas", "seeds", "forcefield", "water_model", "ligand_parameterization", "temperature_k", "pressure_bar", "production_ns"}
        missing_md = sorted(required_md - set(md))
        if missing_md:
            raise ContractError(f"md 缺少字段: {', '.join(missing_md)}")
        if int(md["replicas"]) < 3 or len(set(md["seeds"])) < int(md["replicas"]):
            raise ContractError("MD 至少需要三个独立重复及对应的不同随机种子")
        if float(md["production_ns"]) <= 0:
            raise ContractError("production_ns 必须为正；长度是否充分仍需按体系验证收敛")
    return config


def audit_pdb(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ContractError(f"受体 PDB 不存在: {path}")
    atoms = 0
    residues: set[tuple[str, str, str]] = set()
    chains: set[str] = set()
    altlocs: set[str] = set()
    model_count = 0
    missing_residue_remark = False
    resolution = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atoms += 1
            chains.add(line[21:22].strip() or "_")
            residues.add((line[21:22], line[22:26].strip(), line[17:20].strip()))
            if line[16:17].strip():
                altlocs.add(line[16:17].strip())
        elif line.startswith("MODEL"):
            model_count += 1
        elif line.startswith("REMARK 465"):
            missing_residue_remark = True
        elif line.startswith("REMARK   2 RESOLUTION."):
            match = re.search(r"RESOLUTION\.\s+([0-9.]+)\s+ANGSTROMS", line)
            resolution = float(match.group(1)) if match else None
    if atoms == 0:
        raise ContractError("受体 PDB 没有 ATOM/HETATM 记录")
    return {"atoms": atoms, "residues": len(residues), "chains": sorted(chains), "altlocs": sorted(altlocs), "models": max(model_count, 1), "missing_residue_remark": missing_residue_remark, "resolution_angstrom": resolution}


def audit_ligands(path: Path, base: Path, redocking: bool) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if not path.is_file():
        raise ContractError(f"配体表不存在: {path}")
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = {"ligand_id", "input_structure", "protonation_state", "formal_charge", "charge_method", "role"}
    missing = sorted(required - set(frame))
    if missing:
        raise ContractError(f"配体表缺少列: {', '.join(missing)}")
    if frame["ligand_id"].eq("").any() or frame["ligand_id"].duplicated().any():
        raise ContractError("ligand_id 必须非空且唯一")
    allowed_roles = {"candidate", "native", "positive_control", "decoy"}
    if not set(frame["role"]).issubset(allowed_roles):
        raise ContractError("role 仅允许 candidate/native/positive_control/decoy")
    if redocking and (frame["role"] == "native").sum() != 1:
        raise ContractError("redocking=true 时必须恰有一个 native ligand")
    if not (frame["role"] == "positive_control").any() or not (frame["role"] == "decoy").any():
        raise ContractError("对接验证必须同时包含 positive_control 和 decoy")
    records = []
    resolved_paths = []
    for row in frame.to_dict("records"):
        structure = resolve(base, row["input_structure"])
        if not structure.is_file():
            raise ContractError(f"配体结构不存在: {structure}")
        if not row["protonation_state"] or not row["formal_charge"] or not row["charge_method"]:
            raise ContractError(f"配体 {row['ligand_id']} 缺少质子化/电荷声明")
        resolved_paths.append(str(structure))
        records.append({"ligand_id": row["ligand_id"], "path": str(structure), "sha256": sha256(structure), "role": row["role"], "protonation_state": row["protonation_state"], "formal_charge": row["formal_charge"], "charge_method": row["charge_method"]})
    frame["resolved_structure"] = resolved_paths
    return frame, records


def dependency_audit(engine: str) -> dict[str, bool]:
    return {
        "vina_executable": shutil.which("vina") is not None,
        "vina_python": importlib.util.find_spec("vina") is not None,
        "diffdock_entrypoint": shutil.which("diffdock") is not None,
        "openbabel": shutil.which("obabel") is not None,
        "posebusters": shutil.which("bust") is not None or importlib.util.find_spec("posebusters") is not None,
        "openmm": importlib.util.find_spec("openmm") is not None,
        "mdanalysis": importlib.util.find_spec("MDAnalysis") is not None,
        "selected_engine_available": (shutil.which("vina") is not None or importlib.util.find_spec("vina") is not None) if engine == "vina" else shutil.which("diffdock") is not None,
    }


def stage_status(config: dict[str, Any], deps: dict[str, bool], pdb_audit: dict[str, Any]) -> list[dict[str, str]]:
    receptor_warnings = pdb_audit["missing_residue_remark"] or bool(pdb_audit["altlocs"]) or pdb_audit["models"] > 1
    docking_ready = deps["selected_engine_available"] and deps["openbabel"] and (deps["posebusters"] or not bool(config["docking"].get("require_posebusters", True)))
    md_enabled = bool(config["md"].get("enabled", False))
    md_ready = md_enabled and deps["openmm"] and deps["mdanalysis"] and bool(config["md"].get("ligand_parameterization"))
    return [
        {"stage": "structure_audit", "status": "review_required" if receptor_warnings else "ready", "reason": "missing residues/altlocs/multiple models require explicit resolution" if receptor_warnings else "basic PDB checks passed"},
        {"stage": "docking", "status": "ready" if docking_ready else "blocked", "reason": "engine, Open Babel and pose validation available" if docking_ready else "missing selected engine, Open Babel or PoseBusters"},
        {"stage": "md", "status": "ready" if md_ready else ("not_requested" if not md_enabled else "blocked"), "reason": "OpenMM and MDAnalysis available" if md_ready else "missing OpenMM/MDAnalysis or MD not requested"},
        {"stage": "mechanistic_claim", "status": "blocked", "reason": "requires convergence across replicas plus independent experimental evidence"},
    ]


def write_plan(path: Path, config: dict[str, Any], statuses: list[dict[str, str]]) -> None:
    lines = [
        "# 分子建模执行计划", "",
        "## 固定顺序", "",
        "1. 人工解决结构审计中的 missing residues、altloc、链、辅因子、结晶水与质子化状态。",
        "2. 用 native ligand 做 redocking，并与阳性/诱饵对照一起运行多个 seed。",
        "3. 对 pose 执行几何/碰撞/口袋合理性检查；docking score 只用于同协议内排序。",
        "4. 仅对通过 pose 验证的体系做配体参数化、溶剂化、最小化、NVT/NPT 平衡。",
        f"5. 运行 {config['md'].get('replicas', 0)} 个独立 MD 重复，分别保存 checkpoint、thermodynamic log 和轨迹。",
        "6. 去除平衡段后比较 replica 间 RMSD/RMSF、接触占有率、关键距离和收敛；必要时再设计自由能分析。",
        "7. 不以 docking score 声称亲和力，不以单条轨迹声称稳定性或机制。", "",
        "## 当前状态", "",
    ]
    lines.extend(f"- {row['stage']}: `{row['status']}` — {row['reason']}" for row in statuses)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(config_path: Path) -> Path:
    config_path = config_path.resolve()
    config = load_config(config_path)
    receptor = resolve(config_path.parent, config["receptor_pdb"])
    ligand_table = resolve(config_path.parent, config["ligand_table"])
    out = resolve(config_path.parent, config["output_dir"])
    pdb = audit_pdb(receptor)
    ligands, ligand_audit = audit_ligands(ligand_table, ligand_table.parent, bool(config["docking"]["redocking"]))
    deps = dependency_audit(config["docking"]["engine"])
    statuses = stage_status(config, deps, pdb)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(statuses).to_csv(out / "stage_status.tsv", sep="\t", index=False)
    ligands.to_csv(out / "ligands_validated.csv", index=False)
    write_plan(out / "execution_plan.md", config, statuses)
    audit = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "receptor": str(receptor), "receptor_sha256": sha256(receptor), "pdb_audit": pdb,
        "ligands": ligand_audit, "dependencies": deps, "stages": statuses,
        "docking_seeds": list(config["docking"]["seeds"]),
        "md_seeds": list(config["md"].get("seeds", [])),
        "claim_boundaries": ["docking score is not binding affinity", "a single short MD trajectory is not stability or mechanism evidence"],
        "execution_performed": False,
    }
    (out / "preflight.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "effective_config.yml").write_text(yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    try:
        print(run(args.config))
    except (ContractError, OSError, yaml.YAMLError, pd.errors.ParserError) as exc:
        print(f"[CONTRACT_ERROR] {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

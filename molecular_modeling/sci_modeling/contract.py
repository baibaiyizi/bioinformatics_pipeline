from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


class ContractError(ValueError):
    """A required scientific control or executable input is missing."""


SAFE_ID = re.compile(r"^[A-Za-z0-9_.:-]{1,128}$")
ROLES = {"native", "positive_control", "decoy", "candidate"}


def resolve(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _id(value: Any, field: str) -> str:
    text = str(value).strip()
    if not SAFE_ID.fullmatch(text):
        raise ContractError(f"{field} 含空值、控制字符或未允许字符")
    return text


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    required = {"project_id", "receptor", "pocket", "ligand_table", "output_dir", "docking", "md", "analysis"}
    missing = sorted(required - set(config))
    if missing:
        raise ContractError(f"配置缺少字段: {', '.join(missing)}")
    config["project_id"] = _id(config["project_id"], "project_id")

    receptor = config["receptor"]
    receptor_required = {"input_pdb", "source", "protein_chains", "covalent_ligand", "protonation_method"}
    if not isinstance(receptor, dict) or not receptor_required <= set(receptor):
        raise ContractError("receptor 必须包含 input_pdb/source/protein_chains/covalent_ligand/protonation_method")
    if not isinstance(receptor["source"], dict) or not {"database", "id"} <= set(receptor["source"]):
        raise ContractError("receptor.source 必须包含 database 和 id")
    receptor["source"]["database"] = _id(receptor["source"]["database"], "receptor.source.database")
    receptor["source"]["id"] = _id(receptor["source"]["id"], "receptor.source.id")
    if not receptor["protein_chains"]:
        raise ContractError("必须显式提供 protein_chains")
    if bool(receptor["covalent_ligand"]):
        raise ContractError("共价配体体系禁止自动参数化")
    if receptor["protonation_method"] != "reduce_hbond_network":
        raise ContractError("当前自动流程的 receptor.protonation_method 固定为 reduce_hbond_network")

    pocket = config["pocket"]
    pocket_required = {"basis", "native_ligand_resname", "native_reference_sdf", "center", "size"}
    if not isinstance(pocket, dict) or not pocket_required <= set(pocket):
        raise ContractError("pocket 缺少口袋依据、native ligand 或 box")
    if pocket["basis"] not in {"co_crystal_ligand", "experimentally_validated_site"}:
        raise ContractError("口袋必须基于共晶配体或独立实验验证位点")
    if len(pocket["center"]) != 3 or len(pocket["size"]) != 3 or any(float(x) <= 0 for x in pocket["size"]):
        raise ContractError("pocket.center/size 必须各为三个数且尺寸为正")
    _id(pocket["native_ligand_resname"], "native_ligand_resname")

    docking = config["docking"]
    docking_required = {"engine", "seeds", "exhaustiveness", "num_modes", "require_posebusters", "redocking"}
    if not isinstance(docking, dict) or not docking_required <= set(docking):
        raise ContractError("docking 配置不完整")
    if docking["engine"] != "vina":
        raise ContractError("对接引擎固定为 vina")
    if len(set(int(x) for x in docking["seeds"])) < 3:
        raise ContractError("对接至少需要三个不同 seed")
    redocking = docking["redocking"]
    if not isinstance(redocking, dict) or not {"enabled", "rmsd_threshold_angstrom", "minimum_passing_seeds"} <= set(redocking):
        raise ContractError("redocking 配置不完整")
    if not bool(redocking["enabled"]):
        raise ContractError("必须启用 native ligand redocking")
    if float(redocking["rmsd_threshold_angstrom"]) > 2.0:
        raise ContractError("redocking RMSD 阈值不得高于 2 Å")
    if int(redocking["minimum_passing_seeds"]) < 2:
        raise ContractError("至少两个 seed 的最佳有效 pose 必须通过 redocking")
    if not bool(docking["require_posebusters"]):
        raise ContractError("PoseBusters 检查不得关闭")

    md = config["md"]
    if not isinstance(md, dict) or not {"engine", "ligand_id", "parameterization", "replicas", "seeds", "smoke", "production"} <= set(md):
        raise ContractError("md 配置不完整")
    if md["engine"] != "gromacs":
        raise ContractError("MD 生产引擎固定为 GROMACS")
    parameterization = md["parameterization"]
    fixed = {"protein_forcefield": "ff14SB", "ligand_forcefield": "GAFF2", "charge_method": "AM1-BCC", "water_model": "TIP3P"}
    for field, expected in fixed.items():
        if str(parameterization.get(field, "")).lower() != expected.lower():
            raise ContractError(f"md.parameterization.{field} 必须为 {expected}")
    if int(md["replicas"]) < 3 or len(set(int(x) for x in md["seeds"])) < int(md["replicas"]):
        raise ContractError("MD 至少需要三个独立 replica 及不同 seed")
    production = md["production"]
    if not isinstance(production, dict) or "enabled" not in production:
        raise ContractError("md.production 必须显式声明 enabled")
    if bool(production["enabled"]):
        if "duration_ns" not in production or float(production["duration_ns"]) <= 0:
            raise ContractError("生产 MD 必须显式提供正的 duration_ns；没有通用默认值")
        if not production.get("approved_pose"):
            raise ContractError("生产 MD 必须显式提供 approved_pose")
    elif "duration_ns" in production and production["duration_ns"] not in {None, ""}:
        raise ContractError("未启用 production 时不要填写通用时长")
    return config


def load_ligands(config_path: Path, config: dict[str, Any]) -> pd.DataFrame:
    table_path = resolve(config_path.parent, str(config["ligand_table"]))
    if not table_path.is_file():
        raise ContractError(f"ligand_table 不存在: {table_path}")
    frame = pd.read_csv(table_path, sep="\t" if table_path.suffix.lower() == ".tsv" else ",", dtype=str).fillna("")
    required = {"ligand_id", "input_structure", "role", "protonation_state", "formal_charge", "charge_method"}
    missing = sorted(required - set(frame))
    if missing:
        raise ContractError(f"ligand_table 缺少列: {', '.join(missing)}")
    if frame["ligand_id"].eq("").any() or frame["ligand_id"].duplicated().any():
        raise ContractError("ligand_id 必须非空且唯一")
    if not set(frame["role"]) <= ROLES:
        raise ContractError("role 仅允许 native/positive_control/decoy/candidate")
    if (frame["role"] == "native").sum() != 1:
        raise ContractError("必须恰有一个 native ligand")
    if not (frame["role"] == "positive_control").any() or not (frame["role"] == "decoy").any():
        raise ContractError("必须同时包含阳性对照和 decoy")
    if str(config["md"]["ligand_id"]) not in set(frame["ligand_id"]):
        raise ContractError("md.ligand_id 不在 ligand_table")
    ambiguous = {"", "unknown", "unclear", "unspecified", "not_set", "na", "n/a"}
    for row in frame.to_dict("records"):
        _id(row["ligand_id"], "ligand_id")
        if row["protonation_state"].strip().lower() in ambiguous:
            raise ContractError(f"配体 {row['ligand_id']} 的质子化状态不明确")
        if row["charge_method"].strip().upper() != "AM1-BCC":
            raise ContractError(f"配体 {row['ligand_id']} 的 charge_method 必须为 AM1-BCC")
        try:
            int(row["formal_charge"])
        except ValueError as exc:
            raise ContractError(f"配体 {row['ligand_id']} formal_charge 必须为整数") from exc
        structure = resolve(table_path.parent, row["input_structure"])
        if not structure.is_file():
            raise ContractError(f"配体结构不存在: {structure}")
    frame["resolved_structure"] = [str(resolve(table_path.parent, value)) for value in frame["input_structure"]]
    frame["sha256"] = [sha256(Path(value)) for value in frame["resolved_structure"]]
    return frame

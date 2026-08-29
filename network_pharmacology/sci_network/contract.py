from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml


class ContractError(ValueError):
    """Raised when the scientific or execution contract is incomplete."""


SAFE_ID = re.compile(r"^[A-Za-z0-9_.:-]{1,128}$")
SPECIES = {9606: "homo_sapiens", 10090: "mus_musculus"}


def _identifier(value: Any, field: str) -> str:
    text = str(value).strip()
    if not SAFE_ID.fullmatch(text):
        raise ContractError(f"{field} 含空值、控制字符或未允许字符")
    return text


def resolve(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    required = {"project_id", "species_taxid", "disease_ids", "compounds", "target_seeds", "output_dir"}
    missing = sorted(required - set(config))
    if missing:
        raise ContractError(f"配置缺少字段: {', '.join(missing)}")

    config["project_id"] = _identifier(config["project_id"], "project_id")
    taxid = int(config["species_taxid"])
    if taxid not in SPECIES:
        raise ContractError("species_taxid 仅允许 9606 或 10090")
    config["species_taxid"] = taxid

    if not isinstance(config["disease_ids"], list) or not isinstance(config["compounds"], list) or not isinstance(config["target_seeds"], list):
        raise ContractError("disease_ids、compounds、target_seeds 必须是列表")
    config["disease_ids"] = [_identifier(x, "disease_ids") for x in config["disease_ids"]]

    compound_ids: set[tuple[str, str]] = set()
    for compound in config["compounds"]:
        if not isinstance(compound, dict) or not {"id", "id_type"} <= set(compound):
            raise ContractError("每个 compound 必须包含 id 和 id_type")
        compound["id_type"] = str(compound["id_type"])
        if compound["id_type"] not in {"chembl", "pubchem_cid", "inchikey"}:
            raise ContractError("compound.id_type 仅允许 chembl/pubchem_cid/inchikey")
        compound["id"] = _identifier(compound["id"], "compound.id")
        key = (compound["id_type"], compound["id"])
        if key in compound_ids:
            raise ContractError(f"重复 compound: {key}")
        compound_ids.add(key)

    seed_ids: set[str] = set()
    expected_prefix = "ENSG" if taxid == 9606 else "ENSMUSG"
    for seed in config["target_seeds"]:
        needed = {"ensembl_gene_id", "evidence_source", "evidence_id"}
        if not isinstance(seed, dict) or not needed <= set(seed):
            raise ContractError("每个 target_seed 必须包含 ensembl_gene_id/evidence_source/evidence_id")
        for field in needed:
            seed[field] = _identifier(seed[field], f"target_seed.{field}")
        if not seed["ensembl_gene_id"].startswith(expected_prefix):
            raise ContractError("target_seed 的 Ensembl ID 与 species_taxid 不兼容")
        if seed["ensembl_gene_id"] in seed_ids:
            raise ContractError("target_seed 不得重复")
        seed_ids.add(seed["ensembl_gene_id"])

    if not config["disease_ids"] and not config["compounds"] and not config["target_seeds"]:
        raise ContractError("至少提供一个 disease、compound 或 target seed")

    retrieval = config.setdefault("retrieval", {})
    retrieval.setdefault("page_size", 100)
    retrieval.setdefault("max_targets", 50)
    retrieval.setdefault("max_structure_checks", 10)
    retrieval.setdefault("timeout_seconds", 30)
    retrieval.setdefault("string_required_score", 400)
    for field in ("page_size", "max_targets", "max_structure_checks"):
        if int(retrieval[field]) <= 0:
            raise ContractError(f"retrieval.{field} 必须为正整数")
    if not 0 <= int(retrieval["string_required_score"]) <= 1000:
        raise ContractError("STRING required_score 必须位于 0..1000")

    analysis = config.setdefault("analysis", {})
    analysis.setdefault("string_thresholds", [400, 700, 900])
    analysis.setdefault("reactome_background_uniprot", [])
    thresholds = sorted({int(x) for x in analysis["string_thresholds"]})
    if not thresholds or any(x < 0 or x > 1000 for x in thresholds):
        raise ContractError("analysis.string_thresholds 必须是 0..1000 内的非空列表")
    analysis["string_thresholds"] = thresholds
    background = [_identifier(x, "reactome_background_uniprot") for x in analysis["reactome_background_uniprot"]]
    if len(background) > 80:
        raise ContractError("Reactome 背景逐项查询最多 80 个；更大背景请改用官方 bulk 数据")
    analysis["reactome_background_uniprot"] = list(dict.fromkeys(background))
    return config

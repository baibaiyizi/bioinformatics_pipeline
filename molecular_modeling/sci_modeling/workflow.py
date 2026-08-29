from __future__ import annotations

import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .contract import ContractError, load_config, load_ligands, resolve, sha256
from .runtime import CUDA_GROMACS_MARKER, CommandError, Runner, dependencies, mdrun_used_gpu


STANDARD_RESIDUES = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
WATER = {"HOH", "WAT", "TIP3", "SOL"}
METALS = {"LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG"}


def audit_receptor(path: Path, native_resname: str, chains: list[str]) -> dict[str, Any]:
    if not path.is_file():
        raise ContractError(f"receptor PDB 不存在: {path}")
    atoms, native_atoms = 0, 0
    altlocs: set[str] = set()
    protein_residues: set[str] = set()
    metals: set[str] = set()
    pdb_chains: set[str] = set()
    links: list[str] = []
    missing_residue_remark = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atoms += 1
            record, resname, chain = line[:6].strip(), line[17:20].strip().upper(), line[21:22].strip() or "_"
            element = (line[76:78].strip() or line[12:14].strip()).upper()
            pdb_chains.add(chain)
            if line[16:17].strip():
                altlocs.add(line[16:17].strip())
            if record == "ATOM":
                protein_residues.add(resname)
            if record == "HETATM" and resname == native_resname.upper():
                native_atoms += 1
            if record == "HETATM" and (element in METALS or resname in METALS):
                metals.add(resname or element)
        elif line.startswith("LINK"):
            links.append(line.rstrip())
        elif line.startswith("REMARK 465"):
            missing_residue_remark = True
    if atoms == 0:
        raise ContractError("receptor PDB 没有 ATOM/HETATM")
    missing_chains = sorted(set(chains) - pdb_chains)
    nonstandard = sorted(protein_residues - STANDARD_RESIDUES)
    blockers = []
    if metals:
        blockers.append(f"metal:{','.join(sorted(metals))}")
    if nonstandard:
        blockers.append(f"nonstandard_residue:{','.join(nonstandard)}")
    if links:
        blockers.append("LINK_records_require_covalent_or_cofactor_review")
    if native_atoms == 0:
        blockers.append("native_ligand_missing_from_receptor")
    if missing_chains:
        blockers.append(f"protein_chains_missing:{','.join(missing_chains)}")
    return {"atoms": atoms, "native_ligand_atoms": native_atoms, "pdb_chains": sorted(pdb_chains), "altlocs": sorted(altlocs), "missing_residue_remark": missing_residue_remark, "metals": sorted(metals), "nonstandard_residues": nonstandard, "link_records": links, "automatic_parameterization_blockers": blockers}


def _context(config_path: Path) -> tuple[dict[str, Any], pd.DataFrame, Path, Path, dict[str, Any]]:
    config_path = config_path.resolve()
    config = load_config(config_path)
    ligands = load_ligands(config_path, config)
    receptor = resolve(config_path.parent, str(config["receptor"]["input_pdb"]))
    native_reference = resolve(config_path.parent, str(config["pocket"]["native_reference_sdf"]))
    if not native_reference.is_file():
        raise ContractError(f"native_reference_sdf 不存在: {native_reference}")
    out = resolve(config_path.parent, str(config["output_dir"]))
    audit = audit_receptor(receptor, str(config["pocket"]["native_ligand_resname"]), [str(x) for x in config["receptor"]["protein_chains"]])
    return config, ligands, receptor, out, {"receptor": audit, "native_reference": native_reference}


def _statuses(config: dict[str, Any], receptor_audit: dict[str, Any], deps: dict[str, Any]) -> list[dict[str, str]]:
    structure_blocked = bool(receptor_audit["automatic_parameterization_blockers"])
    review = receptor_audit["missing_residue_remark"] or bool(receptor_audit["altlocs"])
    review_decision = str(config["receptor"].get("structure_review_decision", "")).strip()
    review_resolved = not review or bool(review_decision)
    scientific_ready = not structure_blocked and review_resolved
    docking_deps = all(deps[key] for key in ["vina", "meeko_receptor", "meeko_ligand", "meeko_export", "reduce", "posebusters", "rdkit"])
    amber_deps = all(deps[key] for key in ["antechamber", "parmchk2", "tleap", "parmed", "gmx"])
    return [
        {"stage": "preflight", "status": "blocked" if structure_blocked else ("review_required" if not review_resolved else "ready"), "reason": "|".join(receptor_audit["automatic_parameterization_blockers"]) if structure_blocked else ("missing residues or altlocs require an explicit structure_review_decision" if not review_resolved else (f"structure review accepted: {review_decision}" if review else "scientific inputs passed"))},
        {"stage": "dock", "status": "ready" if docking_deps and scientific_ready else "blocked", "reason": "Vina, Meeko, PoseBusters and RDKit available; structure review resolved" if docking_deps and scientific_ready else ("unresolved structure audit" if not scientific_ready else "missing Vina/Meeko/PoseBusters/RDKit")},
        {"stage": "md-prepare", "status": "ready" if amber_deps and scientific_ready else "blocked", "reason": "AmberTools, ParmEd and GROMACS available; structure review resolved" if amber_deps and scientific_ready else ("unresolved structure audit" if not scientific_ready else "missing AmberTools/ParmEd/GROMACS")},
        {"stage": "md-smoke", "status": "ready" if deps["gmx"] and deps["mdanalysis"] and scientific_ready else "blocked", "reason": "CPU smoke is permitted after validated handoff" if deps["gmx"] and deps["mdanalysis"] and scientific_ready else ("unresolved structure audit" if not scientific_ready else "missing GROMACS or MDAnalysis")},
        {"stage": "md-production", "status": "ready" if deps["gpu_ready"] and bool(config["md"]["production"]["enabled"]) and scientific_ready else "resource-limited", "reason": "real CUDA device, CUDA GROMACS, explicit duration and approved pose passed" if deps["gpu_ready"] and scientific_ready else ("unresolved structure audit" if not scientific_ready else "real CUDA execution unavailable; production is forbidden")},
        {"stage": "mechanistic_claim", "status": "blocked", "reason": "requires converged independent replicas and external experimental evidence"},
    ]


def preflight(config_path: Path) -> Path:
    config, ligands, receptor, out, audits = _context(config_path)
    deps = dependencies()
    statuses = _statuses(config, audits["receptor"], deps)
    out.mkdir(parents=True, exist_ok=True)
    ligands.to_csv(out / "ligands_validated.tsv", sep="\t", index=False)
    pd.DataFrame(statuses).to_csv(out / "stage_status.tsv", sep="\t", index=False)
    payload = {"created_at": datetime.now(timezone.utc).isoformat(), "project_id": config["project_id"], "receptor": str(receptor), "receptor_sha256": sha256(receptor), "native_reference_sdf": str(audits["native_reference"]), "native_reference_sha256": sha256(audits["native_reference"]), "structure_audit": audits["receptor"], "dependencies": deps, "stages": statuses, "claim_boundaries": ["docking score is not binding affinity", "a short trajectory is not stability or mechanism evidence"]}
    (out / "preflight.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "effective_config.yml").write_text(yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return out


def _require_stage_ready(out: Path, stage: str) -> None:
    status_path = out / "stage_status.tsv"
    if not status_path.is_file():
        raise ContractError("必须先运行 preflight")
    table = pd.read_csv(status_path, sep="\t")
    row = table[table["stage"] == stage]
    if row.empty or row.iloc[0]["status"] != "ready":
        reason = "missing status" if row.empty else row.iloc[0]["reason"]
        raise ContractError(f"阶段 {stage} 未通过硬门: {reason}")


def _posebusters_flags(pose_sdf: Path, native_sdf: Path, receptor: Path, output_csv: Path, role: str) -> list[bool]:
    from posebusters import PoseBusters

    if role == "native":
        result = PoseBusters(config="redock").bust(mol_pred=pose_sdf, mol_true=native_sdf, mol_cond=receptor)
    else:
        result = PoseBusters(config="dock").bust(mol_pred=pose_sdf, mol_cond=receptor)
    result.to_csv(output_csv)
    bool_columns = [column for column in result.columns if pd.api.types.is_bool_dtype(result[column])]
    if not bool_columns:
        raise ContractError("PoseBusters 没有返回可判定的布尔检查列")
    return [bool(value) for value in result[bool_columns].all(axis=1).tolist()]


def _pose_rmsds(pose_sdf: Path, native_sdf: Path) -> list[float]:
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    poses = [mol for mol in Chem.SDMolSupplier(str(pose_sdf), removeHs=False) if mol is not None]
    references = [mol for mol in Chem.SDMolSupplier(str(native_sdf), removeHs=False) if mol is not None]
    if not poses or not references:
        raise ContractError("Meeko 输出或 native reference 无法由 RDKit 读取")
    reference = Chem.RemoveHs(references[0])
    values = []
    for pose in poses:
        pose = Chem.RemoveHs(pose)
        if pose.GetNumHeavyAtoms() != reference.GetNumHeavyAtoms():
            values.append(math.nan)
            continue
        try:
            values.append(float(rdMolAlign.GetBestRMS(pose, reference)))
        except RuntimeError:
            values.append(math.nan)
    if not any(math.isfinite(value) for value in values):
        raise ContractError("没有可做重原子对称校正 RMSD 的 pose")
    return values


def _write_selected_pose(source_sdf: Path, pose_index: int, output_sdf: Path) -> None:
    from rdkit import Chem

    poses = [mol for mol in Chem.SDMolSupplier(str(source_sdf), removeHs=False) if mol is not None]
    if pose_index < 0 or pose_index >= len(poses):
        raise ContractError(f"批准的 pose_index 越界: {pose_index}/{len(poses)}")
    writer = Chem.SDWriter(str(output_sdf))
    writer.write(poses[pose_index])
    writer.close()
    if not output_sdf.is_file() or output_sdf.stat().st_size == 0:
        raise ContractError("未能写出批准的单 pose SDF")


def dock(config_path: Path, runner: Runner | None = None) -> Path:
    config, ligands, receptor, out, audits = _context(config_path)
    if not (out / "stage_status.tsv").is_file():
        preflight(config_path)
    _require_stage_ready(out, "dock")
    runner = runner or Runner()
    work = out / "docking"
    work.mkdir(parents=True, exist_ok=True)
    protein = work / "receptor_protein.pdb"
    receptor_h = work / "receptorH.pdb"
    _write_protein_only(receptor, protein, [str(x) for x in config["receptor"]["protein_chains"]])
    reduced = runner.run(["reduce", "-BUILD", str(protein)], work, log=work / "logs/reduce.log")
    receptor_h.write_text(reduced.stdout, encoding="utf-8")
    if not any(line.startswith("ATOM  ") for line in receptor_h.read_text(encoding="utf-8", errors="replace").splitlines()):
        raise ContractError("Reduce 未生成加氢受体")
    prefix = work / "receptor"
    center = [str(float(x)) for x in config["pocket"]["center"]]
    size = [str(float(x)) for x in config["pocket"]["size"]]
    runner.run(["mk_prepare_receptor.py", "-i", str(receptor_h), "-o", str(prefix), "-p", "--box_center", *center, "--box_size", *size], work, log=work / "logs/prepare_receptor.log")
    receptor_pdbqt = prefix.with_suffix(".pdbqt")
    if not receptor_pdbqt.is_file():
        raise ContractError("Meeko 未生成 receptor.pdbqt")
    rows: list[dict[str, Any]] = []
    for ligand in ligands.to_dict("records"):
        ligand_dir = work / str(ligand["ligand_id"])
        ligand_dir.mkdir(parents=True, exist_ok=True)
        ligand_pdbqt = ligand_dir / "ligand.pdbqt"
        runner.run(["mk_prepare_ligand.py", "-i", ligand["resolved_structure"], "-o", str(ligand_pdbqt)], ligand_dir, log=ligand_dir / "prepare.log")
        for seed in [int(x) for x in config["docking"]["seeds"]]:
            seed_dir = ligand_dir / f"seed_{seed}"
            seed_dir.mkdir(parents=True, exist_ok=True)
            pose_pdbqt = seed_dir / "poses.pdbqt"
            runner.run(["vina", "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt), "--center_x", center[0], "--center_y", center[1], "--center_z", center[2], "--size_x", size[0], "--size_y", size[1], "--size_z", size[2], "--exhaustiveness", str(int(config["docking"]["exhaustiveness"])), "--num_modes", str(int(config["docking"]["num_modes"])), "--seed", str(seed), "--out", str(pose_pdbqt)], seed_dir, log=seed_dir / "vina.log")
            pose_sdf = seed_dir / "poses.sdf"
            runner.run(["mk_export.py", str(pose_pdbqt), "-s", str(pose_sdf)], seed_dir, log=seed_dir / "meeko_export.log")
            flags = _posebusters_flags(pose_sdf, audits["native_reference"], protein, seed_dir / "posebusters.csv", str(ligand["role"]))
            rmsds = _pose_rmsds(pose_sdf, audits["native_reference"]) if ligand["role"] == "native" else [math.nan] * len(flags)
            if len(flags) != len(rmsds):
                raise ContractError(f"PoseBusters 与 RDKit pose 数不一致: {len(flags)} != {len(rmsds)}")
            valid_indices = [index for index, valid in enumerate(flags) if valid]
            if ligand["role"] == "native":
                valid_indices = [index for index in valid_indices if math.isfinite(rmsds[index])]
                best_index = min(valid_indices, key=lambda index: rmsds[index]) if valid_indices else None
            else:
                best_index = valid_indices[0] if valid_indices else None
            valid = best_index is not None
            rmsd = rmsds[best_index] if best_index is not None else math.nan
            rows.append({"ligand_id": ligand["ligand_id"], "role": ligand["role"], "seed": seed, "pose_sdf": str(pose_sdf), "best_valid_pose_index": best_index, "posebusters_valid": valid, "symmetry_corrected_heavy_atom_rmsd_angstrom": rmsd, "redocking_pass": bool(ligand["role"] == "native" and valid and rmsd <= float(config["docking"]["redocking"]["rmsd_threshold_angstrom"]))})
    summary = pd.DataFrame(rows)
    summary.to_csv(work / "docking_summary.tsv", sep="\t", index=False)
    native = summary[summary["role"] == "native"]
    passing = native[native["redocking_pass"]]
    required = int(config["docking"]["redocking"]["minimum_passing_seeds"])
    gate = {"passing_native_seeds": int(len(passing)), "required_passing_seeds": required, "passed": bool(len(passing) >= required), "rmsd_threshold_angstrom": float(config["docking"]["redocking"]["rmsd_threshold_angstrom"]), "posebusters_required": True}
    (work / "redocking_gate.json").write_text(json.dumps(gate, ensure_ascii=False, indent=2), encoding="utf-8")
    if not gate["passed"]:
        raise ContractError(f"redocking 硬门未通过: {len(passing)}/{required} seeds")
    md_ligand = str(config["md"]["ligand_id"])
    approved = summary[(summary["ligand_id"] == md_ligand) & summary["posebusters_valid"]].copy()
    if approved.empty:
        raise ContractError("md.ligand_id 没有通过 PoseBusters 的 pose")
    if approved.iloc[0]["role"] == "native":
        approved = approved[approved["redocking_pass"]].sort_values("symmetry_corrected_heavy_atom_rmsd_angstrom")
        if approved.empty:
            raise ContractError("MD native ligand 没有通过 redocking 的单 pose")
    else:
        approved = approved.sort_values(["seed", "best_valid_pose_index"])
    selected = approved.iloc[0]
    approved_sdf = work / f"approved_{md_ligand}.sdf"
    _write_selected_pose(Path(str(selected["pose_sdf"])), int(selected["best_valid_pose_index"]), approved_sdf)
    handoff = {"approved": True, "ligand_id": md_ligand, "pose_sdf": str(approved_sdf), "source_pose_sdf": str(selected["pose_sdf"]), "source_pose_index": int(selected["best_valid_pose_index"]), "seed": int(selected["seed"]), "symmetry_corrected_heavy_atom_rmsd_angstrom": None if pd.isna(selected["symmetry_corrected_heavy_atom_rmsd_angstrom"]) else float(selected["symmetry_corrected_heavy_atom_rmsd_angstrom"]), "redocking_gate": gate, "approved_at": datetime.now(timezone.utc).isoformat()}
    (work / "md_handoff.json").write_text(json.dumps(handoff, ensure_ascii=False, indent=2), encoding="utf-8")
    return work


def _write_protein_only(receptor: Path, output: Path, chains: list[str]) -> None:
    selected = []
    chain_set = set(chains)
    for line in receptor.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM  ") and (line[21:22].strip() or "_") in chain_set:
            selected.append(line)
    if not selected:
        raise ContractError("protein_chains 过滤后没有 ATOM")
    output.write_text("\n".join(selected + ["TER", "END"]) + "\n", encoding="utf-8")


def _write_amber_protein(reduced_receptor: Path, output: Path, chains: list[str]) -> dict[str, str]:
    chain_set = set(chains)
    lines = [line for line in reduced_receptor.read_text(encoding="utf-8", errors="replace").splitlines() if line.startswith("ATOM  ") and (line[21:22].strip() or "_") in chain_set]
    if not lines:
        raise ContractError("Reduce 加氢受体没有可用于 AmberTools 的蛋白 ATOM")
    atoms_by_residue: dict[tuple[str, str, str], set[str]] = {}
    for line in lines:
        if line[17:20].strip().upper() == "HIS":
            key = (line[21:22].strip() or "_", line[22:26].strip(), line[26:27].strip())
            atoms_by_residue.setdefault(key, set()).add(line[12:16].strip().upper())
    assignments: dict[tuple[str, str, str], str] = {}
    for key, atom_names in atoms_by_residue.items():
        hd1, he2 = "HD1" in atom_names, "HE2" in atom_names
        if hd1 and he2:
            assignments[key] = "HIP"
        elif hd1:
            assignments[key] = "HID"
        elif he2:
            assignments[key] = "HIE"
        else:
            raise ContractError(f"Reduce 未明确 HIS 质子化状态: chain={key[0]} resid={key[1]}{key[2]}")
    prepared = []
    for line in lines:
        if line[17:20].strip().upper() == "HIS":
            key = (line[21:22].strip() or "_", line[22:26].strip(), line[26:27].strip())
            line = f"{line[:17]}{assignments[key]:>3}{line[20:]}"
        prepared.append(line)
    output.write_text("\n".join(prepared + ["TER", "END"]) + "\n", encoding="utf-8")
    return {f"{chain}:{resid}{icode}": state for (chain, resid, icode), state in sorted(assignments.items())}


def _write_mdp(path: Path, entries: dict[str, Any]) -> None:
    path.write_text("\n".join(f"{key} = {value}" for key, value in entries.items()) + "\n", encoding="utf-8")


def md_prepare(config_path: Path, runner: Runner | None = None) -> Path:
    config, ligands, receptor, out, audits = _context(config_path)
    _require_stage_ready(out, "md-prepare")
    if audits["receptor"]["automatic_parameterization_blockers"]:
        raise ContractError("结构审计存在自动参数化 blocker")
    handoff_path = out / "docking/md_handoff.json"
    if not handoff_path.is_file():
        raise ContractError("必须先通过 dock 并生成 md_handoff.json")
    handoff = json.loads(handoff_path.read_text(encoding="utf-8"))
    if not handoff.get("approved"):
        raise ContractError("docking pose 未批准")
    runner = runner or Runner()
    work = out / "md_system"
    work.mkdir(parents=True, exist_ok=True)
    protein = work / "protein.pdb"
    reduced_receptor = out / "docking/receptorH.pdb"
    if not reduced_receptor.is_file():
        raise ContractError("缺少通过 docking 阶段生成的 Reduce 加氢受体")
    histidine_states = _write_amber_protein(reduced_receptor, protein, [str(x) for x in config["receptor"]["protein_chains"]])
    ligand_row = ligands[ligands["ligand_id"] == str(config["md"]["ligand_id"])].iloc[0]
    ligand_sdf = Path(str(handoff["pose_sdf"]))
    charge = str(int(ligand_row["formal_charge"]))
    runner.run(["antechamber", "-i", str(ligand_sdf), "-fi", "sdf", "-o", "ligand.mol2", "-fo", "mol2", "-at", "gaff2", "-c", "bcc", "-nc", charge, "-s", "2", "-rn", "LIG"], work, log=work / "logs/antechamber.log")
    runner.run(["parmchk2", "-i", "ligand.mol2", "-f", "mol2", "-o", "ligand.frcmod", "-s", "gaff2"], work, log=work / "logs/parmchk2.log")
    leap = "\n".join(["source leaprc.protein.ff14SB", "source leaprc.gaff2", "source leaprc.water.tip3p", "LIG = loadmol2 ligand.mol2", "loadamberparams ligand.frcmod", "PROT = loadpdb protein.pdb", "COM = combine {PROT LIG}", "check COM", "solvatebox COM TIP3PBOX 10.0", "addionsrand COM Na+ 0", "addionsrand COM Cl- 0", "saveamberparm COM system.prmtop system.inpcrd", "savepdb COM system_amber.pdb", "quit", ""])
    (work / "tleap.in").write_text(leap, encoding="utf-8")
    runner.run(["tleap", "-f", "tleap.in"], work, log=work / "logs/tleap.log")
    parmed_code = "import parmed as p; s=p.load_file('system.prmtop','system.inpcrd'); s.save('system.top',overwrite=True); s.save('system.gro',overwrite=True)"
    runner.run(["python", "-c", parmed_code], work, log=work / "logs/parmed.log")
    required = [work / "system.prmtop", work / "system.inpcrd", work / "system.top", work / "system.gro"]
    if not all(path.is_file() and path.stat().st_size > 0 for path in required):
        raise ContractError("AmberTools/ParmEd 未生成完整 GROMACS 体系")
    manifest = {"created_at": datetime.now(timezone.utc).isoformat(), "protein_forcefield": "ff14SB", "ligand_forcefield": "GAFF2", "charge_method": "AM1-BCC", "water_model": "TIP3P", "receptor_protonation_method": config["receptor"]["protonation_method"], "histidine_states": histidine_states, "ligand_id": config["md"]["ligand_id"], "approved_pose": str(ligand_sdf), "files": {path.name: sha256(path) for path in required}}
    (work / "parameterization_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return work


def _grompp(runner: Runner, gmx: str, cwd: Path, mdp: Path, structure: Path, topology: Path, output: str, log: Path) -> None:
    runner.run([gmx, "grompp", "-f", str(mdp), "-c", str(structure), "-p", str(topology), "-o", output], cwd, log=log)


def md_smoke(config_path: Path, runner: Runner | None = None) -> Path:
    config, _, _, out, audits = _context(config_path)
    _require_stage_ready(out, "md-smoke")
    system = out / "md_system"
    if not (system / "parameterization_manifest.json").is_file():
        raise ContractError("必须先运行 md-prepare")
    runner = runner or Runner()
    deps = dependencies()
    gmx = str(deps["gmx"])
    gpu_args = ["-nb", "gpu", "-pme", "gpu"] if deps["gpu_platform_ready"] else []
    work = out / "md_smoke"
    work.mkdir(parents=True, exist_ok=True)
    steps = int(config["md"]["smoke"].get("steps", 500))
    if steps <= 0 or steps > 50000:
        raise ContractError("md.smoke.steps 必须位于 1..50000，仅用于短能力验收")
    em_mdp = work / "em.mdp"
    _write_mdp(em_mdp, {"integrator": "steep", "nsteps": 500, "emtol": 1000, "cutoff-scheme": "Verlet", "coulombtype": "PME", "rcoulomb": 1.0, "rvdw": 1.0, "constraints": "h-bonds"})
    records = []
    for replica, seed in enumerate([int(x) for x in config["md"]["seeds"][: int(config["md"]["replicas"])]], start=1):
        rep = work / f"replica_{replica}"
        rep.mkdir(parents=True, exist_ok=True)
        _grompp(runner, gmx, rep, em_mdp, system / "system.gro", system / "system.top", "em.tpr", rep / "grompp_em.log")
        runner.run([gmx, "mdrun", "-deffnm", "em", "-ntmpi", "1", "-ntomp", "4", *gpu_args], rep, log=rep / "mdrun_em.log")
        nvt_mdp = rep / "nvt.mdp"
        _write_mdp(nvt_mdp, {"integrator": "md", "dt": 0.002, "nsteps": steps, "continuation": "no", "constraint-algorithm": "lincs", "constraints": "h-bonds", "cutoff-scheme": "Verlet", "nstlist": 10, "rlist": 1.0, "coulombtype": "PME", "rcoulomb": 1.0, "vdwtype": "cut-off", "rvdw": 1.0, "tcoupl": "V-rescale", "tc-grps": "System", "tau_t": 0.1, "ref_t": 300, "pcoupl": "no", "gen_vel": "yes", "gen_temp": 300, "gen_seed": seed, "nstxout-compressed": 10, "nstenergy": 10, "nstlog": 10})
        _grompp(runner, gmx, rep, nvt_mdp, rep / "em.gro", system / "system.top", "nvt.tpr", rep / "grompp_nvt.log")
        first_steps = max(1, steps // 2)
        runner.run([gmx, "mdrun", "-deffnm", "nvt", "-nsteps", str(first_steps), "-cpt", "0.001", "-ntmpi", "1", "-ntomp", "4", *gpu_args], rep, log=rep / "mdrun_nvt_part1.log")
        checkpoint = rep / "nvt.cpt"
        if not checkpoint.is_file():
            raise ContractError("GROMACS smoke 未生成 checkpoint")
        runner.run([gmx, "mdrun", "-deffnm", "nvt", "-cpi", "nvt.cpt", "-append", "-ntmpi", "1", "-ntomp", "4", *gpu_args], rep, log=rep / "mdrun_nvt_resume.log")
        npt_mdp = rep / "npt.mdp"
        _write_mdp(npt_mdp, {"integrator": "md", "dt": 0.002, "nsteps": steps, "continuation": "yes", "constraint-algorithm": "lincs", "constraints": "h-bonds", "cutoff-scheme": "Verlet", "nstlist": 10, "rlist": 1.0, "coulombtype": "PME", "rcoulomb": 1.0, "vdwtype": "cut-off", "rvdw": 1.0, "tcoupl": "V-rescale", "tc-grps": "System", "tau_t": 0.1, "ref_t": 300, "pcoupl": "C-rescale", "pcoupltype": "isotropic", "tau_p": 2.0, "ref_p": 1.0, "compressibility": "4.5e-5", "gen_vel": "no", "nstxout-compressed": 10, "nstenergy": 10, "nstlog": 10})
        _grompp(runner, gmx, rep, npt_mdp, rep / "nvt.gro", system / "system.top", "npt.tpr", rep / "grompp_npt.log")
        runner.run([gmx, "mdrun", "-deffnm", "npt", "-ntmpi", "1", "-ntomp", "4", *gpu_args], rep, log=rep / "mdrun_npt.log")
        log_text = (rep / "mdrun_npt.log").read_text(encoding="utf-8", errors="replace")
        bad = bool(re.search(r"\bNaN\b|segmentation fault|blowing up", log_text, re.IGNORECASE))
        records.append({"replica": replica, "seed": seed, "steps": steps, "checkpoint_resume": True, "trajectory": str(rep / "npt.xtc"), "final_structure": str(rep / "npt.gro"), "no_nan_or_blowup": not bad, "gpu_used": mdrun_used_gpu(log_text)})
    result = pd.DataFrame(records)
    result.to_csv(work / "smoke_summary.tsv", sep="\t", index=False)
    if len(result) < 3 or not result["no_nan_or_blowup"].all():
        raise ContractError("三重复 smoke 或无 NaN/爆炸验收失败")
    if deps["gpu_platform_ready"]:
        if not result["gpu_used"].all():
            raise ContractError("GPU 平台可见，但三条 GROMACS smoke 日志未证明 CUDA 实际执行")
        CUDA_GROMACS_MARKER.parent.mkdir(parents=True, exist_ok=True)
        CUDA_GROMACS_MARKER.write_text(json.dumps({"validated_at": datetime.now(timezone.utc).isoformat(), "gmx": str(Path(gmx).resolve()), "gmx_version": deps["gmx_version"], "replicas": int(len(result))}, ensure_ascii=False, indent=2), encoding="utf-8")
    return work


def md_production(config_path: Path, runner: Runner | None = None) -> Path:
    config, _, _, out, _ = _context(config_path)
    production = config["md"]["production"]
    if not bool(production["enabled"]):
        raise ContractError("production 未显式启用")
    deps = dependencies()
    if not deps["gpu_ready"]:
        raise ContractError("真实 CUDA 设备/GROMACS CUDA 任务未通过，禁止生产 MD")
    approved = resolve(config_path.parent, str(production["approved_pose"]))
    if not approved.is_file():
        raise ContractError("approved_pose 不存在")
    handoff_path = out / "docking/md_handoff.json"
    if not handoff_path.is_file():
        raise ContractError("生产 MD 缺少 docking md_handoff.json")
    handoff = json.loads(handoff_path.read_text(encoding="utf-8"))
    handoff_pose = Path(str(handoff.get("pose_sdf", ""))).resolve()
    if not handoff.get("approved") or approved.resolve() != handoff_pose:
        raise ContractError("production.approved_pose 必须与通过 redocking/PoseBusters 的 handoff pose 一致")
    system = out / "md_system"
    if not (system / "parameterization_manifest.json").is_file():
        raise ContractError("必须先运行 md-prepare")
    runner = runner or Runner()
    work = out / "md_production"
    work.mkdir(parents=True, exist_ok=True)
    duration_ns = float(production["duration_ns"])
    nsteps = int(round(duration_ns * 500000))
    gmx = str(deps["gmx"])
    for replica, seed in enumerate([int(x) for x in config["md"]["seeds"][: int(config["md"]["replicas"])]], start=1):
        rep = work / f"replica_{replica}"
        rep.mkdir(parents=True, exist_ok=True)
        mdp = rep / "production.mdp"
        _write_mdp(mdp, {"integrator": "md", "dt": 0.002, "nsteps": nsteps, "continuation": "yes", "constraints": "h-bonds", "cutoff-scheme": "Verlet", "coulombtype": "PME", "rcoulomb": 1.0, "rvdw": 1.0, "tcoupl": "V-rescale", "tc-grps": "System", "tau_t": 0.1, "ref_t": 300, "pcoupl": "C-rescale", "pcoupltype": "isotropic", "tau_p": 2.0, "ref_p": 1.0, "compressibility": "4.5e-5", "gen_vel": "no", "ld_seed": seed, "nstxout-compressed": 5000, "nstenergy": 1000, "nstlog": 1000})
        start = out / "md_smoke" / f"replica_{replica}" / "npt.gro"
        _grompp(runner, gmx, rep, mdp, start, system / "system.top", "production.tpr", rep / "grompp.log")
        command = [gmx, "mdrun", "-deffnm", "production", "-nb", "gpu", "-pme", "gpu", "-ntmpi", "1", "-ntomp", "8"]
        if (rep / "production.cpt").is_file():
            command += ["-cpi", "production.cpt", "-append"]
        runner.run(command, rep, log=rep / "mdrun.log")
    return work


def analyze(config_path: Path, runner: Runner | None = None) -> Path:
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
    from MDAnalysis.analysis.rms import rmsd as calculate_rmsd
    from MDAnalysis.lib.distances import distance_array

    config, _, _, out, _ = _context(config_path)
    source = out / ("md_production" if bool(config["md"]["production"]["enabled"]) else "md_smoke")
    if not source.is_dir():
        raise ContractError("没有可分析的 MD 轨迹")
    work = out / "analysis"
    work.mkdir(parents=True, exist_ok=True)
    runner = runner or Runner()
    gmx = str(dependencies()["gmx"] or "")
    if not gmx:
        raise ContractError("分析温压/能量 QC 需要 GROMACS")
    summaries = []
    rmsd_rows, rmsf_rows, rg_rows, contact_rows, distance_rows, hbond_rows, thermo_rows = [], [], [], [], [], [], []
    ligand_sel = str(config["analysis"].get("ligand_selection", "resname LIG"))
    analysis_topology = out / "md_system/system.prmtop"
    if not analysis_topology.is_file():
        raise ContractError("缺少可供 MDAnalysis 读取且保留键连接的 Amber prmtop")
    for replica in range(1, int(config["md"]["replicas"]) + 1):
        rep = source / f"replica_{replica}"
        prefix = "production" if source.name == "md_production" else "npt"
        trajectory = rep / f"{prefix}.xtc"
        if not trajectory.is_file():
            raise ContractError(f"replica {replica} 缺少 trajectory")
        energy_xvg = work / f"replica_{replica}_thermodynamics.xvg"
        runner.run([gmx, "energy", "-f", str(rep / f"{prefix}.edr"), "-o", str(energy_xvg)], rep, stdin="Temperature\nPressure\nPotential\nKinetic-En.\nTotal-Energy\n0\n", log=work / f"replica_{replica}_gmx_energy.log")
        legends: dict[int, str] = {}
        for line in energy_xvg.read_text(encoding="utf-8", errors="replace").splitlines():
            legend = re.match(r'@ s(\d+) legend "([^"]+)"', line)
            if legend:
                legends[int(legend.group(1)) + 1] = legend.group(2)
            elif line and not line.startswith(("#", "@")):
                values = [float(value) for value in line.split()]
                for column, value in enumerate(values[1:], start=1):
                    thermo_rows.append({"replica": replica, "time_ps": values[0], "term": legends.get(column, f"column_{column}"), "value": value})
        # MDAnalysis 2.10 does not yet parse GROMACS 2026 TPX v138. The Amber
        # prmtop is the source topology used by ParmEd and preserves the same
        # atom order and bonds as the GROMACS trajectory.
        universe = mda.Universe(str(analysis_topology), str(trajectory))
        protein = universe.select_atoms("protein and backbone")
        ligand = universe.select_atoms(ligand_sel)
        if not protein.n_atoms or not ligand.n_atoms:
            raise ContractError(f"replica {replica} 的蛋白或配体 selection 为空")
        align.AlignTraj(universe, universe, select="protein and backbone", in_memory=True).run()
        universe.trajectory[0]
        protein_reference = protein.positions.copy()
        ligand_reference = ligand.positions.copy()
        hbond_atoms = f"protein or ({ligand_sel})"
        hbond_selector = HydrogenBondAnalysis(universe=universe)
        hydrogen_selection = hbond_selector.guess_hydrogens(hbond_atoms)
        acceptor_selection = hbond_selector.guess_acceptors(hbond_atoms)
        if not hydrogen_selection or not acceptor_selection:
            raise ContractError(f"replica {replica} 无法从拓扑识别蛋白/配体氢键原子")
        hbonds = HydrogenBondAnalysis(
            universe=universe,
            hydrogens_sel=hydrogen_selection,
            acceptors_sel=acceptor_selection,
            between=[["protein", ligand_sel]],
            update_selections=False,
        ).run().results.hbonds
        total_frames = len(universe.trajectory)
        if len(hbonds):
            hbond_frame = pd.DataFrame(hbonds, columns=["frame", "donor_index", "hydrogen_index", "acceptor_index", "distance_angstrom", "angle_degree"])
            for keys, group in hbond_frame.groupby(["donor_index", "hydrogen_index", "acceptor_index"]):
                donor, hydrogen, acceptor = [int(value) for value in keys]
                hbond_rows.append({"replica": replica, "donor_index": donor, "hydrogen_index": hydrogen, "acceptor_index": acceptor, "occupancy": float(group["frame"].nunique() / total_frames), "mean_distance_angstrom": float(group["distance_angstrom"].mean()), "mean_angle_degree": float(group["angle_degree"].mean())})
        ligand_contacts = []
        radii = []
        positions = []
        for frame, ts in enumerate(universe.trajectory):
            time_ps = float(ts.time)
            protein_value = float(calculate_rmsd(protein.positions, protein_reference, center=False, superposition=False))
            ligand_value = float(calculate_rmsd(ligand.positions, ligand_reference, center=False, superposition=False))
            rmsd_rows.extend([{"replica": replica, "time_ps": time_ps, "selection": "protein_backbone", "rmsd_angstrom": protein_value}, {"replica": replica, "time_ps": time_ps, "selection": "ligand_after_protein_fit", "rmsd_angstrom": ligand_value}])
            radii.append(float(protein.radius_of_gyration()))
            distances = distance_array(protein.positions, ligand.positions, box=universe.dimensions)
            ligand_contacts.append(int((distances <= float(config["analysis"].get("contact_cutoff_angstrom", 4.5))).any(axis=1).sum()))
            positions.append(protein.positions.copy())
            for item in config["analysis"].get("key_distances", []):
                a, b = universe.select_atoms(item["selection_a"]), universe.select_atoms(item["selection_b"])
                if len(a) != 1 or len(b) != 1:
                    raise ContractError(f"关键距离 selection 必须各命中一个原子: {item['name']}")
                value = float(distance_array(a.positions, b.positions, box=universe.dimensions)[0, 0])
                distance_rows.append({"replica": replica, "time_ps": time_ps, "distance_name": item["name"], "distance_angstrom": value})
            rg_rows.append({"replica": replica, "time_ps": time_ps, "protein_rg_angstrom": radii[-1]})
            contact_rows.append({"replica": replica, "time_ps": time_ps, "protein_atoms_contacting_ligand": ligand_contacts[-1]})
        coords = np.stack(positions)
        mean_coords = coords.mean(axis=0)
        atom_rmsf = np.sqrt(((coords - mean_coords) ** 2).sum(axis=2).mean(axis=0))
        for residue in protein.residues:
            indices = [list(protein.indices).index(index) for index in residue.atoms.indices if index in set(protein.indices)]
            if indices:
                rmsf_rows.append({"replica": replica, "resid": residue.resid, "resname": residue.resname, "rmsf_angstrom": float(atom_rmsf[indices].mean())})
        replica_rmsd = pd.DataFrame([row for row in rmsd_rows if row["replica"] == replica])
        protein_values = replica_rmsd[replica_rmsd["selection"] == "protein_backbone"]["rmsd_angstrom"].to_numpy()
        ligand_values = replica_rmsd[replica_rmsd["selection"] == "ligand_after_protein_fit"]["rmsd_angstrom"].to_numpy()
        ligand_last = ligand_values[len(ligand_values) // 2 :]
        summaries.append({"replica": replica, "frames": len(universe.trajectory), "duration_ps": float(universe.trajectory.totaltime), "protein_rmsd_last_half_mean": float(protein_values[len(protein_values) // 2 :].mean()), "ligand_rmsd_last_half_mean": float(ligand_last.mean()), "ligand_rmsd_last_half_sd": float(ligand_last.std()), "contact_mean": float(np.mean(ligand_contacts)), "rg_mean": float(np.mean(radii))})
    summary = pd.DataFrame(summaries)
    for name, rows in [("thermodynamics.tsv", thermo_rows), ("rmsd.tsv", rmsd_rows), ("rmsf.tsv", rmsf_rows), ("radius_of_gyration.tsv", rg_rows), ("pocket_contacts.tsv", contact_rows), ("hydrogen_bond_occupancy.tsv", hbond_rows), ("key_distances.tsv", distance_rows)]:
        pd.DataFrame(rows).to_csv(work / name, sep="\t", index=False)
    summary.to_csv(work / "replica_summary.tsv", sep="\t", index=False)
    consistency = {"replicas": len(summary), "ligand_rmsd_between_replica_cv": float(summary["ligand_rmsd_last_half_mean"].std(ddof=1) / summary["ligand_rmsd_last_half_mean"].mean()) if summary["ligand_rmsd_last_half_mean"].mean() else math.nan, "interpretation": "short smoke trajectories validate execution and parsing only; they do not establish stability or mechanism"}
    (work / "replica_consistency.json").write_text(json.dumps(consistency, ensure_ascii=False, indent=2), encoding="utf-8")
    return work

"""One-entrypoint orchestration and evidence audit for the downstream extension."""
from __future__ import annotations

import html
import hashlib
import fcntl
import gzip
import json
import re
import subprocess
import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd

from .config import WorkflowConfig
from .runtime import atomic_json, atomic_text, atomic_tsv, sha256


MODULES = {"beta": "04_beta_states", "enrichment": "05_enrichment",
           "networks": "06_networks", "velocity": "07_velocity"}
STAGES = ("velocity-count", "beta", "enrichment", "networks", "velocity", "report", "audit")


def record_stage(path: Path, stage: str, value: dict) -> None:
    """Preserve other independently running stages when recording completion."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        states = json.loads(path.read_text()) if path.exists() else {}
        states[stage] = value
        atomic_json(path, states)


def archive_metadata_signature(path: Path) -> str:
    listing = subprocess.run(["find", str(path), "-type", "f", "-printf", "%P\t%s\t%T@\n"],
                             check=True, stdout=subprocess.PIPE).stdout
    return hashlib.sha256(b"".join(sorted(listing.splitlines(keepends=True)))).hexdigest()


def verify_archive_guard(run_dir: Path) -> dict:
    path = run_dir / "00_audit/archive_guard.json"
    if not path.exists():
        raise ValueError("Archive baseline was not recorded")
    guard = json.loads(path.read_text())
    if archive_metadata_signature(Path(guard["path"])) != guard["metadata_sha256"]:
        raise ValueError("PFOS1 archive metadata changed during this extension")
    return guard


def verify_execution_inputs(folder: Path) -> dict:
    path = folder / "execution_inputs.json"
    if not path.exists():
        raise ValueError(f"Missing execution input provenance: {path}")
    payload = json.loads(path.read_text())
    if payload.get("verification_status") != "PASS":
        raise ValueError(f"Execution end-of-run verification did not pass: {path}")
    if not payload.get("files"):
        raise ValueError(f"Empty execution input provenance: {path}")
    for item in payload["files"]:
        source = Path(item["path"])
        if "PFOS1" in source.resolve().parts:
            raise ValueError("Archived company data cannot be an execution input")
        if not source.is_file() or sha256(source) != item["sha256"]:
            raise ValueError(f"Execution input changed: {source}")
    return payload


def counting_code_guard(root: Path, folder: Path, reference: Path) -> dict:
    """Bind counter code/environment and verify a resumed decompressed annotation."""
    paths = [root / "src/pfos_scrna/extension_velocity.py",
             root / "environments/velocity/pyproject.toml",
             root / "environments/velocity/uv.lock"]
    files = [{"path": str(path.resolve()), "sha256": sha256(path)} for path in paths]
    receipt = folder / "counting_code_guard.json"
    if receipt.exists() and json.loads(receipt.read_text())["files"] != files:
        raise ValueError("Velocity counter code/environment changed; use a new run-id")
    generated = folder / "cellranger.genes.gtf"
    annotation = None
    if generated.exists():
        source = reference / "genes/genes.gtf.gz"
        digest = hashlib.sha256()
        with gzip.open(source, "rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        if sha256(generated) != digest.hexdigest():
            raise ValueError("Decompressed counting GTF differs from the frozen Cell Ranger annotation")
        annotation = {"path": str(generated.resolve()), "sha256": digest.hexdigest(),
                      "source": str(source.resolve()), "source_sha256": sha256(source)}
    mask = None
    mask_receipt = folder / "repeat_mask_preparation.json"
    if mask_receipt.exists():
        mask = json.loads(mask_receipt.read_text())
        if (mask.get("status") != "PASS" or mask.get("source_contigs") != 61
                or mask.get("normalized_contigs") != 61 or mask.get("gtf_contigs") != 38):
            raise ValueError("Normalized repeat mask lacks complete official contig coverage")
        expected_paths = {"source_mask": reference.parent / "ucsc_repeatmasker/mm39.repeat_mask.gtf",
                          "alias": root / "reference/velocity/mm39.chromAlias.tsv",
                          "fasta_fai": reference / "fasta/genome.fa.fai", "gtf": generated,
                          "normalized_mask": folder / "cellranger.repeat_mask.gtf"}
        for key, expected_path in expected_paths.items():
            if (Path(mask[key]) != expected_path.resolve()
                    or sha256(expected_path) != mask[key + "_sha256"]):
                raise ValueError(f"Normalized repeat mask reference changed: {key}")
    payload = {"files": files, "verified_decompressed_annotation": annotation,
               "verified_repeat_mask": mask}
    if receipt.exists():
        payload["first_recorded_at"] = json.loads(receipt.read_text())["first_recorded_at"]
    else:
        payload["first_recorded_at"] = datetime.now().isoformat()
    atomic_json(receipt, payload)
    return payload


def counting_output_paths() -> set[str]:
    shared = {"counting_status.tsv", "counting_inputs.json", "counting_code_guard.json", "cellranger.genes.gtf",
              "cellranger.repeat_mask.gtf", "repeat_mask_preparation.json"}
    for library in ("PFOS_Lib01", "WT_Lib01"):
        shared.update(f"{library}/{name}" for name in ("beta_barcodes.tsv", "bam_input.json", "filter_completed.json",
                      "beta.coordinate.bam", "cellsorted_beta.coordinate.bam", "count_completed.json", f"{library}.loom"))
    return shared


def verify_counting_seal(folder: Path) -> None:
    path = folder / "counting_output_hashes.tsv"
    if not path.exists():
        raise ValueError("Completed splice counting lacks an intermediate/output SHA-256 seal")
    hashes = pd.read_csv(path, sep="\t")
    if set(hashes.path) != counting_output_paths() or hashes.path.duplicated().any():
        raise ValueError("Incomplete splice counting output seal")
    for row in hashes.itertuples(index=False):
        source = _local_file(folder, row.path)
        if source.stat().st_size != row.bytes or sha256(source) != row.sha256:
            raise ValueError(f"Sealed counting intermediate/output changed: {source}")


def seal_counting_outputs(folder: Path) -> None:
    """Only seal after both counters successfully finish and verify their inputs."""
    if (folder / "counting_output_hashes.tsv").exists():
        verify_counting_seal(folder)
        return
    status = pd.read_csv(folder / "counting_status.tsv", sep="\t")
    if set(status.library_id) != {"PFOS_Lib01", "WT_Lib01"} or not status.status.eq("completed").all():
        raise ValueError("Cannot seal incomplete library splice counts")
    rows = []
    for value in sorted(counting_output_paths()):
        path = _local_file(folder, value)
        rows.append({"path": value, "bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(folder / "counting_output_hashes.tsv", pd.DataFrame(rows))


def completed_counting_or_fresh(folder: Path) -> dict | None:
    """Never resume partial barcode-sorted/filtered BAMs on existence or size alone."""
    if (folder / "counting_status.tsv").exists():
        inputs = json.loads((folder / "counting_inputs.json").read_text())
        if not inputs.get("files"):
            raise ValueError("Missing counting input fingerprints")
        for item in inputs["files"]:
            if sha256(Path(item["path"])) != item["sha256"]:
                raise ValueError(f"Counting input/reference changed: {item['path']}")
        verify_counting_seal(folder)
        return {"status": "completed", "libraries": pd.read_csv(folder / "counting_status.tsv", sep="\t").to_dict("records")}
    if any(folder.rglob("*.bam")) or any(folder.rglob("*.loom")) or any(folder.rglob("filter_completed.json")):
        raise ValueError("Incomplete splice counting contains unsealed intermediates; use a new run-id, not automatic resume")
    return None


def _local_file(folder: Path, value: str) -> Path:
    path = (folder / value).resolve()
    if not path.is_relative_to(folder.resolve()) or not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing, empty or external result: {value}")
    return path


def module_output_paths(manifest: dict) -> set[str]:
    values = set(manifest["tables"] + manifest.get("objects", []) + ["report_manifest.json", "execution_inputs.json"])
    for figure in manifest["figures"]:
        values.update(figure[key] for key in ("png", "pdf", "source"))
    return values


def check_module(folder: Path, *, verify_inputs: bool = True) -> dict:
    """A status-only table or a chapter title never counts as completed analysis."""
    manifest_path = folder / "report_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "completed":
        raise ValueError(f"Module has not completed: {folder}")
    if not manifest.get("method") or not manifest.get("metrics"):
        raise ValueError("Missing method/actual metrics")
    if not manifest.get("tables") or not manifest.get("figures"):
        raise ValueError("Completed module requires real tables and diagnostic figures")
    for value in manifest["tables"]:
        _local_file(folder, value)
    for figure in manifest["figures"]:
        if not figure.get("interpretation"):
            raise ValueError("Each figure requires an interpretation")
        for key in ("png", "pdf", "source"):
            _local_file(folder, figure[key])
        if _local_file(folder, figure["pdf"]).read_bytes()[:5] != b"%PDF-":
            raise ValueError("Invalid PDF")
        if _local_file(folder, figure["png"]).read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
            raise ValueError("Invalid PNG")
    for value in manifest.get("objects", []):
        _local_file(folder, value)
    if verify_inputs:
        verify_execution_inputs(folder)
    seal = folder / "module_output_hashes.tsv"
    if seal.exists():
        hashes = pd.read_csv(seal, sep="\t")
        if hashes.path.duplicated().any() or set(hashes.path) != module_output_paths(manifest):
            raise ValueError("Incomplete or inconsistent module output hash coverage")
        for row in hashes.itertuples(index=False):
            path = _local_file(folder, row.path)
            if path.stat().st_size != row.bytes or sha256(path) != row.sha256:
                raise ValueError(f"Sealed output changed: {path}")
    return manifest


def seal_module(folder: Path) -> None:
    manifest = check_module(folder)
    seal = folder / "module_output_hashes.tsv"
    if seal.exists():
        return
    values = module_output_paths(manifest)
    atomic_tsv(seal, pd.DataFrame([{"path": value, "bytes": (folder / value).stat().st_size,
                                   "sha256": sha256(folder / value)} for value in sorted(values)]))


def chinese_figure_readings(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    readings, current = {}, None
    for line in path.read_text().splitlines():
        if line.startswith("## "):
            current = line[3:].strip()
            readings[current] = {}
        elif current and any(line.startswith(key + "：") for key in ("标题", "分布", "解释", "边界")):
            key, value = line.split("：", 1)
            readings[current][key] = value
    return readings


def freeze_source(config: WorkflowConfig, source_run: Path, run_dir: Path) -> dict:
    source_run = source_run.resolve()
    if not source_run.is_relative_to(config.run_root) or source_run == run_dir.resolve():
        raise ValueError("Source must be a different formal PFOS run")
    audit = json.loads((source_run / "00_audit/audit_status.json").read_text())
    if audit.get("status") != "PASS":
        raise ValueError("Parent audit did not pass")
    h5ad = source_run / "04_annotation/annotated.h5ad"
    parent_hashes = pd.read_csv(source_run / "00_audit/output_hashes.tsv", sep="\t").set_index("path")
    actual_hash = sha256(h5ad)
    if actual_hash != parent_hashes.loc["04_annotation/annotated.h5ad", "sha256"]:
        raise ValueError("Source H5AD no longer matches the completed parent audit")
    files = [config.source, config.samples_path, config.root / "docs/extension_analysis_plan.md",
             source_run / "00_audit/run_manifest.json", source_run / "00_audit/output_hashes.tsv"]
    payload = {"source_run": str(source_run), "input_h5ad": str(h5ad),
               "input_h5ad_sha256": actual_hash,
               "files": [{"path": str(p), "sha256": sha256(p)} for p in files],
               "condition_inference": "NO-GO_n1_per_condition"}
    path = run_dir / "00_audit/source_manifest.json"
    if path.exists():
        if json.loads(path.read_text()) != payload:
            raise ValueError("Extension input/config/reference lineage changed; use a new run-id")
    else:
        atomic_json(path, payload)
    guard_path = run_dir / "00_audit/archive_guard.json"
    if not guard_path.exists():
        archive = config.root.parent / "PFOS1"
        # Signature verified at the start of this authorized extension, before computation.
        expected = "d05cd62a05b166e5159324f302f0c0d600a2084d5a06d5ba18b88561dd2301d2"
        if archive_metadata_signature(archive) != expected:
            raise ValueError("PFOS1 differs from the pre-extension read-only baseline")
        atomic_json(guard_path, {"path": str(archive), "metadata_sha256": expected})
    verify_archive_guard(run_dir)
    return payload


def extension_report(run_dir: Path) -> dict:
    manifests = {name: check_module(run_dir / folder) for name, folder in MODULES.items()}
    root = run_dir.parents[1]
    beta_text = root / "docs/extension_beta_report_zh.md"
    readings = chinese_figure_readings(beta_text)
    rows, claims = [], []
    parts = ["<!doctype html><html lang='zh-CN'><meta charset='utf-8'>",
             "<style>body{max-width:1100px;margin:35px auto;font:16px/1.65 sans-serif}img{max-width:100%}"
             "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:8px}section{margin:40px 0}</style>",
             "<title>PFOS胰岛扩展分析</title><h1>PFOS 胰岛下游补齐报告</h1>",
             "<p>PFOS与WT各一份文库。本报告补充功能注释、β细胞状态、库内通讯与pooled调控候选；"
             "所有条件对照均为两库描述，不能推断动物群体效应。Cluster marker使用log2FC≥0.58、原P&lt;0.05；"
             "marker与基因列表富集P/BH仅用于注释，不能移作PFOS条件显著性。</p>",
             "<p>ORA主展示使用原P&lt;0.05且至少2个命中，完整结果同时保留BH，原P入选不等于BH校正后显著。"
             "通路ES未做置换和标准化，其正负只描述该基因集在本次排名中的分布；不同大小基因集的ES幅度不能直接解释为效应强弱。</p>",
             "<p>本报告由正式自有数据重新计算，公司数据未进入计算。"
             "Mpc2仅作成年代谢调控哨兵/调控重塑定位，不作遗传或因果中介证明。</p>"]
    lineage = json.loads((run_dir / "00_audit/source_manifest.json").read_text())
    parent_report = Path(lineage["source_run"]) / "09_reports/PFOS_islet_scrna_report_zh.html"
    if parent_report.is_file():
        relative = Path("../..") / parent_report.parent.parent.name / "09_reports" / parent_report.name
        parts.append(f"<p><a href='{html.escape(str(relative), quote=True)}'>基础报告：原始输入、QC、全图注释与描述性定位</a>。"
                     "本页只展示本轮新增结果；主计数和既有注释未被下游扩展覆盖。</p>")
    identity_path = run_dir / "04_beta_states/tables/beta_identity_cluster_review.tsv"
    if identity_path.exists():
        review = pd.read_csv(identity_path, sep="\t")
        columns = ["beta_cluster", "n_cells", "median_total_counts", "median_n_genes_by_counts",
                   "median_pct_counts_mt", "median_beta_identity_genes_detected", "cluster_review_status"]
        parts.append("<h2>先读：β候选集的注释不确定性</h2><p>β深入分析继承上游β标签，"
                     "不代表所有细胞都具有同等身份可信度。低复杂度、非同源激素残留与低β身份支持单独标记；"
                     "主分析保留所有上游β标签细胞，排除疑点群只作描述性敏感性，不事后重定主QC或强制命名新亚型。</p>")
        parts.append(review[columns].to_html(index=False, escape=True))
    for module, manifest in manifests.items():
        folder = MODULES[module]
        heading = {"beta": "β候选细胞的亚群、状态与连续体诊断", "enrichment": "GO/KEGG功能富集与描述性通路排序",
                   "networks": "文库内通讯与pooled调控候选", "velocity": "自有BAM剪接计数与逐库动力学诊断"}[module]
        parts.append(f"<section><h2>{heading}</h2><p>{html.escape(str(manifest['method']))}</p>")
        parts.append(f"<p>本模块交付 {len(manifest['figures'])} 张配对图件；每张图下方提供PDF与源数据。</p>"
                     f"<details><summary>实际运行数字与参数</summary><pre>"
                     f"{html.escape(json.dumps(manifest['metrics'], ensure_ascii=False, indent=2))}</pre></details>")
        for i, fig in enumerate(manifest["figures"], 1):
            href = lambda value: html.escape(f"../{folder}/{value}", quote=True)
            title = str(fig.get("title", f"{module}-{i}"))
            interpretation = str(fig["interpretation"])
            chinese = readings.get(Path(fig["png"]).stem) if module == "beta" else None
            if chinese:
                title = chinese["标题"]
                interpretation = "\n\n".join(f"{key}：{chinese[key]}" for key in ("分布", "解释", "边界"))
            parts.append(f"<h3>{html.escape(title)}</h3><img loading='lazy' src='{href(fig['png'])}'>"
                         f"<p style='white-space:pre-line'>{html.escape(interpretation)}</p>"
                         f"<p><a href='{href(fig['pdf'])}'>PDF</a> · <a href='{href(fig['source'])}'>源数据TSV</a></p>")
            rows.append({"module": module, "title": title, **{k: f"{folder}/{fig[k]}" for k in ("png", "pdf", "source")}})
            claims.append({"claim_id": f"{module}-{i:03}", "claim": interpretation,
                           "figure": f"{folder}/{fig['png']}", "result": f"{folder}/{fig['source']}",
                           "scope": "annotation_or_library_descriptive_not_condition_inference"})
        parts.append("<details><summary>完整结果表与数据对象</summary><ul>")
        for table in manifest["tables"] + manifest.get("objects", []):
            href = html.escape(f"../{folder}/{table}", quote=True)
            parts.append(f"<li><a href='{href}'>{html.escape(table)}</a></li>")
        parts.append("</ul></details></section>")
    parts.append("<h2>未开放的结论</h2><p>动物层面的PFOS条件差异、差异通讯、条件特异regulon、"
                 "成年胰岛发育轨迹和普通胰岛inferCNV不具备本项目证据条件。速度与DPT的方向解释以各自诊断门为准。</p></html>")
    reports = run_dir / "09_reports"
    eligibility = [{"analysis": name, "status": "completed", "scope": "annotation_or_library_descriptive",
                    "result": f"{MODULES[name]}/report_manifest.json"} for name in manifests]
    eligibility.extend({"analysis": name, "status": "NO-GO", "scope": "n1_per_condition_not_identifiable",
                        "result": "00_audit/source_manifest.json"}
                       for name in ("condition_expression_inference", "condition_composition_inference",
                                    "differential_communication", "condition_specific_regulon", "condition_velocity_difference"))
    eligibility.append({"analysis": "inferCNV", "status": "NO-GO", "scope": "normal_islets_no_CNV_hypothesis_or_qualified_reference", "result": ""})
    trajectory = json.loads((run_dir / "04_beta_states/trajectory_gate.json").read_text())
    eligibility.append({"analysis": "directional_beta_pseudotime", "status": trajectory["status"],
                        "scope": trajectory["interpretation"], "result": "04_beta_states/trajectory_gate.json"})
    velocity_gate = run_dir / "07_velocity/tables/velocity_eligibility.tsv"
    if velocity_gate.exists():
        for row in pd.read_csv(velocity_gate, sep="\t").itertuples(index=False):
            eligibility.append({"analysis": f"{row.library_id}_local_velocity_projection", "status": row.status,
                                "scope": "within-library beta dynamics screening; not time, fate or condition effect",
                                "result": "07_velocity/tables/velocity_eligibility.tsv"})
    atomic_tsv(run_dir / "00_audit/analysis_eligibility.tsv", pd.DataFrame(eligibility))
    atomic_text(reports / "report.html", "\n".join(parts))
    atomic_tsv(reports / "figure_inventory.tsv", pd.DataFrame(rows))
    atomic_tsv(reports / "claim_result_map.tsv", pd.DataFrame(claims))
    report_inputs = [Path(__file__), beta_text] + [run_dir / folder / "report_manifest.json" for folder in MODULES.values()]
    if identity_path.exists():
        report_inputs.append(identity_path)
    atomic_json(reports / "report_inputs.json", {"files": [{"path": str(p.resolve()), "sha256": sha256(p)} for p in report_inputs]})
    return {"status": "completed", "modules": len(manifests), "figures": len(rows)}


def audit_extension_data(run_dir: Path) -> list[dict]:
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    from .audit import REQUIRED_OBS, REQUIRED_VAR
    beta = ad.read_h5ad(run_dir / "04_beta_states/beta_states.h5ad")
    checks = []
    def check(name: str, passed: bool, detail: str = "") -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})
    check("beta_unique_indices", beta.obs_names.is_unique and beta.var_names.is_unique)
    check("beta_required_obs_var", REQUIRED_OBS <= set(beta.obs) and REQUIRED_VAR <= set(beta.var))
    check("beta_only_source_labels", bool(beta.obs.cell_type.eq("beta").all()))
    check("beta_no_condition_inference", not bool(beta.obs.condition_inference_eligible.any()))
    check("beta_sparse_nonnegative_log_X", sp.issparse(beta.X) and np.isfinite(beta.X.data).all() and (beta.X.data >= 0).all())
    source_path = Path(json.loads((run_dir / "00_audit/source_manifest.json").read_text())["input_h5ad"])
    source = ad.read_h5ad(source_path)
    expected_beta = source.obs_names[source.obs.cell_type.eq("beta")]
    check("beta_exact_parent_cell_set", set(beta.obs_names) == set(expected_beta))
    check("beta_exact_parent_gene_set", set(beta.var_names) == set(source.var_names))
    parent_beta = source[beta.obs_names, beta.var_names]
    for layer in ("counts", "counts_raw"):
        matrix = beta.layers[layer]
        check(f"beta_{layer}_sparse_integer_nonnegative", sp.issparse(matrix) and np.isfinite(matrix.data).all()
              and (matrix.data >= 0).all() and np.equal(matrix.data, np.floor(matrix.data)).all())
        check(f"beta_{layer}_unchanged_from_parent", (matrix != parent_beta.layers[layer]).nnz == 0)
    del parent_beta, source
    for library in sorted(beta.obs.library_id.astype(str).unique()):
        splice = ad.read_h5ad(run_dir / "07_velocity" / f"{library}_splice_counts.h5ad")
        expected_cells = set(beta.obs_names[beta.obs.library_id.eq(library)])
        check(f"{library}_splice_unique_complete_barcodes", splice.obs_names.is_unique and set(splice.obs_names) == expected_cells)
        check(f"{library}_splice_unique_genes", splice.var_names.is_unique and splice.var.gene_ids.is_unique)
        for layer in ("spliced", "unspliced", "ambiguous"):
            matrix = splice.layers[layer]
            check(f"{library}_{layer}_integer_counts", sp.issparse(matrix) and np.isfinite(matrix.data).all()
                  and (matrix.data >= 0).all() and np.equal(matrix.data, np.floor(matrix.data)).all())
        del splice
    eligibility = pd.read_csv(run_dir / "00_audit/analysis_eligibility.tsv", sep="\t")
    required = {"condition_expression_inference", "condition_composition_inference", "differential_communication",
                "condition_specific_regulon", "condition_velocity_difference"}
    selected = eligibility.loc[eligibility.analysis.isin(required)]
    check("n1_all_condition_claims_closed", set(selected.analysis) == required and bool(selected.status.eq("NO-GO").all()))
    return checks


def audit_counting_lineage(run_dir: Path) -> dict:
    folder = run_dir / "07_velocity/counts"
    verify_counting_seal(folder)
    guard = json.loads((folder / "counting_code_guard.json").read_text())
    inputs = json.loads((folder / "counting_inputs.json").read_text())
    for item in guard["files"] + inputs["files"]:
        if sha256(Path(item["path"])) != item["sha256"]:
            raise ValueError(f"Counting code/input changed: {item['path']}")
    annotation = guard.get("verified_decompressed_annotation")
    if not annotation or sha256(Path(annotation["path"])) != annotation["sha256"]:
        raise ValueError("Counting annotation was not verified or changed")
    mask = guard.get("verified_repeat_mask")
    if not mask or mask.get("status") != "PASS":
        raise ValueError("Counting repeat-mask alias preparation was not verified")
    for key in ("source_mask", "alias", "fasta_fai", "gtf", "normalized_mask"):
        if sha256(Path(mask[key])) != mask[key + "_sha256"]:
            raise ValueError(f"Counting repeat-mask reference changed: {key}")
    status = pd.read_csv(folder / "counting_status.tsv", sep="\t")
    if set(status.library_id) != {"PFOS_Lib01", "WT_Lib01"} or not status.status.eq("completed").all():
        raise ValueError("Both library splice counts must complete")
    parent = Path(json.loads((run_dir / "00_audit/source_manifest.json").read_text())["source_run"])
    parent_hashes = pd.read_csv(parent / "00_audit/output_hashes.tsv", sep="\t").set_index("path")
    for row in status.itertuples(index=False):
        lib = folder / row.library_id
        bam = json.loads((lib / "bam_input.json").read_text())
        relative = f"01_cellranger/{row.library_id}/outs/possorted_genome_bam.bam"
        source = parent / relative
        if (Path(bam["path"]) != source.resolve() or bam["sha256"] != parent_hashes.loc[relative, "sha256"]
                or source.stat().st_size != bam["bytes"] or source.stat().st_mtime_ns != bam["mtime_ns"]):
            raise ValueError("Parent BAM no longer matches its actually verified counting receipt")
        completed = json.loads((lib / "count_completed.json").read_text())
        if completed["cells"] != row.beta_cells or sha256(lib / f"{row.library_id}.loom") != completed["loom_sha256"]:
            raise ValueError("Splice completion receipt does not match the actual loom")
    return {"libraries": len(status), "status": "PASS"}


def forbidden_descriptive_columns(columns) -> list[str]:
    """Condition descriptive tables cannot acquire inferential statistics by relabeling."""
    return [str(column) for column in columns if re.search(
        r"(^|_)(p|pval|pvals|pvalue|pvalues|p_value|padj|p_adj|p_adjust|p_adjusted|fdr|qval|qvalue|q_value|nes)($|_)",
        str(column), re.IGNORECASE)]


def audit_expression_network_results(run_dir: Path, modules: dict) -> list[dict]:
    """Audit completed expression modules independently of long-running splice counts."""
    import anndata as ad
    import numpy as np
    source = Path(json.loads((run_dir / "00_audit/source_manifest.json").read_text())["input_h5ad"])
    data = ad.read_h5ad(source, backed="r")
    obs = data.obs.copy()
    data.file.close()
    beta = ad.read_h5ad(run_dir / "04_beta_states/beta_states.h5ad", backed="r")
    beta_groups = set(beta.obs["beta_leiden_0.5"].astype(str))
    beta.file.close()
    libraries = set(obs.library_id.astype(str))
    checks = []
    def check(name: str, passed: bool, detail: str = "") -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})
    enrichment = run_dir / "05_enrichment/tables"
    status = pd.read_csv(enrichment / "ora_status.tsv", sep="\t", dtype={"group": str})
    parent_groups = set(pd.read_csv(source.parent.parent / "03_atlas/cluster_markers.tsv", sep="\t",
                                    usecols=["cluster"], dtype=str).cluster)
    expected_families = {(level, group, database)
                         for level, groups in (("cluster", parent_groups),
                                               ("cell_type", set(obs.cell_type.astype(str)) - {"unknown"}),
                                               ("beta_cluster", beta_groups))
                         for group in groups for database in ("GO_BP", "GO_CC", "GO_MF", "KEGG")}
    actual_families = set(status[["analysis_level", "group", "database"]].itertuples(index=False, name=None))
    check("enrichment_all_groups_databases_accounted", actual_families == expected_families and len(status) == len(expected_families))
    ora = pd.read_csv(enrichment / "ora_full.tsv", sep="\t")
    check("ora_named_terms", bool(ora.term_name.notna().all() and ora.term_name.str.strip().ne("").all()))
    check("ora_valid_p_bh", bool(ora.pvalue.between(0, 1).all() and ora.p_adjust_bh.between(0, 1).all()
                                and (ora.p_adjust_bh + 1e-12 >= ora.pvalue).all()))
    check("ora_nominal_bh_flags_match", bool(ora.passes_nominal_p_0_05.eq(ora.pvalue.lt(.05)).all()
                                           and ora.passes_bh_0_05.eq(ora.p_adjust_bh.lt(.05)).all()))
    marker_summary = pd.read_csv(enrichment / "marker_input_summary.tsv", sep="\t")
    check("marker_user_thresholds", bool(marker_summary.marker_log2fc_threshold.eq(.58).all()
                                        and marker_summary.marker_nominal_pvalue_threshold.eq(.05).all()
                                        and marker_summary.p_adjusted_used_for_marker_selection.eq(False).all()))
    keys = ["analysis_level", "group", "gene_symbol"]
    query = pd.read_csv(enrichment / "marker_query.tsv", sep="\t", dtype=str)
    universe = pd.read_csv(enrichment / "marker_detectable_universe.tsv", sep="\t", dtype=str)
    check("marker_query_within_tested_detectable_background", set(query[keys].itertuples(index=False, name=None))
          <= set(universe[keys].itertuples(index=False, name=None)))
    marker_inputs = {"cluster": (source.parent.parent / "03_atlas/cluster_markers.tsv", "cluster"),
                     "cell_type": (enrichment / "cell_type_markers.tsv", "group"),
                     "beta_cluster": (run_dir / "04_beta_states/tables/beta_cluster_markers_all.tsv", "cluster")}
    for level, (path, group_column) in marker_inputs.items():
        raw = pd.read_csv(path, sep="\t", dtype={group_column: str})
        tested = (np.isfinite(raw.pvals) & np.isfinite(raw.logfoldchanges)
                  & (raw.pct_nz_group.gt(0) | raw.pct_nz_reference.gt(0)))
        positive = tested & raw.logfoldchanges.ge(.58) & raw.pvals.lt(.05)
        for kind, selected, actual in (("background", tested, universe), ("query", positive, query)):
            expected = set(raw.loc[selected, [group_column, "names"]].itertuples(index=False, name=None))
            observed = set(actual.loc[actual.analysis_level.eq(level), ["group", "gene_symbol"]].itertuples(index=False, name=None))
            check(f"{level}_marker_{kind}_recomputed_from_raw_statistics", observed == expected)
    for path in sorted(enrichment.glob("descriptive*.tsv")):
        frame = pd.read_csv(path, sep="\t")
        forbidden = forbidden_descriptive_columns(frame.columns)
        check(f"{path.stem}_no_condition_statistics", not forbidden, ",".join(forbidden))
        if "descriptive_es" in frame:
            check(f"{path.stem}_es_bounds", bool(frame.descriptive_es.between(-1.000001, 1.000001).all()))
    primary_es = pd.read_csv(enrichment / "descriptive_ranked_pathway_es.tsv", sep="\t")
    equal_es = pd.read_csv(enrichment / "descriptive_pathway_es_equal_cell_stability.tsv", sep="\t")
    pathway_key = ["cell_type", "database", "term_id"]
    check("equal_cell_sensitivity_all_primary_pathways", set(primary_es[pathway_key].itertuples(index=False, name=None))
          == set(equal_es[pathway_key].itertuples(index=False, name=None)))
    cells_by_type = obs.groupby(["cell_type", "library_id"], observed=True).size().unstack(fill_value=0)
    expected_balanced = equal_es.cell_type.map(cells_by_type.min(axis=1).clip(upper=1000))
    check("equal_cell_sensitivity_sampling_sizes", bool(equal_es.balanced_cells_per_condition.gt(0).all()
          and equal_es.balanced_cells_per_condition.eq(expected_balanced).all()))
    check("equal_cell_sensitivity_50_attempts_and_valid_counts", bool(equal_es.requested_repeats.eq(50).all()
          and equal_es.valid_repeats.between(0, 50).all()), "valid counts may be <50 when subsampled gene-set coverage fails")
    estimates = ["median_balanced_es", "q10_balanced_es", "q90_balanced_es", "direction_stability",
                 "median_leading_edge_jaccard", "q10_leading_edge_jaccard", "q90_leading_edge_jaccard"]
    unavailable = equal_es.valid_repeats.eq(0)
    check("equal_cell_unavailable_terms_not_called_stable", bool(equal_es.loc[unavailable, estimates].isna().all().all()))
    check("equal_cell_available_stability_has_valid_values", bool(equal_es.loc[~unavailable, estimates].notna().all().all()
          and equal_es.loc[~unavailable, "direction_stability"].between(0, 1).all()
          and equal_es.loc[~unavailable, "median_leading_edge_jaccard"].between(0, 1).all()))
    networks = run_dir / "06_networks/tables"
    lr = pd.read_csv(networks / "communication_lr_full.tsv", sep="\t")
    check("liana_each_library_known_cells", set(lr.library_id) == libraries and "unknown" not in set(lr.source) | set(lr.target))
    check("liana_all_complex_subunits", bool(lr.all_subunits_pass.eq(True).all()))
    check("liana_priority_rule", bool(lr.priority_candidate.eq(lr.magnitude_rank.le(.05)).all()))
    for side in ("ligand", "receptor"):
        values = lr[f"{side}_subunit_expr_props"].astype(str).str.split(";").explode().astype(float)
        check(f"liana_{side}_actual_subunit_detection", bool(values.ge(.1 - 1e-6).all()))
    stable = pd.read_csv(networks / "regulon_seed_stable_edges.tsv", sep="\t")
    seeds = set(stable.seeds.astype(str).str.split(";").explode())
    check("regulon_three_seed_consensus", seeds == {"1907", "2909", "3911"}
          and bool(stable.seed_support.between(2, 3).all())
          and bool(stable.seeds.astype(str).str.split(";").map(lambda x: len(set(x))).eq(stable.seed_support).all()))
    targets = pd.read_csv(networks / "regulon_motif_supported_targets.tsv", sep="\t")
    check("regulon_motif_targets_in_stable_edges", set(targets[["TF", "target"]].itertuples(index=False, name=None))
          <= set(stable[["TF", "target"]].itertuples(index=False, name=None)))
    auc = pd.read_csv(networks / "regulon_auc_all_cells.tsv.gz", sep="\t", index_col=0)
    check("regulon_aucell_all_source_cells", auc.index.is_unique and set(auc.index) == set(obs.index))
    check("regulon_auc_finite_nonnegative", bool(np.isfinite(auc.to_numpy()).all() and (auc.to_numpy() >= 0).all()))
    check("regulon_declared_metrics_match", modules["networks"]["metrics"]["regulon"]["grn_seeds"] == 3
          and modules["networks"]["metrics"]["regulon"]["scored_cells"] == len(auc)
          and modules["networks"]["metrics"]["regulon"]["regulons"] == auc.shape[1])
    return checks


def audit_velocity_results(run_dir: Path) -> list[dict]:
    """A final velocity audit requires actual per-library model diagnostics."""
    libraries = {"PFOS_Lib01", "WT_Lib01"}
    checks = []
    def check(name: str, passed: bool, detail: str = "") -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})
    velocity = pd.read_csv(run_dir / "07_velocity/tables/velocity_eligibility.tsv", sep="\t")
    check("velocity_each_library_actual_screen", set(velocity.library_id) == libraries and len(velocity) == len(libraries)
          and bool(velocity.status.str.startswith(("PASS_exploratory_local_projection", "NO-GO_")).all()))
    from .extension_velocity_model import CRITERIA, projection_eligible
    for row in velocity.to_dict("records"):
        library = row["library_id"]
        sensitivity_path = run_dir / f"07_velocity/tables/{library}_neighbor_sensitivity.tsv"
        if row["status"] in {"PASS_exploratory_local_projection", "NO-GO_directional_projection_screen"}:
            sensitivity = pd.read_csv(sensitivity_path, sep="\t")
            eligible = projection_eligible(row, sensitivity.to_dict("records"))
            check(f"{library}_velocity_direction_gate_recomputed", row["status"].startswith("PASS_") == eligible
                  and set(sensitivity.neighbors) == {15, 50})
        elif row["status"] == "NO-GO_insufficient_splice_coverage":
            check(f"{library}_velocity_coverage_nogo_verified", row["coverage_genes"] < CRITERIA["min_coverage_genes"])
        elif row["status"] == "NO-GO_insufficient_unrelaxed_velocity_genes":
            check(f"{library}_velocity_gene_nogo_verified", row["stochastic_velocity_genes"] < CRITERIA["min_velocity_genes"])
        else:
            check(f"{library}_velocity_status_is_implemented_scientific_gate", False, str(row["status"]))
        projection = run_dir / f"07_velocity/figures/{library}_local_velocity.png"
        check(f"{library}_velocity_projection_only_when_eligible", projection.is_file() == row["status"].startswith("PASS_"))
    return checks


def audit_scientific_results(run_dir: Path, modules: dict) -> list[dict]:
    """Final PASS always includes expression/network and actual velocity results."""
    return audit_expression_network_results(run_dir, modules) + audit_velocity_results(run_dir)


def audit_figure_files(run_dir: Path, manifests: dict) -> list[dict]:
    from PIL import Image
    renderer = shutil.which("gs")
    if renderer is None:
        raise RuntimeError("Ghostscript is required to validate final PDF rendering")
    rows = []
    for module, manifest in manifests.items():
        seen = set()
        folder = run_dir / MODULES[module]
        for figure in manifest["figures"]:
            png, pdf, source = [_local_file(folder, figure[key]) for key in ("png", "pdf", "source")]
            if png.stem in seen or not (png.stem == pdf.stem == source.stem):
                raise ValueError("Duplicate figure or PNG/PDF/source stem mismatch")
            seen.add(png.stem)
            with Image.open(png) as picture:
                width, height = picture.size
                dpi = picture.info.get("dpi", (0, 0))
                picture.verify()
            if min(dpi) < 299 or min(width, height) < 200:
                raise ValueError(f"PNG resolution contract failed: {png}")
            rendered = subprocess.run([renderer, "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=nullpage", str(pdf)],
                                      capture_output=True, text=True)
            if rendered.returncode:
                raise ValueError(f"PDF rendering failed: {pdf}: {rendered.stderr[-500:]}")
            table = pd.read_csv(source, sep="\t", nrows=2)
            if table.empty:
                raise ValueError(f"Figure has empty source table: {source}")
            rows.append({"module": module, "figure_stem": png.stem, "width_px": width, "height_px": height,
                         "dpi_x": dpi[0], "dpi_y": dpi[1], "source_columns": len(table.columns),
                         "pdf_renderer": renderer, "status": "PASS"})
    atomic_tsv(run_dir / "00_audit/figure_contract_checks.tsv", pd.DataFrame(rows))
    return rows


def audit_extension(run_dir: Path) -> dict:
    archive_guard = verify_archive_guard(run_dir)
    counting_lineage = audit_counting_lineage(run_dir)
    for folder in MODULES.values():
        if not (run_dir / folder / "module_output_hashes.tsv").exists():
            raise ValueError(f"Module has not been sealed by the extension entrypoint: {folder}")
    modules = {name: check_module(run_dir / folder) for name, folder in MODULES.items()}
    figure_checks = audit_figure_files(run_dir, modules)
    report = run_dir / "09_reports/report.html"
    if not report.exists():
        raise ValueError("Missing final report")
    for value in re.findall(r"(?:src|href)='([^']+)'", report.read_text()):
        if not (report.parent / value).is_file():
            raise ValueError(f"Broken report link: {value}")
    freeze = json.loads((run_dir / "00_audit/source_manifest.json").read_text())
    if sha256(Path(freeze["input_h5ad"])) != freeze["input_h5ad_sha256"]:
        raise ValueError("Source H5AD changed")
    for item in freeze["files"]:
        if sha256(Path(item["path"])) != item["sha256"]:
            raise ValueError(f"Frozen source/config changed: {item['path']}")
    for item in json.loads((run_dir / "09_reports/report_inputs.json").read_text())["files"]:
        if sha256(Path(item["path"])) != item["sha256"]:
            raise ValueError(f"Report source changed: {item['path']}")
    checks = audit_extension_data(run_dir)
    atomic_tsv(run_dir / "00_audit/data_contract_checks.tsv", pd.DataFrame(checks))
    if any(row["status"] != "PASS" for row in checks):
        raise ValueError("Extension data contract failed; see data_contract_checks.tsv")
    scientific_checks = audit_scientific_results(run_dir, modules)
    atomic_tsv(run_dir / "00_audit/scientific_result_checks.tsv", pd.DataFrame(scientific_checks))
    if any(row["status"] != "PASS" for row in scientific_checks):
        raise ValueError("Scientific result contract failed; see scientific_result_checks.tsv")
    rows = []
    excluded = {"00_audit/output_hashes.tsv", "00_audit/audit_status.json", "00_audit/extension_stages.json"}
    for path in sorted(run_dir.rglob("*")):
        relative = str(path.relative_to(run_dir))
        if path.is_file() and not path.is_symlink() and relative not in excluded:
            # Counting intermediates have their own receipts; hash final looms and tables here.
            if path.suffix == ".bam" or path.name.endswith(".log"):
                continue
            rows.append({"path": relative, "bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_tsv(run_dir / "00_audit/output_hashes.tsv", pd.DataFrame(rows))
    result = {"status": "PASS", "modules": list(modules), "figures": sum(len(x["figures"]) for x in modules.values()),
              "counting_lineage": counting_lineage,
              "figure_contract_passed": len(figure_checks),
              "scientific_result_contract_passed": len(scientific_checks),
              "pfos1_metadata_sha256_unchanged": archive_guard["metadata_sha256"],
              "files_hashed": len(rows), "hash_exclusions": sorted(excluded) + ["intermediate *.bam (sealed separately in counting_output_hashes.tsv)", "runtime *.log"],
              "finished_at": datetime.now().isoformat()}
    atomic_json(run_dir / "00_audit/audit_status.json", result)
    return result


def run_extension(config: WorkflowConfig, source_run: Path, run_id: str, stages: list[str]) -> dict:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", run_id):
        raise ValueError("Invalid run-id")
    if not stages or set(stages) - set(STAGES):
        raise ValueError(f"Stages must be selected from {STAGES}")
    run_dir = config.run_root / run_id
    source_run = config.resolve(source_run)
    lineage = freeze_source(config, source_run, run_dir)
    input_h5ad = Path(lineage["input_h5ad"])
    state_file = run_dir / "00_audit/extension_stages.json"
    states = json.loads(state_file.read_text()) if state_file.exists() else {}
    for stage in stages:
        states[stage] = {"status": "running", "started_at": datetime.now().isoformat()}
        record_stage(state_file, stage, states[stage])
        try:
            folder = run_dir / MODULES[stage] if stage in MODULES else None
            if stage == "velocity":
                audit_counting_lineage(run_dir)
            if folder is not None and (folder / "report_manifest.json").exists():
                result = check_module(folder)
            elif stage == "velocity-count":
                from .extension_velocity import count_velocity
                counts_folder = run_dir / "07_velocity/counts"
                counting_code_guard(config.root, counts_folder, config.reference)
                result = completed_counting_or_fresh(counts_folder)
                if result is None:
                    result = count_velocity(input_h5ad, source_run, counts_folder, config.reference)
                    counting_code_guard(config.root, counts_folder, config.reference)
                    seal_counting_outputs(counts_folder)
            elif stage == "beta":
                from .extension_beta import run_beta
                result = run_beta(input_h5ad, folder)
            elif stage == "enrichment":
                from .extension_enrichment import run_enrichment
                result = run_enrichment(input_h5ad, source_run, folder, config.root / "reference/extension_gene_sets")
            elif stage == "networks":
                from .extension_networks import run_networks
                result = run_networks(input_h5ad, folder, config.reference.parent / "pyscenic_mouse_v10")
            elif stage == "velocity":
                from .extension_velocity_model import model_velocity
                check_module(run_dir / "04_beta_states")
                result = model_velocity(run_dir / "04_beta_states/beta_states.h5ad", folder)
            elif stage == "report":
                result = extension_report(run_dir)
            else:
                result = audit_extension(run_dir)
            if folder is not None:
                seal_module(folder)
            states[stage].update(status="completed", finished_at=datetime.now().isoformat(), result=result)
        except BaseException as exc:
            states[stage].update(status="failed", error=f"{type(exc).__name__}: {exc}")
            record_stage(state_file, stage, states[stage])
            raise
        record_stage(state_file, stage, states[stage])
    return states

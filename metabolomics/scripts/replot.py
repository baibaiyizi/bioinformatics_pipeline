#!/usr/bin/env python3
"""只读冻结结果，在同一文件系统完成重绘后发布；不调用统计模块。"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from plot_lifecycle import PlotRun
from datetime import datetime, timezone


ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def scientific_files(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): digest(p) for p in sorted(root.rglob("*"))
            if p.is_file() and p.suffix.lower() in {".csv", ".tsv", ".rds"}
            and p.name != "figure_manifest.tsv"}


def figure_descriptions(result: Path, comparison: str) -> dict:
    descriptions = {'sample_run_order_total_signal': ['检查真实进样顺序中的总体信号漂移。', '全部登记进样按真实顺序展示总峰面积；这是峰表代理，不是原始 TIC。', '用于技术诊断，不能据此自动剔除样本或确认 pooled QC 身份。'], 'sample_log2_intensity_boxplot': ['比较归一化后各进样的强度分布。', '全部保留特征按样本展示 log2 归一化强度，箱线图保留离群点。', '样本分布用于归一化与异常检查，不代表组间差异检验。'], 'qc_rsd_distribution': ['评估标记 QC 进样中的特征重复性。', '展示 QC RSD 全分布及配置的参考线。', 'QC 来源尚未确认；本次未启用 QC RSD 硬过滤。'], 'sample_rla_summary': ['检查强度相对偏移是否随进样顺序变化。', '按真实顺序展示各进样的 RLA 中位数与零参考线。', '这是技术偏移诊断，不作为暴露效应或自动删样依据。'], 'pca_POS': ['查看正离子模式的总体样本结构。', '同一冻结矩阵中各进样的前两主成分，保留 WT、PFOS、QC。', 'PCA 只描述样本结构，不能代替组间统计或独立验证。'], 'pca_NEG': ['查看负离子模式的总体样本结构。', '负离子特征单独计算 PCA，沿用原矩阵和全部登记进样。', '不同极性分别诊断，不用 PCA 分离程度筛选样本。'], 'sample_spearman_correlation_POS': ['检查正离子模式的样本相似性。', '冻结 log2 强度的样本 Spearman 相关矩阵；聚类仅用于展示。', '相似性用于技术检查，不构成生物学重复或条件检验。'], 'sample_spearman_correlation_NEG': ['检查负离子模式的样本相似性。', '负离子特征的 Spearman 相关矩阵，保留全部样本标记。', '相关性不用于推断因果，也不触发自动删样。'], 'PFOS_vs_WT_volcano': ['展示冻结的 PFOS−WT 代谢特征统计。', '{counts}；横轴 log2FC，纵轴 −log10 limma P。', 'FDR 支持与探索候选分开保留；未知身份条目仅称代谢特征。'], 'PFOS_vs_WT_top_feature_boxplots': ['查看优先特征在两组进样中的强度分布。', '沿用原 candidate_score 排序的前 20 个非 not_sig 特征，分面尺度与所有样本保留。', '展示项按原结果筛选，不能将这 20 项全部称为 FDR 显著代谢物。'], 'PFOS_vs_WT_top_feature_heatmap': ['展示同一批优先特征的联合分布。', '与箱线图相同的 20 个特征，以原 scaled 矩阵聚类；颜色保留标准化强度含义。', '热图是已选择特征的描述，不作为独立分类性能或额外显著性证据。'], 'pathway_quantitative_statistics': ['展示外部通路定量的探索性组间统计。', '{pathways}；效应、P 与 BH/FDR 保留在完整源表中。', '这是外部通路定量结果；未提供本地成员表，不声称本次重跑 ORA 或建立基因—代谢物机制。']}
    descriptions = {name.replace("PFOS_vs_WT", comparison): value for name, value in descriptions.items()}
    with (result / "04_statistics/tables" / f"{comparison}_summary.csv").open() as handle:
        summary = list(csv.DictReader(handle))
    counts = {kind: sum(int(row["n_features"]) for row in summary if row["classification"] == kind)
              for kind in ("robust", "exploratory", "not_sig")}
    total = sum(counts.values())
    descriptions[f"{comparison}_volcano"][1] = descriptions[f"{comparison}_volcano"][1].format(
        counts=f"{total:,} 个特征中，{counts['robust']:,} 个标为 robust，{counts['exploratory']:,} 个标为 exploratory；沿用原统计分类")
    pathway_file = result / "06_pathway/tables/pathway_quantitative_statistics.csv"
    pathways = []
    if pathway_file.is_file():
        with pathway_file.open() as handle:
            pathways = list(csv.DictReader(handle))
    bh = sum(float(row["fdr"]) < .05 for row in pathways if row.get("fdr") not in (None, "", "NA"))
    descriptions["pathway_quantitative_statistics"][1] = descriptions["pathway_quantitative_statistics"][1].format(
        pathways=f"{len(pathways)} 条已有通路全部展示，{bh} 条 BH<0.05")
    descriptions = {name: [text.replace("PFOS−WT", comparison.replace("_vs_", "−")).replace("WT、PFOS、QC", "配置登记的生物组与 QC") for text in values] for name, values in descriptions.items()}
    return descriptions


def replot(result: Path, config: Path, check: bool = False) -> dict:
    with PlotRun(result.parent, result_name=result.name) as run:
        return _replot(result, config, check, run)


def _replot(result: Path, config: Path, check: bool, run: PlotRun) -> dict:
    import yaml
    cfg = yaml.safe_load(config.read_text())
    comparison = f"{cfg['case_group']}_vs_{cfg['control_group']}"
    required = ["03_qc/objects/processed_data.rds",
                f"04_statistics/tables/{comparison}_differential_features.csv"]
    pathway_source = "06_pathway/tables/pathway_quantitative_statistics.csv"
    if (result / pathway_source).exists() or (result / "06_pathway/plots/pathway_quantitative_statistics.pdf").exists():
        required.append(pathway_source)
    missing = [name for name in required if not (result / name).is_file()]
    if missing:
        raise FileNotFoundError(f"缺少冻结绘图来源，未启动统计或覆盖结果：{missing}")
    style_source = (ROOT / "R/figure_style.R").read_text()
    import re
    profile = (ROOT / "R" / re.search(r'PFOS_PROFILE_RELATIVE <- "([^"]+)"', style_source).group(1)).resolve()
    profile_sha = digest(profile)
    before = scientific_files(result)
    expected_figures = {str(p.relative_to(result).with_suffix("")) for p in result.glob("*/plots/*.pdf")}
    code_before = {str(p.relative_to(ROOT)): digest(p) for directory in (ROOT / "R", ROOT / "scripts") for p in sorted(directory.glob("*")) if p.suffix in {".R", ".py"}}
    if check:
        return {"status": "PASS", "mode": "plots-only", "scientific_sources": len(before)}
    work = run.stage(f".{result.name}.plots-")
    stage = work / "result"
    stage.mkdir()
    started = datetime.now(timezone.utc).isoformat()
    receipt = {"status": "STARTED", "mode": "plots-only", "started_at": started,
               "source_signature": before, "config_sha256": digest(config),
               "statistics_recomputed": False, "profile": str(profile), "profile_sha256": profile_sha}
    try:
        shutil.copytree(result, stage, dirs_exist_ok=True,
                        ignore=lambda folder, names: [n for n in names if Path(n).suffix.lower() in {".png", ".pdf"}])
        audit = stage / "01_audit"
        audit.mkdir(exist_ok=True)
        (audit / "replot.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
        logs = stage / "logs"
        logs.mkdir(exist_ok=True)
        with (logs / "plots_only.log").open("w") as handle:
            run.command([os.environ.get("RSCRIPT_BIN", "Rscript"), str(ROOT / "scripts/pipeline.R"),
                            "--config", str(config), "--module", "plots-only", "--result-dir", str(stage)],
                           check=True, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
        after = scientific_files(stage)
        changed = [name for name, sha in before.items() if after.get(name) != sha]
        if changed or scientific_files(result) != before or digest(config) != receipt["config_sha256"] or digest(profile) != profile_sha:
            raise RuntimeError(f"冻结科学来源发生变化，拒绝发布：{changed}")
        png = sorted(stage.glob("*/plots/*.png"))
        pdf = sorted(stage.glob("*/plots/*.pdf"))
        if not png or {p.with_suffix("") for p in png} != {p.with_suffix("") for p in pdf}:
            raise RuntimeError("PNG/PDF 图件缺失或未配对")
        if expected_figures and {str(p.relative_to(stage).with_suffix("")) for p in pdf} != expected_figures:
            raise RuntimeError("现役图件类型覆盖不一致，拒绝发布")
        code_after = {str(p.relative_to(ROOT)): digest(p) for directory in (ROOT / "R", ROOT / "scripts") for p in sorted(directory.glob("*")) if p.suffix in {".R", ".py"}}
        if code_after != code_before:
            raise RuntimeError("绘图源码在执行期间发生变化，拒绝发布")
        descriptions = figure_descriptions(stage, comparison)
        from PIL import Image
        for path in png:
            with Image.open(path) as im:
                if im.mode != "RGB" or abs(im.width - 2161) > 1:
                    raise RuntimeError(f"图件导出规格不一致：{path.name} {im.mode} {im.size}")
        lines = ["figure_id\tpng\tpdf\tsource"]
        gallery = ["<!doctype html><html lang='zh'><meta charset='utf-8'><title>PFOS 肝脏代谢组</title>",
                   "<style>body{font:16px sans-serif;max-width:1100px;margin:auto;padding:30px}img{max-width:100%}figure{margin:36px 0}</style>",
                   f"<h1>{html.escape(str(cfg.get('project_name', '代谢组')))} 完整图件</h1><p>本次为冻结结果重绘；QC 身份与公司方向沿用原配置。完整统计保留 BH/FDR。</p>"]
        for path in png:
            parent = path.parent.parent.name
            source = "03_qc/objects/processed_data.rds" if parent == "03_qc" else required[1] if parent == "04_statistics" else required[2]
            rel = path.relative_to(stage)
            pdf_rel = rel.with_suffix(".pdf")
            lines.append(f"{path.stem}\t{rel}\t{pdf_rel}\t{source}")
            meaning, distribution, conclusion = descriptions[path.stem]
            caption = html.escape(f"生物学意义：{meaning} 数据分布特征：{distribution} 最终结论：{conclusion}")
            gallery.append(f"<figure id='{html.escape(path.stem)}'><h2>{html.escape(path.stem)}</h2><img src='{rel}'><figcaption>{caption} <a href='{pdf_rel}'>PDF</a> · <a href='{source}'>完整来源</a></figcaption></figure>")
        (audit / "tables/figure_manifest.tsv").write_text("\n".join(lines) + "\n")
        (stage / "index.html").write_text("\n".join(gallery) + "</html>\n")
        receipt.update(status="COMPLETED", finished_at=datetime.now(timezone.utc).isoformat(),
                       figure_pairs=len(png), source_unchanged=True,
                       added_display_tables=sorted(set(after) - set(before)),
                       code_signature=code_before,
                       output_signature={str(p.relative_to(stage)): digest(p) for p in [*png, *pdf, stage / "index.html", audit / "tables/figure_manifest.tsv"]})
        (audit / "replot.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
        run.publish(stage, result)
        return receipt
    except BaseException as error:
        receipt.update(status="FAILED", error=str(error))
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, default=Path("result"))
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = args.config if args.config.is_absolute() else ROOT / args.config
    result = args.result_dir if args.result_dir.is_absolute() else ROOT / args.result_dir
    receipt = replot(result.resolve(), config.resolve(), args.check)
    print(json.dumps({name: receipt[name] for name in ("status", "mode", "figure_pairs", "source_unchanged", "scientific_sources") if name in receipt}, ensure_ascii=False))


if __name__ == "__main__":
    main()

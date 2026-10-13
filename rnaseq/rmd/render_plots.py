#!/usr/bin/env python3
"""从已完成结果重绘；完整验收后一次替换 result，不重跑统计分析。"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from PIL import Image
from plot_lifecycle import PlotRun


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modules', default='2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    with PlotRun(root) as run:
        result = root / 'result'
        profile = (root / "../../../config/figure_style.json").resolve() if "projects" in root.parts else root / "config/figure_style.json"
        profile_sha = sha(profile)
        chosen = {int(x) for x in args.modules.replace(' ', ',').split(',') if x}
        if not chosen or not chosen <= (set(range(2, 19)) - {12}):
            parser.error('纯绘图模块为02-11、13-18；12没有已完成WGCNA图。')
        figures = {p.relative_to(result).as_posix(): sha(p) for p in result.rglob('*')
                   if p.suffix.lower() in {'.pdf', '.png'} and p.relative_to(result).parts[0][:2].isdigit()
                   and int(p.relative_to(result).parts[0][:2]) in chosen}
        if not figures:
            raise SystemExit('所选模块没有完成图件，不能声明重绘。')
        frozen = {p.relative_to(root).as_posix(): sha(p) for folder in [result, root / '.cache']
                  for p in folder.rglob('*') if p.is_file() and p.suffix.lower() in {'.csv', '.tsv', '.rds', '.rdata'}
                  and not p.name.endswith('_plot_source.tsv')
                  and p.name not in {'figure_manifest.tsv', 'plot_source_sha256.tsv', 'plot_code_sha256.tsv'}}
        with (root / 'rmd/plot_inputs.tsv').open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                path = root / row['path']
                if not path.is_file():
                    raise SystemExit(f'MISSING: {row["path"]}')
                frozen[row['path']] = sha(path)
        frozen['rmd/plot_samples.tsv'] = sha(root / 'rmd/plot_samples.tsv')
        for folder in [root / 'rawdata', root / '5.post_transcriptional_regulation']:
            for path in folder.rglob('*'):
                if path.is_file() and path.suffix.lower() in {'.csv', '.tsv', '.txt', '.json'}:
                    frozen[path.relative_to(root).as_posix()] = sha(path)
        frozen['result_pca/auto_selection_summary.tsv'] = sha(root / 'result_pca/auto_selection_summary.tsv')
        source_code = {p.relative_to(root).as_posix(): sha(p) for p in (root / 'rmd').rglob('*')
                       if p.is_file() and p.suffix in {'.R', '.Rmd', '.py', '.tsv'}}
        stage = run.stage('.plot-stage-')
        shutil.copytree(result, stage / 'result', dirs_exist_ok=True)
        staged = stage / 'result'
        for relative in figures:
            (staged / relative).unlink()
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', MAGICK_THREAD_LIMIT='1')
        try:
            run.command(['Rscript', str(root / 'rmd/plots_only.R'), str(root), str(staged), args.modules], check=True, env=env)
            missing = [relative for relative in figures if not (staged / relative).is_file()]
            if missing:
                raise RuntimeError('原有图件未全部重绘: ' + ', '.join(missing))
            for relative in figures:
                if relative.endswith('.png'):
                    path = staged / relative
                    with Image.open(path) as image:
                        if image.mode != 'RGB':
                            raise RuntimeError(f'PNG未导出RGB: {relative}')
                        if abs(image.info.get('dpi', (0, 0))[0] - 300) > 0.1:
                            image.save(path, dpi=(300, 300))
            changed = [relative for relative, digest in frozen.items() if sha(root / relative) != digest]
            changed += [relative for relative, digest in frozen.items()
                        if relative.startswith('result/') and sha(stage / relative) != digest]
            if changed:
                raise RuntimeError('冻结数值源发生变化: ' + ', '.join(changed))
            if sha(profile) != profile_sha:
                raise RuntimeError("绘图期间project profile改变")
            changed_code = [relative for relative, digest in source_code.items() if sha(root / relative) != digest]
            if changed_code:
                raise RuntimeError('绘图过程中源码发生变化: ' + ', '.join(changed_code))
            audit = staged / '00_audit/tables'
            audit.mkdir(parents=True, exist_ok=True)
            with (audit / 'plot_source_sha256.tsv').open('w') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow(['path', 'sha256'])
                writer.writerows(sorted(frozen.items()))
            with (audit / 'plot_code_sha256.tsv').open('w') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow(['path', 'sha256'])
                writer.writerows(sorted(source_code.items()))
            invalid_prefixes = ('15_31_', '15_32_', '15_32b_', '15_32c_', '15_36_', '15_37_',
                                '15_38_', '15_39_', '15_54_', '15_55_', '15_56_', '17_13_', '17_14_',
                                '17_21_', '17_22_', '17_23_')
            invalid_files = []
            for path in sorted(staged.rglob('*')):
                if path.is_file() and (path.name.startswith(invalid_prefixes) or path.name == '15_splicing_gene_evidence.csv') and path.suffix.lower() in {'.pdf', '.png', '.csv', '.tsv'}:
                    invalid_files.append({
                        'path': path.relative_to(staged).as_posix(), 'path_base': 'result', 'sha256': sha(path),
                        'scientific_status': 'INVALID_FROZEN_ADAPTER',
                        'role': 'historical_figure' if path.suffix.lower() in {'.pdf', '.png'} else
                                'display_sidecar' if path.name.endswith('_plot_source.tsv') else 'frozen_derived_table',
                        'excluded_from_scientific_conclusions': True,
                    })
            for name in ('gene_method_membership.csv', 'as_gene_evidence_table.csv', 'top_as_candidates.csv'):
                path = root / '5.post_transcriptional_regulation/integrated' / name
                if path.is_file():
                    invalid_files.append({'path': path.relative_to(root).as_posix(), 'path_base': 'analysis_module',
                        'sha256': sha(path), 'scientific_status': 'INVALID_FROZEN_ADAPTER',
                        'role': 'frozen_adapter_source', 'excluded_from_scientific_conclusions': True})
            adapter_audit = {
                'status': 'INVALID_FROZEN_ADAPTER',
                'reason': 'Legacy adapter collapsed unreported methods to FALSE/0; DER all FALSE and DE log2FC all NA in the audited PFOS frozen evidence.',
                'frozen_source_action': 'original_bytes_preserved; no statistical refit; historical display only',
                'replacement': '../../integration/README.md',
                'files': invalid_files,
            }
            (audit / 'frozen_adapter_scientific_status.json').write_text(json.dumps(adapter_audit, ensure_ascii=False, indent=2) + '\n')
            receipt = {
                'completed_at': datetime.now(timezone.utc).isoformat(),
                'status': 'COMPLETED', 'mode': 'plots_only', 'modules': sorted(chosen),
                'contrast': 'PFOS - WT', 'upstream_statistics_rerun': False, 'figure_profile': os.path.relpath(profile,root), 'figure_profile_sha256': profile_sha,
                'display_derivatives': ['PCA', 'correlation', 'heatmap clustering', 'network layout', 'descriptive log-log fit'],
                'source_samples': 'rmd/plot_samples.tsv', 'source_inputs': 'rmd/plot_inputs.tsv',
                'original_figures_preserved': len(figures), 'frozen_files_unchanged': len(frozen),
                'scientific_status_scope': 'FROZEN_SOURCE_UNCHANGED is byte protection, not scientific validation',
                'invalid_frozen_adapter_audit': '00_audit/tables/frozen_adapter_scientific_status.json',
                'invalid_frozen_adapter_figures': [x['path'] for x in invalid_files if x['role'] == 'historical_figure' and x['path'].endswith('.pdf')],
                'invalid_frozen_adapter_sidecars': [x['path'] for x in invalid_files if x['role'] == 'display_sidecar'],
                'invalid_frozen_adapter_tables': [x['path'] for x in invalid_files if x['role'] == 'frozen_derived_table'],
                'invalid_frozen_adapter_inputs': [x['path'] for x in invalid_files if x['role'] == 'frozen_adapter_source'],
            }
            (audit / 'plot_receipt.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')
            run.command(['Rscript', str(root / 'rmd/generate_result_report.R')], check=True,
                           env=dict(env, RNASEQ_ROOT_DIR=str(root), RNASEQ_RESULT_DIR=str(staged)))
            run.publish(staged, result)
            print(json.dumps(receipt, ensure_ascii=False))
        except BaseException:
            print('重绘未发布；必要日志由绘图生命周期收尾保存。')
            raise


if __name__ == '__main__':
    main()

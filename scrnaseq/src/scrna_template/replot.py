"""从已完成运行的签名 TSV 重绘通用图谱；不运行统计。"""
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import tempfile
from .plot_lifecycle import PlotRun

import pandas as pd

from .plotting import categorical_embedding, stacked_composition, sample_palettes, CLUSTER_COLORS, CELL_TYPE_COLORS
from .runtime import atomic_json, atomic_text, sha256


def replot(config, source_run, check=False):
    with PlotRun(config.root) as run:
        return _replot(config, source_run, check, run)


def _replot(config, source_run, check, run):
    source = config.resolve(source_run)
    if json.loads((source / '00_audit/run_result.json').read_text())['status'] != 'PASS':
        raise RuntimeError('来源运行尚未通过执行审计')
    signatures = pd.read_csv(source / '00_audit/output_sha256.tsv', sep='\t').set_index('path').sha256
    stems = sorted(path.stem for path in (source / '08_figures').glob('*.png'))
    if not stems:
        raise RuntimeError('来源运行没有实际图件')
    before = {}
    for stem in stems:
        if not (stem.startswith('umap_') or stem == 'captured_composition'):
            raise RuntimeError(f'通用模板没有该专题图件的纯绘图函数：{stem}')
        relative = f'07_source_data/{stem}.tsv'
        path = source / relative
        if relative not in signatures.index or not path.is_file() or sha256(path) != signatures[relative]:
            raise RuntimeError(f'绘图来源缺失或与原签名不一致：{relative}')
        before[relative] = signatures[relative]
    if check:
        return {'status': 'PASS', 'figure_pairs': len(stems), 'statistics_recomputed': False}
    result = config.root / 'result'
    work = run.stage('.result.plots-')
    stage = work / 'result'
    stage.mkdir()
    rows = []
    conditions, libraries = sample_palettes(config)
    try:
        for stem in stems:
            topic = '04_descriptive' if stem == 'captured_composition' else '03_annotation' if stem == 'umap_cell_type' else '02_atlas'
            plots, tables = stage / topic / 'plots', stage / topic / 'tables'
            tables.mkdir(parents=True, exist_ok=True)
            frame = pd.read_csv(source / f'07_source_data/{stem}.tsv', sep='\t')
            if stem == 'captured_composition':
                stacked_composition(frame, plots / stem)
            else:
                column = 'cell_type' if 'cell_type' in frame else 'label'
                palette = {'umap_condition': conditions, 'umap_library': libraries, 'umap_leiden': CLUSTER_COLORS, 'umap_cell_type': CELL_TYPE_COLORS}[stem]
                categorical_embedding(frame[['UMAP1', 'UMAP2']].to_numpy(), frame[column], stem, plots / stem, palette)
            shutil.copy2(source / f'07_source_data/{stem}.tsv', tables / f'{stem}.tsv')
            rows.append({'figure_id': stem, 'png': f'{topic}/plots/{stem}.png', 'pdf': f'{topic}/plots/{stem}.pdf', 'source': f'{topic}/tables/{stem}.tsv'})
        if any(sha256(source / name) != value for name, value in before.items()):
            raise RuntimeError('来源在绘图期间发生变化，拒绝发布')
        audit = stage / '00_audit'
        (audit / 'tables').mkdir(parents=True)
        pd.DataFrame(rows).to_csv(audit / 'tables/figure_manifest.tsv', sep='\t', index=False)
        receipt = {'status': 'COMPLETED', 'mode': 'plots-only', 'finished_at': datetime.now(timezone.utc).isoformat(), 'source_signature': before, 'source_unchanged': True, 'figure_pairs': len(rows), 'statistics_recomputed': False}
        atomic_json(audit / 'replot.json', receipt)
        atomic_text(stage / 'index.html', '<!doctype html><html lang="zh"><meta charset="utf-8"><h1>单细胞完整重绘图谱</h1><p>图谱与捕获组成保持描述性；样本级推断见原运行。</p>' + ''.join(f'<figure><img style="max-width:100%" src="{r["png"]}"><figcaption>{r["figure_id"]} <a href="{r["pdf"]}">PDF</a> · <a href="{r["source"]}">数值源</a></figcaption></figure>' for r in rows) + '</html>')
        run.publish(stage, result)
        return receipt
    except BaseException as error:
        raise

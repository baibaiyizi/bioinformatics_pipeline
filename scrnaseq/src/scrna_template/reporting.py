from __future__ import annotations
import html
import pandas as pd
from .runtime import atomic_text, atomic_json

def run_report(context):
    samples = context.config.samples()
    output = context.path('09_reports', 'scrna_report_zh.html')
    sections = [f'<h1>{html.escape(str(context.config.data["project"]["name"]))} 单细胞分析</h1>', '<p>细胞用于图谱和注释，推断单位为独立 biological_sample_id；每个混池只计一个样本。同一样本多库先聚合原始整数 counts。</p>']
    sections.append('<p>比较方向：' + html.escape(context.config.data['statistics']['case']) + ' − ' + html.escape(context.config.data['statistics']['control']) + '。默认原始 P 值筛选是探索性结果；BH/FDR 同时保留。</p>')
    sections.extend(['<h2>样本登记</h2>', samples[['library_id','biological_sample_id','unit_type','condition','sample_batch']].to_html(index=False, escape=True)])
    for title, relative in [('QC', '02_qc/qc_summary.tsv'), ('推断资格', '00_audit/analysis_eligibility.tsv'), ('实际推断状态', '06_optional/inference_status.tsv'), ('注释证据', '04_annotation/annotation_evidence.tsv'), ('捕获组成', '05_descriptive/captured_composition.tsv')]:
        path = context.run_dir / relative
        sections.append(f'<h2>{title}</h2>')
        if path.is_file():
            sections.append(pd.read_csv(path, sep='\t', keep_default_na=False).to_html(index=False, escape=True))
        else:
            sections.append('<p>MISSING / not_run</p>')
    for path in sorted(context.path('06_optional').glob('pseudobulk_DE_*.tsv')):
        result = pd.read_csv(path, sep='\t')
        sections.append(f'<h2>{html.escape(str(result.cell_type.iloc[0]))} 的样本级检验</h2><p>检验 {len(result)} 个基因；原始 P 探索性阈值通过 {int(result.passes_raw_p_exploratory.sum())} 个；BH/FDR 阈值通过 {int(result.passes_BH_FDR.sum())} 个。</p>')
        sections.append(f'<a href="../{html.escape(str(path.relative_to(context.run_dir)))}">完整统计表</a>')
    figures = sorted(context.path('08_figures').glob('*.png'))
    sections.append('<h2>图谱与描述</h2>')
    if figures:
        for path in figures:
            sections.append(f'<figure><img style="max-width:800px" src="../08_figures/{path.name}"><figcaption>{html.escape(path.stem)}；描述性图形，不替代样本级推断。</figcaption></figure>')
    else:
        sections.append('<p>本次为仅表报告，图形未生成。</p>')
    sections.append('<h2>能力边界</h2><p>marker score 为注释证据；参考映射、去卷积、条件组成推断、通讯、velocity、GRN 均不在通用核心中执行。原始 FASTQ 路径必须经 MD5/完整成对/稳定性审计。合成数据验证不等于真实数据复现。</p>')
    atomic_text(output, '<!doctype html><html lang="zh"><meta charset="UTF-8"><title>scRNA 报告</title><body>' + '\n'.join(sections) + '</body></html>')
    return output

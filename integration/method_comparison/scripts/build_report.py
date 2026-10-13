#!/usr/bin/env python3
"""核验现役联合结果后更新模块 README；图表直接引用原始结果。"""
import csv
import hashlib
import json
import os
from pathlib import Path


def read(path, delimiter=','):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(os.environ['METHOD_COMPARISON_ROOT']).resolve()
    synthesis = root / 'result/synthesis'
    synthesis.mkdir(parents=True, exist_ok=True)
    receipt_path = synthesis / 'report_receipt.json'
    receipt = {'state': 'running', 'renderer_sha256': sha(Path(__file__))}
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')
    try:
        inputs = []
        for stage in ['describe', 'differential']:
            current = root / 'result' / stage
            status = read(current / 'run_status.tsv', '\t')
            if len(status) != 1 or status[0]['State'] != 'completed':
                raise ValueError('not_run: 联合阶段未完成 ' + stage)
            signature_path = current / 'input_signature.tsv'
            for row in read(signature_path, '\t'):
                path = Path(row['Path'])
                if not path.is_file() or path.stat().st_size != int(row['Size']):
                    raise ValueError('STALE_OUTPUTS: ' + str(path))
                digest = hashlib.md5()
                with path.open('rb') as handle:
                    for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
                        digest.update(chunk)
                if digest.hexdigest() != row['MD5']:
                    raise ValueError('STALE_OUTPUTS: ' + str(path))
            inputs += [current / 'run_status.tsv', signature_path]
        tables = [
            ('方法捕获注释层', 'differential/04_differential_capture/tables/04_annotation/capture_results_method2_vs_method1.csv'),
            ('配对设计', 'differential/04_differential_capture/tables/capture_contrast_manifest.csv'),
            ('外参过程 QC', 'describe/03_spikein_qc/tables/spikein_qc_summary.csv'),
        ]
        annotation = read(root / 'result' / tables[0][1])
        contrast = read(root / 'result' / tables[1][1])[0]
        spike = read(root / 'result' / tables[2][1])[0]
        method2 = sum(x['direction'] == 'up' for x in annotation)
        method1 = sum(x['direction'] == 'down' for x in annotation)
        report = [
            '# 人血清 smallRNA 方法比较', '',
            '## 要回答的问题', '',
            '在同一供体的血清材料中，方法 1 与方法 2 的捕获谱、检出能力和技术偏倚如何变化？方法 3–6 提供单库描述；供体 3 的外参配对单独检查流程表现。', '',
            '## 运行与输入', '',
            '入口不依赖启动目录。先核验单组学已完成的 QC、输入及源码签名，再运行本模块冻结的代码；序列层重建的 180 个类别×文库计数必须与单组学标准计数逐值一致。', '',
            '```bash', 'bash run.sh --check', 'bash run.sh --all',
            'bash run.sh --stage describe', 'bash run.sh --stage differential',
            'bash run.sh --stage report', '```', '',
            '源数据由 [config/project.sh](config/project.sh) 指定；[单组学 QC 与分类](../smallrnaseq/README.md) 提供标准输入。报告只引用本模块结果，不另建图表副本。', '',
            '## 当前结果', '',
            '描述与配对差异阶段均已实际执行；运行状态和输入签名随结果保留。',
            f'注释层共检验 {len(annotation)} 个特征；方法 2 捕获较高 {method2} 个、方法 1 较高 {method1} 个。正 logFC 统一表示方法 2−方法 1。',
            f"显著性指标为 {contrast['SignificanceMetric']}，阈值 {float(contrast['Cutoff']):.3g}；绝对 log2FC 阈值 {float(contrast['AbsoluteLog2FCCutoff']):.3f}。",
            f"供体 3 加 60 fmol cel-miR-39 的文库观察到精确匹配 {spike['ExactReads']} reads（{float(spike['ExactRPM']):.3f} RPM），不能换算提取回收率。", '',
            '## 如何解释', '',
            '方法比较来自供体 2，方法 1/2 只有三个匹配技术分装，采用 ~PairID + Method；DonorN=1。捕获差异限于本供体的技术比较，不代表人群生物学差异，也不支持某方法普遍优于其他方法。方法 3–6 各一库；供体 3 两库不进入方法排名。', '',
            '深度检出曲线沿用 Poisson 期望模型，只在每库实际观察深度以内计算，不能称为精确无放回抽样。供体 3 序列层的 200 次无放回稀释只反映技术抽样波动；注释层分数计数不取整稀释，其等深度抽样保持 not_run。', '',
            '## 完整结果与来源', '',
            '- [描述图册](result/describe/index.html)',
            '- [配对差异图册](result/differential/index.html)',
            '- [描述运行状态](result/describe/run_status.tsv)',
            '- [差异运行状态](result/differential/run_status.tsv)',
            '- [源码快照](SOURCE_SNAPSHOT.tsv)',
            '- [模板来源](template_snapshot.tsv)',
        ]
        for title, relative in tables:
            path = root / 'result' / relative
            inputs.append(path)
            report.append(f'- [{title}](result/{relative})')
        report += ['', '历史输出、旧源码与日志保存在本模块内部 history，保留原执行身份；当前结论以本轮重新执行的结果为准。', '', '分析顺序：01_setup → 02_method_similarity → 03_spikein_qc → 04_differential_capture → 05_capture_bias_diagnostics。外参类别图在 03_spikein_qc/profiles 内按 01–05 连续编号。', '', '## 主要图及解释', '']
        figures = [
            ('方法间 RNA 组成', 'describe/02_method_similarity/plots/02_rna_class_composition/rna_class_composition_by_method', '按文库比较 RNA 类别占比，判断不同方法捕获的 RNA 谱。占比依赖总 reads 分母，不能直接视为绝对 RNA 产量。'),
            ('观察到的检出广度', 'describe/02_method_similarity/plots/03_detection_sensitivity/observed_detection_counts', '展示当前深度下实际检出的特征数。方法 3–6 各只有一库，检出差异同时受测序深度影响，不能据此作普遍方法排名。'),
            ('技术分装的一致性', 'describe/02_method_similarity/plots/06_reproducibility/technical_replicate_agreement', '检查同一方法内部技术分装的一致性；相关性用于评估技术稳定性，不能增加供体数。'),
            ('外参观察信号', 'describe/03_spikein_qc/plots/spikein_detection_overview', '供体 3 的 cel-miR-39 精确匹配与近似匹配分别报告；没有已知全过程回收标定，不能从 reads 推算提取回收率。'),
            ('供体 3 深度敏感性', 'describe/03_spikein_qc/plots/paired_depth_sensitivity', '供体 3 两库单独比较观察深度和等深度抽样敏感性，反映该次实验的技术抽样波动，不参与供体 2 的方法排名。'),
            ('注释层差异捕获', 'differential/04_differential_capture/plots/04_annotation/capture_volcano_method2_vs_method1', '横轴向右表示方法 2 捕获较高，向左表示方法 1 较高；纵轴为显著性。结论限于该供体的技术捕获差异，逐项数值见注释层结果表。'),
            ('类别层差异捕获', 'differential/04_differential_capture/plots/02_category/capture_volcano_method2_vs_method1', '把序列捕获差异汇总到 RNA 大类，分清类别整体变化与少数注释特征变化；同样采用方法 2−方法 1 方向。'),
            ('配对捕获方向的一致性', 'differential/05_capture_bias_diagnostics/plots/02_pair_consistency/capture_pair_logfc_consistency', '每个面板比较两个 PairID 的序列 log2FC。第一、三象限表示方向一致，第二、四象限提示方向不一致；每个点是序列，不是独立供体。'),
        ]
        for number, (title, relative, note) in enumerate(figures, 1):
            for suffix in ['.png', '.pdf']:
                path = root / 'result' / (relative + suffix)
                if not path.is_file():
                    raise ValueError('MISSING: ' + str(path))
                inputs.append(path)
            report += [f'### {number}. {title}', '', note, '',
                       f'![{title}](result/{relative}.png)', '',
                       f'[下载 PDF](result/{relative}.pdf)', '']
        text = '\n'.join(report) + '\n'
        (root / 'README.md').write_text(text)
        import markdown
        html = markdown.markdown(text, extensions=['tables'])
        (synthesis / 'index.html').write_text('<!doctype html><html lang="zh"><head><meta charset="utf-8"><base href="../../"><title>人血清方法比较</title><style>body{font:17px sans-serif;max-width:1100px;margin:auto;padding:30px;line-height:1.7;color:#263445}img{max-width:100%;height:auto}a{color:#2166ac}</style></head><body>' + html + '</body></html>')
        outputs = [root / 'README.md', synthesis / 'index.html']
        receipt.update(state='completed',
                       inputs=[{'path': str(p.relative_to(root)), 'sha256': sha(p)} for p in inputs],
                       outputs={str(p.relative_to(root)): sha(p) for p in outputs})
        print('[REPORT]', root / 'README.md')
    except Exception as error:
        receipt.update(state='failed', detail=str(error))
        raise
    finally:
        receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')


if __name__ == '__main__':
    main()

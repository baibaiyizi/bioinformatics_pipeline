#!/usr/bin/env python3
"""Read-only sample contract; --record explicitly writes execution provenance."""
import argparse
import csv
import hashlib
import json
import os
import re
from pathlib import Path


def audit(root, prefix, require_fastq=False):
    manifest = Path(os.environ[prefix + '_SAMPLES'])
    with manifest.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = ['sample', 'biological_sample_id', 'group', 'stage']
        if os.environ.get('ATAC_INPUT_MODE','bam') != 'counts': required += ['raw_fq1','raw_fq2']
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('样本表缺失或重复列名')
        if not set(required).issubset(reader.fieldnames):
            raise ValueError('缺少样本表列: ' + ', '.join(sorted(set(required)-set(reader.fieldnames))))
        rows = list(reader)
    if not rows:
        raise ValueError('样本表为空')
    reference = os.environ.get(prefix + '_REFERENCE_GROUP', '')
    comparison = os.environ.get(prefix + '_COMPARISON_GROUP', '')
    if not reference or not comparison or reference == comparison:
        raise ValueError('必须显式提供不同的 reference 与 comparison 组')
    if os.environ.get(prefix + '_DESIGN', 'independent_two_group') != 'independent_two_group':
        raise ValueError('仅支持独立两组；配对、重复测量和协变量设计未实现')
    if os.environ.get(prefix+'_ASSEMBLY') != 'GRCm39':
        raise ValueError('本模板仅支持小鼠 GRCm39')
    if prefix == 'ATAC' and os.environ.get('ATAC_LIBRARY_TYPE') != 'paired_end':
        raise ValueError('ATAC 上游仅支持 paired_end')
    seen_samples, seen_bios = set(), set()
    for number, row in enumerate(rows, 2):
        if None in row:
            raise ValueError(f'第 {number} 行多余字段')
        for column in required:
            if not row.get(column) or not row[column].strip():
                raise ValueError(f'第 {number} 行 {column} 为空；不会折叠 TSV 空字段')
            row[column] = row[column].strip()
        for column in ['sample', 'biological_sample_id', 'group']:
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', row[column]):
                raise ValueError(f'第 {number} 行 {column} 仅允许字母、数字、下划线、点、短横线')
        if row['sample'] in seen_samples:
            raise ValueError('重复 sample: ' + row['sample'])
        if row['biological_sample_id'] in seen_bios:
            raise ValueError('重复 biological_sample_id：先在实验单位层面合并技术文库，不能增加 n')
        seen_samples.add(row['sample']); seen_bios.add(row['biological_sample_id'])
        if row['biological_sample_id'].lower() in {'unknown', 'na', 'missing'}:
            raise ValueError('生物学身份未知：不满足两组推断资格')
        for column in ['raw_fq1', 'raw_fq2']:
            if not row.get(column): continue
            p = Path(row[column]); p = p if p.is_absolute() else root / p
            if require_fastq and not p.is_file():
                raise ValueError('缺少 FASTQ: ' + str(p))
        for column in ['pair_id', 'subject_id', 'batch']:
            if row.get(column, '').strip():
                raise ValueError(f'{column} 设计尚未实现；不能忽略样本表中的该字段')
    if {row['group'] for row in rows} != {reference, comparison}:
        raise ValueError('样本表组必须恰好等于显式 reference 和 comparison')
    sizes = {group: sum(row['group'] == group for row in rows) for group in [reference, comparison]}
    if min(sizes.values()) < 2:
        raise ValueError('每组至少需要 2 个独立生物学样本；不等价于功效充分')
    return {'status': 'PASS_INPUT_CONTRACT', 'contrast': comparison+' - '+reference,
            'biological_n': sizes, 'sample_count': len(rows),
            'significance': 'exploratory_raw_pvalue' if os.environ.get(prefix+'_SIGNIFICANCE_MODE', 'pvalue') == 'pvalue' else 'adjusted',
            'analysis_execution': 'not_run', 'samples': rows}


def signature(root, prefix, report):
    declared = set(re.findall(r'^(?:export )?([A-Za-z_][A-Za-z0-9_]*)=', (root/'config/project.env').read_text(), re.M))
    inputs = {'contract': report, 'configuration': {k: v for k, v in os.environ.items()
              if (k.startswith(prefix + '_') or k in declared) and k not in {prefix+'_LOG_DIR',prefix+'_ROOT_DIR',prefix+'_RMD_DIR'}}}
    paths = [row[k] for row in report['samples'] for k in ['raw_fq1','raw_fq2'] if row.get(k)]
    paths += [v for k,v in os.environ.items() if v and (k.startswith(prefix+'_') or k in {'GENOME_FASTA','GENE_GTF','BLACKLIST_BED','TSS_BED'}) and k.endswith(('_FILE','_FASTA','_GTF','_BED'))]
    freshness = {}
    for value in paths:
        path = Path(value); path = path if path.is_absolute() else root/path
        if path.is_file():
            info = path.stat(); freshness[str(path)] = [info.st_size, info.st_mtime_ns]
        else:
            freshness[str(path)] = 'MISSING'
    inputs['input_size_mtime_not_content_hash'] = freshness
    inputs['source_sha256'] = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted([*root.glob('*.sh'), *[p for folder in ['rmd','scripts','config'] for p in (root/folder).rglob('*')]])
        if p.is_file() and p.suffix in {'.R','.Rmd','.sh','.py','.env','.tsv'}}
    return hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--prefix', choices=['WGBS','ATAC'], required=True)
    parser.add_argument('--require-fastq', action='store_true')
    parser.add_argument('--record', choices=['upstream','downstream'])
    parser.add_argument('--verify', choices=['upstream','downstream'])
    args = parser.parse_args()
    try:
        report = audit(args.root, args.prefix, args.require_fastq)
        if args.verify:
            stamp = args.root / 'metadata' / (args.verify+'_configuration.json')
            if not stamp.is_file():
                raise ValueError('UNVERIFIED：没有该输入/配置的来源记录')
            if json.loads(stamp.read_text()).get('sha256') != signature(args.root, args.prefix, report):
                raise ValueError('STALE：当前配置/样本/脚本/输入文件属性与来源记录不同')
        if args.record:
            stamp = args.root / 'metadata' / (args.record+'_configuration.json')
            fingerprint = signature(args.root, args.prefix, report)
            if stamp.exists() and json.loads(stamp.read_text()).get('sha256') != fingerprint:
                raise ValueError('STALE：样本/方向/配置/脚本变更。请使用新的项目输出目录；不自动复用旧结果')
            if not stamp.exists():
                outputs = ['result', '.cache'] if args.record == 'downstream' else (
                    ['02bismark_results_GRCm39','04methykit'] if args.prefix == 'WGBS' else ['02mapping','03peaks','04counts'])
                if any((args.root/f).exists() and any((args.root/f).rglob('*')) for f in outputs):
                    raise ValueError('UNVERIFIED：已有结果缺少模板运行签名；请使用新的项目输出目录')
            stamp.parent.mkdir(parents=True, exist_ok=True)
            stamp.write_text(json.dumps({'sha256': fingerprint, 'input_contract': report,
                           'note':'配置资格记录不代表模块完成；执行状态另行记录'}, ensure_ascii=False, indent=2)+'\n')
        print(json.dumps({k:v for k,v in report.items() if k != 'samples'}, ensure_ascii=False))
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as exc:
        parser.exit(2, 'ERROR: '+str(exc)+'\n')

if __name__ == '__main__':
    main()

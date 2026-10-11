#!/usr/bin/env python3
"""Strict local sample/design contract. Stdlib only; preflight never writes."""
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys

ROOT = Path(os.environ['RNASEQ_ROOT_DIR'])

def path(value):
    p = Path(value)
    return p if p.is_absolute() else ROOT / p

def load():
    if os.environ.get('RNASEQ_DESIGN') != 'independent_two_group':
        raise ValueError('仅支持 independent_two_group；配对、重复测量和协变量模型未实现')
    groups = [os.environ.get('RNASEQ_GROUP1', ''), os.environ.get('RNASEQ_GROUP2', '')]
    if len(set(groups)) != 2 or any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', x) for x in groups):
        raise ValueError('需显式给出两个不同的合法组名 RNASEQ_GROUP1/2')
    if os.environ.get('RNASEQ_INPUT_MODE') not in {'counts', 'fastq'}:
        raise ValueError('RNASEQ_INPUT_MODE 必须是 counts 或 fastq')
    if os.environ.get('RNASEQ_SPECIES') != 'mouse':
        raise ValueError('当前模块已验证的注释实现仅支持 mouse')
    p = path(os.environ['RNASEQ_SAMPLES_TSV'])
    with p.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'sample', 'group', 'biological_id', 'stage', 'eligible', 'raw_fq1', 'raw_fq2'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('样本表缺少列: ' + ','.join(sorted(required - set(reader.fieldnames or []))))
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError('样本表存在重复列名')
        rows = list(reader)
    seen = set(); active = []; bios = set(); reads = set()
    for row in rows:
        if None in row or any(v is None for v in row.values()):
            raise ValueError('样本表存在列数不匹配行')
        row = {k: v.strip() for k, v in row.items()}
        sample = row['sample']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample) or sample in seen:
            raise ValueError('sample 为空、格式无效或重复: ' + sample)
        seen.add(sample)
        if row['eligible'].lower() not in {'true', 'false'}:
            raise ValueError('eligible 必须是 true/false: ' + sample)
        if row['group'] not in groups:
            raise ValueError('未知组别: ' + row['group'])
        if row['eligible'].lower() == 'false':
            continue
        bio = row['biological_id']
        if bio.lower() in {'', 'unknown', 'na', 'missing'} or bio in bios:
            raise ValueError('生物学单位未知或重复；本模板要求每个独立动物/池一列: ' + sample)
        bios.add(bio)
        for key in ('pair_id', 'subject_id', 'batch'):
            if row.get(key, ''):
                raise ValueError('当前两组设计不接受非空 ' + key)
        if os.environ['RNASEQ_INPUT_MODE'] == 'fastq':
            for key in ('raw_fq1', 'raw_fq2'):
                if not row[key]:
                    raise ValueError('FASTQ 模式缺少 ' + key + ': ' + sample)
                fq = str(path(row[key]))
                if '|' in fq or '\n' in fq or fq in reads:
                    raise ValueError('FASTQ 路径重复或含保留字符: ' + fq)
                reads.add(fq)
        active.append(row)
    if any(sum(r['group'] == g for r in active) < 2 for g in groups):
        raise ValueError('每组至少需两个预先合格的独立生物学样本')
    return sorted(active, key=lambda r: (groups.index(r['group']), r['sample']))

def fingerprint(rows):
    keys = ['RNASEQ_INPUT_MODE', 'RNASEQ_DESIGN', 'RNASEQ_GROUP1', 'RNASEQ_GROUP2', 'RNASEQ_SPECIES', 'RNASEQ_LIBRARY_STRANDEDNESS', 'RNASEQ_SEED', 'RNASEQ_DE_USE_PADJ', 'RNASEQ_DE_PVALUE_CUTOFF', 'RNASEQ_DE_LOG2FC_CUTOFF', 'GTF_FILE', 'HISAT2_INDEX', 'GENOME_FASTA']
    config = {k: os.environ.get(k, '') for k in keys}

    sources = {}
    if os.environ['RNASEQ_INPUT_MODE'] == 'counts':
        cp = path(os.environ['RNASEQ_COUNTS_FILE'])
        if cp.is_file():
            sources['counts_sha256'] = hashlib.sha256(cp.read_bytes()).hexdigest()
    else:
        sources['fastq_metadata'] = []
        for row in rows:
            for key in ('raw_fq1', 'raw_fq2'):
                fq = path(row[key])
                stat = fq.stat() if fq.is_file() else None
                sources['fastq_metadata'].append([str(fq), stat.st_size if stat else None, stat.st_mtime_ns if stat else None])
    return hashlib.sha256(json.dumps([rows, config, sources], sort_keys=True).encode()).hexdigest()

def main():
    command = sys.argv[1] if len(sys.argv) > 1 else 'preflight'
    rows = load()
    if command == 'rows':
        for r in rows:
            print('|'.join([r['sample'], r['group'], str(path(r['raw_fq1'])), str(path(r['raw_fq2']))]))
        return
    if command in {'preflight', 'validate', 'prepare'}:
        if os.environ['RNASEQ_INPUT_MODE'] == 'fastq':
            required = [path(r[k]) for r in rows for k in ('raw_fq1', 'raw_fq2')]
            required += [path(os.environ['GTF_FILE'])]
            index = path(os.environ['HISAT2_INDEX'])
            if not list(index.parent.glob(index.name + '.1.ht2*')):
                raise ValueError('MISSING HISAT2_INDEX: ' + str(index))
        else:
            required = [path(os.environ['RNASEQ_COUNTS_FILE'])]
        missing = [str(p) for p in required if not p.is_file() or p.stat().st_size == 0]
        if missing:
            raise ValueError('MISSING: ' + ', '.join(missing))
        if command == 'preflight':
            print('PASS: 配置/身份/独立两组/输入存在性；未执行分析，未创建输出')
            return
    selected = ROOT / '4.DE_analysis/samples.txt'
    if selected.exists():
        with selected.open(newline='') as handle:
            saved = list(csv.reader(handle, delimiter='\t'))
        expected = sorted([[r['group'], r['sample']] for r in rows])
        if sorted(saved) != expected:
            raise ValueError('STALE: samples.txt 身份与全部预先合格样本不一致')
    signature = fingerprint(rows)
    stamp = ROOT / '.cache' / 'input_contract.sha256'
    if not stamp.exists() and any((ROOT / d).exists() for d in ('3.Merge_result', '4.DE_analysis', 'result_pca')):
        raise ValueError('UNVERIFIED: 已有输出没有当前输入签名；请在新复制目录开始，不能自动承认旧缓存')
    if stamp.exists() and stamp.read_text().strip() != signature:
        raise ValueError('STALE: 输入/配置变更。请使用全新输出目录重新运行；旧结果不能证明当前配置已执行')
    if command == 'prepare':
        stamp.parent.mkdir(parents=True, exist_ok=True)
        stamp.write_text(signature + '\n')
        dest = ROOT / '4.DE_analysis'
        dest.mkdir(exist_ok=True)
        (dest / 'samples.txt').write_text(''.join(r['group'] + '\t' + r['sample'] + '\n' for r in rows))
    elif command != 'validate':
        raise ValueError('未知命令: ' + command)

if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        sys.exit(1)

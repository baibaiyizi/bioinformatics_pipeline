#!/usr/bin/env python3
"""只读预检；不创建目录、不载入会写缓存的 R 上下文。"""
import argparse
import csv
import os
from pathlib import Path
import shutil
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True)
    parser.add_argument('--modules', default='')
    args = parser.parse_args()
    root = Path(os.environ['PANDORA_PROJECT_ROOT'])
    species = os.environ['PANDORA_SPECIES']
    endo = os.environ['RUN_ENDOSIRNA'] == '1'
    layout = os.environ['SPORTS_LAYOUT']
    expected_layout = 'nested' if species == 'Mus' else 'flat'
    if layout != expected_layout:
        raise ValueError(f'{species} 必须显式使用 SPORTS_LAYOUT={expected_layout}')
    with Path(os.environ['PANDORA_SAMPLES_FILE']).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'sample', 'group', 'raw_fastq_gz', 'species'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('样本表缺少字段: ' + ','.join(sorted(required - set(reader.fieldnames or []))))
        rows = list(reader)
    samples = [row['sample'].strip() for row in rows]
    if not rows or any(not s or '/' in s or s in {'.', '..'} for s in samples) or len(set(samples)) != len(samples):
        raise ValueError('样本表为空、样本名非法或重复')
    expected_species = {'Mus': {'Mus musculus', 'Mus_musculus'}, 'Homo': {'Homo sapiens', 'Homo_sapiens'}}[species]
    if any(row['species'].strip() not in expected_species for row in rows):
        raise ValueError('样本表 species 与 PROJECT_SPECIES 不一致')
    raw_missing = []
    required_files = []
    for row in rows:
        sample = row['sample'].strip()
        raw = Path(row['raw_fastq_gz'])
        raw = raw if raw.is_absolute() else root / raw
        if not row['raw_fastq_gz'].strip() or not raw.is_file():
            raw_missing.append(str(raw))
        sports_dir = root / 'sportsresult'
        if layout == 'nested':
            sports_dir /= '04_sports_annotation/samples'
        sports_dir /= sample
        if args.stage != 'upstream':
            required_files.extend(sports_dir / (sample + suffix) for suffix in ('_summary.txt', '_length_distribution.txt', '_output.txt'))
        if endo and args.stage != 'upstream':
            for spec in os.environ['ENDOSIRNA_REFS'].split(';'):
                tag, _ = spec.split('=', 1)
                prefix = root / 'sportsresult/05_endosirna/samples' / sample / f'{sample}.{tag}'
                required_files.extend(Path(str(prefix) + suffix) for suffix in ('.stats.tsv', '.sequence_counts.tsv'))
    refs = []
    if endo:
        tags = []
        for spec in os.environ['ENDOSIRNA_REFS'].split(';'):
            tag, prefix = spec.split('=', 1)
            tags.append(tag)
            candidates = [Path(prefix + '.fasta'), Path(prefix + '.fa')]
            found = next((p for p in candidates if p.is_file()), None)
            if found is None:
                raise ValueError('缺少 endo FASTA: ' + prefix)
            refs.append(found)
        if len(tags) != len(set(tags)):
            raise ValueError('ENDOSIRNA_REFS 标签重复')
    modules = set(filter(None, args.modules.replace(' ', ',').split(',')))
    common = {'qc', 'pca', 'sncrna', 'tsrna', 'rsrna', 'mirna', 'pirna'}
    applicable = common.copy()
    if species == 'Homo':
        applicable.add('snrna')
        if os.environ['PANDORA_DESIGN'] == 'method_capture':
            applicable.update({'similarity', 'spikein'})
    if endo:
        applicable.add('endo')
    if args.stage == 'differential':
        applicable = common - {'qc'} if species == 'Mus' else ({'capture', 'capture_bias'} if os.environ['PANDORA_DESIGN'] == 'method_capture' else set())
        if species == 'Mus' and endo:
            applicable.add('endo')
    elif args.stage == 'enrichment':
        applicable = {'mirna'} if species == 'Mus' else set()
    if len(rows) < 3:
        applicable.discard('pca')
    if args.stage != 'upstream' and not applicable:
        raise ValueError('本物种/设计没有可执行的所选阶段')
    if modules - applicable:
        raise ValueError('模块不适用于本阶段/设计: ' + ','.join(sorted(modules - applicable)))
    missing = [str(p) for p in required_files if not p.is_file()]
    if missing:
        raise ValueError('缺少输入文件: ' + '\n'.join(missing))
    if args.stage == 'upstream':
        if raw_missing:
            raise ValueError('原始 FASTQ 缺失: ' + '\n'.join(raw_missing))
        if not Path(os.environ['SPORTS_PL']).is_file():
            raise ValueError('缺少 SPORTS_PL: ' + os.environ['SPORTS_PL'])
        for resource in ('DB_ROOT', 'REFERENCE_ROOT'):
            if not Path(os.environ[resource]).is_dir():
                raise ValueError('缺少资源目录 ' + resource + ': ' + os.environ[resource])
        executables = ['perl', 'bowtie', 'bowtie-build', 'gzip', 'xargs', 'python3', 'sha256sum']
        if endo or os.environ.get('PRETRIM_FASTQ') == '1' or os.environ.get('TRIM_ADAPTER') == '1':
            executables.append('cutadapt')
    else:
        executables = ['Rscript']
    absent = [name for name in executables if shutil.which(name) is None]
    if absent:
        raise ValueError('缺少工具: ' + ','.join(absent))
    print(f'[CHECK] root={root} species={species} RUN_ENDOSIRNA={int(endo)} stage={args.stage} samples={len(rows)}')
    for key in ('DB_ROOT', 'REFERENCE_ROOT', 'SPORTS_PL'):
        print(f'[CHECK] {key}={os.environ[key]}')
    if raw_missing and args.stage != 'upstream':
        print(f'[CHECK] {len(raw_missing)} 个原始 FASTQ 缺失；已有 SPORTS 描述入口可用，原始 upstream 未验证。')
    print('[CHECK] 只读输入检查通过；不代表统计或结果已完成。')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        print('[ERROR] ' + str(error), file=sys.stderr)
        raise SystemExit(2)

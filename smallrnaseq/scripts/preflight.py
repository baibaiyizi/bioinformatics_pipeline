#!/usr/bin/env python3
"""Read-only TSV validation; no pandas or result writes."""
import argparse
import hashlib
import csv
import os
from pathlib import Path
import re
import sys

def samples(mode):
    path=Path(os.environ['SAMPLES_TSV'])
    with path.open() as f:
        reader=csv.DictReader(f,delimiter='\t'); rows=list(reader)
    required={'sample','group','biological_sample_id','unit_type','batch','raw_fq','condition','age_group','species','stage','replicate'}
    if not rows or not required.issubset(reader.fieldnames or []):
        raise ValueError('samples.tsv 为空或缺字段: '+str(required-set(reader.fieldnames or [])))
    names=[row['sample'] for row in rows]
    if len(set(names))!=len(names) or not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',name) for name in names):
        raise ValueError('sample 必须唯一且为安全文件名')
    if any(row.get(k, '') not in {'', 'unknown'} for row in rows for k in ['pair_id','subject_id','timepoint']):
        raise ValueError('paired/repeated designs are not supported by the independent core')
    groups={os.environ['PROJECT_CONTRAST_BASELINE'],os.environ['PROJECT_CONTRAST_CASE']}
    if len(groups)!=2 or {row['group'] for row in rows}!=groups:
        raise ValueError('显式 case/baseline 必须不同且与两组样本完全对应')
    if any(row['condition']!=row['group'] for row in rows):
        raise ValueError('condition 必须与显式 group 一致')
    root=Path(os.environ['SMALLRNA_ROOT_DIR'])
    for row in rows:
        if any(value is None or '\x1f' in value for value in row.values()):raise ValueError('TSV 列数错误或非法分隔符')
        if row['unit_type'] not in {'individual','pool','unknown'}:raise ValueError('unit_type 必须为 individual/pool/unknown')
        if row['species']!=os.environ['PROJECT_SPECIES']:raise ValueError('样本与配置物种不一致')
        if mode=='raw':
            raw=Path(row['raw_fq']); raw=raw if raw.is_absolute() else root/raw
            if not raw.is_file():raise FileNotFoundError('MISSING FASTQ: '+str(raw))
        elif mode=='downstream':
            for suffix in ('summary.txt','length_distribution.txt','output.txt'):
                p=Path(os.environ['SPORTS_RESULT_DIR'])/'04_sports_annotation/samples'/row['sample']/(row['sample']+'_'+suffix)
                if not p.is_file() or p.stat().st_size==0:raise FileNotFoundError('MISSING SPORTS: '+str(p))
    if os.environ['RUN_ENDOSIRNA']=='1':
        reference=os.environ.get('ENDOSIRNA_REFERENCE_TSV','')
        if not reference or not Path(reference).is_file():raise ValueError('RUN_ENDOSIRNA=1 需要显式 ENDOSIRNA_REFERENCE_TSV')
        with open(reference) as f:refs=list(csv.DictReader(f,delimiter='\t'))
        if not refs:raise ValueError('endo-siRNA 参考清单为空')
        for row in refs:
            p=Path(row['ReferenceFile']);p=p if p.is_absolute() else root/p
            if not p.is_file():raise FileNotFoundError('MISSING endo reference: '+str(p))
        if mode=='raw' and not os.environ.get('ENDOSIRNA_REFS',''):raise ValueError('上游 endo-siRNA 需要显式 ENDOSIRNA_REFS')
    return rows

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['raw','downstream','config'],required=True);parser.add_argument('--emit',action='store_true');parser.add_argument('--sha256',action='store_true');args=parser.parse_args()
    try:
        rows=samples(args.mode)
        if args.sha256:
            root=Path(os.environ['SMALLRNA_ROOT_DIR'])
            paths=[Path(os.environ['SMALLRNA_CONFIG']),Path(os.environ['SAMPLES_TSV'])]
            if args.mode=='downstream':
                paths += [Path(os.environ['SPORTS_RESULT_DIR'])/'04_sports_annotation/samples'/row['sample']/(row['sample']+'_'+suffix) for row in rows for suffix in ('summary.txt','length_distribution.txt','output.txt')]
            elif args.mode=='raw':
                paths += [Path(row['raw_fq']) if Path(row['raw_fq']).is_absolute() else root/row['raw_fq'] for row in rows]
            print('path\tsha256')
            for path in paths:
                digest=hashlib.sha256()
                with path.open('rb') as handle:
                    for block in iter(lambda:handle.read(16*1024*1024),b''):digest.update(block)
                print(str(path)+'\t'+digest.hexdigest())
        elif args.emit:
            for row in rows:print('\x1f'.join(row[k] for k in ['group','sample','replicate','raw_fq','condition','age_group','species','stage']))
        else:print(f'PASS {args.mode}: {len(rows)} libraries; design identity reviewed separately')
    except Exception as exc:
        print(f'FAIL: {exc}',file=sys.stderr);sys.exit(1)

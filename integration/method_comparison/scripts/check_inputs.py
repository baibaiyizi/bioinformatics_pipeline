#!/usr/bin/env python3
"""核验单组学当前生产记录和血清实验单位；只读。"""
import csv
import hashlib
import os
from pathlib import Path


def rows(path):
    with path.open(newline='',encoding='utf-8') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def main():
    root=Path(os.environ['PANDORA_PROJECT_ROOT'])
    source=root/'result/describe'
    receipt=rows(source/'run_status.tsv')
    if len(receipt)!=1 or receipt[0]['State']!='completed':raise ValueError('not_run: 单组学描述尚未完成')
    if not any(x['Module']=='qc' and x['Status']=='completed' for x in rows(source/'module_status.tsv')):raise ValueError('not_run: 单组学QC未完成')
    for row in rows(source/'input_signature.tsv'):
        path=Path(row['Path'])
        if not path.is_file() or path.stat().st_size!=int(row['Size']):
            raise ValueError('STALE_OUTPUTS: 单组学输入不存在或大小变化 '+str(path))
        md5=hashlib.md5()
        with path.open('rb') as handle:
            for chunk in iter(lambda:handle.read(8*1024*1024),b''):md5.update(chunk)
        if md5.hexdigest()!=row['MD5']:raise ValueError('STALE_OUTPUTS: 单组学输入或代码改变 '+str(path))
    metadata=rows(root/'config/samples.tsv')
    if len(metadata)!=12 or len({x['sample'] for x in metadata})!=12:raise ValueError('UNVERIFIED: 血清设计必须为12个不同文库')
    comparison=[x for x in metadata if x['analysis_role']=='method_comparison']
    paired=[x for x in comparison if x['method'] in {'1','2'}]
    if len(comparison)!=10 or {x['biological_sample'] for x in comparison}!={'2'}:raise ValueError('UNVERIFIED: 方法比较必须来自供体2的10个文库')
    if len(paired)!=6 or len({x['pair_id'] for x in paired})!=3:raise ValueError('UNVERIFIED: 方法1/2必须为三个匹配PairID和六库')
    for pair in {x['pair_id'] for x in paired}:
        members=[x for x in paired if x['pair_id']==pair]
        if len(members)!=2 or {x['method'] for x in members}!={'1','2'}:raise ValueError('UNVERIFIED: PairID必须各包含一个方法1和方法2文库')
    if not (source/'01_setup/1_expression_matrix/class_reads.tsv').is_file():raise ValueError('MISSING: 标准类别计数矩阵')
    if not (source/'02_single/tables/single_sample_qc_table.csv').is_file():raise ValueError('MISSING: 标准QC表')
    for row in metadata:
        if row['species'] not in {'Homo_sapiens','Homo sapiens'}:raise ValueError('UNVERIFIED: 方法比较只接受人类血清物种配置')
        for suffix in ('_summary.txt','_output.txt','_length_distribution.txt'):
            if not (root/'sportsresult'/row['sample']/(row['sample']+suffix)).is_file():raise ValueError('MISSING: SPORTS标准注释 '+row['sample']+suffix)
    print('[CHECK] 当前单组学签名、12库身份与3个方法1/2技术配对通过；DonorN=1，方法3–6仅描述，donor3外参独立QC。')


if __name__=='__main__':
    try:main()
    except (OSError, ValueError, KeyError) as error:
        import sys
        print('[ERROR] '+str(error),file=sys.stderr)
        raise SystemExit(2)

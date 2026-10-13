#!/usr/bin/env python3
"""显式标准产物驱动的课题联合分析；不猜测配对，不复用旧联合候选表。"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import yaml


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def load_table(path):
    return pd.read_csv(path, sep='\t' if Path(path).suffix in {'.tsv', '.matrix'} else ',', low_memory=False)


def gene_id(values):
    return values.astype('string').str.replace(r'\.[0-9]+$', '', regex=True)


def validate_direction(df, spec, cfg):
    reference, comparison = cfg['reference'], cfg['comparison']
    if spec['kind'] == 'rna':
        for col, expected in [('sampleA', comparison), ('sampleB', reference)]:
            if col not in df or df[col].isna().any() or set(df[col].unique()) != {expected}:
                raise ValueError(f"{spec['id']}: {col} 必须显式为 {expected}")
    elif spec['kind'] == 'wgbs':
        for col, expected in [('group1_label', reference), ('group2_label', comparison)]:
            if col not in df or df[col].isna().any() or set(df[col].unique()) != {expected}:
                raise ValueError(f"{spec['id']}: WGBS方向未知或反向")
        effect = pd.to_numeric(df['diffMethy'], errors='coerce')
        expected = pd.to_numeric(df['meanMethy_group2'], errors='coerce') - pd.to_numeric(df['meanMethy_group1'], errors='coerce')
        if not np.isfinite(effect).all() or not np.isfinite(expected).all() or not np.allclose(effect, expected, rtol=0, atol=1e-7):
            raise ValueError(f"{spec['id']}: diffMethy必须逐行等于PFOS均值减WT均值")
    elif spec['kind'] == 'atac':
        if 'Direction' not in df:
            raise ValueError(f"{spec['id']}: ATAC缺少显式Direction")
        effect = pd.to_numeric(df[spec['effect_column']], errors='coerce')
        if df['Direction'].isna().any() or not np.isfinite(effect).all():
            raise ValueError(f"{spec['id']}: ATAC方向或效应缺失/非有限")
        for label, positive in [(comparison + '-open', True), (reference + '-open', False)]:
            selected = df['Direction'] == label
            if (effect[selected] <= 0 if positive else effect[selected] >= 0).any():
                raise ValueError(f"{spec['id']}: Direction与效应符号冲突")
        if not set(df['Direction'].dropna()).issubset({comparison + '-open', reference + '-open'}):
            raise ValueError(f"{spec['id']}: ATAC含未知方向")


class Integration:
    def __init__(self, root):
        self.root = Path(root).resolve()
        with (self.root / 'config/project.yml').open() as f:
            self.cfg = yaml.safe_load(f)
        if self.cfg['reference'] == self.cfg['comparison']:
            raise ValueError('reference与comparison必须不同')
        self.result = self.root / 'result'
        self.report = self.root / 'README.md'
        self.html_report = self.result / 'index.html'
        self.tables = {}
        self.specs = {s['id']: s for s in self.cfg['inputs']}
        if len(self.specs) != len(self.cfg['inputs']):
            raise ValueError('输入id重复')
        self.audit = []

    def input_path(self, spec):
        return (self.root / spec['path']).resolve()

    def check(self, include_frozen_contract=False):
        self.audit = []
        for spec in self.specs.values():
            path = self.input_path(spec)
            row = dict(id=spec['id'], assay=spec['kind'], context=spec['context'],
                       path=spec['path'], species=spec['species'], genome=spec['genome'],
                       source_status=spec['source_status'], biological_unit=spec.get('biological_unit', 'unknown'),
                       pairing=spec.get('pairing', 'unknown'), quality=spec.get('quality', 'unknown'))
            try:
                if spec['species'] != self.cfg['species'] or spec['genome'] != self.cfg['genome']:
                    raise ValueError('物种或基因组不一致')
                if not path.is_file() or not path.stat().st_size:
                    raise FileNotFoundError('MISSING')
                df = load_table(path)
                required = spec.get('required_columns', [])
                if not set(required).issubset(df):
                    raise ValueError(f"缺列: {set(required) - set(df)}")
                validate_direction(df, spec, self.cfg)
                if spec['kind'] == 'rna':
                    self.check_rna_producer(path, df)
                row.update(status='PASS', rows=len(df), sha256=sha256(path), reason='')
                self.tables[spec['id']] = df
            except (OSError, ValueError, KeyError) as e:
                row.update(status='MISSING' if isinstance(e, FileNotFoundError) else 'INVALID', rows=0, sha256='', reason=str(e))
            self.audit.append(row)
        resources = {'gene_sets': self.cfg['gene_sets'], 'figure_profile': self.cfg['figure_profile']}
        resources.update(self.cfg.get('joint_inputs', {}))
        for context in self.cfg.get('rna_context_order', []):
            for key, template in self.cfg.get('deepening', {}).get('posttranscriptional', {}).items():
                resources[f'deepening:{context}:{key}'] = template.format(context=context)
        for spec in self.specs.values():
            for key in ('expression', 'universe', 'cpg_rdata', 'samples', 'sample_methylation', 'all_sites'):
                if spec.get(key): resources[spec['id'] + ':' + key] = spec[key]
            for group, paths in spec.get('tracks', {}).items():
                for index, relative in enumerate(paths):
                    resources[f"{spec['id']}:track:{group}:{index+1}"] = relative
        for sid, relative in resources.items():
            path = self.root / relative
            exists = path.is_file() and path.stat().st_size > 0
            self.audit.append(dict(id=sid, assay='resource', context='', path=relative,
                                   species=self.cfg['species'], genome=self.cfg['genome'],
                                   source_status='explicit_frozen_resource', biological_unit='not_applicable',
                                   pairing='not_applicable', quality='source_content_only',
                                   status='PASS' if exists else 'MISSING', rows=0,
                                   sha256=sha256(path) if exists else '', reason='' if exists else 'MISSING'))
        if include_frozen_contract:
            required=['00_audit/run.json']+['00_audit/tables/'+n for n in ['source_catalog.tsv','figure_manifest.tsv','membership.tsv','term_id_map.tsv']]
            missing=[str(Path('result')/p) for p in required if not (self.result/p).is_file()]
            status='MISSING' if missing else 'PASS';reason='MISSING: '+', '.join(missing) if missing else ''
            if not missing:
                try:
                    if json.loads((self.result/'00_audit/run.json').read_text())['status']!='COMPLETED':raise ValueError('现役冻结生产不是COMPLETED')
                except (ValueError,KeyError) as e:status='INVALID';reason=str(e)
            self.audit.append(dict(id='frozen_result_contract',assay='contract',context='project',path='result/00_audit',source_status='existing_completed_numeric_sources',status=status,reason=reason,sha256='',rows=0))
        return pd.DataFrame(self.audit)

    def check_rna_producer(self, path, df):
        provenance = path.parent / 'provenance.json'
        p = json.loads(provenance.read_text())
        if not isinstance(p.get('outputs'), dict) or not isinstance(p.get('source'), list):
            raise ValueError('RNA生产凭据格式不完整')
        if {s['role'] for s in p['source']} != {'counts', 'samples', 'annotation', 'producer'}:
            raise ValueError('RNA生产凭据缺少必要来源')
        execution = load_table(path.parent / 'execution.tsv')
        if execution.empty or execution.iloc[-1]['status'] != 'COMPLETED':
            raise ValueError('RNA当前生产状态不是COMPLETED')
        if p.get('status') != 'completed' or p['outputs'].get(path.name) != sha256(path):
            raise ValueError('RNA生产状态或DE签名不匹配')
        samples = load_table(path.parent / 'samples.tsv')
        expression_path = path.parent / 'normalized_expression.tsv'
        if p['outputs'].get('samples.tsv') != sha256(path.parent / 'samples.tsv') or p['outputs'].get('normalized_expression.tsv') != sha256(expression_path):
            raise ValueError('RNA样本表或表达矩阵签名失配')
        expression = pd.read_csv(expression_path, sep='\t', nrows=0)
        if set(expression.columns[1:]) != set(samples['sample']) or set(p['samples']) != set(samples['sample']):
            raise ValueError('RNA DE/表达矩阵/样本范围不一致')
        if samples['sample'].duplicated().any() or df['gene_id'].duplicated().any():
            raise ValueError('RNA样本或基因身份重复')
        module = path.parents[2]
        for source in p['source']:
            # 实际执行源码被冻结；后来入口的修订不能改写已完成模型的来源。
            original = path.parent / 'producer_source.R' if source['role'] == 'producer' else Path(source['path'])
            if not original.is_absolute(): original = module / original
            if not original.is_file() or sha256(original) != source['sha256']:
                raise ValueError('RNA上游输入或producer签名失配: ' + source['role'])

    def require(self, ids):
        failures = [s for s in self.audit if s['id'] in ids and s['status'] != 'PASS']
        if failures:
            raise ValueError('输入资格未通过: ' + '; '.join(s['id'] + ':' + s['reason'] for s in failures))

    def write(self, frame, name):
        path = self.result / name
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False, sep='\t' if path.suffix == '.tsv' else ',')
        return path

    def code_signatures(self):
        files = [self.root / 'run.sh', self.root / 'config/project.yml']
        files.extend(p for p in (self.root / 'scripts').rglob('*')
                     if p.suffix in {'.py', '.R', '.Rmd', '.sh'} and '__pycache__' not in p.parts)
        return {str(p.relative_to(self.root)): sha256(p) for p in sorted(files) if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description='冻结结果驱动的项目级 PFOS 整合与纯绘图')
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--check', action='store_true')
    action.add_argument('--integrate', action='store_true')
    action.add_argument('--plots-only', action='store_true')
    args = parser.parse_args()
    workflow = Integration(args.root)
    if args.check:
        audit = workflow.check(include_frozen_contract=True)
        print(audit[['id', 'status', 'source_status', 'reason']].to_string(index=False))
        return 0 if (audit.status == 'PASS').all() else 1
    from project_integration import run_project
    return run_project(workflow, integrate=args.integrate)


if __name__ == '__main__':
    raise SystemExit(main())

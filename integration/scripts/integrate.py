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

    def check(self):
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

    def rna(self):
        ids = [s['id'] for s in self.specs.values() if s['kind'] == 'rna']
        self.require(ids)
        frames = []
        for sid in ids:
            spec, d = self.specs[sid], self.tables[sid].copy()
            d['gene_id'] = gene_id(d['ensembl_id'])
            d['context'] = spec['context']
            d['nominal_candidate'] = (d.pvalue < self.cfg['p_cutoff']) & (d.log2FoldChange.abs() >= self.cfg['lfc_cutoff'])
            d['fdr_candidate'] = (d.padj < self.cfg['fdr_cutoff']) & (d.log2FoldChange.abs() >= self.cfg['lfc_cutoff'])
            frames.append(d[['context', 'gene_id', 'SYMBOL', 'log2FoldChange', 'pvalue', 'padj', 'nominal_candidate', 'fdr_candidate']])
        long = pd.concat(frames, ignore_index=True)
        self.write(long, 'cross_context/rnaseq/effects.csv')
        effects = long.pivot(index='gene_id', columns='context', values='log2FoldChange')
        self.write(effects.reset_index(), 'cross_context/rnaseq/effect_matrix.csv')
        summary = long.groupby('context').agg(tested=('pvalue', 'count'), nominal=('nominal_candidate', 'sum'), fdr=('fdr_candidate', 'sum')).reset_index()
        self.write(summary, 'cross_context/rnaseq/summary.csv')
        candidates = long[long.nominal_candidate]
        agreement = candidates.groupby('gene_id').agg(contexts=('context', lambda x: ';'.join(x)), n_contexts=('context', 'nunique'), min_effect=('log2FoldChange', 'min'), max_effect=('log2FoldChange', 'max')).reset_index()
        agreement['opposite_directions'] = (agreement.min_effect < 0) & (agreement.max_effect > 0)
        self.write(agreement, 'cross_context/rnaseq/candidate_contexts.csv')
        self.plot_rna(summary)

    def regions(self, sid):
        s, d = self.specs[sid], self.tables[sid].copy()
        chrom = s.get('chrom_column', 'seqnames')
        region = pd.DataFrame({'region_id': d[s['region_column']].astype(str), 'chrom': d[chrom].astype(str), 'start': d['start'], 'end': d['end']})
        if region.region_id.duplicated().any():
            region = region.drop_duplicates()
            if region.region_id.duplicated().any():
                raise ValueError(sid + ': 同一区域id对应冲突坐标')
        return region

    def overlap(self, left, right, name):
        a = self.write(self.regions(left), f'inputs/{left}_regions.csv')
        b = self.write(self.regions(right), f'inputs/{right}_regions.csv')
        out = self.result / name
        out.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(['Rscript', '--vanilla', str(self.root / 'scripts/overlap_regions.R'), str(a), str(b), str(out)], check=True)
        return load_table(out)

    def wgbs(self):
        ids = [s['id'] for s in self.specs.values() if s['kind'] == 'wgbs']
        self.require(ids)
        summary = []
        for i, left in enumerate(ids):
            for right in ids[i + 1:]:
                name = left + '__' + right
                overlaps = self.overlap(left, right, f'cross_context/wgbs/{name}_overlap.csv')
                a, b = self.tables[left], self.tables[right]
                ga = set(gene_id(a[self.specs[left]['gene_column']]).dropna())
                gb = set(gene_id(b[self.specs[right]['gene_column']]).dropna())
                shared = sorted(ga & gb)
                self.write(pd.DataFrame({'gene_id': shared}), f'cross_context/wgbs/{name}_shared_genes.csv')
                quality = '; '.join(sorted({self.specs[left].get('quality', 'unknown'), self.specs[right].get('quality', 'unknown')}))
                summary.append(dict(left=left, right=right, overlap_pairs=len(overlaps), left_overlapping_regions=overlaps.left_region.nunique(), right_overlapping_regions=overlaps.right_region.nunique(), shared_annotated_genes=len(shared), interpretation='跨组织位置/最近基因注释交汇；不等于甲基化遗传', quality=quality))
        self.write(pd.DataFrame(summary), 'cross_context/wgbs/summary.csv')

    def gene_evidence(self, sid):
        s, d = self.specs[sid], self.tables[sid].copy()
        d['gene_id'] = gene_id(d[s['gene_column']])
        d['effect'] = pd.to_numeric(d[s['effect_column']], errors='coerce')
        d['region_feature_annotation'] = d['annotation'].fillna('unknown')
        d['source_region'] = d[s['region_column']].astype(str)
        valid = d.gene_id.notna() & d.gene_id.str.match(r'^ENS(?:MUS)?G[0-9]+$', na=False)
        self.write(d[~valid], f'multiomics/liver/{s["kind"]}_unmapped_rows.csv')
        long = d.loc[valid, ['gene_id', 'source_region', 'region_feature_annotation', 'effect']].copy()
        long['positive'] = long.effect > 0
        long['negative'] = long.effect < 0
        self.write(long, f'multiomics/liver/{s["kind"]}_region_gene_evidence.csv')
        summary = long.groupby('gene_id').agg(regions=('source_region', 'nunique'), positive=('positive', 'sum'), negative=('negative', 'sum'), median_effect=('effect', 'median'), min_effect=('effect', 'min'), max_effect=('effect', 'max')).reset_index()
        summary['mixed_directions'] = (summary.positive > 0) & (summary.negative > 0)
        summary = summary.rename(columns={c: s['kind'] + '_' + c for c in summary.columns if c != 'gene_id'})
        return summary

    def liver(self):
        ids = self.cfg['liver_inputs']
        self.require(ids)
        by_kind = {self.specs[i]['kind']: i for i in ids}
        rna = self.tables[by_kind['rna']].copy()
        rna['gene_id'] = gene_id(rna.ensembl_id)
        evidence = rna[['gene_id', 'SYMBOL', 'log2FoldChange', 'pvalue', 'padj']].rename(columns={'log2FoldChange': 'rna_log2FC', 'pvalue': 'rna_pvalue', 'padj': 'rna_padj'})
        for kind in ('atac', 'wgbs'):
            evidence = evidence.merge(self.gene_evidence(by_kind[kind]), on='gene_id', how='outer', validate='one_to_one')
        evidence['rna_nominal_candidate'] = (evidence.rna_pvalue < self.cfg['p_cutoff']) & (evidence.rna_log2FC.abs() >= self.cfg['lfc_cutoff'])
        evidence['rna_fdr_candidate'] = (evidence.rna_padj < self.cfg['fdr_cutoff']) & (evidence.rna_log2FC.abs() >= self.cfg['lfc_cutoff'])
        evidence['has_atac_region'] = evidence.atac_regions.fillna(0) > 0
        evidence['has_wgbs_region'] = evidence.wgbs_regions.fillna(0) > 0
        evidence['evidence_scope'] = '最近基因链接的区域层探索证据；样本关系及上游来源见输入审计；区域feature注释未必属于最近基因'
        self.write(evidence, 'multiomics/liver/gene_evidence.csv')
        candidates = evidence[evidence.rna_nominal_candidate & (evidence.has_atac_region | evidence.has_wgbs_region)]
        self.write(candidates, 'multiomics/liver/candidates.csv')
        self.overlap(by_kind['atac'], by_kind['wgbs'], 'multiomics/liver/atac_wgbs_overlap.csv')
        metab = self.tables[by_kind['metabolomics']].copy()
        keep = [c for c in ['feature_id', 'name', 'id_kegg', 'confidence_level', 'annotation_ambiguous', 'annotation_match_type', 'log2FC', 'limma_p', 'fdr', 'classification', 'direction'] if c in metab]
        self.write(metab[keep], 'multiomics/liver/metabolic_features.csv')
        mpath = self.cfg.get('metabolic_pathways')
        if mpath:
            path = (self.root / mpath).resolve()
            if path.is_file(): self.write(load_table(path), 'multiomics/liver/metabolic_pathway_effects.csv')
        self.write(pd.DataFrame([dict(connection='metabolite_to_gene', status='MISSING', reason='缺少核实的代谢特征-化合物-反应-基因成员映射；代谢通路效应单独并列，不虚构RNA耦合')]), 'multiomics/liver/connection_eligibility.csv')
        self.plot_liver(evidence)

    def programs(self):
        evidence = load_table(self.result / 'multiomics/liver/gene_evidence.csv')
        selected = evidence[evidence.rna_nominal_candidate | evidence.has_atac_region | evidence.has_wgbs_region].copy()
        genes = self.write(selected[['gene_id']], 'programs/genes.tsv')
        members = self.result / 'programs/go_members.tsv'
        subprocess.run(['Rscript', '--vanilla', str(self.root / 'scripts/annotate_programs.R'), str(genes), self.cfg['orgdb'], str(members)], check=True)
        m = load_table(members).merge(selected, on='gene_id', how='inner', validate='many_to_one')
        counts = m.groupby(['GO', 'term']).agg(rna_candidate_genes=('rna_nominal_candidate', 'sum'), atac_annotated_genes=('has_atac_region', 'sum'), wgbs_annotated_genes=('has_wgbs_region', 'sum'), genes=('gene_id', lambda x: ';'.join(sorted(set(x))))).reset_index()
        counts['interpretation'] = 'GO注释交汇；不是富集检验或通路活性'
        self.write(counts, 'programs/go_evidence.csv')

    def plot_rna(self, summary):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(183 / 25.4, 100 / 25.4))
        x = np.arange(len(summary))
        ax.bar(x - .18, summary.nominal, .36, label=f"Nominal P < {self.cfg['p_cutoff']:g} + effect", color='#e59a9e')
        ax.bar(x + .18, summary.fdr, .36, label=f"BH < {self.cfg['fdr_cutoff']:g} + effect", color='#88abc9')
        ax.set_xticks(x, summary.context, rotation=35, ha='right')
        ax.set_ylabel('RNA candidate genes'); ax.legend(frameon=False, fontsize=8)
        ax.spines[['top', 'right']].set_visible(False)
        self.save_figure(fig, 'cross_context/rnaseq', 'rna_context_summary')

    def plot_liver(self, evidence):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        d = evidence.dropna(subset=['rna_log2FC', 'atac_median_effect'])
        fig, ax = plt.subplots(figsize=(183 / 25.4, 110 / 25.4))
        for present, label, color in [(False, 'ATAC + RNA', '#88abc9'), (True, '+ WGBS region', '#e59a9e')]:
            s = d[d.has_wgbs_region == present]
            ax.scatter(s.rna_log2FC, s.atac_median_effect, s=8, alpha=.6, c=color, label=label, rasterized=True)
        ax.axhline(0, color='#bbbbbb', lw=.6); ax.axvline(0, color='#bbbbbb', lw=.6)
        contrast = self.cfg['comparison'] + ' - ' + self.cfg['reference']
        ax.set(xlabel='RNA log2FC (' + contrast + ')', ylabel='Median linked ATAC effect (' + contrast + ')')
        ax.legend(frameon=False, fontsize=8); ax.spines[['top', 'right']].set_visible(False)
        self.save_figure(fig, 'multiomics/liver', 'liver_gene_effects')

    def save_figure(self, fig, stage_folder, name):
        folder = self.result / stage_folder / 'plots'
        folder.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(folder / (name + '.pdf'))
        fig.savefig(folder / (name + '.png'), dpi=300)
        import matplotlib.pyplot as plt
        plt.close(fig)

    def report_header(self):
        return ['# ' + self.cfg['title'], '', '## 运行与结果', '',
                '```bash', 'bash run.sh --check', 'bash run.sh --all', 'bash run.sh --stage liver', '```', '',
                '输入路径基于本分析目录；从其他目录启动时使用入口脚本的路径。',
                '研究设计见所属课题README；全部图表直接保存在本模块result，不另外复制精选版。',
                '[完整结果](result/) · [浏览报告](result/index.html) · [执行记录](result/run.json)', '']

    def synthesis(self, stages):
        self.result.mkdir(parents=True, exist_ok=True)
        completed = {s['stage'] for s in stages if s['status'] == 'COMPLETED'}
        lines = self.report_header() + [ '先看研究问题与证据，再查各组学的独立报告。比较方向为 **' + self.cfg['comparison'] + ' − ' + self.cfg['reference'] + '**。', '',
                 '本报告的计数与连接来自当前签名通过的联合运行。RNA使用具有实际完成凭据的完整cohort下游产物；各层样本单位、配对、质量和上游来源身份逐项见输入审计。联合计算的结论边界由原始测量层和实验单位决定；来源核验通过不会自动提升上游复现状态。', '', '## 1. 跨组织 RNA：哪些变化共享，哪些依赖组织？', '']
        rna = self.result / 'cross_context/rnaseq/summary.csv'
        if 'rna' in completed and rna.exists():
            d = load_table(rna)
            lines += [f"探索候选使用原始P<{self.cfg['p_cutoff']:g}且|log2FC|≥{self.cfg['lfc_cutoff']:g}；BH支持候选使用FDR<{self.cfg['fdr_cutoff']:g}及相同效应阈值。[关键汇总表](result/cross_context/rnaseq/summary.csv)。", '']
            lines += ['| 组织／阶段 | 有效检验基因 | 探索候选 | BH支持候选 |', '|---|---:|---:|---:|']
            lines += [f'| {r.context} | {r.tested} | {r.nominal} | {r.fdr} |' for r in d.itertuples()]
            lines += ['', '![跨组织RNA候选概览](result/cross_context/rnaseq/plots/rna_context_summary.png)', '', '生物学意义：比较组织特异性与共享响应。数据分布特征：图中每个柱对应上表的实际计数，使用相同效应阈值。最终结论：共享候选不等于跨组织传递，阶段也不是同一个体纵向追踪。', '']
        else: lines += ['本模块未完成。', '']
        lines += ['## 2. 肝脏：开放、甲基化和表达是否落在共同基因／区域？', '']
        path = self.result / 'multiomics/liver/candidates.csv'
        if 'liver' in completed and path.exists():
            d = load_table(path)
            lines += [f'当前共有 **{len(d)}** 个 RNA 探索候选具有通过最近基因注释链接的 ATAC 或 WGBS 区域；其中 **{int((d.has_atac_region & d.has_wgbs_region).sum())}** 个同时具有两类区域证据。完整效应、混合方向及源区域见[候选证据表](result/multiomics/liver/candidates.csv)。', '', '![肝脏基因层效应](result/multiomics/liver/plots/liver_gene_effects.png)', '', '生物学意义：查找不同测量层的共同落点。数据分布特征：每点是基因，ATAC多区域以中位效应展示，完整正负区域及混合方向保留在源表。区域feature注释保存在region_feature_annotation，可能属于另一个基因，不能当作最近基因自身的promoter/intron效应。最终结论：当前是注释连接，配对相关、远端靶向和调控因果均需另外验证。', '']
        else: lines += ['本模块未完成。', '']
        if 'rna' in completed:
            rna_effects = load_table(self.result / 'cross_context/rnaseq/effects.csv')
            symbols = set(rna_effects.loc[rna_effects.nominal_candidate, 'SYMBOL'].dropna())
            for spec in self.specs.values():
                if spec['kind'] == 'celltypes' and spec['id'] in self.tables:
                    cells = self.tables[spec['id']]
                    candidates = cells[cells.gene.isin(symbols)].copy()
                    candidates['count_per_cell'] = candidates['count'] / candidates['cells'].replace(0, np.nan)
                    candidates['interpretation'] = '细胞类型描述背景；单位=' + spec.get('biological_unit', 'unknown') + '；配对=' + spec.get('pairing', 'unknown') + '；来源=' + spec['source_status'] + '；未执行条件效应检验'
                    self.write(candidates, 'synthesis/candidate_celltype_context.csv')
        lines += ['## 3. 代谢与程序证据', '', '代谢特征保留鉴定等级、歧义、原始P和FDR；代谢通路效应与染色质／RNA证据并列阅读。当前缺少经核实的化合物—反应—基因映射，因此相关连接明确为MISSING。GO表描述基因注释交汇，不作为富集显著性或通路活性。', '',
                  ]
        if 'liver' in completed:
            metabolic = load_table(self.result / 'multiomics/liver/metabolic_features.csv')
            nominal = pd.to_numeric(metabolic.limma_p, errors='coerce') < self.cfg['p_cutoff']
            adjusted = pd.to_numeric(metabolic.fdr, errors='coerce') < self.cfg['fdr_cutoff']
            lines += [f"当前纳入 **{len(metabolic)}** 个代谢feature；原始P<{self.cfg['p_cutoff']:g}有 **{int(nominal.sum())}** 个，FDR<{self.cfg['fdr_cutoff']:g}有 **{int(adjusted.sum())}** 个（本段只按P/FDR计数）。feature数量与已确证化合物数量不同，鉴定等级和歧义见[代谢关键表](result/multiomics/liver/metabolic_features.csv)。", '']
        if 'programs' in completed:
            programs = load_table(self.result / 'programs/go_evidence.csv')
            triple = (programs.rna_candidate_genes > 0) & (programs.atac_annotated_genes > 0) & (programs.wgbs_annotated_genes > 0)
            lines += [f'GO注释交汇表包含 **{len(programs)}** 个BP条目，其中 **{int(triple.sum())}** 个条目在三个测量层各有至少一个对应基因。GO条目具有层级依赖；这只是注释覆盖盘点，不能视为独立显著通路。[完整程序证据](result/programs/go_evidence.csv)。', '']
        lines += ['## 4. 跨情境与细胞类型背景', '', '跨组织WGBS计算坐标交集和最近基因注释交汇；不同阶段不视为同一个体纵向追踪。细胞类型背景按自身实验单位解释。', '']
        if 'wgbs' in completed:
            summary = load_table(self.result / 'cross_context/wgbs/summary.csv')
            lines += ['| WGBS情境比较 | 区域交叠对数 | 左侧命中区域 | 右侧命中区域 | 共同最近基因 |', '|---|---:|---:|---:|---:|']
            lines += [f'| {r.left}／{r.right} | {r.overlap_pairs} | {r.left_overlapping_regions} | {r.right_overlapping_regions} | {r.shared_annotated_genes} |' for r in summary.itertuples()]
            lines += ['', '同一区域可能命中多个区域，因此交叠对数不同于独立区域数。[跨情境关键表](result/cross_context/wgbs/summary.csv)。', '']
        lines += [note + '\n' for note in self.cfg.get('context_notes', [])]
        lines += ['## 5. 数据与执行状态', '', '| 输入 | 来源身份 | 核验 | 生物学单位 | 配对 | 质量说明 |', '|---|---|---|---|---|---|']
        lines += [f"| {a['id']} | {a['source_status']} | {a['status']} | {a['biological_unit']} | {a['pairing']} | {a['quality']} |" for a in self.audit]
        lines += ['', '| 联合阶段 | 本次状态 |', '|---|---|']
        lines += [f"| {s['stage']} | {s['status']} |" for s in stages]
        lines += ['', '[完整输入与SHA256](result/input_audit.tsv) · [联合源表](result/) · [分析入口](README.md)', '']
        markdown = '\n'.join(lines)
        self.report.write_text(markdown)
        # 使用已有markdown库渲染；缺依赖时明确失败，不把纯文本伪装成HTML报告。
        import markdown as md
        body = md.markdown(markdown, extensions=['tables'])
        self.html_report.write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><base href="../"><title>' + html.escape(self.cfg['title']) + '</title><style>body{max-width:1000px;margin:40px auto;font:17px/1.7 sans-serif;padding:0 20px;color:#263238}img{max-width:100%}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:6px 12px}a{color:#28658c}</style><body>' + body + '</body></html>')

    def run(self, stage):
        self.result.mkdir(parents=True, exist_ok=True)
        audit = self.check()
        self.write(audit, 'input_audit.tsv')
        stages = []
        old = self.result / 'run.json'
        if old.exists():
            previous = json.loads(old.read_text())
            if previous.get('input_signatures') == {r['id']: r['sha256'] for r in self.audit} and previous.get('code_signatures') == self.code_signatures():
                for s in previous['stages']:
                    outputs = s.get('outputs', {})
                    figure = {'rna': 'cross_context/rnaseq/plots/rna_context_summary', 'liver': 'multiomics/liver/plots/liver_gene_effects'}.get(s['stage'])
                    figure_keys = {figure + '.' + ext for ext in ['png', 'pdf']} if figure else set()
                    if s['stage'] != 'synthesis' and s['status'] == 'COMPLETED' and outputs and figure_keys.issubset(outputs) and all((self.result / p).is_file() and sha256(self.result / p) == sig for p, sig in outputs.items()):
                        stages.append(s)
            archive = self.root / 'history' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
            archive.mkdir(parents=True)
            shutil.copy2(old, archive / 'run.json')
            for path in [self.report, self.html_report]:
                if path.is_file():
                    target = archive / path.relative_to(self.root)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(path, target)
        synthesis_folder = self.result / 'synthesis'
        if synthesis_folder.exists():
            archived = self.root / 'history' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ') / 'synthesis'
            archived.parent.mkdir(parents=True)
            synthesis_folder.rename(archived)
        def persist(status):
            record = dict(time=datetime.now(timezone.utc).isoformat(), status=status, stages=stages,
                          input_signatures={r['id']: r['sha256'] for r in self.audit}, code_signatures=self.code_signatures(),
                          interpretation='本次联合产物；上游来源与生物学单位限制逐项保留')
            old.write_text(json.dumps(record, ensure_ascii=False, indent=2))

        def status_report(status, reason):
            message = '\n'.join(self.report_header()) + '\n当前执行状态：**' + status + '**。\n\n' + reason + '\n'
            self.report.write_text(message)
            self.html_report.write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><body><h1>' + html.escape(self.cfg['title']) + '</h1><p>' + html.escape(status + ': ' + reason) + '</p></body></html>')

        persist('STARTED')
        status_report('STARTED', '本次运行尚未完成。历史报告已保存到本模块内部归档。')
        actions = [('rna', self.rna), ('wgbs', self.wgbs), ('liver', self.liver), ('programs', self.programs)]
        folders = {'rna': 'cross_context/rnaseq', 'wgbs': 'cross_context/wgbs', 'liver': 'multiomics/liver', 'programs': 'programs'}
        for name, action in actions:
            if stage not in ('all', name): continue
            stages = [s for s in stages if s['stage'] != name]
            folder = self.result / folders[name]
            if folder.exists():
                archived = self.root / 'history' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ') / folders[name]
                archived.parent.mkdir(parents=True, exist_ok=True)
                folder.rename(archived)
            try:
                if name == 'programs' and not any(s['stage'] == 'liver' and s['status'] == 'COMPLETED' for s in stages):
                    raise ValueError('需要当前已完成且签名匹配的liver阶段')
                action(); status, reason = 'COMPLETED', ''
            except Exception as e:
                status, reason = 'FAILED', str(e)
            outputs = {str(p.relative_to(self.result)): sha256(p) for p in sorted(folder.rglob('*')) if p.is_file()} if folder.exists() and status == 'COMPLETED' else {}
            stages.append(dict(stage=name, status=status, reason=reason, outputs=outputs))
            persist('STARTED')
            print(name, status, reason, flush=True)
        try:
            self.synthesis(stages)
            paths = [self.report, self.html_report]
            celltypes = self.result / 'synthesis/candidate_celltype_context.csv'
            if celltypes.exists() and 'rna' in {s['stage'] for s in stages if s['status'] == 'COMPLETED'}:
                paths.append(celltypes)
            stages.append(dict(stage='synthesis', status='COMPLETED', reason='', outputs={os.path.relpath(p, self.result): sha256(p) for p in paths}))
        except Exception as e:
            stages.append(dict(stage='synthesis', status='FAILED', reason=str(e), outputs={}))
            status_report('FAILED', '报告生成失败：' + str(e) + '。请查看联合目录result/run.json。')
        finally:
            failed = any(s['status'] == 'FAILED' for s in stages)
            complete = {'rna', 'wgbs', 'liver', 'programs', 'synthesis'}.issubset({s['stage'] for s in stages if s['status'] == 'COMPLETED'})
            persist('FAILED' if failed else ('COMPLETED' if complete else 'PARTIAL'))
        return 1 if failed else 0

    def code_signatures(self):
        files = [self.root / 'run.sh', self.root / 'config/project.yml', *sorted((self.root / 'scripts').glob('*'))]
        return {str(p.relative_to(self.root)): sha256(p) for p in files if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--check', action='store_true')
    action.add_argument('--all', action='store_true')
    action.add_argument('--stage', choices=['rna', 'wgbs', 'liver', 'programs', 'synthesis'])
    args = parser.parse_args()
    workflow = Integration(args.root)
    if args.check:
        audit = workflow.check()
        print(audit[['id', 'status', 'source_status', 'reason']].to_string(index=False))
        return 0 if (audit.status == 'PASS').all() else 1
    return workflow.run(args.stage or 'all')


if __name__ == '__main__':
    raise SystemExit(main())

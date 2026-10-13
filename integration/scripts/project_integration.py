"""冻结数值结果的全项目整合、图件重绘及原子发布。"""
from __future__ import annotations
import hashlib
import html
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import uuid
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from plot_lifecycle import PlotRun


def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def read(p):
    return pd.read_csv(p,sep='\t' if Path(p).suffix=='.tsv' else ',',low_memory=False)


def save(df,stage,relative):
    p=stage/relative;p.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(p,index=False,sep='\t' if p.suffix=='.tsv' else ',')
    return 'result/'+relative


class FrozenSources:
    def __init__(self,stage,catalog):
        self.stage=stage;self.catalog=catalog.set_index('source_key')
    def path(self,key):
        return self.stage/self.catalog.loc[key,'path']
    def table(self,key):
        return read(self.path(key))


def freeze_existing(root,stage):
    """仅接受现役COMPLETED合同；冻结科学源，旧图不进入staging。"""
    current=root/'result'
    if not (current/'00_audit/run.json').is_file():raise ValueError('缺少现役冻结结果合同result/00_audit/run.json；须提供真实COMPLETED数值来源，不能从空result启动')
    receipt=json.loads((current/'00_audit/run.json').read_text())
    if receipt['status']!='COMPLETED':raise ValueError('当前整合结果尚未COMPLETED')
    catalog=read(current/'00_audit/tables/source_catalog.tsv')
    manifest=read(current/'00_audit/tables/figure_manifest.tsv')
    for _,r in catalog.iterrows():
        src=current/r['path']
        if not src.is_file() or digest(src)!=r['sha256']:raise ValueError('冻结源签名失配: '+r['source_key'])
        dst=stage/r['path'];dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    # plots-only 保持项目级数值源内容不变。
    for rel,sig in receipt['derived_tables'].items():
        p=current/rel
        if not p.is_file() or digest(p)!=sig:raise ValueError('项目整合表改变: '+rel)
        dst=stage/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
    for name in ('membership.tsv','term_id_map.tsv'):
        p=current/'00_audit/tables'/name
        dst=stage/'00_audit/tables'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
    return catalog,manifest,receipt


def figure_row(analysis,figure_id,source,title,meaning,distribution,conclusion,kind='heatmap'):
    stem=f'result/{analysis}/plots/{figure_id}'
    return dict(analysis=analysis,figure_id=figure_id,png=stem+'.png',pdf=stem+'.pdf',source_table=source,
                status='COMPLETED',reason='',title=title,biological_meaning=meaning,distribution=distribution,
                conclusion=conclusion,render_kind=kind)


def validate_program_terms(program,terms):
    """冻结标签按原样连接；正式ID缺失只接受显式UNMAPPED记录。"""
    mapped=program[['collection','frozen_term','term_id']].merge(
        terms[['collection','frozen_term','term_id','status']],
        on=['collection','frozen_term'],how='left',validate='many_to_one',suffixes=('_program','_map'))
    actual=mapped.term_id_program.fillna('').astype(str)
    expected=mapped.term_id_map.fillna('').astype(str)
    missing=actual.eq('')
    if mapped.status.isna().any() or actual.ne(expected).any() or (missing&mapped.status.ne('UNMAPPED')).any():
        raise ValueError('冻结基因集标签/正式ID与term_id_map合同不一致；禁止补造ID')
    return dict(program_rows=len(program),unmapped_rows=int(missing.sum()),
                unmapped_frozen_terms=int(mapped.loc[missing,'frozen_term'].nunique()),
                unmapped_labels=sorted(mapped.loc[missing,'frozen_term'].unique().tolist()),
                aggregation_key='逐字一致的collection/frozen_term；UNMAPPED不代表正式GO ID')


def integrate_tables(w,stage,src,run):
    """结果层的集合/注释整合；不重新检验，不合并 P，不把组学方向强制一致。"""
    cfg=w.cfg;new=[];p_cut=cfg['p_cutoff'];f_cut=cfg['fdr_cutoff'];lfc=cfg['lfc_cutoff']
    rows=[]
    def emit(z,relative):return save(z,stage,relative)
    def add(relative,id,z,title,meaning,distribution,conclusion,kind='heatmap'):
        source=emit(z,relative+'/tables/'+id+'.tsv')
        new.append(figure_row(relative,id,source,title,meaning,distribution,conclusion,kind))
    for sid,s in w.specs.items():
        d=w.tables[sid]
        common=dict(input_id=sid,assay=s['kind'],context=s['context'],source_status=s['source_status'],
                    quality=s.get('quality','unknown'),pairing=s.get('pairing','unknown'),
                    generation=s.get('generation','unknown'),source=s['path'])
        if s['kind']=='rna':
            z=pd.DataFrame(dict(entity_id=d.ensembl_id.astype('string').str.replace(r'\.[0-9]+$','',regex=True),symbol=d.SYMBOL,
                                effect=d.log2FoldChange,p=d.pvalue,fdr=d.padj))
            z['detected']=np.isfinite(z.p)&np.isfinite(z.effect)
            z['supported']=z.detected&(z.p<p_cut)&(z.effect.abs()>=lfc)
            z['fdr_supported']=z.detected&(z.fdr<f_cut)&(z.effect.abs()>=lfc)
            z['mapping']='gene_expression';z['unit']='RNA log2FC';z['feature']='gene'
        elif s['kind']=='wgbs':
            cov=src.table(f'result/multiomics/wgbs_rna/{s["context"]}/measured_regions_coverage_eligibility.tsv')
            selected=set(d.gene_id_clean.dropna().str.replace(r'\.[0-9]+$','',regex=True))
            z=cov.rename(columns={'gene_id_clean':'entity_id','diffMethy':'effect','Layer':'feature'}).copy()
            z['symbol']='';z['p']=np.nan;z['fdr']=np.nan
            z['detected']=z.eligible.astype(bool)&np.isfinite(z.effect)
            z['supported']=z.entity_id.isin(selected)
            z['fdr_supported']=pd.NA
            z['mapping']='direct_coverage_gene_or_promoter; support=nearest_gene_of_selected_DMR'
            z['unit']='methylation fraction delta'
            # DMR 映射与可测背景分别保留；未覆盖的 DMR 基因不能消失。
            missing=d.loc[~d.gene_id_clean.isin(set(z.entity_id)),['gene_id_clean','SYMBOL']].drop_duplicates('gene_id_clean')
            if len(missing):
                extra=missing.rename(columns={'gene_id_clean':'entity_id','SYMBOL':'symbol'}).copy()
                extra['effect']=np.nan;extra['p']=np.nan;extra['fdr']=np.nan;extra['feature']='annotation_only'
                extra['detected']=False;extra['supported']=True;extra['fdr_supported']=pd.NA
                extra['mapping']='nearest_gene_of_selected_DMR; direct_coverage_ineligible';extra['unit']='methylation fraction delta'
                z=pd.concat([z,extra],ignore_index=True)
        elif s['kind']=='atac':
            full=src.table('result/multiomics/liver/figures/complete_atac_gene_ranking_source.tsv')
            z=full.rename(columns={'gene_id_clean':'entity_id','Fold':'effect','pvalue':'p','FDR':'fdr'}).copy()
            z['symbol']='';z['detected']=np.isfinite(z.effect)&np.isfinite(z.p)
            selected=set(d.ENSEMBL.dropna().astype(str).str.replace(r'\.[0-9]+$','',regex=True))
            z['supported']=z.entity_id.isin(selected);z['fdr_supported']=z.fdr<f_cut
            z['mapping']='nearest_gene; max_absolute_region_effect; all selected links retained separately'
            z['feature']='consensus_peak_to_gene';z['unit']='ATAC log2FC'
        elif s['kind']=='metabolomics':
            z=d.rename(columns={'feature_id':'entity_id','log2FC':'effect','limma_p':'p'}).copy()
            z['symbol']='';z['detected']=np.isfinite(z.effect)&np.isfinite(z.p)
            z['supported']=z.classification.isin(['robust','exploratory']);z['fdr_supported']=z.classification.eq('robust')
            z['mapping']='no_gene_connection; original classification='+z.classification.astype(str)+'; direction='+z.direction.astype(str);z['feature']='metabolic_feature';z['unit']='metabolic log2FC'
        elif s['kind']=='celltypes':
            # 全景列出第13个测量模块；胰岛没有条件效应检验，不填0效应或伪造阴性。
            z=d.groupby('gene',sort=False)['count'].sum().reset_index().rename(columns={'gene':'entity_id'})
            z['symbol']=z.entity_id;z['effect']=np.nan;z['p']=np.nan;z['fdr']=np.nan
            z['detected']=True;z['supported']=False;z['fdr_supported']=pd.NA
            z['mapping']='islet_expression_symbol; condition_inference=NOT_RUN; original library/celltype rows retained'
            z['feature']='islet_gene_observation';z['unit']='no_condition_effect_estimate';common['assay']='scrnaseq'
        else:continue
        for k,v in common.items():z[k]=v
        keep=list(common)+['entity_id','symbol','feature','effect','unit','p','fdr','detected','supported','fdr_supported','mapping']
        for annotation_field in ('confidence_level','annotation_match_type','annotation_ambiguous'):
            if annotation_field not in z:z[annotation_field]=pd.NA
            keep.append(annotation_field)
        rows.append(z[keep])
    evidence=pd.concat(rows,ignore_index=True)
    evidence['evidence_state']=np.select([evidence.assay.eq('scrnaseq'),~evidence.detected,evidence.supported],['descriptive_only','ineligible_or_unmeasured','supported'],'measured_without_support')
    emit(evidence,'01_project_atlas/01_evidence_map/tables/entity_evidence.tsv')
    entities=evidence.groupby(['assay','context','entity_id'],dropna=False).agg(detected=('detected','max'),supported=('supported','max')).reset_index()
    entities['evidence_state']=np.select([entities.assay.eq('scrnaseq'),entities.supported&~entities.detected,entities.supported,~entities.detected],['descriptive_only','annotation_support_without_direct_gene_coverage','supported','ineligible_or_unmeasured'],'measured_without_support')
    summary=entities.groupby(['assay','context','evidence_state'],dropna=False).size().rename('entities').reset_index()
    add('01_project_atlas/01_evidence_map','measurement_support_atlas',summary,'项目测量与支持全景',
        '区分各层检测资格与按各层原规则选定的实体；RNA和代谢包含raw-P探索，WGBS/ATAC是选定区域到基因的注释，计数不可比较证据强度。','每格为去重实体数；不同组学的实体与单位不同。2-cell WGBS含WT_2cell2低覆盖，为高风险敏感性证据，不可据其稳定性外推。',
        '覆盖多层不等于样本配对或因果链；WGBS支持指选定DMR的基因注释。','atlas')
    gene=evidence[evidence.assay.isin(['rna','wgbs','atac'])].copy()
    # 同组学跨情境：先按 gene/context 聚合，WGBS 混合区域方向保持 mixed。
    keys=['assay','entity_id']
    classes=gene[keys].drop_duplicates().set_index(keys)
    # 区域组学的跨情境方向来自实际选定区域，而非任意另一层的基因覆盖均值。
    direction_rows=[gene[(gene.assay=='rna')&gene.supported][['assay','context','entity_id','effect']]]
    for sid,spec in w.specs.items():
        if spec['kind'] not in {'wgbs','atac'}:continue
        d=w.tables[sid]
        z=pd.DataFrame({'assay':spec['kind'],'context':spec['context'],
            'entity_id':d[spec['gene_column']].astype('string').str.replace(r'\.[0-9]+$','',regex=True),
            'effect':pd.to_numeric(d[spec['effect_column']],errors='coerce')})
        z=z[z.entity_id.str.match(r'^ENSMUSG[0-9]+$',na=False)]
        direction_rows.append(z)
    active=pd.concat(direction_rows,ignore_index=True)
    emit(active,'01_project_atlas/02_context_response/tables/response_direction_basis.tsv')
    active['positive']=active.effect>0;active['negative']=active.effect<0
    grouped=active.groupby(keys).agg(n_supported_contexts=('context','nunique'),
        contexts=('context',lambda v:';'.join(sorted(set(v)))),positive=('positive','max'),negative=('negative','max'))
    classes=classes.join(grouped).reset_index().rename(columns={'entity_id':'gene_id'})
    classes['n_supported_contexts']=classes.n_supported_contexts.fillna(0).astype(int)
    classes['contexts']=classes.contexts.fillna('')
    classes['response_class']=np.select([classes.n_supported_contexts.eq(0),classes.n_supported_contexts.eq(1),
        classes.positive.fillna(False)&classes.negative.fillna(False),classes.positive.fillna(False)|classes.negative.fillna(False)],
        ['no_supported_context','context_limited','repeated_opposite_or_mixed','repeated_same_sign'],'repeated_direction_unknown')
    classes=classes.drop(columns=['positive','negative'])
    emit(classes,'01_project_atlas/02_context_response/tables/gene_response_classes.tsv')
    counts=classes.groupby(['assay','response_class']).size().rename('genes').reset_index()
    add('01_project_atlas/02_context_response','context_response_classes',counts,'跨情境响应类别',
        '辨别情境限定与跨情境重复响应。','方向来自RNA候选或实际选定区域；WGBS不借用其他基因区间的均值，多区域方向冲突保留。','不同阶段和器官是独立情境，不构成个体纵向轨迹。RNA类别按原raw-P+效应探索规则，不代表跨情境BH重复。2-cell WGBS含WT_2cell2低覆盖，为高风险敏感性证据，不可据其稳定性外推。','classes')
    # 肝脏完整测量背景及支持集合；不强行统一甲基化、开放性和表达方向。
    liver=gene[gene.context=='liver']
    flags=liver.groupby(['entity_id','assay']).agg(detected=('detected','max'),supported=('supported','max')).reset_index()
    wide=flags.pivot(index='entity_id',columns='assay',values=['detected','supported']).fillna(False)
    wide.columns=['_'.join(c) for c in wide.columns];wide=wide.reset_index().rename(columns={'entity_id':'gene_id'})
    for assay in ('rna','atac','wgbs'):
        for metric in ('detected','supported'):
            if metric+'_'+assay not in wide:wide[metric+'_'+assay]=False
    wide['n_measured_layers']=wide[['detected_'+a for a in ('rna','atac','wgbs')]].sum(axis=1)
    wide['n_supported_layers']=wide[['supported_'+a for a in ('rna','atac','wgbs')]].sum(axis=1)
    wide['three_layer_support']=wide.n_supported_layers==3
    wide['three_layer_measurable_support']=wide.three_layer_support&(wide.n_measured_layers==3)
    rna_bh=gene[(gene.assay=='rna')&(gene.context=='liver')].groupby('entity_id').fdr_supported.max()
    wide['rna_fdr_supported']=wide.gene_id.map(rna_bh).fillna(False).astype(bool)
    selection_note=f"交集遵循各层原候选规则，RNA包括raw-P探索。三层共同且可测的{int(wide.three_layer_measurable_support.sum())}基因中，RNA有{int((wide.three_layer_measurable_support&wide.rna_fdr_supported).sum())}个满足原BH+效应规则；三层共同不等于三层BH显著。"
    emit(wide,'03_liver_multiomics/01_joint_discovery/tables/liver_gene_layer_support.tsv')
    support_counts=wide.groupby(['n_measured_layers','n_supported_layers']).size().rename('genes').reset_index()
    add('03_liver_multiomics/01_joint_discovery','liver_layer_support',support_counts,'肝脏三层共同基因与检测背景',
        '在共同基因身份上核对表达、开放区域和甲基化证据。','各格为真实基因数；检测与支持分别计数。',
        '同基因多层注释收敛仍不证明调控关系；未使用预设候选评分。组学间样本配对/生物单位关系unknown，结果为未配对的基因/程序注释收敛。'+selection_note,'liver_support')
    supports={a:set(wide.loc[wide['supported_'+a],'gene_id']) for a in ('rna','atac','wgbs')}
    # 冻结 Mouse_KEGG / Mouse_M5 成员、正式 term ID 对照。
    member_path=stage/'00_audit/tables/membership.tsv';term_path=stage/'00_audit/tables/term_id_map.tsv'
    if not member_path.exists():
        run.command(['Rscript','--vanilla',str(w.root/'scripts/export_frozen_memberships.R'),str(w.root/cfg['gene_sets']),str(member_path),str(term_path)],check=True)
    members=read(member_path);terms=read(term_path)
    term_lookup=terms.set_index('frozen_term').to_dict('index')
    program=[]
    def append_program(d,assay,context,method,feature,source,termcol,pcol,fcol,effectcol,genescol=None):
        for _,r in d.iterrows():
            item=term_lookup.get(r[termcol],{})
            program.append(dict(term_id=item.get('term_id',pd.NA),frozen_term=r[termcol],collection=item.get('collection','unknown'),
                assay=assay,context=context,method=method,feature=feature,effect=r.get(effectcol,np.nan),p=r[pcol],fdr=r[fcol],
                supported=pd.notna(r[fcol]) and r[fcol]<f_cut,nominal_supported=pd.notna(r[pcol]) and r[pcol]<p_cut,
                leading_genes=r.get(genescol,'') if genescol else '',source='result/'+str(src.path(source).relative_to(stage)),
                interpretation='original single-layer statistic; no combined P; layer effects retain distinct meaning'))
    rnaora=src.table('result/cross_context/rnaseq/kegg_ora.csv')
    for context,d in rnaora.groupby('context'):
        append_program(d,'rna',context,'frozen_ORA','gene_candidates','result/cross_context/rnaseq/kegg_ora.csv','pathway','p','fdr','fold_enrichment','genes')
    for key in src.catalog.index:
        if not key.endswith('_gsea.tsv'):continue
        if 'effect_product' in key or 'coupling' in key:continue # 含两层效应乘积的排序不能再次作为独立测量层计数。
        name=Path(key).stem
        if '/wgbs_rna/' in key:
            context=key.split('/wgbs_rna/')[1].split('/')[0];assay='wgbs';feature=name.rsplit('_',3)[0]
        elif 'complete_atac_' in name:context='liver';assay='atac';feature='consensus_peak_rank'
        elif 'complete_rna_' in name:context='liver';assay='rna';feature='full_cohort_statistic_rank'
        else:continue
        d=src.table(key)
        append_program(d,assay,context,'frozen_GSEA',feature,key,'pathway','pval','padj','NES','leadingEdge')
    metabolic=w.tables['metabolic_pathway_effects']
    kegg=terms[terms.collection=='KEGG'].copy()
    mappings=[]
    for _,r in metabolic.iterrows():
        hit=kegg[kegg.description==r.pathway]
        if len(hit)!=1:raise ValueError('代谢通路无法唯一对照冻结 KEGG: '+str(r.pathway))
        h=hit.iloc[0]
        mappings.append(dict(metabolic_pathway=r.pathway,term_id=h.term_id,frozen_term=h.frozen_term,status='MATCHED_EXACT_FROZEN_DESCRIPTION',membership_source=cfg['gene_sets']))
        program.append(dict(term_id=h.term_id,frozen_term=h.frozen_term,collection='KEGG',assay='metabolomics',context='liver',method='frozen_pathway_statistic',feature='metabolic_pathway',effect=r.log2FC,p=r.p,fdr=r.fdr,supported=r.fdr<f_cut,nominal_supported=r.p<p_cut,leading_genes='',source=w.specs['metabolic_pathway_effects']['path'],interpretation='pathway-level parallel evidence; metabolite-to-gene mapping MISSING'))
    emit(pd.DataFrame(mappings),'02_program_integration/01_program_map/tables/metabolic_kegg_mapping.tsv')
    program=pd.DataFrame(program)
    mapping_audit=validate_program_terms(program,terms)
    (stage/'00_audit/program_term_mapping_audit.json').write_text(json.dumps(mapping_audit,ensure_ascii=False,indent=2))
    emit(program,'02_program_integration/01_program_map/tables/program_evidence.tsv')
    # 完整单层结果都保存，图只呈现预定排名前 30，并明确未测格。
    selected=(program.groupby('frozen_term').agg(assays=('assay','nunique'),best_fdr=('fdr','min')).reset_index()
              .sort_values(['assays','best_fdr','frozen_term'],ascending=[False,True,True]).head(30).frozen_term)
    program_view=program[program.frozen_term.isin(selected)].copy()
    add('02_program_integration/01_program_map','program_evidence_atlas',program_view,'共同生物学程序的多层证据',
        '在同一冻结基因集标签上对照全部情境和测量层。','保留方法、效应、原始P与BH；灰格表示该层无此结果，浅蓝表示已有结果但BH接近1。正式term_id缺失的标签仍仅按逐字一致frozen_term聚合，数量和名单见00_audit/program_term_mapping_audit.json。',
        'ORA富集幅度、GSEA NES与代谢效应分别解释，不合并成通路活性或总P。2-cell WGBS含WT_2cell2低覆盖，为高风险敏感性证据，不可据其稳定性外推。组学间样本配对/生物单位关系unknown，结果为未配对的基因/程序注释收敛。','program')
    convergence=[]
    liver_programs={term:z for term,z in program[program.context=='liver'].groupby('frozen_term',sort=False)}
    for (collection,term_id,term),z in members.groupby(['collection','term_id','frozen_term'],dropna=False,sort=False):
        ids=set(z.gene_id);hits={a:ids&supports[a] for a in supports}
        present=[a for a in hits if hits[a]];shared=set.intersection(*hits.values())
        pairwise={a+'_'+b:hits[a]&hits[b] for a,b in [('rna','atac'),('rna','wgbs'),('atac','wgbs')]}
        pe=liver_programs.get(term,program.iloc[:0])
        supported_assays=sorted(pe.loc[pe.supported,'assay'].unique())
        convergence.append(dict(collection=collection,term_id=term_id,frozen_term=term,member_genes=len(ids),
            rna_supported_genes=len(hits['rna']),atac_supported_genes=len(hits['atac']),wgbs_supported_genes=len(hits['wgbs']),
            **{name+'_same_gene_count':len(ids) for name,ids in pairwise.items()},
            **{name+'_same_gene_genes':';'.join(sorted(ids)) for name,ids in pairwise.items()},
            n_gene_support_layers=len(present),same_gene_three_layer_count=len(shared),same_gene_three_layer_genes=';'.join(sorted(shared)),
            pathway_bh_supported_assays=';'.join(supported_assays),n_pathway_bh_layers=len(supported_assays),
            metabolic_pathway_measured='metabolomics' in set(pe.assay),metabolic_pathway_bh_supported='metabolomics' in supported_assays,
            convergence_class='same_gene_three_layers' if shared else 'same_gene_two_layers' if any(pairwise.values()) else 'same_program_different_genes' if len(present)>1 else 'single_layer_or_none'))
    convergence=pd.DataFrame(convergence)
    emit(convergence,'02_program_integration/02_gene_convergence/tables/program_gene_convergence.tsv')
    view=convergence.sort_values(['same_gene_three_layer_count','n_pathway_bh_layers','frozen_term'],ascending=[False,False,True]).head(30)
    add('02_program_integration/02_gene_convergence','same_program_and_same_gene',view,'同一程序与同一基因的收敛区别',
        '区分同基因三层、同基因两层及仅同程序不同基因；表中保留三组pairwise共同基因明细。','显示按三层共同基因数、单层BH支持层数排序的30条；完整表保留全部。',
        '基因注释收敛不等于每层通路显著；统计支持和成员支持分开报告。组学间样本配对/生物单位关系unknown，结果为未配对的基因/程序注释收敛。'+selection_note,'convergence')
    fourth=convergence[convergence.metabolic_pathway_measured].copy()
    metabolic_stats=program[program.assay.eq('metabolomics')][['frozen_term','effect','p','fdr','source']].rename(columns={'effect':'metabolic_log2FC','p':'metabolic_p','fdr':'metabolic_fdr','source':'metabolic_source'})
    fourth=fourth.merge(metabolic_stats,on='frozen_term',how='left',validate='one_to_one')
    add('03_liver_multiomics/01_joint_discovery','liver_four_layer_programs',fourth,'肝脏三层基因与第四层代谢程序',
        '把已有的冻结代谢通路统计量与三层基因背景并列。','左侧为6条冻结代谢通路统计对应的三层注释基因数，右侧为冻结代谢log2FC并标原BH；代谢源方向仍未独立核验。',
        '第四层是通路并列证据；代谢物到具体基因的连接保持MISSING。组学间样本配对/生物单位关系unknown，结果为未配对的基因/程序注释收敛。'+selection_note,'convergence')
    # 留一组学检查：集合依赖描述，绝不作为独立验证或重新推断显著性。
    ablations=[]
    for omitted in ('none','rna','atac','wgbs'):
        remaining=[a for a in supports if a!=omitted]
        common=set.intersection(*(supports[a] for a in remaining))
        for g in sorted(common):ablations.append(dict(omitted_assay=omitted,gene_id=g,remaining_layers=';'.join(remaining),support_class='all_remaining_layers',independent_validation=False))
    dep=pd.DataFrame(ablations);emit(dep,'06_synthesis/01_layer_dependency/tables/gene_leave_one_assay_out.tsv')
    dep_summary=dep.groupby('omitted_assay').size().rename('genes').reset_index()
    add('06_synthesis/01_layer_dependency','leave_one_assay_out',dep_summary,'候选集合对各测量层的依赖',
        '观察移除一种组学后共同支持集合如何改变。','固定所有原始统计结果，只重新计算集合交汇。','这是证据依赖性描述，不是独立队列验证。组学间样本配对/生物单位关系unknown，结果为未配对的基因/程序注释收敛。'+selection_note,'dependency')
    program_dep=[]
    for omitted in ('none','rna','atac','wgbs','metabolomics'):
        for term,z in program[program.context=='liver'].groupby('frozen_term',sort=False):
            measured=set(z.assay)-{omitted};supported=set(z.loc[z.supported,'assay'])-{omitted}
            program_dep.append(dict(omitted_assay=omitted,frozen_term=term,measured_layers=';'.join(sorted(measured)),bh_supported_layers=';'.join(sorted(supported)),n_measured=len(measured),n_bh_supported=len(supported)))
    emit(pd.DataFrame(program_dep),'06_synthesis/01_layer_dependency/tables/program_leave_one_assay_out.tsv')
    # 胰岛文库的描述性定位：count/cells 与每库细胞型CPM，均不做处理组检验。
    cell=w.tables['islet_celltypes'].copy();cell['count_per_cell']=cell['count']/cell.cells.replace(0,np.nan)
    totals=cell.groupby(['library_id','cell_type'])['count'].transform('sum')
    cell['cpm']=1e6*cell['count']/totals.replace(0,np.nan)
    symbol_map=gene[gene.assay=='rna'][['entity_id','symbol']].dropna().drop_duplicates()
    candidates=set(wide.loc[wide.n_supported_layers>=2,'gene_id'])
    top_programs=set(convergence.sort_values(['n_pathway_bh_layers','same_gene_three_layer_count'],ascending=False).head(20).frozen_term)
    relevant=set(members.loc[members.frozen_term.isin(top_programs),'gene_id'])|candidates
    mapped=symbol_map[symbol_map.entity_id.isin(relevant)].groupby('symbol').entity_id.agg(lambda v:';'.join(sorted(set(v)))).rename('gene_ids')
    cell=cell.merge(mapped,left_on='gene',right_index=True,how='inner',validate='many_to_one')
    cell['interpretation']='one_library_per_condition; islet descriptive context, not liver cell origin'
    emit(cell,'04_celltype_context/01_global_candidates/tables/islet_candidate_expression.tsv')
    top_symbols=(cell.groupby('gene').cpm.max().sort_values(ascending=False).head(40).index)
    add('04_celltype_context/01_global_candidates','islet_candidate_context',cell[cell.gene.isin(top_symbols)],'候选在胰岛细胞类型中的描述定位',
        '补充候选与共同程序的器官细胞类型背景。','两个文库分别展示细胞型CPM，显示最高表达的40个候选或程序成员；完整表保留其余。',
        '每组一库，不能做处理组推断，也不能把胰岛来源当作肝脏细胞来源。','celltype')
    return pd.DataFrame(new)


def plot_python(stage,manifest,cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import textwrap
    from figure_style import apply_style,save_figure,COLORS,PROFILE,DIVERGING
    apply_style()
    assay_colors=PROFILE['palettes']['assay']
    python_kinds={'atlas','liver_support','program','celltype','convergence','classes','dependency'}
    for _,r in manifest[(manifest.status=='COMPLETED')&(manifest.get('render_kind','').isin(python_kinds)|manifest.figure_id.isin(['rna_context_summary','liver_gene_effects']))].iterrows():
        if pd.isna(r.get('render_kind')) and r.figure_id not in {'rna_context_summary','liver_gene_effects'}:continue
        d=read(stage/str(r.source_table).removeprefix('result/'))
        kind=r.get('render_kind','');fig,ax=plt.subplots();height=145
        if r.figure_id=='rna_context_summary':
            accent=PROFILE['palettes']['accent'];x=np.arange(len(d));ax.bar(x-.18,d.nominal,.36,color=accent['warm'],label='Raw P + effect');ax.bar(x+.18,d.fdr,.36,color=accent['purple'],label='BH + effect')
            ax.set_xticks(x,d.context,rotation=35,ha='right');ax.set_ylabel('RNA candidate genes');ax.legend(frameon=False)
            fig.subplots_adjust(left=.13,right=.97,bottom=.24,top=.93)
        elif r.figure_id=='liver_gene_effects':
            z=d.dropna(subset=['rna_log2FC','atac_median_effect'])
            for present,label,c in [(False,'ATAC + RNA',PROFILE['palettes']['neutral']['dark']),(True,'+ WGBS region',assay_colors['wgbs'])]:
                t=z[z.has_wgbs_region==present];ax.scatter(t.rna_log2FC,t.atac_median_effect,s=8,alpha=.6,c=c,label=label,rasterized=True)
            ax.axhline(0,color=COLORS['NS'],lw=.6);ax.axvline(0,color=COLORS['NS'],lw=.6)
            ax.set(xlabel='RNA log2FC (PFOS - WT)',ylabel='Median annotated ATAC effect (PFOS - WT)');ax.legend(frameon=False)
            fig.subplots_adjust(left=.16,right=.97,bottom=.16,top=.96)
        elif kind in {'atlas','liver_support','program','celltype'}:
            if kind=='atlas':
                d['measurement']=d.assay+' / '+d.context
                m=d.pivot(index='measurement',columns='evidence_state',values='entities').fillna(0)
                values=np.log10(m+1);height=max(150,len(m)*7+55);label='log10(entities + 1)'
            elif kind=='liver_support':
                m=d.pivot(index='n_measured_layers',columns='n_supported_layers',values='genes').fillna(0)
                values=np.log10(m+1);label='log10(genes + 1)'
            elif kind=='program':
                d['measurement']=d.context+' / '+d.assay+' / '+d.feature+' / '+d.method
                d['support']=-np.log10(pd.to_numeric(d.fdr,errors='coerce').clip(lower=1e-10)).clip(upper=10)
                signed=d.method.eq('frozen_GSEA')|d.feature.eq('metabolic_pathway')
                d['direction']=np.where(~signed,'Direction not defined for ORA',np.where(d.effect>0,'Positive',np.where(d.effect<0,'Negative','Zero')))
                terms=list(dict.fromkeys(d.frozen_term));measurements=list(dict.fromkeys(d.measurement));height=240
                x=d.measurement.map({v:i for i,v in enumerate(measurements)});y=d.frozen_term.map({v:i for i,v in enumerate(terms)})
                colors=d.direction.map({'Positive':PROFILE['palettes']['de']['up'],'Negative':PROFILE['palettes']['de']['down'],'Zero':PROFILE['palettes']['de']['ns'],'Direction not defined for ORA':PROFILE['palettes']['neutral']['dark']})
                ax.scatter(x,y,s=5+4*d.support.fillna(0),c=colors,alpha=np.where(d.support.fillna(0)>-np.log10(.05),.9,.35),linewidths=0)
                ax.set_xticks(range(len(measurements)),[str(v).replace(' / rna / ',' | RNA | ').replace(' / atac / ',' | ATAC | ').replace(' / wgbs / ',' | WGBS | ').replace(' / metabolomics / ',' | Metabolism | ').replace('gene_candidates / frozen_ORA','ORA').replace('consensus_peak_rank / frozen_GSEA','peak GSEA').replace('full_cohort_statistic_rank / frozen_GSEA','rank GSEA').replace('cgi_related / frozen_GSEA','CGI GSEA').replace('gene_body / frozen_GSEA','body GSEA').replace('promoter / frozen_GSEA','promoter GSEA').replace('metabolic_pathway / frozen_pathway_statistic','pathway statistic') for v in measurements],rotation=90)
                ax.set_yticks(range(len(terms)),[textwrap.fill(str(v).replace('GOBP_','').replace('_',' '),45) for v in terms]);ax.invert_yaxis()
                from matplotlib.lines import Line2D
                legend=[Line2D([0],[0],marker='o',linestyle='',color=c,label=label,markersize=5) for label,c in [('Positive',PROFILE['palettes']['de']['up']),('Negative',PROFILE['palettes']['de']['down']),('ORA (direction undefined)',PROFILE['palettes']['neutral']['dark'])]]
                ax.legend(handles=legend,frameon=False,loc='upper left',bbox_to_anchor=(1.01,1));fig.subplots_adjust(left=.36,right=.73,bottom=.34,top=.98)
                save_figure(fig,stage/str(r.pdf).removeprefix('result/'),height_mm=min(240,height));plt.close(fig);continue
            else:
                d['measurement']=d.library_id+' / '+d.cell_type
                m=d.pivot_table(index='gene',columns='measurement',values='cpm',aggfunc='sum')
                values=np.log10(m+1);label='log10(cell-type CPM + 1)';height=240
            from matplotlib.colors import LinearSegmentedColormap
            cmap=LinearSegmentedColormap.from_list('pfos_sequential',PROFILE['palettes']['sequential']);cmap.set_bad(COLORS['NS'])
            im=ax.imshow(np.ma.masked_invalid(values.to_numpy(float)),aspect='auto',interpolation='nearest',cmap=cmap)
            labels=list(m.columns)
            if kind=='atlas':
                labels=[{'annotation_support_without_direct_gene_coverage':'Annotation only','descriptive_only':'Descriptive only','ineligible_or_unmeasured':'Ineligible / unmeasured','measured_without_support':'Measured, unsupported','supported':'Layer-selected'}[v] for v in m.columns]
            elif kind=='program':
                labels=[str(v).replace(' / rna / ',' | RNA | ').replace(' / atac / ',' | ATAC | ').replace(' / wgbs / ',' | WGBS | ').replace(' / metabolomics / ',' | Metabolism | ').replace('gene_candidates / frozen_ORA','ORA').replace('consensus_peak_rank / frozen_GSEA','peak GSEA').replace('full_cohort_statistic_rank / frozen_GSEA','rank GSEA').replace('cgi_related / frozen_GSEA','CGI GSEA').replace('gene_body / frozen_GSEA','body GSEA').replace('promoter / frozen_GSEA','promoter GSEA').replace('metabolic_pathway / frozen_pathway_statistic','pathway statistic') for v in m.columns]
            ax.set_xticks(range(len(m.columns)),labels,rotation=90)
            ax.set_yticks(range(len(m)),[textwrap.fill(str(v).replace('GOBP_','').replace('_',' '),45) for v in m.index])
            if kind in {'atlas','liver_support'}:
                for i in range(len(m)):
                    for j in range(len(m.columns)):ax.text(j,i,str(int(m.iloc[i,j])),ha='center',va='center',fontsize=7,color='white' if values.iloc[i,j]>.55*np.nanmax(values.to_numpy()) else PROFILE['palettes']['neutral']['dark'])
            fig.colorbar(im,ax=ax,shrink=.55,label=label,pad=.025)
            if kind=='liver_support':ax.set(xlabel='Layer-selected layers',ylabel='Measured layers')
            left=.36 if kind in {'atlas','program'} else .24 if kind=='celltype' else .17
            fig.subplots_adjust(left=left,right=.91,bottom=.34 if kind=='program' else .23,top=.98)
        elif kind=='convergence':
            z=d.iloc[::-1];y=np.arange(len(z));fourth=r.figure_id=='liver_four_layer_programs'
            if fourth:
                fig.clear();ax=fig.add_axes([.45,.17,.27,.74]);meta=fig.add_axes([.80,.17,.17,.74],sharey=ax)
            for k,(a,c) in enumerate([(a,assay_colors[a]) for a in ('rna','atac','wgbs')]):
                ax.barh(y+(k-1)*.23,z[a+'_supported_genes'],height=.22,color=c,label=a.upper())
            ax.scatter(z.same_gene_three_layer_count,y,marker='D',s=15,color=PROFILE['palettes']['neutral']['dark'],label='same gene / 3 layers',zorder=3)
            ax.set_yticks(y,[textwrap.fill(t.replace('GOBP_','').replace('_',' '),44) for t in z.frozen_term]);ax.set_xlabel('Layer-selected gene annotations')
            height=min(240,max(165,len(z)*7+50))
            if not fourth:
                for pair,marker,offset in [('rna_atac','o',-.15),('rna_wgbs','s',0),('atac_wgbs','^',.15)]:
                    ax.scatter(z[pair+'_same_gene_count'],y+offset,marker=marker,s=23,facecolors='none',edgecolors=PROFILE['palettes']['neutral']['dark'],linewidths=.7,label='same gene '+pair.upper().replace('_',' - '),zorder=4)
            handles,labels=ax.get_legend_handles_labels()
            if fourth:
                bh=z.metabolic_fdr<cfg['fdr_cutoff']
                meta.axvline(0,color=COLORS['NS'],lw=.6)
                meta.scatter(z.metabolic_log2FC,y,s=22,c=[PROFILE['palettes']['accent']['purple'] if v else 'white' for v in bh],edgecolors=PROFILE['palettes']['neutral']['dark'],linewidths=.6)
                for effect,yy,fdr in zip(z.metabolic_log2FC,y,z.metabolic_fdr):meta.annotate('BH '+format(fdr,'.2g'),(effect,yy),xytext=(0,6),textcoords='offset points',ha='center',fontsize=7)
                meta.tick_params(axis='y',left=False,labelleft=False);meta.spines[['left','top','right']].set_visible(False)
                meta.set_xlabel('Metabolic log2FC\n(source unverified)');meta.margins(x=.25)
            else:fig.subplots_adjust(left=.49,right=.97,bottom=.13,top=.93)
            fig.legend(handles,labels,frameon=False,loc='upper right',bbox_to_anchor=(.98,.995),ncol=2)
        elif kind=='classes':
            m=d.pivot(index='response_class',columns='assay',values='genes').fillna(0);y=np.arange(len(m));left=np.zeros(len(m))
            for a in m.columns:ax.barh(y,m[a],left=left,color=assay_colors[a],label=a.upper());left+=m[a].to_numpy()
            ax.set_yticks(y,[v.replace('_',' ') for v in m.index]);ax.set_xlabel('Genes');ax.legend(frameon=False);fig.subplots_adjust(left=.38,right=.97,bottom=.14,top=.95)
        elif kind=='dependency':
            ax.bar(d.omitted_assay,d.genes,color=[PROFILE['palettes']['neutral']['dark'] if v=='none' else assay_colors[v] for v in d.omitted_assay]);ax.set(xlabel='Omitted assay (none = all three)',ylabel='Genes layer-selected in all remaining layers')
            for i,v in enumerate(d.genes):ax.text(i,v,str(v),ha='center',va='bottom',fontsize=7)
            fig.subplots_adjust(left=.17,right=.97,bottom=.17,top=.92)
        else:raise ValueError('未声明绘图类型: '+str(kind))
        save_figure(fig,stage/str(r.pdf).removeprefix('result/'),height_mm=min(240,height));plt.close(fig)


def render_frozen(w,stage,manifest,run):
    from concurrent.futures import ThreadPoolExecutor
    env=os.environ.copy();env.update(PFOS_PLOTS_ONLY='1',FROZEN_RESULT_DIR=str(stage),MAGICK_THREAD_LIMIT='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    commands=[('rna',['Rscript','--vanilla',str(w.root/'scripts/rebuild_rna.R'),str(w.root)]),
              ('liver',['Rscript','--vanilla',str(w.root/'scripts/rebuild_liver_multiomics.R'),str(w.root)])]
    commands += [('wgbs_rna_'+c,['Rscript','--vanilla',str(w.root/'scripts/rebuild_wgbs_rna.R'),str(w.root),c]) for c in w.cfg['joint_rna_context']]
    commands += [('cross_'+c,['Rscript','--vanilla',str(w.root/'scripts/replot_cross_wgbs.R'),str(w.root),c]) for c in w.cfg['cross_wgbs_pairs']]
    if manifest.get('render_kind',pd.Series('',index=manifest.index)).astype(str).str.startswith('deep_').any():
        commands += [('deepening',['Rscript','--vanilla',str(w.root/'scripts/deepening_analysis.R'),str(w.root),str(stage),'plot'])]
    def job(item):
        name,cmd=item
        with (stage/'logs'/f'{name}.log').open('w') as f:
            result=run.command(cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
        if result.returncode:raise RuntimeError(name+'重绘失败，见 '+str(stage/'logs'/f'{name}.log'))
        print('RENDER COMPLETED',name,flush=True)
    with ThreadPoolExecutor(max_workers=8) as pool:
        for future in [pool.submit(job,item) for item in commands]:future.result()
    plot_python(stage,manifest,w.cfg)
    scratch=stage/'00_audit/cross_render'
    if scratch.exists():shutil.rmtree(scratch)


def verify_media(stage,manifest,catalog):
    from PIL import Image
    completed=manifest[manifest.status=='COMPLETED']
    declared=set()
    for _,r in completed.iterrows():
        for col in ('png','pdf','source_table'):
            p=stage/str(r[col]).removeprefix('result/')
            if not p.is_file() or not p.stat().st_size:raise ValueError('图/源缺失: '+str(p))
            if col!='source_table':declared.add(str(p.relative_to(stage)))
        png=stage/str(r.png).removeprefix('result/')
        with Image.open(png) as im:
            dpi=im.info.get('dpi',(0,0))
            if im.mode!='RGB' or abs(im.width-2161)>1 or min(dpi)<299:raise ValueError('成图规格未通过: '+str(png))
            height_mm=im.height/dpi[1]*25.4
            if height_mm>240.5:raise ValueError(f'成图高度超过240 mm: {png} ({height_mm:.1f} mm)')
        if Path(r.png).stem!=Path(r.pdf).stem:raise ValueError('PNG/PDF不同名')
    actual={str(p.relative_to(stage)) for p in stage.rglob('*') if p.suffix in {'.png','.pdf'}}
    if actual!=declared:raise ValueError('媒体清单覆盖失配: '+str(actual^declared))
    for _,r in catalog.iterrows():
        if digest(stage/r.path)!=r.sha256:raise ValueError('科学表格内容改变: '+r.source_key)
    return dict(figure_pairs=len(completed),statuses=manifest.status.value_counts().to_dict(),source_tables_unchanged=len(catalog),unlisted_media=0)


def semantic_qa(root,stage,manifest,cfg):
    profile_path=(root/cfg['figure_profile']).resolve()
    profile=json.loads(profile_path.read_text())
    checks=[]
    def check(name,passed,evidence):
        checks.append(dict(check=name,status='PASS' if passed else 'FAIL',evidence=evidence))
        if not passed:raise ValueError('语义图形QA未通过: '+name+'; '+evidence)
    check('group_palette',profile['palettes']['group']=={'WT':'#5F8FB3','PFOS':'#CE7C8E'},json.dumps(profile['palettes']['group'],ensure_ascii=False))
    check('figure_profile_identity',digest(profile_path)=='e22153228d724956351f4964f014e3720877d8591dbf4b5ecf61ab8a8183c227',digest(profile_path))
    check('fixed_semantic_scales',profile['scales']['row_zscore']==[-2.5,2.5] and profile['scales']['nes']==[-3,3] and profile['scales']['ora_minus_log10_p']==[0,10],json.dumps(profile['scales'],ensure_ascii=False))
    completed=manifest[manifest.status=='COMPLETED']
    check('three_element_interpretation',completed[['biological_meaning','distribution','conclusion']].astype(str).apply(lambda x:x.str.len().gt(0)).all().all(),'all completed rows contain biology, distribution, conclusion')
    deep=read(stage/'00_audit/tables/deepening_manifest.tsv')
    deep_keys=set(zip(deep.analysis.astype(str),deep.figure_id.astype(str)))
    original=sum((str(r.analysis),str(r.figure_id)) not in deep_keys for _,r in completed.iterrows())
    check('original_plot_coverage',original==916 and len(completed)==916+len(deep),f'{original} original and {len(deep)} deepening pairs; {len(completed)} total')
    check('deepening_manifest',len(deep)>0 and deep.status.eq('COMPLETED').all() and len(deep_keys)==len(deep),f'{len(deep)} declared, unique, completed new figure pairs')
    tests=read(stage/'03_liver_multiomics/03_background_calibration/tables/stratified_permutation_tests.tsv')
    check('permutation_contract',len(tests)==3 and tests.B.eq(10000).all() and tests.seed.eq(1).all() and tests.empirical_p.between(1/10001,1).all(),'three prespecified tests; B=10000; seed=1')
    elig=read(stage/'03_liver_multiomics/03_background_calibration/tables/layer_inference_eligibility.tsv')
    check('wgbs_inference_boundary',elig.loc[elig.layer.eq('WGBS'),'eligible'].astype(str).str.lower().eq('false').all(),'WGBS selected-DMR/coverage eligibility only; no three-layer P')
    gsea=read(stage/'02_program_integration/03_context_programs/tables/full_statistic_fgsea.tsv')
    check('gsea_family_rows',set(gsea.collection)=={'Mouse_MH','Mouse_KEGG','Mouse_Reactome'} and gsea.context.nunique()==7,f'{len(gsea)} rows; 7 contexts x 3 collections')
    der=read(stage/'06_synthesis/02_posttranscriptional_patterns/tables/der_adapter_reconstruction_audit.tsv')
    check('der_adapter_rebuilt',(der.legacy_adapter_selected_genes==0).all() and (der.reconstructed_nominal_candidate_genes>0).all() and (der.reconstructed_candidate_bh_genes<=der.reconstructed_nominal_candidate_genes).all(),'legacy DER adapter all zero; current nominal-region candidates rebuilt; BH-supported candidates are a subset')
    cohort=read(stage/'06_synthesis/02_posttranscriptional_patterns/tables/cohort_signature_audit.tsv')
    check('cohort_signature_boundary',(cohort.overlap_allowed.astype(str).str.lower()==cohort.same_sample_signature.astype(str).str.lower()).all(),'overlap_allowed exactly follows sample signature equality')
    method_cohort=read(stage/'06_synthesis/02_posttranscriptional_patterns/tables/method_cohort_signature_audit.tsv')
    check('method_cohort_signatures',len(method_cohort)==35 and method_cohort.status.eq('PASS').all() and method_cohort.same_selected_signature.astype(str).str.lower().eq('true').all(),'35 original DAS/DEU/APA/DER/IS metadata files match the selected cohort and are independently signed')
    method_summary=read(stage/'06_synthesis/02_posttranscriptional_patterns/tables/method_availability_summary.tsv')
    check('posttx_statistical_scope',{'candidate_genes','candidate_bh_genes','reported_fdr_scope'}<=set(method_summary.columns) and method_summary.candidate_bh_genes.fillna(0).le(method_summary.candidate_genes).all() and method_summary.loc[method_summary.method.isin(['DAS','APA','DER']),'reported_fdr_scope'].eq('event/region FDR summary; not gene FDR').all(),'candidate and BH-supported counts separated; event/region minima remain summaries rather than gene inference')
    loci=read(stage/'03_liver_multiomics/04_shared_loci/tables/five_candidate_locus_evidence.tsv')
    per_gene=loci.groupby('gene_id').assay.agg(lambda x:set(x))
    check('five_shared_loci',len(per_gene)==5 and per_gene.map(lambda x:{'RNA','ATAC','WGBS'}<=x).all(),f'{len(per_gene)} genes; each contains RNA/ATAC/WGBS')
    save(pd.DataFrame(checks),stage,'00_audit/tables/semantic_qa.tsv')
    return {r['check']:r['status'] for r in checks}


def reports(stage,manifest,stats):
    program=read(stage/'02_program_integration/01_program_map/tables/program_evidence.tsv')
    mapping_audit=validate_program_terms(program,read(stage/'00_audit/tables/term_id_map.tsv'))
    (stage/'00_audit/program_term_mapping_audit.json').write_text(json.dumps(mapping_audit,ensure_ascii=False,indent=2))
    names={'01_project_atlas':'项目全景','02_program_integration':'共同生物学程序','03_liver_multiomics':'肝脏多层证据','04_celltype_context':'胰岛细胞类型背景','05_pairwise_support':'两两比较支撑','06_synthesis':'证据依赖与综合'}
    css='body{font-family:sans-serif;max-width:1200px;margin:auto;padding:28px;color:#263238}img{max-width:100%;height:auto}figure{padding:20px;border-bottom:1px solid #ddd}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:5px}a{color:#24618b}'
    import fitz
    manifest['pdf_pages']=0
    manifest['png_preview_page']=''
    for index,row in manifest[manifest.status=='COMPLETED'].iterrows():
        with fitz.open(stage/str(row.pdf).removeprefix('result/')) as document:manifest.at[index,'pdf_pages']=len(document)
        manifest.at[index,'png_preview_page']=1
    manifest.to_csv(stage/'00_audit/tables/figure_manifest.tsv',sep='\t',index=False)
    nav=[];readme=['# PFOS mouse 项目级联合分析','','阅读顺序：项目全景 → 共同生物学程序 → 肝脏多层证据 → 胰岛细胞背景 → 两两比较支撑 → 证据依赖。',
        '',f'当前共 {stats["figure_pairs"]} 对图件；原有 916 对图件全部重绘并保留，另新增 {sum(manifest.get("render_kind",pd.Series("",index=manifest.index)).astype(str).str.startswith("deep_"))} 对深化图。771 个既有科学源表保持冻结；深化分析从已完成的单组学输入生成独立新表，包括预设 fgsea、方法特异 ORA 与分层置换，不重跑单组学差异模型。',
        '', '多页位点PDF保留全部页面，图册标明页数；PNG仅为第1页预览，完整位点请打开PDF。', '', '**边界：**样本配对 unknown；2cell WGBS 保留低覆盖风险；testis 为研究者确认的 F0；胰岛每组一库，仅描述。通路支持与同基因多层支持分别报告，代谢物到基因连接保持 MISSING。',
        '', '新项目须先准备真实COMPLETED数值来源及00_audit冻结合同；不支持从空result生成统计全流程。fixture与CONTROL/TREATMENT烟测只验证合同和绘图调用，不等于真实数据整合。',
        '', 'Layer-selected表示按各层原规则选定的候选，RNA与代谢包含raw-P探索；不等于各层BH显著，也不表示证据强度相等。',
        '', f'冻结基因集按原样标签聚合；{mapping_audit["unmapped_frozen_terms"]}个GO_BP标签、{mapping_audit["unmapped_rows"]}行未映射当前正式term_id，保留UNMAPPED，不补造ID。具体标签见[映射审计](result/00_audit/program_term_mapping_audit.json)。',
        '', '`metabolic_pathway_measured`仅表示上游冻结通路统计量存在；不表示通路被直接测量，也不是本轮本地富集计算。', '', '运行：`bash run.sh --check` 核验输入；`bash run.sh --integrate` 生成整合与全部图件；`bash run.sh --plots-only` 从冻结数值源重绘。路径基于本模块根，启动目录不影响解析。','']
    for folder,title in names.items():
        rows=manifest[manifest.png.fillna('').str.startswith('result/'+folder+'/') | manifest.source_table.fillna('').str.startswith('result/'+folder+'/')]
        pieces=[f'<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>{title}</title><style>{css}</style><body><a href="../index.html">返回全项目</a><h1>{title}</h1>']
        for _,r in rows.iterrows():
            pieces.append('<figure id="'+html.escape(str(r.analysis)+'_'+str(r.figure_id))+'"><h2>'+html.escape(str(r.title))+'</h2>')
            if r.status=='COMPLETED':
                media=os.path.relpath(stage/str(r.png).removeprefix('result/'),stage/folder)
                pdf=os.path.relpath(stage/str(r.pdf).removeprefix('result/'),stage/folder)
                pieces.append(f'<a href="{html.escape(pdf)}">PDF</a><img loading="lazy" src="{html.escape(media)}" alt="{html.escape(str(r.title))}">')
                if int(r.pdf_pages)>1:pieces.append(f'<p>PDF含{int(r.pdf_pages)}页，打开查看完整位点图；PNG仅为第1页预览。</p>')
            else:pieces.append('<p>'+html.escape(str(r.status)+': '+str(r.reason))+'</p>')
            for label,col in [('生物学意义','biological_meaning'),('数据分布特征','distribution'),('最终结论','conclusion')]:
                pieces.append('<p><b>'+label+'：</b>'+html.escape(str(r[col]))+'</p>')
            if str(r.source_table):pieces.append('<a href="'+html.escape(os.path.relpath(stage/str(r.source_table).removeprefix('result/'),stage/folder))+'">数值来源</a>')
            pieces.append('</figure>')
        pieces.append('</body></html>');dest=stage/folder/'index.html';dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text('\n'.join(pieces))
        nav.append(f'<li><a href="{folder}/index.html">{title}</a>：{sum(rows.status=="COMPLETED")} 对图件</li>')
        readme += [f'## {title}', '',f'[阅读全部图表](result/{folder}/index.html)。','']
        # 单一 README 覆盖全部图件解释，不建立精选图副本。
        for _,r in rows[(rows.status=='COMPLETED') & rows.get('render_kind',pd.Series('',index=rows.index)).fillna('').ne('')].iterrows():
            readme += [f'### {r.title}', '',f'[PNG]({r.png}) · [PDF]({r.pdf}) · [源表]({r.source_table})',
                       '',f'生物学意义：{r.biological_meaning}',f'数据分布特征：{r.distribution}',f'最终结论：{r.conclusion}','']
    previews=''.join('<figure><h2>'+html.escape(str(r.title))+'</h2><a href="'+html.escape(str(r.pdf).removeprefix('result/'))+'"><img src="'+html.escape(str(r.png).removeprefix('result/'))+'" alt="'+html.escape(str(r.title))+'"></a><p>'+html.escape(str(r.conclusion))+'</p></figure>' for _,r in manifest[manifest.get('render_kind',pd.Series('',index=manifest.index)).fillna('').ne('')].iterrows())
    index=f'<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>PFOS 项目联合分析</title><style>{css}</style><body><h1>PFOS：全局发现与肝脏深入</h1><p>现有科学结果冻结；样本未配对，不作跨组学样本相关或因果推断。</p><ol>'+''.join(nav)+'</ol>'+previews+'<p><a href="00_audit/tables/figure_manifest.tsv">全部图件与状态</a> · <a href="00_audit/run.json">本次执行记录</a></p></body></html>'
    (stage/'index.html').write_text(index)
    return '\n'.join(readme)


def run_project(w,integrate):
    with PlotRun(w.root) as run:
        return _run_project(w, integrate, run)


def _run_project(w,integrate,run):
    root=w.root;current=root/'result';work=run.stage('.result-');stage=work/'result';stage.mkdir()
    (stage/'00_audit/tables').mkdir(parents=True);(stage/'logs').mkdir()
    receipt={'time':datetime.now(timezone.utc).isoformat(),'status':'STARTED','action':'integrate' if integrate else 'plots-only',
             'interpretation':'冻结结果的项目级整合与重绘；不是单组学模型或富集检验重跑'}
    audit_path=stage/'00_audit/run.json'
    try:
        audit=w.check()
        if not (audit.status=='PASS').all():raise ValueError('输入资格未通过: '+audit.loc[audit.status!='PASS',['id','reason']].to_json(orient='records',force_ascii=False))
        signatures={r['id']:r['sha256'] for r in w.audit};code=w.code_signatures()
        print('INPUT IDENTITY PASS; freezing existing numeric sources',flush=True)
        catalog,manifest,parent=freeze_existing(root,stage)
        previous_signatures=parent.get('input_signatures',{})
        if any(signatures.get(key)!=value for key,value in previous_signatures.items()):raise ValueError('冻结数值源与当前原始输入签名不一致，禁止混用陈旧统计')
        print('FROZEN SOURCES',len(catalog),'existing figure pairs',sum(manifest.status=='COMPLETED'),flush=True)
        # 保留原生产记录为独立历史执行证据；不修改其状态、路径或 hash。
        history=root.parents[3]/'maintenance/pfos_unification_20261012/integration'
        history.mkdir(parents=True,exist_ok=True)
        parent_bytes=json.dumps(parent,ensure_ascii=False,indent=2)
        parent_hash=hashlib.sha256(parent_bytes.encode()).hexdigest()
        historical=history/('producer_receipt_'+parent_hash+'.json')
        if not historical.exists():historical.write_text(parent_bytes)
        previous_evidence=history/'previous_execution_evidence'
        for p in current.rglob('*'):
            if p.is_file() and (p.suffix=='.txt' or p.name in {'analysis.html','input_audit.tsv'}):
                dest=previous_evidence/p.relative_to(current);dest.parent.mkdir(parents=True,exist_ok=True)
                if not dest.exists():shutil.copy2(p,dest)
        receipt.update(parent_receipt_sha256=digest(historical),input_signatures=signatures,code_signatures=code)
        save(catalog,stage,'00_audit/tables/source_catalog.tsv');save(audit,stage,'00_audit/tables/input_audit.tsv')
        old=manifest[manifest.analysis.isin(['cross_context/rnaseq','cross_context_rna','liver_atac_rna_wgbs','multiomics/liver','wgbs_2cell_liver','wgbs_2cell_sperm','wgbs_sperm_liver','wgbs_rna_2cell','wgbs_rna_liver','wgbs_rna_sperm'])].copy()
        if integrate:
            new=integrate_tables(w,stage,FrozenSources(stage,catalog),run);manifest=pd.concat([old,new],ignore_index=True)
            print('PROJECT INTEGRATION TABLES READY; new figure families',len(new),flush=True)
            deep_log=stage/'logs/deepening_analyse.log'
            deep_env=os.environ.copy();deep_env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
            with deep_log.open('w') as stream:
                deep=run.command(['Rscript','--vanilla',str(root/'scripts/deepening_analysis.R'),str(root),str(stage),'analyse'],env=deep_env,stdout=stream,stderr=subprocess.STDOUT)
            if deep.returncode:raise RuntimeError('深化分析失败，见 '+str(deep_log))
            deep_manifest=read(stage/'00_audit/tables/deepening_manifest.tsv')
            manifest=pd.concat([manifest,deep_manifest],ignore_index=True)
            profile_path=(root/w.cfg['figure_profile']).resolve()
            deep_receipt={'status':'COMPLETED','seed':w.cfg['deepening']['seed'],'permutations':w.cfg['deepening']['permutations'],
                'program_collections':w.cfg['deepening']['program_collections'],'fgsea':{'eps':0,'parallel':'BiocParallel::SerialParam','bh_family':'context_x_collection'},
                'rna_atac_permutation':{'B':w.cfg['deepening']['permutations'],'strata':['RNA expression quintile','ATAC peak/link count','region category']},
                'wgbs_inference':'selected_DMR_and_coverage_eligibility_only; no_three_layer_P',
                'cohort_rule':'selected6 and full10 compared only when exact sample signature matches; 35 DAS/DEU/APA/DER/IS method metadata files must match selected6 and are content-signed',
                'statistical_scope':'DAS/APA candidates use upstream nominal event P plus effect rule; DER candidates use upstream nominal region P; event/region minimum P/FDR are gene-mapping summaries and never relabeled as gene P/FDR; method-specific ORA remains exploratory',
                'figure_profile_sha256':digest(profile_path),
                'contracts':{'analysis':{'version':'sci-analysis-constructive-2026-10-11-v6','sha256':'49ed9faf6900a1865b64f02cc573bef2f9064671323378d562ac8c4759084ead'},
                    'figure':{'version':'sci-source-data-composition-2026-10-12-v6','sha256':'4d32daedab1f82b731a83c7ab66b9d46d86c8a8f587dfb13de664e4f0a7663cc'}}}
            (stage/'00_audit/deepening_receipt.json').write_text(json.dumps(deep_receipt,ensure_ascii=False,indent=2))
        manifest=manifest.fillna('')
        for col in ('conclusion','distribution','biological_meaning'):
            manifest[col]=manifest[col].str.replace('睾丸代际unknown','睾丸F0',regex=False).str.replace('睾丸代际未知','睾丸F0',regex=False)
        save(manifest,stage,'00_audit/tables/figure_manifest.tsv')
        audit_path.write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
        render_frozen(w,stage,manifest,run)
        stats=verify_media(stage,manifest,catalog)
        stats['semantic_qa']=semantic_qa(root,stage,manifest,w.cfg)
        w.check()
        if {r['id']:r['sha256'] for r in w.audit}!=signatures or w.code_signatures()!=code:raise ValueError('运行期间输入或现役源码改变，禁止发布')
        readme=reports(stage,manifest,stats)
        derived={str(p.relative_to(stage)):digest(p) for p in stage.rglob('*.tsv') if str(p.relative_to(stage)) not in set(catalog.path) and p.name not in {'figure_manifest.tsv','input_audit.tsv','source_catalog.tsv'} and not p.name.startswith('render_')}
        receipt.update(status='COMPLETED',finished_at=datetime.now(timezone.utc).isoformat(),verification=stats,derived_tables=derived,
                       outputs={str(p.relative_to(stage)):digest(p) for p in stage.rglob('*') if p.is_file() and p!=audit_path})
        audit_path.write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
        run.publish(stage, current)
        w.report.write_text(readme)
        (history/'verification.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2))
        refinement=root.parents[3]/'maintenance/pfos_visual_refinement_20261012/integration'
        refinement.mkdir(parents=True,exist_ok=True)
        for source,name in [(current/'00_audit/deepening_receipt.json','deepening_receipt.json'),
                            (current/'00_audit/tables/semantic_qa.tsv','semantic_qa.tsv'),
                            (current/'00_audit/run.json','run.json')]:
            if source.is_file():shutil.copy2(source,refinement/name)
        (refinement/'figure_profile_sha256.txt').write_text(digest((root/w.cfg['figure_profile']).resolve())+'  projects/pfos/config/figure_style.json\n')
        print(json.dumps(dict(status='COMPLETED',**stats),ensure_ascii=False),flush=True)
        return 0
    except BaseException as e:
        receipt.update(status='FAILED',reason=str(e),finished_at=datetime.now(timezone.utc).isoformat())
        run._failure(e)
        print('FAILED; 当前 result 保留。'+str(e),file=__import__('sys').stderr,flush=True)
        return 1

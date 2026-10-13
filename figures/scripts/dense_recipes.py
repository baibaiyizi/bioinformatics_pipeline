#!/usr/bin/env python3
"""高级证据图：仅渲染完整、已冻结的统计来源表。"""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, textwrap
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[1] / '.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd
from figure_style import PROFILE, COLORS, DIVERGING, save_figure
from render import require, checked_numeric


def plot_program_evidence_matrix(data):
    require(data,['context','collection','pathway','NES','padj_context_collection'],['context','collection','pathway'])
    checked_numeric(data,['NES','padj_context_collection'])
    cmap=LinearSegmentedColormap.from_list('effect',DIVERGING);cmap.set_bad(PROFILE['scales']['missing_color'])
    contexts=[c for c in PROFILE['palettes']['context'] if c in set(data.context)]
    if set(contexts)!=set(data.context):raise ValueError('context必须在项目profile登记')
    collections=list(dict.fromkeys(data.collection));row_counts=[len(data.loc[data.collection.eq(c),'pathway'].unique()) for c in collections]
    fig,axes=plt.subplots(len(collections),1,figsize=(183/25.4,(45+9*sum(row_counts))/25.4),squeeze=False,gridspec_kw={'height_ratios':row_counts})
    fig.subplots_adjust(left=.45,right=.87,bottom=.07,top=.93,hspace=.3)
    for ax,collection in zip(axes[:,0],collections):
        sub=data.loc[data.collection.eq(collection)];matrix=sub.pivot(index='pathway',columns='context',values='NES').reindex(columns=contexts)
        if len(matrix)>30:raise ValueError('每集合超过30行时请先生成有记录的科学问题子集；本函数不默默删行')
        artist=ax.imshow(matrix,cmap=cmap,norm=TwoSlopeNorm(0,*PROFILE['scales']['nes']),aspect='auto')
        q=sub.pivot(index='pathway',columns='context',values='padj_context_collection').reindex(index=matrix.index,columns=contexts)
        row,col=np.where(q.to_numpy()<.05);ax.scatter(col,row,s=25,facecolors='none',edgecolors='#374151',linewidths=.7)
        ax.set_yticks(range(len(matrix)),[textwrap.fill(str(v).replace('_',' '),34) for v in matrix.index]);ax.set_xticks(range(len(contexts)),contexts,rotation=35,ha='right');ax.set_title(str(collection))
        if ax is not axes[-1,0]:ax.tick_params(labelbottom=False)
    cb=fig.colorbar(artist,cax=fig.add_axes([.91,.15,.025,.65]),label='NES',ticks=[-3,0,3]);cb.ax.set_yticklabels(['≤−3','0','≥3'])
    fig.suptitle('Frozen preranked programs; circles = BH < 0.05\nBH family = context × collection; declared display subset',fontsize=8,y=.986)
    return fig


def plot_background_calibration(null,tests):
    require(null,['permutation','test','count'],['permutation','test']);require(tests,['test','observed','empirical_p','B','seed'],['test'])
    checked_numeric(null,['count']);checked_numeric(tests,['observed','empirical_p','B','seed'])
    colors={'Open-Up':COLORS['PFOS'],'Closed-Down':COLORS['WT'],'Concordant':PROFILE['palettes']['accent']['purple']}
    fig,axes=plt.subplots(len(tests),1,figsize=(183/25.4,190/25.4),layout='constrained',squeeze=False)
    for ax,row in zip(axes[:,0],tests.itertuples()):
        if row.test not in colors:raise ValueError(f'未登记置换问题: {row.test}')
        values=null.loc[null.test.eq(row.test),'count']
        if len(values)!=int(row.B):raise ValueError('null行数与冻结B不符')
        ax.hist(values,bins=np.arange(values.min(),values.max()+2)-.5,color=colors[row.test],alpha=.75,label='Stratified-label null',rwidth=.9)
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.axvline(row.observed,color='#374151',linestyle='--',linewidth=1,label=f'Observed = {int(row.observed)}')
        ax.set(xlabel='Gene count',ylabel='Permutations',title=f'{row.test}; empirical P={row.empirical_p:.3g}; B={int(row.B)}; seed={int(row.seed)}');ax.legend(frameon=False)
    fig.suptitle('Conditional background calibration; no cross-omics pairing assumed',fontsize=8)
    return fig


def plot_shared_locus_tracks(data):
    require(data,['gene_id','symbol','assay','chr','start','end','effect','unit','direction'])
    checked_numeric(data,['start','end','effect'])
    if (data.start>data.end).any():raise ValueError('start>end')
    genes=list(dict.fromkeys(data.gene_id));height=max(140,min(240,45*len(genes)))
    fig,axes=plt.subplots(len(genes),1,figsize=(183/25.4,height/25.4),layout='constrained',squeeze=False)
    lanes={'WGBS':0,'ATAC':1,'RNA':2};colors={'positive':COLORS['PFOS'],'negative':COLORS['WT'],'zero':PROFILE['palettes']['de']['ns']}
    for ax,gene in zip(axes[:,0],genes):
        sub=data.loc[data.gene_id.eq(gene)]
        if sub.chr.nunique()!=1:raise ValueError('同一个gene涉及多个chr，需先核对映射')
        lane_labels=[]
        for assay in lanes:
            units=sub.loc[sub.assay.eq(assay),'unit'].unique()
            if len(units)>1:raise ValueError('同一locus/assay包含多量纲，需独立分面')
            unit=str(units[0]).removeprefix(assay+' ') if len(units) else 'not measured'
            lane_labels.append(assay+'\n'+textwrap.fill(unit,17))
        for row in sub.itertuples():
            if row.assay not in lanes or row.direction not in colors:raise ValueError('未登记assay/direction')
            y=lanes[row.assay];color=colors[row.direction]
            ax.plot([row.start,row.end],[y,y],color=color,linewidth=4,solid_capstyle='round')
            if row.start==row.end:ax.scatter(row.start,y,s=25,c=color)
            ax.annotate(f'{row.effect:+.3g}',((row.start+row.end)/2,y),xytext=(0,6),textcoords='offset points',ha='center',fontsize=7)
        ax.set(yticks=list(lanes.values()),yticklabels=lane_labels,ylim=(-.6,2.7),xlabel='GRCm39 coordinate',title=f'{sub.symbol.iloc[0]} | {sub.chr.iloc[0]}');ax.margins(x=.1)
        ax.ticklabel_format(axis='x',style='plain',useOffset=False)
    fig.suptitle('Shared coordinates; layer-specific units; selected evidence ≠ causal regulation',fontsize=8)
    return fig


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['program_nes','calibration','locus'],required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--tests',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    data=pd.read_csv(a.input,sep='\t');secondary=None
    if a.kind=='calibration':
        if a.tests is None:p.error('calibration要求--tests，不能重算或猜测经验P')
        secondary=pd.read_csv(a.tests,sep='\t');fig=plot_background_calibration(data,secondary)
    elif a.kind=='program_nes':fig=plot_program_evidence_matrix(data)
    else:fig=plot_shared_locus_tracks(data)
    a.output.mkdir(parents=True,exist_ok=True);outputs=save_figure(fig,a.output/'figure.pdf',height_mm=fig.get_figheight()*25.4);plt.close(fig)
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    shutil.copy2(a.input,a.output/'source.tsv');inputs={'source.tsv':sha(a.input)}
    if a.tests:shutil.copy2(a.tests,a.output/'tests.tsv');inputs['tests.tsv']=sha(a.tests)
    record={'status':'COMPLETED','kind':a.kind,'input_sha256':inputs,'rows':len(data),'script_sha256':sha(__file__),'profile_sha256':sha(Path(__file__).resolve().parent/'../config/figure_style.json'),'statistics_recomputed':False,'outputs':{k:sha(v) for k,v in outputs.items()}}
    (a.output/'figure.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'status':'COMPLETED','kind':a.kind,'rows':len(data)}))
if __name__=='__main__':main()

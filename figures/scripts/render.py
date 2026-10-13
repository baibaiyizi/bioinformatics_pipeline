#!/usr/bin/env python3
"""使用冻结表绘图；本脚本不拟合模型、不检验组差异。"""
from __future__ import annotations
import argparse, hashlib, json, os, shutil
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[1] / '.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from figure_style import PROFILE, COLORS, DIVERGING, configure_profile, save_figure, apply_style


def require(frame, columns, keys=()):
    missing=set(columns)-set(frame)
    if missing: raise ValueError(f'源表缺列: {sorted(missing)}')
    if keys and frame.duplicated(list(keys)).any(): raise ValueError(f'绘图键重复: {keys}')


def checked_numeric(frame, columns):
    for column in columns:
        value=pd.to_numeric(frame[column], errors='raise')
        if np.isinf(value).any(): raise ValueError(f'{column} 含Inf')
        frame[column]=value


def category_color(label):
    label=str(label)
    for name in ('group','cell_type','clusters','databases','context','assay'):
        mapping=PROFILE['palettes'][name]
        if label in mapping: return mapping[label]
    condition=label.split('_',1)[0]
    if condition in COLORS:return COLORS[condition]
    raise ValueError(f'未注册类别颜色: {label}；请在项目profile显式登记')


def draw(frame, kind, value_label):
    diverging=LinearSegmentedColormap.from_list('semantic_diverging',DIVERGING)
    diverging.set_bad(PROFILE['scales']['missing_color'])
    if kind in ('ora','marker'):
        fig=plt.figure(figsize=(183/25.4,(190 if kind=='ora' else 160)/25.4))
        grid=fig.add_gridspec(2,2,width_ratios=[1,.045],height_ratios=[1,.10])
        ax=fig.add_subplot(grid[0,0]);colorbar_ax=fig.add_subplot(grid[0,1]);point_legend_ax=fig.add_subplot(grid[1,:]);point_legend_ax.set_axis_off()
    else:
        fig, ax=plt.subplots(figsize=(183/25.4,140/25.4),layout='constrained')
    if kind=='volcano':
        require(frame,['feature','effect','pvalue','evidence'],['feature']);checked_numeric(frame,['effect','pvalue'])
        if not frame.pvalue.between(0,1).all():raise ValueError('P必须在[0,1]')
        if not set(frame.evidence)<= {'robust','exploratory','not_sig'}:raise ValueError('evidence必须为上游冻结的robust/exploratory/not_sig')
        for evidence,marker in [('not_sig','o'),('exploratory','^'),('robust','o')]:
            sub=frame.loc[frame.evidence.eq(evidence)];colors=np.where(sub.effect.lt(0),COLORS['WT'],COLORS['PFOS']) if evidence!='not_sig' else '#B7BEC6'
            ax.scatter(sub.effect,-np.log10(sub.pvalue.clip(lower=np.finfo(float).tiny)),s=5,c=colors,marker=marker,alpha=.75,linewidths=0,rasterized=True)
        from matplotlib.lines import Line2D
        handles=[Line2D([],[],linestyle='',marker='o',color=c,label=l,markersize=4) for l,c in [('WT / negative',COLORS['WT']),('PFOS / positive',COLORS['PFOS']),('not_sig','#B7BEC6')]]
        handles += [Line2D([],[],linestyle='',marker=m,color='#374151',label=l,markersize=4) for l,m in [('BH + FC supported','o'),('Nominal-P candidate','^')]]
        ax.axvline(0,color='#374151',linewidth=.6);ax.set(xlabel=value_label,ylabel='−log10(nominal P)');ax.legend(handles=handles,frameon=False)
    elif kind=='ora':
        require(frame,['term','count','query_size','pvalue'],['term']);checked_numeric(frame,['count','query_size','pvalue'])
        if not frame.query_size.gt(0).all() or not frame['count'].between(0,frame.query_size).all():raise ValueError('Count/query_size分母不合法')
        if not frame.pvalue.between(0,1).all():raise ValueError('P必须在[0,1]')
        color=-np.log10(frame.pvalue.clip(lower=np.finfo(float).tiny));norm=Normalize(*PROFILE['scales']['ora_minus_log10_p'],clip=True)
        artist=ax.scatter(frame['count']/frame.query_size,np.arange(len(frame)),s=2*frame['count'],c=color,cmap='viridis',norm=norm,edgecolors='#374151',linewidths=.3)
        ax.set_yticks(range(len(frame)),frame.term);ax.set(xlabel='Gene ratio (Count / tested query)',ylabel='');ax.margins(x=.12,y=.04);fig.colorbar(artist,cax=colorbar_ax,label='−log10(nominal P); ≥10 capped')
        point_legend_ax.legend(handles=[ax.scatter([],[],s=2*n,color='#B7BEC6',label=str(n)) for n in (5,20,50)],title='Overlap genes (linear point area)',frameon=False,ncol=3,loc='center')
    elif kind=='marker':
        require(frame,['feature','group','fraction','value'] if 'feature' in frame else ['gene','group','fraction','value'])
        if 'feature' not in frame:frame=frame.rename(columns={'gene':'feature'})
        require(frame,[],['feature','group']);checked_numeric(frame,['fraction','value'])
        if not frame.fraction.between(-1e-12,1+1e-12).all():raise ValueError('检出比例必须在[0,1]；只容忍1e-12数值误差')
        frame["fraction"]=frame.fraction.clip(0,1)
        features=list(dict.fromkeys(frame.feature));groups=sorted(set(frame.group.astype(str)),key=lambda x:(not x.isdigit(),int(x) if x.isdigit() else x))
        xmap={v:i for i,v in enumerate(features)};ymap={v:i for i,v in enumerate(groups)}
        artist=ax.scatter(frame.feature.map(xmap),frame.group.astype(str).map(ymap),s=200*frame.fraction,c=frame.value,cmap=diverging,norm=TwoSlopeNorm(0,*PROFILE['scales']['marker_zscore']),linewidths=.25,edgecolors='#374151')
        ax.set_xticks(range(len(features)),features,rotation=60,ha='right');ax.set_yticks(range(len(groups)),groups);ax.set(xlabel='Marker gene',ylabel='Cluster');cb=fig.colorbar(artist,cax=colorbar_ax,label=value_label,ticks=[-2.5,0,2.5]);cb.ax.set_yticklabels(['≤−2.5','0','≥2.5'])
        point_legend_ax.legend(handles=[ax.scatter([],[],s=200*f,color='#B7BEC6',label=f'{100*f:g}%') for f in (.25,.5,1)],title='Detected cells (linear point area)',frameon=False,ncol=3,loc='center')
    elif kind in ('heatmap','evidence_matrix'):
        if kind=='evidence_matrix':
            require(frame,['feature','library_id','cell_type','condition','value'],['feature','library_id','cell_type'])
            frame=frame.copy();frame['sample']=frame.library_id.astype(str)+' / '+frame.cell_type.astype(str)
            conditions=frame.set_index('sample').condition.to_dict();celltypes=frame.set_index('sample').cell_type.to_dict()
            order=sorted(frame['sample'].unique(),key=lambda x:(list(COLORS).index(conditions[x]),list(PROFILE['palettes']['cell_type']).index(celltypes[x])))
        else:require(frame,['feature','sample','value'],['feature','sample']);order=list(dict.fromkeys(frame['sample']))
        checked_numeric(frame,['value']);matrix=frame.pivot(index='feature',columns='sample',values='value').reindex(columns=order)
        if kind=='evidence_matrix':
            plt.close(fig)
            fig=plt.figure(figsize=(183/25.4,175/25.4),layout='constrained')
            gs=fig.add_gridspec(4,1,height_ratios=[.38,.14,.14,4])
            header=fig.add_subplot(gs[0]);header.axis('off')
            header.text(.0,.95,'Frozen descriptive programs; one library per condition',va='top',transform=header.transAxes,fontsize=8)
            header.legend(handles=[Patch(facecolor=COLORS[k],label=k) for k in PROFILE['palettes']['group'] if k in set(conditions.values())],frameon=False,ncol=2,loc='lower left')
            from matplotlib.colors import to_rgb
            for row,name,mapping in [(1,'Condition',conditions),(2,'Cell type',celltypes)]:
                strip=fig.add_subplot(gs[row]);strip.imshow([[to_rgb(category_color(mapping[label])) for label in order]],aspect='auto');strip.set_xticks([]);strip.set_yticks([]);strip.set_ylabel(name,rotation=0,ha='right',va='center')
                for spine in strip.spines.values():spine.set_visible(False)
            ax=fig.add_subplot(gs[3])
        artist=ax.imshow(matrix.to_numpy(),cmap=diverging,norm=TwoSlopeNorm(0,*PROFILE['scales']['row_zscore']),aspect='auto')
        ax.set_yticks(range(len(matrix)),matrix.index);ax.set_xticks(range(len(order)),order,rotation=65,ha='right');cb=fig.colorbar(artist,ax=ax,label=value_label,ticks=[-2.5,0,2.5]);cb.ax.set_yticklabels(['≤−2.5','0','≥2.5'])
        if kind=='heatmap':ax.set_title('Frozen descriptive values; missing cells remain masked')
    elif kind=='effect':
        require(frame,['feature','effect'],['feature']);checked_numeric(frame,['effect'])
        y=np.arange(len(frame));ax.axvline(0,color='#374151',linewidth=.6)
        colors=np.where(frame.effect.ge(0),COLORS['PFOS'],COLORS['WT']);ax.hlines(y,0,frame.effect,color=colors,linewidth=.8);ax.scatter(frame.effect,y,c=colors,s=20)
        if {'lower','upper'} <= set(frame):
            checked_numeric(frame,['lower','upper']);ax.hlines(y,frame.lower,frame.upper,color='#374151',linewidth=.8)
        ax.set_yticks(y,frame.feature);ax.set_xlabel(value_label)
    elif kind=='embedding':
        require(frame,['x','y','group']);checked_numeric(frame,['x','y'])
        for group,sub in frame.groupby('group',sort=True):ax.scatter(sub.x,sub.y,s=3,c=category_color(group),label=group,rasterized=True)
        ax.set(xlabel='Embedding 1',ylabel='Embedding 2');ax.legend(frameon=False)
    elif kind=='composition':
        require(frame,['library_id','cell_type','proportion'],['library_id','cell_type']);checked_numeric(frame,['proportion'])
        matrix=frame.pivot(index='library_id',columns='cell_type',values='proportion')
        if matrix.isna().any().any():raise ValueError('组成网格缺项，需在源表明确zero或missing')
        if not np.allclose(matrix.sum(axis=1),1):raise ValueError('各库比例之和必须为1')
        order=sorted(matrix.index,key=lambda x:(list(COLORS).index(str(x).split('_',1)[0]),str(x)))
        matrix.loc[order].plot.bar(stacked=True,ax=ax,color={c:category_color(c) for c in matrix},width=.75);ax.set(ylabel=value_label,ylim=(0,1),xlabel='Library');ax.legend(frameon=False,bbox_to_anchor=(1.02,1))
    else:raise ValueError(kind)
    if kind in ('ora','marker'):fig.tight_layout(pad=1.3)
    return fig


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--kind',choices=['volcano','ora','marker','heatmap','effect','embedding','composition','evidence_matrix'],required=True);p.add_argument('--value-label',required=True);p.add_argument('--profile',type=Path);a=p.parse_args()
    if a.profile:configure_profile(a.profile)
    apply_style()
    frame=pd.read_csv(a.input,sep='\t');fig=draw(frame,a.kind,a.value_label);a.output.mkdir(parents=True,exist_ok=True)
    outputs=save_figure(fig,a.output/'figure.pdf',height_mm=fig.get_figheight()*25.4);plt.close(fig);shutil.copy2(a.input,a.output/'source.tsv')
    sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
    receipt={'status':'COMPLETED','kind':a.kind,'rows':len(frame),'input_sha256':sha(a.input),'script_sha256':sha(__file__),'profile_id':PROFILE['profile_id'],'profile_sha256':sha(a.profile or Path(__file__).resolve().parent/'../config/figure_style.json'),'value_label':a.value_label,'statistics_recomputed':False,'outputs':{k:{'path':Path(v).name,'sha256':sha(v)} for k,v in outputs.items()}}
    (a.output/'figure.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'status':'COMPLETED','kind':a.kind,'rows':len(frame)}))
if __name__=='__main__':main()

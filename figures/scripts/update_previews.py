#!/usr/bin/env python3
"""从本轮真实成图回执生成离线图册，并核验原生PDF和PNG。"""
from pathlib import Path
import hashlib
import html
import json
import fitz
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
ENTRIES = [
 ('r_violin','无轮廓小提琴：中位数与四分位摘要','比较甲基化方向及区域类别对应的RNA效应分布。','冻结3173行；填色透明度0.30，无轮廓或箱体；点为中位数，0.25mm线为Q1–Q3。','描述基因层分布；四分位线不是置信区间或动物级检验。'),
 ('python_volcano','代谢特征效应与统计分布','沿用已冻结代谢特征分类，观察效应与统计支持。','26,359特征全部显示；方向用蓝粉、not_sig灰，robust圆形、exploratory三角；P与BH来自原模型。','未知身份仍称代谢特征；不新增统计检验。'),
 ('python_ora','Marker ORA：Python原生','比较标注基因集的命中比例与现有统计支持。','完整已显示24条；GeneRatio有真实query分母，面积线性Count、viridis固定0–10。','用于marker注释，不是PFOS条件推断。'),
 ('r_ora','Marker ORA：R原生','同一源表的R版式候选。','与Python同一真实ORA表和变量语义；项目选一种实现冻结。','候选版式执行不代表论文上游复现。'),
 ('python_marker','Marker表达：Python原生','联合呈现表达强度与捕获细胞检出比例。','预声明marker白名单实际存在的9基因×13群；颜色为原标准化表达，面积为真实检出比例。','群注释证据，不补条件显著性。'),
 ('r_marker','Marker表达：R原生','同一marker证据的R版式候选。','固定名称群序0–12、零中心±2.5色阶、线性检出面积。','与ORA点图区分统计语义。'),
 ('python_program_matrix','程序×文库×细胞类型矩阵','在一张图中对齐多程序和细胞背景。','9程序×2文库×7类型；现有visual_zscore不重算，注释条固定身份色。','每条件一文库，仅描述定位。'),
 ('r_program_matrix','ComplexHeatmap多层证据矩阵','增加捕获细胞数，帮助判断程序位置和信息量。','同一126行程序源表，另带实际捕获细胞数柱和两层注释。','细胞数不是独立动物数，不创建组差异检验。'),
 ('python_composition','捕获细胞组成','检查各文库捕获的细胞类型分布。','实际14个library×celltype比例；类别固定，完整组成网格不补零。','捕获组成不是动物级组织丰度。'),
 ('python_program_nes','跨背景排序程序矩阵','对齐七背景、三个基因集集合的既有排序富集方向。','实际PFOS整合dry-run：完整9777行保留；每集合按跨context最小既有BH取10 pathways、全部contexts；NES±3，圈示族内BH<0.05。','显示子集不代表全局无偏发现；不同背景独立估计，无跨阶段同一动物配对。'),
 ('python_calibration','完整测试背景的分层置换校准','将三类观察重叠放到条件背景零分布中评估。','实际PFOS整合dry-run冻结30000条null、3项既有经验P与B=10000；不重算置换或P。','校准现有背景偏倚，不证明跨组学调控或因果。'),
 ('python_locus','候选位点共享坐标多轨','核对RNA、ATAC和WGBS是否确实在同一坐标邻域提供证据。','实际PFOS整合dry-run17条、5候选；共享GRCm39坐标，图上逐条保留本层效应单位。','selected evidence，未检验层与零效应不混同；不是样本配对或因果轨迹。'),
]

def main():
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    figures=[];issues=[]
    lines=["<!doctype html><html lang='zh'><meta charset='utf-8'><title>生信绘图模板真实预览</title>","<style>body{font:16px sans-serif;max-width:1100px;margin:auto;padding:28px;color:#263445}img{max-width:100%;border:1px solid #e4e7eb}figure{margin:40px 0}p{line-height:1.6}</style>","<h1>生信绘图模板：真实来源预览</h1><p>全部来自冻结PFOS源表；跨背景程序、条件背景校准、共享坐标三类示例来自新整合分析dry-run。这里验证绘图接口，不代表上游论文复现或新增生物学结论。一个项目同一类图固定一个变体。<a href='../README.md'>使用说明</a></p>"]
    for key,title,meaning,distribution,conclusion in ENTRIES:
        folder=ROOT/'previews'/key;receipt=json.loads((folder/'figure.json').read_text())
        if receipt['status']!='COMPLETED':raise ValueError(f'{key}未完成')
        png,pdf=folder/'figure.png',folder/'figure.pdf'
        with Image.open(png) as im:
            image={'mode':im.mode,'size':list(im.size),'dpi':list(im.info.get('dpi',()))}
            if im.mode!='RGB' or abs(im.width-183/25.4*300)>2:issues.append(f'{key}:PNG规格')
        with fitz.open(pdf) as doc:
            font_records=[]
            for font in doc.get_page_fonts(0):
                extracted=doc.extract_font(font[0]);font_records.append({'name':font[3],'type':font[2],'embedded':bool(extracted[3])})
            text=doc[0].get_text();width_mm=doc[0].rect.width*25.4/72
            # Cairo rounds its PDF media box to whole points (518pt = 182.74mm).
            if len(doc)!=1 or abs(width_mm-183)>.4 or not font_records or not all(f['embedded'] for f in font_records) or len(text)<30:issues.append(f'{key}:PDF规格/字体/文本')
        dense=key in {'python_program_nes','python_calibration','python_locus'}
        figure={'recipe':key,'title':title,'validation_state':'example_executed','data_kind':'real_PFOS_frozen_dry_run_source' if dense else 'real_PFOS_frozen_source','png':str(png.relative_to(ROOT/'previews')),'pdf':str(pdf.relative_to(ROOT/'previews')),'png_sha256':sha(png),'pdf_sha256':sha(pdf),'profile_sha256':receipt['profile_sha256'],'statistics_recomputed':False,'media':image,'pdf_width_mm':width_mm,'pdf_fonts':font_records,'pdf_text_characters':len(text)}
        figures.append(figure)
        lines.append(f"<figure id='{key}'><h2>{html.escape(title)}</h2><img src='{key}/figure.png' alt='{html.escape(title)}'><figcaption><p>生物学意义：{html.escape(meaning)}</p><p>数据分布特征：{html.escape(distribution)}</p><p>最终结论：{html.escape(conclusion)}</p><a href='{key}/figure.pdf'>原生PDF</a> · <a href='{key}/figure.json'>执行回执</a></figcaption></figure>")
    (ROOT/'previews/index.html').write_text('\n'.join(lines)+"</html>\n")
    manifest={'status':'COMPLETED' if not issues else 'FAILED','figure_pairs':len(figures),'figures':figures}
    (ROOT/'previews/manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    (ROOT/'previews/validation.json').write_text(json.dumps({'status':'PASS' if not issues else 'FAIL','figure_pairs':len(figures),'issues':issues,'scientific_model_rerun':False,'manual_visual_review':'另见workspace维护报告；自动规格通过不等同视觉或投稿通过'},ensure_ascii=False,indent=2)+'\n')
    if issues:raise ValueError(issues)
    print(f'PASS {len(figures)} preview pairs; RGB300dpi, 183mm PDF, embedded fonts and text.')

if __name__=='__main__':main()

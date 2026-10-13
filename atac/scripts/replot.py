#!/usr/bin/env python3
"""从冻结数值重绘；临时结果通过检查后替换，历史执行签名不改写。"""
import argparse, csv, hashlib, html, json, os, re, shutil, subprocess, sys, tempfile
from pathlib import Path
from datetime import datetime, timezone
from plot_lifecycle import PlotRun

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path);p.add_argument('--assay', choices=('wgbs','atac'),required=True)
    p.add_argument('--steps', default='');p.add_argument('--check', action='store_true')
    a=p.parse_args(); root=a.root.resolve(); result=root/'result'
    with PlotRun(root) as run:
        script=Path(__file__).resolve().with_name('replot.R')
        if not result.is_dir() or not script.is_file(): raise SystemExit('MISSING: result 或重绘源码不存在')
        tables=sorted(x for x in result.rglob('*') if x.is_file() and x.relative_to(result).parts[0] not in ('00_audit','logs') and x.suffix.lower() not in ('.pdf','.png','.html','.md'))
        pdfs=sorted(result.rglob('*.pdf'))
        if not tables or not pdfs: raise SystemExit('MISSING: 没有可追溯数值或原图清单')
        samples_file=Path(os.environ['WGBS_SAMPLES' if a.assay=='wgbs' else 'ATAC_SAMPLES'])
        extra=[samples_file]
        if a.assay=='wgbs':
            extra.append(root/'04methykit/CpG.RData')
            extra.extend(x for folder in ('01rawdata/qc_results_GRCm39','04methykit') for x in (root/folder).rglob('*') if x.is_file() and x.suffix in ('.csv','.json','.txt','.tsv'))
        else:
            extra.extend(sorted((root/'.cache/rmd').glob('*.rds')))
            extra.extend(x for folder in ('01rawdata/qc_results','02mapping/qc','03peaks','04counts','05bigwig','06qc_report','metadata') for x in (root/folder).rglob('*') if x.is_file() and x.suffix.lower() in ('.csv','.tsv','.txt','.json','.log','.narrowpeak','.bw','.gz'))
            with samples_file.open() as handle:
                rows=list(csv.DictReader(handle,delimiter='\t'))
            first={}
            for row in rows:first.setdefault(row['group'],row['sample'])
            for sample in first.values():
                for ext in ('.final.bam','.final.bam.bai'):
                    f=root/'02mapping/bam'/(sample+ext)
                    if f.is_file():extra.append(f)
        extra=sorted(set(extra))
        if any(not f.is_file() for f in extra):raise SystemExit('MISSING: 样本表或保存对象')
        signatures={str(x.relative_to(root)):sha(x) for x in [*tables,*extra]}
        if a.check:
            print(json.dumps(dict(status='READY',module=root.name,assay=a.assay,figures=len(pdfs),numerical_inputs=len(signatures)),ensure_ascii=False));return
        stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        work=run.stage('.plots-only-'); stage=work/'result'; stage.mkdir()
        source_paths=sorted(x for base in (root/'rmd', script.parent) for x in base.glob('*') if x.is_file() and x.suffix in ('.R','.Rmd','.py','.sh'))
        code_signatures={str(x.relative_to(root)) if x.is_relative_to(root) else str(Path('..')/x.relative_to(root.parent)):sha(x) for x in source_paths}
        profile = (root/'config/figure_style.json').resolve()
        profile_relative = os.path.relpath(profile, root)
        code_signatures[profile_relative] = sha(profile)
        receipt=dict(figure_profile=dict(path=profile_relative,sha256=sha(profile),profile_id=json.loads(profile.read_text())['profile_id']),code_signatures=code_signatures,status='STARTED',purpose='统一语义配色与冻结数值重绘',started_at=stamp,assay=a.assay,upstream='not_run',differential='not_run',enrichment_tests='not_run',inputs=signatures)
        try:
            for f in sorted(result.rglob('*')):
                if not f.is_file():continue
                target=stage/f.relative_to(result);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,target)
            audit=stage/'00_audit';(audit/'tables').mkdir(parents=True,exist_ok=True)
            (audit/'plots_only_run.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
            baseline=[]
            for f in pdfs:
                rel=f.relative_to(result)
                baseline.append(dict(figure_id=rel.as_posix(),old_sha256=sha(f)))
            with (audit/'tables/figure_baseline.tsv').open('w') as f:
                w=csv.DictWriter(f,fieldnames=['figure_id','old_sha256'],delimiter='\t');w.writeheader();w.writerows(baseline)
            env=os.environ.copy();env.update(MAGICK_THREAD_LIMIT='1',PFOS_PLOT_SOURCE_ROOT=str(root),PFOS_PLOT_RESULT_ROOT=str(stage),PFOS_PLOT_STEPS=a.steps)
            log=stage/'logs/plots_only.log';log.parent.mkdir(exist_ok=True)
            with log.open('w') as f: run.command(['Rscript',str(script),str(root),str(stage)],env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
            # Existing scientific bytes must be unchanged even if their addresses change.
            for f in tables:
                target=stage/f.relative_to(result)
                if not target.is_file() or sha(target)!=signatures[str(f.relative_to(root))]:raise RuntimeError('科学结果被改写: '+str(f))
            for rel,digest in signatures.items():
                if sha(root/rel)!=digest:raise RuntimeError('输入在绘图期间改变: '+rel)
            for rel,digest in code_signatures.items():
                if sha(root/rel)!=digest:raise RuntimeError('绘图源码在运行期间改变: '+rel)
            from PIL import Image
            import fitz
            statuses=[]
            for row in baseline:
                f=stage/row['figure_id'];png=f.with_suffix('.png')
                if not f.is_file() or not png.is_file() or f.stat().st_size==0 or png.stat().st_size==0:raise RuntimeError('图件配对缺失: '+row['figure_id'])
                changed=sha(f)!=row['old_sha256']
                if changed:
                    with Image.open(png) as im:
                        if im.mode!='RGB' or abs(im.width-round(183/25.4*300))>2:raise RuntimeError('PNG 输出规格错误: '+str(png))
                        dpi=im.info.get('dpi',(0,0))
                        if any(abs(x-300)>.2 for x in dpi):raise RuntimeError('PNG DPI 错误: '+str(png))
                        im.verify()
                    with fitz.open(f) as doc:
                        if abs(doc[0].rect.width-183/25.4*72)>1:raise RuntimeError('PDF 宽度错误: '+str(f))
                if not changed and not a.steps:
                    raise RuntimeError('完整纯绘图遗漏原图: '+row['figure_id'])
                reason='' if changed else ('MISSING_FROZEN_MODEL' if re.search(r'_(pmd|umr_lmr)_model[.]pdf$',str(f)) else 'NOT_SELECTED_OR_DISPLAY_SOURCE_UNVERIFIED')
                statuses.append(dict(figure_id=row['figure_id'],status='RERENDERED' if changed else 'PRESERVED_NOT_RERENDERED',reason=reason,pdf_sha256=sha(f),png_sha256=sha(png)))
            with (audit/'tables/figure_manifest.tsv').open('w') as f:
                w=csv.DictWriter(f,fieldnames=list(statuses[0]),delimiter='\t');w.writeheader();w.writerows(statuses)
            receipt.update(status='COMPLETED',finished_at=datetime.now(timezone.utc).isoformat(),figure_count=len(statuses),rerendered=sum(x['status']=='RERENDERED' for x in statuses),preserved_not_rerendered=sum(x['status']!='RERENDERED' for x in statuses))
            (audit/'plots_only_run.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
            # Full gallery: all real figures, no selected delivery copy.
            cards=[]
            for x in statuses:
                pdf=html.escape(x['figure_id'],quote=True);png=pdf[:-4]+'.png';title=html.escape(Path(x['figure_id']).stem)
                cards.append(f'<article><h3>{title}</h3><p>{x["status"]}</p><a href="{pdf}">PDF</a> · <a href="{png}">PNG</a><img loading="lazy" src="{png}" alt="{title}"></article>')
            (stage/'index.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>完整结果图册</title><style>body{font-family:sans-serif;max-width:1100px;margin:auto}article{border-bottom:1px solid #ddd;padding:20px}img{display:block;max-width:100%;height:auto}</style><h1>完整结果图册</h1><p>本轮统一语义配色与重绘；原统计结果和来源限制保持原样。</p>'+''.join(cards)+'</html>')
            for obsolete in stage.rglob('README.md'): obsolete.unlink()
            for folder in sorted((x for x in stage.rglob('*') if x.is_dir()),key=lambda x:len(x.parts),reverse=True):
                if not any(folder.iterdir()):folder.rmdir()
            run.publish(stage, result)
            print(json.dumps({k:v for k,v in receipt.items() if k!='inputs'},ensure_ascii=False))
        except BaseException:
            print('FAILED；正式结果保持原样；必要诊断由生命周期收尾保存。',file=sys.stderr)
            raise

if __name__=='__main__':main()

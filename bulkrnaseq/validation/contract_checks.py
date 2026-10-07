#!/usr/bin/env python3
"""在新临时副本重放已执行的 22 项输入和入口检查；不写模板目录。"""
from pathlib import Path
import shutil, tempfile, subprocess, csv, json
src=Path(__file__).resolve().parents[1]
def copy_template(destination):
 destination.mkdir()
 for name in ('rmd','scripts','config'):
  shutil.copytree(src/name,destination/name,ignore=shutil.ignore_patterns('__pycache__'))
 for name in ('rnaseq_first.sh','rnaseq_second.sh','rnaseq_third.sh','run_rmd.sh'):
  shutil.copyfile(src/name,destination/name)
 (destination/'rawdata').mkdir()
base=Path(tempfile.mkdtemp(prefix='bulk contract cases ')); d=base/'copy';copy_template(d)
config=d/'config/project.env';config.write_text(config.read_text().replace('RNASEQ_INPUT_MODE=fastq','RNASEQ_INPUT_MODE=counts'))
(d/'rawdata/counts.tsv').write_text('gene_id\tC1\tC2\tT1\tT2\ngene1\t1\t2\t3\t4\n')
sample=d/'config/samples.tsv'; original=sample.read_text(); checks=[]
def run(script,*args): return subprocess.run(['bash',str(d/script),*args],capture_output=True,text=True)
def record(name,condition,detail=''):
 checks.append({'case':name,'pass':bool(condition),'detail':detail}); assert condition,(name,detail)
before={str(p.relative_to(d)) for p in d.rglob('*')}
for script in ['rnaseq_first.sh','rnaseq_second.sh','rnaseq_third.sh','run_rmd.sh']:
 for flag in ['--help','--list','--preflight']:
  r=run(script,flag);record(script+flag,r.returncode==0,r.stderr[-200:])
after={str(p.relative_to(d)) for p in d.rglob('*')};record('readonly_help_list_preflight',before==after)
rows=list(csv.DictReader(original.splitlines(),delimiter='\t'))
for label, mutate in [('duplicate_sample',lambda rs:rs[1].update(sample=rs[0]['sample'])),('duplicate_biological_unit',lambda rs:rs[1].update(biological_id=rs[0]['biological_id'])),('unknown_biological_unit',lambda rs:rs[0].update(biological_id='')),('third_group',lambda rs:rs[0].update(group='other'))]:
 rs=[dict(x) for x in rows];mutate(rs)
 with sample.open('w') as h:
  w=csv.DictWriter(h,fieldnames=rows[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
 r=run('run_rmd.sh','--preflight');record(label,r.returncode!=0,r.stderr.strip())
sample.write_text(original)
# Row permutation does not change prepared group order or direction.
command='ROOT_DIR="$1"; source "$ROOT_DIR/scripts/config.sh"; "$PYTHON_BIN" "$ROOT_DIR/scripts/input_contract.py" prepare'
r=subprocess.run(['bash','-c',command,'bash',str(d)],capture_output=True,text=True);record('prepare_original',r.returncode==0,r.stderr)
selected=(d/'4.DE_analysis/samples.txt').read_text()
sample.write_text(original.splitlines()[0]+'\n'+'\n'.join(reversed(original.splitlines()[1:]))+'\n')
r=subprocess.run(['bash','-c',command,'bash',str(d)],capture_output=True,text=True);record('row_order_invariant',r.returncode==0 and selected==(d/'4.DE_analysis/samples.txt').read_text(),r.stderr)
(d/'rawdata/counts.tsv').write_text('changed')
r=subprocess.run(['bash','-c',command,'bash',str(d)],capture_output=True,text=True);record('stale_counts_rejected',r.returncode!=0 and 'STALE' in r.stderr,r.stderr)
# FASTQ basename differs from sample; dry run must use canonical sample clean names.
d2=base/'fastq copy';copy_template(d2)
(d2/'reference').mkdir();(d2/'reference/genes.gtf').write_text('fixture');(d2/'reference/genome.1.ht2').write_text('fixture')
for row in rows:
 for key in ('raw_fq1','raw_fq2'):(d2/row[key]).write_text('fixture')
before=set(d2.rglob('*'));r=subprocess.run(['bash',str(d2/'rnaseq_first.sh'),'--dry-run','--steps','1'],capture_output=True,text=True)
record('fastq_name_contract',r.returncode==0 and '/fastqc/C1_1.fastq.gz' in r.stdout and 'runA_R1.fastq.gz' in r.stdout,r.stderr)
record('fastq_dry_run_readonly',set(d2.rglob('*'))==before)
(base/'contract_results.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2)+'\n');print(json.dumps(checks,ensure_ascii=False,indent=2));print('Evidence directory:', base)

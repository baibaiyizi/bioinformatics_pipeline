#!/usr/bin/env python3
"""Create a standalone copy with synthetic SPORTS tables; reads no real project data."""
import argparse
import csv
from pathlib import Path
import random
import shutil
parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);parser.add_argument('--per-group',type=int,default=2);args=parser.parse_args()
source=Path(__file__).resolve().parents[1];target=Path(args.output).resolve()
if target.exists():raise FileExistsError(target)
if args.per_group<1:raise ValueError('per-group must be positive')
shutil.copytree(source,target,ignore=shutil.ignore_patterns('__pycache__','result','sportsresult','.cache'))
rng=random.Random(42)
rows=[]
classes=[('miRBase-miRNA_Match_Genome','mmu-miR-example','ACGTACGTACGTACGTACGTAC','miRBase-miRNA'),('GtRNAdb-mature-tRNA_Match_Genome','mature-tRNA-Gly-GCC','GGGGACGTACGTACGTACGTAC','mature-tRNA-Gly-GCC'),('rRNAdb-rRNA_Match_Genome','28S-rRNA','TTTTACGTACGTACGTACGTAC','28S-rRNA'),('piRNAdb-piRNA_Match_Genome','piR-example','CCCCACGTACGTACGTACGTAC','piR-example'),('Unannotated_Match_Genome','-','AAAACGTACGTACGTACGTACG','NO_Annotation')]
for group in ['control','case']:
 for index in range(args.per_group):
  sample=f'{group}_library{index+1}'
  rows.append([sample,group,f'{group}_animal{index+1}','individual','batch1',f'inputs/{sample}.fastq.gz',group,'','Mus musculus','synthetic_only',str(index+1)])
  folder=target/'sportsresult/04_sports_annotation/samples'/sample;folder.mkdir(parents=True)
  values=[rng.randint(100,300) for _ in classes]
  total=sum(values)
  with (folder/(sample+'_summary.txt')).open('w') as f:
   writer=csv.writer(f,delimiter='\t');writer.writerow(['Class','Sub_Class','Reads']);writer.writerows([['Clean_Reads','-',total],['Match_Genome','-',total]])
   for (cls,feature,seq,annotation),reads in zip(classes,values):
    writer.writerow([cls,'-',reads])
    if feature!='-':writer.writerow([cls,feature,reads])
  with (folder/(sample+'_length_distribution.txt')).open('w') as f:
   writer=csv.writer(f,delimiter='\t');writer.writerow(['Class','Length','Reads']);writer.writerows([['Clean_Reads',22,total],['Match_Genome',22,total]])
   for (cls,feature,seq,annotation),reads in zip(classes,values):writer.writerow([cls,22,reads])
  with (folder/(sample+'_output.txt')).open('w') as f:
   writer=csv.writer(f,delimiter='\t');writer.writerow(['ID','Sequence','Length','Reads','Match_Genome','Annotation'])
   for idx,((cls,feature,seq,annotation),reads) in enumerate(zip(classes,values)):writer.writerow([f'seq{idx+1}',seq,22,reads,'Yes',annotation])
with (target/'config/samples.tsv').open('w') as f:
 writer=csv.writer(f,delimiter='\t');writer.writerow(['sample','group','biological_sample_id','unit_type','batch','raw_fq','condition','age_group','species','stage','replicate']);writer.writerows(rows)
p=target/'config/project.env';s=p.read_text().replace("PROJECT_STAGE='configured_tissue'","PROJECT_STAGE='synthetic_only'");p.write_text(s)
print(target)

#!/usr/bin/env python3
"""独立合成 SPORTS/endo 表；不读取真实项目、不调用上游测序工具。"""
import argparse
import csv
from pathlib import Path
import shutil


def table(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(fields)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--species', required=True, choices=['Mus', 'Homo'])
    parser.add_argument('--endo', required=True, choices=['0', '1'])
    args = parser.parse_args()
    target = args.output.resolve()
    if target.exists():
        raise FileExistsError(target)
    source = Path(__file__).resolve().parents[1]
    target.mkdir(parents=True)
    shutil.copy2(source / 'run.sh', target / 'run.sh')
    shutil.copytree(source / 'shared', target / 'shared', ignore=shutil.ignore_patterns('__pycache__'))
    sample = 'synthetic_library_1'
    species = {'Mus': 'Mus musculus', 'Homo': 'Homo sapiens'}[args.species]
    prefix = 'mmu' if args.species == 'Mus' else 'hsa'
    candidate = 'ACGTACGTACGTACGTACGTAC'
    conflict = 'CCGTACGTACGTACGTACGTAC'
    unmatched = 'TCGTACGTACGTACGTACGTAC'
    short = 'GGGTACGTACGTACGTACGT'
    records = [(candidate, 30, 'Yes', 'NO_Annotation'), (conflict, 20, 'Yes', prefix+'-miR-fixture'), (unmatched, 10, 'No', 'NO_Annotation'), (short, 5, 'Yes', 'NO_Annotation')]
    (target / 'config').mkdir()
    (target / 'config/project.sh').write_text(f'''PANDORA_PROJECT_NAME=synthetic_fixture
PANDORA_DESIGN=single_library
PROJECT_SPECIES='{species}'
PROJECT_STAGE=synthetic
PROJECT_SAMPLE_TYPE=synthetic
PROJECT_MATRIX=synthetic
PROJECT_CONDITION=fixture
SPORTS_LAYOUT={'nested' if args.species == 'Mus' else 'flat'}
CONTRAST_BASELINE=''
CONTRAST_CASE=''
RUN_ENDOSIRNA={args.endo}
DB_ROOT=resources/smallrna
REFERENCE_ROOT=resources/genome
SPORTS_PL=resources/sports/source/sports.pl
ENDOSIRNA_REFS='fixture=resources/endo/fixture'
PRETRIM_FASTQ={'1' if args.species == 'Mus' else '0'}
TRIM_ADAPTER={'0' if args.species == 'Mus' else '1'}
CUTADAPT_NEXTSEQ_TRIM={'20' if args.species == 'Mus' else '0'}
''')
    fields = ['sample','group','condition','age_group','species','stage','sample_type','matrix','method','biological_sample','replicate','batch','analysis_role','pair_id','spikein_sequence','spikein_input_fmol','raw_fastq_gz']
    table(target / 'config/samples.tsv', fields, [[sample,'fixture','fixture','',species,'synthetic','synthetic','synthetic','1','synthetic_1','1','fixture','descriptive','','','0','01rawdata/not_generated.fq.gz']])
    sports = target / 'sportsresult'
    if args.species == 'Mus':
        sports /= '04_sports_annotation/samples'
    sports /= sample
    mirna = 'miRBase-miRNA_Match_Genome'
    annotated = [('Clean_Reads','-',65),('Match_Genome','-',55),('Unannotated_Match_Genome','-',35),('Unannotated_Unmatch_Genome','-',10),(mirna,'-',20),(mirna,prefix+'-miR-fixture',20)]
    table(sports / (sample+'_summary.txt'), ['Class','Sub_Class','Reads'], annotated)
    lengths = [('Clean_Reads',22,60),('Clean_Reads',20,5),('Match_Genome',22,50),('Match_Genome',20,5),('Unannotated_Match_Genome',22,30),('Unannotated_Match_Genome',20,5),('Unannotated_Unmatch_Genome',22,10),(mirna,22,20)]
    table(sports / (sample+'_length_distribution.txt'), ['Class','Length','Reads'], lengths)
    table(sports / (sample+'_output.txt'), ['ID','Sequence','Length','Reads','Match_Genome','Annotation'], [[f'seq_{i}',seq,len(seq),count,host,annotation] for i,(seq,count,host,annotation) in enumerate(records,1)])
    ref = target / 'resources/endo/fixture.fasta'
    ref.parent.mkdir(parents=True)
    ref.write_text(''.join(f'>synthetic_{i}\n{seq}\n' for i,(seq,_,_,_) in enumerate(records,1)))
    evidence = target / 'sportsresult/05_endosirna/samples' / sample / (sample+'.fixture')
    table(Path(str(evidence)+'.stats.tsv'), ['sample','ref_tag','endosirna_mismatch'], [[sample,'fixture',0]])
    table(Path(str(evidence)+'.sequence_counts.tsv'), ['sample','ref_tag','sequence','reads'], [[sample,'fixture',seq,count] for seq,count,_,_ in records])
    print(target)


if __name__ == '__main__':
    main()

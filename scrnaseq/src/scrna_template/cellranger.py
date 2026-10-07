from __future__ import annotations
import shutil
from .inputs import fastq_pairs
from .runtime import run_command, atomic_tsv
import pandas as pd

def cellranger_outs(context, library_id):
    return context.path('01_cellranger', library_id, 'outs')

def validate_cellranger(config):
    command = str(config.data['tools']['cellranger'])
    if not shutil.which(command):
        raise FileNotFoundError(f'Cell Ranger 不存在: {command}')
    for name in ['reference.json', 'fasta/genome.fa', 'genes/genes.gtf.gz', 'star']:
        if not (config.reference / name).exists():
            raise FileNotFoundError(f'Cell Ranger reference 缺失: {name}')

def run_cellranger(context):
    config = context.config
    validate_cellranger(config)
    outputs = [cellranger_outs(context, row.library_id) / 'filtered_feature_bc_matrix.h5' for row in config.samples().itertuples()]
    with context.stage('cellranger', outputs):
        stage = context.path('00_audit', 'fastqs')
        stage.mkdir()
        mapping = []
        for library, chunk, left, right, index in fastq_pairs(config):
            for read, source in [('R1', left), ('R2', right), ('I1', index)]:
                if source is None:
                    continue
                link = stage / f'{library}_S1_L001_{read}_{chunk:03d}.fastq.gz'
                link.symlink_to(source)
                mapping.append({'library_id': library, 'source': str(source), 'staged_name': link.name})
        atomic_tsv(context.path('00_audit', 'fastq_mapping.tsv'), pd.DataFrame(mapping))
        for row in config.samples().itertuples():
            command = [config.data['tools']['cellranger'], 'count', f'--id={row.library_id}', f'--transcriptome={config.reference}', f'--fastqs={stage}', f'--sample={row.library_id}', f'--chemistry={row.chemistry}', f"--localcores={int(config.data['resources']['cellranger_cores'])}", f"--localmem={int(config.data['resources']['cellranger_memory_gb'])}", '--create-bam=true', '--disable-ui']
            run_command(command, cwd=context.path('01_cellranger'), log_file=context.path('logs', f'cellranger_{row.library_id}.log'))

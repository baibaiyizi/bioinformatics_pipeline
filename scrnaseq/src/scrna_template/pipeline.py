from __future__ import annotations
from .inputs import validate_inputs, input_paths
from .runtime import RunContext, atomic_tsv, atomic_json

def run_pipeline(config, run_id):
    # Expensive imports belong to execution, not --help / --list-steps.
    from .cellranger import run_cellranger, validate_cellranger
    from .qc import run_qc
    from .analysis import run_atlas
    from .annotation import run_annotation
    from .descriptive import run_descriptive
    from .inference import run_inference
    from .reporting import run_report
    from .audit import audit_run
    audit = validate_inputs(config, full=True)
    paths = input_paths(config)
    if config.data['inputs']['mode'] == 'fastq':
        validate_cellranger(config)
        paths.extend([config.reference/'reference.json', config.reference/'fasta/genome.fa', config.reference/'genes/genes.gtf.gz'])
    context = RunContext.create(config, run_id, paths)
    atomic_tsv(context.path('00_audit', 'input_audit.tsv'), audit)
    try:
        if config.data['inputs']['mode'] == 'fastq':
            run_cellranger(context)
        else:
            context.record_status('cellranger', 'not_run', 'explicit raw h5ad input mode')
        run_qc(context)
        atlas = run_atlas(context)
        annotated = run_annotation(context, atlas)
        counts, metadata = run_descriptive(context, annotated)
        run_inference(context, counts, metadata)
        run_report(context)
        if not audit_run(context):
            raise RuntimeError('运行审计失败；查看 inference_status.tsv 和 audit_results.tsv')
    except Exception as exc:
        atomic_json(context.path('00_audit', 'run_result.json'), {'status': 'FAIL', 'error': f'{type(exc).__name__}: {exc}'})
        run_report(context)
        raise
    return context

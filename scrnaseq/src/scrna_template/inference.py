from __future__ import annotations
import pandas as pd
from .eligibility import assess_eligibility, design_matrix
from .qc import integer_counts
from .runtime import atomic_tsv
import numpy as np

def run_replicated_pseudobulk_deseq2(counts, metadata, case, control, alpha=0.05, min_log2fc=0.58, n_cpus=1):
    if not counts.index.equals(metadata.index):
        raise ValueError('counts/metadata 样本索引及顺序必须完全一致')
    if not counts.index.is_unique or not counts.columns.is_unique:
        raise ValueError('计数矩阵样本/基因 ID 必须唯一')
    libraries = metadata.reset_index().rename(columns={metadata.index.name or 'index': 'biological_sample_id'})
    libraries['library_id'] = libraries.biological_sample_id
    check = assess_eligibility(libraries, 3, case, control)
    if check['eligibility'] != 'GO':
        raise RuntimeError(check['reason'])
    integer_counts(counts.to_numpy())
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    design = '~condition' if metadata.sample_batch.nunique() == 1 else '~sample_batch + condition'
    dds = DeseqDataSet(counts=counts.astype(np.int64), metadata=metadata, design=design, refit_cooks=True, n_cpus=n_cpus, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=['condition', case, control], alpha=alpha, n_cpus=n_cpus, quiet=True)
    stats.summary()
    results = stats.results_df.reset_index(names='gene')
    results['passes_raw_p_exploratory'] = results.pvalue.lt(alpha) & results.log2FoldChange.abs().ge(min_log2fc)
    results['passes_BH_FDR'] = results.padj.lt(alpha) & results.log2FoldChange.abs().ge(min_log2fc)
    results['case'] = case
    results['control'] = control
    results['scope'] = 'exploratory raw P screening; BH adjusted P retained'
    results['statistical_unit'] = 'independent biological_sample_id (pool counts once)'
    return results

def run_inference(context, counts_path, metadata_path):
    status_path = context.path('06_optional', 'inference_status.tsv')
    stats = context.config.data['statistics']
    libraries = context.config.samples()
    global_gate = assess_eligibility(libraries, int(stats['min_independent_samples_per_group']), stats['case'], stats['control'])
    atomic_tsv(context.path('00_audit', 'analysis_eligibility.tsv'), pd.DataFrame([global_gate]))
    counts = pd.read_csv(counts_path, sep='\t', dtype={'cell_type': str, 'aggregation_key': str})
    metadata = pd.read_csv(metadata_path, sep='\t', dtype={'cell_type': str, 'aggregation_key': str, 'biological_sample_id': str, 'condition': str, 'sample_batch': str}, keep_default_na=False)
    statuses = []
    with context.stage('inference', [status_path]):
        for idx, (cell_type, block) in enumerate(counts.groupby('cell_type', sort=True)):
            state = {'cell_type': cell_type, 'eligibility': global_gate['eligibility'], 'execution': 'not_run', 'reason': global_gate['reason'], 'result_file': ''}
            statuses.append(state)
            if not stats['inference_enabled']:
                state.update(execution='disabled', reason='config statistics.inference_enabled=false')
            elif global_gate['eligibility'] != 'GO':
                state['execution'] = 'skipped_ineligible'
            elif cell_type == 'unknown':
                state.update(eligibility='NO-GO', execution='skipped_ineligible', reason='unknown 注释不构成明确细胞类型')
            else:
                meta = metadata.loc[(metadata.cell_type == cell_type) & (metadata.cells >= int(stats['min_cells_per_sample_cell_type']))].copy()
                present = set(meta.biological_sample_id)
                gate = assess_eligibility(libraries.loc[libraries.biological_sample_id.isin(present)], int(stats['min_independent_samples_per_group']), stats['case'], stats['control'])
                state.update(eligibility=gate['eligibility'], reason=gate['reason'])
                if gate['eligibility'] != 'GO':
                    state['execution'] = 'skipped_ineligible'
                else:
                    try:
                        state['execution'] = 'running'
                        atomic_tsv(status_path, pd.DataFrame(statuses))
                        matrix = block.set_index('aggregation_key').drop(columns='cell_type').loc[meta.aggregation_key]
                        matrix.index = pd.Index(meta.biological_sample_id, name='biological_sample_id')
                        meta = meta.set_index('biological_sample_id').loc[matrix.index]
                        matrix = matrix.loc[:, matrix.sum(axis=0) >= int(stats['min_gene_total_count'])]
                        if matrix.shape[1] < 2:
                            raise ValueError('可检验基因少于 2')
                        result = run_replicated_pseudobulk_deseq2(matrix, meta, stats['case'], stats['control'], float(stats['alpha']), float(stats['min_log2fc']), int(context.config.data['resources']['inference_cpus']))
                        result['cell_type'] = cell_type
                        target = context.path('06_optional', f'pseudobulk_DE_{idx:03d}.tsv')
                        atomic_tsv(target, result)
                        state.update(execution='completed', result_file=str(target.relative_to(context.run_dir)))
                    except Exception as exc:
                        state.update(execution='failed', reason=f'{type(exc).__name__}: {exc}')
            atomic_tsv(status_path, pd.DataFrame(statuses))
    return pd.DataFrame(statuses)

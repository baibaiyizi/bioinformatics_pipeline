from __future__ import annotations
import anndata as ad
import numpy as np
import pandas as pd
from .config import UNKNOWN
from .qc import integer_counts
from .runtime import atomic_tsv
from .plotting import stacked_composition

def aggregate_pseudobulk(data, libraries):
    """Aggregate raw integers across cells AND technical libraries by actual sample."""
    counts = integer_counts(data.layers['counts_raw'])
    mapping = libraries.set_index('library_id')
    if not set(data.obs.library_id.astype(str)) <= set(mapping.index):
        raise ValueError('矩阵 library_id 不属于样本表')
    keys, metadata = [], {}
    for lib in data.obs.library_id.astype(str):
        row = mapping.loc[lib]
        unit = str(row.biological_sample_id)
        group = libraries.loc[libraries.biological_sample_id == unit]
        conflict = any(group[col].nunique() > 1 for col in ['condition', 'sample_batch', 'unit_type'])
        resolved = unit.lower() not in UNKNOWN and not conflict
        key = 'sample:' + unit if resolved else 'unresolved-library:' + lib
        keys.append(key)
        metadata[key] = {'biological_sample_id': unit if resolved else 'unknown', 'aggregation_key': key, 'condition': row.condition, 'sample_batch': row.sample_batch, 'unit_type': row.unit_type, 'identity_resolved': resolved}
    obs = data.obs.copy()
    obs['aggregation_key'] = keys
    matrices, rows = {}, []
    for (cell_type, key), indices in obs.groupby(['cell_type', 'aggregation_key'], observed=True).indices.items():
        vector = np.asarray(counts[indices].sum(axis=0)).ravel().astype(np.int64)
        matrices[(str(cell_type), key)] = vector
        rows.append(dict(metadata[key], cell_type=str(cell_type), cells=len(indices), libraries=obs.iloc[indices].library_id.nunique()))
    meta = pd.DataFrame(rows)
    table = pd.DataFrame([{'cell_type': cell_type, 'aggregation_key': key, **dict(zip(data.var_names, vector, strict=True))} for (cell_type, key), vector in matrices.items()])
    return table, meta

def run_descriptive(context, annotation_path):
    outputs = [context.path('05_descriptive', name) for name in ('captured_composition.tsv', 'pseudobulk_counts.tsv', 'pseudobulk_metadata.tsv', 'condition_descriptive_gene_effects.tsv')]
    with context.stage('descriptive', outputs):
        data = ad.read_h5ad(annotation_path)
        grouped = data.obs.groupby(['library_id', 'condition', 'cell_type'], observed=True).size().rename('cells').reset_index()
        grouped['proportion'] = grouped.cells / grouped.groupby('library_id', observed=True).cells.transform('sum')
        grouped['scope'] = 'captured cell composition; no population proportion inference'
        atomic_tsv(outputs[0], grouped)
        atomic_tsv(context.path('07_source_data', 'captured_composition.tsv'), grouped)
        if context.config.data['report'].get('figures', True):
            stacked_composition(grouped, context.path('08_figures', 'captured_composition'), int(context.config.data['report']['dpi']))
        counts, metadata = aggregate_pseudobulk(data, context.config.samples())
        atomic_tsv(outputs[1], counts)
        atomic_tsv(outputs[2], metadata)
        stats = context.config.data['statistics']
        results = []
        for cell_type, block in counts.groupby('cell_type', sort=True):
            matrix = block.set_index('aggregation_key').drop(columns='cell_type')
            meta = metadata.loc[metadata.cell_type == cell_type].set_index('aggregation_key').loc[matrix.index]
            cpm = matrix.div(matrix.sum(axis=1).replace(0, np.nan), axis=0) * 1e6
            means = cpm.groupby(meta.condition).mean()
            case = means.loc[stats['case']] if stats['case'] in means.index else pd.Series(np.nan, index=matrix.columns)
            control = means.loc[stats['control']] if stats['control'] in means.index else pd.Series(np.nan, index=matrix.columns)
            result = pd.DataFrame({'gene': matrix.columns, 'cell_type': cell_type, 'case_cpm': case.values, 'control_cpm': control.values, 'log2_ratio_case_minus_control': np.log2((case.values+1)/(control.values+1)), 'case': stats['case'], 'control': stats['control'], 'scope': 'equal sample mean CPM; unresolved identity remains per-library; descriptive no P/FDR'})
            results.append(result)
        atomic_tsv(outputs[3], pd.concat(results, ignore_index=True))
    return outputs[1], outputs[2]

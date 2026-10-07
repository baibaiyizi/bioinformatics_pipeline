from __future__ import annotations
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
import numpy as np
import pandas as pd
from .cellranger import cellranger_outs
from .runtime import atomic_tsv

def integer_counts(matrix):
    matrix = sp.csr_matrix(matrix)
    if matrix.data.size and (not np.isfinite(matrix.data).all() or (matrix.data < 0).any() or not np.equal(matrix.data, np.rint(matrix.data)).all()):
        raise ValueError('原始 counts 必须是有限非负整数；拒绝归一化/对数表达')
    return matrix.astype(np.int64)

def run_qc(context):
    output = context.path('02_qc', 'qc_combined.h5ad')
    summaries, frames, cells = [], [], []
    config = context.config
    with context.stage('qc', [output, context.path('02_qc', 'qc_summary.tsv')]):
        for row in config.samples().itertuples(index=False):
            data = ad.read_h5ad(config.resolve(row.matrix)) if config.data['inputs']['mode'] == 'h5ad' else sc.read_10x_h5(cellranger_outs(context, row.library_id) / 'filtered_feature_bc_matrix.h5', gex_only=True)
            if not data.var_names.is_unique or not data.obs_names.is_unique:
                raise ValueError('输入矩阵 gene/cell ID 必须唯一；拒绝静默重命名')
            # This entry accepts raw matrices. An explicit counts_raw layer takes precedence over X.
            data.X = integer_counts(data.layers.get('counts_raw', data.X))
            data.layers['counts_raw'] = data.X.copy()
            data.layers['counts'] = data.X.copy()
            data.var['mt'] = data.var_names.str.startswith(config.data['qc']['mt_prefix'])
            sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], inplace=True, percent_top=None, log1p=False)
            data.obs['predicted_doublet'] = False
            doublet_status = 'disabled'
            if config.data['qc']['doublets']:
                import scrublet as scr
                runner = scr.Scrublet(data.X, expected_doublet_rate=float(config.data['qc']['expected_doublet_rate']), random_state=int(config.data['project']['seed']))
                scores, predicted = runner.scrub_doublets(n_prin_comps=min(30, data.n_obs-2, data.n_vars-2), verbose=False)
                if predicted is None:
                    raise RuntimeError(f'{row.library_id} Scrublet 阈值估计失败')
                data.obs['doublet_score'] = scores
                data.obs['predicted_doublet'] = predicted
                doublet_status = 'completed'
            q = config.data['qc']
            data.obs['qc_pass'] = (data.obs.total_counts >= q['min_counts']) & (data.obs.n_genes_by_counts >= q['min_genes']) & (data.obs.pct_counts_mt <= q['max_pct_mt']) & ~data.obs.predicted_doublet
            for key, value in row._asdict().items():
                data.obs[key] = str(value)
            data.obs_names = pd.Index([f'{row.library_id}|{barcode}' for barcode in data.obs_names], name='cell_id')
            before, after = data.n_obs, int(data.obs.qc_pass.sum())
            summaries.append({'library_id': row.library_id, 'condition': row.condition, 'biological_sample_id': row.biological_sample_id, 'called_cells': before, 'retained_cells': after, 'doublet_status': doublet_status, 'ambient_correction': 'not_run_raw_counts_retained'})
            cells.append(data.obs.reset_index()[['cell_id', 'library_id', 'total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'qc_pass']])
            frames.append(data[data.obs.qc_pass].copy())
        if any(set(frame.var_names) != set(frames[0].var_names) for frame in frames):
            raise ValueError('不同文库的基因集合不同；请使用同一参考和完整特征矩阵')
        combined = ad.concat(frames, join='inner', merge='same', index_unique=None)
        if combined.n_obs < 4:
            raise ValueError('QC 后细胞数不足')
        combined.write_h5ad(output, compression='gzip')
        atomic_tsv(context.path('02_qc', 'qc_summary.tsv'), pd.DataFrame(summaries))
        atomic_tsv(context.path('02_qc', 'cell_qc.tsv'), pd.concat(cells, ignore_index=True))
    return output

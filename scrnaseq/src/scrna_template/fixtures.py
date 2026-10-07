from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import yaml

def make_fixture(output, per_group=3, technical_libraries=1, cells_per_library=60):
    root = Path(output).resolve()
    if root.exists():
        raise FileExistsError(root)
    if min(per_group, technical_libraries, cells_per_library) < 1:
        raise ValueError('合成数据计数必须为正')
    (root/'config').mkdir(parents=True)
    (root/'inputs').mkdir()
    # Self-contained contract: never reads production or workspace config.
    config = {'project': {'root': '..', 'name': 'synthetic_only', 'run_root': 'runs', 'seed': 42}, 'samples': 'config/samples.tsv', 'markers': 'config/markers.yml', 'inputs': {'mode': 'h5ad', 'stability_seconds': 300}, 'resources': {'inference_cpus': 1}, 'qc': {'min_counts': 20, 'min_genes': 10, 'max_pct_mt': 100, 'mt_prefix': 'mt-', 'min_cells_per_gene': 1, 'doublets': False, 'expected_doublet_rate': 0.06}, 'preprocessing': {'target_sum': 10000, 'n_hvg': 100, 'n_pcs': 12, 'n_neighbors': 12, 'preferred_resolution': 0.5, 'leiden_resolutions': [0.5], 'exclude_hvg_genes': []}, 'annotation': {'min_fraction': 0.5, 'min_margin': 0.1}, 'statistics': {'design': 'independent_two_group', 'case': 'treated', 'control': 'control', 'inference_enabled': True, 'min_independent_samples_per_group': 3, 'min_cells_per_sample_cell_type': 10, 'min_gene_total_count': 10, 'alpha': 0.05, 'min_log2fc': 0.58, 'marker_min_log2fc': 0.58, 'marker_pvalue': 0.05}, 'report': {'figures': False, 'dpi': 80}}
    genes = [f'Gene{i:03d}' for i in range(150)]
    markers = {'cell_types': {'type_A': {'positive': genes[1:6], 'negative': genes[6:11]}, 'type_B': {'positive': genes[6:11], 'negative': genes[1:6]}}}
    rng = np.random.default_rng(42)
    rows = []
    for condition in ('control', 'treated'):
        for unit in range(per_group):
            subject = f'{condition}_sample{unit+1}'
            for tech in range(technical_libraries):
                lib = f'{subject}_lib{tech+1}'
                truth = np.arange(cells_per_library) % 2
                matrix = rng.negative_binomial(4, 0.65, (cells_per_library, len(genes))).astype(np.int64)
                matrix[truth==0, 1:6] += rng.poisson(15, ((truth==0).sum(),5))
                matrix[truth==1, 6:11] += rng.poisson(15, ((truth==1).sum(),5))
                matrix[:,0] += rng.poisson(60 if condition=='treated' else 1, cells_per_library)
                data = ad.AnnData(sp.csr_matrix(matrix), obs=pd.DataFrame({'synthetic_truth': np.where(truth==0,'type_A','type_B')}, index=[f'cell{i}' for i in range(cells_per_library)]), var=pd.DataFrame(index=genes))
                data.write_h5ad(root/'inputs'/f'{lib}.h5ad')
                rows.append({'library_id': lib, 'biological_sample_id': subject, 'unit_type': 'individual', 'condition': condition, 'sample_batch': 'batch1', 'chemistry': 'synthetic', 'matrix': f'inputs/{lib}.h5ad', 'r1': '', 'r2': '', 'i1': ''})
    pd.DataFrame(rows).to_csv(root/'config/samples.tsv', sep='\t', index=False)
    (root/'config/project.yml').write_text(yaml.safe_dump(config, sort_keys=False))
    (root/'config/markers.yml').write_text(yaml.safe_dump(markers, sort_keys=False))
    return root/'config/project.yml'

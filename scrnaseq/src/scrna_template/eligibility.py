from __future__ import annotations
import numpy as np
import pandas as pd
from .config import UNKNOWN

def known(series):
    return ~series.astype(str).str.strip().str.lower().isin(UNKNOWN)

def sample_metadata(libraries):
    """Do not turn unidentified libraries into presumed independent animals."""
    rows = []
    for unit, group in libraries.loc[known(libraries.biological_sample_id)].groupby('biological_sample_id', sort=True):
        rows.append({'biological_sample_id': unit, 'condition': group.condition.iloc[0], 'sample_batch': group.sample_batch.iloc[0], 'unit_type': group.unit_type.iloc[0], 'n_libraries': len(group)})
    return pd.DataFrame(rows, columns=['biological_sample_id', 'condition', 'sample_batch', 'unit_type', 'n_libraries']).set_index('biological_sample_id')

def design_matrix(metadata):
    return np.column_stack([np.ones(len(metadata)), pd.get_dummies(metadata.condition, drop_first=True, dtype=float), pd.get_dummies(metadata.sample_batch, drop_first=True, dtype=float)])

def assess_eligibility(libraries, minimum=3, case=None, control=None):
    reasons = []
    if not known(libraries.biological_sample_id).all():
        reasons.append('存在未知 biological_sample_id')
    if not libraries.unit_type.isin(['individual', 'pool']).all():
        reasons.append('存在未知生物学单位类型')
    if not known(libraries.sample_batch).all():
        reasons.append('sample_batch 缺失；不能评估批次混杂')
    known_libraries = libraries.loc[known(libraries.biological_sample_id)]
    for unit, group in known_libraries.groupby('biological_sample_id'):
        if any(group[c].nunique() != 1 for c in ('condition', 'sample_batch', 'unit_type')):
            reasons.append(f'{unit} 在 condition/sample_batch/unit_type 上冲突')
    metadata = sample_metadata(libraries)
    counts = metadata.groupby('condition').size().to_dict()
    conditions = [case, control] if case is not None else sorted(libraries.condition.unique())
    if len(conditions) != 2 or any(counts.get(c, 0) < minimum for c in conditions):
        reasons.append(f'独立单位数 {counts}，要求每组至少 {minimum}')
    if len(metadata):
        matrix = design_matrix(metadata)
        if np.linalg.matrix_rank(matrix) != matrix.shape[1] or matrix.shape[0] <= matrix.shape[1]:
            reasons.append('condition + sample_batch 设计不满秩或无残差自由度')
    else:
        reasons.append('无已知独立样本')
    return {'analysis': 'pseudobulk_condition_expression', 'eligibility': 'NO-GO' if reasons else 'GO', 'reason': '; '.join(reasons) or '独立两组、样本数和批次设计通过', 'n_units': str(counts), 'statistical_unit': 'biological_sample_id; each pool counts once'}

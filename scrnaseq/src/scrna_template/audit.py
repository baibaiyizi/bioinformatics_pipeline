from __future__ import annotations
import pandas as pd
from .runtime import atomic_tsv, atomic_json, sha256

def audit_run(context):
    status_path = context.path('06_optional', 'inference_status.tsv')
    statuses = pd.read_csv(status_path, sep='\t', keep_default_na=False)
    rows = []
    for row in statuses.itertuples():
        valid = row.execution in {'completed', 'disabled', 'skipped_ineligible'}
        if row.execution == 'completed':
            valid &= row.eligibility == 'GO' and bool(row.result_file) and (context.run_dir / row.result_file).is_file()
        rows.append({'check': f'inference:{row.cell_type}', 'status': 'PASS' if valid else 'FAIL', 'detail': f'eligibility={row.eligibility}; execution={row.execution}'})
    effects = pd.read_csv(context.path('05_descriptive', 'condition_descriptive_gene_effects.tsv'), sep='\t')
    rows.append({'check': 'descriptive_has_no_condition_pvalues', 'status': 'PASS' if not any('pval' in c.lower() or 'fdr' in c.lower() for c in effects) else 'FAIL', 'detail': 'sample/library-level descriptive CPM only'})
    artifacts = [p for p in context.run_dir.rglob('*') if p.is_file() and p.suffix in {'.tsv', '.h5ad', '.html', '.png', '.pdf'}]
    atomic_tsv(context.path('00_audit', 'output_sha256.tsv'), pd.DataFrame([{'path': str(p.relative_to(context.run_dir)), 'sha256': sha256(p)} for p in artifacts]))
    result = pd.DataFrame(rows)
    atomic_tsv(context.path('00_audit', 'audit_results.tsv'), result)
    atomic_json(context.path('00_audit', 'run_result.json'), {'status': 'PASS' if result.status.eq('PASS').all() else 'FAIL', 'scope': 'execution_integrity; not independent biological replication'})
    return bool(result.status.eq('PASS').all())

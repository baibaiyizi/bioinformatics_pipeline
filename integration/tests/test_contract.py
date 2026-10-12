from pathlib import Path
import importlib.util
import json
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
import yaml

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/integrate.py'
spec = importlib.util.spec_from_file_location('integration', SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class Contract(unittest.TestCase):
    def test_rna_reverse_direction_rejected(self):
        d = pd.DataFrame({'sampleA': ['WT'], 'sampleB': ['PFOS']})
        with self.assertRaises(ValueError):
            module.validate_direction(d, {'id': 'rna', 'kind': 'rna'}, {'reference': 'WT', 'comparison': 'PFOS'})

    def test_wgbs_unknown_direction_rejected(self):
        d = pd.DataFrame({'group1_label': ['WT'], 'group2_label': [None]})
        with self.assertRaises(ValueError):
            module.validate_direction(d, {'id': 'wgbs', 'kind': 'wgbs'}, {'reference': 'WT', 'comparison': 'PFOS'})

    def test_atac_sign_conflict_rejected(self):
        d = pd.DataFrame({'Direction': ['PFOS-open'], 'Fold': [-1.]})
        with self.assertRaises(ValueError):
            module.validate_direction(d, {'id': 'atac', 'kind': 'atac', 'effect_column': 'Fold'}, {'reference': 'WT', 'comparison': 'PFOS'})

    def test_atac_missing_direction_or_nonfinite_effect_rejected(self):
        for direction, effect in [(None, 1.), ('PFOS-open', float('nan')), ('PFOS-open', float('inf'))]:
            with self.subTest(direction=direction, effect=effect), self.assertRaises(ValueError):
                module.validate_direction(pd.DataFrame({'Direction': [direction], 'Fold': [effect]}),
                                          {'id': 'atac', 'kind': 'atac', 'effect_column': 'Fold'}, {'reference': 'WT', 'comparison': 'PFOS'})

    def fixture(self, folder):
        root = Path(folder)
        (root / 'config').mkdir()
        (root / 'scripts').mkdir()
        out = root / 'rnaseq/liver/result/full_cohort'
        out.mkdir(parents=True)
        de = pd.DataFrame({'gene_id': ['ENSMUSG00000000001'], 'sampleA': ['PFOS'], 'sampleB': ['WT']})
        de.to_csv(out / 'DE_results.csv', index=False)
        pd.DataFrame({'sample': ['WT1', 'PFOS1']}).to_csv(out / 'samples.tsv', sep='\t', index=False)
        pd.DataFrame({'gene_id': ['ENSMUSG00000000001'], 'WT1': [1.], 'PFOS1': [2.]}).to_csv(out / 'normalized_expression.tsv', sep='\t', index=False)
        sources=[]
        for role in ['counts','samples','annotation','producer']:
            p=out.parent.parent / (role + '.txt'); p.write_text(role)
            sources.append({'role':role,'path':p.name,'sha256':module.sha256(p)})
        (out/'producer_source.R').write_text('producer')
        p = dict(status='completed', samples=['WT1', 'PFOS1'], source=sources,
                 outputs={name: module.sha256(out / name) for name in ['DE_results.csv', 'samples.tsv', 'normalized_expression.tsv']})
        (out / 'provenance.json').write_text(json.dumps(p))
        pd.DataFrame({'status':['COMPLETED']}).to_csv(out/'execution.tsv',sep='\t',index=False)
        cfg = dict(title='测试', species='mouse', genome='GRCm39', reference='WT', comparison='PFOS', inputs=[dict(id='liver', kind='rna', context='liver', species='mouse', genome='GRCm39', path=str((out / 'DE_results.csv').relative_to(root)), source_status='completed')])
        (root / 'config/project.yml').write_text(yaml.safe_dump(cfg))
        return root, out

    def test_completed_rna_contract_passes(self):
        with tempfile.TemporaryDirectory(prefix='integration space ') as tmp:
            root, _ = self.fixture(tmp)
            self.assertEqual(module.Integration(root).check().iloc[0].status, 'PASS')

    def test_expression_signature_change_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = self.fixture(tmp)
            (out / 'normalized_expression.tsv').write_text('gene_id\tWT1\nENSMUSG1\t2\n')
            self.assertEqual(module.Integration(root).check().iloc[0].status, 'INVALID')

    def test_old_six_library_matrix_cannot_match_ten_library_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, out = self.fixture(tmp)
            d=json.loads((out/'provenance.json').read_text());d['samples'].append('WT2')
            (out/'provenance.json').write_text(json.dumps(d))
            self.assertEqual(module.Integration(root).check().iloc[0].status, 'INVALID')

    def test_upstream_counts_changed_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,out=self.fixture(tmp)
            (out.parent.parent/'counts.txt').write_text('changed')
            self.assertEqual(module.Integration(root).check().iloc[0].status,'INVALID')

    def test_running_producer_cannot_reuse_previous_completed_receipt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,out=self.fixture(tmp)
            pd.DataFrame({'status':['STARTED']}).to_csv(out/'execution.tsv',sep='\t',index=False)
            self.assertEqual(module.Integration(root).check().iloc[0].status,'INVALID')

    def test_wrong_species_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,_=self.fixture(tmp)
            cfg=yaml.safe_load((root/'config/project.yml').read_text());cfg['inputs'][0]['species']='human'
            (root/'config/project.yml').write_text(yaml.safe_dump(cfg))
            self.assertEqual(module.Integration(root).check().iloc[0].status,'INVALID')

    def test_failed_stages_do_not_publish_old_candidate_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,_=self.fixture(tmp)
            old=root/'result/multiomics/liver';old.mkdir(parents=True)
            pd.DataFrame({'has_atac_region':[True],'has_wgbs_region':[True]}).to_csv(old/'candidates.csv',index=False)
            workflow=module.Integration(root);workflow.check();workflow.synthesis([{'stage':'liver','status':'FAILED'}])
            report=(root/'README.md').read_text()
            self.assertNotIn('当前共有 **1**',report)

    def test_render_failure_invalidates_previous_completed_receipt_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, _ = self.fixture(tmp)
            workflow = module.Integration(root)
            workflow.result.mkdir()
            workflow.html_report.write_text('OLD COMPLETED REPORT')
            stale = workflow.result / 'synthesis'; stale.mkdir()
            (stale / 'candidate_celltype_context.csv').write_text('OLD CELLTYPES')
            (workflow.result / 'run.json').write_text(json.dumps({'status': 'COMPLETED', 'stages': []}))
            with patch.object(workflow, 'synthesis', side_effect=RuntimeError('render failed')):
                self.assertEqual(workflow.run('rna'), 1)
            current = json.loads((workflow.result / 'run.json').read_text())
            self.assertEqual(current['status'], 'FAILED')
            self.assertEqual(current['stages'][-1]['status'], 'FAILED')
            self.assertNotIn('OLD COMPLETED REPORT', workflow.html_report.read_text())
            self.assertTrue(list((root / 'history').glob('*/result/index.html')))
            self.assertFalse((stale / 'candidate_celltype_context.csv').exists())


if __name__ == '__main__':
    unittest.main()

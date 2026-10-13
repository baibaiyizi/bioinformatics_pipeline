from pathlib import Path
import importlib.util
import json
import tempfile
import unittest
import sys
from unittest.mock import patch

import pandas as pd
import yaml

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/integrate.py'
sys.path.insert(0, str(SOURCE.parent))
spec = importlib.util.spec_from_file_location('integration', SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class Contract(unittest.TestCase):
    def test_check_lists_missing_frozen_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,_=self.fixture(tmp)
            row=module.Integration(root).check(include_frozen_contract=True).set_index('id').loc['frozen_result_contract']
            self.assertEqual(row.status,'MISSING')
            self.assertIn('result/00_audit/run.json',row.reason)
            self.assertIn('term_id_map.tsv',row.reason)
            self.frozen_fixture(root)
            self.assertEqual(module.Integration(root).check(include_frozen_contract=True).set_index('id').loc['frozen_result_contract'].status,'PASS')

    def test_unmapped_terms_keep_explicit_frozen_labels(self):
        from project_integration import validate_program_terms
        terms=pd.DataFrame({'collection':['GO_BP'],'frozen_term':['GOBP_EXAMPLE'],'term_id':[None],'status':['UNMAPPED']})
        program=terms[['collection','frozen_term','term_id']].copy()
        audit=validate_program_terms(program,terms)
        self.assertEqual((audit['unmapped_rows'],audit['unmapped_frozen_terms']),(1,1))
        for altered in [terms.assign(status='MATCHED'),terms.assign(frozen_term='OTHER'),terms.assign(term_id='GO:0000001')]:
            with self.subTest(),self.assertRaises(ValueError):validate_program_terms(program,altered)

    def test_declared_project_integration_input_contract(self):
        cfg=yaml.safe_load((SOURCE.parent.parent/'config/project.yml').read_text())
        inputs={v['id']:v for v in cfg['inputs']}
        self.assertTrue({'metabolic_pathway_effects','islet_celltypes'}.issubset(inputs))
        self.assertEqual(inputs['islet_celltypes']['kind'],'celltypes')
        self.assertTrue({'pathway','log2FC','p','fdr'}.issubset(inputs['metabolic_pathway_effects']['required_columns']))

    def test_condition_colors_use_declared_groups(self):
        from figure_style import PROFILE
        self.assertEqual(PROFILE['palettes']['group'],{'WT':'#5F8FB3','PFOS':'#CE7C8E'})

    def test_current_inputs_cannot_reuse_different_frozen_statistics(self):
        from project_integration import run_project
        with tempfile.TemporaryDirectory() as tmp:
            root,_=self.fixture(tmp);self.frozen_fixture(root)
            receipt=json.loads((root/'result/00_audit/run.json').read_text());receipt['input_signatures']={'liver':'stale'}
            (root/'result/00_audit/run.json').write_text(json.dumps(receipt))
            self.assertEqual(run_project(module.Integration(root),integrate=False),1)
            self.assertTrue((root/'result/00_audit/run.json').exists())
            audit=json.loads(next(root.glob('.result-*/00_audit/run.json')).read_text())
            self.assertEqual(audit['status'],'FAILED')
            self.assertIn('输入签名',audit['reason'])

    def test_rna_reverse_direction_rejected(self):
        d = pd.DataFrame({'sampleA': ['WT'], 'sampleB': ['PFOS']})
        with self.assertRaises(ValueError):
            module.validate_direction(d, {'id': 'rna', 'kind': 'rna'}, {'reference': 'WT', 'comparison': 'PFOS'})

    def test_wgbs_unknown_direction_rejected(self):
        d = pd.DataFrame({'group1_label': ['WT'], 'group2_label': [None]})
        with self.assertRaises(ValueError):
            module.validate_direction(d, {'id': 'wgbs', 'kind': 'wgbs'}, {'reference': 'WT', 'comparison': 'PFOS'})

    def test_wgbs_effect_must_equal_actual_group_means(self):
        spec={'id':'wgbs','kind':'wgbs'};cfg={'reference':'WT','comparison':'PFOS'}
        d=pd.DataFrame({'group1_label':['WT'],'group2_label':['PFOS'],'meanMethy_group1':[.2],'meanMethy_group2':[.5],'diffMethy':[.3]})
        module.validate_direction(d,spec,cfg)
        for value in [-.3,float('nan'),float('inf')]:
            with self.subTest(value=value),self.assertRaises(ValueError):
                module.validate_direction(d.assign(diffMethy=value),spec,cfg)

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
        (root / 'gene_sets.rds').write_bytes(b'frozen gene set fixture')
        (root / 'figure_style.json').write_text('{}')
        cfg = dict(gene_sets='gene_sets.rds', figure_profile='figure_style.json', rna_context_order=[], title='测试', species='mouse', genome='GRCm39', reference='WT', comparison='PFOS', inputs=[dict(id='liver', kind='rna', context='liver', species='mouse', genome='GRCm39', path=str((out / 'DE_results.csv').relative_to(root)), source_status='completed')])
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

    def test_gene_set_change_changes_input_signature(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, _ = self.fixture(tmp)
            workflow = module.Integration(root)
            before = workflow.check().set_index('id').loc['gene_sets', 'sha256']
            (root / 'gene_sets.rds').write_bytes(b'changed static membership')
            after = workflow.check().set_index('id').loc['gene_sets', 'sha256']
            self.assertNotEqual(before, after)


    def frozen_fixture(self,root):
        from project_integration import digest
        out=root/'result/05_pairwise_support/01_rna_context/tables';out.mkdir(parents=True)
        source=out/'summary.csv';source.write_text('context,tested,nominal,fdr\nliver,10,2,1\n')
        audit=root/'result/00_audit/tables';audit.mkdir(parents=True)
        source_key='result/cross_context/rnaseq/summary.csv'
        relative=str(source.relative_to(root/'result'))
        pd.DataFrame([dict(source_key=source_key,path=relative,sha256=digest(source))]).to_csv(audit/'source_catalog.tsv',sep='\t',index=False)
        m=pd.DataFrame([dict(analysis='cross_context_rna',figure_id='rna_context_summary',png='result/05_pairwise_support/01_rna_context/plots/rna_context_summary.png',pdf='result/05_pairwise_support/01_rna_context/plots/rna_context_summary.pdf',source_table='result/'+relative,status='COMPLETED',reason='',title='test',biological_meaning='test',distribution='test',conclusion='test')])
        m.to_csv(audit/'figure_manifest.tsv',sep='\t',index=False)
        for name in ['membership.tsv','term_id_map.tsv']:(audit/name).write_text('frozen_term\nexample\n')
        receipt=dict(status='COMPLETED',derived_tables={})
        (root/'result/00_audit/run.json').write_text(json.dumps(receipt))
        stage=root/'stage';stage.mkdir()
        return source,stage

    def test_frozen_sources_move_without_copying_old_media(self):
        from project_integration import freeze_existing,digest
        with tempfile.TemporaryDirectory(prefix='integration source space ') as tmp:
            root=Path(tmp);source,stage=self.frozen_fixture(root)
            old=source.parent/'plots';old.mkdir();(old/'old.png').write_bytes(b'old image')
            cat,m,_=freeze_existing(root,stage)
            self.assertEqual(digest(stage/cat.iloc[0].path),digest(source))
            self.assertFalse(list(stage.rglob('*.png')))
            self.assertIn('/plots/',m.iloc[0].png)
            self.assertIn('/tables/',m.iloc[0].source_table)

    def test_changed_frozen_source_rejects_publication(self):
        from project_integration import freeze_existing
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source,stage=self.frozen_fixture(root)
            source.write_text('changed')
            with self.assertRaisesRegex(ValueError,'签名失配'):freeze_existing(root,stage)
            self.assertTrue((root/'result/00_audit/run.json').exists())

    def test_pending_producer_rejected_before_plotting(self):
        from project_integration import freeze_existing
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);_,stage=self.frozen_fixture(root)
            (root/'result/00_audit/run.json').write_text(json.dumps(dict(status='STARTED')))
            with self.assertRaisesRegex(ValueError,'COMPLETED'):freeze_existing(root,stage)

    def test_figure_identity_preserved_and_images_paired(self):
        from project_integration import freeze_existing
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);_,stage=self.frozen_fixture(root)
            original=pd.read_csv(root/'result/00_audit/tables/figure_manifest.tsv',sep='\t')
            _,manifest,_=freeze_existing(root,stage)
            pd.testing.assert_frame_equal(manifest,original)
            self.assertEqual(Path(manifest.iloc[0].png).stem,Path(manifest.iloc[0].pdf).stem)

    def test_no_statistical_fallback_in_plot_helper(self):
        helper=(SOURCE.parent/'joint_plot_helpers.R').read_text()
        self.assertNotIn('fgseaMultilevel',helper)
        self.assertNotIn('phyper(',helper)
        self.assertIn('纯绘图禁止运行富集检验',helper)

    def test_cross_context_source_is_signed_in_new_scripts_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,_=self.fixture(tmp)
            branch=root/'scripts/cross_wgbs/pair';branch.mkdir(parents=True)
            p=branch/'com.Rmd';p.write_text('original')
            w=module.Integration(root);old=w.code_signatures();p.write_text('changed')
            self.assertNotEqual(old,w.code_signatures())


if __name__ == '__main__':
    unittest.main()

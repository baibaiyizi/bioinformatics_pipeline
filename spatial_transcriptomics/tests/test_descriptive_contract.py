import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

spec = importlib.util.spec_from_file_location('spatial_pipeline', Path(__file__).parents[1] / 'scripts/spatial_pipeline.py')
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)

class DescriptiveContractTests(unittest.TestCase):
    def fixture(self, root):
        ctx = m.Context(root, copy.deepcopy(m.DEFAULT_CONFIG), [])
        rng = np.random.default_rng(9)
        counts = rng.poisson(4, (16, 12))
        obs = pd.DataFrame({'sample': ['one'] * 8 + ['two'] * 8, 'group': ['A'] * 8 + ['B'] * 8,
                            'score_type': rng.normal(size=16)}, index=[f's{i}' for i in range(16)])
        a = ad.AnnData(-np.ones(counts.shape), obs=obs, var=pd.DataFrame(index=[f'g{i}' for i in range(12)]))
        a.layers['counts'] = sparse.csr_matrix(counts)
        m.write_adata(a, ctx.integrated_dir() / 'all_samples.annotated.h5ad')
        for stage in ('09_deconvolution', '11_programs_descriptive', '14_report_index'):
            (ctx.step_dir(stage) / 'tables').mkdir(parents=True)
            (ctx.step_dir(stage) / 'plots').mkdir()
        return ctx, counts

    def test_sample_summary_uses_counts_and_no_inference(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx, counts = self.fixture(Path(tmp))
            m.step_11_programs_descriptive(ctx, '11_programs_descriptive')
            base = ctx.step_dir('11_programs_descriptive') / 'tables'
            out = pd.read_csv(base / 'sample_expression_descriptive.tsv', sep='\t')
            self.assertEqual(len(out), 24)
            self.assertEqual(out.query('sample == "one" and gene == "g0"').sum_raw_counts.iloc[0], counts[:8, 0].sum())
            self.assertTrue((out.mean_log1p_normalized > 0).all())
            self.assertFalse(any('pvalue' in x.lower() or 'logfc' in x.lower() for x in out.columns))
            status = json.loads((base / 'inference_status.json').read_text())
            self.assertFalse(status['pvalues_produced'])
            m.step_14_report_index(ctx, '14_report_index')
            self.assertTrue((ctx.result_dir / 'README.md').is_file())

    def test_marker_score_is_not_abundance(self):
        with tempfile.TemporaryDirectory() as tmp:
            ctx, _ = self.fixture(Path(tmp))
            m.step_09_deconvolution(ctx, '09_deconvolution')
            base = ctx.step_dir('09_deconvolution') / 'tables'
            self.assertTrue((base / 'marker_scores_by_spot.tsv').is_file())
            self.assertEqual(json.loads((base / 'method_status.json').read_text())['output_kind'], 'marker_score_not_abundance')

    def test_group_conflict_and_missing_counts_rejected(self):
        for missing in (False, True):
            with tempfile.TemporaryDirectory() as tmp:
                ctx, _ = self.fixture(Path(tmp))
                path = ctx.integrated_dir() / 'all_samples.annotated.h5ad'
                a = ad.read_h5ad(path)
                if missing:
                    del a.layers['counts']
                else:
                    a.obs.loc['s0', 'group'] = 'B'
                a.write_h5ad(path)
                with self.assertRaises(ValueError):
                    m.step_11_programs_descriptive(ctx, '11_programs_descriptive')

class InputBindingTests(unittest.TestCase):
    def raw_fixture(self, root):
        cfg=copy.deepcopy(m.DEFAULT_CONFIG)
        ctx=m.Context(root,cfg,[{'sample':'A','group':'control','organism':'mouse','input_type':'h5ad','input_path':'input.h5ad'}])
        a=ad.AnnData(np.full((8,12),0.25),obs=pd.DataFrame(index=[f's{i}' for i in range(8)]),var=pd.DataFrame(index=[f'g{i}' for i in range(12)]))
        a.layers['counts']=sparse.csr_matrix(np.arange(96).reshape(8,12)+1)
        a.obsm['spatial']=np.column_stack([np.arange(8),np.arange(8)])
        a.write_h5ad(root/'input.h5ad')
        return ctx

    def test_empty_paths_rejected_without_outputs(self):
        for mode in ('h5ad','matrix_coords'):
            with tempfile.TemporaryDirectory() as tmp:
                ctx=self.raw_fixture(Path(tmp));ctx.manifest[0].update(input_type=mode,input_path='',count_file='',spatial_file='')
                with self.assertRaises(SystemExit):m.bind_run(ctx,write=False)
                self.assertFalse(ctx.result_dir.exists())

    def test_counts_restored_before_qc_and_changed_input_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ctx=self.raw_fixture(root)
            m.bind_run(ctx,write=False);self.assertFalse(ctx.result_dir.exists())
            m.bind_run(ctx,write=True);m.execute_step(ctx,'02_build_objects')
            result=ad.read_h5ad(ctx.processed_dir()/'all_samples.raw.h5ad')
            np.testing.assert_array_equal(result.X.toarray(),result.layers['counts'].toarray())
            m.bind_run(ctx,write=False)
            a=ad.read_h5ad(root/'input.h5ad');a.layers['counts'][0,0]=200;a.write_h5ad(root/'input.h5ad')
            with self.assertRaisesRegex(ValueError,'STALE'):m.bind_run(ctx,write=False)

    def test_missing_config_and_unsigned_outputs_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ctx=self.raw_fixture(root)
            with self.assertRaises(FileNotFoundError):m.load_config(root/'missing.yml')
            ctx.processed_dir().mkdir();(ctx.processed_dir()/'old.h5ad').write_text('old')
            with self.assertRaisesRegex(ValueError,'UNVERIFIED'):m.bind_run(ctx,write=False)

    def test_custom_manifest_and_result_sequential_shell(self):
        import shutil, subprocess, os, yaml
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/"copy space"
            shutil.copytree(Path(__file__).parents[1],root,ignore=shutil.ignore_patterns('__pycache__','result'))
            ctx=self.raw_fixture(root)
            ctx.config['paths']['manifest']='custom.tsv'
            ctx.config['paths']['result_dir']='custom_result'
            pd.DataFrame(ctx.manifest).to_csv(root/'custom.tsv',sep='\t',index=False)
            (root/'config/spatial_config.yaml').write_text(yaml.safe_dump(ctx.config))
            env=os.environ.copy();env['SPATIAL_PYTHON_BIN']=sys.executable
            for args in [['--check'],['--steps','1,2'],['--steps','14']]:
                result=subprocess.run(['bash',str(root/'run_spatial.sh'),*args],env=env,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            states=json.loads((root/'custom_result/step_status.json').read_text())
            self.assertEqual(states['14_report_index']['status'],'completed')
            self.assertTrue((root/'custom_result/README.md').is_file())

    def test_failed_upstream_cannot_supply_old_h5ad(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ctx=self.raw_fixture(root)
            m.bind_run(ctx,write=True);m.execute_step(ctx,'02_build_objects')
            m.write_json(ctx.result_dir/'step_status.json',{'02_build_objects':{'status':'failed'}})
            with self.assertRaisesRegex(ValueError,'上游失败'):
                m.execute_step(ctx,'03_qc_filter')
            m.execute_step(ctx,'14_report_index')
            self.assertIn('failed',(ctx.result_dir/'README.md').read_text())

if __name__ == '__main__':
    unittest.main()

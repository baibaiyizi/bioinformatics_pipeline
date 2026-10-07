from pathlib import Path
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml

from scrna_template.config import WorkflowConfig
from scrna_template.eligibility import assess_eligibility
from scrna_template.fixtures import make_fixture
from scrna_template.descriptive import aggregate_pseudobulk
from scrna_template.inference import run_replicated_pseudobulk_deseq2
from scrna_template.inputs import full_pair_check, parse_md5_manifest, validate_inputs, InputGateError
from scrna_template.qc import integer_counts


def libraries(n=3, technical=1):
    return pd.DataFrame([{'library_id': f'{c}_{i}_{t}', 'biological_sample_id': f'{c}_{i}', 'condition': c, 'sample_batch': 'b1', 'unit_type': 'individual'} for c in ['case','control'] for i in range(n) for t in range(technical)])


def test_eligibility_counts_units_not_libraries_or_pool_animals():
    frame = libraries(1, 3)
    assert assess_eligibility(frame)['eligibility'] == 'NO-GO'
    frame = libraries(3,2)
    frame['unit_type'] = 'pool'
    frame['n_animals_in_pool'] = 100
    gate = assess_eligibility(frame)
    assert gate['eligibility'] == 'GO'
    assert "'case': 3" in gate['n_units']


@pytest.mark.parametrize('kind', ['unknown_unit', 'confounded_batch', 'missing_batch', 'unit_conflict', 'unit_type_conflict'])
def test_disqualifying_metadata_is_descriptive_only(kind):
    frame = libraries(3,2)
    if kind == 'unknown_unit':
        frame.loc[0,'biological_sample_id'] = 'unknown'
    elif kind == 'confounded_batch':
        frame['sample_batch'] = frame.condition
    elif kind == 'missing_batch':
        frame.loc[0,'sample_batch'] = ''
    elif kind == 'unit_conflict':
        frame.loc[0,'condition'] = 'control'
    else:
        frame.loc[0,'unit_type'] = 'pool'
    assert assess_eligibility(frame)['eligibility'] == 'NO-GO'


def test_technical_libraries_aggregate_integer_raw_counts():
    libs = libraries(1,2)
    data = ad.AnnData(sp.csr_matrix(np.full((4,3),100.0)), obs=pd.DataFrame({'library_id':libs.library_id.tolist(),'cell_type':['A']*4}, index=['c0','c1','c2','c3']), var=pd.DataFrame(index=['g0','g1','g2']))
    raw = np.arange(12).reshape(4,3)
    data.layers['counts_raw'] = sp.csr_matrix(raw)
    counts, meta = aggregate_pseudobulk(data, libs)
    assert len(counts) == 2
    observed = counts.set_index('aggregation_key').loc['sample:case_0', ['g0','g1','g2']].to_numpy(dtype=int)
    np.testing.assert_array_equal(observed, raw[:2].sum(axis=0))
    assert set(meta.libraries) == {2}


def test_unknown_identity_never_merged_across_libraries():
    libs = libraries(1,2)
    libs['biological_sample_id'] = 'unknown'
    data = ad.AnnData(sp.csr_matrix(np.ones((4,2),dtype=int)), obs=pd.DataFrame({'library_id':libs.library_id.tolist(),'cell_type':['A']*4},index=['a','b','c','d']), var=pd.DataFrame(index=['g1','g2']))
    data.layers['counts_raw'] = data.X.copy()
    counts, meta = aggregate_pseudobulk(data,libs)
    assert len(counts) == 4
    assert not meta.identity_resolved.any()


@pytest.mark.parametrize('bad', [np.array([[0.5]]), np.array([[-1]]), np.array([[np.nan]])])
def test_raw_integer_gate(bad):
    with pytest.raises(ValueError):
        integer_counts(bad)


def test_replicated_planted_effect_and_direction():
    libs = libraries()
    metadata = libs.set_index('biological_sample_id')[['condition','sample_batch','unit_type']]
    rng = np.random.default_rng(12)
    counts = pd.DataFrame(rng.negative_binomial(10,0.5,(6,50))+10,index=metadata.index,columns=[f'Gene{i}' for i in range(50)])
    counts.loc[metadata.condition=='case','Gene0'] += 1000
    result = run_replicated_pseudobulk_deseq2(counts,metadata,'case','control').set_index('gene')
    assert result.loc['Gene0','log2FoldChange'] > 0.58
    assert result.loc['Gene0','pvalue'] < 0.05
    assert 'padj' in result
    reverse = run_replicated_pseudobulk_deseq2(counts.iloc[::-1],metadata.iloc[::-1],'control','case').set_index('gene')
    assert reverse.loc['Gene0','log2FoldChange'] < -0.58
    assert np.isclose(result.loc['Gene0','log2FoldChange'], -reverse.loc['Gene0','log2FoldChange'])


def test_config_retains_empty_optional_fields_and_rejects_duplicate(tmp_path):
    path = make_fixture(tmp_path/'fixture',1)
    config = WorkflowConfig.load(path)
    assert config.samples().i1.eq('').all()
    samples = config.samples()
    samples.loc[1,'library_id'] = samples.loc[0,'library_id']
    samples.to_csv(config.samples_path,sep='\t',index=False)
    with pytest.raises(ValueError, match='library_id'):
        WorkflowConfig.load(path)


def test_paired_design_rejected(tmp_path):
    path = make_fixture(tmp_path/'fixture',1)
    config = yaml.safe_load(path.read_text())
    config['statistics']['design'] = 'paired'
    path.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError,match='配对'):
        WorkflowConfig.load(path)


def fastq(path, names):
    with gzip.open(path,'wt') as handle:
        for name in names:
            handle.write(f'@{name}\nACGT\n+\nIIII\n')


def test_full_fastq_pair_optional_i1(tmp_path):
    left,right,index = [tmp_path/name for name in ['raw_R1.fastq.gz','raw_R2.fastq.gz','raw_I1.fastq.gz']]
    fastq(left,['a/1','b/1']); fastq(right,['a/2','b/2']); fastq(index,['a','b'])
    assert full_pair_check(left,right) == 2
    assert full_pair_check(left,right,index) == 2
    fastq(right,['a/2','wrong/2'])
    with pytest.raises(InputGateError,match='配对'):
        full_pair_check(left,right)


def test_full_input_audit_preserves_names_and_stability(tmp_path):
    config_path = make_fixture(tmp_path/'fixture',1)
    config = yaml.safe_load(config_path.read_text())
    config['inputs'].update(mode='fastq',md5_manifest='inputs/MD5.txt',stability_seconds=30,require_no_open_writer=False)
    config_path.write_text(yaml.safe_dump(config))
    samples = pd.read_csv(config_path.parent/'samples.tsv',sep='\t',keep_default_na=False)
    root = config_path.parent.parent
    manifest = []
    for idx in range(2):
        for read in ['r1','r2']:
            rel = f'inputs/original_vendor_{idx}_{read}.fastq.gz'
            fastq(root/rel,['a','b'])
            samples.loc[idx,read]=rel
            manifest.append(f'{hashlib.md5((root/rel).read_bytes()).hexdigest()}  {Path(rel).name}')
    samples.to_csv(config_path.parent/'samples.tsv',sep='\t',index=False)
    (root/'inputs/MD5.txt').write_text('\n'.join(manifest))
    cfg = WorkflowConfig.load(config_path)
    with pytest.raises(InputGateError, match='稳定'):
        validate_inputs(cfg,full=True)
    for p in (root/'inputs').iterdir():
        os.utime(p,(time.time()-120,time.time()-120))
    result = validate_inputs(cfg,full=True)
    assert result.status.eq('PASS').all()
    assert not (root/'runs').exists()
    assert len(list((root/'inputs').glob('original_vendor*'))) == 4


def test_copied_template_e2e_n1_and_help_no_artifacts(tmp_path):
    source = Path(__file__).resolve().parents[1]
    target = tmp_path/'copied_template'
    shutil.copytree(source,target,ignore=shutil.ignore_patterns('__pycache__','.pytest_cache','runs'))
    env = dict(os.environ,SCRNA_PYTHON=sys.executable,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
    subprocess.run([str(target/'run_scrna.sh'),'--help'],check=True,cwd=tmp_path,env=env,capture_output=True)
    subprocess.run([str(target/'run_scrna.sh'),'--list-steps'],check=True,cwd=tmp_path,env=env,capture_output=True)
    assert not (target/'runs').exists()
    fixture = tmp_path/'standalone'
    subprocess.run([str(target/'run_scrna.sh'),'fixture','--output',str(fixture),'--per-group','1'],check=True,cwd=tmp_path,env=env,capture_output=True)
    config = fixture/'config/project.yml'
    subprocess.run([str(target/'run_scrna.sh'),'--config',str(config),'validate'],check=True,cwd=tmp_path,env=env,capture_output=True)
    assert not (fixture/'runs').exists()
    subprocess.run([str(target/'run_scrna.sh'),'--config',str(config),'run','--run-id','smoke'],check=True,cwd=tmp_path,env=env,capture_output=True)
    run = fixture/'runs/smoke'
    statuses = pd.read_csv(run/'06_optional/inference_status.tsv',sep='\t')
    assert statuses.execution.eq('skipped_ineligible').all()
    assert json.loads((run/'00_audit/run_result.json').read_text())['status']=='PASS'
    assert '仅表报告' in (run/'09_reports/scrna_report_zh.html').read_text()
    assert not list((run/'06_optional').glob('pseudobulk_DE*'))


def test_inference_failure_records_failed_execution(tmp_path,monkeypatch):
    from scrna_template.runtime import RunContext
    from scrna_template.inference import run_inference
    path=make_fixture(tmp_path/'fixture',3)
    config=WorkflowConfig.load(path)
    run=tmp_path/'run'
    run.mkdir()
    context=RunContext(config,'test',run,'signature')
    samples=config.samples()
    keys=['sample:'+s for s in samples.biological_sample_id]
    counts=pd.DataFrame({'cell_type':['A']*6,'aggregation_key':keys,'g1':range(10,16),'g2':range(20,26)})
    meta=samples[['biological_sample_id','condition','sample_batch','unit_type']].copy()
    meta['cell_type']='A'; meta['aggregation_key']=keys; meta['cells']=50
    counts.to_csv(tmp_path/'counts.tsv',sep='\t',index=False);meta.to_csv(tmp_path/'meta.tsv',sep='\t',index=False)
    def fail(*args,**kwargs):
        raise RuntimeError('injected fit failure')
    monkeypatch.setattr('scrna_template.inference.run_replicated_pseudobulk_deseq2',fail)
    result=run_inference(context,tmp_path/'counts.tsv',tmp_path/'meta.tsv')
    assert result.eligibility.eq('GO').all()
    assert result.execution.eq('failed').all()
    assert result.reason.str.contains('injected fit failure').all()

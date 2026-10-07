#!/usr/bin/env python3
"""在临时副本中重放峰表核心、方向和生物学单位验证；不运行真实XCMS。"""
from pathlib import Path
import tempfile,shutil,subprocess,json
import numpy as np,pandas as pd,yaml
src=Path(__file__).resolve().parents[1] / 'metabolomics';root=Path(tempfile.mkdtemp(prefix='model-lcms-validation-'));log=(root/'run.log').open('w');cases=[]
for label,detected,reverse in [('missing_prefill',False,False),('with_prefill',True,False),('reverse',False,True)]:
 dst=root/(label+' space');shutil.copytree(src,dst);(dst/'inputs').mkdir();rng=np.random.default_rng(29);samples=['C1','T1','C2','T2','C3','T3'];a=rng.lognormal(6,.25,(60,6));a[:5,[1,3,5]]*=5
 frame=pd.DataFrame(a,index=[f'f{i}' for i in range(60)],columns=samples).rename_axis('feature_id');frame.to_csv(dst/'inputs/filled.csv');pd.DataFrame({'feature_id':frame.index,'mz':np.arange(60)+100,'rt':np.arange(60)+30,'polarity':'POS'}).to_csv(dst/'inputs/features.csv',index=False)
 cfg=yaml.safe_load((dst/'config/project.yml').read_text())
 if detected:
  d=frame.copy();d.iloc[10,:]=np.nan;d.to_csv(dst/'inputs/detected.csv');cfg['input']['detected_table']='inputs/detected.csv'
 if reverse:cfg['case_group'],cfg['control_group']=cfg['control_group'],cfg['case_group']
 (dst/'config/project.yml').write_text(yaml.safe_dump(cfg,sort_keys=False))
 subprocess.run(['bash',str(dst/'run.sh'),'check'],stdout=log,stderr=log,check=True);assert not (dst/'result').exists()
 subprocess.run(['Rscript',str(dst/'pipeline.R'),'--module','all'],stdout=log,stderr=log,check=True)
 result=pd.read_csv(dst/'result/04_statistics'/f"{cfg['case_group']}_vs_{cfg['control_group']}_differential_features.csv").set_index('feature_id').sort_index()
 assert (dst/'result/report.md').exists();q=pd.read_csv(dst/'result/03_qc/feature_filter_decision.csv')
 if not detected:assert q.detection_rate.isna().all() and (q.detection_status=='not_available_no_prefill_matrix').all()
 else:assert not bool(q.set_index('feature_id').loc['f10','keep'])
 if label=='missing_prefill':base=result
 if reverse:
  np.testing.assert_allclose(base.log2FC,-result.log2FC,atol=1e-9);np.testing.assert_allclose(base.limma_p,result.limma_p,rtol=1e-8)
 cases.append(label)
# duplicate biological units rejected, no outputs created
neg=root/'bad_units';shutil.copytree(src,neg);(neg/'inputs').mkdir();shutil.copy2(root/'missing_prefill space/inputs/filled.csv',neg/'inputs/filled.csv');shutil.copy2(root/'missing_prefill space/inputs/features.csv',neg/'inputs/features.csv')
si=pd.read_csv(neg/'config/samples.tsv',sep='\t');si.loc[1,'biological_sample_id']=si.loc[0,'biological_sample_id'];si.to_csv(neg/'config/samples.tsv',sep='\t',index=False)
p=subprocess.run(['bash',str(neg/'run.sh'),'check'],stdout=log,stderr=log);assert p.returncode!=0 and not (neg/'result').exists();cases.append('duplicate_unit_rejected')
# Failed producers cannot supply an old successful matrix.
failed=root/'failed_producer';shutil.copytree(src,failed);(failed/'inputs').mkdir()
for name in ['filled.csv','features.csv']:
 shutil.copy2(root/'missing_prefill space/inputs'/name,failed/'inputs'/name)
subprocess.run(['Rscript',str(failed/'pipeline.R'),'--module','xcms'],stdout=log,stderr=log,check=True)
state=failed/'result/module_status.tsv';state.write_text(state.read_text()+'xcms\tfailed\tsynthetic_failure\n')
blocked=subprocess.run(['Rscript',str(failed/'pipeline.R'),'--module','preprocess'],capture_output=True,text=True)
assert blocked.returncode!=0 and '上游失败或未完成' in blocked.stderr
cases.append('failed_producer_rejected')
log.close();print(json.dumps({'root':str(root),'passed':cases},ensure_ascii=False));(root / 'validation_summary.json').write_text(json.dumps({'root':str(root),'passed':cases},ensure_ascii=False))

#!/usr/bin/env python3
"""用具备Scanpy的解释器运行，在临时副本完成空间counts至描述报告。"""
from pathlib import Path
import tempfile,shutil,subprocess,os,json,sys
import numpy as np,pandas as pd,anndata as ad,yaml
from scipy import sparse
src=Path(__file__).resolve().parents[1] / 'spatial_transcriptomics'; root=Path(tempfile.mkdtemp(prefix='spatial-template-e2e-'))/'copy with space';shutil.copytree(src,root,ignore=shutil.ignore_patterns('__pycache__', 'result'));(root/'inputs').mkdir();rng=np.random.default_rng(29)
rows=[]
for i in range(2):
 counts=rng.poisson(rng.lognormal(.9,.7,(1,150)),size=(40,150));counts[:20,:8]+=10
 a=ad.AnnData(np.full(counts.shape,0.25),obs=pd.DataFrame(index=[f's{j}' for j in range(40)]),var=pd.DataFrame(index=[f'g{j}' for j in range(150)]));a.layers['counts']=sparse.csr_matrix(counts);a.obsm['spatial']=np.column_stack([np.arange(40)%8,np.arange(40)//8]);a.write_h5ad(root/f'inputs/A{i}.h5ad');rows.append({'sample':f'A{i}','group':['control','treated'][i],'organism':'mouse','input_type':'h5ad','input_path':f'inputs/A{i}.h5ad','exclude':'0'})
pd.DataFrame(rows).to_csv(root/'01rawdata/spatial_info.tsv',sep='\t',index=False)
pd.DataFrame({'tissue':['generic']*4,'cell_type':['typeA','typeA','typeB','typeB'],'gene':['g0','g1','g10','g11']}).to_csv(root/'config/marker_genes.tsv',sep='\t',index=False)
cfg=yaml.safe_load((root/'config/spatial_config.yaml').read_text());cfg['qc'].update(min_counts=1,min_genes=1,min_cells_per_gene=1);cfg['normalization'].update(n_top_genes=80,n_pcs=8,regress_out=[]);cfg['programs']['n_components']=3;cfg['integration']['n_neighbors']=6;(root/'config/spatial_config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
env=os.environ.copy();env.update(SPATIAL_PYTHON_BIN=sys.executable,NUMBA_NUM_THREADS='2',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MPLBACKEND='Agg')
with (root.parent/'execution.log').open('w') as log:
 subprocess.run(['bash',str(root/'run_spatial.sh'),'--check'],env=env,stdout=log,stderr=log,check=True);assert not (root/'result').exists()
 subprocess.run(['bash',str(root/'run_spatial.sh'),'--steps','1,2,3,4,5,7,8,9,11,13,14'],env=env,stdout=log,stderr=log,check=True)
out=pd.read_csv(root/'result/11_programs_descriptive/tables/sample_expression_descriptive.tsv',sep='\t');assert len(out)==300 and not any('pvalue' in c.lower() or 'logfc' in c.lower() for c in out)
status=json.loads((root/'result/step_status.json').read_text());assert len(status)==11 and all(x['status']=='completed' for x in status.values())
result={'path':str(root),'spots':80,'genes':150,'steps_completed':len(status),'png':len(list((root/'result').rglob('*.png'))),'pdf':len(list((root/'result').rglob('*.pdf'))),'all_pass':True};(root.parent / 'validation_summary.json').write_text(json.dumps(result));print(json.dumps(result))

#!/usr/bin/env python3
"""只读验证模板语法、来源基线和用户预存修改；结果写到指定审计 JSON。"""
from pathlib import Path
import argparse
import ast
import csv
import hashlib
import json
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ('bulkrnaseq','wgbsseq','bulkatacseq','scrnaseq','smallrnaseq','metabolomics','spatial_transcriptomics','singlecell_multiome','spatial_metabolomics')
EXCLUDED = {'__pycache__','.venv','.cache','.pytest_cache','result','runs'}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'template_audit/syntax_and_integrity.json')
    args = parser.parse_args()
    checks = []
    rfiles = []
    for name in TEMPLATES:
        directory = ROOT/name
        if not directory.is_dir():
            checks.append({'kind':'directory','path':name,'status':'FAIL_MISSING'})
            continue
        for path in sorted(directory.rglob('*')):
            if not path.is_file() or any(x in EXCLUDED for x in path.relative_to(directory).parts): continue
            suffix = path.suffix.lower()
            if suffix not in {'.py','.r','.rmd','.sh'}: continue
            record = {'kind':suffix,'path':str(path.relative_to(ROOT)),'status':'PASS'}
            try:
                if suffix == '.py': ast.parse(path.read_text(),filename=str(path))
                elif suffix == '.sh':
                    p = subprocess.run(['bash','-n',str(path)],capture_output=True,text=True)
                    if p.returncode: raise ValueError(p.stderr)
                else:
                    rfiles.append(path)
                    continue
            except (SyntaxError,ValueError) as exc: record.update(status='FAIL',error=str(exc))
            checks.append(record)
    with tempfile.TemporaryDirectory(prefix='template-parse-') as tmp:
        source = Path(tmp)/'files.txt'; source.write_text('\n'.join(map(str,rfiles))+'\n')
        code = r'''
args <- commandArgs(trailingOnly=TRUE)
files <- readLines(args[[1]])
for (f in files) {
  err <- tryCatch({
    if (grepl("[.]Rmd$",f,ignore.case=TRUE)) {
      lines <- readLines(f,warn=FALSE); selected <- character(); inside <- FALSE
      for (line in lines) {
        if (grepl("^```\\{r([ ,}]|$)",line)) {inside<-TRUE;next}
        if (inside && grepl("^```\\s*$",line)) {inside<-FALSE;next}
        if (inside) selected <- c(selected,line)
      }
      parse(text=selected)
    } else parse(f)
    ""
  },error=function(e)conditionMessage(e))
  cat(if(nzchar(err))"FAIL" else "PASS",f,gsub("[\\t\\n]"," ",err),sep="\t");cat("\n")
}
'''
        rscript = Path(tmp)/'parse.R'; rscript.write_text(code)
        p = subprocess.run(['Rscript',str(rscript),str(source)],capture_output=True,text=True)
        if p.returncode: checks.append({'kind':'Rscript','status':'FAIL','error':p.stderr})
        for line in p.stdout.splitlines():
            status, file, error = line.split('\t',2)
            checks.append({'kind':'R/Rmd','path':str(Path(file).relative_to(ROOT)),'status':status,'error':error})
    baseline = json.loads((ROOT/'template_audit/baseline.json').read_text())
    external = [x for x in baseline['sources'] if not Path(x['path']).is_relative_to(ROOT)]
    changed = [x['path'] for x in external if not Path(x['path']).is_file() or hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()!=x['sha256']]
    provenance_errors = []
    provenance_external = {}
    for name in TEMPLATES[:7]:
        manifest = ROOT/name/'provenance.tsv'
        if not manifest.is_file():
            provenance_errors.append(f"MISSING {name}/provenance.tsv")
            continue
        for row in csv.DictReader(manifest.open(), delimiter='\t'):
            path = Path(row['source_path'])
            if not path.is_absolute(): path = ROOT.parent/path
            if path.is_relative_to(ROOT): continue  # Original template hashes intentionally change.
            provenance_external[str(path)] = row['sha256']
    for raw, expected in provenance_external.items():
        path = Path(raw)
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            provenance_errors.append(raw)
    dirty = baseline['preexisting_dirty_file']
    preserved = hashlib.sha256(Path(dirty['path']).read_bytes()).hexdigest()==dirty['sha256']
    diff = subprocess.run(['git','diff','--check'],cwd=ROOT,capture_output=True,text=True)
    payload = {'checks':checks,'syntax_pass':sum(x['status']=='PASS' for x in checks),'syntax_fail':sum(x['status']!='PASS' for x in checks),'external_source_files_checked':len(external),'external_source_changes':changed,'preexisting_dirty_file_preserved':preserved,'provenance_external_files_checked':len(provenance_external),'provenance_errors':provenance_errors,'git_diff_check':diff.returncode,'git_diff_errors':diff.stdout+diff.stderr}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in payload.items() if k!='checks'},ensure_ascii=False))
    return int(bool(payload['syntax_fail'] or changed or not preserved or provenance_errors or diff.returncode))

if __name__=='__main__':
    raise SystemExit(main())

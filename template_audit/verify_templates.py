#!/usr/bin/env python3
"""验证当前模板源码语法；旧提取来源及旧生产验收不作当前运行证明。"""
from pathlib import Path
import argparse
import ast
import csv
import hashlib
import json
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = ('rnaseq','wgbs','atac','scrnaseq','smallrnaseq','metabolomics','spatial_transcriptomics','singlecell_multiome','spatial_metabolomics','integration')
EXCLUDED = {'__pycache__','.venv','.cache','.pytest_cache','.maintenance','result','runs','history','.archive','archive'}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'template_audit/current_syntax.json')
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
    diff = subprocess.run(['git','diff','--check'],cwd=ROOT,capture_output=True,text=True)
    signatures = {x['path']:hashlib.sha256((ROOT/x['path']).read_bytes()).hexdigest()
                  for x in checks if 'path' in x and (ROOT/x['path']).is_file()}
    payload = {'checks':checks,'syntax_pass':sum(x['status']=='PASS' for x in checks),
               'syntax_fail':sum(x['status']!='PASS' for x in checks),
               'source_sha256':signatures,'git_diff_check':diff.returncode,
               'git_diff_errors':diff.stdout+diff.stderr,
               'boundary':'当前语法检查；不证明生产数据重跑。旧来源校验见 history/20260929。'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in payload.items() if k not in {'checks','source_sha256'}},ensure_ascii=False))
    return int(bool(payload['syntax_fail'] or diff.returncode))

if __name__=='__main__':
    raise SystemExit(main())

#!/usr/bin/env python3
"""生成四个物种/endo合成案例，实际验证路径与 reads 分类守恒。"""
import argparse
import csv
import hashlib
import os
from pathlib import Path
import subprocess
import sys


def inventory(root):
    return {p.relative_to(root).as_posix(): (p.stat().st_size, p.stat().st_mtime_ns) for p in root.rglob('*') if p.is_file()}


def rows(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    source = Path(__file__).resolve().parents[1]
    checks = []
    for i, (species, endo) in enumerate([('Mus','0'),('Mus','1'),('Homo','0'),('Homo','1')]):
        root = output / ('root A with spaces' if i % 2 else 'root_B') / f'{species}_{endo}'
        subprocess.run([sys.executable, str(source/'scripts/make_fixture.py'), '--output', str(root), '--species', species, '--endo', endo], check=True)
        run = ['bash', str(root/'run.sh'), '--project', '.']
        env = dict(os.environ, RUN_ENDOSIRNA='1' if endo == '0' else '0', PANDORA_SPECIES='INVALID')
        before = inventory(root)
        commands = [('help', run+['--help'], 0), ('list', run+['--list'], 0), ('describe_check', run+['--stage','describe','--check'], 0), ('upstream_missing_fastq',run+['--stage','upstream','--check'], 2), ('invalid_module',run+['--modules','not_a_module','--check'], 2)]
        for label, command, expected in commands:
            result = subprocess.run(command, cwd=output, env=env, capture_output=True, text=True)
            (output/f'{species}_{endo}_{label}.log').write_text(result.stdout+result.stderr)
            if result.returncode != expected or inventory(root) != before:
                raise AssertionError(f'{species}_{endo}: {label} exit={result.returncode}, readonly={inventory(root)==before}')
            checks.append([species,endo,label,'PASS'])
        modules = 'sncrna,endo' if endo == '1' else 'sncrna'
        # 绝对数据集路径由中心入口调用；实现定位于入口自身，与数据集/调用cwd无关。
        command = ['bash',str(root/'run.sh'),'--project',str(root),'--stage','describe','--modules',modules]
        with (output/f'{species}_{endo}_execution.log').open('w') as log:
            subprocess.run(command,cwd=output,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        stage = root/'result/describe'
        setup = stage / ('01_setup' if species == 'Homo' else '1_setup')
        receipt = rows(stage/'run_status.tsv')[0]
        if receipt['State'] != 'completed' or receipt['Project'] != 'synthetic_fixture':
            raise AssertionError('stage not completed or explicit project identity lost')
        for table in ['composition_integrity_check','length_composition_integrity_check','output_integrity_check']:
            if any(row['Pass'] != 'TRUE' for row in rows(setup/f'{table}.tsv')):
                raise AssertionError(table)
        counts = rows(setup/'1_expression_matrix/class_reads.tsv')
        by_class = {row['MainCategory']: float(row['synthetic_library_1']) for row in counts}
        if sum(by_class.values()) != 65 or by_class.get('endo-siRNA',0) != (30 if endo == '1' else 0):
            raise AssertionError(f'classification boundary failed: {by_class}')
        if by_class.get('miRNA') != 20:
            raise AssertionError('explicit miRNA conflict was reclassified')
        for entry in rows(stage/'input_signature.tsv'):
            path = Path(entry['Path'])
            if path.stat().st_size != int(entry['Size']) or hashlib.md5(path.read_bytes()).hexdigest() != entry['MD5']:
                raise AssertionError('signature mismatch: '+str(path))
        checks.extend([[species,endo,'actual_sncrna_endo_execution','PASS'],[species,endo,'exclusive_reads_and_conflict_boundary','PASS'],[species,endo,'current_source_signature','PASS']])
        print(species,endo,'completed',flush=True)
    with (output/'checks.tsv').open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['Species','RUN_ENDOSIRNA','Check','Status']);writer.writerows(checks)


if __name__ == '__main__':
    main()

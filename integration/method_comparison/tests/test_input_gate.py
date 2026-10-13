"""小型输入门fixture；不声称真实数据复现。"""
import contextlib
import csv
import hashlib
import importlib.util
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

METHOD_ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('method_gate',METHOD_ROOT/'scripts/check_inputs.py')
gate=importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)
report_spec=importlib.util.spec_from_file_location('method_report',METHOD_ROOT/'scripts/build_report.py')
report=importlib.util.module_from_spec(report_spec)
report_spec.loader.exec_module(report)


class InputGate(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='method input gate ')
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        source=METHOD_ROOT/'../smallrnaseq/config/samples.tsv'
        if source.is_file():
            with source.open() as f:self.metadata=list(csv.DictReader(f,delimiter='\t'))
        else:
            self.metadata=[{'sample':f'library_{i}','species':'Homo_sapiens','analysis_role':'method_comparison' if i<10 else 'spikein_pair','biological_sample':'2' if i<10 else '3','method':str(1 if i<3 else 2 if i<6 else i-3),'pair_id':f'pair_{i%3}' if i<6 else ''} for i in range(12)]
        stage=self.root/'result/describe'
        (stage/'01_setup/1_expression_matrix').mkdir(parents=True)
        (stage/'02_single/tables').mkdir(parents=True)
        (stage/'run_status.tsv').write_text('State\ncompleted\n')
        (stage/'module_status.tsv').write_text('Module\tStatus\nqc\tcompleted\n')
        (stage/'01_setup/1_expression_matrix/class_reads.tsv').write_text('fixture\n')
        (stage/'02_single/tables/single_sample_qc_table.csv').write_text('fixture\n')
        (self.root/'config').mkdir()
        self.write_metadata()
        for row in self.metadata:
            sample=self.root/'sportsresult'/row['sample'];sample.mkdir(parents=True)
            for suffix in ('_summary.txt','_output.txt','_length_distribution.txt'):(sample/(row['sample']+suffix)).write_text('fixture\n')

    def write_metadata(self):
        path=self.root/'config/samples.tsv'
        with path.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(self.metadata[0]),delimiter='\t');writer.writeheader();writer.writerows(self.metadata)
        (self.root/'result/describe/input_signature.tsv').write_text(f'Path\tSize\tMD5\n{path}\t{path.stat().st_size}\t{hashlib.md5(path.read_bytes()).hexdigest()}\n')

    def run_gate(self):
        with patch.dict(os.environ,{'PANDORA_PROJECT_ROOT':str(self.root)}),contextlib.redirect_stdout(io.StringIO()):gate.main()

    def test_current_input_passes_without_writing(self):
        before={p.relative_to(self.root):(p.stat().st_size,p.stat().st_mtime_ns) for p in self.root.rglob('*') if p.is_file()}
        self.run_gate()
        after={p.relative_to(self.root):(p.stat().st_size,p.stat().st_mtime_ns) for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(before,after)

    def test_incomplete_producer_is_rejected(self):
        (self.root/"result/describe/run_status.tsv").write_text("State\nrunning\n")
        with self.assertRaisesRegex(ValueError,"not_run"):self.run_gate()

    def test_changed_producer_input_is_stale(self):
        with (self.root/'config/samples.tsv').open('a') as f:f.write('\n')
        with self.assertRaisesRegex(ValueError,'STALE_OUTPUTS'):self.run_gate()

    def test_wrong_species_is_rejected(self):
        self.metadata[0]['species']='Mus_musculus';self.write_metadata()
        with self.assertRaises(ValueError):self.run_gate()

    def test_duplicate_pair_member_is_rejected(self):
        self.metadata[0]['pair_id']=self.metadata[1]['pair_id'];self.write_metadata()
        with self.assertRaises(ValueError):self.run_gate()

    def test_report_rejects_running_analysis_and_records_failure(self):
        (self.root/'result/describe/run_status.tsv').write_text('State\nrunning\n')
        with patch.dict(os.environ,{'METHOD_COMPARISON_ROOT':str(self.root)}):
            with self.assertRaisesRegex(ValueError,'not_run'):report.main()
        import json
        receipt=json.loads((self.root/'result/synthesis/report_receipt.json').read_text())
        self.assertEqual(receipt['state'],'failed')
        self.assertFalse((self.root/'README.md').exists())

    def test_report_rejects_changed_signed_input(self):
        with (self.root/'config/samples.tsv').open('a') as handle:handle.write('\n')
        with patch.dict(os.environ,{'METHOD_COMPARISON_ROOT':str(self.root)}):
            with self.assertRaisesRegex(ValueError,'STALE_OUTPUTS'):report.main()
        import json
        receipt=json.loads((self.root/'result/synthesis/report_receipt.json').read_text())
        self.assertEqual(receipt['state'],'failed')
        self.assertFalse((self.root/'README.md').exists())


if __name__=='__main__':unittest.main()

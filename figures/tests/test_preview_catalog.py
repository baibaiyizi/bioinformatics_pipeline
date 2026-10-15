"""Test catalog-driven navigation using copies of existing frozen receipts."""
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('update_previews', ROOT / 'scripts/update_previews.py')
updater = importlib.util.module_from_spec(spec)
spec.loader.exec_module(updater)


class PreviewCatalogTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='sci-preview-', dir='/tmp')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'figures'
        shutil.copytree(ROOT, self.root, ignore=shutil.ignore_patterns('__pycache__', '.cache', 'tests'))

    def rows(self):
        return updater.read_catalog(self.root)

    def write_rows(self, rows):
        with (self.root / 'recipe_catalog.tsv').open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)

    def navigation(self):
        return {name: (self.root / 'previews' / name).read_bytes()
                for name in ('index.html', 'manifest.json', 'validation.json')}

    def test_arbitrary_cwd_idempotent_and_preserves_evidence(self):
        preserved = {p: updater.sha(p) for p in self.root.rglob('*') if p.is_file()
                     and (p.suffix in {'.png', '.pdf'} or p.name == 'figure.json'
                          or p.parent.name in {'scripts', 'examples', 'config'})}
        command = [sys.executable, str(self.root / 'scripts/update_previews.py')]
        subprocess.run(command, cwd=self.temp.name, check=True, capture_output=True)
        first = self.navigation()
        subprocess.run(command, cwd='/', check=True, capture_output=True)
        self.assertEqual(first, self.navigation())
        self.assertEqual(preserved, {p: updater.sha(p) for p in preserved})
        manifest = json.loads(first['manifest.json'])
        self.assertEqual(manifest['recipe_count'], 19)
        self.assertEqual(manifest['figure_pairs'], 16)
        self.assertEqual({r['family'] for r in manifest['unrendered_recipes']},
                         {'effect', 'effect_heatmap', 'embedding'})
        self.assertTrue(all(r['execution_state'] == 'not_run' for r in manifest['unrendered_recipes']))

    def test_catalog_text_drives_html_and_is_escaped(self):
        rows = self.rows()
        rows[0]['preview_title'] = '<script>alert("x")</script>'
        rows[0]['question'] = "Question & 'quoted'"
        rows[0]['preview_distribution'] = 'Changed source distribution'
        rows[0]['limitations'] = 'Changed evidence boundary'
        self.write_rows(rows)
        manifest = updater.build(self.root)
        page = (self.root / 'previews/index.html').read_text()
        self.assertNotIn('<script>', page)
        self.assertIn('&lt;script&gt;', page)
        self.assertIn('Question &amp; &#x27;quoted&#x27;', page)
        self.assertIn('Changed source distribution', page)
        self.assertIn('Changed evidence boundary', page)
        self.assertEqual(manifest['figures'][0]['title'], rows[0]['preview_title'])

    def test_not_run_cannot_borrow_preview(self):
        rows = self.rows()
        rows[0]['execution_state'] = 'not_run'
        self.write_rows(rows)
        before = self.navigation()
        with self.assertRaisesRegex(ValueError, '未执行条目不能附带预览'):
            updater.build(self.root)
        self.assertEqual(before, self.navigation())

    def test_stale_entrypoint_receipt_rejected_before_navigation_write(self):
        script = self.root / 'scripts/render.py'
        script.write_text(script.read_text() + '\n# changed producer\n')
        before = self.navigation()
        with self.assertRaisesRegex(ValueError, '回执 hash 不符'):
            updater.build(self.root)
        self.assertEqual(before, self.navigation())

    def test_stale_source_receipt_rejected(self):
        source = self.root / 'previews/python_volcano/source.tsv'
        source.write_text(source.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, '回执 hash 不符'):
            updater.build(self.root)

    def test_preview_path_cannot_escape_template(self):
        rows = self.rows()
        rows[0]['preview'] = '../../other/figure.png'
        self.write_rows(rows)
        with self.assertRaisesRegex(ValueError, '路径越界'):
            updater.build(self.root)

    def test_method_links_and_frozen_recipe_inputs(self):
        for row in updater.read_catalog(ROOT):
            filename, anchor = row['method_ref'].split('#', 1)
            target = (ROOT / filename).resolve()
            self.assertTrue(target.is_file())
            self.assertIn(f'<a id="{anchor}"></a>', target.read_text())
        rows = self.rows()
        locus = next(r for r in rows if r['family'] == 'locus')
        self.assertIn('GRCm39', locus['prerequisites'])
        self.assertIn('人源', locus['limitations'])
        composition = next(r for r in rows if r['family'] == 'composition')
        self.assertEqual(composition['inputs'], 'library_id,cell_type,proportion')


if __name__ == '__main__':
    unittest.main()

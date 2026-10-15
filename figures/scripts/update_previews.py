#!/usr/bin/env python3
"""从 recipe_catalog.tsv 和既有成图回执生成图册；不执行绘图或统计。"""
from pathlib import Path
import csv
import hashlib
import html
import json
import shlex

import fitz
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local_path(root, value):
    """Catalog artifacts must remain inside this template, including symlinks."""
    path = (root / value).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f'模板路径越界：{value}')
    return path


def read_catalog(root):
    with (root / 'recipe_catalog.tsv').open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'family', 'variant', 'entrypoint', 'inputs', 'meaning', 'question',
                    'prerequisites', 'limitations', 'method_ref', 'validation_state',
                    'execution_state', 'preview', 'preview_title',
                    'preview_distribution', 'preview_data_kind'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('recipe_catalog 缺少必要字段')
        rows = list(reader)
    if not rows:
        raise ValueError('recipe_catalog 不能为空')
    seen, previews = set(), set()
    for row in rows:
        key = (row['family'], row['variant'])
        if key in seen or any(not row[name] for name in required - {
                'preview', 'preview_title', 'preview_distribution', 'preview_data_kind'}):
            raise ValueError(f'重复或不完整 recipe：{key}')
        seen.add(key)
        if not local_path(root, shlex.split(row['entrypoint'])[0]).is_file():
            raise ValueError(f'入口不存在：{key}')
        if row['preview']:
            if row['execution_state'] != 'COMPLETED' or row['validation_state'] not in {
                    'example_executed', 'project_validated'}:
                raise ValueError(f'未执行条目不能附带预览：{key}')
            preview = local_path(root, row['preview'])
            if preview in previews or preview.name != 'figure.png' or not preview.is_relative_to(root / 'previews'):
                raise ValueError(f'重复或无效预览路径：{key}')
            previews.add(preview)
            if not all(row[name] for name in ('preview_title', 'preview_distribution', 'preview_data_kind')):
                raise ValueError(f'预览缺少说明或来源：{key}')
        elif any(row[name] for name in ('preview_title', 'preview_distribution', 'preview_data_kind')):
            raise ValueError(f'无图条目不能附带预览说明：{key}')
    return rows


def check_receipt(root, row, folder):
    receipt = json.loads((folder / 'figure.json').read_text(encoding='utf-8'))
    if receipt['status'] != 'COMPLETED' or receipt['statistics_recomputed'] is not False:
        raise ValueError(f'{folder.name}：回执未完成或重算了统计')
    checks = {
        local_path(root, shlex.split(row['entrypoint'])[0]): receipt['script_sha256'],
        root / 'config/figure_style.json': receipt['profile_sha256'],
    }
    # These are the three existing producer receipt formats, not migration fallbacks.
    if 'output_sha256' in receipt:
        checks.update({folder / name: receipt['output_sha256'][name]
                       for name in ('figure.pdf', 'figure.png')})
        checks[root / 'scripts/render_examples.R'] = receipt['runner_sha256']
        checks[root / 'scripts/figure_style.R'] = receipt['style_sha256']
        checks.update({local_path(root, 'examples/' + name): digest
                       for name, digest in receipt['input_sha256'].items()})
    else:
        for suffix in ('pdf', 'png'):
            output = receipt['outputs'][suffix]
            checks[folder / f'figure.{suffix}'] = output['sha256'] if isinstance(output, dict) else output
        inputs = receipt['input_sha256']
        if isinstance(inputs, str):
            checks[folder / 'source.tsv'] = inputs
        else:
            checks.update({local_path(folder, name): digest for name, digest in inputs.items()})
    for path, expected in checks.items():
        if not path.is_file() or sha(path) != expected:
            raise ValueError(f'{folder.name}：回执 hash 不符：{path.name}')
    if 'data_kind' in receipt and receipt['data_kind'] != row['preview_data_kind']:
        raise ValueError(f'{folder.name}：来源分类与回执不符')
    return receipt


def build(root=ROOT):
    root = root.resolve()
    rows = read_catalog(root)
    figures, unrendered, issues = [], [], []
    escape = html.escape
    lines = ["<!doctype html><html lang='zh'><meta charset='utf-8'><title>生信绘图模板预览</title>",
             "<style>body{font:16px sans-serif;max-width:1100px;margin:auto;padding:28px;color:#263445}img{max-width:100%;border:1px solid #e4e7eb}figure{margin:40px 0}p{line-height:1.6}</style>",
             "<h1>生信绘图模板：来源与执行状态</h1><p>目录由 recipe_catalog.tsv 生成；只有目录与真实回执一致的已执行图件显示预览。示例执行不代表上游论文复现或项目科学验证。<a href='../README.md'>使用说明</a> · <a href='../recipe_catalog.tsv'>完整目录与输入前提</a></p>"]
    for row in rows:
        if not row['preview']:
            unrendered.append({name: row[name] for name in (
                'family', 'variant', 'entrypoint', 'validation_state', 'execution_state')})
            continue
        png = local_path(root, row['preview'])
        folder, pdf = png.parent, png.with_suffix('.pdf')
        receipt = check_receipt(root, row, folder)
        key = folder.name
        with Image.open(png) as im:
            media = {'mode': im.mode, 'size': list(im.size), 'dpi': list(im.info.get('dpi', ()))}
            if im.mode != 'RGB' or abs(im.width - 183 / 25.4 * 300) > 2 or len(media['dpi']) != 2 or any(abs(dpi - 300) > 1 for dpi in media['dpi']):
                issues.append(f'{key}:PNG规格')
        with fitz.open(pdf) as doc:
            font_records = []
            for font in doc.get_page_fonts(0):
                extracted = doc.extract_font(font[0])
                font_records.append({'name': font[3], 'type': font[2], 'embedded': bool(extracted[3])})
            text, width_mm = doc[0].get_text(), doc[0].rect.width * 25.4 / 72
            # Cairo rounds its PDF media box to whole points (518pt = 182.74mm).
            if len(doc) != 1 or abs(width_mm - 183) > .4 or not font_records or not all(f['embedded'] for f in font_records) or len(text) < 30:
                issues.append(f'{key}:PDF规格/字体/文本')
        figures.append({'recipe': key, 'title': row['preview_title'],
                        'validation_state': row['validation_state'], 'data_kind': row['preview_data_kind'],
                        'png': str(png.relative_to(root / 'previews')), 'pdf': str(pdf.relative_to(root / 'previews')),
                        'png_sha256': sha(png), 'pdf_sha256': sha(pdf), 'profile_sha256': receipt['profile_sha256'],
                        'statistics_recomputed': False, 'media': media, 'pdf_width_mm': width_mm,
                        'pdf_fonts': font_records, 'pdf_text_characters': len(text)})
        lines.append(f"<figure id='{escape(key)}'><h2>{escape(row['preview_title'])}</h2>"
                     f"<img src='{escape(key)}/figure.png' alt='{escape(row['preview_title'])}'><figcaption>"
                     f"<p>生物学意义：{escape(row['question'])}</p>"
                     f"<p>数据分布特征：{escape(row['preview_distribution'])}</p>"
                     f"<p>最终结论：{escape(row['limitations'])}</p>"
                     f"<p>输入前提：{escape(row['prerequisites'])}</p>"
                     f"<p>来源：{escape(row['preview_data_kind'])}；验证状态：{escape(row['validation_state'])}</p>"
                     f"<a href='{escape(key)}/figure.pdf'>原生PDF</a> · <a href='{escape(key)}/figure.json'>执行回执</a></figcaption></figure>")
    lines.append('<h2>尚无预览的条目</h2><ul>')
    for row in unrendered:
        label = ' / '.join(row[name] for name in ('family', 'variant', 'validation_state', 'execution_state'))
        lines.append(f'<li>{escape(label)}；未附预览</li>')
    lines.append('</ul></html>')
    manifest = {'status': 'COMPLETED' if not issues else 'FAILED', 'recipe_count': len(rows),
                'figure_pairs': len(figures), 'figures': figures, 'unrendered_recipes': unrendered}
    validation = {'status': 'PASS' if not issues else 'FAIL', 'figure_pairs': len(figures),
                  'issues': issues, 'receipt_hashes_checked': True, 'scientific_model_rerun': False,
                  'manual_visual_review': 'not_run_by_updater；自动规格通过不等同视觉或投稿通过'}
    # All receipts are checked before replacing any derived navigation files.
    (root / 'previews/index.html').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    for filename, record in (('manifest.json', manifest), ('validation.json', validation)):
        (root / 'previews' / filename).write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    if issues:
        raise ValueError(issues)
    return manifest


if __name__ == '__main__':
    result = build()
    print(f"PASS {result['recipe_count']} recipes; {result['figure_pairs']} preview pairs; receipts, RGB300dpi, 183mm PDF, embedded fonts and text.")

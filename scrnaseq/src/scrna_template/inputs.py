from __future__ import annotations
import gzip
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import pandas as pd

class InputGateError(RuntimeError):
    pass

def parse_md5_manifest(path):
    entries = {}
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        match = re.match(r'^([0-9a-fA-F]{32})\s+[* ]?(.+?)\s*$', line)
        bsd = re.match(r'^MD5\s*\((.+)\)\s*=\s*([0-9a-fA-F]{32})\s*$', line)
        if match:
            digest, name = match.groups()
        elif bsd:
            name, digest = bsd.groups()
        else:
            raise InputGateError('MD5 格式无效')
        basename = Path(name.replace('\\', '/')).name
        if basename in entries:
            raise InputGateError('MD5 basename 重复，无法唯一映射')
        entries[basename] = digest.lower()
    if not entries:
        raise InputGateError('MD5 清单为空')
    return entries

def md5sum(path):
    digest = hashlib.md5(usedforsecurity=False)
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(16*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()

def fastq_records(path):
    with gzip.open(path, 'rt') as stream:
        while True:
            header = stream.readline()
            if not header:
                return
            sequence, plus, quality = [stream.readline() for _ in range(3)]
            if not header.startswith('@') or not plus.startswith('+') or not quality or len(sequence.rstrip()) != len(quality.rstrip()):
                raise InputGateError(f'FASTQ 结构损坏: {path.name}')
            yield re.sub(r'/[12]$', '', header.split()[0][1:])

def full_pair_check(r1, r2, i1=None):
    from itertools import zip_longest
    streams = [fastq_records(r1), fastq_records(r2)]
    if i1 is not None:
        streams.append(fastq_records(i1))
    n = 0
    for names in zip_longest(*streams):
        if None in names or len(set(names)) != 1:
            raise InputGateError(f'FASTQ 配对/记录数不一致: {r1.name}')
        n += 1
    if n == 0:
        raise InputGateError('FASTQ 为空')
    return n

def fastq_pairs(config):
    pairs = []
    for row in config.samples().itertuples():
        paths = {key: [config.resolve(p.strip()) for p in getattr(row, key).split(';') if p.strip()] for key in ('r1', 'r2', 'i1')}
        if len(paths['r1']) != len(paths['r2']) or (paths['i1'] and len(paths['i1']) != len(paths['r1'])):
            raise InputGateError(f'{row.library_id}: r1/r2/i1 分号分隔路径数量不匹配')
        for idx, (left, right) in enumerate(zip(paths['r1'], paths['r2'], strict=True)):
            pairs.append((row.library_id, idx+1, left, right, paths['i1'][idx] if paths['i1'] else None))
    all_paths = [p for *_, left, right, index in pairs for p in (left, right, index) if p is not None]
    if len(set(all_paths)) != len(all_paths):
        raise InputGateError('同一 FASTQ 被多个文库/读段重复指定')
    if len({p.name for p in all_paths}) != len(all_paths):
        raise InputGateError('FASTQ basename 重复，MD5 清单无法唯一对应')
    return pairs

def input_paths(config):
    if config.data['inputs']['mode'] == 'h5ad':
        paths = [config.resolve(p) for p in config.samples().matrix]
        if len(set(paths)) != len(paths):
            raise InputGateError('多个 library_id 指向同一个矩阵，可能重复计数')
        return paths
    pairs = fastq_pairs(config)
    return [p for _, _, left, right, index in pairs for p in (left, right, index) if p is not None] + [config.resolve(config.data['inputs']['md5_manifest'])]

def validate_inputs(config, full=False):
    """Read-only, including full checks. No audit/result directories are created."""
    paths = input_paths(config)
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise InputGateError(f'MISSING/空输入: {path}')
    config.markers()
    if config.data['inputs']['mode'] == 'h5ad':
        return pd.DataFrame([{'path': str(p), 'status': 'PRESENT', 'scope': 'integer counts checked during QC'} for p in paths])
    manifest = paths[-1]
    entries = parse_md5_manifest(manifest)
    if set(entries) != {p.name for p in paths[:-1]}:
        raise InputGateError('MD5 条目必须恰好覆盖配置的 FASTQ，包括可选 I1')
    minimum = float(config.data['inputs']['stability_seconds'])
    if minimum < 0:
        raise InputGateError('stability_seconds 不得为负')
    states = {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in paths}
    if any(time.time() - p.stat().st_mtime < minimum for p in paths):
        raise InputGateError('FASTQ/MD5 最近仍变化，尚未达到稳定时间门')
    if config.data['inputs'].get('require_no_open_writer', True):
        if shutil.which('lsof') is None:
            raise InputGateError('缺少 lsof，无法检查活动写句柄')
        for path in paths:
            process = subprocess.run(['lsof', '-F', 'fa', '--', str(path)], text=True, capture_output=True)
            if process.returncode not in (0,1) or process.stderr.strip():
                raise InputGateError(f'lsof 检查不能确认成功: {process.stderr}')
            if any(line in {'aw', 'au'} for line in process.stdout.splitlines()):
                raise InputGateError(f'存在活动写句柄: {path}')
    if not full:
        return pd.DataFrame([{'path': str(p), 'status': 'PRESENT', 'scope': 'MD5/pair audit not_run'} for p in paths])
    rows = []
    for path in paths[:-1]:
        digest = md5sum(path)
        if digest != entries[path.name]:
            raise InputGateError(f'MD5 FAIL: {path}')
        rows.append({'path': str(path), 'md5': digest, 'status': 'PASS'})
    for library, _, left, right, index in fastq_pairs(config):
        n = full_pair_check(left, right, index)
        rows.append({'path': f'{library}:{left.name}', 'read_pairs': n, 'status': 'PASS'})
    if any(states[p] != (p.stat().st_size, p.stat().st_mtime_ns) for p in paths):
        raise InputGateError('校验期间输入大小/mtime 改变')
    return pd.DataFrame(rows)

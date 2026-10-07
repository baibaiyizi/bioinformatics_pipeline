from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re
import pandas as pd
import yaml

UNKNOWN = {'', 'unknown', 'na', 'nan', 'none'}
REQUIRED_SAMPLE_COLUMNS = ('library_id', 'biological_sample_id', 'unit_type', 'condition', 'sample_batch', 'chemistry', 'matrix', 'r1', 'r2', 'i1')

@dataclass(frozen=True)
class WorkflowConfig:
    source: Path
    data: dict

    @classmethod
    def load(cls, source):
        path = Path(source).expanduser().resolve()
        data = yaml.safe_load(path.read_text())
        if not isinstance(data, dict):
            raise ValueError('project.yml 必须是映射')
        cfg = cls(path, data)
        for key in ('project', 'samples', 'inputs', 'qc', 'preprocessing', 'statistics', 'annotation', 'report'):
            if key not in data:
                raise ValueError(f'缺配置 {key}')
        if data['inputs']['mode'] not in {'h5ad', 'fastq'}:
            raise ValueError('inputs.mode 只能是 h5ad 或 fastq')
        if data['statistics']['design'] != 'independent_two_group':
            raise ValueError('仅支持独立两组；配对、重复测量设计未实现')
        if data['statistics']['case'] == data['statistics']['control']:
            raise ValueError('case 与 control 必须不同')
        if int(data['statistics']['min_independent_samples_per_group']) < 3:
            raise ValueError('推断每组至少 3 个独立生物学样本')
        cfg.samples()
        return cfg

    @property
    def root(self):
        p = Path(self.data['project']['root']).expanduser()
        return p.resolve() if p.is_absolute() else (self.source.parent / p).resolve()

    def resolve(self, value):
        p = Path(value).expanduser()
        return p.resolve() if p.is_absolute() else (self.root / p).resolve()

    @property
    def run_root(self):
        return self.resolve(self.data['project'].get('run_root', 'runs'))

    @property
    def samples_path(self):
        return self.resolve(self.data['samples'])

    @property
    def markers_path(self):
        return self.resolve(self.data['markers'])

    @property
    def reference(self):
        return self.resolve(self.data['reference']['transcriptome'])

    def markers(self):
        value = yaml.safe_load(self.markers_path.read_text())
        if not isinstance(value, dict) or 'cell_types' not in value:
            raise ValueError('markers.yml 必须包含 cell_types 映射')
        return value

    def samples(self):
        frame = pd.read_csv(self.samples_path, sep='\t', dtype=str, keep_default_na=False).apply(lambda s: s.str.strip())
        missing = set(REQUIRED_SAMPLE_COLUMNS) - set(frame.columns)
        if missing or frame.empty:
            raise ValueError(f'samples.tsv 为空或缺字段: {sorted(missing)}')
        if frame.library_id.duplicated().any() or not frame.library_id.map(lambda s: bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', s))).all():
            raise ValueError('library_id 必须唯一且只能包含字母数字下划线点和短横线')
        expected = {str(self.data['statistics']['case']), str(self.data['statistics']['control'])}
        if set(frame.condition) != expected:
            raise ValueError('样本 condition 必须恰好覆盖显式 case/control')
        if not set(frame.unit_type) <= {'individual', 'pool', 'unknown', ''}:
            raise ValueError('unit_type 只能是 individual、pool、unknown；混池计作一个单位')
        for column in ('subject_id', 'pair_id', 'timepoint'):
            if column in frame and (~frame[column].str.lower().isin(UNKNOWN)).any():
                raise ValueError(f'{column} 提示配对/重复测量；本模板未实现该设计')
        mode = self.data['inputs']['mode']
        for column in (('matrix',) if mode == 'h5ad' else ('r1', 'r2', 'chemistry')):
            if frame[column].str.lower().isin(UNKNOWN).any():
                raise ValueError(f'{mode} 输入必须提供 {column}')
        return frame

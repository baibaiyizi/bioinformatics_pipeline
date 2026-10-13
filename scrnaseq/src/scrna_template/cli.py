from __future__ import annotations
import argparse
import json

def main():
    parser = argparse.ArgumentParser(description='可复制的 scRNA-seq 独立两组模板')
    parser.add_argument('--config', default='config/project.yml')
    parser.add_argument('--list-steps', action='store_true')
    sub = parser.add_subparsers(dest='command')
    validate = sub.add_parser('validate', help='只读输入和设计预检查')
    validate.add_argument('--full-input-audit', action='store_true')
    run = sub.add_parser('run', help='执行核心流程；不接受已有 run-id')
    run.add_argument('--run-id', required=True)
    plots = sub.add_parser('plots', help='从已签名数值源重绘完整图谱，不调用统计')
    plots.add_argument('--source-run', required=True)
    plots.add_argument('--check', action='store_true')
    fixture = sub.add_parser('fixture', help='生成独立合成输入/配置，不能写入已存在目录')
    fixture.add_argument('--output', required=True)
    fixture.add_argument('--per-group', type=int, default=3)
    fixture.add_argument('--technical-libraries', type=int, default=1)
    args = parser.parse_args()
    if args.list_steps:
        print('inputs -> cellranger (FASTQ only) -> qc -> atlas -> annotation -> descriptive -> optional pseudobulk inference -> report -> audit')
        return 0
    if args.command == 'fixture':
        from .fixtures import make_fixture
        print(make_fixture(args.output, args.per_group, args.technical_libraries))
        return 0
    if args.command is None:
        parser.print_help()
        return 0
    from .config import WorkflowConfig
    from .inputs import validate_inputs
    from .eligibility import assess_eligibility
    config = WorkflowConfig.load(args.config)
    if args.command == 'validate':
        audit = validate_inputs(config, full=args.full_input_audit)
        stats = config.data['statistics']
        print(audit.to_string(index=False))
        print(json.dumps(assess_eligibility(config.samples(), stats['min_independent_samples_per_group'], stats['case'], stats['control']), ensure_ascii=False, indent=2))
        return 0
    if args.command == 'plots':
        from .replot import replot
        print(json.dumps(replot(config, args.source_run, args.check), ensure_ascii=False))
        return 0
    from .pipeline import run_pipeline
    print(run_pipeline(config, args.run_id).run_dir)
    return 0

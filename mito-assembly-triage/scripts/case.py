#!/usr/bin/env python3
"""任务内证据记录与报告；不执行分析命令，不访问跨任务经验库。

init DIR --case-id ID --issue TEXT --input ROLE PATH [--input ROLE PATH ...]
input DIR --role ROLE --path PATH
event DIR --record event.json
candidate DIR --record candidate.json
candidate-check DIR --record check.json
conclusion DIR --record conclusion.json [--update] [--expected-revision N]
validate DIR [--verify-files]
report DIR [--output report.md]

退出码：0 完成；1 输入/格式/关联/文件校验失败；3 缺少 jsonschema。
版本 3 与历史 case schema 分开，不隐式迁移。
候选是已有工具生成的文件；candidate/candidate-check 只登记关系与复核结果。
"""
from pathlib import Path
import sys

from _case_store import CaseStore, MissingDependency, read_json, snapshot, files_in
from _case_report import render
from _evidence_io import InputArgumentParser, atomic_write_text, protect_output


def main(argv=None):
    parser = InputArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    init = commands.add_parser('init')
    init.add_argument('directory')
    init.add_argument('--case-id', required=True)
    init.add_argument('--issue', required=True)
    init.add_argument('--taxon', default='未确定')
    init.add_argument('--case-type', choices=('abnormal_case', 'normal_validation_case', 'tool_failure_case'),
                      default='abnormal_case')
    init.add_argument('--input', nargs=2, action='append', required=True, metavar=('ROLE', 'PATH'))
    add_input = commands.add_parser('input')
    add_input.add_argument('directory')
    add_input.add_argument('--role', required=True)
    add_input.add_argument('--path', required=True)
    add_input.add_argument('--expected-revision', type=int)
    for action in ('event', 'conclusion', 'candidate', 'candidate-check'):
        command = commands.add_parser(action)
        command.add_argument('directory')
        command.add_argument('--record', required=True, type=Path)
        command.add_argument('--expected-revision', type=int)
        if action == 'conclusion':
            command.add_argument('--update', action='store_true')
    check = commands.add_parser('validate')
    check.add_argument('directory')
    check.add_argument('--verify-files', action='store_true')
    report = commands.add_parser('report')
    report.add_argument('directory')
    report.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    store = CaseStore(args.directory)
    try:
        from _deps import require_stage
        require_stage('case_records', __file__)
        if args.action == 'init':
            case = store.init(args.case_id, args.issue, args.taxon, args.case_type, args.input)
        elif args.action in ('input', 'event', 'conclusion', 'candidate', 'candidate-check'):
            record = dict(role=args.role, path=args.path) if args.action == 'input' else read_json(args.record)
            case = store.mutate(args.action, record, update=getattr(args, 'update', False),
                                expected_revision=args.expected_revision)
        else:
            with store.locked():
                case = store.load(verify_files=args.action == 'report' or args.verify_files)
                if args.action == 'report':
                    output = args.output if args.output is not None else store.directory/'report.md'
                    protect_output(output, [store.path, store.lockpath] + [item['path'] for item in files_in(case)])
                    case_hash = snapshot(store.path, 'case', 'ledger')['sha256']
                    atomic_write_text(output, render(case, case_hash))
                    print('报告已生成: %s' % output)
                else:
                    print('FORMAT/BUSINESS: VALID；文件校验: %s；不证明科学结论。' %
                          ('已完成' if args.verify_files else '未执行（加 --verify-files）'))
        print('案例 %s 修订 %d；输入 %d，事件 %d，候选 %d，结论 %d' %
              (case['case_id'], case['revision'], len(case['inputs']), len(case['events']),
               len(case.get('candidates', [])), len(case['conclusions'])))
        if args.action in ('init', 'input'):
            for item in case['inputs']:
                print('%s\t%s\t%s\t%s' % (item['id'], item['role'], item['path'], item['sha256']))
        return 0
    except MissingDependency as exc:
        print('环境错误: %s' % exc, file=sys.stderr)
        return 3
    except (OSError, ValueError, TypeError) as exc:
        print('ERROR: %s' % exc, file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

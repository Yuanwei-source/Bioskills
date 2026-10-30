#!/usr/bin/env python3
"""本地跨样本经验库；只保存经验与案例引用，不复制样本数据或修改 Skill。

propose   从版本 3 案例和 JSON 草稿生成 candidate（须 --authorize-local）
review    审核当前修订（须 --authorize-local）
search    先按状态/领域/范围过滤，再用 FTS5 检索摘要
show      按经验 ID 读取完整修订
supersede 追加新修订并保留历史（须 --authorize-local）
rebuild   从 JSON 事实源重建派生索引（须 --authorize-local）

经验库由 --library 或 MITO_KNOWLEDGE_DIR 显式指定。
退出码：0 完成；1 输入/校验/状态错误；3 缺少依赖。
"""
import json
import os
from pathlib import Path
import sqlite3
import sys

from _case_store import MissingDependency
from _evidence_io import InputArgumentParser
from _lesson_store import LessonStore


def parser():
    result = InputArgumentParser(description=__doc__)
    result.add_argument('--library', default=os.environ.get('MITO_KNOWLEDGE_DIR'),
                        help='本地经验库目录；也可由 MITO_KNOWLEDGE_DIR 指定')
    commands = result.add_subparsers(dest='action', required=True)
    for action in ('propose', 'supersede'):
        command = commands.add_parser(action)
        if action == 'supersede':
            command.add_argument('lesson_id')
        command.add_argument('--case', required=True, help='版本 3 案例目录')
        command.add_argument('--record', required=True, type=Path, help='经验草稿 JSON')
        command.add_argument('--allow-unresolved', action='store_true')
        command.add_argument('--authorize-local', action='store_true')
    review = commands.add_parser('review')
    review.add_argument('lesson_id')
    review.add_argument('--status', required=True,
                        choices=('verified', 'rejected', 'deprecated', 'withdrawn'))
    review.add_argument('--reviewer', required=True)
    review.add_argument('--reason', required=True)
    review.add_argument('--authorize-local', action='store_true')
    search = commands.add_parser('search')
    search.add_argument('--query', default='')
    search.add_argument('--status', default='verified',
                        choices=('verified', 'candidate', 'rejected', 'deprecated', 'withdrawn', 'all'))
    search.add_argument('--domain', choices=('tool', 'annotation', 'biology'))
    search.add_argument('--taxon')
    search.add_argument('--data-type')
    search.add_argument('--limit', type=int, default=5)
    show = commands.add_parser('show')
    show.add_argument('lesson_id')
    show.add_argument('--revision', type=int)
    rebuild = commands.add_parser('rebuild')
    rebuild.add_argument('--authorize-local', action='store_true')
    return result


def print_cards(cards):
    if not cards:
        print('未找到匹配经验；这不证明没有相关案例。')
        return
    for card in cards:
        print('%s rev%d [%s/%s] %s' % (card['lesson_id'], card['revision'],
                                      card['status'], card['domain'], card['title']))
        print('  范围: taxa=%s; data_types=%s; source=%s' % (
            ','.join(card['taxa']), ','.join(card['data_types']), card['source_case_type']))
        print('  线索: %s' % ' | '.join(card['clues']))
        print('  下一步: %s' % card['next_test'])
        print('  用 show %s 查看完整证据引用；经验不替代当前样本检查。' % card['lesson_id'])


def main(argv=None):
    args = parser().parse_args(argv)
    if not args.library:
        parser().error('须显式提供 --library 或设置 MITO_KNOWLEDGE_DIR')
    try:
        from _deps import require_stage
        require_stage('lesson_memory', __file__)
        store = LessonStore(args.library)
        if args.action == 'propose':
            revision = store.propose(args.case, args.record, args.authorize_local, args.allow_unresolved)
            print('经验候选已写入: %s rev%d [candidate]；尚未审核。' % (
                revision['lesson_id'], revision['revision']))
        elif args.action == 'review':
            revision = store.review(args.lesson_id, args.status, args.reviewer,
                                    args.reason, args.authorize_local)
            print('经验状态已更新: %s rev%d [%s]；审核记录不等于独立科学验证。' % (
                revision['lesson_id'], revision['revision'], revision['validation_status']))
        elif args.action == 'search':
            print_cards(store.search(args.query, args.status, args.domain,
                                     args.taxon, args.data_type, args.limit))
        elif args.action == 'show':
            revision = store.show(args.lesson_id, args.revision)
            print(json.dumps(revision, ensure_ascii=False, indent=2, allow_nan=False))
        elif args.action == 'supersede':
            revision = store.supersede(args.lesson_id, args.case, args.record,
                                       args.authorize_local, args.allow_unresolved)
            print('新修订已写入: %s rev%d [candidate]；旧修订保留为 superseded。' % (
                revision['lesson_id'], revision['revision']))
        elif args.action == 'rebuild':
            count = store.rebuild(args.authorize_local)
            print('索引已从 JSON 事实源重建: %d 条当前经验。' % count)
        return 0
    except MissingDependency as exc:
        print('环境错误: %s' % exc, file=sys.stderr)
        return 3
    except (OSError, ValueError, TypeError, sqlite3.Error) as exc:
        print('ERROR: %s' % exc, file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

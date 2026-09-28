#!/usr/bin/env python3
"""Inspect source-attributed taxon profiles without applying them to analyses."""
import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / 'config' / 'taxon_profiles.json'
SCHEMA = ROOT / 'schemas' / 'taxon-profiles.schema.json'


def load_profiles():
    try:
        import jsonschema
    except ImportError as exc:
        raise RuntimeError('缺少 jsonschema；请使用环境检查所列 Python 环境') from exc
    try:
        data = json.loads(REGISTRY.read_text(encoding='utf-8'))
        schema = json.loads(SCHEMA.read_text(encoding='utf-8'))
        jsonschema.validate(data, schema, format_checker=jsonschema.FormatChecker())
    except (OSError, ValueError, jsonschema.ValidationError) as exc:
        raise ValueError('类群档案或 schema 无效: %s' % exc) from exc
    names = [profile['scientific_name'].casefold() for profile in data['profiles']]
    if len(names) != len(set(names)):
        raise ValueError('类群档案含重复 scientific_name')
    return data['profiles']


def show_profile(profile):
    code = profile['mitochondrial_code']
    identifier = (' (NCBI TaxID %d)' % profile['ncbi_taxon_id']
                  if 'ncbi_taxon_id' in profile else ' (未登记种级 NCBI TaxID)')
    print('物种: %s%s' % (profile['scientific_name'], identifier))
    print('分类: %s' % ' > '.join(profile['lineage']))
    print('线粒体遗传密码表建议: table %d [%s]' % (code['table'], code['status']))
    print('物种特异证实: %s' % ('是' if code['species_specific_confirmed'] else '否'))
    print('自动选择: 否；分析命令仍须显式传入 --table')
    print('证据:')
    for item in code['evidence']:
        print('  - %s (%s)' % (item['claim'], item['url']))
        print('    适用范围: %s' % item['scope'])
    print('限制:')
    for item in code['limitations']:
        print('  - %s' % item)
    if profile['engineering_thresholds']:
        print('工程阈值:')
        for item in profile['engineering_thresholds']:
            print('  - %s=%s [%s]; %s' % (item['name'], item['value'], item['status'],
                                          item['source']['url']))
    else:
        print('工程阈值: 无类群特异配置')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest='command', required=True)
    subparsers.add_parser('list', help='列出已登记类群')
    show = subparsers.add_parser('show', help='显示一个类群档案')
    show.add_argument('scientific_name', help='完整学名')
    args = parser.parse_args(argv)

    try:
        profiles = load_profiles()
    except (RuntimeError, ValueError) as exc:
        print('ERROR: %s' % exc, file=sys.stderr)
        return 3 if isinstance(exc, RuntimeError) else 1

    if args.command == 'list':
        for profile in sorted(profiles, key=lambda item: item['scientific_name'].casefold()):
            code = profile['mitochondrial_code']
            taxid = 'NCBI TaxID=%d' % profile['ncbi_taxon_id'] if 'ncbi_taxon_id' in profile else 'NCBI TaxID=未登记'
            print('%s\t%s\ttable %d [%s]' % (
                profile['scientific_name'], taxid, code['table'], code['status']))
        return 0

    query = args.scientific_name.strip().casefold()
    matches = [profile for profile in profiles
               if profile['scientific_name'].casefold() == query]
    if not matches:
        print('未找到类群档案: %s（不根据近缘类群自动推断）' % args.scientific_name,
              file=sys.stderr)
        return 2
    show_profile(matches[0])
    return 0


if __name__ == '__main__':
    sys.exit(main())

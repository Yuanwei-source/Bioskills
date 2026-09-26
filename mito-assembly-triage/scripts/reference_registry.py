#!/usr/bin/env python3
"""公共参考序列的登记与获取（reference registry）。

为什么需要它：reads 与 FASTA 都有 SHA-256，但**诊断所比较的公共参考**从来没有被固定下来。
于是"和近缘物种比较"这种说法无法复现——参考库更新一次（`NC_060773.1` → `.2`），
同一个案例就可能得出不同结论，而证据链上没有任何痕迹。

本脚本把参考当成**分析输入**来管理，落实 `references/reference-policy.md` 的五条原则：

1. 公共参考可下载（不需要额外授权）；**样本序列不得上传**（那才需要显式许可）。
2. 每条参考必须登记 `source` / `accession`(**含版本**) / `hash` / `purposes`。
3. 参考版本必须固定：只接受 `ACCESSION.VERSION`；同一 accession 内容变化即报错（drift）。
4. 参考等级（L1 同种 … L5 远缘）**限制用途**：L5 只能做定位/筛查，不得判读 tRNA 丢失或重排。
5. 参考相似性只是线索，不是独立证据。

登记落在 `$MITO_KNOWLEDGE_DIR/references/registry.json`（**不写进 skill 安装目录**）。

用法：
    reference_registry.py register --file ref.gb --accession NC_060773.1 --source NCBI \
        --purposes gene_order_comparison --level L3 [--registry DIR]
    reference_registry.py acquire --accession NC_060773.1 --purposes boundary_validation \
        [--out DIR] [--dry-run]
    reference_registry.py record-database --name 'MITOS2 refseq invertebrate' \
        --path /refs/invertebrate --version refseq89m --purposes taxon_screen
    reference_registry.py list | verify
"""
import argparse
import datetime
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DATA_ROOT = os.environ.get('XDG_DATA_HOME') or os.path.join(os.path.expanduser('~'), '.local', 'share')
KNOWLEDGE = os.environ.get('MITO_KNOWLEDGE_DIR') or os.path.join(DEFAULT_DATA_ROOT, 'mito-assembly-triage')

REGISTRY_FORMAT = 'mito-reference-registry-1'
REGISTRY_FILE = 'registry.json'
DEFAULT_REGISTRY_DIR = os.path.join(KNOWLEDGE, 'references')

REFERENCE_KINDS = ('sequence', 'database')
REFERENCE_PURPOSES = ('gene_order_comparison', 'annotation_comparison', 'boundary_validation',
                      'identity_check', 'taxon_screen', 'contamination_screen', 'structure_check')
REFERENCE_LEVELS = ('L1', 'L2', 'L3', 'L4', 'L5')
# 等级限制用途：远缘参考只能定位/筛查，不得用于判读边界、顺序或结构
LEVEL_PURPOSES = {
    'L1': set(REFERENCE_PURPOSES),
    'L2': set(REFERENCE_PURPOSES),
    'L3': {'gene_order_comparison', 'annotation_comparison', 'identity_check',
           'taxon_screen', 'contamination_screen'},
    'L4': {'identity_check', 'taxon_screen', 'contamination_screen'},
    'L5': {'taxon_screen'},
}
# accession 必须带版本号：NC_060773.1 而不是 NC_060773
ACCESSION_VERSION = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*\.\d+$')
HEX64 = re.compile(r'^[0-9a-f]{64}$')


class AcquisitionError(RuntimeError):
    """获取/校验参考失败（网络、格式或内容不一致）。"""


def registry_path(directory):
    return pathlib.Path(directory) / REGISTRY_FILE


def load_registry(directory):
    path = registry_path(directory)
    if not path.exists():
        return {'format': REGISTRY_FORMAT, 'references': []}
    try:
        data = json.loads(path.read_text(encoding='utf-8'))
    except (ValueError, OSError) as exc:
        raise ValueError('registry 无法解析: %s (%s)' % (path, exc)) from exc
    if not isinstance(data, dict) or data.get('format') != REGISTRY_FORMAT:
        raise ValueError('registry 格式不兼容（期望 %s）: %s' % (REGISTRY_FORMAT, path))
    if not isinstance(data.get('references'), list):
        raise ValueError('registry.references 必须是数组: %s' % path)
    return data


def save_registry(directory, data):
    directory = pathlib.Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    tmp = registry_path(directory).with_suffix('.json.tmp')
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    os.replace(tmp, registry_path(directory))


def next_id(data):
    used = []
    for record in data.get('references', []):
        match = re.match(r'^ref-(\d+)$', str(record.get('id', '')))
        if match:
            used.append(int(match.group(1)))
    return 'ref-%03d' % ((max(used) + 1) if used else 1)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def directory_sha256(path):
    """Cheap, order-stable fingerprint of a database directory (name/size/mtime)."""
    path = pathlib.Path(path)
    entries = []
    if path.is_file():
        return file_sha256(path)
    for item in sorted(path.rglob('*')):
        if item.is_file():
            stat = item.stat()
            entries.append('%s\t%d\t%d' % (item.relative_to(path), stat.st_size, int(stat.st_mtime)))
    return hashlib.sha256('\n'.join(entries).encode('utf-8')).hexdigest()


def extract_sequence(text):
    """Sequence letters from a GenBank or FASTA payload (uppercased, no whitespace)."""
    if 'ORIGIN' in text and '//' in text:
        body = text.split('ORIGIN', 1)[1].split('//', 1)[0]
    elif text.lstrip().startswith('>'):
        body = '\n'.join(line for line in text.splitlines() if not line.startswith('>'))
    else:
        body = text
    return re.sub(r'[^A-Za-z]', '', body).upper()


def declared_length(text):
    match = re.search(r'^LOCUS\s+\S+\s+(\d+)\s+bp', text, re.M)
    return int(match.group(1)) if match else None


def declared_version(text):
    match = re.search(r'^VERSION\s+(\S+)', text, re.M)
    return match.group(1) if match else None


def validate_record(record):
    """Fail-closed validation of one registry entry. Returns a list of errors."""
    errors = []
    kind = record.get('kind')
    if kind not in REFERENCE_KINDS:
        return ['kind 取值非法: %s (允许: %s)' % (kind, '/'.join(REFERENCE_KINDS))]
    if not str(record.get('source') or '').strip():
        errors.append('必须记录 source（参考来自哪里，否则无法回溯）')
    purposes = record.get('purposes')
    if not isinstance(purposes, list) or not purposes:
        errors.append('必须声明 purposes（该参考被允许用于什么判断）')
        purposes = []
    for purpose in purposes:
        if purpose not in REFERENCE_PURPOSES:
            errors.append('purposes 取值非法: %s (允许: %s)' % (purpose, '/'.join(REFERENCE_PURPOSES)))
    if not str(record.get('retrieved') or '').strip():
        errors.append('必须记录 retrieved（获取日期）')
    if kind == 'sequence':
        accession = str(record.get('accession') or '')
        if not accession:
            errors.append('sequence 记录必须有 accession')
        elif not ACCESSION_VERSION.match(accession):
            errors.append('accession 必须固定版本（形如 NC_060773.1），当前: %s '
                          '—— 不带版本号无法复现' % accession)
        if not HEX64.match(str(record.get('file_sha256') or '')):
            errors.append('sequence 记录必须有 file_sha256（64 位十六进制）')
        if not HEX64.match(str(record.get('sequence_sha256') or '')):
            errors.append('sequence 记录必须有 sequence_sha256')
        if not isinstance(record.get('length_bp'), int) or record.get('length_bp', 0) <= 0:
            errors.append('sequence 记录必须有正的 length_bp')
        level = str(record.get('level') or '')
        if level not in REFERENCE_LEVELS:
            errors.append('sequence 记录必须有 level（%s）' % '/'.join(REFERENCE_LEVELS))
        else:
            allowed = LEVEL_PURPOSES[level]
            rejected = sorted(set(purposes) - allowed)
            if rejected:
                errors.append('level=%s（%s）不得用于 %s；该等级只允许 %s'
                              % (level, LEVEL_LEVEL_NAMES[level], '/'.join(rejected),
                                 '/'.join(sorted(allowed))))
    else:  # database
        if not str(record.get('name') or '').strip():
            errors.append('database 记录必须有 name')
        if not str(record.get('path') or '').strip():
            errors.append('database 记录必须有 path')
        if not str(record.get('version') or '').strip():
            errors.append('database 记录必须有 version（库版本/构建时间戳）'
                          '—— 否则今天与半年后的 top hit 可能不同')
        if not HEX64.match(str(record.get('file_sha256') or '')):
            errors.append('database 记录必须有 file_sha256（目录指纹或库文件 hash）')
    return errors


LEVEL_LEVEL_NAMES = {'L1': '同种', 'L2': '同属', 'L3': '同科', 'L4': '同目', 'L5': '远缘'}


def register_record(directory, record, force=False):
    """Validate, deduplicate/idempotence-check and persist one record."""
    errors = validate_record(record)
    if errors:
        raise ValueError('参考记录非法: %s' % '; '.join(errors))
    data = load_registry(directory)
    for existing in data.get('references', []):
        if existing.get('accession') and existing.get('accession') == record.get('accession'):
            if existing.get('file_sha256') == record.get('file_sha256'):
                return existing  # 同一文件，幂等
            if not force:
                raise ValueError(
                    '同一 accession(%s) 的内容 hash 与已登记记录不同：%s != %s —— '
                    '公共库可能已更新该记录（版本漂移）；不得静默替换，'
                    '请确认后带 --force 重新登记或改用新版本 accession'
                    % (record.get('accession'), record.get('file_sha256')[:12],
                       str(existing.get('file_sha256'))[:12]))
    record.setdefault('id', next_id(data))
    record.setdefault('registered', datetime.date.today().isoformat())
    data.setdefault('references', []).append(record)
    save_registry(directory, data)
    return record


def register_file(directory, path, accession, source, purposes, organism=None, taxon=None,
                  level=None, retrieved=None, version=None, force=False):
    path = pathlib.Path(path).resolve()
    if not path.is_file():
        raise ValueError('参考文件不存在: %s' % path)
    text = path.read_text(encoding='utf-8', errors='replace')
    if '//' not in text and not text.lstrip().startswith('>'):
        raise AcquisitionError('不是 GenBank/FASTA 格式（无法作为参考使用）: %s' % path)
    sequence = extract_sequence(text)
    declared = declared_length(text)
    if declared is not None and declared != len(sequence):
        raise AcquisitionError('序列长度(%d) 与 LOCUS 声明(%d) 不一致，文件可能被截断: %s'
                               % (len(sequence), declared, path))
    record = {'kind': 'sequence', 'accession': accession, 'organism': organism, 'taxon': taxon,
              'source': source, 'retrieved': retrieved or datetime.date.today().isoformat(),
              'version': version, 'level': level, 'purposes': list(purposes),
              'path': str(path), 'file_sha256': file_sha256(path),
              'sequence_sha256': hashlib.sha256(sequence.encode('ascii')).hexdigest(),
              'length_bp': declared or len(sequence)}
    return register_record(directory, record, force=force)


def record_database(directory, name, path, version, purposes, source=None, taxon_filter=None,
                    retrieved=None, force=False):
    path = pathlib.Path(path)
    if not path.exists():
        raise ValueError('数据库路径不存在: %s' % path)
    record = {'kind': 'database', 'name': name, 'path': str(path.resolve()),
              'version': version, 'source': source or 'local', 'taxon_filter': taxon_filter,
              'retrieved': retrieved or datetime.date.today().isoformat(),
              'purposes': list(purposes), 'file_sha256': directory_sha256(path)}
    return register_record(directory, record, force=force)


def plan_acquisition(accession, purposes, outdir, source='NCBI Nucleotide',
                     level=None, fetcher='efetch', database='nucleotide'):
    """Return the exact command + intended record WITHOUT downloading anything."""
    record = {'kind': 'sequence', 'accession': accession, 'source': source,
              'retrieved': datetime.date.today().isoformat(), 'level': level,
              'purposes': list(purposes)}
    errors = [e for e in validate_record(dict(record, file_sha256='0' * 64,
                                             sequence_sha256='0' * 64, length_bp=1))
              if 'sha256' not in e and 'length_bp' not in e]
    if errors:
        raise ValueError('无法规划获取: %s' % '; '.join(errors))
    outdir = pathlib.Path(outdir).resolve()
    target = outdir / (accession + '.gb')
    command = '%s -db %s -id %s -format gbwithparts > %s' % (fetcher, database, accession, target)
    return {'accession': accession, 'command': command, 'target': str(target),
            'record': record, 'outdir': str(outdir)}


def acquire(directory, outdir, accession, purposes, source='NCBI Nucleotide', level=None,
            database='nucleotide', force=False):
    plan = plan_acquisition(accession, purposes, outdir, source=source, level=level,
                            database=database)
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fetcher = 'efetch'
    if not any(os.access(os.path.join(entry, fetcher), os.X_OK)
               for entry in os.environ.get('PATH', '').split(os.pathsep) if entry):
        raise AcquisitionError('找不到 %s（NCBI edirect）；公共参考获取需要它，'
                               '或先用 register 登记手工获取的文件' % fetcher)
    completed = subprocess.run([fetcher, '-db', database, '-id', accession,
                                '-format', 'gbwithparts'],
                               capture_output=True, timeout=300)
    payload = completed.stdout.decode('utf-8', 'replace')
    if completed.returncode != 0 or not payload.strip():
        raise AcquisitionError('%s 获取失败（exit=%d）: %s'
                               % (accession, completed.returncode,
                                  completed.stderr.decode('utf-8', 'replace')[:200]))
    if '//' not in payload:
        raise AcquisitionError('下载内容不是完整 GenBank 记录（缺少 "//" 结束符），'
                               '拒绝登记以免把截断文件当参考')
    actual = declared_version(payload)
    if actual and actual != accession:
        raise AcquisitionError('下载内容的 VERSION(%s) 与请求的 accession(%s) 不一致，拒绝登记'
                               % (actual, accession))
    target = pathlib.Path(plan['target'])
    target.write_text(payload, encoding='utf-8')
    return register_file(directory, target, accession=accession, source=source,
                         purposes=purposes, level=level, retrieved=plan['record']['retrieved'],
                         version=database, force=force)


def verify_registry(directory):
    """Re-check every registered file: missing, or content changed."""
    problems = []
    data = load_registry(directory)
    for record in data.get('references', []):
        path = record.get('path')
        if not path:
            problems.append('%s: 记录没有 path' % record.get('id'))
            continue
        if not os.path.exists(path):
            problems.append('%s: 文件不存在/无法访问: %s' % (record.get('id'), path))
            continue
        current = directory_sha256(path) if record.get('kind') == 'database' else file_sha256(path)
        if record.get('file_sha256') and current != record.get('file_sha256'):
            problems.append('%s: 内容 hash 已变化（登记 %s…，当前 %s…）: %s'
                            % (record.get('id'), str(record['file_sha256'])[:12],
                               current[:12], path))
    return problems


def main(argv=None):
    parser = argparse.ArgumentParser(description='公共参考序列的登记与获取')
    sub = parser.add_subparsers(dest='command', required=True)

    def with_registry(target):
        target.add_argument('--registry', default=argparse.SUPPRESS,
                            help='registry 目录（默认 $MITO_KNOWLEDGE_DIR/references）')
        return target

    register = with_registry(sub.add_parser('register', help='登记已获取的参考文件'))
    register.add_argument('--file', required=True)
    register.add_argument('--accession', required=True)
    register.add_argument('--source', required=True)
    register.add_argument('--purposes', nargs='+', required=True)
    register.add_argument('--organism')
    register.add_argument('--taxon')
    register.add_argument('--level', choices=list(REFERENCE_LEVELS))
    register.add_argument('--version')
    register.add_argument('--retrieved')
    register.add_argument('--force', action='store_true')

    acquisition = with_registry(sub.add_parser('acquire', help='下载公共参考并登记（默认先看 --dry-run）'))
    acquisition.add_argument('--accession', required=True)
    acquisition.add_argument('--purposes', nargs='+', required=True)
    acquisition.add_argument('--source', default='NCBI Nucleotide')
    acquisition.add_argument('--level', choices=list(REFERENCE_LEVELS))
    acquisition.add_argument('--database', default='nucleotide')
    acquisition.add_argument('--out', default=None)
    acquisition.add_argument('--dry-run', action='store_true')
    acquisition.add_argument('--force', action='store_true')

    database = with_registry(sub.add_parser('record-database', help='登记本地数据库版本（BLAST/refseq 库）'))
    database.add_argument('--name', required=True)
    database.add_argument('--path', required=True)
    database.add_argument('--version', required=True)
    database.add_argument('--purposes', nargs='+', required=True)
    database.add_argument('--source')
    database.add_argument('--taxon-filter')
    database.add_argument('--force', action='store_true')

    with_registry(sub.add_parser('list', help='列出已登记参考'))
    with_registry(sub.add_parser('verify', help='校验已登记文件是否缺失或变化'))

    args = parser.parse_args(argv)
    if not hasattr(args, 'registry'):
        args.registry = DEFAULT_REGISTRY_DIR
    try:
        if args.command == 'register':
            record = register_file(args.registry, args.file, accession=args.accession,
                                   source=args.source, purposes=args.purposes,
                                   organism=args.organism, taxon=args.taxon, level=args.level,
                                   retrieved=args.retrieved, version=args.version, force=args.force)
            print('已登记 %s: %s (%s, %d bp)' % (record['id'], record['accession'],
                                                record['file_sha256'][:12], record['length_bp']))
        elif args.command == 'acquire':
            outdir = args.out or args.registry
            if args.dry_run:
                plan = plan_acquisition(args.accession, args.purposes, outdir,
                                        source=args.source, level=args.level,
                                        database=args.database)
                print('dry-run（未下载）: %s' % plan['command'])
                print('将登记到: %s' % registry_path(args.registry))
            else:
                outdir = args.out or args.registry
                record = acquire(args.registry, outdir, args.accession, args.purposes,
                                 source=args.source, level=args.level, database=args.database,
                                 force=args.force)
                print('已获取并登记 %s: %s (%s, %d bp) -> %s'
                      % (record['id'], record['accession'], record['file_sha256'][:12],
                         record['length_bp'], record['path']))
        elif args.command == 'record-database':
            record = record_database(args.registry, args.name, args.path, args.version,
                                     args.purposes, source=args.source,
                                     taxon_filter=args.taxon_filter, force=args.force)
            print('已登记数据库 %s: %s (%s, %s)' % (record['id'], record['name'],
                                                  record['version'], record['file_sha256'][:12]))
        elif args.command == 'list':
            data = load_registry(args.registry)
            if not data['references']:
                print('registry 为空: %s' % registry_path(args.registry))
            for record in data['references']:
                print('%s  %-14s %-9s %s' % (record.get('id'),
                                             record.get('accession') or record.get('name'),
                                             record.get('level') or record.get('kind'),
                                             ', '.join(record.get('purposes') or [])))
        elif args.command == 'verify':
            problems = verify_registry(args.registry)
            for problem in problems:
                print('⚠ %s' % problem)
            if problems:
                return 1
            print('全部参考文件与登记一致')
    except (ValueError, AcquisitionError, OSError, subprocess.SubprocessError) as exc:
        print('✗ %s' % exc, file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

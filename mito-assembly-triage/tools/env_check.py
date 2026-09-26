#!/usr/bin/env python3
"""环境就绪检查：一次部署，多次使用（并保证"环境变了会被发现"）。

三种模式，对应两种使用场景：

* `--setup`（第一次 / 换机 / 升级后）：按 `config/dependencies.json` 全量盘点，
  按 essential / extended / optional 分组报告，缺失项给出**用途 + 安装方式**；
  essential 齐全时写入 environment lock（默认 `$MITO_KNOWLEDGE_DIR/environment.lock.json`）。
* `--daily`（日常）：读 lock 做**轻量**检查——只探测"本次需要的那一组"是否还在。
  不再让用户面对三十个工具，也不会因为无关工具缺失而拦住任何事。
* `--stage NAME`（单步门禁）：只检查某个 capability 的依赖，供脚本前置调用（缺则退出 3）。

**不变式（重要）**：lock 是**缓存**，不是信任凭证。日常模式仍然真实探测所需工具；
lock 只用于避免重复全量盘点、并提供"上次验证版本"用于漂移对比。若只信 lock，
就退化成"第一次通过、以后永远相信"——正是本 skill 反复清除的静默降级。

退出码：`0` 就绪 / `2` essential 缺失（与 `scripts/check_env.sh` 一致）/ `3` 指定 stage 的依赖缺失 /
`1` 用法或清单错误。
"""
import argparse
import datetime
import json
import os
import pathlib
import re
import shutil
import socket
import subprocess
import sys

SKILL_DIR = pathlib.Path(__file__).resolve().parent.parent
MANIFEST = SKILL_DIR / 'config' / 'dependencies.json'
DEFAULT_DATA_ROOT = os.environ.get('XDG_DATA_HOME') or os.path.join(os.path.expanduser('~'), '.local', 'share')
KNOWLEDGE = os.environ.get('MITO_KNOWLEDGE_DIR') or os.path.join(DEFAULT_DATA_ROOT, 'mito-assembly-triage')
DEFAULT_LOCK = os.path.join(KNOWLEDGE, 'environment.lock.json')
LOCK_FORMAT = 'mito-environment-lock-1'
TIERS = ('essential', 'extended', 'optional')
VERSION_RE = re.compile(r'\d+(?:\.\d+)+')
EXIT_OK, EXIT_USAGE, EXIT_ESSENTIAL, EXIT_STAGE = 0, 1, 2, 3


def load_manifest(path=None):
    path = pathlib.Path(path or MANIFEST)
    try:
        data = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, ValueError) as exc:
        raise ValueError('依赖清单无法解析: %s (%s)' % (path, exc)) from exc
    if data.get('format') != 'mito-dependency-manifest-1':
        raise ValueError('依赖清单格式不兼容: %s' % path)
    return data


def manifest_entries(manifest):
    """[(id, tier, probe, purpose, install)]，python 库与外部工具并列。"""
    entries = []
    for key in ('python_libs', 'tools'):
        for item in manifest.get(key) or []:
            entries.append((item['id'], item['tier'], item.get('probe') or {},
                            item.get('purpose', ''), item.get('install', '')))
    return entries


def _first_version(text):
    match = VERSION_RE.search(text or '')
    return match.group(0) if match else None


def _resolve_in_env_value(name, value):
    """env 变量可能指向：可执行文件、目录、或 PATH 列表。返回命中的可执行路径或 None。"""
    if not value:
        return None
    if os.path.isfile(value) and os.access(value, os.X_OK):
        return value
    for entry in value.split(os.pathsep):
        if not entry:
            continue
        candidate = os.path.join(entry, name)
        if os.access(candidate, os.X_OK):
            return candidate
    return None


def probe_command(name, search_dir=None, env_value=None):
    """(可执行路径, 版本) 或 (None, None)。

    `search_dir`（CLI 的 `--path`）表示**只在该目录内查找**——用于"检查某个 env 的 bin 是否齐"，
    因此不再回退 PATH；否则回退 PATH。
    """
    found = _resolve_in_env_value(name, env_value)
    if found is None and search_dir:
        candidate = os.path.join(search_dir, name)
        found = candidate if os.access(candidate, os.X_OK) else None
    if found is None and not search_dir:
        found = shutil.which(name)
    if found is None:
        return None, None
    version = None
    try:
        proc = subprocess.run([found, '--version'], capture_output=True, text=True, timeout=20)
        version = _first_version((proc.stdout or '') + (proc.stderr or ''))
    except (OSError, subprocess.SubprocessError):
        pass
    return found, version


def probe_python_module(module, interpreter=None, search_dir=None, distribution=None):
    """(解释器, 版本) 或 (None, None)。

    可用性只取决于**导入是否成功**；版本号按发行版名查，查不到就不报版本
    （导入名与发行版名并不总是相同，例如模块 `Bio` 的发行版名是 `biopython`）。
    """
    interpreter = interpreter or sys.executable
    if search_dir:  # 允许测试/多环境场景指定解释器目录
        candidate = os.path.join(search_dir, os.path.basename(interpreter))
        if os.access(candidate, os.X_OK):
            interpreter = candidate
    code = ("import importlib\n"
            "importlib.import_module(%r)\n" % module +
            "try:\n"
            "    import importlib.metadata as _m\n"
            "    print(_m.version(%r))\n" % (distribution or module) +
            "except Exception:\n"
            "    print('installed')\n")
    try:
        proc = subprocess.run([interpreter, '-c', code], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None, None
    if proc.returncode != 0:
        return None, None
    out = (proc.stdout or '').strip()
    version = _first_version(out)
    return interpreter, version


def probe_network(host, timeout=4):
    try:
        with socket.create_connection((host, 443), timeout=timeout):
            return host, None
    except OSError:
        return None, None


def probe(entry, search_dir=None, allow_network=False):
    """按 probe 规格探测一个依赖，返回 (available: bool, detail: str|None, version: str|None)。"""
    kind = (entry.get('probe') or {}).get('type')
    env_var = (entry.get('probe') or {}).get('env')
    env_path = os.environ.get(env_var) if env_var else None
    if kind == 'python-interpreter':
        # "本 skill 用的 Python"就是运行本工具的解释器；--path 指定 env 时优先用该目录下的同名解释器
        interpreter = sys.executable
        if search_dir:
            candidate = os.path.join(search_dir, os.path.basename(interpreter))
            if os.access(candidate, os.X_OK):
                interpreter = candidate
        try:
            proc = subprocess.run([interpreter, '--version'], capture_output=True, text=True, timeout=20)
            version = _first_version((proc.stdout or '') + (proc.stderr or ''))
        except (OSError, subprocess.SubprocessError):
            return False, interpreter, None
        return True, interpreter, version
    if kind == 'command':
        path, version = probe_command(entry['probe']['name'], search_dir, env_path)
        return bool(path), path, version
    if kind == 'python-module':
        found, version = probe_python_module(entry['probe']['module'], env_path, search_dir,
                                            entry['probe'].get('distribution'))
        return bool(found), found, version
    if kind == 'python-import':
        found, version = probe_python_module(entry['probe']['module'], None, search_dir,
                                            entry['probe'].get('distribution'))
        return bool(found), found, version
    if kind == 'env-dir':
        base = os.environ.get(entry['probe']['env']) or ''
        sub = os.environ.get(entry['probe'].get('sub') or '') or ''
        target = os.path.join(base, sub) if sub else base
        return bool(base and os.path.isdir(target)), target or '(未设置)', None
    if kind == 'network':
        if not allow_network:
            return None, '未探测（需要 --network）', None
        host, _ = probe_network(entry['probe']['host'])
        return bool(host), host, None
    return False, '未知 probe 类型', None


def inventory(manifest, search_dir=None, allow_network=False):
    report = []
    for dep_id, tier, spec, purpose, install in manifest_entries(manifest):
        available, detail, version = probe({'probe': spec}, search_dir, allow_network)
        report.append({'id': dep_id, 'tier': tier, 'available': available, 'detail': detail,
                       'version': version, 'purpose': purpose, 'install': install})
    return report


def group_by_tier(report):
    grouped = {tier: [] for tier in TIERS}
    for item in report:
        grouped.setdefault(item['tier'], []).append(item)
    return grouped


def stage_requirements(manifest, stage):
    stages = manifest.get('stages') or {}
    if stage not in stages:
        raise KeyError(stage)
    return stages[stage]


def dependencies_for(manifest, ids):
    wanted = set(ids)
    return [item for item in manifest_entries(manifest) if item[0] in wanted]


def write_lock(path, manifest, report, mode):
    payload = {
        'format': LOCK_FORMAT,
        'checked_at': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
        'mode': mode,
        'skill_dir': str(SKILL_DIR),
        'env': {key: os.environ.get(key) for key in
                ('CONDA_ROOT', 'MITOS2_PY', 'MITOS2_REFDIR', 'MITOS2_REFSEQVER',
                 'MITOS2_EXTRA_PATH', 'MINIMAP2', 'PLOT_PY')},
        'available': {item['id']: {'version': item['version'], 'detail': item['detail']}
                      for item in report if item['available']},
        'missing': {item['id']: {'tier': item['tier'], 'purpose': item['purpose'],
                                 'install': item['install']}
                    for item in report if not item['available']},
    }
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    os.replace(tmp, path)
    return payload


def read_lock(path):
    path = pathlib.Path(path)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return None
    return data if data.get('format') == LOCK_FORMAT else None


def print_tier_report(grouped):
    labels = {'essential': '核心（essential）', 'extended': '增强（extended）', 'optional': '按需（optional）'}
    for tier in TIERS:
        print('\n=== %s ===' % labels[tier])
        for item in grouped.get(tier, []):
            mark = '✓' if item['available'] else ('-' if item['available'] is None else '✗')
            version = ' %s' % item['version'] if item['version'] else ''
            print('  %s %s%s' % (mark, item['id'], version))
            if not item['available']:
                print('      %s' % item['purpose'])
                print('      安装: %s' % item['install'])


def main(argv=None):
    parser = argparse.ArgumentParser(description='环境就绪检查（部署一次 / 日常轻量）')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--setup', action='store_true', help='全量盘点；essential 齐全时写 lock')
    mode.add_argument('--daily', action='store_true', help='读 lock 做轻量检查（默认模式）')
    mode.add_argument('--stage', metavar='NAME', help='只检查某个 capability 的依赖（缺则退出 3）')
    parser.add_argument('--json', action='store_true')
    parser.add_argument('--strict', action='store_true', help='extended 也必须齐全才算就绪')
    parser.add_argument('--network', action='store_true', help='额外探测联网类依赖（较慢）')
    parser.add_argument('--path', default=None, help='只在指定目录内查找可执行文件（检查某个 env 的 bin）')
    parser.add_argument('--lock', default=DEFAULT_LOCK, help='lock 路径（默认 $MITO_KNOWLEDGE_DIR/environment.lock.json）')
    parser.add_argument('--manifest', default=str(MANIFEST))
    args = parser.parse_args(argv)

    try:
        manifest = load_manifest(args.manifest)
    except ValueError as exc:
        print('✗ %s' % exc, file=sys.stderr)
        return EXIT_USAGE

    if args.stage:
        try:
            spec = stage_requirements(manifest, args.stage)
        except KeyError:
            print('✗ 未知 stage: %s（可选：%s）'
                  % (args.stage, ', '.join(sorted(manifest.get('stages') or {}))), file=sys.stderr)
            return EXIT_USAGE
        entries = dependencies_for(manifest, spec.get('requires') or [])
        report = [dict(item=dep_id, tier=tier, purpose=purpose, install=install,
                       **dict(zip(('available', 'detail', 'version'),
                                  probe({'probe': spec2}, args.path, args.network))))
                  for dep_id, tier, spec2, purpose, install in entries]
        missing = [item for item in report if not item['available']]
        if args.json:
            print(json.dumps({'stage': args.stage, 'capability': spec.get('capability'),
                              'ready': not missing, 'dependencies': report}, ensure_ascii=False, indent=1))
        else:
            print('stage: %s（%s）' % (args.stage, spec.get('capability')))
            for item in report:
                mark = '✓' if item['available'] else '-'
                print('  %s %s' % (mark, item['item']))
            if missing:
                print('  该步骤缺少依赖，无法执行：')
                for item in missing:
                    print('    ✗ %s —— %s\n      安装: %s' % (item['item'], item['purpose'], item['install']))
                print('  → 本步骤未执行；其余不依赖这些工具的步骤不受影响。')
        return EXIT_OK if not missing else EXIT_STAGE

    report = inventory(manifest, args.path, args.network)
    grouped = group_by_tier(report)
    essential_missing = [i for i in grouped['essential'] if not i['available']]
    extended_missing = [i for i in grouped['extended'] if not i['available']]
    ready = not essential_missing and not (args.strict and extended_missing)

    if args.setup:
        if args.json:
            print(json.dumps({'mode': 'setup', 'ready': ready, 'dependencies': report},
                             ensure_ascii=False, indent=1))
        else:
            print('Mito skill 环境盘点（--setup）')
            print_tier_report(grouped)
            if essential_missing:
                print('\n核心环境: 未就绪 —— 缺 %d 项（上面已给用途与安装方式）' % len(essential_missing))
                print('补齐后重跑: bash scripts/check_env.sh --setup')
            else:
                print('\n核心环境: READY')
                if extended_missing:
                    print('增强环境: 缺 %d 项 —— 基础诊断可用，全能力需补齐' % len(extended_missing))
                if ready:
                    lock = write_lock(args.lock, manifest, report, 'setup')
                    print('已记录环境 lock: %s（checked_at=%s）' % (args.lock, lock['checked_at']))
        return EXIT_OK if ready else EXIT_ESSENTIAL

    # 日常模式：lock 是缓存，仍需真实探测
    lock = read_lock(args.lock)
    if lock is None:
        report = inventory(manifest, args.path, args.network)
        essential_missing = [i for i in report if i['tier'] == 'essential' and not i['available']]
        if args.json:
            print(json.dumps({'mode': 'daily', 'lock': None, 'ready': not essential_missing,
                              'hint': '尚无 lock，建议运行 --setup 记录环境',
                              'dependencies': report}, ensure_ascii=False, indent=1))
        else:
            print('环境: 尚无 lock（%s）' % args.lock)
            if essential_missing:
                print('核心环境未就绪，缺:')
                for item in essential_missing:
                    print('  ✗ %s —— %s\n    安装: %s' % (item['id'], item['purpose'], item['install']))
            else:
                print('核心环境就绪；建议运行一次 --setup 记录 lock（之后检查更快）。')
        return EXIT_OK if not essential_missing else EXIT_ESSENTIAL

    report = inventory(manifest, args.path, args.network)
    by_id = {item['id']: item for item in report}
    relevant = [i for i in by_id.values() if i['tier'] == 'essential'] + \
               ([i for i in by_id.values() if i['tier'] == 'extended'] if args.strict else [])
    missing = [i for i in relevant if not i['available']]
    drifted = []
    for item in relevant:
        if not item['available']:
            continue
        before = (lock.get('available') or {}).get(item['id']) or {}
        if before.get('version') and item['version'] and before['version'] != item['version']:
            drifted.append((item['id'], before['version'], item['version']))
    if args.json:
        print(json.dumps({'mode': 'daily', 'lock': args.lock, 'ready': not missing,
                          'missing': [i['id'] for i in missing],
                          'version_drift': drifted}, ensure_ascii=False, indent=1))
    else:
        print('环境状态（lock: %s，上次验证 %s）' % (args.lock, lock.get('checked_at')))
        if missing:
            print('✗ 环境已变化 —— 以下工具本次需要但探测不到：')
            for item in missing:
                print('    ✗ %s —— %s\n      安装: %s' % (item['id'], item['purpose'], item['install']))
            print('  → 相关步骤停止；补齐后重跑 --setup 更新 lock。')
        else:
            print('✓ 就绪（本次需要 %d 项，均可用）' % len(relevant))
        if drifted:
            print('⚠ 版本与 lock 不同（可能是有意升级；若影响复现请记录）:')
            for dep_id, old, new in drifted:
                print('    %s: %s → %s' % (dep_id, old, new))
    return EXIT_OK if not missing else EXIT_ESSENTIAL


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""
进化引擎: 管理 skill 的知识层 (knowledge/)

命令:
  search --query "关键词"           # 检索历史案例 (按信号/物种/内容)
  add-case --sample X --species Y   # 追加案例 (写入 MITO_KNOWLEDGE_DIR/cases/X.md)
  update-signals --signal S --judgment J --ref 案例   # 追加信号
  update-pitfalls --pitfall P --fix F                # 追加坑
  suggest-promotions                # 模式提炼: 统计重复信号/动作, 建议固化
  stats                             # 显示知识库统计
"""
import sys, os, re, argparse, datetime, hashlib, json, pathlib, shutil, tempfile, urllib.request

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DATA_ROOT = os.environ.get('XDG_DATA_HOME') or os.path.join(os.path.expanduser('~'), '.local', 'share')
KNOWLEDGE = os.environ.get('MITO_KNOWLEDGE_DIR') or os.path.join(DEFAULT_DATA_ROOT, 'mito-assembly-triage')
CASES = os.path.join(KNOWLEDGE, 'cases')
SIGNALS = os.path.join(KNOWLEDGE, 'signals.md')
PITFALLS = os.path.join(KNOWLEDGE, 'pitfalls.md')
STATS = os.path.join(KNOWLEDGE, 'stats.md')
LESSON_CANDIDATES = os.path.join(KNOWLEDGE, 'lessons', 'candidates')
LESSON_VERIFIED = os.path.join(KNOWLEDGE, 'lessons', 'verified')
PUBLIC_KNOWLEDGE = os.path.join(KNOWLEDGE, 'public')


def case_files():
    if not os.path.isdir(CASES):
        return []
    return sorted(fn for fn in os.listdir(CASES) if fn.endswith('.md'))

def structured_case_files():
    if not os.path.isdir(CASES): return []
    return sorted(str(p) for p in pathlib.Path(CASES).glob('**/case.json'))

def lesson_dirs():
    os.makedirs(LESSON_CANDIDATES, exist_ok=True)
    os.makedirs(LESSON_VERIFIED, exist_ok=True)

def _json(path):
    with open(path, encoding='utf-8') as fh: return json.load(fh)

def _write_json(path, value):
    path = pathlib.Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    with open(tmp, 'w', encoding='utf-8') as fh: json.dump(value, fh, ensure_ascii=False, indent=2); fh.write('\n')
    os.replace(tmp, path)

def _tokens(value):
    return set(re.findall(r'[\w-]{2,}', json.dumps(value, ensure_ascii=False).lower()))

def structured_search(query):
    terms = set(query.lower().split()); results=[]
    for path in structured_case_files():
        try: case = _json(path)
        except (OSError, ValueError): continue
        words = _tokens(case)
        score = len(terms & words)
        if score: results.append((score, case.get('decision', {}).get('status', 'UNKNOWN'), path, case))
    for path in pathlib.Path(LESSON_CANDIDATES).glob('*.json') if os.path.isdir(LESSON_CANDIDATES) else []:
        try: lesson = _json(path)
        except (OSError, ValueError): continue
        score = len(terms & _tokens(lesson))
        if score: results.append((score, lesson.get('validation_status', 'candidate'), str(path), lesson))
    for score, status, path, _ in sorted(results, key=lambda x: (-x[0], x[2])):
        print('%d\t%s\t%s' % (score, status, path))
    if not results: print('未找到结构化经验: %s' % query)
    return results

def detect_lesson_conflicts(lesson, directory=None):
    directory = directory or LESSON_CANDIDATES
    conflicts=[]
    for path in pathlib.Path(directory).glob('*.json') if os.path.isdir(directory) else []:
        try: other = _json(path)
        except (OSError, ValueError): continue
        same_scope = _tokens(lesson.get('applicable_when', [])) & _tokens(other.get('applicable_when', []))
        same_clue = _tokens(lesson.get('diagnostic_clues', [])) & _tokens(other.get('diagnostic_clues', []))
        if same_scope and same_clue and lesson.get('suggested_next_test') != other.get('suggested_next_test'):
            conflicts.append(str(path))
    return conflicts

def propose_lesson(args):
    lesson_dirs(); case = _json(pathlib.Path(args.case) / 'case.json')
    if case.get('decision', {}).get('status') == 'UNRESOLVED' and not args.allow_unresolved:
        raise ValueError('未解决案例不能直接提炼；使用 --allow-unresolved 仅生成候选并保留限制')
    lesson = {'lesson_id': args.lesson_id or case.get('case_id'), 'applicable_when': [case.get('issue', {}).get('type')],
              'not_applicable_when': [], 'diagnostic_clues': [case.get('issue', {}).get('user_observation', '')],
              'suggested_next_test': args.next_test, 'supporting_case_ids': [case.get('case_id')],
              'counterexample_case_ids': [], 'sources': [{'case_id': case.get('case_id'), 'path': str(args.case)}],
              'validation_status': 'candidate', 'version': '1.0.0', 'last_reviewed': None}
    conflicts = detect_lesson_conflicts(lesson)
    if conflicts: lesson['conflicts'] = conflicts
    path = os.path.join(LESSON_CANDIDATES, lesson['lesson_id'] + '.json')
    _write_json(path, lesson); print('候选经验写入: %s' % path)
    if conflicts: print('冲突候选: %s' % ', '.join(conflicts))

def export_contribution(args):
    if not args.authorize: raise ValueError('必须显式指定 --authorize；默认不导出/上传')
    case = _json(pathlib.Path(args.case) / 'case.json')
    payload = json.dumps(case, ensure_ascii=False)
    sensitive = re.compile(r'(fastq|bam|token|password|secret|api[_-]?key|/home/|/root/)', re.I)
    if sensitive.search(payload): raise ValueError('检测到敏感路径或凭据字段，拒绝生成贡献文件')
    contribution = {'format': 'mito-experience-contribution-1', 'generated_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'case': case, 'review_status': 'candidate'}
    _write_json(args.output, contribution); print('贡献文件已生成（未上传）: %s' % args.output)

def _read_source(source):
    if re.match(r'^https?://', source):
        with urllib.request.urlopen(source, timeout=30) as response: return response.read()
    with open(source, 'rb') as fh: return fh.read()

def sync_public(args):
    raw = _read_source(args.manifest)
    manifest = json.loads(raw.decode('utf-8'))
    if manifest.get('format') != 'mito-public-knowledge-1': raise ValueError('manifest 格式不兼容')
    staging = pathlib.Path(tempfile.mkdtemp(prefix='mito-public-', dir=KNOWLEDGE))
    destination = pathlib.Path(PUBLIC_KNOWLEDGE)
    backup = pathlib.Path(str(destination) + '.previous')
    moved_old = False
    try:
        for item in manifest.get('items', []):
            if item.get('status') in ('revoked', 'withdrawn'): continue
            rel = pathlib.PurePosixPath(item.get('path', ''))
            if not rel.parts or rel.is_absolute() or '..' in rel.parts: raise ValueError('非法公共知识路径')
            data = _read_source(item.get('url') or os.path.join(os.path.dirname(args.manifest), str(rel)))
            if hashlib.sha256(data).hexdigest() != item.get('sha256'): raise ValueError('内容校验失败: %s' % rel)
            target = staging.joinpath(*rel.parts); target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes(data)
        manifest_path = staging / 'manifest.json'; manifest_path.write_bytes(raw)
        if destination.exists():
            if backup.exists(): shutil.rmtree(backup)
            os.replace(destination, backup)
            moved_old = True
        os.replace(staging, destination)
        print('公共知识同步完成（本地案例未覆盖，远程内容不执行）: %s' % destination)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True); raise
    finally:
        if moved_old and not destination.exists() and backup.exists():
            os.replace(backup, destination)


def ensure_cases_dir():
    os.makedirs(CASES, exist_ok=True)


def ensure_knowledge_files():
    ensure_cases_dir()
    defaults = [
        (SIGNALS, '# 信号映射表\n\n| 信号 | 判定 | 处理 | 来源 |\n|---|---|---|---|\n'),
        (PITFALLS, '# 常见工具坑\n\n| 坑 | 现象 | 解决 | 来源 |\n|---|---|---|---|\n'),
        (STATS, '# 样品处理统计\n\n暂无样品记录。\n'),
    ]
    for path, content in defaults:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if not os.path.exists(path):
            with open(path, 'x', encoding='utf-8') as fh:
                fh.write(content)

def search(query):
    print('=== 检索案例 (关键词: %s) ===' % query)
    kw = [k.lower() for k in query.split()]
    found = 0
    for fn in case_files():
        content = open(os.path.join(CASES, fn)).read()
        if all(k in content.lower() for k in kw):
            found += 1
            print('\n▶ %s' % fn)
            # 打印 frontmatter 摘要
            m = re.search(r'^---\n(.*?)\n---', content, re.S)
            if m:
                for line in m.group(1).split('\n'):
                    if ':' in line and not line.startswith(' '):
                        print('   %s' % line.strip())
            # 打印信号/判定/新知识小节
            for sec in ['信号', '判定', '新知识']:
                mm = re.search(r'## %s.*?\n(.*?)(?=\n## |\Z)' % sec, content, re.S)
                if mm:
                    lines = [l.strip() for l in mm.group(1).split('\n') if l.strip()][:3]
                    print('   [%s] %s' % (sec, ' | '.join(lines)))
    if not found:
        print('未找到匹配案例 — 这是新的问题类型! 处理完请 add-case 沉淀')
    return found

def case_path(sample):
    if not sample or not isinstance(sample, str) or '/' in sample or '\\' in sample or '\x00' in sample or os.path.basename(sample) != sample or sample in ('.', '..'):
        raise ValueError('sample must be a single file name')
    path = os.path.join(CASES, '%s.md' % sample)
    if os.path.islink(path):
        raise ValueError('case path must not be a symlink')
    return path

def add_case(args):
    ensure_knowledge_files()
    fn = case_path(args.sample)
    force = getattr(args, 'force', False)
    if os.path.exists(fn) and not force:
        raise FileExistsError(fn)
    content = """---
sample: %(sample)s
species: %(species)s
date: %(date)s
tools: %(tools)s
status: 待审核
signals: [%(signals)s]
---

# %(sample)s 线粒体组装诊断案例

## 信号 (Signal)
%(signals)s

## 判定 (Diagnosis)
%(diagnosis)s

## 动作 (Actions)
%(actions)s

## 新知识 (Lessons)
%(lessons)s
""" % {
        'sample': args.sample, 'species': args.species or '待鉴定',
        'date': datetime.date.today().isoformat(), 'tools': args.tools or '待记录',
        'signals': args.signals or '待补充', 'diagnosis': args.diagnosis or '待补充',
        'actions': args.actions or '待补充', 'lessons': args.lessons or '待补充',
    }
    with open(fn, 'w' if force else 'x') as fh:
        fh.write(content)
    print('✓ 案例写入: %s (状态=待审核；审核后再改为已确认)' % fn)
    print('  记得同时更新 signals.md / pitfalls.md (有新知识时) 和 stats.md')


def update_signals(args):
    ensure_knowledge_files()
    if not os.path.exists(SIGNALS):
        print('✗ 找不到 signals.md'); return
    line = '| %s | %s | %s | %s |' % (args.signal, args.judgment, '待处理', args.ref or '手动')
    with open(SIGNALS, 'a') as fh:
        fh.write(line + '\n')
    print('✓ 已追加信号: %s' % args.signal)

def update_pitfalls(args):
    ensure_knowledge_files()
    if not os.path.exists(PITFALLS):
        print('✗ 找不到 pitfalls.md'); return
    with open(PITFALLS, 'a') as fh:
        fh.write('| %s | %s | %s | %s |\n' % (args.pitfall, '见案例', args.fix or '待解决', datetime.date.today().isoformat()))
    print('✓ 已追加坑: %s' % args.pitfall)

def suggest_promotions():
    print('=== 模式提炼 (统计 signals.md 出现频次) ===')
    # 统计案例数
    n_cases = len(case_files())
    print('案例数: %d' % n_cases)
    # 从案例的 frontmatter 统计 signals 出现频次
    sig_count = {}
    for fn in case_files():
        content = open(os.path.join(CASES, fn)).read()
        m = re.search(r'^signals: \[(.*?)\]', content, re.M)
        if m:
            for s in m.group(1).split(','):
                s = s.strip()
                sig_count[s] = sig_count.get(s, 0) + 1
    if sig_count:
        print('\n信号频次:')
        for s, c in sorted(sig_count.items(), key=lambda x: -x[1]):
            mark = '  ★ 出现%d次 → 仅审核通过案例才可建议提升到 SKILL.md!' % c if c >= 3 else ''
            print('  %s: %d%s' % (s, c, mark))
    # 脚本固化建议
    print('\n建议: 当某手工操作重复 3 次以上, 固化为 scripts/ 下的脚本;')
    print('当某信号重复 3 次以上, 提升为 SKILL.md 的固定检查项。')

def stats():
    files = case_files()
    print('知识库位置: %s' % KNOWLEDGE)
    print('案例数: %d' % len(files))
    for f in files:
        c = open(os.path.join(CASES, f)).read()
        m = re.search(r'^status: (.*)', c, re.M)
        print('  - %s [%s]' % (f.replace('.md',''), (m.group(1) if m else '?')))
    print('\n文件:')
    for d in [CASES, SIGNALS, PITFALLS, STATS]:
        print('  %s' % d)

# Structured case helpers are deliberately integrated here rather than exposed
# as a second diagnostic engine. Legacy Markdown cases remain readable.
def v2_case_root(raw):
    root = pathlib.Path(raw).resolve()
    root.mkdir(parents=True, exist_ok=True)
    return root

def file_sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def case_init(args):
    root = v2_case_root(args.directory)
    inputs = []
    for role, raw in args.input:
        path = pathlib.Path(raw).resolve()
        if not path.is_file(): raise ValueError('输入不存在: %s' % raw)
        inputs.append({'role': role, 'path': str(path), 'sha256': file_sha256(path)})
    case = {'schema_version': '2.0', 'case_id': args.case_id or root.name,
            'taxon': {'name': args.taxon, 'taxid': None, 'genetic_code': None},
            'inputs': inputs, 'issue': {'type': args.issue, 'user_observation': args.observation},
            'hypotheses': [], 'events_file': 'events.jsonl',
            'decision': {'status': 'UNRESOLVED', 'confidence': 'not_assessable', 'rationale': '尚未完成足够检查'},
            'modifications': [], 'validation': [], 'lessons_proposed': []}
    tmp = root / 'case.json.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh: json.dump(case, fh, ensure_ascii=False, indent=2); fh.write('\n')
    os.replace(tmp, root / 'case.json')
    (root / 'events.jsonl').touch()
    print('case initialized: %s' % (root / 'case.json'))

def case_event(args):
    root = pathlib.Path(args.directory).resolve()
    with open(root / 'case.json', encoding='utf-8') as fh: case = json.load(fh)
    record = {'timestamp': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'action': args.action,
              'command': args.command, 'tool_versions': args.tool_version, 'result': args.result,
              'motivation': args.motivation, 'impact': args.impact}
    with open(root / case.get('events_file', 'events.jsonl'), 'a', encoding='utf-8') as fh:
        fh.write(json.dumps(record, ensure_ascii=False) + '\n')
    print('event recorded')

def case_validate(args):
    root = pathlib.Path(args.directory).resolve()
    with open(root / 'case.json', encoding='utf-8') as fh: case = json.load(fh)
    errors = []
    if case.get('schema_version') != '2.0': errors.append('schema_version must be 2.0')
    if case.get('decision', {}).get('status') not in {'RESOLVED', 'NO_CHANGE', 'UNRESOLVED'}: errors.append('invalid decision.status')
    if case.get('decision', {}).get('confidence') not in {'high', 'moderate', 'low', 'not_assessable'}: errors.append('invalid decision.confidence')
    if not (root / case.get('events_file', 'events.jsonl')).is_file(): errors.append('events file missing')
    if errors:
        print('INVALID'); [print('- ' + e) for e in errors]; return 1
    print('VALID'); return 0

def case_report(args):
    root = pathlib.Path(args.directory).resolve()
    with open(root / 'case.json', encoding='utf-8') as fh: case = json.load(fh)
    lines = ['# 诊断案例 %s' % case['case_id'], '', '## 事件', '']
    with open(root / case.get('events_file', 'events.jsonl'), encoding='utf-8') as fh:
        for line in fh:
            e = json.loads(line)
            lines.append('- `%s`：%s；影响：%s' % (e.get('action'), e.get('result'), e.get('impact') or '未记录'))
    d = case['decision']
    lines += ['', '## 当前判定', '', '- 状态：`%s`' % d['status'], '- 置信等级：`%s`' % d['confidence'], '- 理由：%s' % d['rationale'], '', '## 证据边界', '', '缺少 reads 时不得写成 raw-read-supported。']
    (root / 'case.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('report written')

def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(0)
    cmd = sys.argv[1]
    if cmd == 'search':
        if '--query' not in sys.argv:
            print('需要 --query'); sys.exit(1)
        search(sys.argv[sys.argv.index('--query')+1])
    elif cmd == 'add-case':
        ap = argparse.ArgumentParser()
        for a in ['sample','species','tools','signals','diagnosis','actions','lessons']:
            ap.add_argument('--%s' % a)
        ap.add_argument('--force', action='store_true')
        try:
            add_case(ap.parse_args(sys.argv[2:]))
        except (ValueError, FileExistsError) as exc:
            print('✗ %s' % exc, file=sys.stderr)
            sys.exit(1)
    elif cmd == 'update-signals':
        ap = argparse.ArgumentParser()
        ap.add_argument('--signal', required=True); ap.add_argument('--judgment', required=True)
        ap.add_argument('--ref')
        update_signals(ap.parse_args(sys.argv[2:]))
    elif cmd == 'update-pitfalls':
        ap = argparse.ArgumentParser()
        ap.add_argument('--pitfall', required=True); ap.add_argument('--fix')
        update_pitfalls(ap.parse_args(sys.argv[2:]))
    elif cmd == 'suggest-promotions':
        suggest_promotions()
    elif cmd == 'stats':
        stats()
    elif cmd == 'search-structured':
        ap = argparse.ArgumentParser(); ap.add_argument('--query', required=True); structured_search(ap.parse_args(sys.argv[2:]).query)
    elif cmd == 'propose-lesson':
        ap = argparse.ArgumentParser(); ap.add_argument('--case', required=True); ap.add_argument('--next-test', required=True); ap.add_argument('--lesson-id'); ap.add_argument('--allow-unresolved', action='store_true')
        try: propose_lesson(ap.parse_args(sys.argv[2:]))
        except (ValueError, OSError, json.JSONDecodeError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
    elif cmd == 'export-contribution':
        ap = argparse.ArgumentParser(); ap.add_argument('--case', required=True); ap.add_argument('--output', required=True); ap.add_argument('--authorize', action='store_true')
        try: export_contribution(ap.parse_args(sys.argv[2:]))
        except (ValueError, OSError, json.JSONDecodeError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
    elif cmd == 'sync-public':
        ap = argparse.ArgumentParser(); ap.add_argument('--manifest', required=True)
        try: sync_public(ap.parse_args(sys.argv[2:]))
        except (ValueError, OSError, json.JSONDecodeError, urllib.error.URLError) as exc: print('✗ 公共知识未同步: %s' % exc, file=sys.stderr); sys.exit(1)
    elif cmd in ('case-init', 'case-event', 'case-validate', 'case-report'):
        ap = argparse.ArgumentParser()
        if cmd == 'case-init':
            ap.add_argument('directory'); ap.add_argument('--case-id'); ap.add_argument('--issue', required=True); ap.add_argument('--observation', required=True); ap.add_argument('--taxon'); ap.add_argument('--input', nargs=2, action='append', default=[], metavar=('ROLE', 'PATH'))
            try: case_init(ap.parse_args(sys.argv[2:]))
            except (ValueError, OSError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
        elif cmd == 'case-event':
            ap.add_argument('directory'); ap.add_argument('--action', required=True); ap.add_argument('--result', required=True); ap.add_argument('--impact', default=''); ap.add_argument('--command', default=''); ap.add_argument('--tool-version', default=''); ap.add_argument('--motivation', default='')
            case_event(ap.parse_args(sys.argv[2:]))
        elif cmd == 'case-validate':
            ap.add_argument('directory'); sys.exit(case_validate(ap.parse_args(sys.argv[2:])))
        else:
            ap.add_argument('directory'); case_report(ap.parse_args(sys.argv[2:]))
    else:
        print('未知命令: %s' % cmd); print(__doc__)

if __name__ == '__main__':
    main()

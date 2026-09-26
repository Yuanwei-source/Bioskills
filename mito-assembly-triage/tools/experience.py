#!/usr/bin/env python3
"""
进化引擎: 管理 skill 的知识层 (knowledge/)

命令:
  search --query "关键词"           # 检索历史案例 (按信号/物种/内容)
  add-case --sample X --species Y   # 追加案例 (写入 MITO_KNOWLEDGE_DIR/cases/X.md)
  case-anomaly <dir> --id A1 --claim '...' --status UNRESOLVED --confidence low \
      [--reads-support NOT_ASSESSED] [--event-action <ACTION> ...] [--update]
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
MAX_PUBLIC_ITEM_BYTES = 2 * 1024 * 1024
LESSON_STATUSES = {'candidate', 'verified', 'rejected', 'deprecated', 'withdrawn', 'superseded'}
# candidate 以外的状态都是终态/审核结论; 只有 verified 可进入公共同步
LESSON_REVIEW_STATUSES = LESSON_STATUSES - {'candidate'}
LESSON_RETIRED_STATUSES = {'rejected', 'deprecated', 'withdrawn', 'superseded', 'revoked'}
CASE_SCHEMA = os.path.join(SKILL_DIR, 'schemas', 'case.schema.json')
SAFE_ID = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$')


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

def safe_identifier(value, label='identifier'):
    if not isinstance(value, str) or not SAFE_ID.fullmatch(value) or value in ('.', '..'):
        raise ValueError('%s 只能包含字母、数字、点、下划线和短横线' % label)
    return value

def safe_case_member(root, value, label='case member'):
    if not isinstance(value, str) or pathlib.PurePath(value).name != value or '\\' in value or value in ('.', '..') or '\x00' in value:
        raise ValueError('%s 必须是案例目录内的单一文件名' % label)
    root = pathlib.Path(root).resolve(); target = root / value
    if target.is_symlink(): raise ValueError('%s 不能是符号链接' % label)
    return target

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
    sources = ((LESSON_CANDIDATES, 'local-candidate'), (LESSON_VERIFIED, 'local-verified'),
               (PUBLIC_KNOWLEDGE, 'public'))
    for directory, source in sources:
        for path in pathlib.Path(directory).glob('**/*.json') if os.path.isdir(directory) else []:
            if path.name == 'manifest.json': continue
            try: lesson = _json(path)
            except (OSError, ValueError): continue
            if lesson.get('validation_status') in LESSON_RETIRED_STATUSES: continue
            score = len(terms & _tokens(lesson))
            if score:
                item = dict(lesson); item['knowledge_source'] = source
                results.append((score, lesson.get('validation_status', 'candidate'), str(path), item))
    for score, status, path, _ in sorted(results, key=lambda x: (-x[0], x[2])):
        print('%d\t%s\t%s' % (score, status, path))
    if not results: print('未找到结构化经验: %s' % query)
    return results

def detect_lesson_conflicts(lesson, directory=None):
    conflicts=[]
    directories = [directory] if directory else [LESSON_CANDIDATES, LESSON_VERIFIED]
    paths = []
    for current in directories:
        if current and os.path.isdir(current): paths.extend(pathlib.Path(current).glob('*.json'))
    for path in paths:
        try: other = _json(path)
        except (OSError, ValueError): continue
        if other.get('lesson_id') == lesson.get('lesson_id'): continue
        same_scope = _tokens(lesson.get('applicable_when', [])) & _tokens(other.get('applicable_when', []))
        same_clue = _tokens(lesson.get('diagnostic_clues', [])) & _tokens(other.get('diagnostic_clues', []))
        decisions = {lesson.get('decision_status'), other.get('decision_status')} - {None, 'UNRESOLVED'}
        opposite_scope = (_tokens(lesson.get('applicable_when', [])) & _tokens(other.get('not_applicable_when', []))) or (_tokens(other.get('applicable_when', [])) & _tokens(lesson.get('not_applicable_when', [])))
        if same_scope and same_clue and (len(decisions) > 1 or opposite_scope):
            conflicts.append(str(path))
    return conflicts

def propose_lesson(args):
    lesson_dirs(); case = _json(pathlib.Path(args.case) / 'case.json')
    if case.get('decision', {}).get('status') == 'UNRESOLVED' and not args.allow_unresolved:
        raise ValueError('未解决案例不能直接提炼；使用 --allow-unresolved 仅生成候选并保留限制')
    lesson_id = safe_identifier(args.lesson_id or case.get('case_id'), 'lesson_id')
    lesson = {'lesson_id': lesson_id, 'applicable_when': [case.get('issue', {}).get('type')],
              'not_applicable_when': [], 'diagnostic_clues': [case.get('issue', {}).get('user_observation', '')],
              'suggested_next_test': args.next_test, 'supporting_case_ids': [case.get('case_id')],
              'counterexample_case_ids': [], 'sources': [{'case_id': case.get('case_id'), 'path': str(args.case)}],
              'decision_status': case.get('decision', {}).get('status'),
              'validation_status': 'candidate', 'version': '1.0.0', 'last_reviewed': None}
    conflicts = detect_lesson_conflicts(lesson)
    if conflicts: lesson['conflicts'] = conflicts
    path = os.path.join(LESSON_CANDIDATES, lesson['lesson_id'] + '.json')
    _write_json(path, lesson); print('候选经验写入: %s' % path)
    if conflicts: print('冲突候选: %s' % ', '.join(conflicts))

def review_lesson(args):
    safe_identifier(args.lesson_id, 'lesson_id')
    if args.status not in LESSON_REVIEW_STATUSES: raise ValueError('非法经验状态')
    source = pathlib.Path(LESSON_CANDIDATES) / (args.lesson_id + '.json')
    if not source.exists(): source = pathlib.Path(LESSON_VERIFIED) / (args.lesson_id + '.json')
    if not source.exists(): raise ValueError('找不到 lesson: %s' % args.lesson_id)
    lesson = _json(source); previous = lesson.get('validation_status', 'candidate')
    lesson['validation_status'] = args.status; lesson['last_reviewed'] = datetime.date.today().isoformat()
    lesson.setdefault('review_history', []).append({'from': previous, 'to': args.status, 'reviewer': args.reviewer, 'reason': args.reason, 'date': lesson['last_reviewed']})
    target = pathlib.Path(LESSON_VERIFIED) / source.name if args.status == 'verified' else pathlib.Path(LESSON_CANDIDATES) / source.name
    _write_json(target, lesson)
    if target != source: source.unlink()
    print('经验状态已更新: %s -> %s' % (previous, args.status))

def export_contribution(args):
    if not args.authorize: raise ValueError('必须显式指定 --authorize；默认不导出/上传')
    case = _json(pathlib.Path(args.case) / 'case.json')
    safe_case = {'schema_version': case.get('schema_version'),
                 'case_id': hashlib.sha256(str(case.get('case_id', '')).encode()).hexdigest()[:16],
                 'taxon': case.get('taxon', {}),
                 'inputs': [{'role': item.get('role'), 'sha256': item.get('sha256')} for item in case.get('inputs', [])],
                 'issue': {'type': case.get('issue', {}).get('type')},
                 'hypotheses': case.get('hypotheses', []), 'decision': case.get('decision', {}),
                 'validation': case.get('validation', []), 'modifications': case.get('modifications', [])}
    safe_case = sanitize_for_share(safe_case)
    payload = json.dumps(safe_case, ensure_ascii=False)
    sensitive = re.compile(r'(token|password|secret|api[_-]?key|BEGIN (?:RSA|OPENSSH) PRIVATE KEY)', re.I)
    if sensitive.search(payload): raise ValueError('检测到凭据字段，拒绝生成贡献文件')
    contribution = {'format': 'mito-experience-contribution-1', 'generated_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'case': safe_case, 'review_status': 'candidate'}
    _write_json(args.output, contribution); print('贡献文件已生成（未上传）: %s' % args.output)

def sanitize_for_share(value):
    if isinstance(value, dict):
        blocked = re.compile(r'(?:^|_)(?:path|command|token|password|secret|api[_-]?key|account)(?:$|_)', re.I)
        return {key: sanitize_for_share(item) for key, item in value.items() if not blocked.search(str(key))}
    if isinstance(value, list): return [sanitize_for_share(item) for item in value]
    if isinstance(value, str):
        value = re.sub(r'(?<!\w)(?:/[\w. -]+){2,}', '[local-path-redacted]', value)
        value = re.sub(r'\b[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}\b', '[email-redacted]', value)
    return value

def _read_source(source, max_bytes=MAX_PUBLIC_ITEM_BYTES):
    if re.match(r'^https?://', source):
        with urllib.request.urlopen(source, timeout=30) as response: data = response.read(max_bytes + 1)
    else:
        with open(source, 'rb') as fh: data = fh.read(max_bytes + 1)
    if len(data) > max_bytes: raise ValueError('公共知识文件超过大小限制')
    return data

def validate_public_lesson(data, expected_path):
    try: lesson = json.loads(data.decode('utf-8'))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc: raise ValueError('公共知识不是有效 JSON: %s' % expected_path) from exc
    required = {'lesson_id', 'applicable_when', 'not_applicable_when', 'diagnostic_clues', 'suggested_next_test', 'supporting_case_ids', 'counterexample_case_ids', 'sources', 'validation_status', 'version', 'last_reviewed'}
    missing = sorted(required - set(lesson))
    if missing: raise ValueError('公共知识缺少字段 %s: %s' % (','.join(missing), expected_path))
    safe_identifier(lesson['lesson_id'], 'lesson_id')
    if lesson['validation_status'] not in LESSON_STATUSES: raise ValueError('公共知识状态不兼容: %s' % expected_path)
    if lesson['validation_status'] != 'verified': raise ValueError('公共知识必须为 verified: %s' % expected_path)
    return lesson

def sync_public(args):
    os.makedirs(KNOWLEDGE, exist_ok=True)
    raw = _read_source(args.manifest)
    manifest = json.loads(raw.decode('utf-8'))
    if manifest.get('format') != 'mito-public-knowledge-1': raise ValueError('manifest 格式不兼容')
    staging = pathlib.Path(tempfile.mkdtemp(prefix='mito-public-', dir=KNOWLEDGE))
    destination = pathlib.Path(PUBLIC_KNOWLEDGE)
    backup = pathlib.Path(str(destination) + '.previous')
    moved_old = False
    try:
        seen = set()
        for item in manifest.get('items', []):
            if item.get('status') in LESSON_RETIRED_STATUSES: continue
            rel = pathlib.PurePosixPath(item.get('path', ''))
            if not rel.parts or rel.is_absolute() or '..' in rel.parts: raise ValueError('非法公共知识路径')
            if str(rel) in seen: raise ValueError('manifest 包含重复路径: %s' % rel)
            seen.add(str(rel))
            data = _read_source(item.get('url') or os.path.join(os.path.dirname(args.manifest), str(rel)))
            if hashlib.sha256(data).hexdigest() != item.get('sha256'): raise ValueError('内容校验失败: %s' % rel)
            validate_public_lesson(data, rel)
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
    requested = pathlib.Path(raw)
    if requested.is_symlink(): raise ValueError('案例目录不能是符号链接')
    root = requested.resolve()
    root.mkdir(parents=True, exist_ok=True)
    return root

def file_sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def _hypothesis_statements(raw):
    """[(explicit_id|None, explanation)] from --hypothesis values.

    A caller who writes `--hypothesis "H1: boundary error"` has already supplied the
    id; prepending our own would produce `id=H1` with an explanation starting with
    `H1:` (the F-05 duplication).  A leading `H<n>:` / `H<n>、` prefix is therefore
    stripped from the explanation, and when EVERY statement carries an explicit id
    those ids are honoured as given.
    """
    statements = [str(text).strip() for text in (raw or []) if str(text).strip()]
    parsed, explicit = [], []
    for text in statements:
        match = re.match(r'^\s*H(\d+)\s*[:：.、)\]]\s*(.*)$', text)
        if match:
            explicit.append('H%s' % match.group(1))
            parsed.append(match.group(2).strip() or text)
        else:
            explicit.append(None)
            parsed.append(text)
    if all(explicit):
        if len(set(explicit)) != len(explicit):
            raise ValueError('--hypothesis 的显式编号重复: %s' % ','.join(sorted(explicit)))
        return list(zip(explicit, parsed))
    return [(None, text) for text in parsed]


def case_init(args):
    root = v2_case_root(args.directory)
    inputs = []
    for role, raw in args.input:
        path = pathlib.Path(raw).resolve()
        if not path.is_file(): raise ValueError('输入不存在: %s' % raw)
        inputs.append({'role': role, 'path': str(path), 'sha256': file_sha256(path)})
    pairs = _hypothesis_statements(getattr(args, 'hypothesis', None))
    if not pairs:
        raise ValueError('case-init 至少要给出 1 条候选解释: --hypothesis "..." '
                         '(SKILL.md 的 HYPOTHESIZE 步骤建议先列出多条竞争解释)')
    hypotheses = []
    for index, (explicit, text) in enumerate(pairs):
        hypotheses.append({'id': explicit or 'H%d' % (index + 1), 'explanation': text,
                           'support': [], 'against': [], 'unknown': []})
    case = {'schema_version': '2.0', 'case_id': args.case_id or root.name,
            'taxon': {'name': args.taxon, 'taxid': None, 'genetic_code': None},
            'inputs': inputs, 'issue': {'type': args.issue, 'user_observation': args.observation},
            'hypotheses': hypotheses, 'events_file': 'events.jsonl',
            'decision': {'status': 'UNRESOLVED', 'confidence': 'not_assessable', 'rationale': '尚未完成足够检查'},
            'modifications': [], 'validation': [], 'lessons_proposed': []}
    tmp = root / 'case.json.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh: json.dump(case, fh, ensure_ascii=False, indent=2); fh.write('\n')
    os.replace(tmp, root / 'case.json')
    events_path = safe_case_member(root, 'events.jsonl', 'events_file')
    events_path.touch()
    print('case initialized: %s' % (root / 'case.json'))

def case_event(args):
    root = pathlib.Path(args.directory).resolve()
    with open(root / 'case.json', encoding='utf-8') as fh: case = json.load(fh)
    record = {'timestamp': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'action': args.action,
              'command': args.command, 'tool_versions': args.tool_version, 'result': args.result,
              'motivation': args.motivation, 'impact': args.impact}
    events_path = safe_case_member(root, case.get('events_file', 'events.jsonl'), 'events_file')
    with open(events_path, 'a', encoding='utf-8') as fh:
        fh.write(json.dumps(record, ensure_ascii=False) + '\n')
    print('event recorded')

DECISION_STATUSES = ('RESOLVED', 'NO_CHANGE', 'UNRESOLVED')
CONFIDENCE_LEVELS = ('high', 'moderate', 'low', 'not_assessable')
READS_SUPPORT_LEVELS = ('NOT_ASSESSED', 'READS_CONSISTENT', 'READS_DISCRIMINATING')


def _case_errors_without_jsonschema(case):
    """Explicit equivalent of schemas/case.schema.json.

    This is what runs when ``jsonschema`` is unavailable.  CI installs jsonschema,
    so the primary path runs there and the two cross-comparison tests are skipped;
    this fallback is therefore covered by the tests that call it DIRECTLY
    (``test_fallback_verdicts_are_fixed``, ``SchemaKeywordCoverageTests`` and the
    Schema-keyword x JSON-basic-type matrix).  It **must not be weaker than the
    schema**: a weaker fallback would silently skip the ``anomalies[]`` enum checks
    while the docs claim they are validated, and it must never raise on a legal
    JSON value.  That is exactly the failure mode this project keeps hitting.
    """
    if not isinstance(case, dict):
        return ['case.json 必须是 JSON 对象']
    errors = []
    required = ('schema_version', 'case_id', 'inputs', 'issue', 'hypotheses', 'decision', 'events_file')
    for key in required:
        if key not in case:
            errors.append('缺少必需字段: %s' % key)
    if case.get('schema_version') != '2.0':
        errors.append('schema_version 必须是 "2.0"')
    if 'case_id' in case and (not isinstance(case['case_id'], str) or not case['case_id']):
        errors.append('case_id 必须是非空字符串 (minLength: 1)')
    if 'events_file' in case and not isinstance(case['events_file'], str):
        errors.append('events_file 必须是字符串')
    for key in ('inputs', 'hypotheses', 'anomalies', 'modifications',
                'validation', 'lessons_proposed'):
        if key in case and not isinstance(case[key], list):
            errors.append('%s 必须是数组' % key)
    hypotheses = case.get('hypotheses')
    if isinstance(hypotheses, list):
        if not hypotheses:
            errors.append('hypotheses 不能为空 (至少 1 条候选解释)')
        for index, item in enumerate(hypotheses):
            if not isinstance(item, dict) or \
                    not {'id', 'explanation', 'support', 'against', 'unknown'} <= set(item):
                errors.append('hypotheses[%d] 缺少必需字段 '
                              '(id/explanation/support/against/unknown)' % index)
    inputs = case.get('inputs')
    if isinstance(inputs, list):
        for index, item in enumerate(inputs):
            if not isinstance(item, dict) or not {'role', 'path', 'sha256'} <= set(item):
                errors.append('inputs[%d] 缺少必需字段 (role/path/sha256)' % index)
    decision = case.get('decision')
    if 'decision' in case:
        if not isinstance(decision, dict):
            # explicit null/str is NOT the same as an absent key
            errors.append('decision 必须是对象')
        else:
            if not {'status', 'confidence', 'rationale'} <= set(decision):
                errors.append('decision 缺少必需字段 (status/confidence/rationale)')
            if decision.get('status') not in DECISION_STATUSES:
                errors.append('decision.status 取值非法: %s' % decision.get('status'))
            if decision.get('confidence') not in CONFIDENCE_LEVELS:
                errors.append('decision.confidence 取值非法: %s' % decision.get('confidence'))
    if 'issue' in case:
        issue = case['issue']
        if not isinstance(issue, dict):
            errors.append('issue 必须是对象')
        elif 'type' not in issue:
            errors.append('issue 缺少必需字段 (type)')
    if 'taxon' in case and not isinstance(case['taxon'], dict):
        errors.append('taxon 必须是对象')
    anomalies = case.get('anomalies')
    if isinstance(anomalies, list):
        for index, item in enumerate(anomalies):
            if not isinstance(item, dict):
                errors.append('anomalies[%d] 必须是对象' % index)
                continue
            if not {'id', 'claim', 'status', 'confidence'} <= set(item):
                errors.append('anomalies[%d] 缺少必需字段 (id/claim/status/confidence)' % index)
            if item.get('status') not in DECISION_STATUSES:
                errors.append('anomalies[%d].status 取值非法: %s' % (index, item.get('status')))
            if item.get('confidence') not in CONFIDENCE_LEVELS:
                errors.append('anomalies[%d].confidence 取值非法: %s' % (index, item.get('confidence')))
            if 'reads_support' in item and item['reads_support'] not in READS_SUPPORT_LEVELS:
                errors.append('anomalies[%d].reads_support 取值非法: %s'
                              % (index, item['reads_support']))
    return errors


def case_anomaly(args):
    """Write one anomaly into case.json's anomalies[] -- the format the schema owns.

    Deliberately NOT a second record format: the entry goes into the same
    schema-validated case.json, and the resulting case is re-validated before the
    file is replaced.  `--event-action` must name an existing event action, so the
    anomaly's link to the original tool output is real rather than decorative.
    """
    root = pathlib.Path(args.directory).resolve()
    case_path = root / 'case.json'
    if not case_path.is_file():
        raise ValueError('案例不存在: %s' % case_path)
    with open(case_path, encoding='utf-8') as fh:
        case = json.load(fh)
    if not isinstance(case, dict):
        raise ValueError('case.json 必须是 JSON 对象')
    anomalies = case.setdefault('anomalies', [])
    if not isinstance(anomalies, list):
        raise ValueError('case.json 的 anomalies 必须是数组')
    anomaly_id = str(args.id).strip()
    claim = str(args.claim).strip()
    if not anomaly_id:
        raise ValueError('--id 不能为空')
    if not claim:
        raise ValueError('--claim 不能为空')
    # enums are read from the single set of constants the schema path also uses
    if args.status not in DECISION_STATUSES:
        raise ValueError('--status 取值非法: %s (允许: %s)'
                         % (args.status, '/'.join(DECISION_STATUSES)))
    if args.confidence not in CONFIDENCE_LEVELS:
        raise ValueError('--confidence 取值非法: %s (允许: %s)'
                         % (args.confidence, '/'.join(CONFIDENCE_LEVELS)))
    if args.reads_support and args.reads_support not in READS_SUPPORT_LEVELS:
        raise ValueError('--reads-support 取值非法: %s (允许: %s)'
                         % (args.reads_support, '/'.join(READS_SUPPORT_LEVELS)))

    actions = [str(a).strip() for a in (args.event_action or []) if str(a).strip()]
    if actions:
        events_path = safe_case_member(root, case.get('events_file', 'events.jsonl'),
                                       'events_file')
        known = set()
        if events_path.is_file():
            for line in events_path.read_text(encoding='utf-8').splitlines():
                if not line.strip():
                    continue
                try:
                    record = json.loads(line)
                except ValueError:
                    continue
                if isinstance(record, dict) and record.get('action'):
                    known.add(str(record['action']))
        unknown = [a for a in actions if a not in known]
        if unknown:
            raise ValueError('--event-action 引用了不存在的事件: %s '
                             '(关联必须真实, 不能凭空挂接)' % ', '.join(unknown))

    entry = {'id': anomaly_id, 'claim': claim, 'status': args.status,
             'confidence': args.confidence}
    if args.reads_support:
        entry['reads_support'] = args.reads_support
    if actions:
        entry['evidence_events'] = list(dict.fromkeys(actions))

    matches = [index for index, item in enumerate(anomalies)
               if isinstance(item, dict) and item.get('id') == anomaly_id]
    if matches and not args.update:
        raise ValueError('anomaly id 已存在: %s (如需替换请加 --update)' % anomaly_id)
    if matches:
        anomalies[matches[0]] = entry
    else:
        anomalies.append(entry)

    errors = case_schema_errors(case)
    if errors:
        raise ValueError('写回后案例不符合 schema, 已放弃写入: %s' % '; '.join(errors))

    tmp = root / 'case.json.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(case, fh, ensure_ascii=False, indent=2)
        fh.write('\n')
    os.replace(tmp, case_path)
    print('anomaly recorded: %s (%s/%s%s)'
          % (anomaly_id, entry['status'], entry['confidence'],
             ', reads_support=%s' % entry['reads_support'] if args.reads_support else ''))


def case_schema_errors(case):
    """Enforce schemas/case.schema.json; fall back to an equivalent explicit check."""
    schema_path = CASE_SCHEMA
    try:
        import jsonschema
    except ImportError:
        jsonschema = None
    if jsonschema is not None and os.path.exists(schema_path):
        try:
            jsonschema.validate(case, _json(schema_path))
        except jsonschema.ValidationError as exc:
            location = '/'.join(str(part) for part in exc.absolute_path) or '<root>'
            return ['schema 校验失败: %s (%s)' % (exc.message, location)]
        return []
    return _case_errors_without_jsonschema(case)


def case_validate(args):
    root = pathlib.Path(args.directory).resolve()
    with open(root / 'case.json', encoding='utf-8') as fh: case = json.load(fh)
    if not isinstance(case, dict):
        print('INVALID'); print('- case.json 必须是 JSON 对象'); return 1
    errors = case_schema_errors(case)
    if case.get('schema_version') != '2.0': errors.append('schema_version must be 2.0')
    # An explicit "decision": null is not the same as an absent key: the default in
    # case.get('decision', {}) does NOT apply, so the value must be type-checked
    # before any field is read.  These redundant business checks run after the
    # schema so the CLI still reports INVALID instead of raising AttributeError.
    decision = case.get('decision')
    if not isinstance(decision, dict):
        decision = {}
        errors.append('decision 必须是对象')
    if decision.get('status') not in {'RESOLVED', 'NO_CHANGE', 'UNRESOLVED'}: errors.append('invalid decision.status')
    if decision.get('confidence') not in {'high', 'moderate', 'low', 'not_assessable'}: errors.append('invalid decision.confidence')
    try: events_path = safe_case_member(root, case.get('events_file', 'events.jsonl'), 'events_file')
    except ValueError as exc: errors.append(str(exc)); events_path = None
    if events_path is not None and not events_path.is_file(): errors.append('events file missing')
    if errors:
        print('INVALID'); [print('- ' + e) for e in errors]; return 1
    print('VALID'); return 0

def case_report(args):
    root = pathlib.Path(args.directory).resolve()
    with open(root / 'case.json', encoding='utf-8') as fh: case = json.load(fh)
    lines = ['# 诊断案例 %s' % case['case_id'], '', '## 事件', '']
    with open(safe_case_member(root, case.get('events_file', 'events.jsonl'), 'events_file'), encoding='utf-8') as fh:
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
    elif cmd == 'review-lesson':
        ap = argparse.ArgumentParser(); ap.add_argument('--lesson-id', required=True); ap.add_argument('--status', required=True); ap.add_argument('--reviewer', required=True); ap.add_argument('--reason', required=True)
        try: review_lesson(ap.parse_args(sys.argv[2:]))
        except (ValueError, OSError, json.JSONDecodeError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
    elif cmd == 'sync-public':
        ap = argparse.ArgumentParser(); ap.add_argument('--manifest', required=True)
        try: sync_public(ap.parse_args(sys.argv[2:]))
        except (ValueError, OSError, json.JSONDecodeError, urllib.error.URLError) as exc: print('✗ 公共知识未同步: %s' % exc, file=sys.stderr); sys.exit(1)
    elif cmd in ('case-init', 'case-event', 'case-anomaly', 'case-validate', 'case-report'):
        ap = argparse.ArgumentParser()
        if cmd == 'case-init':
            ap.add_argument('directory'); ap.add_argument('--case-id'); ap.add_argument('--issue', required=True); ap.add_argument('--observation', required=True); ap.add_argument('--taxon'); ap.add_argument('--input', nargs=2, action='append', default=[], metavar=('ROLE', 'PATH')); ap.add_argument('--hypothesis', action='append', default=[], metavar='TEXT')
            try: case_init(ap.parse_args(sys.argv[2:]))
            except (ValueError, OSError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
        elif cmd == 'case-event':
            ap.add_argument('directory'); ap.add_argument('--action', required=True); ap.add_argument('--result', required=True); ap.add_argument('--impact', default=''); ap.add_argument('--command', default=''); ap.add_argument('--tool-version', default=''); ap.add_argument('--motivation', default='')
            case_event(ap.parse_args(sys.argv[2:]))
        elif cmd == 'case-anomaly':
            ap.add_argument('directory'); ap.add_argument('--id', required=True); ap.add_argument('--claim', required=True); ap.add_argument('--status', required=True); ap.add_argument('--confidence', required=True); ap.add_argument('--reads-support', dest='reads_support', default=None); ap.add_argument('--event-action', dest='event_action', action='append', default=[], metavar='ACTION'); ap.add_argument('--update', action='store_true')
            try: case_anomaly(ap.parse_args(sys.argv[2:]))
            except (ValueError, OSError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
        elif cmd == 'case-validate':
            ap.add_argument('directory'); sys.exit(case_validate(ap.parse_args(sys.argv[2:])))
        else:
            ap.add_argument('directory'); case_report(ap.parse_args(sys.argv[2:]))
    else:
        print('未知命令: %s' % cmd); print(__doc__)

if __name__ == '__main__':
    main()

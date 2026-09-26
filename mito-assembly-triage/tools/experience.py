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
# 公共参考登记表：与案例同属证据链，但落在知识目录（不进 skill 安装目录）
REFERENCE_REGISTRY_DIR = os.path.join(KNOWLEDGE, 'references')
REFERENCE_REGISTRY_FILE = 'registry.json'
MAX_PUBLIC_ITEM_BYTES = 2 * 1024 * 1024
LESSON_STATUSES = {'candidate', 'verified', 'rejected', 'deprecated', 'withdrawn', 'superseded'}
# candidate 以外的状态都是终态/审核结论; 只有 verified 可进入公共同步
LESSON_REVIEW_STATUSES = LESSON_STATUSES - {'candidate'}
LESSON_RETIRED_STATUSES = {'rejected', 'deprecated', 'withdrawn', 'superseded', 'revoked'}
LESSON_DOMAINS = ('tool', 'annotation', 'biology')
# `tool_failure_case` 只能是关于工具/环境的 lesson，不得升为生物学或样本质量结论。
CASE_TYPE_LESSON_DOMAINS = {'abnormal_case': ('biology', 'annotation', 'tool'),
                            'normal_validation_case': ('biology', 'annotation', 'tool'),
                            'tool_failure_case': ('tool',)}
CASE_SCHEMA = os.path.join(SKILL_DIR, 'schemas', 'case.schema.json')
# 案例类型：只记"异常"会让经验系统偏向"线粒体一定有问题"（特异性偏差），
# 因此正常验证案例与工具故障案例都必须可记录。
CASE_TYPES = ('abnormal_case', 'normal_validation_case', 'tool_failure_case')
# lesson 的推广上限：默认 fail-closed（single_case / none），防止"一次案例 → 规则"。
GENERALIZATION_SCOPES = ('single_case', 'species', 'genus', 'family', 'order', 'multi_taxon')
TRANSFERABILITY_LEVELS = ('none', 'low', 'moderate', 'high')
# 支持的**独立**案例数下限（独立 = case_id 不同且记录的类群不同；无法证明独立的按不独立处理）。
SCOPE_MIN_CASES = {'single_case': 1, 'species': 1, 'genus': 2, 'family': 2, 'order': 3, 'multi_taxon': 3}
TRANSFERABILITY_MIN_CASES = {'none': 0, 'low': 1, 'moderate': 3, 'high': 5}
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

def _case_taxon_name(case):
    """Normalised taxon name of a case, or None when it is not recorded."""
    taxon = case.get('taxon')
    if isinstance(taxon, dict):
        name = taxon.get('name')
        if isinstance(name, str) and name.strip():
            return name.strip().lower()
    return None


def _load_related_cases(paths, label):
    """[(case_dict, path)] for each --supporting-case / --counterexample-case value."""
    loaded = []
    for raw in (paths or []):
        case_path = pathlib.Path(raw) / 'case.json'
        if not case_path.is_file():
            raise ValueError('%s 案例不存在: %s' % (label, raw))
        loaded.append((_json(case_path), str(raw)))
    return loaded


def _lesson_scope_errors(scope, transferability, support, counter):
    """Fail-closed rules tying a lesson's claim to its INDEPENDENT support.

    Independence is deliberately hard to establish: two cases count as independent
    only when both record a taxon and the taxa differ.  A case with no recorded taxon
    contributes nothing, because "we cannot show these are independent" must not be
    rounded up to "independent" - that rounding is exactly how one case becomes a rule.
    """
    errors = []
    if scope not in GENERALIZATION_SCOPES:
        errors.append('generalization_scope 取值非法: %s (允许: %s)'
                      % (scope, '/'.join(GENERALIZATION_SCOPES)))
    if transferability not in TRANSFERABILITY_LEVELS:
        errors.append('transferability 取值非法: %s (允许: %s)'
                      % (transferability, '/'.join(TRANSFERABILITY_LEVELS)))
    if errors:
        return errors
    ids = [c.get('case_id') for c in support if c.get('case_id')]
    duplicated = sorted({item for item in ids if ids.count(item) > 1})
    if duplicated:
        errors.append('支持案例的 case_id 重复: %s（同一案例不得重复计数；'
                      '独立性必须同时满足 case_id 不同与类群不同）' % ','.join(duplicated))
    total = len(set(ids)) or 1
    if errors:
        return errors
    independent = len({_case_taxon_name(c) for c in support if _case_taxon_name(c)})
    if scope == 'single_case' and transferability != 'none':
        errors.append('generalization_scope=single_case 只能配 transferability=none')
    if counter and transferability != 'none':
        errors.append('存在反例（counterexample）时不得提高 transferability')
    if scope != 'single_case':
        if total < SCOPE_MIN_CASES[scope]:
            errors.append('generalization_scope=%s 需要至少 %d 个支持案例（当前 %d）'
                          % (scope, SCOPE_MIN_CASES[scope], total))
        if independent < SCOPE_MIN_CASES[scope]:
            errors.append('generalization_scope=%s 需要至少 %d 个不同类群的独立支持案例'
                          '（无法证明独立的不计数，当前 %d）'
                          % (scope, SCOPE_MIN_CASES[scope], independent))
    if independent < TRANSFERABILITY_MIN_CASES[transferability]:
        errors.append('transferability=%s 需要至少 %d 个不同类群的独立支持案例（当前 %d）'
                      % (transferability, TRANSFERABILITY_MIN_CASES[transferability], independent))
    return errors


def propose_lesson(args):
    lesson_dirs(); case = _json(pathlib.Path(args.case) / 'case.json')
    if case.get('decision', {}).get('status') == 'UNRESOLVED' and not args.allow_unresolved:
        raise ValueError('未解决案例不能直接提炼；使用 --allow-unresolved 仅生成候选并保留限制')
    lesson_id = safe_identifier(args.lesson_id or case.get('case_id'), 'lesson_id')
    source_case_type = case.get('case_type') or 'abnormal_case'
    allowed_domains = CASE_TYPE_LESSON_DOMAINS.get(source_case_type)
    if allowed_domains is None:
        raise ValueError('案例的 case_type 非法: %s' % source_case_type)
    lesson_domain = getattr(args, 'lesson_domain', None) or allowed_domains[0]
    if lesson_domain not in LESSON_DOMAINS:
        raise ValueError('lesson_domain 取值非法: %s (允许: %s)'
                         % (lesson_domain, '/'.join(LESSON_DOMAINS)))
    if lesson_domain not in allowed_domains:
        raise ValueError('case_type=%s 的案例只能产出 %s 类 lesson（当前请求 %s）：'
                         '工具/环境故障不得升为生物学或样本质量结论'
                         % (source_case_type, '/'.join(allowed_domains), lesson_domain))
    scope = getattr(args, 'generalization_scope', None) or 'single_case'
    transferability = getattr(args, 'transferability', None) or 'none'
    extra_support = _load_related_cases(getattr(args, 'supporting_case', None), '--supporting-case')
    counter = _load_related_cases(getattr(args, 'counterexample_case', None), '--counterexample-case')
    support = [(case, str(args.case))] + extra_support
    support_ids = [item.get('case_id') for item, _ in support if item.get('case_id')]
    counter_ids = [item.get('case_id') for item, _ in counter if item.get('case_id')]
    shared = sorted(set(support_ids) & set(counter_ids))
    if shared:
        raise ValueError('同一案例不能同时作为支持与反例: %s' % ', '.join(shared))
    errors = _lesson_scope_errors(scope, transferability, [item for item, _ in support],
                                 [item for item, _ in counter])
    if errors:
        raise ValueError('; '.join(errors))
    lesson = {'lesson_id': lesson_id, 'applicable_when': [case.get('issue', {}).get('type')],
              'not_applicable_when': [], 'diagnostic_clues': [case.get('issue', {}).get('user_observation', '')],
              'suggested_next_test': args.next_test, 'supporting_case_ids': support_ids,
              'counterexample_case_ids': counter_ids,
              'sources': [{'case_id': item.get('case_id'), 'path': path} for item, path in support],
              'generalization_scope': scope, 'transferability': transferability,
              'lesson_domain': lesson_domain, 'source_case_type': source_case_type,
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
                 'case_type': case.get('case_type', 'abnormal_case'),
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
    # fail-closed：缺失的推广上限按最保守值补齐，而不是当成"无限制"
    lesson.setdefault('generalization_scope', 'single_case')
    lesson.setdefault('transferability', 'none')
    if lesson['generalization_scope'] not in GENERALIZATION_SCOPES:
        raise ValueError('公共知识 generalization_scope 非法: %s (%s)'
                         % (lesson['generalization_scope'], expected_path))
    if lesson['transferability'] not in TRANSFERABILITY_LEVELS:
        raise ValueError('公共知识 transferability 非法: %s (%s)'
                         % (lesson['transferability'], expected_path))
    if lesson['generalization_scope'] == 'single_case' and lesson['transferability'] != 'none':
        raise ValueError('公共知识 single_case 只能配 transferability=none: %s' % expected_path)
    if lesson['transferability'] != 'none':
        case_ids = {item for item in (lesson.get('supporting_case_ids') or []) if item}
        needed = max(SCOPE_MIN_CASES[lesson['generalization_scope']],
                     TRANSFERABILITY_MIN_CASES[lesson['transferability']])
        if len(case_ids) < needed:
            raise ValueError('公共知识声明了 transferability=%s，需要至少 %d 个可区分支持案例'
                             '（当前 %d）: %s'
                             % (lesson['transferability'], needed, len(case_ids), expected_path))
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
    case_type = getattr(args, 'case_type', None) or 'abnormal_case'
    if case_type not in CASE_TYPES:
        raise ValueError('--case-type 取值非法: %s (允许: %s)' % (case_type, '/'.join(CASE_TYPES)))
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
            'case_type': case_type,
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


def case_anomaly(args):
    """Write one anomaly into case.json's anomalies[] -- the format the schema owns.

    Deliberately NOT a second record format: the entry goes into the same
    schema-validated case.json, and the resulting case is re-validated before the
    file is replaced.  `--event-action` must name an existing event action, so the
    anomaly's link to the original tool output is real rather than decorative.
    """
    # 缺依赖即失败：校验不可用时不得改动案例（宁可失败也不产生未校验写入）
    _require_jsonschema()
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

    # 业务矛盾不阻断写入（没有"改 decision"的命令，阻断会把人卡死），但必须告警：
    # 格式错误 = 记录不可用 → 拒绝；业务矛盾 = 记录可用但自相矛盾 → 显式提示。
    for finding in case_business_errors(case, root):
        print('⚠ [%s] %s' % (finding['code'], finding['message']))

    # 业务矛盾不阻断写入（没有"改 decision"的命令，阻断会把人卡死），但必须告警：
    # 格式错误 = 记录不可用 → 拒绝；业务矛盾 = 记录可用但自相矛盾 → 显式提示。
    for finding in case_business_errors(case, root):
        print('⚠ [%s] %s' % (finding['code'], finding['message']))

    tmp = root / 'case.json.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(case, fh, ensure_ascii=False, indent=2)
        fh.write('\n')
    os.replace(tmp, case_path)
    print('anomaly recorded: %s (%s/%s%s)'
          % (anomaly_id, entry['status'], entry['confidence'],
             ', reads_support=%s' % entry['reads_support'] if args.reads_support else ''))


class MissingDependency(RuntimeError):
    """必需依赖缺失。

    宁可失败也不产生弱化结论：同一结论必须对应同一证据路径。因此本 skill 不再为案例校验
    提供第二份实现——缺库时报错并给出安装命令，由 CLI 以退出码 3 返回（环境/依赖故障，
    与 cox1_id.py 的 3 同一含义）。
    """

    EXIT_CODE = 3


def _require_jsonschema():
    try:
        import jsonschema
    except ImportError as exc:
        raise MissingDependency(
            '缺少必需依赖 jsonschema（案例校验的唯一实现）：pip install jsonschema') from exc
    return jsonschema


def case_schema_errors(case):
    """按 `schemas/case.schema.json` 校验案例——**只有 jsonschema 这一条路径**。

    早期版本另有一份手写的“等价”实现，库缺失时静默切换。两份实现都对同一案例说
    `VALID`，却需要人肉同步，且已偏离过两次（`decision: null`、`case_id` 的
    `minLength`）。现在缺库即失败，不再降级。
    """
    schema_path = CASE_SCHEMA
    if not os.path.exists(schema_path):
        raise MissingDependency('缺少 schema 文件: %s（skill 安装不完整）' % schema_path)
    jsonschema = _require_jsonschema()
    try:
        jsonschema.validate(case, _json(schema_path))
    except jsonschema.ValidationError as exc:
        location = '/'.join(str(part) for part in exc.absolute_path) or '<root>'
        return ['schema 校验失败: %s (%s)' % (exc.message, location)]
    return []


# 跨字段业务规则（schema 之外）。每条规则有稳定错误码；FORMAT 与 BUSINESS 分开输出，
# 使 `case-validate` 的 VALID 只表示"格式合法 + 已定义业务规则无矛盾"，绝不隐含科学正确。
BUSINESS_RULES = {
    'DECISION_ANOMALY_CONFLICT':
        '案例级 decision.status 为 RESOLVED，却存在 status 为 UNRESOLVED 的异常（自相矛盾）',
    'ANOMALY_ID_DUPLICATE': 'anomalies[].id 重复，逐异常判定不可区分',
    'NORMAL_CASE_HAS_UNRESOLVED_ANOMALY':
        'case_type=normal_validation_case（已核实无异常）却挂着未解决异常，分类与内容矛盾',
    'ANOMALY_EVENT_UNKNOWN':
        'anomalies[].evidence_events 引用了 events.jsonl 中不存在的事件（证据链断开）',
    'EVENT_HYPOTHESIS_UNKNOWN': '事件 impact 引用了 hypotheses[] 中不存在的假设编号',
    'INPUT_FILE_MISSING': 'inputs[].path 在磁盘上不存在（环境性；仅提示，不改判读）',
    'INPUT_SHA256_MISMATCH': 'inputs[].sha256 与磁盘文件不一致（输入被替换或记录有误）',
    'EVENTS_UNPARSABLE': 'events.jsonl 存在无法解析为 JSON 对象的行',
}


def case_business_errors(case, root=None, verify_inputs=False):
    """跨字段业务矛盾（不是格式问题）。返回 [{'code','message','rule'}]。

    与格式校验分开的理由：schema 只表达结构与类型；"案例判 RESOLVED 却挂着 UNRESOLVED 异常"
    这类关系属业务自洽性，混进 schema 会把"格式合法"与"科学自洽"混为一谈。规则集合见
    `references/developer-contract.md` §1，措辞与 `references/evidence-standard.md` §7 一一对应。
    """
    findings = []

    def add(code, message, severity='error'):
        findings.append({'code': code, 'message': message, 'rule': BUSINESS_RULES[code],
                         'severity': severity})

    anomalies = [a for a in (case.get('anomalies') or []) if isinstance(a, dict)]
    decision = case.get('decision') if isinstance(case.get('decision'), dict) else {}
    unresolved = [str(a.get('id')) for a in anomalies if a.get('status') == 'UNRESOLVED']

    if decision.get('status') == 'RESOLVED' and unresolved:
        add('DECISION_ANOMALY_CONFLICT',
            '案例级 decision.status=RESOLVED，但异常 %s 为 UNRESOLVED' % ', '.join(unresolved))
    ids = [a.get('id') for a in anomalies if a.get('id')]
    duplicates = sorted({i for i in ids if ids.count(i) > 1})
    if duplicates:
        add('ANOMALY_ID_DUPLICATE', 'anomalies[].id 重复: %s' % ', '.join(duplicates))
    if case.get('case_type') == 'normal_validation_case' and unresolved:
        add('NORMAL_CASE_HAS_UNRESOLVED_ANOMALY',
            'case_type=normal_validation_case 但异常 %s 为 UNRESOLVED' % ', '.join(unresolved))

    # events.jsonl：可解析性、事件→异常、事件→假设的引用完整性
    event_actions, events = [], []
    if root is not None:
        try:
            events_path = safe_case_member(pathlib.Path(root),
                                           case.get('events_file', 'events.jsonl'), 'events_file')
        except ValueError:
            events_path = None
        if events_path is not None and events_path.is_file():
            for number, line in enumerate(events_path.read_text(encoding='utf-8').splitlines(), 1):
                if not line.strip():
                    continue
                try:
                    record = json.loads(line)
                except ValueError:
                    add('EVENTS_UNPARSABLE', 'events.jsonl 第 %d 行不是合法 JSON' % number)
                    continue
                if not isinstance(record, dict):
                    add('EVENTS_UNPARSABLE', 'events.jsonl 第 %d 行不是 JSON 对象' % number)
                    continue
                events.append(record)
                if record.get('action'):
                    event_actions.append(str(record['action']))
    known_actions = set(event_actions)
    for anomaly in anomalies:
        for action in (anomaly.get('evidence_events') or []):
            if str(action) not in known_actions:
                add('ANOMALY_EVENT_UNKNOWN', '异常 %s 引用的事件不存在: %s'
                    % (anomaly.get('id'), action))
    hypothesis_ids = {str(h.get('id')) for h in (case.get('hypotheses') or [])
                      if isinstance(h, dict) and h.get('id')}
    for record in events:
        for token in re.findall(r'\bH\d+\b', str(record.get('impact') or '')):
            if token not in hypothesis_ids:
                add('EVENT_HYPOTHESIS_UNKNOWN', '事件 %s 的 impact 引用不存在的假设: %s'
                    % (record.get('action'), token))

    # 输入文件：存在性总是校验；哈希只在 --verify-inputs 时比对（大文件重新哈希代价高）
    for item in (case.get('inputs') or []):
        if not isinstance(item, dict) or not item.get('path'):
            continue
        path = pathlib.Path(str(item['path']))
        if not path.is_file():
            # 文件不在这台机器上不等于记录自相矛盾（归档/换机很常见）→ note，不改判读
            add('INPUT_FILE_MISSING', '输入文件不存在: %s (%s)' % (item['path'], item.get('role')),
                severity='note')
        elif verify_inputs and item.get('sha256'):
            actual = file_sha256(path)
            if actual != item['sha256']:
                add('INPUT_SHA256_MISMATCH', '输入 %s 的 sha256 与记录不一致（记录 %s…，实际 %s…）'
                    % (item.get('role'), str(item['sha256'])[:12], actual[:12]))
    return findings


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
    business = case_business_errors(case, root, verify_inputs=getattr(args, 'verify_inputs', False))
    fatal = [f for f in business if f.get('severity', 'error') == 'error']
    print('FORMAT: %s' % ('OK' if not errors else 'INVALID (%d 项)' % len(errors)))
    for error in errors:
        print('- [SCHEMA] %s' % error)
    print('BUSINESS: %s' % ('无矛盾' if not business
                            else '%d 项矛盾（其中 %d 项为 error 级）' % (len(business), len(fatal))))
    for finding in business:
        print('- [%s/%s] %s' % (finding['code'], finding.get('severity', 'error'), finding['message']))
    if errors or fatal:
        # 判读词本身必须保持独立字面量：它被 doc-contract 审计保护，也便于机器判读。
        print('INVALID')
        print('  cause: %s' % ('format' if errors else 'business'))
        return 1
    print('VALID')
    print('  格式合法 + 已定义业务规则无矛盾；不证明科学结论')
    return 0

def case_reference(args):
    """Link a registered public reference into the case's evidence chain.

    The case must cite the registry entry (id + accession with version + file hash),
    not "a close relative": a bare citation cannot be reproduced once GenBank updates
    the record.  An unknown id, or a purpose the record never declared, fails closed so
    dangling citations cannot be written.
    """
    # 缺依赖即失败：校验不可用时不得改动案例（宁可失败也不产生未校验写入）
    _require_jsonschema()
    root = pathlib.Path(args.directory).resolve()
    case_path = root / 'case.json'
    if not case_path.is_file():
        raise ValueError('案例不存在: %s' % case_path)
    with open(case_path, encoding='utf-8') as fh:
        case = json.load(fh)
    if not isinstance(case, dict):
        raise ValueError('case.json 必须是 JSON 对象')
    registry_dir = pathlib.Path(getattr(args, 'registry', None) or REFERENCE_REGISTRY_DIR)
    registry_file = registry_dir / REFERENCE_REGISTRY_FILE
    if not registry_file.is_file():
        raise ValueError('参考登记表不存在: %s（先用 scripts/reference_registry.py register/acquire 登记）'
                         % registry_file)
    with open(registry_file, encoding='utf-8') as fh:
        registry = json.load(fh)
    record = next((item for item in (registry.get('references') or [])
                   if isinstance(item, dict) and item.get('id') == args.reference_id), None)
    if record is None:
        raise ValueError('参考登记表中没有该 id: %s（不允许引用未登记的参考，否则证据链断开）'
                         % args.reference_id)
    purpose = str(args.purpose).strip()
    declared = list(record.get('purposes') or [])
    if purpose not in declared:
        raise ValueError('该参考登记时未声明用途 %s（已声明: %s）—— 等级限制用途，'
                         '不得事后扩大使用范围'
                         % (purpose, '/'.join(declared) or '无'))
    entry = {'reference_id': record.get('id'), 'accession': record.get('accession'),
             'level': record.get('level'), 'purpose': purpose,
             'file_sha256': record.get('file_sha256'), 'source': record.get('source')}
    entries = case.setdefault('references', [])
    if not isinstance(entries, list):
        raise ValueError('case.json 的 references 必须是数组')
    matches = [index for index, item in enumerate(entries)
               if isinstance(item, dict) and item.get('reference_id') == entry['reference_id']
               and item.get('purpose') == purpose]
    if matches and not getattr(args, 'update', False):
        raise ValueError('该参考与用途已关联: %s/%s (如需替换请加 --update)'
                         % (entry['reference_id'], purpose))
    if matches:
        entries[matches[0]] = entry
    else:
        entries.append(entry)
    errors = case_schema_errors(case)
    if errors:
        raise ValueError('写回后案例不符合 schema, 已放弃写入: %s' % '; '.join(errors))
    tmp = root / 'case.json.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(case, fh, ensure_ascii=False, indent=2)
        fh.write('\n')
    os.replace(tmp, case_path)
    print('参考已关联: %s (%s, level=%s, purpose=%s, hash=%s)'
          % (entry['reference_id'], entry['accession'], entry['level'], purpose,
             str(entry['file_sha256'])[:12]))


def case_report(args):
    root = pathlib.Path(args.directory).resolve()
    with open(root / 'case.json', encoding='utf-8') as fh: case = json.load(fh)
    events = []
    events_path = safe_case_member(root, case.get('events_file', 'events.jsonl'), 'events_file')
    with open(events_path, encoding='utf-8') as fh:
        for line in fh:
            if line.strip():
                events.append(json.loads(line))
    anomalies = case.get('anomalies') or []
    supported = [a for a in anomalies if a.get('status') in ('RESOLVED', 'NO_CHANGE')]
    unsupported = [a for a in anomalies if a.get('status') == 'UNRESOLVED']
    decision = case['decision']

    lines = ['# 诊断案例 %s' % case['case_id'], '',
             '- 案例类型：`%s`' % case.get('case_type', 'abnormal_case'), '']
    lines += ['## 1. 观察事实（事件）', '']
    if events:
        for event in events:
            lines.append('- `%s`：%s；影响：%s' % (event.get('action'), event.get('result'),
                                            event.get('impact') or '未记录'))
    else:
        lines.append('- 尚无事件记录')

    lines += ['', '## 2. 证据矩阵（逐异常结论）', '']
    if anomalies:
        lines += ['| 异常 | claim | 状态 | 置信 | reads 支持 | 证据事件 |',
                  '|---|---|---|---|---|---|']
        for item in anomalies:
            lines.append('| `%s` | %s | `%s` | `%s` | `%s` | %s |'
                         % (item.get('id'), item.get('claim'), item.get('status'),
                            item.get('confidence'), item.get('reads_support') or 'NOT_ASSESSED',
                            ', '.join('`%s`' % action
                                      for action in (item.get('evidence_events') or [])) or '—'))
    else:
        lines.append('尚无逐异常判定：请用 `case-anomaly` 逐条写入 `anomalies[]`；'
                     '案例级 `decision` 不能代替逐异常结论。')

    lines += ['', '## 3. 有证据支持的结论', '']
    if supported:
        for item in supported:
            lines.append('- `%s`：%s（`%s`/`%s`，reads=`%s`）'
                         % (item.get('id'), item.get('claim'), item.get('status'),
                            item.get('confidence'), item.get('reads_support') or 'NOT_ASSESSED'))
    else:
        lines.append('- 目前没有达到 `RESOLVED`/`NO_CHANGE` 的逐异常结论。')

    lines += ['', '## 4. 无证据支持的声明（不得当作结论使用）', '']
    if unsupported:
        for item in unsupported:
            lines.append('- `%s`：%s（`UNRESOLVED`/`%s`）'
                         % (item.get('id'), item.get('claim'), item.get('confidence')))
        lines += ['', '以上条目**尚缺**能改变判定或置信档位的证据；不得在报告或结论中写成\u201c已排除\u201d。']
    else:
        lines.append('- 无（当前没有 `UNRESOLVED` 的逐异常记录）')

    lines += ['', '## 5. 下一步最小实验', '']
    if unsupported:
        for item in unsupported:
            lines.append('- `%s`：需要能改变该结论或置信档位的最小检查'
                         '（缺什么输入见 `diagnostic-decision-tree.md` §3/§6）'
                         % item.get('id'))
    elif not anomalies:
        lines.append('- 先完成逐异常判定（`case-anomaly`）后再确定下一步。')
    else:
        lines.append('- 无：当前异常均已判定。')

    lines += ['', '## 6. 当前判定（案例级汇总，不代替逐异常结论）', '',
              '- 状态：`%s`' % decision['status'],
              '- 置信等级：`%s`' % decision['confidence'],
              '- 理由：%s' % decision['rationale']]
    hypotheses = case.get('hypotheses') or []
    if hypotheses:
        lines += ['', '## 7. 假设与未测项', '']
        for item in hypotheses:
            lines.append('- `%s`：%s（支持 %d / 反证 %d / 未测 %d）'
                         % (item.get('id'), item.get('explanation'),
                            len(item.get('support') or []), len(item.get('against') or []),
                            len(item.get('unknown') or [])))
            for unknown in (item.get('unknown') or []):
                lines.append('  - 未测：%s' % unknown)

    lines += ['', '## 8. 参考（已登记：版本 + hash 固定）', '']
    references = case.get('references') or []
    if references:
        lines += ['| 参考 | accession(含版本) | 等级 | 用途 | file_sha256 |',
                  '|---|---|---|---|---|']
        for item in references:
            lines.append('| `%s` | %s | %s | %s | `%s` |'
                         % (item.get('reference_id'), item.get('accession'),
                            item.get('level') or '—', item.get('purpose'),
                            str(item.get('file_sha256') or '')[:12]))
    else:
        lines.append('- 本案例未引用已登记参考（如引用过参考，必须用 `case-reference` 登记，'
                     '否则无法复现：见 `references/reference-policy.md`）')

    lines += ['', '## 9. 证据边界', '',
              '缺少 reads 时不得写成 raw-read-supported；`case-validate` 只校验记录格式，'
              '不证明科学结论；阴性结果的强度取决于覆盖与读长'
              '（见 `evidence-standard.md` §3.2）。']
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
        ap = argparse.ArgumentParser(); ap.add_argument('--case', required=True); ap.add_argument('--next-test', required=True); ap.add_argument('--lesson-id'); ap.add_argument('--allow-unresolved', action='store_true'); ap.add_argument('--lesson-domain', dest='lesson_domain', default=None, choices=list(LESSON_DOMAINS)); ap.add_argument('--supporting-case', dest='supporting_case', action='append', default=[], metavar='DIR'); ap.add_argument('--counterexample-case', dest='counterexample_case', action='append', default=[], metavar='DIR'); ap.add_argument('--generalization-scope', dest='generalization_scope', default='single_case', choices=list(GENERALIZATION_SCOPES)); ap.add_argument('--transferability', default='none', choices=list(TRANSFERABILITY_LEVELS))
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
    elif cmd in ('case-init', 'case-event', 'case-anomaly', 'case-reference', 'case-validate', 'case-report'):
        ap = argparse.ArgumentParser()
        if cmd == 'case-init':
            ap.add_argument('directory'); ap.add_argument('--case-id'); ap.add_argument('--case-type', dest='case_type', default='abnormal_case', choices=list(CASE_TYPES)); ap.add_argument('--issue', required=True); ap.add_argument('--observation', required=True); ap.add_argument('--taxon'); ap.add_argument('--input', nargs=2, action='append', default=[], metavar=('ROLE', 'PATH')); ap.add_argument('--hypothesis', action='append', default=[], metavar='TEXT')
            try: case_init(ap.parse_args(sys.argv[2:]))
            except (ValueError, OSError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
        elif cmd == 'case-event':
            ap.add_argument('directory'); ap.add_argument('--action', required=True); ap.add_argument('--result', required=True); ap.add_argument('--impact', default=''); ap.add_argument('--command', default=''); ap.add_argument('--tool-version', default=''); ap.add_argument('--motivation', default='')
            case_event(ap.parse_args(sys.argv[2:]))
        elif cmd == 'case-anomaly':
            ap.add_argument('directory'); ap.add_argument('--id', required=True); ap.add_argument('--claim', required=True); ap.add_argument('--status', required=True); ap.add_argument('--confidence', required=True); ap.add_argument('--reads-support', dest='reads_support', default=None); ap.add_argument('--event-action', dest='event_action', action='append', default=[], metavar='ACTION'); ap.add_argument('--update', action='store_true')
            try: case_anomaly(ap.parse_args(sys.argv[2:]))
            except MissingDependency as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(MissingDependency.EXIT_CODE)
            except (ValueError, OSError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
        elif cmd == 'case-reference':
            ap.add_argument('directory'); ap.add_argument('--reference-id', dest='reference_id', required=True); ap.add_argument('--purpose', required=True); ap.add_argument('--registry', default=None); ap.add_argument('--update', action='store_true')
            try: case_reference(ap.parse_args(sys.argv[2:]))
            except MissingDependency as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(MissingDependency.EXIT_CODE)
            except (ValueError, OSError, json.JSONDecodeError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
        elif cmd == 'case-validate':
            ap.add_argument('directory'); ap.add_argument('--verify-inputs', dest='verify_inputs', action='store_true')
            try: sys.exit(case_validate(ap.parse_args(sys.argv[2:])))
            except MissingDependency as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(MissingDependency.EXIT_CODE)
            except (ValueError, OSError, json.JSONDecodeError) as exc: print('✗ %s' % exc, file=sys.stderr); sys.exit(1)
        else:
            ap.add_argument('directory'); case_report(ap.parse_args(sys.argv[2:]))
    else:
        print('未知命令: %s' % cmd); print(__doc__)

if __name__ == '__main__':
    main()

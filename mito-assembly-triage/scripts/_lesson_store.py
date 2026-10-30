"""Authorized local lesson memory; JSON is authoritative and FTS5 is derived."""
from contextlib import closing, contextmanager
from copy import deepcopy
from functools import lru_cache
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import stat

from _case_store import CaseStore, MissingDependency, now, read_json, snapshot
from _evidence_io import atomic_write_text

SCHEMA = Path(__file__).resolve().parents[1] / 'schemas/lesson.schema.json'
STATUSES = ('candidate', 'verified', 'rejected', 'deprecated', 'withdrawn', 'superseded')
DRAFT_REQUIRED = {
    'title', 'claim', 'lesson_domain', 'diagnostic_clues', 'suggested_next_test',
    'applicable_when', 'not_applicable_when',
}
DRAFT_OPTIONAL = {'counterexample_case_ids', 'conflicts'}
TRANSITIONS = {
    'candidate': {'verified', 'rejected'},
    'verified': {'deprecated', 'withdrawn'},
    'rejected': set(),
    'deprecated': set(),
    'withdrawn': set(),
    'superseded': set(),
}
CONTENT_FIELDS = (
    'title', 'claim', 'lesson_domain', 'source_case_type', 'diagnostic_clues',
    'suggested_next_test', 'applicable_when', 'not_applicable_when',
    'supporting_case_ids', 'counterexample_case_ids', 'sources', 'conflicts',
    'generalization_scope', 'transferability',
)
PRIVATE_PATTERNS = (
    re.compile(r'\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b', re.I),
    re.compile(r'(?i)\b(?:api[_-]?key|token|password|secret)\s*[:=]\s*\S+'),
    re.compile(r'(?<![:\w/])/(?:[^/\s]+/)+[^/\s]+'),
    re.compile(r'(?m)^\s*>[^\n]+\n\s*[ACGTRYSWKMBDHVNacgtryswkmbdhvn\s]{20,}'),
)


@lru_cache(maxsize=1)
def validator():
    try:
        from jsonschema import Draft202012Validator
    except ImportError as exc:
        raise MissingDependency('缺少 jsonschema；安装: python -m pip install jsonschema') from exc
    schema = read_json(SCHEMA)
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)


@lru_cache(maxsize=1)
def ensure_fts5():
    try:
        with closing(sqlite3.connect(':memory:')) as connection:
            connection.execute('CREATE VIRTUAL TABLE fts_probe USING fts5(value)')
    except sqlite3.Error as exc:
        raise MissingDependency('Python SQLite 未提供 FTS5；请使用带 FTS5 的 Python 3') from exc


def content_hash(revision):
    content = {field: revision[field] for field in CONTENT_FIELDS}
    encoded = json.dumps(content, ensure_ascii=False, sort_keys=True, separators=(',', ':'),
                         allow_nan=False).encode('utf-8')
    return hashlib.sha256(encoded).hexdigest()


def privacy_check(value, path='record'):
    if isinstance(value, dict):
        for key, item in value.items():
            privacy_check(item, '%s.%s' % (path, key))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            privacy_check(item, '%s[%d]' % (path, index))
    elif isinstance(value, str):
        if any(pattern.search(value) for pattern in PRIVATE_PATTERNS) or re.fullmatch(
                r'[ACGTNacgtn\s]{40,}', value):
            raise ValueError('经验草稿的 %s 可能含路径、凭据或样本序列；请改为泛化描述' % path)


def validate_revision(revision):
    errors = sorted(validator().iter_errors(revision), key=lambda error: str(list(error.path)))
    if errors:
        raise ValueError('LESSON FORMAT: ' + '; '.join(
            '%s: %s' % ('/'.join(map(str, error.path)), error.message) for error in errors[:8]))
    if len(revision['sources']) != 1 or revision['supporting_case_ids'] != [revision['sources'][0]['case_id']]:
        raise ValueError('MVP 每条修订必须且只能关联一个来源案例')
    if revision['generalization_scope'] != 'single_case' or revision['transferability'] != 'none':
        raise ValueError('MVP 单案例经验只允许 single_case/none')
    if revision['source_case_type'] == 'tool_failure_case' and revision['lesson_domain'] != 'tool':
        raise ValueError('tool_failure_case 只能产生 tool 领域经验')
    if revision['lesson_id'] in revision['conflicts']:
        raise ValueError('经验不能与自身冲突')
    if revision['content_hash'] != content_hash(revision):
        raise ValueError('经验内容哈希不匹配')
    if revision['revision'] == 1 and revision['supersedes'] is not None:
        raise ValueError('首修订不得声明 supersedes')
    if revision['revision'] > 1:
        expected = {'lesson_id': revision['lesson_id'], 'revision': revision['revision'] - 1}
        if revision['supersedes'] != expected:
            raise ValueError('修订链必须指向前一修订')
    if revision['review_history'] and revision['last_reviewed'] is None:
        raise ValueError('有审核历史时必须记录 last_reviewed')
    return revision


def validate_document(document):
    if not isinstance(document, dict) or set(document) != {
            'format', 'lesson_id', 'current_revision', 'revisions'}:
        raise ValueError('经验文档须仅含 format/lesson_id/current_revision/revisions')
    if document['format'] != 'mito-lesson-document-1':
        raise ValueError('不支持的经验文档格式')
    revisions = document['revisions']
    if not isinstance(revisions, list) or not revisions:
        raise ValueError('经验文档至少需要一个修订')
    for index, revision in enumerate(revisions, 1):
        validate_revision(revision)
        if revision['lesson_id'] != document['lesson_id'] or revision['revision'] != index:
            raise ValueError('经验 ID 或修订序号不连续')
        if index < len(revisions) and revision['validation_status'] != 'superseded':
            raise ValueError('历史修订必须标记 superseded')
    if document['current_revision'] != len(revisions):
        raise ValueError('current_revision 与修订数量不符')
    if revisions[-1]['validation_status'] == 'superseded':
        raise ValueError('当前修订不得标记 superseded')
    return document


def load_draft(path):
    draft = read_json(path)
    if not isinstance(draft, dict):
        raise ValueError('经验草稿必须为 JSON object')
    missing = DRAFT_REQUIRED - set(draft)
    unknown = set(draft) - DRAFT_REQUIRED - DRAFT_OPTIONAL
    if missing or unknown:
        raise ValueError('经验草稿缺字段: %s；未知字段: %s' % (
            ','.join(sorted(missing)) or '无', ','.join(sorted(unknown)) or '无'))
    privacy_check(draft)
    result = deepcopy(draft)
    result.setdefault('counterexample_case_ids', [])
    result.setdefault('conflicts', [])
    return result


def source_case(directory):
    store = CaseStore(directory)
    case = store.load(verify_files=False)
    if not case['conclusions']:
        raise ValueError('来源案例至少需要一条结论')
    digest = snapshot(store.path, 'case', 'ledger')
    return case, digest


def check_unresolved(case, allow):
    if not allow and any(item['status'] == 'UNRESOLVED' for item in case['conclusions']):
        raise ValueError('案例含 UNRESOLVED 结论；仅生成候选时也须显式 --allow-unresolved')


def next_lesson_id(documents):
    numbers = [int(match.group(1)) for document in documents
               for match in [re.fullmatch(r'L(\d{6})', document['lesson_id'])] if match]
    return 'L%06d' % (max(numbers, default=0) + 1)


def build_revision(draft, case, digest, lesson_id, number, previous=None):
    timestamp = now()
    revision = dict(
        lesson_id=lesson_id,
        revision=number,
        validation_status='candidate',
        lesson_domain=draft['lesson_domain'],
        source_case_type=case['case_type'],
        title=draft['title'],
        claim=draft['claim'],
        diagnostic_clues=draft['diagnostic_clues'],
        suggested_next_test=draft['suggested_next_test'],
        applicable_when=draft['applicable_when'],
        not_applicable_when=draft['not_applicable_when'],
        supporting_case_ids=[case['case_id']],
        counterexample_case_ids=draft['counterexample_case_ids'],
        sources=[{'case_id': case['case_id'], 'path': digest['path'], 'case_sha256': digest['sha256']}],
        conflicts=draft['conflicts'],
        generalization_scope='single_case',
        transferability='none',
        review_history=[],
        created_at=previous['created_at'] if previous else timestamp,
        updated_at=timestamp,
        last_reviewed=None,
        supersedes={'lesson_id': lesson_id, 'revision': number - 1} if previous else None,
        content_hash='',
    )
    revision['content_hash'] = content_hash(revision)
    return validate_revision(revision)


class LessonStore:
    def __init__(self, directory):
        requested = Path(directory)
        if requested.is_symlink():
            raise ValueError('经验库目录不得为符号链接')
        self.directory = requested.resolve()
        self.lessons = self.directory / 'lessons'
        self.lockpath = self.directory / '.memory.lock'
        self.index = self.directory / 'index.sqlite3'

    @contextmanager
    def locked(self, create=False, authorized=False):
        if not authorized:
            raise ValueError('写入本地经验库须显式 --authorize-local')
        validator()
        ensure_fts5()
        if create:
            self.directory.mkdir(parents=True, exist_ok=True, mode=0o700)
        if not self.directory.is_dir():
            raise ValueError('经验库不存在: %s' % self.directory)
        fd = os.open(self.lockpath, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
        try:
            info = os.fstat(fd)
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise ValueError('经验库锁必须为独立普通文件')
            fcntl.flock(fd, fcntl.LOCK_EX)
            if self.lessons.is_symlink() or self.index.is_symlink():
                raise ValueError('经验目录和索引不得为符号链接')
            self.lessons.mkdir(exist_ok=True, mode=0o700)
            yield
        finally:
            os.close(fd)

    def lesson_path(self, lesson_id):
        if not isinstance(lesson_id, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', lesson_id):
            raise ValueError('非法经验 ID')
        return self.lessons / (lesson_id + '.json')

    def source_paths(self):
        if not self.lessons.is_dir():
            return []
        paths = sorted(self.lessons.glob('*.json'))
        for path in paths:
            if path.is_symlink() or not path.is_file():
                raise ValueError('经验源文件必须为普通文件: %s' % path)
        return paths

    def fingerprint(self):
        digest = hashlib.sha256()
        for path in self.source_paths():
            digest.update(path.name.encode('utf-8'))
            digest.update(b'\0')
            with path.open('rb') as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b''):
                    digest.update(block)
            digest.update(b'\0')
        return digest.hexdigest()

    def documents(self):
        return [validate_document(read_json(path)) for path in self.source_paths()]

    def load(self, lesson_id):
        path = self.lesson_path(lesson_id)
        if not path.exists():
            raise ValueError('经验不存在: %s' % lesson_id)
        document = validate_document(read_json(path))
        if document['lesson_id'] != lesson_id:
            raise ValueError('文件名与经验 ID 不一致')
        return document

    def save(self, document):
        validate_document(document)
        path = self.lesson_path(document['lesson_id'])
        if path.is_symlink():
            raise ValueError('经验源文件不得为符号链接')
        atomic_write_text(path, json.dumps(document, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

    def check_duplicate(self, revision, exclude=None):
        for document in self.documents():
            if document['lesson_id'] == exclude:
                continue
            if document['revisions'][-1]['content_hash'] == revision['content_hash']:
                raise ValueError('相同经验已存在: %s' % document['lesson_id'])

    def propose(self, case_directory, draft_path, authorized=False, allow_unresolved=False):
        draft = load_draft(draft_path)
        case, digest = source_case(case_directory)
        check_unresolved(case, allow_unresolved)
        build_revision(draft, case, digest, next_lesson_id(self.documents()), 1)
        with self.locked(create=True, authorized=authorized):
            existing = self.documents()
            lesson_id = next_lesson_id(existing)
            revision = build_revision(draft, case, digest, lesson_id, 1)
            self.check_duplicate(revision)
            document = dict(format='mito-lesson-document-1', lesson_id=lesson_id,
                            current_revision=1, revisions=[revision])
            self.save(document)
            self.rebuild_unlocked()
            return revision

    def review(self, lesson_id, status, reviewer, reason, authorized=False):
        if status not in TRANSITIONS or status == 'superseded':
            raise ValueError('非法审核状态')
        if not reviewer.strip() or not reason.strip():
            raise ValueError('审核人和审核理由不能为空')
        with self.locked(authorized=authorized):
            document = self.load(lesson_id)
            revision = document['revisions'][-1]
            old_status = revision['validation_status']
            if status not in TRANSITIONS[old_status]:
                raise ValueError('不允许的状态转换: %s → %s' % (old_status, status))
            if status == 'verified':
                self.verify_sources(revision)
            timestamp = now()
            revision['review_history'].append(dict(from_status=old_status, to_status=status,
                reviewer=reviewer, reason=reason, reviewed_at=timestamp))
            revision['validation_status'] = status
            revision['last_reviewed'] = timestamp
            revision['updated_at'] = timestamp
            self.save(document)
            self.rebuild_unlocked()
            return revision

    def verify_sources(self, revision):
        for source in revision['sources']:
            case_path = Path(source['path'])
            store = CaseStore(case_path.parent)
            with store.locked():
                digest = snapshot(case_path, 'case', 'ledger')
                if digest['sha256'] != source['case_sha256']:
                    raise ValueError('来源案例自提案后已变化；请生成新修订')
                case = store.load(verify_files=False)
                if case['case_id'] != source['case_id']:
                    raise ValueError('来源案例 ID 不匹配')
                if any(item['status'] == 'UNRESOLVED' for item in case['conclusions']):
                    raise ValueError('含 UNRESOLVED 结论的案例不能审核为 verified')

    def supersede(self, lesson_id, case_directory, draft_path, authorized=False, allow_unresolved=False):
        draft = load_draft(draft_path)
        case, digest = source_case(case_directory)
        check_unresolved(case, allow_unresolved)
        with self.locked(authorized=authorized):
            document = self.load(lesson_id)
            previous = document['revisions'][-1]
            previous_source = previous['sources'][0]
            if (previous['supporting_case_ids'] != [case['case_id']] or
                    Path(previous_source['path']).resolve() != Path(digest['path']).resolve()):
                raise ValueError('MVP 新修订必须沿用同一来源案例文件；新样本请单独 propose')
            revision = build_revision(draft, case, digest, lesson_id, previous['revision'] + 1, previous)
            if revision['content_hash'] == previous['content_hash']:
                raise ValueError('新修订的来源和科学内容与上一版相同')
            self.check_duplicate(revision, exclude=lesson_id)
            previous['validation_status'] = 'superseded'
            previous['updated_at'] = now()
            document['revisions'].append(revision)
            document['current_revision'] = revision['revision']
            self.save(document)
            self.rebuild_unlocked()
            return revision

    def show(self, lesson_id, revision_number=None):
        document = self.load(lesson_id)
        number = document['current_revision'] if revision_number is None else revision_number
        if number < 1 or number > len(document['revisions']):
            raise ValueError('修订不存在: %s' % number)
        return document['revisions'][number - 1]

    def rebuild(self, authorized=False):
        with self.locked(authorized=authorized):
            return self.rebuild_unlocked()

    def rebuild_unlocked(self):
        documents = self.documents()
        fingerprint = self.fingerprint()
        with closing(sqlite3.connect(self.index)) as connection:
            connection.execute('CREATE TABLE IF NOT EXISTS index_state '
                               '(id INTEGER PRIMARY KEY CHECK(id=1), fingerprint TEXT NOT NULL)')
            connection.execute('CREATE VIRTUAL TABLE IF NOT EXISTS lesson_index USING fts5('
                               'lesson_id UNINDEXED, revision UNINDEXED, status UNINDEXED, '
                               'domain UNINDEXED, source_case_type UNINDEXED, taxa UNINDEXED, '
                               'data_types UNINDEXED, title, claim, clues, next_test, '
                               'content_hash UNINDEXED, tokenize="unicode61")')
            connection.execute('DELETE FROM lesson_index')
            connection.execute('DELETE FROM index_state')
            for document in documents:
                revision = document['revisions'][-1]
                connection.execute('INSERT INTO lesson_index VALUES (?,?,?,?,?,?,?,?,?,?,?,?)', (
                    revision['lesson_id'], revision['revision'], revision['validation_status'],
                    revision['lesson_domain'], revision['source_case_type'],
                    json.dumps(revision['applicable_when']['taxa'], ensure_ascii=False),
                    json.dumps(revision['applicable_when']['data_types'], ensure_ascii=False),
                    revision['title'], revision['claim'], '\n'.join(revision['diagnostic_clues']),
                    revision['suggested_next_test'], revision['content_hash']))
            connection.execute('INSERT INTO index_state VALUES (1,?)', (fingerprint,))
            connection.commit()
        os.chmod(self.index, 0o600)
        return len(documents)

    def search(self, query='', status='verified', domain=None, taxon=None, data_type=None, limit=5):
        ensure_fts5()
        if limit < 1:
            raise ValueError('limit 必须大于 0')
        if not self.index.exists():
            raise ValueError('经验索引不存在；请先 propose 或 rebuild')
        if self.index.is_symlink():
            raise ValueError('经验索引不得为符号链接')
        fingerprint = self.fingerprint()
        tokens = re.findall(r'\w+', query, flags=re.UNICODE)
        if query.strip() and not tokens:
            return []
        match = ' AND '.join('"' + token.replace('"', '""') + '"' for token in tokens)
        where, parameters = [], []
        if match:
            where.append('lesson_index MATCH ?')
            parameters.append(match)
        if status != 'all':
            if status not in STATUSES or status == 'superseded':
                raise ValueError('非法检索状态')
            where.append('status = ?')
            parameters.append(status)
        if domain:
            where.append('domain = ?')
            parameters.append(domain)
        sql = 'SELECT lesson_id,revision,status,domain,source_case_type,taxa,data_types,title,claim,clues,next_test,bm25(lesson_index) FROM lesson_index'
        if where:
            sql += ' WHERE ' + ' AND '.join(where)
        sql += ' ORDER BY bm25(lesson_index),lesson_id'
        uri = self.index.resolve().as_uri() + '?mode=ro'
        with closing(sqlite3.connect(uri, uri=True)) as connection:
            recorded = connection.execute('SELECT fingerprint FROM index_state WHERE id=1').fetchone()
            if recorded is None or recorded[0] != fingerprint:
                raise ValueError('经验索引已过期；请运行 rebuild')
            rows = connection.execute(sql, parameters).fetchall()
        results = []
        wanted_taxon = taxon.casefold() if taxon else None
        wanted_data = data_type.casefold() if data_type else None
        for row in rows:
            taxa = json.loads(row[5])
            data_types = json.loads(row[6])
            if wanted_taxon and wanted_taxon not in {item.casefold() for item in taxa}:
                continue
            if wanted_data and wanted_data not in {item.casefold() for item in data_types}:
                continue
            results.append(dict(lesson_id=row[0], revision=row[1], status=row[2], domain=row[3],
                source_case_type=row[4], taxa=taxa, data_types=data_types, title=row[7], claim=row[8],
                clues=row[9].splitlines(), next_test=row[10], score=row[11]))
            if len(results) == limit:
                break
        return results

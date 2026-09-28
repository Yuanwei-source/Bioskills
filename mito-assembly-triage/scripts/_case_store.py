"""Task-local, versioned evidence ledger. No execution, upload or knowledge store.

One atomic case.json contains immutable events and conclusion revision history.
The lock serializes cooperating writers. Hashes establish content identity, not
the truth of a scientific claim or the authenticity of imported execution metadata.
"""
from contextlib import contextmanager
from copy import deepcopy
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import stat

from _evidence_io import atomic_write_text, protect_output

SCHEMA = Path(__file__).resolve().parents[1] / 'schemas/task-case.schema.json'


class MissingDependency(RuntimeError):
    pass


def now():
    return datetime.now(timezone.utc).isoformat()


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('重复 JSON 字段: %s' % key)
            result[key] = value
        return result

    def constant(value):
        raise ValueError('非法 JSON 数值: %s' % value)
    with open(path, encoding='utf-8') as handle:
        return json.load(handle, object_pairs_hook=pairs, parse_constant=constant)


def validator():
    try:
        from jsonschema import Draft202012Validator
    except ImportError as exc:
        raise MissingDependency('缺少 jsonschema；安装: python -m pip install jsonschema') from exc
    schema = read_json(SCHEMA)
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema)


def snapshot(path, identifier, role):
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError('只接受普通文件: %s' % path)
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError('只接受普通文件: %s' % path)
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            digest.update(chunk)
        after = os.fstat(handle.fileno())
    if (before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (after.st_size, after.st_mtime_ns, after.st_ctime_ns):
        raise ValueError('文件在计算哈希期间发生变化: %s' % path)
    return dict(id=identifier, role=role, path=str(path), sha256=digest.hexdigest(), size=after.st_size)


def files_in(case):
    yield from case['inputs']
    for event in case['events']:
        yield from event['artifacts']


def unique(items, label):
    result = {}
    for item in items:
        if item['id'] in result:
            raise ValueError('%s ID 重复: %s' % (label, item['id']))
        result[item['id']] = item
    return result


def candidate_status(results):
    if results and all(result == 'pass' for result in results):
        return 'verified'
    if results and all(result == 'fail' for result in results):
        return 'rejected'
    if results:
        return 'inconclusive'
    return 'proposed'


def validate(case, verify_files=False):
    errors = sorted(validator().iter_errors(case), key=lambda error: str(list(error.path)))
    if errors:
        raise ValueError('FORMAT: ' + '; '.join('%s: %s' % ('/'.join(map(str, e.path)), e.message)
                                              for e in errors[:8]))
    inputs = unique(case['inputs'], 'input')
    events = unique(case['events'], 'event')
    unique(case['conclusions'], 'conclusion')
    candidates = unique(case.get('candidates', []), 'candidate')
    for event in events.values():
        unique(event['artifacts'], 'artifact in ' + event['id'])
        if not event['command'][0].strip():
            raise ValueError('event command[0] 不能为空')
        if not set(event['input_ids']) <= inputs.keys():
            raise ValueError('event 引用了未登记输入: %s' % event['id'])
        if event['execution_status'] == 'not_run':
            if event['exit_code'] is not None:
                raise ValueError('not_run 不得登记退出码')
        elif event['exit_code'] is None:
            raise ValueError('已执行事件必须登记实际退出码')
        if event['execution_status'] == 'completed' and (not event['artifacts'] or not event['input_ids']):
            raise ValueError('completed 事件需要输入 ID 和真实产物')

    def conclusion_rules(conclusion, historical=False):
        if not historical and case['case_type'] == 'normal_validation_case' and conclusion['status'] == 'UNRESOLVED':
            raise ValueError('normal_validation_case 不得包含未解决结论；请使用 abnormal_case')
        if conclusion['status'] == 'RESOLVED' and not conclusion['evidence_for']:
            raise ValueError('RESOLVED 必须关联支持证据')
        if conclusion['confidence'] == 'high' and not conclusion['evidence_for']:
            raise ValueError('high 置信度不能没有支持证据')
        if conclusion['status'] == 'RESOLVED' and conclusion['confidence'] == 'not_assessable':
            raise ValueError('RESOLVED 与 not_assessable 矛盾')
        if conclusion['reads_support'] != 'NOT_ASSESSED':
            supporting_events = [events[e['event_id']] for e in conclusion['evidence_for'] if e['event_id'] in events]
            if not any(any(inputs[i]['role'] in ('reads', 'bam') for i in event['input_ids'])
                       for event in supporting_events):
                raise ValueError('reads 支持声明需要关联使用 reads/bam 角色输入的支持事件')
        competitors = set(conclusion['competitive_reference_ids'])
        if not competitors <= inputs.keys():
            raise ValueError('竞争参考必须先登记为输入')
        if conclusion['reads_support'] == 'READS_DISCRIMINATING':
            if not competitors:
                raise ValueError('READS_DISCRIMINATING 必须登记竞争参考')
            if any(inputs[i]['role'] != 'competitive_reference' for i in competitors):
                raise ValueError('竞争参考输入须明确标为 competitive_reference')
            support_inputs = set().union(*(set(events[e['event_id']]['input_ids'])
                for e in conclusion['evidence_for'] if e['event_id'] in events))
            if not competitors <= support_inputs:
                raise ValueError('竞争参考未用于所引用的支持事件')
        review = conclusion['review']
        if review['status'] in ('reviewed', 'disputed') and not review['reviewer']:
            raise ValueError('已复核/有分歧必须记录实际复核人')
        for evidence in conclusion['evidence_for'] + conclusion['evidence_against']:
            event = events.get(evidence['event_id'])
            if event is None:
                raise ValueError('证据引用了未知 event: %s' % evidence['event_id'])
            if evidence['artifact_id'] not in {a['id'] for a in event['artifacts']}:
                raise ValueError('证据引用了未知 artifact: %s' % evidence['artifact_id'])
            if event['execution_status'] != 'completed' and conclusion['conclusion_type'] != 'tool_execution':
                raise ValueError('失败/未执行事件不能作为科学支持或反证')

    for candidate in candidates.values():
        base = inputs.get(candidate['base_input_id'])
        derived = inputs.get(candidate['candidate_input_id'])
        producer = events.get(candidate['producer_event_id'])
        if base is None or derived is None or derived['role'] != 'candidate':
            raise ValueError('候选必须关联已登记的基础输入与 candidate 输入')
        if producer is None or producer['execution_status'] != 'completed':
            raise ValueError('候选必须由已完成事件产生')
        artifact = next((a for a in producer['artifacts'] if a['id'] == candidate['artifact_id']), None)
        if artifact is None or (artifact['path'], artifact['sha256']) != (derived['path'], derived['sha256']):
            raise ValueError('候选输入必须与生成事件的指定产物完全一致')
        if candidate['base_input_id'] not in producer['input_ids']:
            raise ValueError('候选生成事件未登记其基础输入')
        if base['path'] == derived['path'] or (base['sha256'] == derived['sha256'] and base['size'] == derived['size']):
            raise ValueError('候选必须与基础输入分离且内容有差异')
        results = []
        for check in candidate['checks']:
            event = events.get(check['event_id'])
            if event is None or event['execution_status'] != 'completed':
                raise ValueError('候选验证必须关联已完成事件')
            if candidate['candidate_input_id'] not in event['input_ids']:
                raise ValueError('候选验证事件必须实际登记候选输入')
            if check['result'] == 'pass' and candidate['change_type'] in ('sequence', 'structure') and not any(
                    inputs[input_id]['role'] in ('reads', 'bam') for input_id in event['input_ids']):
                raise ValueError('碱基/结构候选通过验证必须关联 reads 或 bam 输入')
            results.append(check['result'])
        expected_status = candidate_status(results)
        if candidate['status'] != expected_status:
            raise ValueError('候选状态须与登记的验证结果一致: %s' % expected_status)

    for conclusion in case['conclusions']:
        conclusion_rules(conclusion)
    for entry in case['conclusion_history']:
        if entry['revision'] >= case['revision']:
            raise ValueError('历史修订号必须早于当前版本')
        if entry['conclusion']['id'] not in {c['id'] for c in case['conclusions']}:
            raise ValueError('历史结论缺少对应的当前结论')
        conclusion_rules(entry['conclusion'], historical=True)
    if verify_files:
        checked = {}
        for item in files_in(case):
            path = item['path']
            if path not in checked:
                checked[path] = snapshot(path, item['id'], item['role'])
            actual = checked[path]
            if (item['sha256'], item['size']) != (actual['sha256'], actual['size']):
                raise ValueError('文件内容漂移: %s' % path)
    return case


class CaseStore:
    def __init__(self, directory):
        self.directory = Path(directory).resolve()
        self.path = self.directory / 'case.json'
        self.lockpath = self.directory / '.case.lock'

    @contextmanager
    def locked(self, create=False):
        # Dependency failure must not leave an apparent initialized task.
        validator()
        if create:
            self.directory.mkdir(parents=True, exist_ok=True)
        fd = os.open(self.lockpath, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
        try:
            info = os.fstat(fd)
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise ValueError('任务锁必须为独立普通文件')
            fcntl.flock(fd, fcntl.LOCK_EX)
            if self.path.is_symlink():
                raise ValueError('case.json 不得为符号链接')
            yield
        finally:
            os.close(fd)

    def reserved(self, paths):
        for path in paths:
            protect_output(path, [self.path, self.lockpath])

    def load(self, verify_files=False):
        case = validate(read_json(self.path), verify_files)
        self.reserved(item['path'] for item in files_in(case))
        return case

    def save(self, case):
        validate(case, verify_files=True)
        self.reserved(item['path'] for item in files_in(case))
        atomic_write_text(self.path, json.dumps(case, ensure_ascii=False, indent=2, allow_nan=False)+'\n')

    def init(self, identifier, issue, taxon, case_type, inputs):
        with self.locked(create=True):
            if self.path.exists():
                raise ValueError('案例已存在，拒绝覆盖')
            timestamp = now()
            case = dict(schema_version='3.0', case_id=identifier, issue=issue, taxon=taxon,
                        case_type=case_type, revision=1, created_at=timestamp, updated_at=timestamp,
                        inputs=[snapshot(path, 'I%d' % (i+1), role) for i, (role, path) in enumerate(inputs)],
                        events=[], conclusions=[], conclusion_history=[], candidates=[])
            self.save(case)
            return case

    def mutate(self, action, record, update=False, expected_revision=None):
        with self.locked():
            case = self.load(verify_files=True)
            if expected_revision is not None and expected_revision != case['revision']:
                raise ValueError('案例版本已变化，请重新读取后提交')
            revision = case['revision']
            timestamp = now()
            if action == 'input':
                identifier = 'I%d' % (len(case['inputs'])+1)
                if any(item['id'] == identifier for item in case['inputs']):
                    raise ValueError('输入 ID 序列不连续，拒绝自动分配')
                case['inputs'].append(snapshot(record['path'], identifier, record['role']))
            elif action == 'event':
                event = deepcopy(record)
                if not isinstance(event, dict):
                    raise ValueError('事件必须为 JSON object')
                if 'recorded_at' in event:
                    raise ValueError('recorded_at 由程序写入')
                event['recorded_at'] = timestamp
                # Validate raw artifact shapes before replacing their computed fields.
                artifacts = event.get('artifacts', [])
                if not isinstance(artifacts, list):
                    raise ValueError('artifacts 必须为数组')
                converted = []
                for artifact in artifacts:
                    if not isinstance(artifact, dict) or set(artifact) != {'id', 'role', 'path'}:
                        raise ValueError('artifact 仅接受 id/role/path；hash 由程序计算')
                    if not isinstance(artifact['path'], str) or not Path(artifact['path']).is_absolute():
                        raise ValueError('artifact path 必须为绝对路径')
                    converted.append(snapshot(artifact['path'], artifact['id'], artifact['role']))
                event['artifacts'] = converted
                case['events'].append(event)
            elif action == 'candidate':
                proposal = deepcopy(record)
                required = {'id', 'change_type', 'base_input_id', 'producer_event_id', 'artifact_id',
                            'rationale', 'expected_effect', 'risks'}
                if not isinstance(proposal, dict) or set(proposal) != required:
                    raise ValueError('候选字段须为 ' + ', '.join(sorted(required)))
                if proposal['id'] in {item['id'] for item in case.get('candidates', [])}:
                    raise ValueError('候选 ID 已存在')
                inputs = unique(case['inputs'], 'input')
                events = unique(case['events'], 'event')
                base = inputs.get(proposal['base_input_id'])
                event = events.get(proposal['producer_event_id'])
                if base is None or event is None or event['execution_status'] != 'completed':
                    raise ValueError('候选需要已登记的基础输入和已完成生成事件')
                if proposal['base_input_id'] not in event['input_ids']:
                    raise ValueError('候选生成事件未使用声明的基础输入')
                artifact = next((a for a in event['artifacts'] if a['id'] == proposal['artifact_id']), None)
                if artifact is None:
                    raise ValueError('候选 artifact 不属于指定生成事件')
                protect_output(artifact['path'], [item['path'] for item in case['inputs']])
                candidate_id = 'I%d' % (len(case['inputs']) + 1)
                case['inputs'].append(snapshot(artifact['path'], candidate_id, 'candidate'))
                candidate = dict(proposal, candidate_input_id=candidate_id, created_at=timestamp,
                                 checks=[], status='proposed')
                case.setdefault('candidates', []).append(candidate)
            elif action == 'candidate-check':
                check = deepcopy(record)
                required = {'candidate_id', 'event_id', 'result', 'criteria', 'notes'}
                if not isinstance(check, dict) or set(check) != required:
                    raise ValueError('验证记录字段须为 ' + ', '.join(sorted(required)))
                candidate = next((c for c in case.get('candidates', []) if c['id'] == check['candidate_id']), None)
                event = next((e for e in case['events'] if e['id'] == check['event_id']), None)
                if candidate is None or event is None or event['execution_status'] != 'completed':
                    raise ValueError('验证必须引用已登记候选和已完成事件')
                if candidate['candidate_input_id'] not in event['input_ids']:
                    raise ValueError('验证事件未使用该候选输入')
                if check['result'] == 'pass' and event['exit_code'] != 0:
                    raise ValueError('非零退出码的验证事件不能登记为通过')
                if check['result'] == 'pass' and candidate['change_type'] in ('sequence', 'structure'):
                    inputs = {item['id']: item for item in case['inputs']}
                    if not any(inputs[input_id]['role'] in ('reads', 'bam') for input_id in event['input_ids']):
                        raise ValueError('碱基/结构候选通过验证必须关联 reads 或 bam 输入')
                if any(item['event_id'] == check['event_id'] for item in candidate['checks']):
                    raise ValueError('同一事件不能重复作为候选验证')
                check.pop('candidate_id')
                check['recorded_at'] = timestamp
                candidate['checks'].append(check)
                results = [item['result'] for item in candidate['checks']]
                candidate['status'] = candidate_status(results)
            elif action == 'conclusion':
                conclusion = deepcopy(record)
                if not isinstance(conclusion, dict) or 'recorded_at' in conclusion:
                    raise ValueError('结论须为 object，recorded_at 由程序写入')
                conclusion['recorded_at'] = timestamp
                old = next((c for c in case['conclusions'] if c['id'] == conclusion.get('id')), None)
                if old is not None:
                    if not update:
                        raise ValueError('结论 ID 已存在，显式 --update 才可修订')
                    case['conclusion_history'].append(dict(revision=revision, replaced_at=timestamp,
                                                           conclusion=deepcopy(old)))
                    case['conclusions'][case['conclusions'].index(old)] = conclusion
                elif update:
                    raise ValueError('--update 的结论不存在')
                else:
                    case['conclusions'].append(conclusion)
            else:
                raise ValueError('未知 mutation')
            case['revision'] += 1
            case['updated_at'] = timestamp
            self.save(case)
            return case

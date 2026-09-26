#!/usr/bin/env python3
"""#4：JSON Schema 之外的跨字段业务规则。

为什么单独一层：schema 只表达结构与类型；"案例判 RESOLVED 却挂着 UNRESOLVED 异常"这类关系属于
**业务自洽性**，混进 schema 会把"格式合法"与"科学自洽"混为一谈。因此：

* 每条规则有**稳定错误码**（见 `BUSINESS_RULES`），正例与反例都进测试；
* `case-validate` 把 `FORMAT` 与 `BUSINESS` 分开打印，`VALID` 只表示"格式合法 + 已定义业务规则
  无矛盾"，不隐含科学结论；
* 写入口（`case-anomaly` / `case-reference`）**不因业务矛盾阻断写入**（没有"改 decision"的命令，
  阻断会把人卡死），但必须告警。
"""
import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
EXPERIENCE = ROOT / 'tools' / 'experience.py'


def load_module(name='experience_business'):
    spec = importlib.util.spec_from_file_location(name, EXPERIENCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _Base(unittest.TestCase):
    def setUp(self):
        self.module = load_module(self.__class__.__name__)
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = pathlib.Path(temp.name) / 'case-001'
        self.dir.mkdir(parents=True)
        (self.dir / 'events.jsonl').write_text('', encoding='utf-8')

    def write_case(self, **overrides):
        case = {
            'schema_version': '2.0', 'case_id': 'case-001', 'case_type': 'abnormal_case',
            'inputs': [], 'issue': {'type': 'internal_stop', 'user_observation': 'synthetic'},
            'hypotheses': [{'id': 'H1', 'explanation': 'e', 'support': [], 'against': [],
                            'unknown': []}],
            'events_file': 'events.jsonl', 'anomalies': [],
            'decision': {'status': 'UNRESOLVED', 'confidence': 'low', 'rationale': 'r'},
        }
        case.update(overrides)
        (self.dir / 'case.json').write_text(json.dumps(case, ensure_ascii=False, indent=2) + '\n',
                                            encoding='utf-8')
        return case

    def business(self, case, **kwargs):
        return self.module.case_business_errors(case, self.dir, **kwargs)

    def codes(self, findings):
        return sorted(f['code'] for f in findings)


class RuleTests(_Base):
    def test_every_rule_has_a_stable_code_and_description(self):
        rules = self.module.BUSINESS_RULES
        self.assertTrue(rules)
        for code, description in rules.items():
            with self.subTest(code=code):
                self.assertRegex(code, r'^[A-Z][A-Z0-9_]+$')
                self.assertTrue(description.strip())

    # 1) 案例级 RESOLVED 与异常级 UNRESOLVED 矛盾
    def test_resolved_case_with_unresolved_anomaly_is_a_conflict(self):
        case = self.write_case(decision={'status': 'RESOLVED', 'confidence': 'moderate',
                                         'rationale': 'ok'},
                               anomalies=[{'id': 'A1', 'claim': 'x', 'status': 'UNRESOLVED',
                                           'confidence': 'low'}])
        self.assertIn('DECISION_ANOMALY_CONFLICT', self.codes(self.business(case)))

    def test_resolved_case_with_only_resolved_anomalies_is_clean(self):
        case = self.write_case(decision={'status': 'RESOLVED', 'confidence': 'moderate',
                                         'rationale': 'ok'},
                               anomalies=[{'id': 'A1', 'claim': 'x', 'status': 'RESOLVED',
                                           'confidence': 'moderate'}])
        self.assertEqual(self.business(case), [])

    def test_unresolved_case_with_unresolved_anomaly_is_clean(self):
        case = self.write_case(anomalies=[{'id': 'A1', 'claim': 'x', 'status': 'UNRESOLVED',
                                           'confidence': 'low'}])
        self.assertEqual(self.business(case), [])

    # 2) 异常 ID 唯一性
    def test_duplicate_anomaly_ids_are_flagged(self):
        case = self.write_case(anomalies=[
            {'id': 'A1', 'claim': 'x', 'status': 'UNRESOLVED', 'confidence': 'low'},
            {'id': 'A1', 'claim': 'y', 'status': 'RESOLVED', 'confidence': 'moderate'}])
        self.assertIn('ANOMALY_ID_DUPLICATE', self.codes(self.business(case)))

    # 3) 正常验证案例不得携带未解决异常
    def test_normal_validation_case_cannot_carry_unresolved_anomalies(self):
        case = self.write_case(case_type='normal_validation_case',
                               anomalies=[{'id': 'A1', 'claim': 'x', 'status': 'UNRESOLVED',
                                           'confidence': 'low'}])
        self.assertIn('NORMAL_CASE_HAS_UNRESOLVED_ANOMALY', self.codes(self.business(case)))
        ok = self.write_case(case_type='normal_validation_case', anomalies=[
            {'id': 'A1', 'claim': 'checked', 'status': 'RESOLVED', 'confidence': 'moderate'}])
        self.assertEqual(self.business(ok), [])

    # 4) 异常引用的证据事件必须真实存在
    def test_anomaly_referencing_an_unknown_event_is_flagged(self):
        (self.dir / 'events.jsonl').write_text(
            json.dumps({'action': 'annot_check', 'result': 'x', 'impact': 'H1:against'}) + '\n',
            encoding='utf-8')
        bad = self.write_case(anomalies=[{'id': 'A1', 'claim': 'x', 'status': 'UNRESOLVED',
                                          'confidence': 'low',
                                          'evidence_events': ['depth_analysis']}])
        self.assertIn('ANOMALY_EVENT_UNKNOWN', self.codes(self.business(bad)))
        good = self.write_case(anomalies=[{'id': 'A1', 'claim': 'x', 'status': 'UNRESOLVED',
                                           'confidence': 'low',
                                           'evidence_events': ['annot_check']}])
        self.assertEqual(self.business(good), [])

    # 5) 事件 impact 引用的假设必须存在
    def test_event_impact_referencing_an_unknown_hypothesis_is_flagged(self):
        (self.dir / 'events.jsonl').write_text(
            json.dumps({'action': 'annot_check', 'result': 'x', 'impact': 'H9:against'}) + '\n',
            encoding='utf-8')
        case = self.write_case()
        self.assertIn('EVENT_HYPOTHESIS_UNKNOWN', self.codes(self.business(case)))
        (self.dir / 'events.jsonl').write_text(
            json.dumps({'action': 'annot_check', 'result': 'x', 'impact': 'H1:against'}) + '\n',
            encoding='utf-8')
        self.assertEqual(self.business(case), [])

    # 6) 输入文件必须存在；哈希只在显式要求时比对
    def test_missing_input_file_is_flagged(self):
        case = self.write_case(inputs=[{'role': 'assembly_fasta',
                                        'path': str(self.dir / 'nope.fa'),
                                        'sha256': 'a' * 64}])
        self.assertIn('INPUT_FILE_MISSING', self.codes(self.business(case)))

    def test_sha256_mismatch_is_only_checked_when_asked(self):
        target = self.dir / 'asm.fa'
        target.write_text('>x\nACGT\n', encoding='utf-8')
        case = self.write_case(inputs=[{'role': 'assembly_fasta', 'path': str(target),
                                        'sha256': 'b' * 64}])
        self.assertEqual(self.codes(self.business(case)), [],
                         '默认不做重哈希（大 FASTQ 代价高）；只校验存在性')
        self.assertIn('INPUT_SHA256_MISMATCH',
                      self.codes(self.business(case, verify_inputs=True)))

    def test_correct_sha256_passes_when_verified(self):
        target = self.dir / 'asm.fa'
        target.write_text('>x\nACGT\n', encoding='utf-8')
        digest = self.module.file_sha256(target)
        case = self.write_case(inputs=[{'role': 'assembly_fasta', 'path': str(target),
                                        'sha256': digest}])
        self.assertEqual(self.business(case, verify_inputs=True), [])

    # 7) events.jsonl 必须可解析
    def test_unparsable_events_line_is_flagged(self):
        (self.dir / 'events.jsonl').write_text('not json at all\n', encoding='utf-8')
        self.assertIn('EVENTS_UNPARSABLE', self.codes(self.business(self.write_case())))


class CliOutputTests(_Base):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(EXPERIENCE), 'case-validate', str(self.dir), *args],
                              capture_output=True, text=True, timeout=180)

    def test_clean_case_is_valid_and_separates_the_two_layers(self):
        self.write_case(anomalies=[{'id': 'A1', 'claim': 'x', 'status': 'UNRESOLVED',
                                    'confidence': 'low'}])
        result = self.run_cli()
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, blob)
        self.assertIn('FORMAT: OK', blob)
        self.assertIn('BUSINESS: 无矛盾', blob)
        self.assertIn('VALID', blob)

    def test_business_contradiction_is_invalid_and_named(self):
        self.write_case(decision={'status': 'RESOLVED', 'confidence': 'moderate', 'rationale': 'ok'},
                        anomalies=[{'id': 'A1', 'claim': 'x', 'status': 'UNRESOLVED',
                                    'confidence': 'low'}])
        result = self.run_cli()
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 1, blob)
        self.assertIn('FORMAT: OK', blob, '格式没问题，不能报成格式错误')
        self.assertIn('DECISION_ANOMALY_CONFLICT/error', blob)
        self.assertIn('cause: business', blob)
        self.assertFalse(any(line.strip() == 'VALID' for line in blob.splitlines()),
                         '有业务矛盾时不得给出 VALID')

    def test_verify_inputs_flag_changes_the_verdict(self):
        target = self.dir / 'asm.fa'
        target.write_text('>x\nACGT\n', encoding='utf-8')
        self.write_case(inputs=[{'role': 'assembly_fasta', 'path': str(target),
                                 'sha256': 'c' * 64}])
        self.assertEqual(self.run_cli().returncode, 0)
        strict = self.run_cli('--verify-inputs')
        self.assertEqual(strict.returncode, 1, strict.stdout + strict.stderr)
        self.assertIn('INPUT_SHA256_MISMATCH', strict.stdout + strict.stderr)


class WritePathTests(_Base):
    """业务矛盾不阻断写入（没有改 decision 的命令），但必须告警。"""

    def write_via_cli(self, *args):
        return subprocess.run([sys.executable, str(EXPERIENCE), *args],
                              capture_output=True, text=True, timeout=180)

    def test_adding_an_unresolved_anomaly_warns_but_writes(self):
        self.write_case(decision={'status': 'RESOLVED', 'confidence': 'moderate', 'rationale': 'ok'})
        result = self.write_via_cli('case-anomaly', str(self.dir), '--id', 'A1', '--claim', 'x',
                                    '--status', 'UNRESOLVED', '--confidence', 'low')
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, blob)
        self.assertIn('DECISION_ANOMALY_CONFLICT', blob, '必须显式告警')
        self.assertIn('⚠', blob)
        data = json.loads((self.dir / 'case.json').read_text(encoding='utf-8'))
        self.assertEqual([a['id'] for a in data['anomalies']], ['A1'],
                         '告警不得变成静默阻断（没有改 decision 的命令，阻断会卡死用户）')

    def test_format_error_still_blocks_the_write(self):
        self.write_case()
        data = json.loads((self.dir / 'case.json').read_text(encoding='utf-8'))
        data['anomalies'] = 'not-an-array'
        (self.dir / 'case.json').write_text(json.dumps(data, ensure_ascii=False), encoding='utf-8')
        before = (self.dir / 'case.json').read_text(encoding='utf-8')
        result = self.write_via_cli('case-anomaly', str(self.dir), '--id', 'A1', '--claim', 'x',
                                    '--status', 'UNRESOLVED', '--confidence', 'low')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((self.dir / 'case.json').read_text(encoding='utf-8'), before,
                         '格式错误仍然必须拒绝写入')


if __name__ == '__main__':
    unittest.main()

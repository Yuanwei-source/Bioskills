#!/usr/bin/env python3
"""缺少必需依赖时必须响亮失败，而不是静默降级。

事故形状（要消除的）：`case-validate` 曾有两份实现——`jsonschema` 读 schema，以及
`tools/experience.py` 里一份 87 行手写"等价"实现。没有库时静默切到第二份，
两者都会对同一案例说 "VALID"，却需要人肉同步（已偏过两次）。同一结论必须对应同一证据路径，
因此那份实现被删除；缺库时 `case-validate` 必须**拒绝工作**并给出安装命令。

这些用例用自带的 blocker（临时目录里放一个 `raise ImportError` 的 `jsonschema.py`）
在**无 jsonschema 环境**下运行真实 CLI。
"""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPERIENCE = ROOT / 'tools' / 'experience.py'

SYNTHETIC_CASE = {
    "schema_version": "2.0", "case_id": "case-001", "case_type": "abnormal_case",
    "inputs": [], "issue": {"type": "internal_stop", "user_observation": "synthetic"},
    "hypotheses": [{"id": "H1", "explanation": "boundary error", "support": [], "against": [],
                    "unknown": []}],
    "events_file": "events.jsonl", "anomalies": [],
    "decision": {"status": "UNRESOLVED", "confidence": "low", "rationale": "synthetic"},
}


def load_module(name='experience_missing_dep'):
    spec = importlib.util.spec_from_file_location(name, EXPERIENCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MissingDependencyTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        # 伪造"没装 jsonschema"的环境
        self.blocker = self.root / 'blocker'
        self.blocker.mkdir()
        (self.blocker / 'jsonschema.py').write_text(
            'raise ImportError("jsonschema intentionally absent")\n', encoding='utf-8')
        self.env = dict(os.environ)
        self.env['PYTHONPATH'] = str(self.blocker)
        self.case = self.root / 'case-001'
        self.case.mkdir()
        (self.case / 'case.json').write_text(json.dumps(SYNTHETIC_CASE, ensure_ascii=False, indent=2),
                                             encoding='utf-8')
        (self.case / 'events.jsonl').touch()

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(EXPERIENCE), *[str(a) for a in args]],
                              capture_output=True, text=True, timeout=180, env=self.env)

    def case_text(self):
        return (self.case / 'case.json').read_text(encoding='utf-8')

    # --- 1) 校验入口：拒绝工作并给出安装命令 -------------------------------
    def test_case_validate_without_jsonschema_exits_3(self):
        result = self.run_cli('case-validate', self.case)
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 3,
                         '缺少必需依赖必须以 3 退出（环境/依赖故障），实际: %s\n%s'
                         % (result.returncode, blob))
        self.assertIn('jsonschema', blob)
        self.assertIn('pip install', blob, '必须给出可执行的安装命令')
        self.assertNotIn('VALID', blob, '不得在缺依赖时给出任何裁定')

    # --- 2) 写入口：不得在未校验的情况下落盘 -------------------------------
    def test_case_anomaly_without_jsonschema_refuses_to_write(self):
        before = self.case_text()
        result = self.run_cli('case-anomaly', self.case, '--id', 'A1', '--claim', 'x',
                              '--status', 'UNRESOLVED', '--confidence', 'low')
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertEqual(self.case_text(), before, '校验不可用时不得写入案例')

    def test_case_reference_without_jsonschema_refuses_to_write(self):
        before = self.case_text()
        result = self.run_cli('case-reference', self.case, '--reference-id', 'ref-001',
                              '--purpose', 'gene_order_comparison', '--registry',
                              str(self.root / 'nope'))
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertEqual(self.case_text(), before)

    # --- 3) 双实现必须真的被删掉 -------------------------------------------
    def test_duplicate_validator_is_gone(self):
        module = load_module('experience_duplicate_check')
        self.assertFalse(hasattr(module, '_case_errors_without_jsonschema'),
                         '第二份案例校验实现必须删除：同一结论只能有一条证据路径')
        self.assertTrue(hasattr(module, 'MissingDependency'),
                        '需要有一个可区分的"缺依赖"异常，供 CLI 以退出码 3 处理')


if __name__ == '__main__':
    unittest.main()

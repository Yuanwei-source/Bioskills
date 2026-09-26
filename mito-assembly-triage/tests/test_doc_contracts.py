#!/usr/bin/env python3
"""代码 ↔ 文档一致性护栏的测试。

背景（真实事故）：一次文档压缩把 6 个 CLI flag 的**唯一**文档位置删掉了，其中 4 个是
`required=True`。单元测试全绿，因为代码没变——漂移发生在"文档侧"。
`tools/audit_code_docs.py` 就是为这一类失败加的常驻检查；本文件既验证它在当前仓库上是干净的，
也用**合成漂移**验证它真的会拦下来（否则这个护栏本身不可信）。
"""
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / 'tools' / 'audit_code_docs.py'


def load_audit():
    spec = importlib.util.spec_from_file_location('audit_code_docs', AUDIT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _Fixture(unittest.TestCase):
    """把被审计的最小文件集复制到临时根目录，便于合成漂移。"""

    def setUp(self):
        self.module = load_audit()
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name) / 'skill'
        (self.root / 'references').mkdir(parents=True)
        (self.root / 'tools').mkdir(parents=True)
        for rel in self.module.CODE_FILES:
            src = ROOT / rel
            if src.is_file():
                dst = self.root / rel
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
        shutil.copy2(ROOT / 'SKILL.md', self.root / 'SKILL.md')
        for doc in sorted((ROOT / 'references').glob('*.md')):
            shutil.copy2(doc, self.root / 'references' / doc.name)
        shutil.copy2(ROOT / 'tools' / 'doc_contract_baseline.json',
                     self.root / 'tools' / 'doc_contract_baseline.json')

    def audit(self):
        baseline = self.module.load_baseline(self.root)
        return self.module.audit(self.root, baseline)

    def kinds(self, result):
        return {kind for kind, _, _ in result['problems']}

    def strip_from_docs(self, token):
        """从**所有**文档里删掉某个 token（模拟压缩时把它冲掉）。"""
        touched = 0
        for path in [self.root / 'SKILL.md'] + sorted((self.root / 'references').glob('*.md')):
            text = path.read_text(encoding='utf-8')
            if token in text:
                path.write_text(text.replace(token, 'TOKEN_REMOVED'), encoding='utf-8')
                touched += 1
        if not touched:
            self.fail('文档里本来就找不到 %s' % token)
        return touched

    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(AUDIT), '--root', str(self.root), *args],
                              capture_output=True, text=True, timeout=180)


class RepoIsCleanTests(unittest.TestCase):
    """当前仓库必须干净：否则这个护栏会淹没在误报里。"""

    def test_audit_reports_no_problems(self):
        module = load_audit()
        result = module.audit(ROOT, module.load_baseline(ROOT))
        self.assertEqual(result['problems'], [],
                         '代码↔文档审计发现漂移：%s' % result['problems'][:5])

    def test_cli_exits_zero(self):
        proc = subprocess.run([sys.executable, str(AUDIT)], capture_output=True, text=True, timeout=180)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)

    def test_baseline_matches_current_state(self):
        proc = subprocess.run([sys.executable, str(AUDIT), '--check-baseline'],
                              capture_output=True, text=True, timeout=180)
        self.assertEqual(proc.returncode, 0,
                         '基线文件与当前状态不一致，需 --update-baseline 明确更新:\n' + proc.stdout)

    def test_interface_flags_are_on_the_protected_list(self):
        """接口语义类 flag 必须"有文档"，不能被豁免掉。"""
        baseline = json.loads((ROOT / 'tools' / 'doc_contract_baseline.json').read_text(encoding='utf-8'))
        documented = set(baseline['documented'])
        for token in ('--input', '--next-test', '--lesson-id', '--reason', '--manifest', '--query',
                      '--reference-id', '--event-action', '--case-type', '--lesson-domain',
                      '--authorize', '--allow-public-upload'):
            self.assertIn(token, documented, '%s 必须保持有文档' % token)

    def test_json_output_is_machine_readable(self):
        proc = subprocess.run([sys.executable, str(AUDIT), '--json'],
                              capture_output=True, text=True, timeout=180)
        payload = json.loads(proc.stdout)
        self.assertIn('problems', payload)
        self.assertEqual(payload['problems'], [])


class DriftDetectionTests(_Fixture):
    """合成漂移必须被拦下——这就是那次真实事故的形状。"""

    def test_documentation_loss_is_detected(self):
        self.strip_from_docs('--next-test')
        result = self.audit()
        self.assertIn('documentation-lost', self.kinds(result),
                      'documented 列表里的 token 消失必须报 documentation-lost')

    def test_new_undocumented_token_is_detected(self):
        """新增一个没写文档的 flag：审计必须要求"写文档或显式豁免"。"""
        target = self.root / 'scripts' / 'annot_check.py'
        text = target.read_text(encoding='utf-8')
        target.write_text(text + "\n# ap.add_argument('--brand-new-flag')\n", encoding='utf-8')
        result = self.audit()
        self.assertIn('new-token-undocumented', self.kinds(result))
        subjects = {s for k, s, _ in result['problems'] if k == 'new-token-undocumented'}
        self.assertIn('--brand-new-flag', subjects)

    def test_token_removed_from_code_is_detected(self):
        target = self.root / 'tools' / 'experience.py'
        text = target.read_text(encoding='utf-8')
        self.assertIn('--input', text)
        target.write_text(text.replace('--input', '--input-renamed'), encoding='utf-8')
        result = self.audit()
        self.assertIn('token-removed-from-code', self.kinds(result))

    def test_rule_evidence_must_exist(self):
        target = self.root / 'scripts' / 'annot_check.py'
        text = target.read_text(encoding='utf-8')
        self.assertIn('CIRCULAR_DECLARATION_CHECK', text)
        target.write_text(text.replace('CIRCULAR_DECLARATION_CHECK', 'SOMETHING_ELSE'), encoding='utf-8')
        result = self.audit()
        self.assertIn('rule-evidence-missing', self.kinds(result))

    def test_rule_must_stay_documented(self):
        for token in ('CIRCULAR_DECLARATION_CHECK', 'TRANSL_EXCEPT_MATCHED', 'GENE_SET_DIFF'):
            with self.subTest(token=token):
                self.strip_from_docs(token)
        result = self.audit()
        self.assertIn('rule-undocumented', self.kinds(result))

    def test_clean_fixture_has_no_problems(self):
        """先证明夹具本身是干净的，否则上面的断言没有意义。"""
        self.assertEqual(self.audit()['problems'], [])

    def test_cli_fails_on_drift(self):
        self.strip_from_docs('--next-test')
        proc = self.run_cli()
        self.assertNotEqual(proc.returncode, 0, '漂移时 CLI 必须以非 0 退出')
        self.assertIn('documentation-lost', proc.stdout + proc.stderr)


if __name__ == '__main__':
    unittest.main()

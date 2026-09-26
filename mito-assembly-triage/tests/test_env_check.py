#!/usr/bin/env python3
"""环境就绪检查的测试：两模式 + 依赖清单 + lock 语义。

最重要的不变式：**lock 是缓存，不是信任凭证**。日常模式必须真实探测本次需要的工具；
若只信 lock，就退化成"第一次通过、以后永远相信"——这类静默降级正是本 skill 反复清除的缺陷。
"""
import importlib.util
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENV_CHECK = ROOT / 'tools' / 'env_check.py'
MANIFEST = ROOT / 'config' / 'dependencies.json'
TOOL_CHECK = ROOT / 'references' / 'tool_check.md'
ESSENTIAL_COMMANDS = ('blastn', 'makeblastdb', 'samtools', 'minimap2')


def load_module(name='env_check'):
    spec = importlib.util.spec_from_file_location(name, ENV_CHECK)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ManifestIntegrityTests(unittest.TestCase):
    """清单是依赖的唯一来源，它自身必须完整且与文档一致。"""

    def setUp(self):
        self.manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
        self.ids = {item['id'] for key in ('python_libs', 'tools')
                    for item in self.manifest.get(key) or []}

    def test_format_and_required_fields(self):
        self.assertEqual(self.manifest['format'], 'mito-dependency-manifest-1')
        for key in ('python_libs', 'tools'):
            for item in self.manifest[key]:
                with self.subTest(dep=item.get('id')):
                    for field in ('id', 'tier', 'probe', 'purpose', 'install'):
                        self.assertIn(field, item)
                    self.assertIn(item['tier'], ('essential', 'extended', 'optional'))
                    self.assertIn(item['probe'].get('type'),
                                  ('command', 'python-module', 'python-import', 'python-interpreter', 'env-dir', 'network'))

    def test_stages_reference_known_dependencies(self):
        for stage, spec in self.manifest['stages'].items():
            with self.subTest(stage=stage):
                self.assertTrue(spec.get('requires'), 'stage 必须有 requires')
                self.assertTrue(spec.get('capability'))
                unknown = sorted(set(spec['requires']) - self.ids)
                self.assertEqual(unknown, [], 'stage %s 引用了清单里不存在的依赖' % stage)

    def test_only_three_tiers_and_all_used(self):
        tiers = {item['tier'] for key in ('python_libs', 'tools') for item in self.manifest[key]}
        self.assertEqual(tiers, {'essential', 'extended', 'optional'})

    def test_check_env_tool_list_is_a_subset_of_the_manifest(self):
        """`check_env.sh` 里的工具名不得出现清单外的条目——防止第二份清单悄悄漂移。"""
        text = (ROOT / 'scripts' / 'check_env.sh').read_text(encoding='utf-8')
        match = re.search(r'for t in ([^\n;]+);', text)
        self.assertIsNotNone(match, 'check_env.sh 的工具循环未找到')
        unknown = [name for name in match.group(1).split() if name not in self.ids]
        self.assertEqual(unknown, [], 'check_env.sh 列了清单外的工具: %s' % unknown)

    def test_every_manifest_dependency_is_documented(self):
        """清单里的每个依赖都必须能在文档里查到（单一来源不能脱离文档）。"""
        text = TOOL_CHECK.read_text(encoding='utf-8')
        for dep_id in sorted(self.ids):
            with self.subTest(dep=dep_id):
                self.assertIn(dep_id, text, '%s 未写进 tool_check.md' % dep_id)


class _Stubs(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        self.lock = self.root / 'environment.lock.json'
        for name in ESSENTIAL_COMMANDS:
            self.make_stub(name)
        self.env = {k: v for k, v in os.environ.items()
                    if k not in ('MINIMAP2', 'MITOS2_PY', 'MITOS2_REFDIR', 'MITOS2_REFSEQVER',
                                 'MITOS2_EXTRA_PATH', 'PLOT_PY')}
        self.env['MITO_KNOWLEDGE_DIR'] = str(self.root / 'knowledge')

    def make_stub(self, name, version='1.2.3'):
        path = self.bin / name
        path.write_text('#!/bin/sh\necho "%s %s"\n' % (name, version), encoding='utf-8')
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        return path

    def run_env_check(self, *args, env=None):
        return subprocess.run([sys.executable, str(ENV_CHECK), '--path', str(self.bin), *args],
                              capture_output=True, text=True, timeout=180,
                              env=env or self.env)


class StageGateTests(_Stubs):
    def test_ready_stage_exits_zero(self):
        result = self.run_env_check('--stage', 'gene_locating')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('gene_locating', result.stdout)

    def test_missing_dependency_exits_three_with_install_hint(self):
        (self.bin / 'blastn').unlink()
        result = self.run_env_check('--stage', 'gene_locating')
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 3, blob)
        self.assertIn('blastn', blob)
        self.assertIn('安装', blob)
        self.assertIn('本步骤未执行', blob)

    def test_unknown_stage_is_a_usage_error(self):
        result = self.run_env_check('--stage', 'no_such_stage')
        self.assertEqual(result.returncode, 1)
        self.assertIn('未知 stage', result.stdout + result.stderr)

    def test_stage_requires_match_manifest(self):
        module = load_module('env_check_stage')
        manifest = module.load_manifest()
        deps = module.dependencies_for(manifest, manifest['stages']['read_evidence']['requires'])
        self.assertEqual(sorted(d[0] for d in deps), ['samtools'])


class SetupModeTests(_Stubs):
    def test_missing_essential_do_not_write_a_lock(self):
        (self.bin / 'samtools').unlink()
        result = self.run_env_check('--setup', '--lock', str(self.lock))
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, blob)
        self.assertFalse(self.lock.exists(), '核心不全时不得写 lock')
        self.assertIn('未就绪', blob)

    def test_complete_essentials_write_a_lock(self):
        result = self.run_env_check('--setup', '--lock', str(self.lock))
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 0, blob)
        self.assertTrue(self.lock.exists())
        data = json.loads(self.lock.read_text(encoding='utf-8'))
        self.assertEqual(data['format'], 'mito-environment-lock-1')
        self.assertTrue(data['checked_at'])
        self.assertIn('samtools', data['available'])
        self.assertNotIn('samtools', data.get('missing') or {})

    def test_extended_missing_does_not_block_without_strict(self):
        result = self.run_env_check('--setup', '--lock', str(self.lock))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        blob = result.stdout
        self.assertIn('增强（extended）', blob)
        self.assertTrue(self.lock.exists())

    def test_strict_requires_extended(self):
        result = self.run_env_check('--setup', '--strict', '--lock', str(self.lock))
        self.assertEqual(result.returncode, 2, '缺 extended 时 --strict 应不就绪')
        self.assertFalse(self.lock.exists())


class DailyModeTests(_Stubs):
    """lock 只是缓存：环境变了必须被发现。"""

    def test_without_lock_reports_status_and_suggests_setup(self):
        result = self.run_env_check('--daily', '--lock', str(self.lock))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('尚无 lock', result.stdout)
        self.assertFalse(self.lock.exists(), '日常模式不得偷偷写 lock')

    def test_ready_environment_uses_the_lock(self):
        self.run_env_check('--setup', '--lock', str(self.lock))
        result = self.run_env_check('--daily', '--lock', str(self.lock))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('就绪', result.stdout)

    def test_tool_disappearing_after_setup_is_detected(self):
        """核心不变式：删掉一个工具后，日常模式必须报错，而不是相信 lock。"""
        self.run_env_check('--setup', '--lock', str(self.lock))
        (self.bin / 'samtools').unlink()
        result = self.run_env_check('--daily', '--lock', str(self.lock))
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 2, blob)
        self.assertIn('环境已变化', blob)
        self.assertIn('samtools', blob)

    def test_version_drift_is_reported_without_blocking(self):
        self.run_env_check('--setup', '--lock', str(self.lock))
        self.make_stub('samtools', version='9.9.9')
        result = self.run_env_check('--daily', '--lock', str(self.lock))
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('版本与 lock 不同', result.stdout)

    def test_json_output_is_machine_readable(self):
        self.run_env_check('--setup', '--lock', str(self.lock))
        result = self.run_env_check('--daily', '--json', '--lock', str(self.lock))
        payload = json.loads(result.stdout)
        self.assertTrue(payload['ready'])
        self.assertEqual(payload['missing'], [])


class LockLocationTests(_Stubs):
    def test_default_lock_lives_in_the_knowledge_dir(self):
        expected = Path(self.env['MITO_KNOWLEDGE_DIR']) / 'environment.lock.json'
        result = self.run_env_check('--setup')       # 不传 --lock
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(expected.is_file(), '默认 lock 必须落在 $MITO_KNOWLEDGE_DIR')
        self.assertNotIn(str(ROOT), str(expected), 'lock 不得写进 skill 安装目录')

    def test_broken_manifest_is_a_usage_error(self):
        bad = self.root / 'bad.json'
        bad.write_text('{"format": "nope"}', encoding='utf-8')
        result = subprocess.run([sys.executable, str(ENV_CHECK), '--manifest', str(bad)],
                                capture_output=True, text=True, timeout=60, env=self.env)
        self.assertEqual(result.returncode, 1)
        self.assertIn('清单', result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()

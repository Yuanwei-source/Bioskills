#!/usr/bin/env python3
"""脚本前置门禁（issue #15 的三件事）。

1. `annot_check.py` 不再有 "缺 Biopython 就 return None" 的静默降级；
2. 缺依赖的失败契约统一为**退出码 3 + 用途 + 安装方式 + "该步骤未执行"**（不再有
   静默降级 / 退出 1 / 裸回溯三种行为）；
3. 每个脚本按 stage 声明自己的边界，依赖清单只有一份（config/dependencies.json）。

以及一条容易忘的边界：**只阻断该步骤**——缺 samtools 不该拦住只读 FASTA 的检查。
"""
import importlib.util
import os
import pathlib
import re
import stat
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts'

# 脚本 → 它应当声明的 stage
STAGE_OF = {
    'seq_stats.py': 'seq_stats',
    'annot_check.py': 'annot_check',
    'blast_genes.py': 'gene_locating',
    'circularize.py': 'circularize',
    'depth_analysis.py': 'read_evidence',
    'mitos2_to_genbank.py': 'mitos2_bridge',
}
# 刻意不按 stage 门禁的脚本：它们的网络类依赖只在具体动作里需要（且已自有退出码契约）
UNGATED = {'cox1_id.py': '网络依赖只在查询时需要；失败已用退出码 3 表达',
           'reference_registry.py': 'register/list/verify 可离线工作，只有 acquire 需要网络',
           '_deps.py': '门禁实现自身'}

TINY_FASTA = '>seq1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n'


class _Env(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = pathlib.Path(temp.name)
        self.base_env = {k: v for k, v in os.environ.items()
                         if k not in ('MINIMAP2', 'MITOS2_PY', 'MITOS2_REFDIR',
                                      'MITOS2_REFSEQVER', 'MITOS2_EXTRA_PATH', 'PLOT_PY')}
        self.fasta = self.root / 'tiny.fasta'
        self.fasta.write_text(TINY_FASTA, encoding='utf-8')

    # 造一个"缺 Biopython"的环境
    def env_without_biopython(self):
        blocker = self.root / 'nobio'
        blocker.mkdir(exist_ok=True)
        (blocker / 'Bio.py').write_text('raise ImportError("Bio blocked for test")\n',
                                        encoding='utf-8')
        env = dict(self.base_env)
        env['PYTHONPATH'] = str(blocker)
        return env

    # 造一个"缺外部命令"的环境（只剩 python3）
    def env_with_minimal_path(self):
        bindir = self.root / 'minbin'
        bindir.mkdir(exist_ok=True)
        link = bindir / 'python3'
        if not link.exists():
            link.symlink_to(sys.executable)
        env = dict(self.base_env)
        env['PATH'] = str(bindir)
        return env

    def run_script(self, name, *args, env=None):
        return subprocess.run([sys.executable, str(SCRIPTS / name), *[str(a) for a in args]],
                              capture_output=True, text=True, timeout=180, env=env or self.base_env)


class DeclarationTests(_Env):
    def test_every_script_declares_its_stage(self):
        for name, stage in STAGE_OF.items():
            with self.subTest(script=name):
                text = (SCRIPTS / name).read_text(encoding='utf-8')
                self.assertIn("require_stage('%s'" % stage, text,
                              '%s 未声明 stage=%s' % (name, stage))

    def test_ungated_scripts_have_a_documented_reason(self):
        for name, reason in UNGATED.items():
            with self.subTest(script=name):
                self.assertTrue(reason, name)

    def test_no_script_has_a_silent_import_fallback(self):
        """杀掉那类"缺库就换行为"的写法（annot_check 曾 return None）。"""
        offenders = []
        for path in sorted(SCRIPTS.glob('*.py')):
            text = path.read_text(encoding='utf-8')
            if 'except ImportError' in text:
                offenders.append(path.name)
        self.assertEqual(offenders, [], '仍存在 except ImportError 降级路径: %s' % offenders)

    def test_shell_scripts_gate_on_their_stage(self):
        for name, stage in (('run_mitos2.sh', 'annot_independent'),
                            ('run_circular_map.sh', 'circular_plot')):
            with self.subTest(script=name):
                text = (SCRIPTS / name).read_text(encoding='utf-8')
                self.assertIn('env_check.py" --stage %s' % stage, text)


class MissingDependencyContractTests(_Env):
    def _assert_gate_output(self, result, script):
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 3,
                         '%s 缺依赖应以 3 退出，实际 %s\n%s' % (script, result.returncode, blob))
        self.assertNotIn('Traceback', blob, '%s 不得抛回溯' % script)
        self.assertIn('未执行', blob)
        self.assertIn('安装', blob)

    def test_biopython_missing_blocks_biopython_scripts(self):
        env = self.env_without_biopython()
        for name in ('seq_stats.py', 'annot_check.py', 'mitos2_to_genbank.py',
                     'blast_genes.py', 'circularize.py'):
            with self.subTest(script=name):
                self._assert_gate_output(self.run_script(name, self.fasta, self.fasta, env=env), name)

    def test_annot_check_reports_no_verdict_without_biopython(self):
        result = self.run_script('annot_check.py', self.fasta, env=self.env_without_biopython())
        blob = result.stdout + result.stderr
        self.assertEqual(result.returncode, 3)
        for verdict in ('VALID', 'INVALID', 'ERROR'):
            self.assertNotIn(verdict, blob, '缺依赖时不得给出任何注释结论')

    def test_missing_samtools_blocks_only_the_read_step(self):
        env = self.env_with_minimal_path()
        blocked = self.run_script('depth_analysis.py', self.fasta, self.fasta, env=env)
        self.assertEqual(blocked.returncode, 3, blocked.stdout + blocked.stderr)
        self.assertIn('samtools', blocked.stdout + blocked.stderr)
        # 同一环境下，只依赖 Biopython 的步骤仍应正常工作（按阶段阻断，不搞全局门）
        allowed = self.run_script('seq_stats.py', self.fasta, env=env)
        self.assertEqual(allowed.returncode, 0, allowed.stdout + allowed.stderr)
        self.assertIn('seq1', allowed.stdout)

    def test_help_works_when_dependencies_are_present(self):
        result = self.run_script('annot_check.py', '--help')        # 默认环境依赖齐全
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('注释质检', result.stdout + result.stderr)

    def test_help_is_not_a_loophole_around_the_gate(self):
        """模块级 import 会先执行，所以 help 也不该在缺依赖时假装可用。"""
        result = self.run_script('annot_check.py', '--help', env=self.env_without_biopython())
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertIn('未执行', result.stdout + result.stderr)

    def test_unknown_stage_fails_closed_without_traceback(self):
        code = ("import sys, os; sys.path.insert(0, %r);"
                "from _deps import require_stage;"
                "require_stage('no_such_stage', 'x.py')" % str(SCRIPTS))
        result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True,
                                timeout=120, env=self.base_env)
        self.assertEqual(result.returncode, 3)
        self.assertNotIn('Traceback', result.stdout + result.stderr)
        self.assertIn('依赖清单', result.stdout + result.stderr)

    def test_gate_uses_the_manifest_and_does_not_reimplement_it(self):
        """门禁必须从 env_check 走；脚本里不应出现自己写的依赖列表。"""
        text = (SCRIPTS / '_deps.py').read_text(encoding='utf-8')
        self.assertIn('env_check.py', text)
        self.assertIn('stage_requirements', text)
        self.assertNotIn('shutil.which', text, '探测逻辑必须只有一份（在 env_check 里）')
        for name in STAGE_OF:
            body = (SCRIPTS / name).read_text(encoding='utf-8')
            self.assertNotIn('shutil.which', body, '%s 不应自己探测依赖' % name)


class StageBoundaryTests(_Env):
    """stage→requires 必须与实际用法一致（否则门禁会拦错或放错）。"""

    def setUp(self):
        super().setUp()
        spec = importlib.util.spec_from_file_location('mito_env_check_boundary',
                                                      ROOT / 'tools' / 'env_check.py')
        self.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.module)
        self.manifest = self.module.load_manifest()

    def test_read_evidence_needs_only_samtools(self):
        self.assertEqual(sorted(self.manifest['stages']['read_evidence']['requires']), ['samtools'])

    def test_alignment_and_annotation_bridge_stages_exist(self):
        for stage, expected in (('read_alignment', ['minimap2']),
                                ('mitos2_bridge', ['biopython', 'python3']),
                                ('annot_check', ['biopython', 'python3'])):
            with self.subTest(stage=stage):
                self.assertEqual(sorted(self.manifest['stages'][stage]['requires']), expected)

    def test_every_gated_stage_exists_in_the_manifest(self):
        stages = set(self.manifest['stages'])
        for name, stage in STAGE_OF.items():
            with self.subTest(script=name):
                self.assertIn(stage, stages)


if __name__ == '__main__':
    unittest.main()

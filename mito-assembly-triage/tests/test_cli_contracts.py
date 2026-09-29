import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


class CliContracts(unittest.TestCase):
    def test_dependency_probe_finds_tools_in_separate_conda_envs(self):
        sys.path.insert(0, str(ROOT/'tools'))
        from env_check import probe_command
        with tempfile.TemporaryDirectory() as tmp:
            tool = Path(tmp) / 'envs' / 'assembly' / 'bin' / 'example-tool'
            tool.parent.mkdir(parents=True)
            tool.write_text('#!/bin/sh\necho example-tool 1.2.3\n')
            tool.chmod(0o755)
            with patch.dict(os.environ, {'CONDA_ROOT': tmp}), patch('env_check.shutil.which', return_value=None):
                found, version = probe_command('example-tool')
            self.assertEqual(found, str(tool))
            self.assertEqual(version, '1.2.3')

    def test_r_dependency_probe_selects_conda_rscript_with_required_package(self):
        sys.path.insert(0, str(ROOT/'tools'))
        from env_check import probe_r_package
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp)/'envs'/'a'/'bin'/'Rscript'
            good = Path(tmp)/'envs'/'z'/'bin'/'Rscript'
            for path, success in ((bad, 'false'), (good, 'true')):
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('#!/bin/sh\n'
                                'if [ "$1" = "--version" ]; then echo "R scripting front-end version 4.2.0"; exit 0; fi\n'
                                'if [ "$1" = "-e" ] && [ "' + success + '" = "true" ]; then echo 1.5.3; exit 0; fi\n'
                                'exit 1\n')
                path.chmod(0o755)
            with patch.dict(os.environ, {'CONDA_ROOT': tmp}), patch('env_check.shutil.which', return_value=None):
                found, version = probe_r_package('reshape2', search_conda=True)
            self.assertEqual(found, str(good))
            self.assertEqual(version, '1.5.3')

    def test_run_mitos2_uses_verified_rscript_without_shadowing_mitos_python(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pybin, helpers = root/'mitos-python'/'bin', root/'helpers'
            pybin.mkdir(parents=True)
            helpers.mkdir()
            interpreter = pybin/'python'
            interpreter.write_text('#!/bin/sh\n'
                                   'if [ "$1" = "--version" ]; then echo "Python 3.12.0"; exit 0; fi\n'
                                   'echo "$PATH"\n'
                                   'first="${PATH%%:*}"\n'
                                   'readlink "$first/Rscript"\n')
            interpreter.chmod(0o755)
            rscript = helpers/'Rscript'
            rscript.write_text('#!/bin/sh\nexit 0\n')
            rscript.chmod(0o755)
            outdir = root/'out'
            env = dict(os.environ, PATH='/usr/bin:/bin', CONDA_ROOT=str(root/'no-conda'),
                       MITOS2_PY=str(interpreter), MITOS2_REFDIR=str(root),
                       MITOS2_RSCRIPT=str(rscript), MITOS2_EXTRA_PATH=str(helpers),
                       MITO_ENV_GATE='off')
            result = subprocess.run(['bash', str(ROOT/'scripts/run_mitos2.sh'),
                                     '-i', str(root/'input.fa'), '-o', str(outdir), '-c', '5'],
                                    env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            output_lines = result.stdout.strip().splitlines()
            effective_path = output_lines[-2].split(':')
            runtime_dir = Path(effective_path[0])
            self.assertEqual(runtime_dir.parent, outdir)
            self.assertTrue(runtime_dir.name.startswith('.mitos2-runtime.'))
            self.assertEqual(Path(output_lines[-1]), rscript)
            self.assertEqual(effective_path[1], str(pybin))
            self.assertEqual(effective_path[2], str(helpers))
            self.assertFalse(runtime_dir.exists())  # temporary shim is cleaned up

    def test_help_does_not_require_external_tools(self):
        for name in ('depth_analysis', 'blast_genes', 'circularize', 'annot_check', 'mitos2_to_genbank',
                     'run_illumina_candidates', 'run_candidate_annotations', 'fastq_qc'):
            with self.subTest(tool=name):
                result = subprocess.run([sys.executable, str(ROOT/'scripts'/(name+'.py')), '--help'],
                    capture_output=True, text=True, env=dict(os.environ, PATH='/nonexistent', MITO_ENV_GATE='on'))
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_annotation_requires_explicit_table_before_loading_input(self):
        result = subprocess.run([sys.executable, str(ROOT/'scripts/annot_check.py'), 'missing.gb'],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn('--table', result.stdout)

    def test_invalid_window_and_thresholds(self):
        for extra in (['--window', '0'], ['--min-depth', '-1'], ['--min-mapq', '255']):
            result = subprocess.run([sys.executable, str(ROOT/'scripts/depth_analysis.py'),
                'sample.bam', 'seq.fa', *extra], capture_output=True, text=True)
            self.assertEqual(result.returncode, 1)
            self.assertNotIn('Traceback', result.stderr)


if __name__ == '__main__':
    unittest.main()

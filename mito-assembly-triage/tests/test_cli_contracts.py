import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class CliContracts(unittest.TestCase):
    def test_help_does_not_require_external_tools(self):
        for name in ('depth_analysis', 'blast_genes', 'circularize', 'annot_check', 'mitos2_to_genbank'):
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

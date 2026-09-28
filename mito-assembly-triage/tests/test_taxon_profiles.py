import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT/'scripts/taxon_profiles.py'


class TaxonProfileTests(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True)

    def test_show_cemus_profile_keeps_code_provisional_and_advisory(self):
        result = self.run_cli('show', 'Cemus macaoensis')
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('Cemus macaoensis', result.stdout)
        self.assertIn('Hemiptera', result.stdout)
        self.assertIn('Delphacidae', result.stdout)
        self.assertIn('table 5', result.stdout.lower())
        self.assertIn('provisional', result.stdout.lower())
        self.assertIn('自动选择: 否', result.stdout)

    def test_qjxh_species_profile_does_not_invent_species_taxid_or_code(self):
        result = self.run_cli('show', 'Neomohunia sinuatipenis')
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('Neomohunia sinuatipenis', result.stdout)
        self.assertIn('Mukariini', result.stdout)
        self.assertIn('未登记种级 NCBI TaxID', result.stdout)
        self.assertIn('species-specific', result.stdout)
        self.assertIn('自动选择: 否', result.stdout)

    def test_unknown_taxon_is_not_silently_guessed(self):
        result = self.run_cli('show', 'QJXH')
        self.assertEqual(result.returncode, 2)
        self.assertIn('未找到', result.stderr)

    def test_list_has_only_registered_taxa(self):
        result = self.run_cli('list')
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('Cemus macaoensis', result.stdout)
        self.assertIn('Neomohunia sinuatipenis', result.stdout)


if __name__ == '__main__':
    unittest.main()

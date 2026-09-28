import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANNOT_CHECK = ROOT / 'scripts' / 'annot_check.py'


@unittest.skipUnless(importlib.util.find_spec('Bio'), 'Biopython is required for GenBank fixture')
class PartialStopCompletionTests(unittest.TestCase):
    def run_check(self, transl_except):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqFeature import FeatureLocation, SeqFeature
        from Bio.SeqRecord import SeqRecord

        sequence = Seq('ATG' + 'AAA' * 9 + 'T')
        record = SeqRecord(sequence, id='partial_stop', name='partial_stop', description='fixture')
        record.annotations.update(molecule_type='DNA', topology='linear', data_file_division='INV')
        record.features = [
            SeqFeature(FeatureLocation(0, len(sequence)), type='source', qualifiers={
                'organism': ['Test insect'], 'mol_type': ['genomic DNA']}),
            SeqFeature(FeatureLocation(0, len(sequence)), type='CDS', qualifiers={
                'gene': ['cox1'], 'product': ['cytochrome c oxidase subunit I'],
                'codon_start': ['1'], 'transl_table': ['5'],
                'transl_except': [transl_except],
                'note': ["TAA stop codon is completed by the addition of 3' A residues to the mRNA"]}),
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            genbank = Path(temp_dir) / 'partial_stop.gb'
            SeqIO.write(record, genbank, 'genbank')
            return subprocess.run(
                [sys.executable, str(ANNOT_CHECK), str(genbank), '--table', '5'],
                text=True, capture_output=True, check=False)

    def test_terminal_term_exception_is_validated_at_exact_partial_stop(self):
        result = self.run_check('(pos:31,aa:TERM)')
        self.assertIn('PARTIAL_STOP_COMPLETION_VALIDATED: cox1', result.stdout)
        self.assertNotIn('TRANSL_EXCEPT_UNEXPLAINED: cox1', result.stdout)
        self.assertIn('未独立验证 RNA 加尾', result.stdout)

    def test_wrong_terminal_position_does_not_validate_partial_stop(self):
        result = self.run_check('(pos:28,aa:TERM)')
        self.assertNotIn('PARTIAL_STOP_COMPLETION_VALIDATED: cox1', result.stdout)
        self.assertIn('TRANSL_EXCEPT_UNEXPLAINED: cox1', result.stdout)


if __name__ == '__main__':
    unittest.main()

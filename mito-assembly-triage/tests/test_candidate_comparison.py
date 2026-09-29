import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('compare_assembly_candidates', ROOT / 'scripts' / 'compare_assembly_candidates.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CandidateComparisonTests(unittest.TestCase):
    def test_fasta_parser_preserves_records_and_normalizes_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'candidate.fa'
            path.write_text('>one description\nacgt N\n>two\nTTAA\n', encoding='ascii')
            self.assertEqual(MODULE.parse_fasta(path), [('one', 'ACGTN'), ('two', 'TTAA')])

    def test_fasta_parser_rejects_sequence_before_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'bad.fa'
            path.write_text('ACGT\n', encoding='ascii')
            with self.assertRaisesRegex(ValueError, '首个 FASTA 标识前'):
                MODULE.parse_fasta(path)

    def test_coverage_parser_reads_samtools_coverage_columns(self):
        result = MODULE.parse_coverage(
            '#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n'
            't0001_001\t1\t1000\t42\t900\t90.00\t12.50\t35.00\t58.00\n')
        self.assertEqual(result['t0001_001']['mapped_reads'], 42)
        self.assertEqual(result['t0001_001']['covered_bases'], 900)
        self.assertEqual(result['t0001_001']['breadth_percent'], 90.0)
        self.assertEqual(result['t0001_001']['mean_mapping_quality'], 58.0)


if __name__ == '__main__':
    unittest.main()

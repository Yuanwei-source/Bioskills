import importlib.util
import json
from pathlib import Path
import subprocess
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

    def test_pair_summary_rejects_duplicate_primary_qname_mate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            samtools = root/'samtools'
            rows = [
                'q1\t99\tt1\t1\t60\t4M\t=\t1\t0\tACGT\tIIII',
                'q1\t147\tt1\t1\t60\t4M\t=\t1\t0\tACGT\tIIII',
                'q2\t0\tt1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII',
            ]
            samtools.write_text('#!/bin/sh\nprintf "%s\\n" ' + ' '.join("'%s'" % row for row in rows) + '\n')
            samtools.chmod(0o755)
            counts, provenance = MODULE.summarize_primary_pairs(
                'input.bam', str(samtools), {'t1': 'candidate_001'}, 20, root/'view.log')
            self.assertEqual(counts['candidate_001']['primary_mapped_fragments'], 2)
            self.assertEqual(counts['candidate_001']['proper_pair_fragments_same_candidate'], 1)
            self.assertEqual(provenance['command'][2:4], ['-F', '3840'])
            rows.append(rows[0])
            samtools.write_text('#!/bin/sh\nprintf "%s\\n" ' + ' '.join("'%s'" % row for row in rows) + '\n')
            with self.assertRaisesRegex(ValueError, '重复 QNAME/mate'):
                MODULE.summarize_primary_pairs(
                    'input.bam', str(samtools), {'t1': 'candidate_001'}, 20, root/'duplicate.log')

    def test_coverage_command_excludes_all_nonprimary_and_duplicate_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidate = root/'candidate.fa'
            candidate.write_text('>mt\nACGTACGT\n')
            r1, r2 = root/'r1.fq', root/'r2.fq'
            r1.write_text('@r/1\nACGT\n+\nIIII\n')
            r2.write_text('@r/2\nACGT\n+\nIIII\n')
            manifest = root/'manifest.json'
            manifest.write_text(json.dumps({
                'format': 'mito-illumina-assembly-candidates-1',
                'candidates': [{'path': str(candidate), 'producer': 'synthetic'}]}))
            bwa = root/'bwa'
            bwa.write_text('#!/bin/sh\nexit 0\n')
            bwa.chmod(0o755)
            coverage_args = root/'coverage.args'
            samtools = root/'samtools'
            samtools.write_text(
                '#!/bin/sh\n'
                'if [ "$1" = "--version" ]; then echo "samtools 1.0"; exit 0; fi\n'
                'if [ "$1" = "sort" ]; then\n'
                '  out=; while [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
                '  cat > /dev/null; : > "$out"; exit 0\n'
                'fi\n'
                'if [ "$1" = "coverage" ]; then\n'
                '  printf "%s\\n" "$@" > "' + str(coverage_args) + '"\n'
                '  printf "#rname\\tstartpos\\tendpos\\tnumreads\\tcovbases\\tcoverage\\tmeandepth\\tmeanbaseq\\tmeanmapq\\n"\n'
                '  printf "t0001_001\\t1\\t8\\t0\\t0\\t0\\t0\\t0\\t0\\n"; exit 0\n'
                'fi\n'
                'exit 0\n')
            samtools.chmod(0o755)
            outdir = root/'evidence'
            result = subprocess.run([
                __import__('sys').executable, str(ROOT/'scripts/compare_assembly_candidates.py'),
                '--assembly-manifest', str(manifest), '--r1', str(r1), '--r2', str(r2),
                '--outdir', str(outdir), '--bwa', str(bwa), '--samtools', str(samtools)],
                capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertEqual(coverage_args.read_text().splitlines(), [
                'coverage', '-q', '20', '-Q', '20', '--ff',
                'UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY',
                str(outdir/'competitive_candidates.sorted.bam')])


if __name__ == '__main__':
    unittest.main()

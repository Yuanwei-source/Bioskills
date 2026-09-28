"""Real samtools/BLAST integration on generated public-free fixtures."""
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]


class ToolIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def bam(self, length=100, covered=100, count=10, chrom='mt', missing_quality=False, mapq=60):
        if not shutil.which('samtools'):
            self.skipTest('samtools not installed')
        sam = self.root/'input.sam'
        lines = ['@HD\tVN:1.6\tSO:coordinate', '@SQ\tSN:%s\tLN:%d' % (chrom, length)]
        for i in range(count):
            lines.append('\t'.join([str(i), str(16*(i%2)), chrom, '1', str(mapq), str(covered)+'M',
                                    '*', '0', '0', 'A'*covered, '*' if missing_quality else 'I'*covered]))
        sam.write_text('\n'.join(lines)+'\n')
        bam = self.root/'test.bam'
        subprocess.run(['samtools', 'sort', '-o', str(bam), str(sam)], check=True, capture_output=True)
        subprocess.run(['samtools', 'index', str(bam)], check=True, capture_output=True)
        return bam

    def run_depth(self, bam, length=100, chrom='mt', extra=()):
        fasta = self.root/'genome.fa'
        fasta.write_text('>'+chrom+'\n'+'A'*length+'\n')
        return subprocess.run([sys.executable, str(ROOT/'scripts/depth_analysis.py'), str(bam),
                               str(fasta), '--output-json', str(self.root/'depth.json'), *extra],
                              capture_output=True, text=True)

    def test_real_zero_coverage_returns_review(self):
        result = self.run_depth(self.bam(count=0))
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        report = json.loads((self.root/'depth.json').read_text())
        self.assertEqual(report['global_mean'], 0)
        self.assertEqual(report['zero_coverage_bases'], 100)

    def test_real_good_coverage(self):
        result = self.run_depth(self.bam())
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        report = json.loads((self.root/'depth.json').read_text())
        self.assertEqual(report['global_mean'], 10)

    def test_real_missing_quality_is_not_high_quality(self):
        result = self.run_depth(self.bam(missing_quality=True))
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        report = json.loads((self.root/'depth.json').read_text())
        self.assertEqual(report['global_mean'], 0)

    def test_real_unknown_mapq_excluded(self):
        result = self.run_depth(self.bam(mapq=255))
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        self.assertEqual(json.loads((self.root/'depth.json').read_text())['global_mean'], 0)

    def test_real_junction_does_not_count_gap_spans(self):
        if not shutil.which('samtools'):
            self.skipTest('samtools not installed')
        sys.path.insert(0, str(ROOT/'scripts'))
        from circularize import junction_evidence
        lines = ['@HD\tVN:1.6\tSO:coordinate', '@SQ\tSN:mt\tLN:300']
        for name, pos, cigar in [('deletion', 90, '10M100D10M'), ('skip', 90, '10M100N10M'),
                                 ('continuous', 95, '20M')]:
            lines.append('\t'.join([name, '0', 'mt', str(pos), '60', cigar,
                                     '*', '0', '0', 'A'*20, 'I'*20]))
        sam, bam = self.root/'junction.sam', self.root/'junction.bam'
        sam.write_text('\n'.join(lines)+'\n')
        subprocess.run(['samtools', 'sort', '-o', str(bam), str(sam)], check=True, capture_output=True)
        subprocess.run(['samtools', 'index', str(bam)], check=True, capture_output=True)
        evidence = junction_evidence(str(bam), 'mt:100-110')
        self.assertEqual(evidence['support'], 1)
        self.assertEqual(evidence['names'], ['continuous'])

    def test_real_reference_mismatch_is_error(self):
        result = self.run_depth(self.bam(chrom='other'))
        self.assertEqual(result.returncode, 1, result.stdout+result.stderr)
        self.assertFalse((self.root/'depth.json').exists())

    def test_real_last_window_weighting(self):
        result = self.run_depth(self.bam(length=501, covered=500), length=501)
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        report = json.loads((self.root/'depth.json').read_text())
        self.assertAlmostEqual(report['global_mean'], 5000/501)

    def test_real_missing_index_does_not_become_negative_evidence(self):
        bam = self.bam(count=0)
        Path(str(bam)+'.bai').unlink()
        result = self.run_depth(bam)
        self.assertEqual(result.returncode, 1, result.stdout+result.stderr)
        self.assertFalse((self.root/'depth.json').exists())

    def test_blast_two_copies_are_ambiguous_and_one_copy_works(self):
        if not all(shutil.which(tool) for tool in ('blastn', 'makeblastdb')):
            self.skipTest('BLAST+ not installed')
        rng = random.Random(713)
        sequence = ''.join(rng.choice('ACGT') for _ in range(300))
        record = SeqRecord(Seq(sequence), id='ref', annotations={'molecule_type': 'DNA'})
        record.features = [SeqFeature(FeatureLocation(0, len(sequence), strand=1),
                                      type='CDS', qualifiers={'gene': ['cox1']})]
        ref = self.root/'ref.gb'
        SeqIO.write(record, ref, 'genbank')
        reverse = str(Seq(sequence).reverse_complement())
        for name, sequences in [('forward', [sequence]), ('reverse', [reverse]),
                                ('two_forward', [sequence, sequence]), ('both', [sequence, reverse])]:
            with self.subTest(case=name):
                target = self.root/'target.fa'
                target.write_text('>mt\n'+('N'*100).join(sequences)+'\n')
                output = self.root/('order_%s.tsv' % name)
                result = subprocess.run([sys.executable, str(ROOT/'scripts/blast_genes.py'),
                    str(ref), str(target), '--out', str(output)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0 if len(sequences) == 1 else 1, result.stdout+result.stderr)
                self.assertEqual(output.exists(), len(sequences) == 1)
                if name == 'reverse':
                    self.assertEqual(output.read_text().split('\t')[4], '-')

    def test_blast_repeated_reference_names_use_coordinate_labels(self):
        if not all(shutil.which(tool) for tool in ('blastn', 'makeblastdb')):
            self.skipTest('BLAST+ not installed')
        rng = random.Random(9182)
        sequence = ''.join(rng.choice('ACGT') for _ in range(900))
        record = SeqRecord(Seq(sequence), id='ref', annotations={'molecule_type': 'DNA'})
        record.features = [
            SeqFeature(FeatureLocation(100, 250, strand=1), type='tRNA',
                       qualifiers={'product': ['tRNA-Leu']}),
            SeqFeature(FeatureLocation(500, 650, strand=-1), type='tRNA',
                       qualifiers={'product': ['tRNA-Leu']}),
        ]
        ref = self.root/'duplicate_reference.gb'
        SeqIO.write(record, ref, 'genbank')
        target = self.root/'target.fa'
        target.write_text('>mt\n'+sequence+'\n')
        output = self.root/'order.tsv'
        result = subprocess.run([sys.executable, str(ROOT/'scripts/blast_genes.py'),
            str(ref), str(target), '--out', str(output)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        rows = [line.split('\t') for line in output.read_text().splitlines()]
        self.assertTrue(all(len(row) == 8 for row in rows))
        self.assertTrue(all(row[7] == 'mt' for row in rows))
        self.assertEqual([row[0] for row in rows], ['tRNA-Leu@101..250:+', 'tRNA-Leu@501..650:-'])
        self.assertEqual([(row[2], row[3], row[4]) for row in rows],
                         [('101', '250', '+'), ('501', '650', '-')])


if __name__ == '__main__':
    unittest.main()

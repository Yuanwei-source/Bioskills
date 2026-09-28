"""Synthetic scientific counterexamples; no private sample data required."""
import contextlib
import importlib
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
with patch.dict(os.environ, {'MITO_ENV_GATE': 'off'}):
    depth = importlib.import_module('depth_analysis')
    circular = importlib.import_module('circularize')
    blast = importlib.import_module('blast_genes')


def sam(name='r', flag=0, pos=1, cigar='4M', seq='ACGA', qual='IIII', mapq=60, rg=None):
    row = [name, str(flag), 'mt', str(pos), str(mapq), cigar, '=', '1', '0', seq, qual]
    return row + (['RG:Z:' + rg] if rg else [])


def hsp(qlo, qhi, lo, hi, identity=100, strand='+'):
    return dict(qstart=qlo, qend=qhi, start=lo, end=hi, strand=strand,
                identity=identity, evalue=1e-10, query_len=100, aln_len=qhi-qlo+1,
                target='mt', type='CDS')


class BaseEvidenceTests(unittest.TestCase):
    def test_reverse_sam_already_oriented(self):
        self.assertEqual(depth.read_base_and_quality_at_reference(
            sam(flag=16, qual='I!5?'), 1), ('A', 40))
        self.assertEqual(depth.read_base_and_quality_at_reference(
            sam(flag=16, qual='I!5?'), 2), ('C', 0))

    def test_cigar_offsets(self):
        row = sam(cigar='2S2M1I2M2D2M', seq='TTACAGTCA', qual='I'*9)
        self.assertEqual(depth.read_base_at_reference(row, 3), 'G')
        self.assertIsNone(depth.read_base_at_reference(row, 5))
        self.assertEqual(depth.read_base_at_reference(row, 7), 'C')

    def test_missing_quality_never_callable(self):
        rows = [sam(str(i), flag=16*(i%2), seq='ATAT', qual='*') for i in range(8)]
        result = depth.base_support(rows, 2)
        self.assertFalse(result['callable'])
        self.assertEqual(result['depth'], 0)

    def test_overlapping_mates_count_once(self):
        rows = [sam(str(i), flag=flag, seq='ATAT') for i in range(3) for flag in (99, 147)]
        result = depth.base_support(rows, 1)
        self.assertEqual(result['depth'], 3)
        self.assertFalse(result['callable'])

    def test_conflicting_mates_are_not_a_consensus(self):
        result = depth.base_support([sam('pair', 99), sam('pair', 147, seq='TCGA')], 1)
        self.assertEqual(result['depth'], 0)

    def test_good_distinct_templates_remain_usable(self):
        rows = [sam(str(i), flag=16*(i%2), seq='ATAT') for i in range(6)]
        result = depth.base_support(rows, 1)
        self.assertTrue(result['callable'])
        self.assertEqual(result['depth'], 6)

    def test_unknown_mapq_and_ambiguous_bases_excluded(self):
        rows = [sam(str(i), mapq=255, seq='ATAT') for i in range(6)]
        self.assertEqual(depth.base_support(rows, 1)['depth'], 0)
        rows = [sam(str(i), seq='NNNN') for i in range(6)]
        self.assertEqual(depth.base_support(rows, 1)['depth'], 0)

    def test_same_name_in_different_read_groups(self):
        result = depth.base_support([sam(rg='a'), sam(rg='b')], 1)
        self.assertEqual(result['depth'], 2)


class JunctionTests(unittest.TestCase):
    def test_overlap_must_be_spanned_with_external_flanks(self):
        self.assertEqual(circular.join_anchor_region(100, 30, 180, 10),
                         'mitogenome_candidate:61-110')
        with self.assertRaises(ValueError):
            circular.join_anchor_region(100, 100, 180, 10)

    def evidence(self, rows, region='mt:100-110'):
        with patch.object(circular.subprocess, 'run', return_value=
                          subprocess.CompletedProcess([], 0, '\n'.join('\t'.join(r) for r in rows), '')):
            return circular.junction_evidence('test.bam', region)

    def test_deletion_and_skip_cannot_support_junction(self):
        for op in ('D', 'N'):
            with self.subTest(op=op):
                self.assertEqual(self.evidence([sam(pos=90, cigar='10M100'+op+'10M',
                                                   seq='A'*20, qual='I'*20)])['support'], 0)

    def test_exclusive_end_is_not_a_covered_base(self):
        self.assertEqual(self.evidence([sam(pos=100, cigar='10M', seq='A'*10, qual='I'*10)])['support'], 0)

    def test_contiguous_anchor_supported(self):
        self.assertEqual(self.evidence([sam(pos=95, cigar='20M', seq='A'*20, qual='I'*20)])['support'], 1)

    def test_insertion_at_junction_is_not_exact_adjacency(self):
        self.assertEqual(self.evidence([sam(pos=95, cigar='10M2I10M', seq='A'*22,
                                           qual='I'*22)])['support'], 0)

    def test_read_and_template_counts_are_distinct(self):
        rows = [sam('pair', f, 95, '20M', 'A'*20, 'I'*20) for f in (99, 147)]
        result = self.evidence(rows)
        self.assertEqual(result['support'], 1)
        self.assertEqual(sum(result['strand_counts'].values()), 1)


class BlastTests(unittest.TestCase):
    def test_two_loci_remain_two_candidates(self):
        result = blast.locus_candidates([hsp(1, 100, 1, 100), hsp(1, 100, 1001, 1100)])
        self.assertEqual(len(result), 2)
        selected, errors = blast.select_gene_hits({'cox1': result})
        self.assertFalse(selected)
        self.assertTrue(errors)

    def test_distinct_partial_chains_keep_separate_loci(self):
        result = blast.locus_candidates([hsp(1, 50, 1, 50), hsp(51, 100, 51, 100),
                                        hsp(1, 50, 1001, 1050), hsp(51, 100, 1051, 1100)])
        self.assertEqual(len(result), 2)
        self.assertTrue(all(hit['aln_len'] == 100 for hit in result))

    def test_two_complete_loci_cannot_be_merged(self):
        self.assertIsNone(blast.merge_hsps([hsp(1, 100, 1, 100), hsp(1, 100, 1001, 1100)]))

    def test_weighted_identity(self):
        result = blast.merge_hsps([hsp(1, 10, 1, 10, 100), hsp(11, 100, 11, 100, 60)])
        self.assertAlmostEqual(result['identity'], 64)
        selected, errors = blast.select_gene_hits({'cox1': [result]})
        self.assertFalse(selected)
        self.assertTrue(errors)

    def test_target_bases_cannot_be_reused(self):
        self.assertIsNone(blast.merge_hsps([hsp(1, 50, 1, 50), hsp(51, 100, 1, 50)]))

    def test_reverse_collinear_fragments_merge(self):
        result = blast.merge_hsps([hsp(1, 50, 151, 200, strand='-'),
                                  hsp(51, 100, 101, 150, strand='-')])
        self.assertIsNotNone(result)
        self.assertEqual(result['aln_len'], 100)


if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python3
"""Failing reproductions for the expert review findings in the Pi lane:
`annot_check.py` biological semantics and `cox1_id.py` structured BLAST output.

Regressions covered (expert T-matrix items in this lane):
  T11  GenBank partial CDS / codon_start contradiction / transl_except mismatch
  T12  circular rotation and whole-genome reverse complement must not read as
       a gene rearrangement; gene-set loss must not be conflated with order
  T13  /product naming, bare trnL/trnS, unknown rRNA name, ambiguous overlap exception
  T14  structured BLAST parsing (multi-HSP coverage, wrapped descriptions,
       no-hit, duplicate accessions, network vs no-match)

`depth_analysis.py`, `blast_genes.py` and `circularize.py` are deliberately out of
scope here (other lane).
"""
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gb_fixtures import (  # noqa: E402
    ORDER, biopython_available, build_synthetic_ref, cds_block, load_module,
    run_annot_check, write_gb, write_gb_raw,
)

REASON = "合成测试: 单基因记录, 基因集差异已逐项确认"


def _cds_sequence(length=102, seed=1, terminal=None):
    """A clean table-5 CDS; `terminal` overrides the last 3 bases (e.g. '' or 'TA')."""
    block = cds_block(length, seed=seed)
    if terminal is not None:
        block = block[:-len(terminal)] + terminal if terminal else block[:-3]
    return block


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class CdSStartAndPartialTests(unittest.TestCase):
    """T11a/b: partiality comes from the location, not from /codon_start alone."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def test_location_partial_helper_reads_both_strands(self):
        module = load_module("annot_partial", Path("scripts") / "annot_check.py")
        from Bio import SeqIO
        path = self.dir / "partial.gb"
        write_gb_raw(path, "A" * 400, [
            {"location": "<1..100", "type": "CDS", "gene": "p5"},
            {"location": "200..>300", "type": "CDS", "gene": "p3"},
            {"location": "complement(<250..>350)", "type": "CDS", "gene": "p5rc"},
            {"location": "1..100", "type": "CDS", "gene": "full"},
        ])
        record = SeqIO.read(path, "genbank")
        by_gene = {f.qualifiers["gene"][0]: f for f in record.features if f.type == "CDS"}
        self.assertEqual(module.feature_partial(by_gene["p5"]), (True, False))
        self.assertEqual(module.feature_partial(by_gene["p3"]), (False, True))
        self.assertEqual(module.feature_partial(by_gene["p5rc"]), (True, True))
        self.assertEqual(module.feature_partial(by_gene["full"]), (False, False))

    def test_five_prime_partial_does_not_require_a_start_codon(self):
        # no valid ATN start at all, but the location says the 5' end is missing
        sequence = "TTT" + "AAA" * 20 + "TAA"
        path = self.dir / "p5.gb"
        write_gb_raw(path, sequence + "A" * 50, [
            {"location": "<1..%d" % len(sequence), "type": "CDS", "gene": "cox1"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertNotIn("[ERROR]", result.stdout)
        self.assertIn("PARTIAL_CDS_5P", result.stdout)

    def test_complete_location_with_codon_start_two_is_a_contradiction(self):
        sequence = "GG" + "ATG" + "AAA" * 20 + "TAA"
        path = self.dir / "contradiction.gb"
        write_gb_raw(path, sequence, [
            {"location": "1..%d" % len(sequence), "type": "CDS", "gene": "cox1",
             "codon_start": 2},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("CODON_START_CONFLICT", result.stdout)
        self.assertIn("[WARN ]", result.stdout)

    def test_three_prime_partial_missing_stop_is_not_an_error(self):
        sequence = "ATG" + "AAA" * 20  # no stop codon at all
        path = self.dir / "p3.gb"
        write_gb_raw(path, sequence + "A" * 50, [
            {"location": "1..>%d" % len(sequence), "type": "CDS", "gene": "cox1"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertNotIn("[ERROR]", result.stdout)
        self.assertIn("PARTIAL_CDS_3P", result.stdout)

    def test_internal_stop_is_a_conflict_not_a_partial(self):
        sequence = "ATG" + "AAA" * 5 + "TAA" + "AAA" * 5 + "TAA"
        path = self.dir / "internal.gb"
        write_gb_raw(path, sequence + "A" * 50, [
            {"location": "1..%d" % len(sequence), "type": "CDS", "gene": "cox1"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("[ERROR]", result.stdout)
        self.assertIn("内部终止", result.stdout)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class TranslExceptTests(unittest.TestCase):
    """T11c: /transl_except must be matched against the stop it claims to explain."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def _gb(self, name, transl_except):
        # ATG AAA AAA TAA AAA AAA TAA -> one internal stop at codon 4 (record 10..12)
        sequence = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA"
        spec = {"location": "1..%d" % len(sequence), "type": "CDS", "gene": "cox1"}
        if transl_except:
            spec["transl_except"] = transl_except
        path = self.dir / name
        write_gb_raw(path, sequence + "A" * 40, [spec])
        return path

    def test_matching_exception_explains_the_internal_stop(self):
        result = run_annot_check(self._gb("match.gb", "(pos:10..12,aa:Trp)"),
                                 "--allow-atypical", REASON)
        self.assertIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertNotIn("[ERROR]", result.stdout)

    def test_non_matching_exception_does_not_excuse_the_stop(self):
        # codon 4 is the only internal stop; declaring it at codon 2 explains nothing
        result = run_annot_check(self._gb("nomatch.gb", "(pos:4..6,aa:Trp)"),
                                 "--allow-atypical", REASON)
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)
        self.assertIn("内部终止", result.stdout)

    def test_exception_without_any_internal_stop_is_flagged(self):
        sequence = "ATG" + "AAA" * 20 + "TAA"   # clean CDS, nothing to explain
        path = self.dir / "clean.gb"
        write_gb_raw(path, sequence + "A" * 40, [
            {"location": "1..%d" % len(sequence), "type": "CDS", "gene": "cox1",
             "transl_except": "(pos:10..12,aa:Trp)"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class GeneOrderCycleTests(unittest.TestCase):
    """T12: rotation and whole-genome reverse complement are representation changes."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)
        self.module = load_module("annot_cycle", Path("scripts") / "annot_check.py")

    def _cycle(self, genes):
        return [(name, strand) for name, _kind, _length, strand in genes]

    def test_cycle_helpers_treat_rotation_and_reverse_as_equivalent(self):
        base = self._cycle(ORDER)
        rotated = base[3:] + base[:3]
        reversed_cycle = [(name, "-" if strand == "+" else "+")
                          for name, strand in reversed(base)]
        self.assertTrue(self.module.cycles_equivalent(base, rotated))
        self.assertTrue(self.module.cycles_equivalent(base, reversed_cycle))
        swapped = base[:]
        swapped[4], swapped[5] = swapped[5], swapped[4]
        self.assertFalse(self.module.cycles_equivalent(base, swapped))

    def test_adjacency_diff_reports_only_real_arrangement_changes(self):
        base = self._cycle(ORDER)
        rotated = base[7:] + base[:7]
        self.assertEqual(self.module.adjacency_diff(base, rotated), [])
        swapped = base[:]
        swapped[4], swapped[5] = swapped[5], swapped[4]
        self.assertTrue(self.module.adjacency_diff(base, swapped))

    def test_rotated_query_origin_is_not_a_rearrangement(self):
        ref, sequence = build_synthetic_ref(self.dir)
        shift = ORDER[0][2]        # rotate by exactly one gene -> no feature wraps
        total = len(sequence)
        features = []
        for name, kind, length, strand in ORDER:
            start = sum(item[2] for item in ORDER[:len(features)])
            features.append({"gene": name, "type": kind, "strand": strand,
                             "start": (start - shift) % total, "end": (start - shift) % total + length})
        rotated = self.dir / "rotated.gb"
        write_gb(rotated, sequence, features)
        result = run_annot_check(rotated, "--ref", ref)
        self.assertIn("邻接关系一致", result.stdout)
        self.assertNotIn("ARRANGEMENT_DIFF", result.stdout)

    def test_missing_gene_is_reported_separately_from_order(self):
        ref, sequence = build_synthetic_ref(self.dir)
        from Bio import SeqIO
        record = SeqIO.read(ref, "genbank")
        record.features = [f for f in record.features
                           if f.qualifiers.get("gene", [""])[0] != "nad3"]
        missing = self.dir / "missing.gb"
        SeqIO.write(record, missing, "genbank")
        result = run_annot_check(missing, "--ref", ref)
        self.assertIn("GENE_SET_DIFF", result.stdout)
        self.assertNotIn("ARRANGEMENT_DIFF", result.stdout)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class ProductNamingTests(unittest.TestCase):
    """T13: natural-language /product names and undetermined identities."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)
        self.module = load_module("annot_names", Path("scripts") / "annot_check.py")

    def test_product_names_map_to_standard_identities(self):
        cases = {
            "NADH dehydrogenase subunit 5": "nd5",
            "NADH dehydrogenase subunit 4L": "nd4l",
            "cytochrome c oxidase subunit 1": "cox1",
            "cytochrome c oxidase subunit III": "cox3",
            "cytochrome b": "cytb",
            "ATP synthase F0 subunit 6": "atp6",
            "16S ribosomal RNA": "rrnl",
            "12S ribosomal RNA": "rrns",
        }
        for product, expected in cases.items():
            with self.subTest(product=product):
                self.assertEqual(self.module._canonical_key(product), expected)

    def test_trna_product_does_not_guess_leu_or_ser_subtype(self):
        self.assertEqual(self.module._canonical_key("tRNA-Leu"), "trnl")
        self.assertEqual(self.module._canonical_key("tRNA-Ser"), "trns")
        self.assertEqual(self.module._canonical_key("tRNA-Ala"), "trna")

    def test_unknown_rrna_name_is_undetermined_not_rrns(self):
        self.assertEqual(self.module._canonical_key("ribosomal RNA"), "rrna")

    def test_unknown_rrna_does_not_borrow_the_rrns_length_band(self):
        path = self.dir / "unknown-rrna.gb"
        write_gb_raw(path, "A" * 1500, [
            {"location": "1..1300", "type": "rRNA", "gene": "rrna"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertNotIn("超出范围 600-850", result.stdout)
        self.assertIn("UNDETERMINED_RRNA", result.stdout)

    def test_bare_trnl_from_product_is_undetermined(self):
        path = self.dir / "product-trna.gb"
        write_gb_raw(path, "A" * 300, [
            {"location": "1..70", "type": "tRNA", "product": "tRNA-Leu"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("UNDETERMINED_TRNA", result.stdout)

    def test_ambiguous_overlap_exception_warns_it_cannot_be_located(self):
        # two genes share the identity "cox1" -> a name-only exception is ambiguous
        block = cds_block(102, seed=41)
        path = self.dir / "dup.gb"
        write_gb_raw(path, block * 2 + "A" * 200, [
            {"location": "1..102", "type": "CDS", "gene": "cox1"},
            {"location": "62..164", "type": "CDS", "gene": "cox1"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON,
                                 "--tolerate-overlap", "cox1,cox1")
        self.assertIn("TOLERATE_OVERLAP_AMBIGUOUS", result.stdout)


class ContractWordingTests(unittest.TestCase):
    def test_require_circular_states_it_only_checks_the_declaration(self):
        result = run_annot_check("/nonexistent.gb", "--require-circular")
        self.assertIn("CIRCULAR_DECLARATION_CHECK", result.stdout + result.stderr)

    def test_unregistered_start_exception_carries_a_limitation(self):
        if not biopython_available():
            self.skipTest("Biopython is not installed")
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        directory = Path(temp.name)
        sequence = "CGA" + "AAA" * 20 + "TAA"
        path = directory / "cga.gb"
        write_gb_raw(path, sequence + "A" * 40, [
            {"location": "1..%d" % len(sequence), "type": "CDS", "gene": "cox1"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON,
                                 "--tolerate-start", "cox1:CGA")
        self.assertIn("EXCEPTION_NOT_REGISTERED", result.stdout)


class BlastXmlParsingTests(unittest.TestCase):
    """T14: structured BLAST output instead of text regexes."""

    def setUp(self):
        self.module = load_module("cox1_structured", Path("scripts") / "cox1_id.py")

    @staticmethod
    def _hsp(q_from, q_to, identity, align_len, bitscore=500):
        return ("<Hsp><Hsp_bit-score>%d</Hsp_bit-score>"
                "<Hsp_query-from>%d</Hsp_query-from><Hsp_query-to>%d</Hsp_query-to>"
                "<Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>%d</Hsp_hit-to>"
                "<Hsp_identity>%d</Hsp_identity><Hsp_align-len>%d</Hsp_align-len></Hsp>"
                % (bitscore, q_from, q_to, align_len, identity, align_len))

    def _xml(self, hits, query_len=1000):
        body = "".join(hits)
        return ("<?xml version=\"1.0\"?><BlastOutput><BlastOutput_db>nt</BlastOutput_db>"
                "<BlastOutput_query-len>%d</BlastOutput_query-len><BlastOutput_iterations>"
                "<Iteration><Iteration_query-len>%d</Iteration_query-len>"
                "<Iteration_hits>%s</Iteration_hits></Iteration>"
                "</BlastOutput_iterations></BlastOutput>" % (query_len, query_len, body))

    def test_multi_hsp_hit_coverage_is_the_query_union(self):
        hit = ("<Hit><Hit_accession>NC_1</Hit_accession>"
               "<Hit_def>Species A mitochondrion</Hit_def><Hit_len>15000</Hit_len><Hit_hsps>"
               + self._hsp(1, 400, 400, 400) + self._hsp(601, 1000, 400, 400)
               + "</Hit_hsps></Hit>")
        parsed = self.module.parse_blast_xml(self._xml([hit]))
        self.assertEqual(parsed["query_len"], 1000)
        self.assertEqual(len(parsed["hits"]), 1)
        best = parsed["hits"][0]
        self.assertAlmostEqual(best["coverage"], 0.8, places=6)
        self.assertAlmostEqual(best["identity"], 100.0, places=6)
        self.assertEqual(best["hsp_count"], 2)

    def test_overlapping_hsps_are_not_double_counted(self):
        hit = ("<Hit><Hit_accession>NC_2</Hit_accession><Hit_def>x</Hit_def><Hit_len>15000</Hit_len>"
               "<Hit_hsps>" + self._hsp(1, 500, 480, 500) + self._hsp(400, 800, 380, 400)
               + "</Hit_hsps></Hit>")
        best = self.module.parse_blast_xml(self._xml([hit]))["hits"][0]
        self.assertAlmostEqual(best["coverage"], 0.8, places=6)

    def test_wrapped_description_is_normalised(self):
        hit = ("<Hit><Hit_accession>NC_3</Hit_accession>"
               "<Hit_def>Species B\n mitochondrion,\n complete genome</Hit_def>"
               "<Hit_len>15000</Hit_len><Hit_hsps>" + self._hsp(1, 1000, 990, 1000)
               + "</Hit_hsps></Hit>")
        best = self.module.parse_blast_xml(self._xml([hit]))["hits"][0]
        self.assertNotIn("\n", best["description"])
        self.assertIn("complete genome", best["description"])

    def test_no_hits_is_a_status_not_a_failure(self):
        parsed = self.module.parse_blast_xml(self._xml([]))
        self.assertEqual(parsed["hits"], [])

    def test_duplicate_accessions_are_one_candidate(self):
        def hit(accession, identity):
            return ("<Hit><Hit_accession>%s</Hit_accession><Hit_def>same species</Hit_def>"
                    "<Hit_len>15000</Hit_len><Hit_hsps>%s</Hit_hsps></Hit>"
                    % (accession, self._hsp(1, 1000, identity, 1000)))
        parsed = self.module.parse_blast_xml(self._xml([hit("NC_9", 995), hit("NC_9", 990)]))
        distinct = self.module.distinct_candidates(parsed["hits"])
        self.assertEqual(len(distinct), 1)

    def test_search_info_status_is_parsed(self):
        waiting = "<BlastSearchInfo><Status>WAITING</Status><RID>ABC123</RID></BlastSearchInfo>"
        ready = "<BlastSearchInfo><Status>READY</Status><RID>ABC123</RID></BlastSearchInfo>"
        self.assertEqual(self.module.parse_search_info(waiting), ("WAITING", "ABC123"))
        self.assertEqual(self.module.parse_search_info(ready), ("READY", "ABC123"))
        self.assertEqual(self.module.parse_search_info("not xml"), (None, None))


if __name__ == "__main__":
    unittest.main()

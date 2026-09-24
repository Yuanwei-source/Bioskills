#!/usr/bin/env python3
"""Pre-merge boundary checks requested by the reviewer (Pi lane).

1. feature_partial() must read Biopython's BeforePosition/AfterPosition objects,
   not the location string; must handle minus strand, compound/cross-origin
   locations and both-ends partial.
2. /transl_except must match the actual codon and reading frame, including on the
   minus strand, and must never explain a second unexplained internal stop.
3. COX1 XML2 multi-HSP aggregation must not double count overlapping query bases
   into a precise-looking composite identity; conflicting HSPs are AMBIGUOUS and
   cross-accession candidates are retained.
"""
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gb_fixtures import (  # noqa: E402
    ROOT, biopython_available, load_module, run_annot_check, write_gb_raw,
)

REASON = "合成测试: 单基因记录, 基因集差异已逐项确认"


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class PartialPositionObjectTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module("annot_partial_obj", Path("scripts") / "annot_check.py")
        from Bio.SeqFeature import (AfterPosition, BeforePosition, CompoundLocation,
                                    ExactPosition, FeatureLocation, SeqFeature)
        self.Before = BeforePosition
        self.After = AfterPosition
        self.Exact = ExactPosition
        self.FeatureLocation = FeatureLocation
        self.CompoundLocation = CompoundLocation
        self.SeqFeature = SeqFeature

    def test_object_result_is_authoritative_and_conflict_is_reported(self):
        class FakeLocation:
            def __init__(self, parts, strand):
                self.parts = parts
                self.strand = strand

            def __str__(self):
                return "[0:100](+)"          # string says "complete"

        feature = SimpleNamespace(location=FakeLocation(
            [SimpleNamespace(start=self.Before(0), end=self.Exact(100))], 1))
        detail = self.module.feature_partial_detail(feature)
        self.assertEqual(detail["source"], "position_objects")
        self.assertEqual((detail["five"], detail["three"]), (True, False))
        self.assertTrue(detail["conflict"])
        self.assertEqual(detail["string_result"], (False, False))
        # the object result wins: never silently downgraded to "complete"
        self.assertEqual(self.module.feature_partial(feature), (True, False))

    def test_no_conflict_when_objects_and_string_agree(self):
        location = self.CompoundLocation([
            self.FeatureLocation(self.Exact(0), self.After(10), strand=-1),
            self.FeatureLocation(self.Before(749), self.Exact(780), strand=-1),
        ])
        detail = self.module.feature_partial_detail(self._feature(location))
        self.assertEqual(detail["source"], "position_objects")
        self.assertFalse(detail["conflict"])
        self.assertEqual((detail["five"], detail["three"]), (True, True))

    def _feature(self, location):
        return self.SeqFeature(location, type="CDS", qualifiers={"gene": ["x"]})

    def test_plus_strand_partials_from_position_objects(self):
        five = self._feature(self.FeatureLocation(self.Before(0), self.Exact(100), strand=1))
        three = self._feature(self.FeatureLocation(self.Exact(100), self.After(200), strand=1))
        both = self._feature(self.FeatureLocation(self.Before(100), self.After(200), strand=1))
        plain = self._feature(self.FeatureLocation(0, 100, strand=1))
        self.assertEqual(self.module.feature_partial(five), (True, False))
        self.assertEqual(self.module.feature_partial(three), (False, True))
        self.assertEqual(self.module.feature_partial(both), (True, True))
        self.assertEqual(self.module.feature_partial(plain), (False, False))

    def test_minus_strand_partials_follow_the_transcription_direction(self):
        # 5' end of a minus-strand gene is the high coordinate -> AfterPosition on end
        five = self._feature(self.FeatureLocation(self.Exact(200), self.After(300), strand=-1))
        three = self._feature(self.FeatureLocation(self.Before(200), self.Exact(300), strand=-1))
        both = self._feature(self.FeatureLocation(self.Before(200), self.After(300), strand=-1))
        self.assertEqual(self.module.feature_partial(five), (True, False))
        self.assertEqual(self.module.feature_partial(three), (False, True))
        self.assertEqual(self.module.feature_partial(both), (True, True))

    def test_cross_origin_compound_location_both_ends_partial(self):
        # parts are in transcription order: first part carries the 5' end
        location = self.CompoundLocation([
            self.FeatureLocation(self.Exact(0), self.After(10), strand=-1),
            self.FeatureLocation(self.Before(749), self.Exact(780), strand=-1),
        ])
        self.assertEqual(self.module.feature_partial(self._feature(location)), (True, True))
        plain = self.CompoundLocation([
            self.FeatureLocation(self.Exact(349), self.Exact(400), strand=1),
            self.FeatureLocation(self.Exact(0), self.Exact(50), strand=1),
        ])
        self.assertEqual(self.module.feature_partial(self._feature(plain)), (False, False))

    def test_position_objects_are_used_even_when_the_string_shows_no_markers(self):
        class FakeLocation:
            def __init__(self, parts, strand):
                self.parts = parts
                self.strand = strand

            def __str__(self):
                return "join{[0:100](+)}"      # deliberately marker-free

        feature = SimpleNamespace(
            location=FakeLocation([SimpleNamespace(start=self.Before(0), end=self.Exact(100))], 1))
        self.assertEqual(self.module.feature_partial(feature), (True, False))

    def test_location_string_is_only_a_fallback(self):
        class FakeLocation:
            def __init__(self, parts, strand):
                self.parts = parts
                self.strand = strand

            def __str__(self):
                return "[<0:100](+)"

        # plain ints: no position objects at all, so the string fallback must fire
        feature = SimpleNamespace(
            location=FakeLocation([SimpleNamespace(start=0, end=100)], 1))
        self.assertEqual(self.module.feature_partial(feature), (True, False))


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class TranslExceptReadingFrameTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    @staticmethod
    def _revcomp(sequence):
        from Bio.Seq import Seq
        return str(Seq(sequence).reverse_complement())

    def test_minus_strand_exception_is_matched_at_the_correct_codon(self):
        coding = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA"   # internal stop = codon 4
        stored = self._revcomp(coding)
        path = self.dir / "minus-match.gb"
        write_gb_raw(path, stored + "A" * 40, [
            {"location": "complement(1..%d)" % len(coding), "type": "CDS", "gene": "cox1",
             "transl_except": "(pos:10..12,aa:Trp)"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertNotIn("[ERROR]", result.stdout)

    def test_minus_strand_exception_at_the_wrong_codon_stays_unexplained(self):
        coding = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA"
        stored = self._revcomp(coding)
        path = self.dir / "minus-nomatch.gb"
        write_gb_raw(path, stored + "A" * 40, [
            {"location": "complement(1..%d)" % len(coding), "type": "CDS", "gene": "cox1",
             "transl_except": "(pos:4..6,aa:Trp)"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_one_exception_cannot_excuse_two_internal_stops(self):
        # two internal stops: codon 4 (bases 10..12) and codon 7 (bases 19..21)
        coding = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA"
        path = self.dir / "two-stops.gb"
        write_gb_raw(path, coding + "A" * 40, [
            {"location": "1..%d" % len(coding), "type": "CDS", "gene": "cox1",
             "transl_except": "(pos:10..12,aa:Trp)"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)
        self.assertIn("1 个内部终止", result.stdout)

    def test_non_codon_position_range_does_not_explain_anything(self):
        coding = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA"
        path = self.dir / "bad-range.gb"
        write_gb_raw(path, coding + "A" * 40, [
            {"location": "1..%d" % len(coding), "type": "CDS", "gene": "cox1",
             "transl_except": "(pos:8..13,aa:Trp)"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class ExceptionRegistryTests(unittest.TestCase):
    """Unregistered or non-matching exceptions must stay REVIEW, never validated."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def _cga(self):
        sequence = "CGA" + "AAA" * 20 + "TAA"
        path = self.dir / "cga.gb"
        write_gb_raw(path, sequence + "A" * 40, [
            {"location": "1..%d" % len(sequence), "type": "CDS", "gene": "cox1"},
        ])
        return path

    def test_registry_entry_with_missing_fields_stays_review(self):
        import json
        registry = self.dir / "registry.json"
        registry.write_text(json.dumps({"exceptions": [
            {"gene": "cox1", "codon": "CGA"}]}), encoding="utf-8")
        result = run_annot_check(self._cga(), "--allow-atypical", REASON,
                                 "--tolerate-start", "cox1:CGA",
                                 "--exception-registry", str(registry))
        self.assertIn("EXCEPTION_RECORD_INCOMPLETE", result.stdout)
        self.assertNotIn("[ERROR]", result.stdout)

    def test_taxon_mismatch_is_reported(self):
        import json
        registry = self.dir / "registry2.json"
        registry.write_text(json.dumps({"exceptions": [
            {"gene": "cox1", "codon": "CGA", "taxon": "Lepidoptera",
             "source": "DOI 10.1000/x", "rationale": "documented CGA start"}]}),
            encoding="utf-8")
        result = run_annot_check(self._cga(), "--allow-atypical", REASON,
                                 "--tolerate-start", "cox1:CGA",
                                 "--exception-registry", str(registry),
                                 "--taxon", "Diptera")
        self.assertIn("EXCEPTION_TAXON_MISMATCH", result.stdout)

    def test_matching_registry_entry_is_quoted(self):
        import json
        registry = self.dir / "registry3.json"
        registry.write_text(json.dumps({"exceptions": [
            {"gene": "cox1", "codon": "CGA", "taxon": "Lepidoptera",
             "source": "DOI 10.1000/x", "rationale": "documented CGA start"}]}),
            encoding="utf-8")
        result = run_annot_check(self._cga(), "--allow-atypical", REASON,
                                 "--tolerate-start", "cox1:CGA",
                                 "--exception-registry", str(registry),
                                 "--taxon", "Lepidoptera")
        self.assertIn("taxon=Lepidoptera", result.stdout)
        self.assertIn("NONCANONICAL_START_REVIEW", result.stdout)
        self.assertNotIn("EXCEPTION_TAXON_MISMATCH", result.stdout)
        self.assertNotIn("EXCEPTION_RECORD_INCOMPLETE", result.stdout)


class HspCollinearityTests(unittest.TestCase):
    """Non-overlapping HSPs must still be collinear on the target to be summari­sable."""

    def setUp(self):
        self.module = load_module("cox1_collinear", Path("scripts") / "cox1_id.py")

    @staticmethod
    def _hsp(query_from, query_to, hit_from, hit_to, identity, align_len, bitscore=500):
        return ("<Hsp><Hsp_bit-score>%d</Hsp_bit-score>"
                "<Hsp_query-from>%d</Hsp_query-from><Hsp_query-to>%d</Hsp_query-to>"
                "<Hsp_hit-from>%d</Hsp_hit-from><Hsp_hit-to>%d</Hsp_hit-to>"
                "<Hsp_identity>%d</Hsp_identity><Hsp_align-len>%d</Hsp_align-len></Hsp>"
                % (bitscore, query_from, query_to, hit_from, hit_to, identity, align_len))

    def _hit(self, accession, hsps, query_len=1000):
        xml = ("<?xml version=\"1.0\"?><BlastOutput><BlastOutput_db>nt</BlastOutput_db>"
               "<BlastOutput_iterations><Iteration><Iteration_query-len>%d</Iteration_query-len>"
               "<Iteration_hits><Hit><Hit_accession>%s</Hit_accession><Hit_def>species</Hit_def>"
               "<Hit_len>15000</Hit_len><Hit_hsps>%s</Hit_hsps></Hit>"
               "</Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>"
               % (query_len, accession, hsps))
        return self.module.parse_blast_xml(xml)["hits"][0]

    def test_collinear_non_overlapping_hsps_are_aggregated(self):
        hit = self._hit("NC_A", self._hsp(1, 400, 1, 400, 400, 400)
                        + self._hsp(601, 1000, 500, 899, 380, 400))
        self.assertTrue(hit["collinearity"]["collinear"])
        self.assertFalse(hit["ambiguous_alignment"])
        self.assertIsNotNone(hit["identity"])
        self.assertAlmostEqual(hit["coverage"], 0.8, places=6)

    def test_scattered_hsps_are_not_treated_as_one_alignment(self):
        # query coverage is high, but the two halves hit far-apart target positions
        hit = self._hit("NC_B", self._hsp(1, 400, 1, 400, 400, 400)
                        + self._hsp(601, 1000, 9000, 9399, 390, 400))
        self.assertFalse(hit["collinearity"]["collinear"])
        self.assertTrue(hit["ambiguous_alignment"])
        self.assertTrue(hit["conflicting_alignment"])
        self.assertIn("non_collinear_hsps", hit["ambiguity_reasons"])
        self.assertIsNone(hit["identity"])
        # high coverage must not make it auto-selectable
        self.assertIsNone(self.module.select_supported_hit([hit], min_coverage=0.8,
                                                           query_len=1000))

    def test_opposite_target_direction_is_conflicting(self):
        hit = self._hit("NC_C", self._hsp(1, 400, 1, 400, 400, 400)
                        + self._hsp(601, 1000, 899, 500, 390, 400))
        self.assertFalse(hit["collinearity"]["consistent_direction"])
        self.assertTrue(hit["conflicting_alignment"])
        self.assertIsNone(hit["identity"])

    def test_exact_duplicate_hsp_is_overlap_not_conflict(self):
        hsp = self._hsp(1, 500, 1, 500, 495, 500)
        hit = self._hit("NC_D", hsp + hsp)
        self.assertTrue(hit["ambiguous_alignment"])          # not aggregated (conservative)
        self.assertFalse(hit["conflicting_alignment"])        # a duplicate is not a contradiction
        self.assertEqual(hit["ambiguity_reasons"], ["overlapping_hsps"])
        self.assertEqual(hit["redundant_bases"], 500)

    def test_blocked_hit_keeps_its_raw_evidence(self):
        hit = self._hit("NC_E", self._hsp(1, 400, 1, 400, 400, 400)
                        + self._hsp(601, 1000, 9000, 9399, 390, 400))
        self.assertEqual(len(hit["hsps"]), 2)
        self.assertEqual(hit["hsps"][1]["hit_from"], 9000)
        self.assertIsNotNone(hit["identity_range"])
        self.assertGreater(hit["bitscore"], 0)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class PartialSilentFailureRegression(unittest.TestCase):
    """The original partial-detection bug, kept as a permanent regression case.

    The first implementation parsed ``str(location)`` with a regex written against
    the GenBank FILE syntax (``<1..100``).  Biopython's ``str()`` emits
    ``[<0:100](+)`` (colon, 0-based), so the regex matched nothing and
    ``feature_partial()`` returned ``(False, False)`` for every feature: the check
    ran, printed nothing, and silently turned every partial CDS into a complete
    one.  Guards kept in place:
      1. the position objects are authoritative (the string cannot cause this);
      2. the string fallback accepts BOTH syntaxes;
      3. a marker-free location must never invent partiality.
    """

    def setUp(self):
        self.module = load_module("annot_partial_regression", Path("scripts") / "annot_check.py")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def test_real_genbank_round_trip_detects_both_ends(self):
        from Bio import SeqIO
        path = self.dir / "partial.gb"
        write_gb_raw(path, "A" * 400, [
            {"location": "<1..100", "type": "CDS", "gene": "p5"},
            {"location": "200..>300", "type": "CDS", "gene": "p3"},
            {"location": "complement(<250..>350)", "type": "CDS", "gene": "both_rc"},
        ])
        record = SeqIO.read(path, "genbank")
        by_gene = {f.qualifiers["gene"][0]: f for f in record.features if f.type == "CDS"}
        self.assertEqual(self.module.feature_partial(by_gene["p5"]), (True, False))
        self.assertEqual(self.module.feature_partial(by_gene["p3"]), (False, True))
        self.assertEqual(self.module.feature_partial(by_gene["both_rc"]), (True, True))

    def test_fallback_accepts_both_print_and_file_syntax(self):
        class FakeLocation:
            def __init__(self, text, strand=1):
                self.parts = [SimpleNamespace(start=0, end=100)]   # no position objects
                self.strand = strand
                self._text = text

            def __str__(self):
                return self._text

        for text, expected in (("[<0:100](+)", (True, False)),   # Biopython print form (colon)
                               ("<1..100", (True, False)),        # GenBank file form (dots)
                               ("[<0:>100](+)", (True, True)),    # both ends
                               ("200..>300", (False, True))):
            with self.subTest(text=text):
                feature = SimpleNamespace(location=FakeLocation(text))
                self.assertEqual(self.module.feature_partial(feature), expected)

    def test_fallback_does_not_invent_partiality(self):
        class FakeLocation:
            def __init__(self, text):
                self.parts = [SimpleNamespace(start=0, end=100)]
                self.strand = 1
                self._text = text

            def __str__(self):
                return self._text

        feature = SimpleNamespace(location=FakeLocation("[0:100](+)"))
        self.assertEqual(self.module.feature_partial(feature), (False, False))

    def test_partiality_is_visible_in_the_report(self):
        sequence = "TTT" + "AAA" * 20 + "TAA"
        path = self.dir / "visible.gb"
        write_gb_raw(path, sequence + "A" * 40, [
            {"location": "<1..%d" % len(sequence), "type": "CDS", "gene": "cox1"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("PARTIAL_CDS_5P", result.stdout)


class BlastIdentityAggregationTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module("cox1_agg", Path("scripts") / "cox1_id.py")

    @staticmethod
    def _hsp(query_from, query_to, identity, align_len, bitscore=500):
        return ("<Hsp><Hsp_bit-score>%d</Hsp_bit-score>"
                "<Hsp_query-from>%d</Hsp_query-from><Hsp_query-to>%d</Hsp_query-to>"
                "<Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>%d</Hsp_hit-to>"
                "<Hsp_identity>%d</Hsp_identity><Hsp_align-len>%d</Hsp_align-len></Hsp>"
                % (bitscore, query_from, query_to, align_len, identity, align_len))

    def _xml(self, hits, query_len=1000):
        return ("<?xml version=\"1.0\"?><BlastOutput><BlastOutput_db>nt</BlastOutput_db>"
                "<BlastOutput_iterations><Iteration><Iteration_query-len>%d</Iteration_query-len>"
                "<Iteration_hits>%s</Iteration_hits></Iteration></BlastOutput_iterations>"
                "</BlastOutput>" % (query_len, "".join(hits)))

    @staticmethod
    def _hit(accession, hsps):
        return ("<Hit><Hit_accession>%s</Hit_accession><Hit_def>species %s</Hit_def>"
                "<Hit_len>15000</Hit_len><Hit_hsps>%s</Hit_hsps></Hit>"
                % (accession, accession, hsps))

    def test_overlapping_hsps_are_ambiguous_and_not_aggregated(self):
        hit = self._hit("NC_1", self._hsp(1, 500, 480, 500) + self._hsp(400, 800, 380, 400))
        best = self.module.parse_blast_xml(self._xml([hit]))["hits"][0]
        self.assertAlmostEqual(best["coverage"], 0.8, places=6)
        self.assertTrue(best["ambiguous_alignment"])
        self.assertIsNone(best["identity"])
        self.assertGreater(best["redundant_bases"], 0)
        self.assertIn("overlap", best["identity_method"])

    def test_non_overlapping_hsps_keep_a_weighted_identity(self):
        hit = self._hit("NC_2", self._hsp(1, 400, 400, 400) + self._hsp(601, 1000, 400, 400))
        best = self.module.parse_blast_xml(self._xml([hit]))["hits"][0]
        self.assertAlmostEqual(best["identity"], 100.0, places=6)
        self.assertFalse(best["ambiguous_alignment"])

    def test_raw_hsps_are_preserved(self):
        hit = self._hit("NC_3", self._hsp(1, 400, 400, 400) + self._hsp(601, 1000, 390, 400))
        best = self.module.parse_blast_xml(self._xml([hit]))["hits"][0]
        self.assertEqual(len(best["hsps"]), 2)
        self.assertEqual(best["hsps"][0]["query_from"], 1)
        self.assertEqual(best["hsps"][1]["query_to"], 1000)
        self.assertIn("bitscore", best["hsps"][0])

    def test_conflicting_overlaps_are_flagged(self):
        # same query span, different identity -> cannot be summarised
        hit = self._hit("NC_4", self._hsp(1, 500, 495, 500) + self._hsp(1, 500, 400, 500))
        best = self.module.parse_blast_xml(self._xml([hit]))["hits"][0]
        self.assertTrue(best["ambiguous_alignment"])
        self.assertIsNone(best["identity"])

    def test_cross_accession_candidates_are_retained(self):
        hits = ([self._hit("NC_5", self._hsp(1, 1000, 990, 1000))]
                + [self._hit("NC_6", self._hsp(1, 1000, 970, 1000))])
        parsed = self.module.parse_blast_xml(self._xml(hits))
        distinct = self.module.distinct_candidates(parsed["hits"])
        self.assertEqual(len(distinct), 2)
        self.assertEqual({item["accession"] for item in distinct}, {"NC_5", "NC_6"})

    def test_ambiguous_hit_is_not_selected_as_a_confident_candidate(self):
        hit = self._hit("NC_7", self._hsp(1, 500, 495, 500) + self._hsp(400, 800, 380, 400))
        best = self.module.parse_blast_xml(self._xml([hit]))["hits"][0]
        self.assertIsNone(self.module.select_supported_hit([best], min_coverage=0.8,
                                                           query_len=1000))


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Permanent failure reproductions for the PR1 independent review (DO NOT MERGE).

Every case below failed on `d3ca6f6` and is kept as a regression guard:

  P1-1  SearchInfo: NCBI answers ``FORMAT_OBJECT=SearchInfo`` with the QBlastInfo
        PLAIN TEXT block (``RID = ...`` / ``Status=WAITING``).  The parser only
        matched XML tags, so the remote flow could never obtain a RID and could
        never see READY.
  P1-2  The tool requests ``FORMAT_TYPE=XML2`` but the parser only understood the
        legacy ``-outfmt 5`` dialect: a real XML2 result containing hits was
        reported as ``no_match`` (exit 1) instead of a format failure (exit 3).
  P1-3  Query-disjoint but TARGET-overlapping HSPs were aggregated into one
        precise identity and auto-selected as the best candidate.
  P1-4  ``/transl_except`` accepted a wrong-frame 3 nt window and any amino-acid
        string, silently excusing a real internal stop.
  P1-5  The no-jsonschema fallback (the one CI runs) accepted cases the schema
        rejects: empty ``case_id`` and non-array ``modifications`` /
        ``validation`` / ``lessons_proposed``.
  P2-1  ``pos:complement(...)`` and compound-CDS exceptions were dropped; the
        former silently, without even a TRANSL_EXCEPT_UNEXPLAINED finding.
  P2-2  The CLI claimed the raw HSPs were "retained" but printed/persisted none.
  P2-3  A registry whose top level is a JSON array raised AttributeError.

The XML fixtures under tests/fixtures/ are verbatim ``blastn`` 2.12.0+ output, so
the production dialect is pinned by a versioned artefact rather than by a
hand-written sample that can drift from the real format.
"""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gb_fixtures import (  # noqa: E402
    ROOT, biopython_available, have, load_module, run_annot_check, run_script, write_gb_raw,
)

REASON = "合成测试: 单基因记录, 基因集差异已逐项确认"
FIXTURES = Path(__file__).resolve().parent / "fixtures"
XML2_FIXTURE = FIXTURES / "blast_xml2_outfmt16.xml"
XML5_FIXTURE = FIXTURES / "blast_xml_outfmt5.xml"

# ATG AAA AAA TAA AAA AAA TAA -> the only internal stop is codon 4 (record 10..12)
INTERNAL_STOP_CDS = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA"


def _revcomp(sequence):
    from Bio.Seq import Seq
    return str(Seq(sequence).reverse_complement())


def _hsp(query_from, query_to, hit_from, hit_to, identity, align_len, bitscore=500):
    return ("<Hsp><Hsp_bit-score>%d</Hsp_bit-score>"
            "<Hsp_query-from>%d</Hsp_query-from><Hsp_query-to>%d</Hsp_query-to>"
            "<Hsp_hit-from>%d</Hsp_hit-from><Hsp_hit-to>%d</Hsp_hit-to>"
            "<Hsp_identity>%d</Hsp_identity><Hsp_align-len>%d</Hsp_align-len></Hsp>"
            % (bitscore, query_from, query_to, hit_from, hit_to, identity, align_len))


def _blast_xml(accession, hsps, query_len=1000, hit_len=15000, description="species"):
    return ("<?xml version=\"1.0\"?><BlastOutput><BlastOutput_db>nt</BlastOutput_db>"
            "<BlastOutput_iterations><Iteration><Iteration_query-len>%d</Iteration_query-len>"
            "<Iteration_hits><Hit><Hit_accession>%s</Hit_accession><Hit_def>%s</Hit_def>"
            "<Hit_len>%d</Hit_len><Hit_hsps>%s</Hit_hsps></Hit></Iteration_hits>"
            "</Iteration></BlastOutput_iterations></BlastOutput>"
            % (query_len, accession, description, hit_len, hsps))


class SearchInfoProtocolTests(unittest.TestCase):
    """P1-1: the format NCBI actually answers with."""

    def setUp(self):
        self.module = load_module("cox1_searchinfo", Path("scripts") / "cox1_id.py")

    def test_qblast_plain_text_is_parsed(self):
        text = ("QBlastInfoBegin\n"
                "    RID = ABC123DEF\n"
                "    RTOE = 17\n"
                "Status=WAITING\n"
                "QBlastInfoEnd")
        self.assertEqual(self.module.parse_search_info(text), ("WAITING", "ABC123DEF"))

    def test_ready_status_in_plain_text(self):
        self.assertEqual(self.module.parse_search_info("RID = XYZ789\nStatus=READY"),
                         ("READY", "XYZ789"))

    def test_spaces_around_the_equals_sign(self):
        self.assertEqual(self.module.parse_search_info("RID = XYZ789\nStatus = READY"),
                         ("READY", "XYZ789"))

    def test_put_response_without_a_status_still_yields_the_rid(self):
        self.assertEqual(self.module.parse_search_info("RID = ABC123\nRTOE = 11"),
                         (None, "ABC123"))

    def test_xml_search_info_is_still_accepted(self):
        text = "<BlastSearchInfo><Status>READY</Status><RID>ABC123</RID></BlastSearchInfo>"
        self.assertEqual(self.module.parse_search_info(text), ("READY", "ABC123"))

    def test_unknown_payload_is_not_mistaken_for_a_waiting_search(self):
        payload = "<html><body>Error: service temporarily unavailable</body></html>"
        self.assertEqual(self.module.parse_search_info(payload), (None, None))


class BlastXmlDialectTests(unittest.TestCase):
    """P1-2: the dialect the tool requests (XML2) must be readable."""

    def setUp(self):
        self.module = load_module("cox1_dialect", Path("scripts") / "cox1_id.py")

    def test_real_xml2_outfmt16_matches_are_parsed(self):
        parsed = self.module.parse_blast_xml(XML2_FIXTURE.read_text(encoding="utf-8"))
        self.assertEqual(parsed["query_len"], 800)
        self.assertEqual(len(parsed["hits"]), 1)
        hit = parsed["hits"][0]
        self.assertEqual(hit["hsp_count"], 1)
        self.assertEqual(hit["hit_len"], 1200)
        self.assertAlmostEqual(hit["coverage"], 1.0, places=6)
        self.assertAlmostEqual(hit["identity"], 100.0, places=6)
        self.assertTrue(hit["accession"])
        self.assertIn("mitochondrion", hit["description"])
        self.assertIsNotNone(hit["hsps"][0]["hit_from"])

    def test_legacy_outfmt5_fixture_still_parses(self):
        parsed = self.module.parse_blast_xml(XML5_FIXTURE.read_text(encoding="utf-8"))
        self.assertEqual(parsed["query_len"], 800)
        self.assertEqual(len(parsed["hits"]), 1)
        self.assertEqual(parsed["hits"][0]["hsp_count"], 1)

    def test_unrecognised_root_is_a_format_failure(self):
        with self.assertRaises(ValueError):
            self.module.parse_blast_xml("<?xml version='1.0'?><foo><bar>1</bar></foo>")

    def test_recognised_root_without_known_fields_is_a_format_failure(self):
        with self.assertRaises(ValueError):
            self.module.parse_blast_xml(
                "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\"/>")

    def test_recognised_root_with_unreadable_hits_is_a_format_failure(self):
        xml = ("<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\">"
               "<query-len>10</query-len><Hit><num>1</num>"
               "<description><HitDescr><id>x</id></HitDescr></description></Hit></BlastXML2>")
        with self.assertRaises(ValueError):
            self.module.parse_blast_xml(xml)

    def test_legitimate_no_hit_is_not_a_format_failure(self):
        xml = ("<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\"><query-len>500</query-len>"
               "<stat><Statistics><db-len>1</db-len></Statistics></stat></BlastXML2>")
        parsed = self.module.parse_blast_xml(xml)
        self.assertEqual(parsed["hits"], [])
        self.assertEqual(parsed["query_len"], 500)


class Cox1CliExitBoundaryTests(unittest.TestCase):
    """P1-2 / P2-2: exit-code boundary and evidence persistence at the CLI."""

    def setUp(self):
        self.module = load_module("cox1_cli", Path("scripts") / "cox1_id.py")
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name)
        self.fasta = self.dir / "genome.fa"
        self.fasta.write_text(">g\n%s\n" % ("ACGT" * 300), encoding="utf-8")

    def _run(self, payload, extra=()):
        module = self.module
        saved = (module.request_search, module.poll_search_info,
                 module.fetch_results_xml, module.time)
        module.request_search = lambda query, max_results=10: "RID = RID1\nRTOE = 0"
        module.poll_search_info = lambda rid: "Status=READY"
        module.fetch_results_xml = lambda rid, max_results=10: payload
        module.time = SimpleNamespace(sleep=lambda seconds: None)
        argv = sys.argv
        sys.argv = ["cox1_id.py", str(self.fasta), "--allow-public-upload",
                    "--coords", "1,1200"] + list(extra)
        try:
            with contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    module.main()
            return caught.exception.code
        finally:
            (module.request_search, module.poll_search_info,
             module.fetch_results_xml, module.time) = saved
            sys.argv = argv

    def test_format_failure_is_exit_3_not_no_match(self):
        self.assertEqual(self._run("<html>not a BLAST report</html>"), 3)

    def test_real_xml2_result_reaches_a_verdict(self):
        self.assertEqual(self._run(XML2_FIXTURE.read_text(encoding="utf-8")), 0)

    def test_legitimate_no_hit_stays_exit_1(self):
        payload = ("<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\">"
                   "<query-len>1200</query-len></BlastXML2>")
        self.assertEqual(self._run(payload), 1)

    def test_output_json_persists_every_hsp(self):
        out = self.dir / "results.json"
        code = self._run(XML2_FIXTURE.read_text(encoding="utf-8"),
                         extra=["--output-json", str(out)])
        self.assertEqual(code, 0)
        data = json.loads(out.read_text(encoding="utf-8"))
        self.assertEqual(data["verdict"], "provisional_candidate")
        candidate = data["candidates"][0]
        self.assertEqual(data["selection"], candidate["accession"])
        self.assertEqual(candidate["hsps"][0]["query_from"], 1)
        self.assertIn("bitscore", candidate["hsps"][0])

    def test_blocked_candidate_prints_its_raw_hsp_table(self):
        payload = _blast_xml("NC_BLOCK", _hsp(1, 400, 1000, 1399, 400, 400)
                             + _hsp(401, 800, 1200, 1599, 400, 400))
        module = self.module
        saved = (module.request_search, module.poll_search_info,
                 module.fetch_results_xml, module.time)
        module.request_search = lambda query, max_results=10: "RID = RID1"
        module.poll_search_info = lambda rid: "Status=READY"
        module.fetch_results_xml = lambda rid, max_results=10: payload
        module.time = SimpleNamespace(sleep=lambda seconds: None)
        argv = sys.argv
        sys.argv = ["cox1_id.py", str(self.fasta), "--allow-public-upload", "--coords", "1,1200"]
        captured = io.StringIO()
        try:
            with contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(captured):
                with self.assertRaises(SystemExit) as caught:
                    module.main()
        finally:
            (module.request_search, module.poll_search_info,
             module.fetch_results_xml, module.time) = saved
            sys.argv = argv
        self.assertEqual(caught.exception.code, 1)
        self.assertIn("subject_overlap", captured.getvalue())
        self.assertIn("HSP q1..400 -> s1000..1399", captured.getvalue())


class SubjectOverlapTests(unittest.TestCase):
    """P1-3: reused target positions are not two independent alignments."""

    def setUp(self):
        self.module = load_module("cox1_subject_overlap", Path("scripts") / "cox1_id.py")

    def _hit(self, accession, hsps):
        return self.module.parse_blast_xml(_blast_xml(accession, hsps))["hits"][0]

    def test_plus_strand_subject_overlap_is_a_blocker(self):
        hit = self._hit("NC_OVERLAP", _hsp(1, 400, 1000, 1399, 400, 400)
                        + _hsp(401, 800, 1200, 1599, 400, 400))
        self.assertTrue(hit["subject_overlap"])
        self.assertEqual(hit["subject_overlap_bases"], 200)
        self.assertIn("subject_overlap", hit["ambiguity_reasons"])
        self.assertIsNone(hit["identity"])
        self.assertIsNone(self.module.select_supported_hit([hit], min_coverage=0.8,
                                                           query_len=1000))

    def test_minus_strand_subject_overlap_is_detected(self):
        # monotonic in the target's own (descending) direction, yet the two
        # subject intervals still share 200 bp
        hit = self._hit("NC_MINUS_OVERLAP", _hsp(1, 400, 1599, 1200, 400, 400)
                        + _hsp(401, 800, 1399, 1000, 400, 400))
        self.assertEqual(hit["collinearity"]["subject_strand"], "-")
        self.assertTrue(hit["collinearity"]["collinear"])
        self.assertTrue(hit["subject_overlap"])
        self.assertEqual(hit["subject_overlap_bases"], 200)
        self.assertIsNone(hit["identity"])

    def test_duplicate_subject_interval_is_a_blocker(self):
        hit = self._hit("NC_DUP", _hsp(1, 400, 500, 899, 400, 400)
                        + _hsp(601, 1000, 500, 899, 390, 400))
        self.assertEqual(hit["redundant_bases"], 0)          # query side is clean
        self.assertTrue(hit["subject_overlap"])
        self.assertEqual(hit["subject_overlap_bases"], 400)
        self.assertIsNone(hit["identity"])

    def test_adjacent_target_blocks_are_not_overlap(self):
        # the real contiguous minus-strand shape, locked by BlastCoordinateConventionTests
        hit = self._hit("NC_MINUS_OK", _hsp(1, 400, 1800, 1401, 400, 400)
                        + _hsp(501, 900, 1400, 1001, 400, 400))
        self.assertFalse(hit["subject_overlap"])
        self.assertEqual(hit["subject_overlap_bases"], 0)
        self.assertTrue(hit["collinearity"]["collinear"])
        self.assertFalse(hit["ambiguous_alignment"])
        self.assertIsNotNone(hit["identity"])

    def test_touching_target_blocks_are_not_overlap(self):
        hit = self._hit("NC_TOUCH", _hsp(1, 400, 1, 400, 400, 400)
                        + _hsp(401, 800, 401, 800, 400, 400))
        self.assertEqual(hit["subject_overlap_bases"], 0)
        self.assertIsNotNone(hit["identity"])
        self.assertEqual(hit["coverage"], 0.8)

    def test_exact_duplicate_hsp_is_both_query_and_target_overlap(self):
        hsp = _hsp(1, 500, 1, 500, 495, 500)
        hit = self._hit("NC_EXACT", hsp + hsp)
        self.assertIn("overlapping_hsps", hit["ambiguity_reasons"])
        self.assertIn("subject_overlap", hit["ambiguity_reasons"])
        self.assertFalse(hit["conflicting_alignment"])


class SubjectOverlapViaRealBlastnTests(unittest.TestCase):
    """P1-3 on real output: a duplicated target region must not be aggregated."""

    @unittest.skipUnless(have("blastn"), "blastn is required")
    def test_real_duplicated_subject_region_is_blocked(self):
        import random
        import subprocess
        module = load_module("cox1_real_overlap", Path("scripts") / "cox1_id.py")
        rng = random.Random(17)
        subject = "".join(rng.choice("ACGT") for _ in range(2000))
        # two query halves that both target the SAME 400 bp stretch
        query = subject[600:1000] + "N" * 100 + subject[600:1000]
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            (directory / "s.fa").write_text(">s\n%s\n" % subject, encoding="utf-8")
            (directory / "q.fa").write_text(">q\n%s\n" % query, encoding="utf-8")
            result = subprocess.run(
                ["blastn", "-query", str(directory / "q.fa"),
                 "-subject", str(directory / "s.fa"), "-outfmt", "16", "-dust", "no"],
                capture_output=True, text=True, check=True, timeout=120)
        parsed = module.parse_blast_xml(result.stdout)
        self.assertTrue(parsed["hits"], "expected a hit")
        hit = parsed["hits"][0]
        self.assertTrue(hit["subject_overlap"], hit["hsps"])
        self.assertIsNone(hit["identity"])
        self.assertIn("subject_overlap", hit["ambiguity_reasons"])


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class TranslExceptExactCodonTests(unittest.TestCase):
    """P1-4 / P2-1: the exception must name the exact codon of a real internal stop."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name)

    def _run(self, name, sequence, transl_except, location="1..%d", registry=None,
             taxon="Lepidoptera", gene="cox1"):
        path = self.dir / name
        if "%d" in location:
            location = location % len(sequence)
        spec = {"location": location, "type": "CDS"}
        if gene:
            spec["gene"] = gene
        if transl_except is not None:
            spec["transl_except"] = transl_except
        write_gb_raw(path, sequence + "A" * 40, [spec])
        extra = []
        if registry:
            extra += ["--exception-registry", _write_registry(self.dir, name, *registry)]
        if taxon:
            extra += ["--taxon", taxon]
        return run_annot_check(path, "--allow-atypical", REASON, *extra)

    def test_exact_codon_with_evidence_is_validated(self):
        result = self._run("ok.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp)",
                           registry=[_transl_record()])
        self.assertIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("TRANSL_EXCEPT_VALIDATED", result.stdout)
        self.assertNotIn("[ERROR]", result.stdout)

    def test_exact_codon_without_evidence_stays_unverified(self):
        # a position/reading-frame match is a syntax fact, not a biological claim:
        # without an audited record the internal stop keeps its ERROR
        result = self._run("ok-bare.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp)")
        self.assertIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("TRANSL_EXCEPT_DECLARED_UNVERIFIED", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_VALIDATED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_missing_taxon_flag_keeps_the_stop_as_error(self):
        # P1-1: a record whose taxon cannot be tied to this sample proves nothing
        result = self._run("no-taxon.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp)",
                           registry=[_transl_record()], taxon=None)
        self.assertIn("TRANSL_EXCEPT_DECLARED_UNVERIFIED", result.stdout)
        self.assertIn("--taxon", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_VALIDATED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_cds_without_a_gene_identity_is_never_validated(self):
        # the record itself is valid, but a CDS whose /gene and /product both
        # canonicalise to '' must not have gene-keyed evidence applied to it
        result = self._run("unnamed.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp)",
                           registry=[_transl_record()], gene=None)
        self.assertIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("TRANSL_EXCEPT_DECLARED_UNVERIFIED", result.stdout)
        self.assertIn("gene", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_VALIDATED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_wrong_frame_three_nt_window_is_not_accepted(self):
        result = self._run("wrong-frame.gb", INTERNAL_STOP_CDS, "(pos:11..13,aa:Trp)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)
        self.assertIn("内部终止", result.stdout)

    def test_offset_window_ending_on_the_stop_is_not_accepted(self):
        result = self._run("shifted.gb", INTERNAL_STOP_CDS, "(pos:9..11,aa:Trp)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_invalid_amino_acid_is_not_accepted(self):
        result = self._run("bad-aa.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Foo)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("aa:Foo", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_minus_strand_complement_position_is_validated(self):
        result = self._run("minus-ok.gb", _revcomp(INTERNAL_STOP_CDS),
                           "(pos:complement(10..12),aa:Trp)", location="complement(1..%d)",
                           registry=[_transl_record()])
        self.assertIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("TRANSL_EXCEPT_VALIDATED", result.stdout)
        self.assertNotIn("[ERROR]", result.stdout)

    def test_minus_strand_complement_at_the_wrong_codon_is_not_accepted(self):
        result = self._run("minus-bad.gb", _revcomp(INTERNAL_STOP_CDS),
                           "(pos:complement(11..13),aa:Trp)", location="complement(1..%d)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_compound_cds_exception_inside_a_part_is_validated(self):
        # CDS = join(1..6, 57..62): coding order is 1..6 then 57..62, so the codons
        # are ATG (1..3), AAA (4..6), TAA (57..59, internal stop) and TAA (60..62,
        # terminal stop).
        sequence = list("A" * 200)
        sequence[0:6] = list("ATGAAA")
        sequence[56:62] = list("TAATAA")
        sequence = "".join(sequence)
        result = self._run("join-ok.gb", sequence, "(pos:57..59,aa:Trp)",
                           location="join(1..6,57..62)",
                           registry=[_transl_record(pos="57..59")])
        self.assertIn("TRANSL_EXCEPT_VALIDATED", result.stdout)
        self.assertNotIn("[ERROR]", result.stdout)

    def test_compound_cds_exception_outside_its_codon_is_not_accepted(self):
        sequence = list("A" * 200)
        sequence[0:6] = list("ATGAAA")
        sequence[56:62] = list("TAATAA")
        sequence = "".join(sequence)
        result = self._run("join-bad.gb", sequence, "(pos:58..60,aa:Trp)",
                           location="join(1..6,57..62)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_unparsable_qualifier_is_reported_not_dropped(self):
        result = self._run("unparsable.gb", INTERNAL_STOP_CDS, "(pos:exon1,aa:Trp)")
        self.assertIn("TRANSL_EXCEPT_UNPARSED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_coding_position_list_matches_feature_extract(self):
        """The codon map must agree with how Biopython actually extracts the CDS."""
        import random
        from Bio import SeqIO
        from Bio.Seq import Seq
        module = load_module("annot_coding_positions", Path("scripts") / "annot_check.py")
        rng = random.Random(23)
        sequence = "".join(rng.choice("ACGT") for _ in range(400))
        cases = {
            "plus.gb": "1..30",
            "minus.gb": "complement(1..30)",
            "join.gb": "join(1..12,101..121)",
            "join_minus.gb": "complement(join(1..12,101..121))",
        }
        for name, location in cases.items():
            with self.subTest(location=location):
                path = self.dir / name
                write_gb_raw(path, sequence, [{"location": location, "type": "CDS"}])
                record = SeqIO.read(path, "genbank")
                feature = record.features[0]
                positions = module.coding_position_list(feature)
                extracted = str(feature.extract(sequence))
                translated = "".join(str(Seq(extracted[i:i + 3]).translate())
                                     for i in range(0, len(extracted) - 2, 3))
                complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
                strand = feature.location.strand or 0
                coding = "".join(complement.get(sequence[p - 1], sequence[p - 1])
                                  if strand < 0 else sequence[p - 1]
                                  for p in positions)
                rebuilt = "".join(str(Seq(coding[i:i + 3]).translate())
                                  for i in range(0, len(coding) - 2, 3))
                self.assertEqual(len(positions), len(extracted))
                self.assertEqual(rebuilt, translated)


class ExceptionRegistryTypeTests(unittest.TestCase):
    """P2-3: malformed registries fail through the CLI contract, not a traceback."""

    def setUp(self):
        self.module = load_module("annot_registry_types", Path("scripts") / "annot_check.py")
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name)

    def _load(self, payload):
        path = self.dir / "registry.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                self.module.load_exception_registry(str(path))
        return caught.exception.code

    def test_top_level_array_is_a_controlled_error(self):
        self.assertEqual(self._load([]), 1)

    def test_top_level_string_is_a_controlled_error(self):
        self.assertEqual(self._load("not-a-registry"), 1)

    def test_exceptions_not_an_array_is_a_controlled_error(self):
        self.assertEqual(self._load({"exceptions": {"gene": "cox1"}}), 1)

    def test_entry_not_an_object_is_a_controlled_error(self):
        self.assertEqual(self._load({"exceptions": ["cox1"]}), 1)

    def test_valid_registry_still_loads(self):
        path = self.dir / "ok.json"
        path.write_text(json.dumps({"exceptions": [
            {"gene": "cox1", "codon": "CGA", "taxon": "Lepidoptera",
             "source": "DOI 10.1000/x", "rationale": "documented"}]}), encoding="utf-8")
        registry = self.module.load_exception_registry(str(path))
        self.assertIn(("cox1", "CGA"), registry['start'])
        self.assertEqual(registry['transl_except'], {})

    def test_valid_transl_except_record_lands_in_its_own_section(self):
        path = self.dir / "transl.json"
        path.write_text(json.dumps({"exceptions": [
            {"gene": "cox1", "codon": "TAA", "amino_acid": "Trp", "transl_table": 5,
             "taxon": "Lepidoptera", "source": "DOI 10.1000/x", "rationale": "r",
             "pos": "10..12"}]}), encoding="utf-8")
        registry = self.module.load_exception_registry(str(path))
        self.assertEqual(registry['start'], {})
        self.assertIn(("cox1", "TAA", "trp"), registry['transl_except'])


class SchemaKeywordCoverageTests(unittest.TestCase):
    """P1-5: every declared top-level type/const keyword is enforced by the fallback."""

    def setUp(self):
        self.module = load_module("experience_keyword_coverage", Path("tools") / "experience.py")

    @staticmethod
    def _base():
        return {
            "schema_version": "2.0", "case_id": "case-1",
            "inputs": [{"role": "assembly_fasta", "path": "/tmp/x.fa", "sha256": "a" * 64}],
            "issue": {"type": "internal_stop", "user_observation": "x"},
            "hypotheses": [{"id": "H1", "explanation": "e", "support": [], "against": [],
                            "unknown": []}],
            "decision": {"status": "RESOLVED", "confidence": "moderate", "rationale": "r"},
            "anomalies": [], "events_file": "events.jsonl",
        }

    def test_every_top_level_keyword_is_enforced(self):
        violations = {
            "schema_version": "9.9", "case_id": 5, "case_id_empty": "", "taxon": "Diptera",
            "inputs": "x", "issue": "x", "hypotheses": "x", "events_file": 1,
            "decision": "x", "anomalies": "x", "modifications": "x",
            "validation": "x", "lessons_proposed": "x", "case_type": "anomalous",
        }
        for key, bad in violations.items():
            with self.subTest(field=key):
                case = self._base()
                case[key.replace("_empty", "")] = bad
                self.assertTrue(self.module._case_errors_without_jsonschema(case),
                                "%s=%r was accepted by the fallback" % (key, bad))

    def test_a_legal_case_still_passes_the_fallback(self):
        self.assertEqual(self.module._case_errors_without_jsonschema(self._base()), [])

    def test_legal_array_valued_optionals_are_accepted(self):
        case = self._base()
        case.update({"modifications": [{"x": 1}], "validation": ["a"],
                     "lessons_proposed": ["b"]})
        self.assertEqual(self.module._case_errors_without_jsonschema(case), [])


def _hsp_missing(fields):
    """An <Hsp> with the named fields deliberately omitted."""
    values = {'query-from': 1, 'query-to': 1000, 'hit-from': 100, 'hit-to': 1099,
              'identity': 990, 'align-len': 1000, 'bit-score': 900.0}
    for field in fields:
        values.pop(field)
    return ("<Hsp><Hsp_query-from>%s</Hsp_query-from><Hsp_query-to>%s</Hsp_query-to>"
            "<Hsp_hit-from>%s</Hsp_hit-from><Hsp_hit-to>%s</Hsp_hit-to>"
            "<Hsp_identity>%s</Hsp_identity><Hsp_align-len>%s</Hsp_align-len>"
            "<Hsp_bit-score>%s</Hsp_bit-score></Hsp>"
            % tuple('' if key not in values else values[key]
                    for key in ('query-from', 'query-to', 'hit-from', 'hit-to',
                                'identity', 'align-len', 'bit-score')))


def _write_registry(directory, name, *entries):
    """Write a temporary --exception-registry file and return its path."""
    path = Path(directory) / (name + ".registry.json")
    path.write_text(json.dumps({"exceptions": list(entries)}, ensure_ascii=False),
                    encoding="utf-8")
    return str(path)


def _transl_record(pos="10..12", codon="TAA", amino_acid="Trp", transl_table=5,
                   taxon="Lepidoptera", source="DOI 10.1000/example",
                   rationale="documented exception", **extra):
    """A complete, auditable, SITE-BOUND /transl_except registry entry.

    ``pos`` is the genomic position set of the single stop codon this record covers
    (``None`` omits it, which is only valid together with ``scope="gene_wide"``).
    """
    record = {"gene": "cox1", "codon": codon, "amino_acid": amino_acid,
              "transl_table": transl_table, "taxon": taxon, "source": source,
              "rationale": rationale}
    if pos is not None:
        record["pos"] = pos
    record.update(extra)
    return record


class HspFieldValidationTests(unittest.TestCase):
    """P1-1: a half-readable HSP must be a format failure, not a candidate.

    Without the subject coordinates the direction, collinearity and target-reuse
    checks cannot run.  The old code skipped unusable HSPs and let the remaining
    "single" HSP take the fast path with collinear=True and no blocker, so a
    truncated response produced identity=99% and exit code 0.
    """

    def setUp(self):
        self.module = load_module("cox1_hsp_fields", Path("scripts") / "cox1_id.py")

    def _parse(self, hsp):
        return self.module.parse_blast_xml(_blast_xml("NC_X", hsp))

    def test_valid_hsp_still_parses(self):
        hit = self._parse(_hsp(1, 1000, 100, 1099, 990, 1000))["hits"][0]
        self.assertEqual(hit["hsp_count"], 1)
        self.assertEqual(hit["hsps"][0]["hit_from"], 100)
        self.assertIsNotNone(hit["identity"])

    def test_missing_subject_coordinates_is_a_format_failure(self):
        for field in ("hit-from", "hit-to"):
            with self.subTest(field=field):
                with self.assertRaises(ValueError) as caught:
                    self._parse(_hsp_missing([field]))
                self.assertIn(field, str(caught.exception))

    def test_missing_query_coordinates_is_a_format_failure(self):
        with self.assertRaises(ValueError):
            self._parse(_hsp_missing(["query-to"]))

    def test_missing_counts_is_a_format_failure(self):
        for field in ("identity", "align-len", "bit-score"):
            with self.subTest(field=field):
                with self.assertRaises(ValueError):
                    self._parse(_hsp_missing([field]))

    def test_non_numeric_field_is_a_format_failure(self):
        hsp = ("<Hsp><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>1000</Hsp_query-to>"
               "<Hsp_hit-from>100</Hsp_hit-from><Hsp_hit-to>1099</Hsp_hit-to>"
               "<Hsp_identity>990</Hsp_identity><Hsp_align-len>ABC</Hsp_align-len>"
               "<Hsp_bit-score>900</Hsp_bit-score></Hsp>")
        with self.assertRaises(ValueError):
            self._parse(hsp)

    def test_zero_coordinate_is_a_format_failure(self):
        with self.assertRaises(ValueError):
            self._parse(_hsp(0, 1000, 100, 1099, 990, 1000))

    def test_identity_above_align_len_is_a_format_failure(self):
        with self.assertRaises(ValueError):
            self._parse(_hsp(1, 500, 100, 599, 900, 500))

    def test_zero_align_len_is_a_format_failure(self):
        with self.assertRaises(ValueError):
            self._parse(_hsp(1, 500, 100, 599, 0, 0))

    def test_half_readable_hit_never_reaches_auto_selection(self):
        with self.assertRaises(ValueError):
            self.module.parse_blast_xml(
                _blast_xml("NC_MISSING", _hsp_missing(["hit-from", "hit-to"])))


class HspFieldValidationCliTests(unittest.TestCase):
    """P1-1 at the CLI: the half-readable response must exit 3, never 0."""

    def setUp(self):
        self.module = load_module("cox1_hsp_fields_cli", Path("scripts") / "cox1_id.py")
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.fasta = Path(temp.name) / "genome.fa"
        self.fasta.write_text(">g\n%s\n" % ("ACGT" * 300), encoding="utf-8")

    def _run(self, payload):
        module = self.module
        saved = (module.request_search, module.poll_search_info,
                 module.fetch_results_xml, module.time)
        module.request_search = lambda query, max_results=10: "RID = RID1"
        module.poll_search_info = lambda rid: "Status=READY"
        module.fetch_results_xml = lambda rid, max_results=10: payload
        module.time = SimpleNamespace(sleep=lambda seconds: None)
        argv = sys.argv
        sys.argv = ["cox1_id.py", str(self.fasta), "--allow-public-upload", "--coords", "1,1200"]
        try:
            with contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    module.main()
            return caught.exception.code
        finally:
            (module.request_search, module.poll_search_info,
             module.fetch_results_xml, module.time) = saved
            sys.argv = argv

    def test_missing_subject_coordinates_exits_3(self):
        payload = _blast_xml("NC_MISSING", _hsp_missing(["hit-from", "hit-to"]))
        self.assertEqual(self._run(payload), 3)

    def test_missing_bit_score_exits_3(self):
        payload = _blast_xml("NC_NOSCORE", _hsp_missing(["bit-score"]))
        self.assertEqual(self._run(payload), 3)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class TranslExceptSemanticsTests(unittest.TestCase):
    """P1-2 / P2-1: the declaration's strand, aa token and grammar must all hold."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name)

    def _run(self, name, sequence, transl_except, location="1..%d", registry=None, taxon="Lepidoptera"):
        path = self.dir / name
        if "%d" in location:
            location = location % len(sequence)
        write_gb_raw(path, sequence + "A" * 40, [
            {"location": location, "type": "CDS", "gene": "cox1",
             "transl_except": transl_except}])
        extra = []
        if registry:
            extra += ["--exception-registry", _write_registry(self.dir, name, *registry)]
        if taxon:
            extra += ["--taxon", taxon]
        return run_annot_check(path, "--allow-atypical", REASON, *extra)

    def test_aa_term_does_not_explain_an_internal_stop(self):
        result = self._run("term.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:TERM)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_bare_position_on_a_minus_strand_cds_is_not_accepted(self):
        result = self._run("bare-minus.gb", _revcomp(INTERNAL_STOP_CDS),
                           "(pos:10..12,aa:Trp)", location="complement(1..%d)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_complement_position_on_a_plus_strand_cds_is_not_accepted(self):
        result = self._run("comp-plus.gb", INTERNAL_STOP_CDS,
                           "(pos:complement(10..12),aa:Trp)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_trailing_comma_is_unparsed_and_does_not_explain(self):
        result = self._run("comma.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp),")
        self.assertIn("TRANSL_EXCEPT_UNPARSED", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_mixed_strand_compound_position_is_not_accepted(self):
        # join(complement(10..11),12) OR-ed its leaves into "minus", hiding the
        # plus-strand leaf and letting the declaration clear a real stop
        result = self._run("mixed.gb", _revcomp(INTERNAL_STOP_CDS),
                           "(pos:join(complement(10..11),12),aa:Trp)",
                           location="complement(1..%d)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("不同链方向", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_mixed_strand_compound_position_ignores_a_valid_registry(self):
        result = self._run("mixed-reg.gb", _revcomp(INTERNAL_STOP_CDS),
                           "(pos:join(complement(10..11),12),aa:Trp)",
                           location="complement(1..%d)", registry=[_transl_record()])
        self.assertIn("[ERROR]", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_VALIDATED", result.stdout)

    def test_trailing_garbage_is_unparsed_and_does_not_explain(self):
        result = self._run("junk.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp),BROKEN")
        self.assertIn("TRANSL_EXCEPT_UNPARSED", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_missing_pos_token_is_reported_as_unparsed(self):
        result = self._run("nopos.gb", INTERNAL_STOP_CDS, "(position:10..12,aa:Trp)")
        self.assertIn("TRANSL_EXCEPT_UNPARSED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_unbalanced_parentheses_are_unparsed(self):
        result = self._run("unbalanced.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp")
        self.assertIn("TRANSL_EXCEPT_UNPARSED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_fuzzy_location_is_unparsed(self):
        result = self._run("fuzzy.gb", INTERNAL_STOP_CDS, "(pos:<10..12,aa:Trp)")
        self.assertIn("TRANSL_EXCEPT_UNPARSED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_cross_part_codon_is_validated_with_evidence(self):
        # CDS join(1..5,101..107): codons are {1,2,3} {4,5,101} {102,103,104} {105,106,107}
        sequence = list("A" * 200)
        sequence[0:5] = list("ATGTA")
        sequence[100:107] = list("AATATAA")
        sequence = "".join(sequence)
        result = self._run("join-codon.gb", sequence,
                           "(pos:join(4..5,101),aa:Leu)", location="join(1..5,101..107)",
                           registry=[_transl_record(pos="join(4..5,101)", amino_acid="Leu")])
        self.assertIn("TRANSL_EXCEPT_VALIDATED", result.stdout)
        self.assertNotIn("[ERROR]", result.stdout)

    def test_cross_part_codon_without_evidence_stays_unverified(self):
        sequence = list("A" * 200)
        sequence[0:5] = list("ATGTA")
        sequence[100:107] = list("AATATAA")
        sequence = "".join(sequence)
        result = self._run("join-bare.gb", sequence,
                           "(pos:join(4..5,101),aa:Leu)", location="join(1..5,101..107)")
        self.assertIn("TRANSL_EXCEPT_DECLARED_UNVERIFIED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_cross_part_codon_with_wrong_positions_is_not_accepted(self):
        sequence = list("A" * 200)
        sequence[0:5] = list("ATGTA")
        sequence[100:107] = list("AATATAA")
        sequence = "".join(sequence)
        result = self._run("join-wrong.gb", sequence,
                           "(pos:join(4..5,102),aa:Leu)", location="join(1..5,101..107)")
        self.assertIn("TRANSL_EXCEPT_UNEXPLAINED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_one_record_covers_only_its_own_site(self):
        # P1-2: two identical TAA stops, one site-bound record -> only codon 4 is
        # validated; the second stop keeps its ERROR
        sequence = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA"
        result = self._run("one-site.gb", sequence,
                           "(pos:10..12,aa:Trp),(pos:19..21,aa:Trp)",
                           registry=[_transl_record(pos="10..12")])
        self.assertEqual(result.stdout.count("TRANSL_EXCEPT_VALIDATED"), 1)
        self.assertEqual(result.stdout.count("TRANSL_EXCEPT_DECLARED_UNVERIFIED"), 1)
        self.assertIn("[ERROR]", result.stdout)

    def test_codon_index_binding_covers_only_that_codon(self):
        sequence = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA"
        result = self._run("one-index.gb", sequence,
                           "(pos:10..12,aa:Trp),(pos:19..21,aa:Trp)",
                           registry=[_transl_record(pos=None, codon_index=4)])
        self.assertEqual(result.stdout.count("TRANSL_EXCEPT_VALIDATED"), 1)
        self.assertEqual(result.stdout.count("TRANSL_EXCEPT_DECLARED_UNVERIFIED"), 1)
        self.assertIn("[ERROR]", result.stdout)

    def test_explicit_gene_wide_scope_covers_multiple_sites(self):
        # reusing one record across sites requires an EXPLICIT, auditable scope
        sequence = "ATG" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA" + "AAA" * 2 + "TAA"
        result = self._run("gene-wide.gb", sequence,
                           "(pos:10..12,aa:Trp),(pos:19..21,aa:Trp)",
                           registry=[_transl_record(pos=None, scope="gene_wide")])
        self.assertEqual(result.stdout.count("TRANSL_EXCEPT_VALIDATED"), 2)
        self.assertIn("scope=gene_wide", result.stdout)
        self.assertNotIn("[ERROR]", result.stdout)

    def test_matching_position_without_evidence_is_not_accepted(self):
        # A legal but unproven amino acid (TAA -> Gln is not a known reassignment)
        # must not be promoted to an accepted exception: position/reading-frame
        # agreement is a syntax fact only, so the internal stop keeps its ERROR.
        result = self._run("wrong-aa.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Gln)")
        self.assertIn("TRANSL_EXCEPT_MATCHED", result.stdout)
        self.assertIn("TRANSL_EXCEPT_DECLARED_UNVERIFIED", result.stdout)
        self.assertNotIn("TRANSL_EXCEPT_VALIDATED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_incomplete_registry_record_stays_unverified(self):
        record = _transl_record()
        record.pop("source")
        result = self._run("incomplete.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp)",
                           registry=[record])
        self.assertIn("TRANSL_EXCEPT_DECLARED_UNVERIFIED", result.stdout)
        self.assertIn("source", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_registry_table_mismatch_stays_unverified(self):
        result = self._run("table-mismatch.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp)",
                           registry=[_transl_record(transl_table=2)])
        self.assertIn("TRANSL_EXCEPT_DECLARED_UNVERIFIED", result.stdout)
        self.assertIn("transl_table", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_registry_taxon_mismatch_stays_unverified(self):
        result = self._run("taxon-mismatch.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Trp)",
                           registry=[_transl_record(taxon="Lepidoptera")], taxon="Diptera")
        self.assertIn("TRANSL_EXCEPT_DECLARED_UNVERIFIED", result.stdout)
        self.assertIn("taxon", result.stdout)
        self.assertIn("[ERROR]", result.stdout)

    def test_registry_amino_acid_mismatch_stays_unverified(self):
        # the record exists for TAA->Trp, but the declaration claims Gln
        result = self._run("aa-mismatch.gb", INTERNAL_STOP_CDS, "(pos:10..12,aa:Gln)",
                           registry=[_transl_record()])
        self.assertIn("TRANSL_EXCEPT_DECLARED_UNVERIFIED", result.stdout)
        self.assertIn("[ERROR]", result.stdout)


class TranslExceptRegistryFormatTests(unittest.TestCase):
    """P1-2 / P2-1: a malformed transl_except registry fails closed at LOAD time.

    A record that does not say which stop it covers would silently become a
    gene-wide rule, and evidence fields with the wrong JSON type would be
    stringified into a non-empty "audited" value.
    """

    def setUp(self):
        self.module = load_module("annot_transl_registry", Path("scripts") / "annot_check.py")
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name)

    def _load(self, *entries):
        path = self.dir / "registry.json"
        path.write_text(json.dumps({"exceptions": list(entries)}, ensure_ascii=False),
                        encoding="utf-8")
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                self.module.load_exception_registry(str(path))
        return caught.exception.code

    @staticmethod
    def _record(**overrides):
        record = {"gene": "cox1", "codon": "TAA", "amino_acid": "Trp", "transl_table": 5,
                  "taxon": "Lepidoptera", "source": "DOI 10.1000/example",
                  "rationale": "documented", "pos": "10..12"}
        record.update(overrides)
        return record

    def test_record_without_a_site_binding_is_rejected(self):
        record = self._record()
        record.pop("pos")
        self.assertEqual(self._load(record), 1)

    def test_codon_index_and_pos_together_are_rejected(self):
        self.assertEqual(self._load(self._record(codon_index=4)), 1)

    def test_unknown_scope_is_rejected(self):
        record = self._record()
        record.pop("pos")
        record["scope"] = "gene"
        self.assertEqual(self._load(record), 1)

    def test_scope_with_a_site_binding_is_rejected(self):
        self.assertEqual(self._load(self._record(scope="gene_wide")), 1)

    def test_duplicate_site_records_are_rejected(self):
        self.assertEqual(self._load(self._record(), self._record()), 1)

    def _load_ok(self, *entries):
        path = self.dir / "ok.json"
        path.write_text(json.dumps({"exceptions": list(entries)}, ensure_ascii=False),
                        encoding="utf-8")
        return self.module.load_exception_registry(str(path))

    def test_distinct_sites_are_kept(self):
        registry = self._load_ok(self._record(), self._record(pos="19..21"))
        self.assertEqual(len(registry["transl_except"][("cox1", "TAA", "trp")]), 2)

    def test_same_site_for_different_taxa_is_kept(self):
        # P2-1: the runtime SELECTS on taxon, so the duplicate identity must not
        # collapse records that differ only by taxon
        registry = self._load_ok(self._record(taxon="Lepidoptera"),
                                 self._record(taxon="Diptera"))
        self.assertEqual(len(registry["transl_except"][("cox1", "TAA", "trp")]), 2)

    def test_same_site_for_different_tables_is_kept(self):
        registry = self._load_ok(self._record(transl_table=5), self._record(transl_table=2))
        self.assertEqual(len(registry["transl_except"][("cox1", "TAA", "trp")]), 2)

    def test_same_site_same_taxon_same_table_is_rejected(self):
        self.assertEqual(self._load(self._record(), self._record()), 1)

    def test_empty_or_punctuation_gene_selector_is_rejected(self):
        # The bug: an empty selector matches an empty selector, so a CDS with no
        # /gene or /product would have gene-keyed evidence applied to it.
        for gene in ("", "?", "---", "   ", "(CUN)"):
            with self.subTest(gene=gene):
                self.assertEqual(self._load(self._record(gene=gene)), 1)

    def test_empty_amino_acid_selector_is_rejected(self):
        self.assertEqual(self._load(self._record(amino_acid="")), 1)

    def test_explicit_null_amino_acid_is_a_malformed_transl_except_record(self):
        # P2-1: the KEY PRESENCE is what classifies a record; null must not silently
        # downgrade a malformed transl_except record into a start-codon record
        self.assertEqual(self._load(self._record(amino_acid=None)), 1)

    def test_explicit_null_gene_or_codon_is_rejected(self):
        for field in ("gene", "codon"):
            with self.subTest(field=field):
                self.assertEqual(self._load(self._record(**{field: None})), 1)

    def test_selector_classification_three_way(self):
        # absent key -> start record; "" -> malformed transl_except; null -> malformed
        start = {"gene": "cox1", "codon": "CGA", "taxon": "Lepidoptera",
                 "source": "DOI", "rationale": "r"}
        registry = self._load_ok(start)
        self.assertIn(("cox1", "CGA"), registry["start"])
        self.assertEqual(registry["transl_except"], {})
        self.assertEqual(self._load(self._record(amino_acid="")), 1)
        self.assertEqual(self._load(self._record(amino_acid=None)), 1)

    def test_invalid_codon_selector_is_rejected(self):
        for codon in ("", "TA", "TAXA", "TXX"):
            with self.subTest(codon=codon):
                self.assertEqual(self._load(self._record(codon=codon)), 1)

    def test_invalid_amino_acid_token_is_rejected(self):
        self.assertEqual(self._load(self._record(amino_acid="Zzz")), 1)

    def test_gene_wide_scope_record_is_kept(self):
        record = self._record()
        record.pop("pos")
        record["scope"] = "gene_wide"
        registry = self._load_ok(record)
        self.assertEqual(registry["transl_except"][("cox1", "TAA", "trp")][0]["site"],
                         ("gene_wide", None))

    def test_wrong_evidence_field_types_are_rejected(self):
        cases = (("taxon", ["Coleoptera"]), ("taxon", 3), ("source", {"doi": "x"}),
                 ("rationale", 1), ("transl_table", "5"), ("codon_index", "4"),
                 ("gene", 1), ("amino_acid", 1))
        for field, value in cases:
            with self.subTest(field=field, value=value):
                self.assertEqual(self._load(self._record(**{field: value})), 1)

    def test_fuzzy_position_is_rejected(self):
        self.assertEqual(self._load(self._record(pos="<10..12")), 1)

    def test_mixed_strand_position_is_rejected(self):
        self.assertEqual(self._load(self._record(pos="join(complement(10..11),12)")), 1)

    def test_start_records_are_unaffected(self):
        path = self.dir / "start.json"
        path.write_text(json.dumps({"exceptions": [
            {"gene": "cox1", "codon": "CGA", "taxon": "Lepidoptera",
             "source": "DOI", "rationale": "r"}]}), encoding="utf-8")
        registry = self.module.load_exception_registry(str(path))
        self.assertIn(("cox1", "CGA"), registry["start"])
        self.assertEqual(registry["transl_except"], {})


class StartRegistrySelectorTests(unittest.TestCase):
    """Follow-up P2: start records/CLI must not accept an empty canonical gene selector.

    An empty canonical gene matches the empty canonical gene of a CDS that has
    neither /gene nor /product, so an unnamed CDS could be reported as "accepted
    per audited record" -- evidence bound to no gene identity at all.
    """

    def setUp(self):
        self.module = load_module("annot_start_selector", Path("scripts") / "annot_check.py")
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name)
        # the argument validation runs after the GenBank file is read, so the CLI
        # tests need a readable record rather than a nonexistent path
        self.genome = self.dir / "genome.gb"
        code = "CGA" + "AAA" * 20 + "TAA"
        write_gb_raw(self.genome, code + "A" * 40,
                     [{"location": "1..%d" % len(code), "type": "CDS", "gene": "cox1"}])

    def _load(self, *entries):
        path = self.dir / "start.json"
        path.write_text(json.dumps({"exceptions": list(entries)}, ensure_ascii=False),
                        encoding="utf-8")
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                self.module.load_exception_registry(str(path))
        return caught.exception.code

    def _load_ok(self, *entries):
        path = self.dir / "ok.json"
        path.write_text(json.dumps({"exceptions": list(entries)}, ensure_ascii=False),
                        encoding="utf-8")
        return self.module.load_exception_registry(str(path))

    @staticmethod
    def _record(**overrides):
        record = {"gene": "cox1", "codon": "CGA", "taxon": "Lepidoptera",
                  "source": "DOI 10.1000/example", "rationale": "documented"}
        record.update(overrides)
        return record

    def test_empty_canonical_gene_is_rejected(self):
        for gene in ("", "?", "---", "   ", "(CUN)"):
            with self.subTest(gene=gene):
                self.assertEqual(self._load(self._record(gene=gene)), 1)

    def test_empty_codon_is_rejected(self):
        self.assertEqual(self._load(self._record(codon="")), 1)

    def test_valid_start_record_still_loads(self):
        registry = self._load_ok(self._record())
        self.assertIn(("cox1", "CGA"), registry["start"])

    def test_codon_must_be_three_iupac_bases(self):
        for codon in ("", "   ", "CG", "XXXX", "C1A"):
            with self.subTest(codon=codon):
                self.assertEqual(self._load(self._record(codon=codon)), 1)

    def test_codon_is_stripped_and_upper_cased(self):
        for codon in ("CGA", "cga", " CGA ", "\tcga\n"):
            with self.subTest(codon=codon):
                registry = self._load_ok(self._record(codon=codon))
                self.assertIn(("cox1", "CGA"), registry["start"])

    def test_duplicate_gene_codon_is_rejected(self):
        # Two records for the same (gene, codon) cannot both be represented: the
        # second silently OVERWROTE the first, so this fails closed now.  The old
        # test only asserted the key existed and therefore passed on the overwrite.
        self.assertEqual(self._load(self._record(taxon="Lepidoptera"),
                                    self._record(taxon="Diptera")), 1)
        self.assertEqual(self._load(self._record(), self._record()), 1)

    def test_distinct_gene_or_codon_still_coexist(self):
        registry = self._load_ok(self._record(codon="CGA"), self._record(codon="CGG"))
        self.assertEqual(sorted(registry["start"]), [("cox1", "CGA"), ("cox1", "CGG")])
        registry = self._load_ok(self._record(codon="CGA"),
                                 self._record(gene="cox2", codon="CGA"))
        self.assertEqual(sorted(registry["start"]), [("cox1", "CGA"), ("cox2", "CGA")])

    @unittest.skipUnless(biopython_available(), "Biopython is not installed")
    def test_cli_rejects_a_malformed_codon(self):
        for token in ("cox1:CG", "cox1:XXXX", "cox1:   "):
            with self.subTest(token=token):
                result = run_annot_check(self.genome, "--tolerate-start", token)
                self.assertEqual(result.returncode, 1)
                self.assertIn("三个 IUPAC 碱基", result.stdout)

    @unittest.skipUnless(biopython_available(), "Biopython is not installed")
    def test_cli_accepts_a_stripped_upper_cased_codon(self):
        result = run_annot_check(self.genome, "--allow-atypical", REASON,
                                 "--tolerate-start", "cox1:cga")
        self.assertNotIn("三个 IUPAC 碱基", result.stdout)
        self.assertIn("已按 --tolerate-start 接受", result.stdout)

    @unittest.skipUnless(biopython_available(), "Biopython is not installed")
    def test_cli_rejects_an_empty_gene_selector(self):
        result = run_annot_check(self.genome, "--tolerate-start", ":CGA")
        self.assertEqual(result.returncode, 1)
        self.assertIn("基因名归一化后为空", result.stdout)

    @unittest.skipUnless(biopython_available(), "Biopython is not installed")
    def test_cli_rejects_a_punctuation_only_gene_selector(self):
        result = run_annot_check(self.genome, "--tolerate-start", "?:CGA")
        self.assertEqual(result.returncode, 1)
        self.assertIn("基因名归一化后为空", result.stdout)

    @unittest.skipUnless(biopython_available(), "Biopython is not installed")
    def test_cli_still_rejects_a_missing_codon(self):
        result = run_annot_check(self.genome, "--tolerate-start", "cox1")
        self.assertEqual(result.returncode, 1)
        self.assertIn("基因:密码子", result.stdout)

    @unittest.skipUnless(biopython_available(), "Biopython is not installed")
    def test_unnamed_cds_is_not_reported_as_accepted(self):
        code = "CGA" + "AAA" * 20 + "TAA"
        genome = self.dir / "unnamed.gb"
        write_gb_raw(genome, code + "A" * 40,
                     [{"location": "1..%d" % len(code), "type": "CDS"}])
        registry = self.dir / "start.json"
        registry.write_text(json.dumps({"exceptions": [{k: v for k, v in self._record().items()}]}),
                            encoding="utf-8")
        # the registry cannot even load with an empty canonical gene, so the CLI fails
        # closed instead of describing an unidentified CDS as an audited exception
        result = run_annot_check(genome, "--allow-atypical", REASON,
                                 "--tolerate-start", ":CGA",
                                 "--exception-registry", str(registry))
        self.assertEqual(result.returncode, 1)
        self.assertIn("基因名归一化后为空", result.stdout)
        self.assertNotIn("按已审计例外记录接受", result.stdout)

    @unittest.skipUnless(biopython_available(), "Biopython is not installed")
    def test_named_cds_with_a_named_record_is_still_accepted_as_review(self):
        code = "CGA" + "AAA" * 20 + "TAA"
        genome = self.dir / "named.gb"
        write_gb_raw(genome, code + "A" * 40,
                     [{"location": "1..%d" % len(code), "type": "CDS", "gene": "cox1"}])
        registry = self.dir / "named.json"
        registry.write_text(json.dumps({"exceptions": [self._record()]}), encoding="utf-8")
        result = run_annot_check(genome, "--allow-atypical", REASON,
                                 "--tolerate-start", "cox1:CGA",
                                 "--exception-registry", str(registry), "--taxon", "Lepidoptera")
        self.assertIn("按已审计例外记录接受", result.stdout)
        self.assertNotIn("EXCEPTION_NOT_REGISTERED", result.stdout)


class NullSelectorMessageTests(unittest.TestCase):
    """Follow-up P3: the null selector message must not misdescribe the record type.

    Only the amino_acid KEY decides that a record is a transl_except record;
    gene/codon belong to start records too, so their null message must not claim
    to be about transl_except classification.
    """

    def setUp(self):
        self.module = load_module("annot_null_selector", Path("scripts") / "annot_check.py")
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name)

    def _message(self, **overrides):
        record = {"gene": "cox1", "codon": "CGA", "taxon": "Lepidoptera",
                  "source": "DOI 10.1000/example", "rationale": "documented"}
        record.update(overrides)
        path = self.dir / "registry.json"
        path.write_text(json.dumps({"exceptions": [record]}), encoding="utf-8")
        captured = io.StringIO()
        try:
            with contextlib.redirect_stdout(captured):
                self.module.load_exception_registry(str(path))
        except SystemExit:
            return captured.getvalue()
        self.fail("expected the registry load to fail")

    def test_gene_null_message_does_not_claim_transl_except_classification(self):
        message = self._message(gene=None)
        self.assertIn("gene 键存在但为 null", message)
        self.assertNotIn("transl_except 记录", message)

    def test_codon_null_message_does_not_claim_transl_except_classification(self):
        message = self._message(codon=None)
        self.assertIn("codon 键存在但为 null", message)
        self.assertNotIn("transl_except 记录", message)

    def test_amino_acid_null_message_does_mention_the_record_type(self):
        message = self._message(amino_acid=None)
        self.assertIn("amino_acid 键存在但为 null", message)
        self.assertIn("transl_except 记录", message)


class SchemaTypeMatrixTests(unittest.TestCase):
    """P1-3: every legal JSON value must yield a verdict, never a traceback.

    The previous fixed data set only used the string "x" for array-type errors.
    Strings are iterable, which hid `enumerate(1)` raising TypeError.  This is
    the Schema-keyword x JSON-basic-type matrix instead of hand-picked cases.
    """

    VALUES = (1, 0, -1, "x", "", [], [1], {}, {"a": 1}, None, True, 1.5)
    FIELDS = ("schema_version", "case_id", "case_type", "taxon", "inputs", "issue",
              "hypotheses", "events_file", "decision", "anomalies", "modifications",
              "validation", "lessons_proposed")

    def setUp(self):
        self.module = load_module("experience_type_matrix", Path("tools") / "experience.py")
        self.schema = json.loads((ROOT / "schemas" / "case.schema.json").read_text(encoding="utf-8"))

    @staticmethod
    def _base():
        return {
            "schema_version": "2.0", "case_id": "case-1",
            "inputs": [{"role": "assembly_fasta", "path": "/tmp/x.fa", "sha256": "a" * 64}],
            "issue": {"type": "internal_stop"},
            "hypotheses": [{"id": "H1", "explanation": "e", "support": [], "against": [],
                            "unknown": []}],
            "decision": {"status": "RESOLVED", "confidence": "moderate", "rationale": "r"},
            "anomalies": [], "events_file": "events.jsonl",
        }

    def test_no_json_value_crashes_the_fallback(self):
        for field in self.FIELDS:
            for value in self.VALUES:
                with self.subTest(field=field, value=value):
                    case = self._base()
                    case[field] = value
                    try:
                        errors = self.module._case_errors_without_jsonschema(case)
                    except Exception as exc:  # noqa: BLE001 - this IS the assertion
                        self.fail("%s=%r raised %s: %s" % (field, value, type(exc).__name__, exc))
                    self.assertIsInstance(errors, list)

    def test_no_json_value_crashes_the_public_validator(self):
        for field in self.FIELDS:
            for value in self.VALUES:
                with self.subTest(field=field, value=value):
                    case = self._base()
                    case[field] = value
                    try:
                        self.module.case_schema_errors(case)
                    except Exception as exc:  # noqa: BLE001
                        self.fail("case_schema_errors %s=%r raised %s: %s"
                                  % (field, value, type(exc).__name__, exc))

    def test_matrix_agrees_with_jsonschema(self):
        try:
            import jsonschema
        except ImportError:
            self.skipTest("jsonschema is not installed")
        disagreements = []
        for field in self.FIELDS:
            for value in self.VALUES:
                case = self._base()
                case[field] = value
                try:
                    jsonschema.validate(case, self.schema)
                    schema_ok = True
                except jsonschema.ValidationError:
                    schema_ok = False
                fallback_ok = not self.module._case_errors_without_jsonschema(case)
                if schema_ok != fallback_ok:
                    disagreements.append('%s=%r jsonschema=%s fallback=%s'
                                         % (field, value, schema_ok, fallback_ok))
        self.assertEqual(disagreements, [])

    def test_a_legal_case_passes_every_path(self):
        case = self._base()
        self.assertEqual(self.module._case_errors_without_jsonschema(case), [])
        self.assertEqual(self.module.case_schema_errors(case), [])


class CaseValidateEntryPointTests(unittest.TestCase):
    """P1-4: the public CLI must report INVALID, never leak a traceback."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name) / "case"
        self.dir.mkdir()
        (self.dir / "events.jsonl").touch()

    def _payload(self, decision_value="__absent__"):
        case = SchemaTypeMatrixTests._base()
        if decision_value == "__absent__":
            case.pop("decision", None)
        else:
            case["decision"] = decision_value
        return case

    def _run_cli(self, payload):
        (self.dir / "case.json").write_text(json.dumps(payload), encoding="utf-8")
        return run_script(Path("tools") / "experience.py", "case-validate", self.dir)

    def test_explicit_null_decision_is_invalid_without_traceback(self):
        result = self._run_cli(self._payload(None))
        self.assertEqual(result.returncode, 1)
        self.assertIn("INVALID", result.stdout)
        self.assertNotIn("Traceback", result.stderr)
        self.assertNotIn("AttributeError", result.stderr)

    def test_wrong_decision_type_is_invalid_without_traceback(self):
        for value in ("RESOLVED", 1, []):
            with self.subTest(value=value):
                result = self._run_cli(self._payload(value))
                self.assertEqual(result.returncode, 1)
                self.assertIn("INVALID", result.stdout)
                self.assertNotIn("Traceback", result.stderr)

    def test_top_level_array_case_is_invalid_without_traceback(self):
        result = self._run_cli([])
        self.assertEqual(result.returncode, 1)
        self.assertIn("INVALID", result.stdout)
        self.assertNotIn("Traceback", result.stderr)

    def test_valid_case_reports_valid(self):
        for absent_or_null in (None, "__absent__"):
            with self.subTest(decision=absent_or_null):
                result = self._run_cli(self._payload(absent_or_null))
                self.assertEqual(result.returncode, 1)      # null / missing is invalid
        result = self._run_cli(self._payload("KEEP"))
        self.assertEqual(result.returncode, 1)
        result = self._run_cli(SchemaTypeMatrixTests._base())
        self.assertEqual(result.returncode, 0)
        self.assertIn("VALID", result.stdout)


class FixtureHygieneTests(unittest.TestCase):
    """P3-1: `git diff --check` must stay clean (new blank line at EOF)."""

    def test_xml_fixtures_end_with_exactly_one_newline(self):
        for path in sorted(FIXTURES.glob("*.xml")):
            with self.subTest(fixture=path.name):
                text = path.read_text(encoding="utf-8")
                self.assertTrue(text.endswith("\n"), path.name)
                self.assertFalse(text.endswith("\n\n"), "%s ends with a blank line" % path.name)
                for number, line in enumerate(text.split("\n"), start=1):
                    self.assertEqual(line, line.rstrip(),
                                     "%s:%d has trailing whitespace" % (path.name, number))


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Failing reproductions for the annotation-check defects (P1-P7) plus the two
expert-added defects (bare trnL/trnS normalisation, whole-genome orientation).

These are written against the *intended* interface; they fail on the pre-fix
script and pass afterwards.
"""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gb_fixtures import (  # noqa: E402
    ROOT, ORDER, biopython_available, build_synthetic_ref, cds_block, load_module,
    run_annot_check, write_gb,
)

REASON = "合成测试: 单基因记录, 基因集差异已逐项确认"


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class NoncanonicalStartTests(unittest.TestCase):
    """P1: a non-canonical start must be reviewable, not auto-rejected."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def _cga_cox1(self):
        sequence = "CGA" + "AAA" * 20 + "TAA" + "A" * 50
        path = self.dir / "cga.gb"
        write_gb(path, sequence, [
            {"gene": "cox1", "type": "CDS", "start": 0, "end": 66, "strand": "+"},
        ])
        return path

    def test_noncanonical_start_is_error_without_declared_exception(self):
        result = run_annot_check(self._cga_cox1(), "--allow-atypical", REASON)
        self.assertIn("[ERROR]", result.stdout)
        self.assertIn("NONCANONICAL_START_REVIEW", result.stdout)
        self.assertIn("CGA", result.stdout)
        self.assertEqual(result.returncode, 1)

    def test_noncanonical_start_reviews_with_gene_codon_exception(self):
        result = run_annot_check(self._cga_cox1(), "--allow-atypical", REASON,
                                 "--tolerate-start", "cox1:CGA")
        self.assertNotIn("[ERROR]", result.stdout)
        self.assertIn("NONCANONICAL_START_REVIEW", result.stdout)
        self.assertNotEqual(result.returncode, 1)

    def test_exception_is_gene_and_codon_specific(self):
        # cox1:TTG must not excuse cox1:CGA, and cox2:CGA must not excuse cox1:CGA
        for bad in ("cox1:TTG", "cox2:CGA"):
            with self.subTest(exception=bad):
                result = run_annot_check(self._cga_cox1(), "--allow-atypical", REASON,
                                         "--tolerate-start", bad)
                self.assertIn("[ERROR]", result.stdout)
                self.assertIn("未匹配任何起始密码子", result.stdout)

    def test_real_partial_cds_codon_start_is_not_a_start_error(self):
        # 5' partial: first complete codon sits at offset 2, so it is not a start codon
        sequence = "GG" + "ATG" + "AAA" * 20 + "TAA"
        path = self.dir / "partial.gb"
        write_gb(path, sequence, [
            {"gene": "cox1", "type": "CDS", "start": 0, "end": len(sequence),
             "strand": "+", "codon_start": 2},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertNotIn("[ERROR]", result.stdout)
        self.assertIn("5' 端不完整", result.stdout)
        self.assertNotEqual(result.returncode, 1)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class IncompleteStopTests(unittest.TestCase):
    """P6: a note must not excuse a non-T/TA terminal codon."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def test_sense_codon_terminal_is_error_even_with_note(self):
        sequence = "ATG" + "AAA" * 19 + "GAT" + "A" * 60
        path = self.dir / "sense-terminal.gb"
        write_gb(path, sequence, [
            {"gene": "cox1", "type": "CDS", "start": 0, "end": 63, "strand": "+",
             "note": "ok"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("[ERROR]", result.stdout)
        self.assertIn("终止密码子", result.stdout)
        self.assertEqual(result.returncode, 1)

    def test_genuine_t_prefix_terminal_is_review_only(self):
        sequence = "ATG" + "AAA" * 20 + "TA" + "A" * 60  # length 65 -> 2 bases past the last codon
        path = self.dir / "ta-terminal.gb"
        write_gb(path, sequence, [
            {"gene": "cox1", "type": "CDS", "start": 0, "end": 65, "strand": "+"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertNotIn("[ERROR]", result.stdout)
        self.assertIn("T/TA", result.stdout)
        self.assertNotEqual(result.returncode, 1)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class CrossOriginAndOverlapTests(unittest.TestCase):
    """P3/P7 and the new default overlap policy."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def test_cross_origin_join_uses_per_segment_length_and_overlap(self):
        length = 2000
        bases = ["A"] * length
        bases[1949:1952] = list("ATG")
        bases[48:51] = list("TAA")
        sequence = "".join(bases)
        path = self.dir / "join.gb"
        write_gb(path, sequence, [
            {"gene": "cox1", "type": "CDS", "strand": "+",
             "segments": [(1949, 2000), (0, 51)]},
            {"gene": "trnL2", "type": "tRNA", "start": 60, "end": 130, "strand": "+"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("102bp", result.stdout)
        # 60..130 does not overlap either segment of the join
        self.assertNotIn("cox1", self._overlap_lines(result.stdout))

    @staticmethod
    def _overlap_lines(stdout):
        return "\n".join(line for line in stdout.splitlines() if "重叠" in line)

    def test_long_cds_trna_overlap_defaults_to_review_not_error(self):
        path = self.dir / "long-overlap.gb"
        write_gb(path, "A" * 500, [
            {"gene": "atp8", "type": "CDS", "start": 0, "end": 100, "strand": "+"},
            {"gene": "trnQ", "type": "tRNA", "start": 90, "end": 170, "strand": "+"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertNotIn("[ERROR]", result.stdout)
        self.assertIn("[WARN ]", result.stdout)
        self.assertIn("10bp", result.stdout)
        self.assertNotEqual(result.returncode, 1)

    def test_short_overlap_is_recorded_at_info_level(self):
        path = self.dir / "short-overlap.gb"
        write_gb(path, "A" * 500, [
            {"gene": "atp8", "type": "CDS", "start": 0, "end": 100, "strand": "+"},
            {"gene": "atp6", "type": "CDS", "start": 95, "end": 195, "strand": "+"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON)
        self.assertIn("[INFO ]", result.stdout)
        self.assertIn("5bp", result.stdout)

    def test_tolerate_overlap_is_case_insensitive_and_reports_unused(self):
        path = self.dir / "tolerate.gb"
        write_gb(path, "A" * 500, [
            {"gene": "nad5", "type": "CDS", "start": 0, "end": 100, "strand": "-"},
            {"gene": "trnH", "type": "tRNA", "start": 90, "end": 170, "strand": "-"},
        ])
        result = run_annot_check(path, "--allow-atypical", REASON,
                                 "--tolerate-overlap", "ND5,trnH",
                                 "--tolerate-overlap", "atp6,trnA")
        self.assertIn("ACCEPT", result.stdout)
        self.assertIn("未匹配任何重叠对", result.stdout)
        self.assertNotIn("[ERROR]", result.stdout)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class TrnaIdentityTests(unittest.TestCase):
    """New defect: an unclassified trnL/trnS must stay undetermined."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)
        self.module = load_module("annot_check_identity", Path("scripts") / "annot_check.py")

    def _bare(self):
        path = self.dir / "bare-trna.gb"
        write_gb(path, "A" * 400, [
            {"gene": "trnL", "type": "tRNA", "start": 0, "end": 70, "strand": "+"},
            {"gene": "trnS", "type": "tRNA", "start": 100, "end": 170, "strand": "+"},
        ])
        return path

    def test_bare_trnl_trns_are_not_autonormalized(self):
        result = run_annot_check(self._bare())
        self.assertIn("UNDETERMINED_TRNA", result.stdout)
        self.assertNotIn("缺少基因身份: trnl1", result.stdout)
        self.assertNotIn("缺少基因身份: trnl2", result.stdout)
        self.assertNotIn("缺少基因身份: trns1", result.stdout)
        self.assertNotIn("缺少基因身份: trns2", result.stdout)
        self.assertNotIn("非标准/未知基因身份: trnl", result.stdout)

    def test_canonical_gene_keeps_bare_trnl_trns_undetermined(self):
        self.assertEqual(self.module.canonical_gene(SimpleNamespace(qualifiers={"gene": ["trnL"]})), "trnl")
        self.assertEqual(self.module.canonical_gene(SimpleNamespace(qualifiers={"gene": ["trnS"]})), "trns")

    def test_anticodon_resolves_trna_type(self):
        cases = {"tag": "trnl1", "taa": "trnl2", "gct": "trns1", "tga": "trns2"}
        for anticodon, expected in cases.items():
            with self.subTest(anticodon=anticodon):
                feature = SimpleNamespace(qualifiers={
                    "gene": ["trnL"],
                    "anticodon": ["(pos:1..70,aa:Leu,seq:%s)" % anticodon],
                })
                self.assertEqual(self.module.resolve_trna_identity(feature), expected)

    def test_parenthesised_gene_name_is_accepted(self):
        feature = SimpleNamespace(qualifiers={"gene": ["trnL2(UUR)"]})
        self.assertEqual(self.module.canonical_gene(feature), "trnl2")
        feature = SimpleNamespace(qualifiers={"gene": ["trnS1(AGN)"]})
        self.assertEqual(self.module.canonical_gene(feature), "trns1")


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
class OrientationAndAliasTests(unittest.TestCase):
    """New defect: whole-genome revcomp; P2/P5 alias normalisation."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def test_reverse_complemented_genome_does_not_warn_on_strand_distribution(self):
        from Bio.Seq import Seq

        ref, sequence = build_synthetic_ref(self.dir)
        reverse = str(Seq(sequence).reverse_complement())
        total = len(sequence)
        features = []
        for name, kind, length, strand in ORDER:
            start = sum(item[2] for item in ORDER[:len(features)])
            features.append({
                "gene": name, "type": kind, "strand": "-" if strand == "+" else "+",
                "start": total - (start + length), "end": total - start,
            })
        path = self.dir / "revcomp.gb"
        write_gb(path, reverse, features)
        result = run_annot_check(path, "--table", "5")
        self.assertNotIn("正链 CDS", result.stdout)
        self.assertIn("反向互补", result.stdout)

    def test_cross_origin_feature_helpers(self):
        module = load_module("annot_check_segments", Path("scripts") / "annot_check.py")
        feature = SimpleNamespace(location=SimpleNamespace(
            parts=[SimpleNamespace(start=1950, end=2000), SimpleNamespace(start=0, end=50)]))
        self.assertEqual(module.feature_segments(feature), [(1950, 2000), (0, 50)])
        self.assertEqual(module.feature_length(feature), 100)
        self.assertEqual(module.segment_overlap([(1950, 2000), (0, 50)], [(60, 130)]), 0)
        self.assertEqual(module.segment_overlap([(1950, 2000), (0, 50)], [(30, 70)]), 20)

    def test_reference_gene_aliases_do_not_report_false_order_difference(self):
        ref, _ = build_synthetic_ref(self.dir)
        from Bio import SeqIO
        alias_map = {"nad1": "nd1", "nad2": "nd2", "nad3": "nd3", "nad4": "nd4",
                     "nad4l": "nd4l", "nad5": "nd5", "nad6": "nd6", "cob": "cytb"}
        record = SeqIO.read(ref, "genbank")
        for feature in record.features:
            if feature.type != "CDS":
                continue
            name = feature.qualifiers.get("gene", [""])[0].lower()
            if name in alias_map:
                feature.qualifiers["gene"] = [alias_map[name]]
        aliased = self.dir / "ref-aliased.gb"
        SeqIO.write(record, aliased, "genbank")

        result = run_annot_check(ref, "--ref", aliased)
        self.assertIn("基因顺序匹配(CDS): 一致", result.stdout)
        self.assertNotIn("基因顺序与参考不同", result.stdout)

    def test_rrna_aliases_use_the_correct_length_range(self):
        for name, length, forbidden in (("lrrna", 1300, "600-850"), ("srrna", 750, "1100-1500")):
            with self.subTest(gene=name):
                path = self.dir / ("%s.gb" % name)
                write_gb(path, "A" * 1500, [
                    {"gene": name, "type": "rRNA", "start": 0, "end": length, "strand": "-"},
                ])
                result = run_annot_check(path, "--allow-atypical", REASON)
                self.assertNotIn("超出范围 %s" % forbidden, result.stdout)


class ArgumentContractTests(unittest.TestCase):
    def test_invalid_overlap_severity_exits_1(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "annot_check.py"),
             "/nonexistent.gb", "--overlap-severity", "bogus"],
            capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 1)

    def test_allow_atypical_requires_a_reason(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "annot_check.py"),
             "/nonexistent.gb", "--allow-atypical"],
            capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 1)
        self.assertIn("--allow-atypical", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
import importlib.util
import shutil
import subprocess
import tempfile
import time
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]


def load_experience():
    path = ROOT / "tools" / "experience.py"
    spec = importlib.util.spec_from_file_location("experience", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_script(name, relative_path):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ExperienceSecurityTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.module = load_experience()
        self.knowledge = Path(self.temp_dir.name) / "knowledge"
        self.cases = self.knowledge / "cases"
        self.cases.mkdir(parents=True)
        self.module.KNOWLEDGE = str(self.knowledge)
        self.module.CASES = str(self.cases)
        self.module.SIGNALS = str(self.knowledge / "signals.md")
        self.module.PITFALLS = str(self.knowledge / "pitfalls.md")
        self.module.STATS = str(self.knowledge / "stats.md")

    def args(self, sample, **overrides):
        values = {
            "sample": sample,
            "species": "test",
            "tools": "test",
            "signals": "signal",
            "diagnosis": "diagnosis",
            "actions": "actions",
            "lessons": "lessons",
            "force": False,
        }
        values.update(overrides)
        return SimpleNamespace(**values)

    def test_rejects_sample_path_traversal(self):
        with self.assertRaises(ValueError):
            self.module.add_case(self.args("../escape"))
        self.assertFalse((Path(self.temp_dir.name) / "escape.md").exists())

    def test_does_not_overwrite_existing_case(self):
        path = self.cases / "sample.md"
        path.write_text("original", encoding="utf-8")
        with self.assertRaises(FileExistsError):
            self.module.add_case(self.args("sample"))
        self.assertEqual(path.read_text(encoding="utf-8"), "original")

    def test_force_replaces_existing_case(self):
        path = self.cases / "sample.md"
        path.write_text("original", encoding="utf-8")
        self.module.add_case(self.args("sample", force=True))
        self.assertIn("sample: sample", path.read_text(encoding="utf-8"))

    def test_force_does_not_follow_case_symlink(self):
        outside = Path(self.temp_dir.name) / "outside.md"
        outside.write_text("outside", encoding="utf-8")
        link = self.cases / "sample.md"
        link.symlink_to(outside)
        with self.assertRaises(ValueError):
            self.module.add_case(self.args("sample", force=True))
        self.assertEqual(outside.read_text(encoding="utf-8"), "outside")

    def test_missing_case_directory_is_safe_for_reads(self):
        self.module.CASES = str(Path(self.temp_dir.name) / "missing")
        self.assertEqual(self.module.case_files(), [])
        self.module.stats()
        self.module.suggest_promotions()

    def test_structured_case_lifecycle_is_integrated_with_experience(self):
        assembly = Path(self.temp_dir.name) / "assembly.fasta"
        assembly.write_text(">ctg1\nACGTN\n", encoding="utf-8")
        case = Path(self.temp_dir.name) / "case"
        self.module.case_init(SimpleNamespace(directory=str(case), case_id=None, issue="internal_stop", observation="nad5 stop", taxon=None, input=[["assembly_fasta", str(assembly)]]))
        self.module.case_event(SimpleNamespace(directory=str(case), action="annot_check", result="stop remains", impact="H1:against", command="", tool_version="", motivation=""))
        self.assertEqual(self.module.case_validate(SimpleNamespace(directory=str(case))), 0)
        self.module.case_report(SimpleNamespace(directory=str(case)))
        self.assertIn("annot_check", (case / "case.md").read_text(encoding="utf-8"))


class ExternalUploadTests(unittest.TestCase):
    def test_cox1_requires_explicit_public_upload_opt_in(self):
        missing = Path(tempfile.gettempdir()) / "mito-cox1-missing.fasta"
        result = subprocess.run(
            ["python3", str(ROOT / "scripts" / "cox1_id.py"), str(missing), "--coords", "1,4"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--allow-public-upload", result.stdout + result.stderr)

    def test_cox1_rejects_low_identity_and_tied_hits(self):
        module = load_script("cox1_id", Path("scripts") / "cox1_id.py")
        self.assertIsNone(module.select_supported_hit([(80.0, "low")], 85.0))
        self.assertIsNone(module.select_supported_hit([(95.0, "a"), (95.0, "b")], 85.0))
        self.assertIsNone(module.select_supported_hit([(95.0, "a"), (94.5, "b")], 85.0))
        self.assertEqual(module.select_supported_hit([(95.0, "a"), (90.0, "b")], 85.0), (95.0, "a"))


class BackgroundScriptSecurityTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.root = Path(self.temp_dir.name)
        (self.root / "scripts").mkdir()
        for name in ("run_bg.sh", "check_bg.sh"):
            shutil.copy2(ROOT / "scripts" / name, self.root / "scripts" / name)

    def run_script(self, name, *args):
        return subprocess.run(
            ["bash", str(self.root / "scripts" / name), *args],
            capture_output=True,
            text=True,
            timeout=10,
        )

    def test_rejects_task_name_path_traversal(self):
        marker = self.root / "marker"
        result = self.run_script("run_bg.sh", "../escape", "--", "touch", str(marker))
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(marker.exists())
        self.assertFalse((self.root / "escape.status").exists())

    def test_default_mode_does_not_parse_command_string(self):
        marker = self.root / "marker"
        result = self.run_script("run_bg.sh", "safe", "--", f"touch {marker}")
        time.sleep(0.3)
        self.assertFalse(marker.exists())

    def test_explicit_trusted_shell_mode_executes_shell_command(self):
        marker = self.root / "marker"
        result = self.run_script("run_bg.sh", "trusted", "--trusted-shell", f"touch {marker}")
        self.assertEqual(result.returncode, 0)
        for _ in range(20):
            if marker.exists():
                break
            time.sleep(0.1)
        self.assertTrue(marker.exists())

    def test_check_rejects_task_name_path_traversal(self):
        status = self.root / "escape.status"
        status.write_text("done\n", encoding="utf-8")
        result = self.run_script("check_bg.sh", "../escape")
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("已完成", result.stdout)

    def test_launcher_does_not_hold_output_pipes(self):
        try:
            result = subprocess.run(
                ["bash", str(self.root / "scripts" / "run_bg.sh"), "detach", "--", "sleep", "1"],
                capture_output=True,
                text=True,
                timeout=2,
            )
        except subprocess.TimeoutExpired:
            self.fail("background supervisor kept stdout/stderr open")
        self.assertEqual(result.returncode, 0)
        status_file = self.root / "logs" / "detach.status"
        for _ in range(70):
            if status_file.exists() and status_file.read_text(encoding="utf-8").strip() == "done":
                break
            time.sleep(0.1)
        self.assertEqual(status_file.read_text(encoding="utf-8").strip(), "done")


class GeneHitSelectionTests(unittest.TestCase):
    def setUp(self):
        self.module = load_script("blast_genes", Path("scripts") / "blast_genes.py")

    def test_rejects_partial_gene_hits(self):
        hits = {"COX1": [{"evalue": 1e-40, "identity": 98.0, "aln_len": 50, "query_len": 100, "start": 1, "end": 50, "strand": "+", "type": "CDS"}]}
        selected, errors = self.module.select_gene_hits(hits)
        self.assertEqual(selected, {})
        self.assertTrue(any("COX1" in error for error in errors))

    def test_rejects_ambiguous_gene_loci(self):
        base = {"evalue": 1e-40, "identity": 98.0, "aln_len": 100, "query_len": 100, "start": 1, "end": 100, "strand": "+", "type": "CDS"}
        hits = {"COX1": [base, dict(base, start=500, end=600)]}
        selected, errors = self.module.select_gene_hits(hits)
        self.assertEqual(selected, {})
        self.assertTrue(any("COX1" in error and "多个" in error for error in errors))

    def test_merges_overlapping_hsps_once(self):
        h = {"evalue": 1e-40, "identity": 98.0, "aln_len": 500,
             "query_len": 1000, "qstart": 1, "qend": 900,
             "start": 1, "end": 900, "strand": "+", "type": "CDS"}
        merged = self.module.merge_hsps([h, dict(h, qstart=500, qend=1000, start=500, end=1000)])
        self.assertEqual(merged["aln_len"], 1000)


class CircularizationSelectionTests(unittest.TestCase):
    def setUp(self):
        self.module = load_script("circularize", Path("scripts") / "circularize.py")

    def test_rejects_equal_top_candidates(self):
        self.assertIsNone(self.module.select_unique_candidate([(10, "a"), (10, "b")]))

    def test_selects_unique_top_candidate(self):
        self.assertEqual(self.module.select_unique_candidate([(10, "a"), (9, "b")]), (10, "a"))

    def test_validates_circular_gene_order(self):
        self.assertTrue(self.module.circular_order_matches(["b", "c", "a"], ["a", "b", "c"]))
        self.assertFalse(self.module.circular_order_matches(["a", "c", "b", "d"], ["a", "b", "c", "d"]))

    def test_rejects_inverted_gene_orientation(self):
        self.assertTrue(self.module.circular_annotation_matches(
            [("b", "+"), ("c", "-"), ("a", "+")],
            [("a", "+"), ("b", "+"), ("c", "-")],
        ))
        self.assertFalse(self.module.circular_annotation_matches(
            [("a", "+"), ("b", "-"), ("c", "-")],
            [("a", "+"), ("b", "+"), ("c", "-")],
        ))

    def test_union_and_scaffold_overlap_are_nonredundant(self):
        self.assertEqual(self.module.interval_union_length([(0, 900), (500, 1500)]), 1500)
        joined, overlap = self.module.join_scaffolds("AAACCC", "CCCGGG", minimum_overlap=3)
        self.assertEqual((joined, overlap), ("AAACCCGGG", 3))


class AnnotationNormalizationTests(unittest.TestCase):
    def setUp(self):
        self.module = load_script("annot_check", Path("scripts") / "annot_check.py")

    def test_rrna_names_follow_biological_16s_12s_mapping(self):
        self.assertEqual(self.module.canonical_gene(SimpleNamespace(qualifiers={"gene": ["16S"]})), "rrnl")
        self.assertEqual(self.module.canonical_gene(SimpleNamespace(qualifiers={"gene": ["12S"]})), "rrns")


class DepthMappingTests(unittest.TestCase):
    def setUp(self):
        self.module = load_script("depth_analysis", Path("scripts") / "depth_analysis.py")

    def test_maps_forward_read_through_cigar(self):
        fields = ["read", "0", "chr", "100", "60", "2S3M1I4M", "*", "0", "0", "ACGTACGTACGT", "="]
        self.assertEqual(self.module.read_base_at_reference(fields, 102), "A")

    def test_maps_reverse_read_using_reference_orientation(self):
        fields = ["read", "16", "chr", "100", "60", "10M", "*", "0", "0", "ACGTACGTAC", "="]
        self.assertEqual(self.module.read_base_at_reference(fields, 103), "C")

    def test_zero_depth_fails_quality_gate(self):
        self.assertFalse(self.module.depth_quality_gate(0.0, []))
        self.assertFalse(self.module.depth_quality_gate(10.0, [1]))
        self.assertTrue(self.module.depth_quality_gate(10.0, []))

    def test_coverage_intervals_include_reference_tail(self):
        self.assertEqual(self.module.coverage_intervals(1000, 400), [(1, 400), (401, 800), (801, 1000)])


class AnnotationQualityTests(unittest.TestCase):
    def test_rejects_duplicate_genes_with_correct_counts(self):
        try:
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            from Bio.SeqFeature import FeatureLocation, SeqFeature
            from Bio import SeqIO
        except ImportError:
            self.skipTest("Biopython is not installed")
        sequence = ""
        features = []
        position = 0
        for feature_type, names in (("CDS", ["COX1"] * 13), ("tRNA", ["trnA"] * 22), ("rRNA", ["rrnL"] * 2)):
            for name in names:
                bases = "ATG" + "AAA" * 18 + "TAA" if feature_type == "CDS" else "A" * 60
                features.append(SeqFeature(
                    FeatureLocation(position, position + len(bases), strand=1),
                    type=feature_type,
                    qualifiers={"gene": [name]},
                ))
                sequence += bases
                position += len(bases)
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "duplicate.gb"
            record = SeqRecord(Seq(sequence), id="test", name="test", description="", annotations={"molecule_type": "DNA"})
            record.features = features
            SeqIO.write(record, path, "genbank")
            result = subprocess.run(
                ["python3", str(ROOT / "scripts" / "annot_check.py"), str(path)],
                capture_output=True,
                text=True,
                timeout=10,
            )
        self.assertEqual(result.returncode, 1)
        self.assertIn("基因身份", result.stdout)


if __name__ == "__main__":
    unittest.main()

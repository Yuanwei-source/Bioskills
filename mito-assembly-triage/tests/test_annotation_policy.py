#!/usr/bin/env python3
"""annotation / diagnostic-policy regression tests.

Encodes the six minimum regression scenarios of the references revision
(`references/annotation_quality.md`, `diagnostic-playbook.md`,
`evidence-standard.md`, `standard_gene_order.md`, `learning-policy.md`):

1. a correct genome with a typical CDS-CDS overlap or a credible incomplete
   T/TA stop must not be auto-"fixed";
2. a gene missed by annotation tools must be reported as an annotation finding,
   not as true sequence loss;
3. an annotation-only correction without FASTQ/BAM is allowed, with
   reads-support recorded as NOT_ASSESSED;
4. insufficient repeat-spanning reads or missing unique junction evidence
   must end as UNRESOLVED, without claiming an unverified closure;
5. known exceptions (truncated tRNA, class-specific trnP/trnT order) must not
   be misjudged by a blanket rule;
6. local paths/hashes/sample names and downloaded lesson text must never leak
   into public exports or be executed.
"""
import hashlib
import importlib.util
import json
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]

#: Typical insect (Drosophila-type) arrangement; the anchor documented in
#: references/standard_gene_order.md.  Not a universal rule.
ORDER = [
    ("trnI", "tRNA", 70, "+"), ("trnQ", "tRNA", 68, "-"), ("trnM", "tRNA", 70, "+"),
    ("nad2", "CDS", 300, "+"), ("trnW", "tRNA", 69, "+"), ("trnC", "tRNA", 66, "-"),
    ("trnY", "tRNA", 68, "-"), ("cox1", "CDS", 300, "+"), ("trnL2", "tRNA", 70, "+"),
    ("cox2", "CDS", 300, "+"), ("trnK", "tRNA", 70, "+"), ("trnD", "tRNA", 68, "+"),
    ("atp8", "CDS", 300, "+"), ("atp6", "CDS", 300, "+"), ("cox3", "CDS", 300, "+"),
    ("trnG", "tRNA", 68, "+"), ("nad3", "CDS", 300, "+"), ("trnA", "tRNA", 70, "+"),
    ("trnR", "tRNA", 68, "+"), ("trnN", "tRNA", 70, "+"), ("trnS1", "tRNA", 68, "+"),
    ("trnE", "tRNA", 70, "+"), ("trnF", "tRNA", 67, "-"), ("nad5", "CDS", 300, "-"),
    ("trnH", "tRNA", 68, "-"), ("nad4", "CDS", 300, "-"), ("nad4l", "CDS", 300, "-"),
    ("trnT", "tRNA", 69, "+"), ("trnP", "tRNA", 68, "-"), ("nad6", "CDS", 300, "+"),
    ("cob", "CDS", 300, "+"), ("trnS2", "tRNA", 70, "+"), ("nad1", "CDS", 300, "-"),
    ("trnL1", "tRNA", 70, "-"), ("rrnL", "rRNA", 1300, "-"), ("trnV", "tRNA", 70, "-"),
    ("rrnS", "rRNA", 750, "-"),
]

#: table-5 sense codons usable inside a CDS (stops excluded)
CODON_POOL = [
    a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
    if a + b + c not in ("TAA", "TAG", "AGA", "AGG")
]


def load_module(name, relative_path):
    path = ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def biopython_available():
    try:
        import Bio  # noqa: F401
    except ImportError:
        return False
    return True


def write_gb(path, sequence, features):
    """features: [{name, type, start, end, strand(+-), note?}] with 0-based half-open coords."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    seq_features = []
    for spec in features:
        qualifiers = {"gene": [spec["name"]]}
        if spec.get("note"):
            qualifiers["note"] = [spec["note"]]
        seq_features.append(SeqFeature(
            FeatureLocation(spec["start"], spec["end"], strand=1 if spec["strand"] == "+" else -1),
            type=spec["type"],
            qualifiers=qualifiers,
        ))
    record = SeqRecord(Seq(sequence), id="test", name="test", description="",
                       annotations={"molecule_type": "DNA", "topology": "circular"})
    record.features = seq_features
    SeqIO.write(record, path, "genbank")


def build_synthetic_ref(directory):
    """Write ref.gb plus two overlapping scaffolds (s1/s2.fa) under `directory`.

    Sequences are random but CDS are stop-free CDS in table 5, so the record is a
    positive control for the documented typical-insect arrangement.
    """
    random.seed(7)
    from Bio.Seq import Seq

    blocks, features, position = [], [], 0
    for name, kind, length, strand in ORDER:
        if kind == "CDS":
            block = "ATG" + "".join(random.choice(CODON_POOL) for _ in range(98)) + "TAA"
            # features are extracted in feature orientation: store the minus-strand
            # genes as reverse complements so the extracted CDS is ATG...TAA.
            if strand == "-":
                block = str(Seq(block).reverse_complement())
        else:
            block = "".join(random.choice("ACGT") for _ in range(length))
        assert len(block) == length
        blocks.append(block)
        features.append({"name": name, "type": kind, "start": position, "end": position + length,
                         "strand": strand})
        position += length
    sequence = "".join(blocks)
    ref = Path(directory) / "ref.gb"
    write_gb(ref, sequence, features)
    (Path(directory) / "s1.fa").write_text(">s1\n%s\n" % sequence[:3000], encoding="utf-8")
    (Path(directory) / "s2.fa").write_text(">s2\n%s\n" % sequence[2980:], encoding="utf-8")
    return ref, sequence


def run_annot_check(path, *extra):
    return subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "annot_check.py"), str(path), *extra],
        capture_output=True, text=True, timeout=60,
    )


class OverlapAndStopPolicyTests(unittest.TestCase):
    """Scenario 1: never auto-fix a correct genome (typical overlap / T/TA stop)."""

    def setUp(self):
        if not biopython_available():
            self.skipTest("Biopython is not installed")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def _two_feature_gb(self, second_type, second_name):
        path = self.dir / "overlap.gb"
        write_gb(path, "A" * 500, [
            {"name": "atp8", "type": "CDS", "start": 0, "end": 100, "strand": "+"},
            {"name": second_name, "type": second_type, "start": 90, "end": 190, "strand": "+"},
        ])
        return path

    def test_cds_cds_overlap_is_warning_not_error(self):
        path = self._two_feature_gb("CDS", "atp6")
        result = run_annot_check(path)
        self.assertRegex(result.stdout, r"WARN :\s+atp8.*atp6.*重叠\s+10bp")
        self.assertNotRegex(result.stdout, r"ERROR:.*atp8.*atp6.*重叠")

    def test_cds_trna_overlap_defaults_to_error_and_is_configurable(self):
        path = self._two_feature_gb("tRNA", "trnQ")
        default = run_annot_check(path)
        self.assertRegex(default.stdout, r"ERROR:.*atp8.*trnQ.*重叠\s+10bp")

        downgraded = run_annot_check(path, "--overlap-severity", "warn")
        self.assertRegex(downgraded.stdout, r"WARN :\s+atp8.*trnQ.*重叠\s+10bp")
        self.assertNotRegex(downgraded.stdout, r"ERROR:.*atp8.*trnQ")

        tolerated = run_annot_check(path, "--tolerate-overlap", "atp8,trnQ")
        self.assertRegex(tolerated.stdout, r"ACCEPT:.*atp8.*trnQ.*重叠\s+10bp")

    def test_overlap_is_never_a_reason_to_rewrite_the_sequence(self):
        path = self._two_feature_gb("CDS", "atp6")
        before = hashlib.sha256(path.read_bytes()).hexdigest()
        run_annot_check(path, "--overlap-severity", "warn")
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), before)

    def test_incomplete_stop_with_note_is_warning_not_error(self):
        sequence = "ATG" + "A" * 27 + "A" * 100  # extracted CDS ends with AAA -> incomplete
        path = self.dir / "incomplete.gb"
        write_gb(path, sequence, [
            {"name": "cox1", "type": "CDS", "start": 0, "end": 30, "strand": "+",
             "note": "3' 端不完整, 转录后多聚腺苷酸化"},
        ])
        with_note = run_annot_check(path)
        self.assertIn("已按 note 记录", with_note.stdout)
        self.assertNotRegex(with_note.stdout, r"\[ERROR\].*缺少完整终止密码子")

        no_note = self.dir / "incomplete-nonote.gb"
        write_gb(no_note, sequence, [
            {"name": "cox1", "type": "CDS", "start": 0, "end": 30, "strand": "+"},
        ])
        without_note = run_annot_check(no_note)
        self.assertRegex(without_note.stdout, r"\[ERROR\].*缺少完整终止密码子")


class MissingGenePolicyTests(unittest.TestCase):
    """Scenario 2: a missed gene is an annotation finding, not sequence loss."""

    def setUp(self):
        if not biopython_available():
            self.skipTest("Biopython is not installed")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "missing.gb"
        sequence = "ATG" + "A" * 87 + "TAA" + "A" * 100
        write_gb(self.path, sequence, [
            {"name": "cox1", "type": "CDS", "start": 0, "end": 93, "strand": "+"},
        ])

    def test_missing_gene_is_reported_as_annotation_identity(self):
        result = run_annot_check(self.path)
        self.assertIn("缺少基因身份", result.stdout)
        for forbidden in ("序列缺失", "真实丢失", "序列已丢失"):
            self.assertNotIn(forbidden, result.stdout)

    def test_allow_atypical_downgrades_gene_set_findings(self):
        result = run_annot_check(self.path, "--allow-atypical")
        self.assertIn("待核查(非典型类群)", result.stdout)
        self.assertNotEqual(result.returncode, 1)


class StatusVocabularyTests(unittest.TestCase):
    """Scenario 3 + learning-policy: annotation-only fix without reads; lesson retirement states."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.module = load_module("experience_policy", Path("tools") / "experience.py")
        self.module.KNOWLEDGE = str(self.root / "knowledge")
        self.module.CASES = str(self.root / "knowledge" / "cases")
        self.module.LESSON_CANDIDATES = str(self.root / "knowledge" / "lessons" / "candidates")
        self.module.LESSON_VERIFIED = str(self.root / "knowledge" / "lessons" / "verified")
        self.module.PUBLIC_KNOWLEDGE = str(self.root / "knowledge" / "public")
        self.module.ensure_cases_dir()

    def _case(self, name, status="RESOLVED"):
        case = self.root / name
        case.mkdir(parents=True)
        data = {"schema_version": "2.0", "case_id": name, "inputs": [],
                "issue": {"type": "boundary_error", "user_observation": "cox1 边界"},
                "hypotheses": [{"id": "H1", "explanation": "边界错误",
                                "support": [], "against": [], "unknown": []}],
                "events_file": "events.jsonl",
                "decision": {"status": status, "confidence": "moderate", "rationale": "同源+ORF"},
                "validation": [{"check": "annot_check", "reads_support": "NOT_ASSESSED"}]}
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8")
        (case / "events.jsonl").touch()
        return case

    def test_annotation_only_correction_without_reads_can_be_recorded(self):
        case = self._case("annotation-only")
        self.assertEqual(self.module.case_validate(SimpleNamespace(directory=str(case))), 0)
        data = json.loads((case / "case.json").read_text(encoding="utf-8"))
        self.assertEqual(data["validation"][0]["reads_support"], "NOT_ASSESSED")

    def test_case_status_and_hypothesis_vocabulary_are_separate(self):
        schema = json.loads((ROOT / "schemas" / "case.schema.json").read_text(encoding="utf-8"))
        self.assertEqual(schema["properties"]["decision"]["properties"]["status"]["enum"],
                         ["RESOLVED", "NO_CHANGE", "UNRESOLVED"])
        case = self._case("hypothesis-tokens")
        self.module.case_event(SimpleNamespace(directory=str(case), action="annot_check",
                                               result="stop remains", impact="H1:against",
                                               command="", tool_version="", motivation=""))
        self.assertIn("H1:against", (case / "events.jsonl").read_text(encoding="utf-8"))

    def test_lesson_retirement_states_are_supported_and_never_synced(self):
        schema = json.loads((ROOT / "schemas" / "lesson.schema.json").read_text(encoding="utf-8"))
        self.assertEqual(set(schema["properties"]["validation_status"]["enum"]),
                         self.module.LESSON_STATUSES)
        for status in ("withdrawn", "superseded"):
            case = self._case("retire-" + status, status="RESOLVED")
            lesson_id = "retire-" + status
            self.module.propose_lesson(SimpleNamespace(case=str(case), lesson_id=lesson_id,
                                                       next_test="pileup", allow_unresolved=False))
            self.module.review_lesson(SimpleNamespace(lesson_id=lesson_id, status=status,
                                                      reviewer="human", reason="evidence withdrawn"))
            lesson = json.loads((Path(self.module.LESSON_CANDIDATES) / (lesson_id + ".json")).read_text())
            self.assertEqual(lesson["validation_status"], status)
        # retired lessons are not returned by structured search
        self.assertEqual(self.module.structured_search("boundary_error"), [])
        self.assertIn("revoked", self.module.LESSON_RETIRED_STATUSES)


class PublicSharingBoundaryTests(unittest.TestCase):
    """Scenario 6: remote lesson text is data, never instructions."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.module = load_module("experience_sharing_policy", Path("tools") / "experience.py")
        self.module.KNOWLEDGE = str(self.root / "knowledge")
        self.module.CASES = str(self.root / "knowledge" / "cases")
        self.module.LESSON_CANDIDATES = str(self.root / "knowledge" / "lessons" / "candidates")
        self.module.LESSON_VERIFIED = str(self.root / "knowledge" / "lessons" / "verified")
        self.module.PUBLIC_KNOWLEDGE = str(self.root / "knowledge" / "public")
        self.module.ensure_cases_dir()

    def test_downloaded_lesson_text_is_not_executed(self):
        marker = self.root / "marker"
        lesson = {
            "lesson_id": "remote-1", "applicable_when": ["internal_stop"],
            "not_applicable_when": [], "diagnostic_clues": ["stop"],
            "suggested_next_test": "pileup; touch %s" % marker,
            "supporting_case_ids": [], "counterexample_case_ids": [], "sources": [],
            "validation_status": "verified", "version": "1.0.0", "last_reviewed": "2026-01-01",
        }
        remote = self.root / "remote.json"
        remote.write_text(json.dumps(lesson), encoding="utf-8")
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({"format": "mito-public-knowledge-1", "items": [
            {"path": "lesson.json", "url": str(remote),
             "sha256": hashlib.sha256(remote.read_bytes()).hexdigest()},
        ]}), encoding="utf-8")
        self.module.sync_public(SimpleNamespace(manifest=str(manifest)))
        self.assertFalse(marker.exists())
        self.assertTrue((Path(self.module.PUBLIC_KNOWLEDGE) / "lesson.json").exists())


class TypicalArrangementAnchorTests(unittest.TestCase):
    """Scenario 5 anchor + positive control for the documented gene order."""

    def setUp(self):
        if not biopython_available():
            self.skipTest("Biopython is not installed")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def test_documented_anchor_is_self_consistent(self):
        counts = {"CDS": 0, "tRNA": 0, "rRNA": 0}
        for _, kind, _, _ in ORDER:
            counts[kind] += 1
        self.assertEqual(counts, {"CDS": 13, "tRNA": 22, "rRNA": 2})
        plus = sum(1 for _, kind, _, strand in ORDER if kind == "CDS" and strand == "+")
        self.assertEqual(plus, 9)

    def test_documented_anchor_has_no_hard_annotation_error(self):
        ref, _ = build_synthetic_ref(self.dir)
        result = run_annot_check(ref, "--require-circular")
        self.assertNotEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertNotIn("[ERROR]", result.stdout)


class TrnaExceptionTests(unittest.TestCase):
    """Scenario 5: a truncated (but real) tRNA must not be hard-failed."""

    def setUp(self):
        if not biopython_available():
            self.skipTest("Biopython is not installed")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def _short_trna(self, note):
        path = self.dir / ("short-%s.gb" % ("note" if note else "plain"))
        spec = {"name": "trnS1", "type": "tRNA", "start": 0, "end": 45, "strand": "+"}
        if note:
            spec["note"] = "trnS1 缺 DHU 臂, 结构与同源证据一致"
        write_gb(path, "A" * 200, [spec])
        return path

    def test_truncated_trna_with_note_is_warning(self):
        result = run_annot_check(self._short_trna(True))
        self.assertIn("长度 45bp <50", result.stdout)
        self.assertNotRegex(result.stdout, r"\[ERROR\].*tRNA trnS1 长度")

    def test_truncated_trna_without_note_is_error(self):
        result = run_annot_check(self._short_trna(False))
        self.assertRegex(result.stdout, r"\[ERROR\].*tRNA trnS1 长度 45bp 异常短")


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
@unittest.skipIf(any(shutil.which(tool) is None for tool in ("blastn", "makeblastdb", "minimap2")),
                 "blastn/makeblastdb/minimap2 are required")
class CircularizationEvidenceTests(unittest.TestCase):
    """Scenario 4: no unique junction evidence -> UNRESOLVED, never a claimed closure."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)
        self.ref, self.sequence = build_synthetic_ref(self.dir)

    def _run(self, outdir, *extra):
        return subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "circularize.py"),
             str(self.dir / "s1.fa"), str(self.dir / "s2.fa"), str(self.ref), str(outdir), *extra],
            capture_output=True, text=True, timeout=600,
        )

    def test_missing_reads_evidence_yields_unresolved(self):
        result = self._run(self.dir / "out")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("UNRESOLVED", result.stdout + result.stderr)
        self.assertNotIn("CANDIDATE_ACCEPTED", result.stdout)

    def test_insufficient_spanning_reads_are_rejected(self):
        if shutil.which("samtools") is None:
            self.skipTest("samtools is not available")
        first = self._run(self.dir / "out-bam")
        candidate = self.dir / "out-bam" / "genome_candidate.fasta"
        if not candidate.exists():
            self.skipTest("circularize did not reach the candidate stage: %s" % first.stderr)
        from Bio import SeqIO
        length = len(next(SeqIO.parse(candidate, "fasta")).seq)
        sam = self.dir / "reads.sam"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:mitogenome_candidate\tLN:%d\n"
            "r1\t0\tmitogenome_candidate\t100\t60\t50M\t*\t0\t0\t%s\t*\n"
            "r2\t0\tmitogenome_candidate\t200\t60\t50M\t*\t0\t0\t%s\t*\n"
            % (length, "A" * 50, "A" * 50), encoding="utf-8",
        )
        bam = self.dir / "reads.bam"
        subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)],
                       capture_output=True, text=True, check=True, timeout=120)
        subprocess.run(["samtools", "index", str(bam)],
                       capture_output=True, text=True, check=True, timeout=120)
        result = self._run(self.dir / "out-bam2", "--bam", str(bam), "--reads-validated")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("需要至少 3 条", result.stderr)
        self.assertNotIn("CANDIDATE_ACCEPTED", result.stdout)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Failing reproductions for the evidence-contract defects:

P8  JSON Schema is declared but never enforced by ``case-validate``
P9  circularize treats a reference order mismatch as a hard failure (exit 1)
P10 read-level evidence (MAPQ / base quality / strand / duplicates / depth) is
    not produced, so the strongest claims are not backed by the documented fields
P11 cox1_id reports a single hit with no query coverage and prints unfounded
    species-level certainty tiers
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gb_fixtures import (  # noqa: E402
    ROOT, biopython_available, cds_block, have, load_module, run_script, write_gb,
)


class CaseSchemaTests(unittest.TestCase):
    """P8: schemas/case.schema.json must actually be enforced."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.module = load_module("experience_schema", Path("tools") / "experience.py")
        self.module.KNOWLEDGE = str(self.root / "knowledge")

    def _case(self, name, **overrides):
        case = self.root / name
        case.mkdir(parents=True)
        data = {
            "schema_version": "2.0", "case_id": name,
            "inputs": [{"role": "assembly_fasta", "path": "/tmp/x.fa", "sha256": "a" * 64}],
            "issue": {"type": "internal_stop", "user_observation": "nad5 stop"},
            "hypotheses": [{"id": "H1", "explanation": "边界错误",
                            "support": [], "against": [], "unknown": []}],
            "events_file": "events.jsonl",
            "decision": {"status": "RESOLVED", "confidence": "moderate", "rationale": "同源"},
        }
        data.update(overrides)
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8")
        (case / "events.jsonl").touch()
        return case

    def validate(self, case):
        return self.module.case_validate(SimpleNamespace(directory=str(case)))

    def test_valid_case_passes(self):
        self.assertEqual(self.validate(self._case("ok")), 0)

    def test_missing_hypotheses_is_invalid(self):
        case = self._case("no-hypotheses")
        data = json.loads((case / "case.json").read_text(encoding="utf-8"))
        del data["hypotheses"]
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8")
        self.assertEqual(self.validate(case), 1)

    def test_empty_hypotheses_is_invalid(self):
        self.assertEqual(self.validate(self._case("empty-hypotheses", hypotheses=[])), 1)

    def test_unknown_decision_status_is_invalid(self):
        case = self._case("bad-status")
        data = json.loads((case / "case.json").read_text(encoding="utf-8"))
        data["decision"]["status"] = "DONE"
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8")
        self.assertEqual(self.validate(case), 1)

    def test_malformed_input_entry_is_invalid(self):
        case = self._case("bad-input", inputs=[{"role": "assembly_fasta"}])
        self.assertEqual(self.validate(case), 1)

    def test_per_anomaly_status_and_confidence_are_validated(self):
        good = self._case("anomalies-ok", anomalies=[
            {"id": "A1", "claim": "cox1 非典型起始", "status": "UNRESOLVED",
             "confidence": "not_assessable", "reads_support": "NOT_ASSESSED"},
            {"id": "A2", "claim": "trnS1 缺 DHU 臂属真实特征", "status": "RESOLVED",
             "confidence": "moderate", "reads_support": "NOT_ASSESSED"},
        ])
        self.assertEqual(self.validate(good), 0)
        bad_status = self._case("anomalies-bad-status", anomalies=[
            {"id": "A1", "claim": "x", "status": "DONE", "confidence": "high"}])
        self.assertEqual(self.validate(bad_status), 1)
        bad_reads = self._case("anomalies-bad-reads", anomalies=[
            {"id": "A1", "claim": "x", "status": "RESOLVED", "confidence": "high",
             "reads_support": "raw-read-supported"}])
        self.assertEqual(self.validate(bad_reads), 1)

    def test_case_init_requires_at_least_one_hypothesis(self):
        with self.assertRaises(ValueError):
            self.module.case_init(SimpleNamespace(
                directory=str(self.root / "init-none"), case_id=None, issue="internal_stop",
                observation="obs", taxon=None, input=[], hypothesis=[]))
        self.module.case_init(SimpleNamespace(
            directory=str(self.root / "init-one"), case_id=None, issue="internal_stop",
            observation="obs", taxon=None, input=[], hypothesis=["H1 边界错误"]))
        self.assertEqual(self.validate(self.root / "init-one"), 0)


class JunctionEvidenceTests(unittest.TestCase):
    """P10: junction support must be molecule/quality evidence, not a bare count."""

    def setUp(self):
        self.module = load_module("circularize_junction", Path("scripts") / "circularize.py")

    @staticmethod
    def _read(name, flag=0, mapq=60, pos=80, cigar="50M"):
        return "\t".join([name, str(flag), "mitogenome_candidate", str(pos), str(mapq), cigar,
                          "*", "0", "0", "A" * 50, "I" * 50])

    def test_duplicates_are_excluded_and_strands_counted(self):
        lines = [
            self._read("fwd", flag=0),
            self._read("rev", flag=16),
            self._read("dup", flag=1024),          # 0x400 PCR/optical duplicate
            self._read("secondary", flag=256),
            self._read("lowmapq", mapq=5),
        ]
        with patch.object(self.module.subprocess, "run",
                          return_value=SimpleNamespace(stdout="\n".join(lines))):
            evidence = self.module.junction_evidence("x.bam", "mitogenome_candidate:90-110",
                                                     min_mapq=20)
        self.assertEqual(evidence["support"], 2)
        self.assertEqual(evidence["duplicates_excluded"], 1)
        self.assertEqual(evidence["strand_counts"]["+"], 1)
        self.assertEqual(evidence["strand_counts"]["-"], 1)


@unittest.skipUnless(biopython_available(), "Biopython is not installed")
@unittest.skipUnless(have("blastn", "makeblastdb", "minimap2"), "blastn/minimap2 required")
class CircularizeReviewTests(unittest.TestCase):
    """P9: a reference order mismatch is REVIEW, not a hard failure."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)

    def _swapped_scaffold_case(self):
        names = ["geneA", "geneB", "geneC", "geneD"]
        blocks = {name: cds_block(300, seed=31 + index) for index, name in enumerate(names)}
        sequence = "".join(blocks[name] for name in names)
        features = [{"gene": name, "type": "CDS", "start": index * 300,
                     "end": (index + 1) * 300, "strand": "+"}
                    for index, name in enumerate(names)]
        ref = write_gb(self.dir / "ref.gb", sequence, features)
        (self.dir / "s1.fa").write_text(">s1\n%s\n" % (blocks["geneA"] + blocks["geneB"]),
                                        encoding="utf-8")
        # C and D are swapped on the second scaffold -> candidate order A,B,D,C
        (self.dir / "s2.fa").write_text(">s2\n%s\n" % (blocks["geneD"] + blocks["geneC"]),
                                        encoding="utf-8")
        return ref

    def test_order_mismatch_returns_review_and_keeps_the_candidate(self):
        ref = self._swapped_scaffold_case()
        outdir = self.dir / "out"
        result = run_script(Path("scripts") / "circularize.py",
                            self.dir / "s1.fa", self.dir / "s2.fa", ref, outdir)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertIn("REVIEW", result.stdout + result.stderr)
        self.assertNotIn("CANDIDATE_ACCEPTED", result.stdout)
        self.assertTrue((outdir / "genome_candidate.fasta").exists())


class BaseSupportTests(unittest.TestCase):
    """P10: single-base replacement needs MAPQ, base quality, strand and depth."""

    def setUp(self):
        self.module = load_module("depth_base_support", Path("scripts") / "depth_analysis.py")

    @staticmethod
    def _read(name, flag=0, mapq=60, seq="AAAAAAAAAA", qual="IIIIIIIIII", pos=100):
        return [name, str(flag), "chr", str(pos), str(mapq), "10M", "*", "0", "0", seq, qual]

    def _support(self, records, **overrides):
        kwargs = {"min_mapq": 20, "min_baseq": 20, "min_depth": 5, "max_strand_fraction": 0.9}
        kwargs.update(overrides)
        return self.module.base_support(records, 103, **kwargs)

    def test_shallow_evidence_is_not_callable(self):
        evidence = self._support([self._read("r1"), self._read("r2")])
        self.assertFalse(evidence["callable"])
        self.assertIn("深度", " ".join(evidence["reasons"]))

    def test_strand_biased_evidence_is_not_callable(self):
        records = [self._read("f%d" % i) for i in range(6)]
        evidence = self._support(records)
        self.assertFalse(evidence["callable"])
        self.assertIn("链向", " ".join(evidence["reasons"]))

    def test_balanced_deep_evidence_is_callable(self):
        records = [self._read("f%d" % i) for i in range(4)]
        records += [self._read("r%d" % i, flag=16, seq="TTTTTTTTTT") for i in range(4)]
        evidence = self._support(records)
        self.assertTrue(evidence["callable"], evidence)
        self.assertEqual(evidence["base"], "A")

    def test_duplicates_and_low_mapq_are_excluded(self):
        records = [self._read("f%d" % i) for i in range(4)]
        records += [self._read("r%d" % i, flag=16, seq="TTTTTTTTTT") for i in range(4)]
        records += [self._read("dup", flag=1024), self._read("low", mapq=3)]
        evidence = self._support(records)
        self.assertEqual(evidence["excluded"]["duplicate"], 1)
        self.assertEqual(evidence["excluded"]["low_mapq"], 1)
        self.assertTrue(evidence["callable"], evidence)


class Cox1VerdictTests(unittest.TestCase):
    """P11: query coverage, ambiguity, and no unfounded species certainty."""

    def setUp(self):
        self.module = load_module("cox1_verdict", Path("scripts") / "cox1_id.py")

    def test_low_coverage_hits_are_rejected(self):
        hits = [(99.0, 400, "short hit"), (95.0, 900, "full hit")]
        best = self.module.select_supported_hit(hits, min_identity=85.0,
                                                min_coverage=0.8, query_len=1000)
        self.assertIsNotNone(best)
        self.assertEqual(best[2], "full hit")

    def test_only_low_coverage_hits_gives_no_verdict(self):
        hits = [(99.0, 400, "short hit")]
        self.assertIsNone(self.module.select_supported_hit(hits, min_identity=85.0,
                                                           min_coverage=0.8, query_len=1000))

    def test_verdict_never_claims_species(self):
        for identity, coverage in ((99.9, 0.99), (97.5, 0.95), (90.0, 0.9), (80.0, 0.5)):
            verdict = self.module.interpret(identity, coverage)
            self.assertIn(verdict["status"],
                          ("provisional_candidate", "ambiguous", "insufficient"))
            self.assertNotIn("同种", verdict["message"])

    def test_unfounded_certainty_tiers_are_gone(self):
        source = (ROOT / "scripts" / "cox1_id.py").read_text(encoding="utf-8")
        self.assertNotIn(">97%", source)
        self.assertNotIn("同属不同种", source)
        self.assertIn("query coverage", source.lower())


def _fixed_case(**overrides):
    """A valid case dict; overrides produce the fixed equivalence data set."""
    case = {
        "schema_version": "2.0", "case_id": "fixed-case",
        "inputs": [{"role": "assembly_fasta", "path": "/tmp/x.fa", "sha256": "a" * 64}],
        "issue": {"type": "internal_stop", "user_observation": "nad5 stop"},
        "hypotheses": [{"id": "H1", "explanation": "frame/boundary",
                        "support": [], "against": [], "unknown": []}],
        "events_file": "events.jsonl",
        "decision": {"status": "RESOLVED", "confidence": "moderate", "rationale": "homology"},
        "anomalies": [{"id": "A1", "claim": "boundary fix", "status": "RESOLVED",
                       "confidence": "moderate", "reads_support": "NOT_ASSESSED"}],
    }
    case.update(overrides)
    return case


class SchemaValidatorEquivalenceTests(unittest.TestCase):
    """The jsonschema path and the no-dependency fallback must follow ONE schema.

    Technical debt being managed here: the two validators are separate code, so a
    schema change can drift them apart silently.  This fixed data set pins the
    verdicts; ``test_both_paths_agree`` additionally requires identical results
    wherever jsonschema is installed.  If a future schema edit only touches one
    implementation, one of these two tests fails.

    Note: a cross-field contradiction (a case-level RESOLVED decision next to an
    UNRESOLVED anomaly) is deliberately NOT rejected by either validator -- the
    schema does not express it, and adding the rule to only one implementation
    would create exactly the divergence this test guards against.
    """

    FIXED_CASES = {
        "legal": (_fixed_case(), True),
        "missing_field": (_fixed_case(decision=None), False),
        "wrong_field_type": (_fixed_case(inputs="not-an-array"), False),
        "illegal_enum": (_fixed_case(decision={"status": "DONE", "confidence": "moderate",
                                                "rationale": "x"}), False),
        "illegal_nested_enum": (_fixed_case(anomalies=[
            {"id": "A1", "claim": "x", "status": "RESOLVED", "confidence": "high",
             "reads_support": "raw-read-supported"}]), False),
        "conflicting_status": (_fixed_case(anomalies=[
            {"id": "A1", "claim": "x", "status": "UNRESOLVED",
             "confidence": "not_assessable"}]), True),
        "nested_object_error": (_fixed_case(anomalies=[
            {"id": "A1", "status": "RESOLVED", "confidence": "high"}]), False),
        "empty_hypotheses": (_fixed_case(hypotheses=[]), False),
        "inputs_missing_sha256": (_fixed_case(inputs=[{"role": "assembly_fasta"}]), False),
    }

    def setUp(self):
        self.module = load_module("experience_schema_equiv", Path("tools") / "experience.py")

    def test_fallback_verdicts_are_fixed(self):
        for name, (case, expected_valid) in self.FIXED_CASES.items():
            with self.subTest(case=name):
                errors = self.module._case_errors_without_jsonschema(case)
                self.assertEqual(not errors, expected_valid,
                                 "%s -> %s" % (name, errors))

    def test_both_paths_agree(self):
        try:
            import jsonschema
        except ImportError:
            self.skipTest("jsonschema is not installed; the fallback is covered by "
                          "test_fallback_verdicts_are_fixed")
        schema = json.loads((ROOT / "schemas" / "case.schema.json").read_text(encoding="utf-8"))
        for name, (case, _expected) in self.FIXED_CASES.items():
            with self.subTest(case=name):
                try:
                    jsonschema.validate(case, schema)
                    schema_ok = True
                except jsonschema.ValidationError:
                    schema_ok = False
                fallback_ok = not self.module._case_errors_without_jsonschema(case)
                self.assertEqual(schema_ok, fallback_ok,
                                 "%s: jsonschema=%s fallback=%s" % (name, schema_ok, fallback_ok))

    def test_public_entry_point_matches_the_fixed_verdicts(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        self.module.KNOWLEDGE = str(root / "knowledge")
        for index, (name, (case, expected_valid)) in enumerate(self.FIXED_CASES.items()):
            with self.subTest(case=name):
                directory = root / ("case-%d" % index)
                directory.mkdir()
                data = dict(case)
                data.pop("decision", None) if case.get("decision") is None else None
                (directory / "case.json").write_text(json.dumps(data), encoding="utf-8")
                (directory / "events.jsonl").touch()
                verdict = self.module.case_validate(SimpleNamespace(directory=str(directory)))
                self.assertEqual(verdict == 0, expected_valid, "%s" % name)


if __name__ == "__main__":
    unittest.main()

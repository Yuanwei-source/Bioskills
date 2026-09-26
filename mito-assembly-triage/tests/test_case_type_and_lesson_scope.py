#!/usr/bin/env python3
"""Failing reproductions for two experience-system gaps found in the expert review.

  case_type   A learning system that only stores "abnormal" cases drifts towards the
              belief that every mitogenome is broken (specificity bias).  There was no
              way to record a case as `normal_validation_case` (checked, nothing wrong)
              or `tool_failure_case` (the failure was in the tool/environment).
  lesson scope  A lesson drawn from a single case could be written without any stated
              limit on how far it transfers, so "nad6 split observed once" can be read
              as "insect nad6 is often split".  `generalization_scope` / `transferability`
              did not exist, and there was no fail-closed rule tying them to the number of
              INDEPENDENT supporting cases.

All data here is synthetic.
"""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
EXPERIENCE = ROOT / "tools" / "experience.py"


def load_module(name="experience_scope"):
    spec = importlib.util.spec_from_file_location(name, EXPERIENCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CaseTypeCliTests(unittest.TestCase):
    """case_type must be recordable and validated, including the normal case."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.tmp = Path(temp.name)
        self.env = dict(os.environ)
        self.env["MITO_KNOWLEDGE_DIR"] = str(self.tmp / "knowledge")

    def _init(self, name, *extra):
        directory = self.tmp / name
        result = subprocess.run(
            [sys.executable, str(EXPERIENCE), "case-init", str(directory),
             "--case-id", name, "--issue", "internal_stop",
             "--observation", "synthetic observation",
             "--hypothesis", "H1: boundary error", *extra],
            capture_output=True, text=True, timeout=180, env=self.env)
        return result, directory

    def _case(self, directory):
        return json.loads((directory / "case.json").read_text(encoding="utf-8"))

    def test_default_case_type_is_abnormal(self):
        result, directory = self._init("c-default")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self._case(directory)["case_type"], "abnormal_case")

    def test_normal_validation_case_is_recorded(self):
        result, directory = self._init("c-normal", "--case-type", "normal_validation_case")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self._case(directory)["case_type"], "normal_validation_case")
        validate = subprocess.run(
            [sys.executable, str(EXPERIENCE), "case-validate", str(directory)],
            capture_output=True, text=True, timeout=180, env=self.env)
        self.assertEqual(validate.returncode, 0, validate.stdout + validate.stderr)

    def test_tool_failure_case_is_recorded(self):
        result, directory = self._init("c-tool", "--case-type", "tool_failure_case")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self._case(directory)["case_type"], "tool_failure_case")

    def test_invalid_case_type_is_rejected_without_writing(self):
        result, directory = self._init("c-bad", "--case-type", "anomalous")
        self.assertNotEqual(result.returncode, 0, "非法 case_type 必须被拒绝")
        self.assertFalse((directory / "case.json").exists(), "失败时不得留下案例")

    def test_both_validation_paths_agree_on_case_type(self):
        module = load_module("experience_case_type")
        schema = json.loads((ROOT / "schemas" / "case.schema.json").read_text(encoding="utf-8"))
        self.assertIn("case_type", schema["properties"], "schema 必须声明 case_type")
        self.assertEqual(set(schema["properties"]["case_type"]["enum"]),
                         set(module.CASE_TYPES))

        def base():
            return {
                "schema_version": "2.0", "case_id": "case-1", "inputs": [],
                "issue": {"type": "internal_stop"},
                "hypotheses": [{"id": "H1", "explanation": "e", "support": [], "against": [],
                                "unknown": []}],
                "decision": {"status": "RESOLVED", "confidence": "moderate", "rationale": "r"},
                "events_file": "events.jsonl",
            }

        bad = base()
        bad["case_type"] = "anomalous"
        self.assertTrue(module.case_schema_errors(bad), "jsonschema 路径必须拒绝非法 case_type")
        self.assertTrue(module._case_errors_without_jsonschema(bad),
                        "兜底路径必须拒绝非法 case_type（两条路径必须等价）")
        for value in module.CASE_TYPES:
            good = base()
            good["case_type"] = value
            self.assertEqual(module.case_schema_errors(good), [], value)
            self.assertEqual(module._case_errors_without_jsonschema(good), [], value)


class LessonScopeTests(unittest.TestCase):
    """generalization_scope / transferability must be fail-closed."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.module = load_module("experience_lesson_scope")
        self.module.KNOWLEDGE = str(self.root / "knowledge")
        self.module.CASES = str(self.root / "knowledge" / "cases")
        self.module.LESSON_CANDIDATES = str(self.root / "knowledge" / "lessons" / "candidates")
        self.module.LESSON_VERIFIED = str(self.root / "knowledge" / "lessons" / "verified")
        self.module.PUBLIC_KNOWLEDGE = str(self.root / "knowledge" / "public")

    def make_case(self, name, taxon=None, status="RESOLVED"):
        case = self.root / name
        case.mkdir(parents=True, exist_ok=True)
        data = {"schema_version": "2.0", "case_id": name, "inputs": [],
                "issue": {"type": "internal_stop", "user_observation": "same clue"},
                "hypotheses": [{"id": "H1", "explanation": "e", "support": [], "against": [],
                                "unknown": []}],
                "events_file": "events.jsonl",
                "taxon": {"name": taxon, "taxid": None, "genetic_code": None},
                "decision": {"status": status, "confidence": "moderate", "rationale": "evidence"}}
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8")
        (case / "events.jsonl").touch()
        return case

    def _propose(self, case, lesson_id, **extra):
        args = SimpleNamespace(case=str(case), lesson_id=lesson_id, next_test="pileup",
                               allow_unresolved=False, **extra)
        self.module.propose_lesson(args)
        return json.loads((Path(self.module.LESSON_CANDIDATES)
                           / (lesson_id + ".json")).read_text(encoding="utf-8"))

    def _lesson_path(self, lesson_id):
        return Path(self.module.LESSON_CANDIDATES) / (lesson_id + ".json")

    def test_default_scope_is_single_case_and_no_transferability(self):
        case = self.make_case("case-a", taxon="Species one")
        lesson = self._propose(case, "lesson-default")
        self.assertEqual(lesson["generalization_scope"], "single_case")
        self.assertEqual(lesson["transferability"], "none")

    def test_single_case_cannot_claim_transferability(self):
        case = self.make_case("case-a", taxon="Species one")
        with self.assertRaises(ValueError):
            self._propose(case, "lesson-bold", generalization_scope="order",
                          transferability="high")
        self.assertFalse(self._lesson_path("lesson-bold").exists(),
                         "被拒绝时不得写出候选经验")

    def test_scope_beyond_single_case_requires_independent_cases(self):
        case = self.make_case("case-a", taxon="Species one")
        with self.assertRaises(ValueError):
            self._propose(case, "lesson-family", generalization_scope="family",
                          transferability="low")
        self.assertFalse(self._lesson_path("lesson-family").exists())

    def test_same_taxon_cases_are_not_independent(self):
        first = self.make_case("case-a", taxon="Species one")
        second = self.make_case("case-b", taxon="Species one")
        with self.assertRaises(ValueError):
            self._propose(first, "lesson-same-taxon", supporting_case=[str(second)],
                          generalization_scope="family", transferability="low")
        self.assertFalse(self._lesson_path("lesson-same-taxon").exists())

    def test_independent_taxa_support_low_transferability(self):
        first = self.make_case("case-a", taxon="Species one")
        second = self.make_case("case-b", taxon="Species two")
        lesson = self._propose(first, "lesson-independent", supporting_case=[str(second)],
                               generalization_scope="family", transferability="low")
        self.assertEqual(lesson["generalization_scope"], "family")
        self.assertEqual(lesson["transferability"], "low")
        self.assertEqual(sorted(lesson["supporting_case_ids"]), ["case-a", "case-b"])

    def test_missing_taxon_cannot_establish_independence(self):
        first = self.make_case("case-a", taxon=None)
        second = self.make_case("case-b", taxon=None)
        with self.assertRaises(ValueError):
            self._propose(first, "lesson-no-taxon", supporting_case=[str(second)],
                          generalization_scope="family", transferability="low")
        self.assertFalse(self._lesson_path("lesson-no-taxon").exists())

    def test_counterexample_blocks_transferability(self):
        first = self.make_case("case-a", taxon="Species one")
        second = self.make_case("case-b", taxon="Species two")
        counter = self.make_case("case-c", taxon="Species three")
        with self.assertRaises(ValueError):
            self._propose(first, "lesson-counter", supporting_case=[str(second)],
                          counterexample_case=[str(counter)],
                          generalization_scope="family", transferability="low")
        lesson = self._propose(first, "lesson-counter-recorded",
                               supporting_case=[str(second)],
                               counterexample_case=[str(counter)])
        self.assertEqual(lesson["transferability"], "none")
        self.assertEqual(lesson["counterexample_case_ids"], ["case-c"])

    def test_enum_values_are_validated(self):
        case = self.make_case("case-a", taxon="Species one")
        with self.assertRaises(ValueError):
            self._propose(case, "lesson-bad-scope", generalization_scope="kingdom")
        with self.assertRaises(ValueError):
            self._propose(case, "lesson-bad-transfer", transferability="certain")
        self.assertFalse(self._lesson_path("lesson-bad-scope").exists())
        self.assertFalse(self._lesson_path("lesson-bad-transfer").exists())

    def test_case_cannot_be_both_support_and_counterexample(self):
        case = self.make_case("case-a", taxon="Species one")
        with self.assertRaises(ValueError):
            self._propose(case, "lesson-self", counterexample_case=[str(case)])
        self.assertFalse(self._lesson_path("lesson-self").exists())

    def test_lesson_schema_declares_the_two_fields(self):
        schema = json.loads((ROOT / "schemas" / "lesson.schema.json").read_text(encoding="utf-8"))
        self.assertEqual(set(schema["properties"]["generalization_scope"]["enum"]),
                         set(self.module.GENERALIZATION_SCOPES))
        self.assertEqual(set(schema["properties"]["transferability"]["enum"]),
                         set(self.module.TRANSFERABILITY_LEVELS))

    def test_proposed_lesson_still_has_the_required_fields(self):
        case = self.make_case("case-a", taxon="Species one")
        lesson = self._propose(case, "lesson-required")
        schema = json.loads((ROOT / "schemas" / "lesson.schema.json").read_text(encoding="utf-8"))
        for key in schema["required"]:
            self.assertIn(key, lesson, "缺少必需字段: %s" % key)


if __name__ == "__main__":
    unittest.main()

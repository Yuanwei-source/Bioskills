#!/usr/bin/env python3
"""Failing reproductions for the Copilot Code Review findings on PR #8.

Findings covered here (the documentation-only ones are fixed in the same PR):

  R-1  `conclusion-report.md` describes a report contract (anomalies[], evidence matrix,
       the mandated four sections) that `case-report` does not produce -- the command
       writes only the event list, the case-level decision and a generic note.
  R-3  `learning-policy.md` §2 restricts `tool_failure_case` to tool/environment lessons,
       but `propose_lesson` accepts such a case exactly like an abnormal one.
  R-5  `generalization_scope` / `transferability` are optional in the lesson schema and
       `validate_public_lesson()` does not require them, so a verified or imported lesson
       can bypass the fail-closed single_case/none contract.
  R-6  The independence count used distinct taxon names only, so the same `case_id`
       supplied twice with different `taxon` values was counted as several independent
       cases.
  R-7  `case_type` is not carried into the lesson, so the origin of a lesson is lost.
  R-8  `export_contribution()` builds `safe_case` from a field allowlist that omits
       `case_type`, so the classification silently disappears at the sharing boundary.

All data here is synthetic.
"""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]


def load_module(name="experience_review_round"):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / "experience.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _Base(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.module = load_module(self.__class__.__name__)
        self.module.KNOWLEDGE = str(self.root / "knowledge")
        self.module.CASES = str(self.root / "knowledge" / "cases")
        self.module.LESSON_CANDIDATES = str(self.root / "knowledge" / "lessons" / "candidates")
        self.module.LESSON_VERIFIED = str(self.root / "knowledge" / "lessons" / "verified")
        self.module.PUBLIC_KNOWLEDGE = str(self.root / "knowledge" / "public")
        for name in ("knowledge", "knowledge/cases", "knowledge/lessons/candidates",
                     "knowledge/lessons/verified", "knowledge/public"):
            (self.root / name).mkdir(parents=True, exist_ok=True)

    def make_case(self, name, taxon="Species one", status="RESOLVED", case_type="abnormal_case",
                  anomalies=None, events=()):
        case = self.root / name
        case.mkdir(parents=True, exist_ok=True)
        data = {"schema_version": "2.0", "case_id": name, "case_type": case_type, "inputs": [],
                "issue": {"type": "internal_stop", "user_observation": "same clue"},
                "hypotheses": [{"id": "H1", "explanation": "boundary error", "support": [],
                                "against": [], "unknown": ["needs pileup"]}],
                "events_file": "events.jsonl",
                "taxon": {"name": taxon, "taxid": None, "genetic_code": None},
                "anomalies": anomalies or [],
                "decision": {"status": status, "confidence": "moderate", "rationale": "evidence"}}
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8")
        with open(case / "events.jsonl", "w", encoding="utf-8") as fh:
            for action, result, impact in events:
                fh.write(json.dumps({"action": action, "result": result, "impact": impact}) + "\n")
        return case

    def lesson(self, lesson_id):
        return json.loads((Path(self.module.LESSON_CANDIDATES) / (lesson_id + ".json")).read_text(encoding="utf-8"))

    def lesson_exists(self, lesson_id):
        return (Path(self.module.LESSON_CANDIDATES) / (lesson_id + ".json")).exists()


class CaseReportLayerTests(_Base):
    """R-1: the documented report layer must be what case-report actually writes."""

    def setUp(self):
        super().setUp()
        self.case = self.make_case(
            "case-report", anomalies=[
                {"id": "A-01", "claim": "49 bp zero-coverage gap", "status": "UNRESOLVED",
                 "confidence": "high", "reads_support": "READS_CONSISTENT",
                 "evidence_events": ["depth_analysis"]},
                {"id": "A-02", "claim": "tRNA trnI not detected", "status": "UNRESOLVED",
                 "confidence": "low"},
                {"id": "A-03", "claim": "gene boundaries agree with reference", "status": "RESOLVED",
                 "confidence": "moderate"},
            ],
            events=[("depth_analysis", "gap confirmed", "H1:support"),
                    ("annot_check", "11 errors", "H1:against")])
        self.module.case_report(SimpleNamespace(directory=str(self.case)))
        self.md = (self.case / "case.md").read_text(encoding="utf-8")

    def test_anomalies_are_rendered(self):
        for token in ("A-01", "A-02", "A-03", "49 bp zero-coverage gap", "trnI not detected"):
            self.assertIn(token, self.md, "case-report 必须渲染 anomalies[]（缺 %s）" % token)

    def test_status_and_confidence_are_rendered_per_anomaly(self):
        self.assertIn("UNRESOLVED", self.md)
        self.assertIn("READS_CONSISTENT", self.md)

    def test_four_mandated_sections_are_present(self):
        for heading in ("观察事实", "有证据支持的结论", "无证据支持的声明", "下一步最小实验"):
            self.assertIn(heading, self.md, "缺少报告段落: %s" % heading)

    def test_evidence_matrix_is_present(self):
        self.assertIn("证据矩阵", self.md)

    def test_unresolved_claims_are_listed_as_unsupported(self):
        unsupported = self.md.split("无证据支持的声明", 1)[1]
        self.assertIn("49 bp zero-coverage gap", unsupported)
        self.assertNotIn("gene boundaries agree with reference", unsupported.split("##")[0])

    def test_case_type_is_reported(self):
        self.assertIn("abnormal_case", self.md)

    def test_empty_anomalies_are_stated_explicitly(self):
        bare = self.make_case("case-bare")
        self.module.case_report(SimpleNamespace(directory=str(bare)))
        text = (bare / "case.md").read_text(encoding="utf-8")
        self.assertIn("尚无逐异常判定", text,
                      "没有逐异常记录时必须明说，而不是静默省略")


class LessonDomainTests(_Base):
    """R-3/R-7: a tool_failure_case must not yield a biological lesson."""

    def test_tool_failure_case_yields_a_tool_lesson_and_records_origin(self):
        case = self.make_case("case-tool", case_type="tool_failure_case")
        self.module.propose_lesson(SimpleNamespace(case=str(case), lesson_id="lesson-tool",
                                                  next_test="check format", allow_unresolved=False))
        lesson = self.lesson("lesson-tool")
        self.assertEqual(lesson.get("lesson_domain"), "tool")
        self.assertEqual(lesson.get("source_case_type"), "tool_failure_case")

    def test_tool_failure_case_cannot_produce_a_biology_lesson(self):
        case = self.make_case("case-tool", case_type="tool_failure_case")
        with self.assertRaises(ValueError):
            self.module.propose_lesson(SimpleNamespace(case=str(case), lesson_id="lesson-bio",
                                                       next_test="pileup", allow_unresolved=False,
                                                       lesson_domain="biology"))
        self.assertFalse(self.lesson_exists("lesson-bio"), "被拒绝时不得写出候选经验")

    def test_abnormal_case_can_still_produce_a_biology_lesson(self):
        case = self.make_case("case-abnormal")
        self.module.propose_lesson(SimpleNamespace(case=str(case), lesson_id="lesson-abn",
                                                  next_test="pileup", allow_unresolved=False,
                                                  lesson_domain="biology"))
        lesson = self.lesson("lesson-abn")
        self.assertEqual(lesson["lesson_domain"], "biology")
        self.assertEqual(lesson["source_case_type"], "abnormal_case")

    def test_invalid_lesson_domain_is_rejected(self):
        case = self.make_case("case-abnormal")
        with self.assertRaises(ValueError):
            self.module.propose_lesson(SimpleNamespace(case=str(case), lesson_id="lesson-bad",
                                                       next_test="pileup", allow_unresolved=False,
                                                       lesson_domain="astrology"))
        self.assertFalse(self.lesson_exists("lesson-bad"))

    def test_lesson_schema_declares_domain_and_origin(self):
        schema = json.loads((ROOT / "schemas" / "lesson.schema.json").read_text(encoding="utf-8"))
        self.assertEqual(set(schema["properties"]["lesson_domain"]["enum"]),
                         set(self.module.LESSON_DOMAINS))
        self.assertIn("source_case_type", schema["properties"])


class PublicLessonScopeTests(_Base):
    """R-5: imported/verified lessons must not bypass the fail-closed scope contract."""

    def _validate(self, lesson):
        return self.module.validate_public_lesson(
            json.dumps(lesson).encode("utf-8"), "lesson.json")

    @staticmethod
    def _public(**extra):
        lesson = {"lesson_id": "public-1", "applicable_when": ["stop"],
                  "not_applicable_when": [], "diagnostic_clues": ["stop"],
                  "suggested_next_test": "pileup", "supporting_case_ids": [],
                  "counterexample_case_ids": [], "sources": [], "validation_status": "verified",
                  "version": "1.0.0", "last_reviewed": "2026-09-24"}
        lesson.update(extra)
        return lesson

    def test_missing_scope_fields_are_normalised_to_the_safe_defaults(self):
        lesson = self._validate(self._public())
        self.assertEqual(lesson["generalization_scope"], "single_case")
        self.assertEqual(lesson["transferability"], "none")

    def test_invalid_scope_enum_is_rejected(self):
        with self.assertRaises(ValueError):
            self._validate(self._public(generalization_scope="kingdom"))
        with self.assertRaises(ValueError):
            self._validate(self._public(transferability="certain"))

    def test_single_case_with_transferability_is_rejected(self):
        with self.assertRaises(ValueError):
            self._validate(self._public(generalization_scope="single_case",
                                        transferability="high"))

    def test_claimed_transferability_needs_multiple_supporting_cases(self):
        with self.assertRaises(ValueError):
            self._validate(self._public(generalization_scope="order",
                                        transferability="moderate",
                                        supporting_case_ids=["case-a"]))
        lesson = self._validate(self._public(generalization_scope="family",
                                            transferability="low",
                                            supporting_case_ids=["case-a", "case-b"]))
        self.assertEqual(lesson["transferability"], "low")


class IndependenceCountTests(_Base):
    """R-6: one case_id must count once, whatever taxon values it is given."""

    def test_duplicate_supporting_case_id_is_rejected(self):
        first = self.make_case("case-a", taxon="Species one")
        dup_one = self.make_case("dup-one", taxon="Species one")
        dup_two = self.make_case("dup-two", taxon="Species two")
        # both extra directories claim the same case_id
        for path in (dup_one, dup_two):
            data = json.loads((path / "case.json").read_text(encoding="utf-8"))
            data["case_id"] = "case-a"
            (path / "case.json").write_text(json.dumps(data), encoding="utf-8")
        third = self.make_case("case-b", taxon="Species three")
        with self.assertRaises(ValueError):
            self.module.propose_lesson(SimpleNamespace(
                case=str(first), lesson_id="lesson-dup", next_test="pileup",
                allow_unresolved=False, supporting_case=[str(dup_one), str(dup_two), str(third)],
                generalization_scope="family", transferability="moderate"))
        self.assertFalse(self.lesson_exists("lesson-dup"))

    def test_genuinely_independent_cases_still_work(self):
        first = self.make_case("case-a", taxon="Species one")
        second = self.make_case("case-b", taxon="Species two")
        third = self.make_case("case-c", taxon="Species three")
        self.module.propose_lesson(SimpleNamespace(
            case=str(first), lesson_id="lesson-ok", next_test="pileup", allow_unresolved=False,
            supporting_case=[str(second), str(third)],
            generalization_scope="order", transferability="moderate"))
        self.assertEqual(self.lesson("lesson-ok")["transferability"], "moderate")


class ContributionCaseTypeTests(_Base):
    """R-8: the sharing boundary must not drop the case classification."""

    def test_export_preserves_case_type(self):
        case = self.make_case("case-normal", case_type="normal_validation_case")
        output = self.root / "contribution.json"
        self.module.export_contribution(SimpleNamespace(case=str(case), output=str(output),
                                                        authorize=True))
        payload = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(payload["case"].get("case_type"), "normal_validation_case")


if __name__ == "__main__":
    unittest.main()

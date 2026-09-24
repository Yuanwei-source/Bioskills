import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]


def load_module():
    spec = importlib.util.spec_from_file_location("experience_sharing", ROOT / "tools" / "experience.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


class ExperienceSharingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name); self.module = load_module()
        self.module.KNOWLEDGE = str(self.root / "knowledge")
        self.module.CASES = str(self.root / "knowledge" / "cases")
        self.module.LESSON_CANDIDATES = str(self.root / "knowledge" / "lessons" / "candidates")
        self.module.LESSON_VERIFIED = str(self.root / "knowledge" / "lessons" / "verified")
        self.module.PUBLIC_KNOWLEDGE = str(self.root / "knowledge" / "public")
        self.module.ensure_cases_dir()

    def make_case(self, name, issue="internal_stop", status="RESOLVED"):
        case = self.root / name; case.mkdir()
        data = {"schema_version": "2.0", "case_id": name, "inputs": [],
                "issue": {"type": issue, "user_observation": "same clue"}, "events_file": "events.jsonl",
                "decision": {"status": status, "confidence": "moderate", "rationale": "evidence"}}
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8"); (case / "events.jsonl").touch()
        return case

    def test_contribution_requires_authorization_and_redacts_local_paths(self):
        case = self.make_case("case-a")
        data = json.loads((case / "case.json").read_text())
        data["inputs"] = [{"role": "assembly_fasta", "path": "/home/user/private.fa", "sha256": "a" * 64}]
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8")
        out = self.root / "contribution.json"
        with self.assertRaises(ValueError):
            self.module.export_contribution(SimpleNamespace(case=str(case), output=str(out), authorize=False))
        self.module.export_contribution(SimpleNamespace(case=str(case), output=str(out), authorize=True))
        self.assertNotIn("/home/user", out.read_text(encoding="utf-8"))

    def test_lesson_id_and_events_file_cannot_escape(self):
        case = self.make_case("case-safe")
        with self.assertRaises(ValueError):
            self.module.propose_lesson(SimpleNamespace(case=str(case), lesson_id="../escape", next_test="pileup", allow_unresolved=False))
        data = json.loads((case / "case.json").read_text()); data["events_file"] = "../outside.jsonl"
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaises(ValueError):
            self.module.case_event(SimpleNamespace(directory=str(case), action="x", result="x", impact="", command="", tool_version="", motivation=""))

    def test_conflicting_lessons_are_reported_not_merged(self):
        case = self.make_case("case-a")
        self.module.propose_lesson(SimpleNamespace(case=str(case), lesson_id="lesson-a", next_test="pileup", allow_unresolved=False))
        case2 = self.make_case("case-b", status="NO_CHANGE")
        self.module.propose_lesson(SimpleNamespace(case=str(case2), lesson_id="lesson-b", next_test="annotation", allow_unresolved=False))
        lesson = json.loads((Path(self.module.LESSON_CANDIDATES) / "lesson-b.json").read_text())
        self.assertTrue(lesson["conflicts"])

    def test_lesson_review_promotes_and_records_history(self):
        case = self.make_case("case-review")
        self.module.propose_lesson(SimpleNamespace(case=str(case), lesson_id="review-me", next_test="pileup", allow_unresolved=False))
        self.module.review_lesson(SimpleNamespace(lesson_id="review-me", status="verified", reviewer="human", reason="independent evidence"))
        verified = Path(self.module.LESSON_VERIFIED) / "review-me.json"
        self.assertTrue(verified.exists())
        self.assertEqual(json.loads(verified.read_text())["review_history"][0]["from"], "candidate")

    def test_public_sync_verifies_hash_and_skips_revoked(self):
        public = self.root / "remote"; public.mkdir(); good = public / "lesson.json"
        lesson = {"lesson_id": "public-1", "applicable_when": ["stop"], "not_applicable_when": [], "diagnostic_clues": ["stop"], "suggested_next_test": "pileup", "supporting_case_ids": [], "counterexample_case_ids": [], "sources": [], "validation_status": "verified", "version": "1.0.0", "last_reviewed": "2026-09-24"}
        good.write_text(json.dumps(lesson), encoding="utf-8")
        manifest = self.root / "manifest.json"; manifest.write_text(json.dumps({"format": "mito-public-knowledge-1", "items": [
            {"path": "lesson.json", "url": str(good), "sha256": hashlib.sha256(good.read_bytes()).hexdigest()},
            {"path": "revoked.json", "url": str(good), "sha256": hashlib.sha256(good.read_bytes()).hexdigest(), "status": "revoked"}]}), encoding="utf-8")
        self.module.sync_public(SimpleNamespace(manifest=str(manifest)))
        self.assertTrue((Path(self.module.PUBLIC_KNOWLEDGE) / "lesson.json").exists())
        self.assertFalse((Path(self.module.PUBLIC_KNOWLEDGE) / "revoked.json").exists())
        self.assertTrue(self.module.structured_search("stop"))

    def test_failed_public_update_keeps_previous_cache(self):
        destination = Path(self.module.PUBLIC_KNOWLEDGE); destination.mkdir(parents=True)
        old = destination / "old.json"; old.write_text('{"old":true}', encoding="utf-8")
        manifest = self.root / "bad-manifest.json"
        manifest.write_text(json.dumps({"format": "mito-public-knowledge-1", "items": [{"path": "new.json", "url": str(old), "sha256": "0" * 64}]}), encoding="utf-8")
        with self.assertRaises(ValueError):
            self.module.sync_public(SimpleNamespace(manifest=str(manifest)))
        self.assertEqual(old.read_text(encoding="utf-8"), '{"old":true}')

    def test_offline_search_uses_local_cases(self):
        case = Path(self.module.CASES) / "case-local"; case.mkdir(parents=True)
        data = {"schema_version": "2.0", "case_id": "case-local", "inputs": [],
                "issue": {"type": "rearrangement", "user_observation": "same clue"}, "events_file": "events.jsonl",
                "decision": {"status": "UNRESOLVED", "confidence": "low", "rationale": "offline"}}
        (case / "case.json").write_text(json.dumps(data), encoding="utf-8"); (case / "events.jsonl").touch()
        result = self.module.structured_search("rearrangement")
        self.assertTrue(any(str(case / "case.json") == item[2] for item in result))


if __name__ == "__main__": unittest.main()

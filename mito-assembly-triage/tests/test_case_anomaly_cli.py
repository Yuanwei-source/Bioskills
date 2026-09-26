#!/usr/bin/env python3
"""Failing reproductions for the two experience-system defects from real-data validation.

  F-04  SKILL.md step 6 requires a per-anomaly decision (`anomalies[]` with
        id/claim/status/confidence and optional reads_support, schema-validated by
        `schemas/case.schema.json`), but `tools/experience.py` provides **no command
        that writes that field** -- only validation reads it.  During the real-data
        validation the anomalies therefore had to be inserted with a manual JSON edit,
        which is exactly the "second record path" the skill forbids.
  F-05  `case-init` prepends its own `H<n>` id without stripping one already present
        in `--hypothesis`, so `--hypothesis "H1: ..."` is stored as `id=H1` with an
        explanation starting with `H1:` ("H1 H1: ...").

These tests are CLI-level and use only synthetic cases.
"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gb_fixtures import ROOT  # noqa: E402

EXPERIENCE = ROOT / "tools" / "experience.py"
SCHEMA = ROOT / "schemas" / "case.schema.json"


def _run(*args, timeout=180):
    return subprocess.run([sys.executable, str(EXPERIENCE), *[str(a) for a in args]],
                          capture_output=True, text=True, timeout=timeout)


class CaseAnomalyCliTests(unittest.TestCase):
    """F-04: anomalies[] must be writable through the tool that owns the format."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name) / "case-001"
        result = _run("case-init", self.dir, "--case-id", "case-001",
                      "--issue", "internal_stop", "--observation", "nad5 internal stop",
                      "--taxon", "Lepidoptera",
                      "--hypothesis", "boundary/frame error",
                      "--hypothesis", "base error")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.case_path = self.dir / "case.json"

    def _case(self):
        return json.loads(self.case_path.read_text(encoding="utf-8"))

    def _add(self, *args):
        return _run("case-anomaly", self.dir, *args)

    def test_command_is_known(self):
        result = self._add("--id", "A1", "--claim", "x", "--status", "UNRESOLVED",
                           "--confidence", "low")
        blob = result.stdout + result.stderr
        self.assertNotIn("未知命令", blob,
                         "F-04: experience.py 必须提供写入 anomalies[] 的命令")

    def test_writes_a_schema_valid_anomaly(self):
        result = self._add("--id", "A1", "--claim", "49 bp zero-coverage gap",
                           "--status", "UNRESOLVED", "--confidence", "high",
                           "--reads-support", "READS_CONSISTENT")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        case = self._case()
        self.assertEqual(len(case["anomalies"]), 1)
        entry = case["anomalies"][0]
        self.assertEqual(entry["id"], "A1")
        self.assertEqual(entry["status"], "UNRESOLVED")
        self.assertEqual(entry["confidence"], "high")
        self.assertEqual(entry["reads_support"], "READS_CONSISTENT")
        validate = _run("case-validate", self.dir)
        self.assertEqual(validate.returncode, 0, validate.stdout + validate.stderr)

    def test_appends_without_disturbing_existing_entries(self):
        self._add("--id", "A1", "--claim", "one", "--status", "RESOLVED",
                  "--confidence", "moderate")
        self._add("--id", "A2", "--claim", "two", "--status", "UNRESOLVED",
                  "--confidence", "low")
        case = self._case()
        self.assertEqual([a["id"] for a in case["anomalies"]], ["A1", "A2"])
        # the rest of the case must survive untouched
        self.assertEqual(len(case["hypotheses"]), 2)
        self.assertEqual(case["issue"]["type"], "internal_stop")

    def test_duplicate_id_is_rejected(self):
        self._add("--id", "A1", "--claim", "one", "--status", "RESOLVED",
                  "--confidence", "moderate")
        before = self.case_path.read_text(encoding="utf-8")
        result = self._add("--id", "A1", "--claim", "again", "--status", "RESOLVED",
                           "--confidence", "moderate")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("A1", result.stdout + result.stderr)
        self.assertEqual(self.case_path.read_text(encoding="utf-8"), before,
                         "ID 冲突时不得改动案例文件")

    def test_update_flag_replaces_an_existing_id(self):
        self._add("--id", "A1", "--claim", "one", "--status", "UNRESOLVED",
                  "--confidence", "low")
        result = self._add("--id", "A1", "--claim", "one (revised)", "--status",
                           "RESOLVED", "--confidence", "moderate", "--update")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        case = self._case()
        self.assertEqual(len(case["anomalies"]), 1)
        self.assertEqual(case["anomalies"][0]["status"], "RESOLVED")

    def test_enum_values_are_validated(self):
        cases = {
            'status': ['--id', 'A10', '--claim', 'x', '--status', 'DONE',
                       '--confidence', 'low'],
            'confidence': ['--id', 'A11', '--claim', 'x', '--status', 'UNRESOLVED',
                           '--confidence', 'certain'],
            'reads_support': ['--id', 'A12', '--claim', 'x', '--status', 'UNRESOLVED',
                              '--confidence', 'low', '--reads-support', 'SUPPORTED'],
        }
        for label, args in cases.items():
            with self.subTest(field=label):
                before = self.case_path.read_text(encoding='utf-8')
                result = self._add(*args)
                self.assertNotEqual(result.returncode, 0, '%s 的非法枚举必须被拒绝' % label)
                self.assertEqual(self.case_path.read_text(encoding='utf-8'), before,
                                 '枚举校验失败时不得改动案例文件')

    def test_required_fields_are_enforced(self):
        result = self._add("--id", "A1", "--claim", "x", "--status", "UNRESOLVED")
        self.assertNotEqual(result.returncode, 0, "--confidence 缺失必须失败")

    def test_event_association_must_reference_a_real_event(self):
        result = self._add("--id", "A1", "--claim", "x", "--status", "UNRESOLVED",
                           "--confidence", "low", "--event-action", "depth_analysis")
        self.assertNotEqual(result.returncode, 0,
                            "引用不存在的事件应当失败，否则关联是假的")
        event = _run("case-event", self.dir, "--action", "depth_analysis",
                     "--result", "coverage profile", "--impact", "H1:support")
        self.assertEqual(event.returncode, 0, event.stdout + event.stderr)
        result = self._add("--id", "A1", "--claim", "x", "--status", "UNRESOLVED",
                           "--confidence", "low", "--event-action", "depth_analysis")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self._case()["anomalies"][0]["evidence_events"],
                         ["depth_analysis"])

    def test_no_second_data_format_is_created(self):
        self._add("--id", "A1", "--claim", "one", "--status", "RESOLVED",
                  "--confidence", "moderate")
        names = sorted(p.name for p in self.dir.iterdir())
        self.assertEqual(names, ["case.json", "events.jsonl"],
                         "异常必须写进既有 case.json，不得另立数据文件")

    def test_failed_write_leaves_no_temporary_files(self):
        self._add("--id", "A1", "--claim", "one", "--status", "RESOLVED",
                  "--confidence", "moderate")
        leftovers = [p.name for p in self.dir.iterdir() if p.name.endswith(".tmp")]
        self.assertEqual(leftovers, [], "写入必须是原子的，不留 .tmp")

    def test_case_json_must_exist(self):
        empty = Path(tempfile.mkdtemp()) / "no-case"
        empty.mkdir(parents=True)
        result = _run("case-anomaly", empty, "--id", "A1", "--claim", "x",
                      "--status", "UNRESOLVED", "--confidence", "low")
        self.assertNotEqual(result.returncode, 0)


class HypothesisIdTests(unittest.TestCase):
    """F-05: case-init must not duplicate an id the caller already wrote."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name) / "case-hyp"

    def _init(self, *hypotheses):
        args = ["case-init", self.dir, "--case-id", "case-hyp", "--issue", "internal_stop",
                "--observation", "obs"]
        for text in hypotheses:
            args += ["--hypothesis", text]
        result = _run(*args)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return json.loads((self.dir / "case.json").read_text(encoding="utf-8"))

    def test_prefixed_hypothesis_is_not_duplicated(self):
        case = self._init("H1: boundary/frame error", "H2: base error")
        ids = [h["id"] for h in case["hypotheses"]]
        texts = [h["explanation"] for h in case["hypotheses"]]
        self.assertEqual(ids, ["H1", "H2"])
        for text in texts:
            self.assertFalse(text.startswith("H1"), "不得出现 H1 H1: 之类重复编号")
            self.assertFalse(text.startswith("H2"))

    def test_plain_hypotheses_still_get_generated_ids(self):
        case = self._init("boundary/frame error", "base error")
        self.assertEqual([h["id"] for h in case["hypotheses"]], ["H1", "H2"])
        self.assertEqual(case["hypotheses"][0]["explanation"], "boundary/frame error")

    def test_explicit_ids_are_honoured(self):
        case = self._init("H7: seventh", "H3: third")
        self.assertEqual([h["id"] for h in case["hypotheses"]], ["H7", "H3"])
        self.assertEqual(case["hypotheses"][0]["explanation"], "seventh")
        self.assertEqual(case["hypotheses"][1]["explanation"], "third")

    def test_duplicate_explicit_ids_fail_closed(self):
        args = ["case-init", self.dir, "--case-id", "case-hyp", "--issue",
                "internal_stop", "--observation", "obs",
                "--hypothesis", "H1: a", "--hypothesis", "H1: b"]
        result = _run(*args)
        self.assertNotEqual(result.returncode, 0, "重复的显式假设编号必须失败")
        self.assertFalse((self.dir / "case.json").exists(),
                         "失败时不得留下半成品 case.json")


if __name__ == "__main__":
    unittest.main()

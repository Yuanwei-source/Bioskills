import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "v2_case.py"


class V2CaseTests(unittest.TestCase):
    def test_case_lifecycle_records_hash_events_and_report(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            fasta = root / "assembly.fasta"
            fasta.write_text(">ctg1\nACGTN\n", encoding="utf-8")
            case = root / "case"
            def run(*args):
                return subprocess.run(["python3", str(SCRIPT), *map(str, args)], capture_output=True, text=True, check=False)
            self.assertEqual(run("init", case, "--issue", "internal_stop", "--observation", "nad5 stop", "--input", "assembly_fasta", fasta).returncode, 0)
            self.assertEqual(run("event", case, "--action", "annot_check", "--result", "stop remains", "--impact", "H1:against").returncode, 0)
            self.assertEqual(run("validate", case).stdout.strip(), "VALID")
            self.assertEqual(run("report", case).returncode, 0)
            data = json.loads((case / "case.json").read_text(encoding="utf-8"))
            self.assertEqual(data["schema_version"], "2.0")
            self.assertEqual(len(data["inputs"][0]["sha256"]), 64)
            self.assertIn("annot_check", (case / "case.md").read_text(encoding="utf-8"))

    def test_missing_input_is_rejected_without_case(self):
        with tempfile.TemporaryDirectory() as td:
            case = Path(td) / "case"
            result = subprocess.run(["python3", str(SCRIPT), "init", str(case), "--issue", "x", "--observation", "x", "--input", "assembly_fasta", str(case / "missing.fa")], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((case / "case.json").exists())


if __name__ == "__main__":
    unittest.main()

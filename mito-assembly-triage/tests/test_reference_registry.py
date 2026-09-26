#!/usr/bin/env python3
"""Failing reproductions for the reference-acquisition governance gap.

Reads and assemblies are pinned by SHA-256, but the PUBLIC REFERENCE a diagnosis
compares against had no policy, no registry and no version pinning:

  * nothing separated "download a public reference" from "upload the user's sample";
  * a reference was used and cited as "a close relative" with no accession, date or hash,
    so the comparison cannot be reproduced (and a silent GenBank update changes the answer);
  * L1-L5 reference levels existed, but nothing stopped an L5 (distant) record from being
    used to judge tRNA loss or rearrangement;
  * taxon screening used a local BLAST database whose version was never recorded;
  * a downloaded reference did not become part of the case's evidence chain.

This module defines the registry that fixes that.  All data here is synthetic.
"""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY_SCRIPT = ROOT / "scripts" / "reference_registry.py"
EXPERIENCE = ROOT / "tools" / "experience.py"

SYNTHETIC_GB = """LOCUS       NC_TEST1                  60 bp    DNA     circular INV 01-JAN-2024
DEFINITION  synthetic test mitogenome.
ACCESSION   NC_TEST1
VERSION     NC_TEST1.1
ORGANISM  synthetic organism
FEATURES             Location/Qualifiers
     source          1..60
     CDS             1..1500
                     /gene="cox1"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""


def load_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReferenceRecordTests(unittest.TestCase):
    """A reference record must pin accession.version, source, purpose and integrity."""

    def setUp(self):
        self.module = load_module("reference_registry_unit", Path("scripts") / "reference_registry.py")

    def _base(self, **extra):
        record = {"kind": "sequence", "accession": "NC_060773.1",
                  "organism": "Mimotettix multispinosus",
                  "taxon": "Hemiptera: Cicadellidae", "source": "NCBI GenBank",
                  "retrieved": "2026-09-26", "version": "GenBank 2023-04",
                  "level": "L3", "purposes": ["gene_order_comparison"],
                  "file_sha256": "a" * 64, "sequence_sha256": "b" * 64, "length_bp": 15754}
        record.update(extra)
        return record

    def test_a_complete_record_passes(self):
        self.assertEqual(self.module.validate_record(self._base()), [])

    def test_bare_accession_without_version_is_rejected(self):
        errors = self.module.validate_record(self._base(accession="NC_060773"))
        self.assertTrue(errors, "未固定版本的 accession 必须被拒绝（否则无法复现）")
        self.assertTrue(any("版本" in e for e in errors), errors)

    def test_source_and_purpose_are_required(self):
        self.assertTrue(self.module.validate_record(self._base(source="")))
        self.assertTrue(self.module.validate_record(self._base(purposes=[])))
        self.assertTrue(self.module.validate_record(self._base(purposes=["just_looking"])))

    def test_integrity_fields_are_required_for_sequences(self):
        self.assertTrue(self.module.validate_record(self._base(file_sha256="")))
        self.assertTrue(self.module.validate_record(self._base(file_sha256="zz")))
        self.assertTrue(self.module.validate_record(self._base(sequence_sha256="")))

    def test_level_restricts_what_the_reference_may_be_used_for(self):
        # L5 = 远缘参考：只能做定位/筛查，不得用来判读 tRNA 丢失或重排
        self.assertTrue(self.module.validate_record(
            self._base(level="L5", purposes=["annotation_comparison"])))
        self.assertTrue(self.module.validate_record(
            self._base(level="L5", purposes=["gene_order_comparison"])))
        self.assertEqual(self.module.validate_record(
            self._base(level="L5", purposes=["taxon_screen"])), [])
        # L3 同科：可以看基因排列，但不得用于边界真值
        self.assertTrue(self.module.validate_record(
            self._base(level="L3", purposes=["boundary_validation"])))
        self.assertEqual(self.module.validate_record(
            self._base(level="L3", purposes=["gene_order_comparison"])), [])
        # L1 同种：全部用途
        self.assertEqual(self.module.validate_record(
            self._base(level="L1", purposes=["boundary_validation", "gene_order_comparison"])), [])

    def test_database_records_need_path_and_version_not_accession(self):
        database = {"kind": "database", "name": "MITOS2 refseq invertebrate",
                    "path": "/refs/invertebrate", "version": "refseq89m",
                    "source": "local RefSeq mirror", "retrieved": "2026-09-26",
                    "taxon_filter": "Metazoa", "purposes": ["taxon_screen"],
                    "file_sha256": "c" * 64}
        self.assertEqual(self.module.validate_record(database), [])
        self.assertTrue(self.module.validate_record(dict(database, version="")))
        self.assertTrue(self.module.validate_record(dict(database, path="")))

    def test_unknown_kind_is_rejected(self):
        self.assertTrue(self.module.validate_record(self._base(kind="magic")))


class RegistryStorageTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = Path(temp.name) / "references"
        self.module = load_module("reference_registry_store", Path("scripts") / "reference_registry.py")
        self.gb = Path(temp.name) / "NC_TEST1.1.gb"
        self.gb.write_text(SYNTHETIC_GB, encoding="utf-8")

    def test_register_computes_hashes_and_assigns_an_id(self):
        record = self.module.register_file(self.dir, self.gb, accession="NC_TEST1.1",
                                          source="NCBI GenBank", purposes=["gene_order_comparison"],
                                          organism="synthetic organism", taxon="Testidae",
                                          level="L3", retrieved="2026-09-26")
        self.assertEqual(record["id"], "ref-001")
        self.assertEqual(len(record["file_sha256"]), 64)
        self.assertEqual(record["length_bp"], 60)
        self.assertTrue(record["sequence_sha256"])
        registry = json.loads((self.dir / "registry.json").read_text(encoding="utf-8"))
        self.assertEqual(registry["format"], "mito-reference-registry-1")
        self.assertEqual(len(registry["references"]), 1)

    def test_registering_the_same_file_twice_is_idempotent(self):
        first = self.module.register_file(self.dir, self.gb, accession="NC_TEST1.1",
                                          source="NCBI GenBank",
                                          purposes=["gene_order_comparison"], level="L3",
                                          retrieved="2026-09-26")
        second = self.module.register_file(self.dir, self.gb, accession="NC_TEST1.1",
                                           source="NCBI GenBank",
                                           purposes=["gene_order_comparison"], level="L3",
                                           retrieved="2026-09-26")
        self.assertEqual(first["id"], second["id"])
        registry = json.loads((self.dir / "registry.json").read_text(encoding="utf-8"))
        self.assertEqual(len(registry["references"]), 1, "同一文件不得重复登记")

    def test_same_accession_with_different_content_is_a_drift_error(self):
        self.module.register_file(self.dir, self.gb, accession="NC_TEST1.1",
                                  source="NCBI GenBank", purposes=["gene_order_comparison"],
                                  level="L3", retrieved="2026-09-26")
        changed = Path(self.gb.parent) / "changed.gb"
        changed.write_text(SYNTHETIC_GB.replace("acgtacgtac", "tttttttttt", 1), encoding="utf-8")
        with self.assertRaises(ValueError) as ctx:
            self.module.register_file(self.dir, changed, accession="NC_TEST1.1",
                                      source="NCBI GenBank", purposes=["gene_order_comparison"],
                                      level="L3", retrieved="2026-09-27")
        self.assertIn("hash", str(ctx.exception).lower())

    def test_verify_detects_missing_and_modified_files(self):
        record = self.module.register_file(self.dir, self.gb, accession="NC_TEST1.1",
                                          source="NCBI GenBank",
                                          purposes=["gene_order_comparison"], level="L3",
                                          retrieved="2026-09-26")
        self.assertEqual(self.module.verify_registry(self.dir), [])
        self.gb.write_text(SYNTHETIC_GB + "\n", encoding="utf-8")
        problems = self.module.verify_registry(self.dir)
        self.assertTrue(any("变化" in p or "hash" in p for p in problems), problems)
        self.gb.unlink()
        problems = self.module.verify_registry(self.dir)
        self.assertTrue(any("缺失" in p or "不存在" in p for p in problems), problems)
        self.assertEqual(record["id"], "ref-001")

    def test_registry_never_lives_in_the_skill_install_dir(self):
        record = self.module.register_file(self.dir, self.gb, accession="NC_TEST1.1",
                                          source="NCBI GenBank",
                                          purposes=["gene_order_comparison"], level="L3",
                                          retrieved="2026-09-26")
        self.assertTrue(str(self.dir).startswith(str(self.dir.parent)))
        self.assertNotIn(str(ROOT), json.dumps(record["path"]))
        self.assertFalse((ROOT / "registry.json").exists(),
                         "运行期登记不得写进 skill 安装目录")


class AcquisitionTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.dir = self.root / "references"
        self.module = load_module("reference_registry_acquire", Path("scripts") / "reference_registry.py")

    def test_dry_run_does_not_download_but_reports_the_exact_command(self):
        plan = self.module.plan_acquisition(accession="NC_060773.1", source="ncbi",
                                           purposes=["gene_order_comparison"], outdir=self.dir,
                                           level="L3")
        self.assertIn("efetch", plan["command"])
        self.assertIn("NC_060773.1", plan["command"])
        self.assertFalse((self.dir / "registry.json").exists())

    def test_bare_accession_cannot_be_planned(self):
        with self.assertRaises(ValueError) as ctx:
            self.module.plan_acquisition(accession="NC_060773", source="ncbi",
                                        purposes=["gene_order_comparison"], outdir=self.dir,
                                        level="L3")
        self.assertIn("accession", str(ctx.exception))

    def test_acquisition_uses_the_fetcher_and_registers_the_result(self):
        # a stub fetcher stands in for efetch: no network, same code path
        stub_dir = self.root / "bin"
        stub_dir.mkdir()
        stub = stub_dir / "efetch"
        stub.write_text(textwrap.dedent("""\
            #!/usr/bin/env bash
            # $* holds -db nucleotide -id NC_TEST1.1 -format gbwithparts
            cat "$SYNTHETIC_GB"
            """), encoding="utf-8")
        stub.chmod(0o755)
        synthetic = self.root / "synthetic.gb"
        synthetic.write_text(SYNTHETIC_GB, encoding="utf-8")
        env = dict(os.environ, PATH="%s:%s" % (stub_dir, os.environ.get("PATH", "")),
                   SYNTHETIC_GB=str(synthetic))
        result = subprocess.run(
            [sys.executable, str(REGISTRY_SCRIPT), "acquire", "--accession", "NC_TEST1.1",
             "--purposes", "gene_order_comparison", "--level", "L3",
             "--registry", str(self.dir), "--out", str(self.dir), "--source", "stub"],
            capture_output=True, text=True, timeout=120, env=env)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        registry = json.loads((self.dir / "registry.json").read_text(encoding="utf-8"))
        self.assertEqual(len(registry["references"]), 1)
        entry = registry["references"][0]
        self.assertEqual(entry["accession"], "NC_TEST1.1")
        self.assertEqual(len(entry["file_sha256"]), 64)
        self.assertTrue(Path(entry["path"]).exists())

    def test_truncated_download_is_rejected_and_not_registered(self):
        stub_dir = self.root / "bin2"
        stub_dir.mkdir()
        stub = stub_dir / "efetch"
        stub.write_text("#!/usr/bin/env bash\nprintf 'LOCUS  broken\\n'\n", encoding="utf-8")
        stub.chmod(0o755)
        env = dict(os.environ, PATH="%s:%s" % (stub_dir, os.environ.get("PATH", "")))
        result = subprocess.run(
            [sys.executable, str(REGISTRY_SCRIPT), "acquire", "--accession", "NC_TEST1.1",
             "--purposes", "gene_order_comparison", "--registry", str(self.dir),
             "--out", str(self.dir), "--source", "stub"],
            capture_output=True, text=True, timeout=120, env=env)
        self.assertNotEqual(result.returncode, 0, "截断/非 GenBank 的下载必须失败")
        self.assertFalse((self.dir / "registry.json").exists(), "失败时不得登记")


class CaseReferenceLinkTests(unittest.TestCase):
    """The case must cite the registry entry, not 'a close relative'."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.knowledge = self.root / "knowledge"
        self.module = load_module("experience_ref_link", Path("tools") / "experience.py")
        self.module.KNOWLEDGE = str(self.knowledge)
        self.module.CASES = str(self.knowledge / "cases")
        self.module.REFERENCE_REGISTRY_DIR = str(self.knowledge / "references")
        self.registry = load_module("reference_registry_for_case",
                                    Path("scripts") / "reference_registry.py")
        gb = self.root / "NC_TEST1.1.gb"
        gb.write_text(SYNTHETIC_GB, encoding="utf-8")
        self.record = self.registry.register_file(
            Path(self.module.REFERENCE_REGISTRY_DIR), gb, accession="NC_TEST1.1",
            source="NCBI GenBank", purposes=["gene_order_comparison", "boundary_validation"],
            organism="synthetic organism", taxon="Testidae", level="L2", retrieved="2026-09-26")
        self.case = self.root / "case-001"
        self.module.case_init(type("A", (), {
            "directory": str(self.case), "case_id": "case-001", "case_type": "abnormal_case",
            "issue": "gene_missing", "observation": "nad6 split", "taxon": "Testidae",
            "input": [], "hypothesis": ["H1: annotation error", "H2: assembly error"]})())

    def _case(self):
        return json.loads((self.case / "case.json").read_text(encoding="utf-8"))

    def _link(self, **kwargs):
        args = {"directory": str(self.case), "reference_id": self.record["id"],
                "purpose": "gene_order_comparison", "update": False}
        args.update(kwargs)
        return self.module.case_reference(type("A", (), args)())

    def test_unknown_reference_id_is_rejected(self):
        with self.assertRaises(ValueError):
            self._link(reference_id="ref-999")
        self.assertEqual(self._case().get("references", []), [])

    def test_linking_records_the_pinned_identity(self):
        self._link()
        entry = self._case()["references"][0]
        self.assertEqual(entry["reference_id"], "ref-001")
        self.assertEqual(entry["accession"], "NC_TEST1.1")
        self.assertEqual(entry["file_sha256"], self.record["file_sha256"])
        self.assertEqual(entry["purpose"], "gene_order_comparison")
        validate = self.module.case_validate(type("A", (), {"directory": str(self.case)})())
        self.assertEqual(validate, 0, "关联参考后案例仍须 schema 合法")

    def test_purpose_must_be_declared_on_the_record(self):
        with self.assertRaises(ValueError):
            self._link(purpose="contamination_screen")
        self.assertEqual(self._case().get("references", []), [])

    def test_duplicate_link_is_rejected_unless_update(self):
        self._link()
        with self.assertRaises(ValueError):
            self._link()
        self._link(purpose="boundary_validation")

    def test_report_cites_accession_and_hash(self):
        self._link()
        self.module.case_report(type("A", (), {"directory": str(self.case)})())
        text = (self.case / "case.md").read_text(encoding="utf-8")
        self.assertIn("NC_TEST1.1", text)
        self.assertIn(self.record["file_sha256"][:12], text)
        self.assertIn("参考", text)

    def test_case_schema_declares_the_references_array(self):
        schema = json.loads((ROOT / "schemas" / "case.schema.json").read_text(encoding="utf-8"))
        self.assertIn("references", schema["properties"])
        self.assertEqual(schema["properties"]["references"]["type"], "array")


class ReferenceRegistrySchemaTests(unittest.TestCase):
    def test_schema_declares_the_registry_contract(self):
        schema = json.loads((ROOT / "schemas" / "reference-registry.schema.json")
                            .read_text(encoding="utf-8"))
        self.assertEqual(schema["properties"]["format"]["const"], "mito-reference-registry-1")
        item = schema["properties"]["references"]["items"]
        for key in ("accession", "source", "purposes", "level", "file_sha256"):
            self.assertIn(key, item["properties"])


if __name__ == "__main__":
    unittest.main()

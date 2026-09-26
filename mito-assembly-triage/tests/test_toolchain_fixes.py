#!/usr/bin/env python3
"""Failing reproductions for the three toolchain defects found by real-data validation.

  F-01  `scripts/check_env.sh` reports "核心依赖缺失" (exit 2) when CONDA_ROOT is
        unset, even though every tool exists -- a misleading false negative that
        sends the operator looking for missing software instead of missing
        environment variables.
  F-02  `scripts/run_mitos2.sh` does not create `--outdir`; MITOS2 aborts with
        `no such directory <outdir>`, and `references/tool_check.md` attributes
        that exact message to MITOS2_REFDIR, which is the wrong cause.
  F-03  MITOS2's `runmitos.py` writes no GenBank, `scripts/annot_check.py` accepts
        only GenBank, and the skill ships no converter; MitoFinder needs a
        reference GenBank as input.  So the annotation -> QC loop cannot be closed
        with the shipped tools.

All fixtures here are SYNTHETIC.  The unpublished validation sample is never used:
its data lives outside the repository and is gitignored (`tests/tests_data`).
"""
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gb_fixtures import ROOT  # noqa: E402

SKILL = ROOT
CHECK_ENV = SKILL / "scripts" / "check_env.sh"
RUN_MITOS2 = SKILL / "scripts" / "run_mitos2.sh"
CONVERTER = SKILL / "scripts" / "mitos2_to_genbank.py"
TOOL_CHECK = SKILL / "references" / "tool_check.md"


def _run(cmd, cwd=None, env=None, timeout=180):
    return subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True,
                          timeout=timeout, shell=isinstance(cmd, str))


class CheckEnvHintTests(unittest.TestCase):
    """F-01: an unset CONDA_ROOT must be reported as such, not as missing software."""

    def _env_without_conda(self):
        env = {k: v for k, v in os.environ.items()
               if k not in ('CONDA_ROOT', 'MITOS2_PY', 'MITOS2_REFDIR',
                            'MITOS2_EXTRA_PATH', 'PLOT_PY', 'MINIMAP2', 'PROFILE')}
        return env

    def test_unset_conda_root_is_named_in_the_output(self):
        result = subprocess.run(["bash", str(CHECK_ENV)], capture_output=True, text=True,
                                env=self._env_without_conda(), timeout=180)
        blob = result.stdout + result.stderr
        self.assertIn("CONDA_ROOT", blob,
                      "当 CONDA_ROOT 未设置时，输出必须点名该变量，而不是只说依赖缺失")
        self.assertIn("假阴性", blob,
                      "必须说明此时工具发现可能不完整（假阴性）")

    def test_unset_env_still_exits_nonzero(self):
        result = subprocess.run(["bash", str(CHECK_ENV)], capture_output=True, text=True,
                                env=self._env_without_conda(), timeout=180)
        self.assertNotEqual(result.returncode, 0)

    def test_configured_env_reports_no_false_negative_hint(self):
        env = dict(os.environ)
        env.update({"CONDA_ROOT": "/home/dell/miniforge3"})
        if not Path(env["CONDA_ROOT"]).is_dir():
            self.skipTest("CONDA_ROOT not present on this machine")
        result = subprocess.run(["bash", str(CHECK_ENV)], capture_output=True, text=True,
                                env=env, timeout=300)
        self.assertIn("CONDA_ROOT=/home/dell/miniforge3", result.stdout)
        self.assertNotIn("未设置 CONDA_ROOT", result.stdout)


class RunMitos2OutdirTests(unittest.TestCase):
    """F-02: run_mitos2.sh must create the output directory it is given."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.tmp = Path(temp.name)

    def test_missing_outdir_is_created_before_invoking_mitos2(self):
        env = dict(os.environ)
        # /bin/true accepts any argument list, so this exercises the wrapper logic
        # (argument parsing + mkdir) without needing a real MITOS2 installation.
        env.update({"MITOS2_PY": "/bin/true", "MITOS2_REFDIR": str(self.tmp / "ref"),
                    "MITOS2_REFSEQVER": "refseq89m", "MITOS2_EXTRA_PATH": "",
                    "PROFILE": "test", "CONDA_ROOT": str(self.tmp),
                    # 本用例考的是包装脚本的 arg/mkdir 逻辑（stub 解释器），显式关闭环境门禁
                    "MITO_ENV_GATE": "off"})
        genome = self.tmp / "g.fasta"
        genome.write_text(">g\nACGTACGTACGT\n", encoding="utf-8")
        outdir = self.tmp / "nested" / "mitos2_out"
        result = subprocess.run(
            ["bash", str(RUN_MITOS2), "-i", str(genome), "-o", str(outdir), "-c", "5"],
            capture_output=True, text=True, env=env, timeout=180)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(outdir.is_dir(),
                        "run_mitos2.sh 必须自己创建 --outdir（MITOS2 不会创建）")

    def test_long_outdir_flag_is_also_created(self):
        env = dict(os.environ)
        env.update({"MITOS2_PY": "/bin/true", "MITOS2_REFDIR": str(self.tmp / "ref"),
                    "MITOS2_REFSEQVER": "refseq89m", "MITOS2_EXTRA_PATH": "",
                    "PROFILE": "test", "CONDA_ROOT": str(self.tmp),
                    # 本用例考的是包装脚本的 arg/mkdir 逻辑（stub 解释器），显式关闭环境门禁
                    "MITO_ENV_GATE": "off"})
        genome = self.tmp / "g.fasta"
        genome.write_text(">g\nACGTACGTACGT\n", encoding="utf-8")
        outdir = self.tmp / "long_flag_out"
        result = subprocess.run(
            ["bash", str(RUN_MITOS2), "-i", str(genome), "--outdir", str(outdir), "-c", "5"],
            capture_output=True, text=True, env=env, timeout=180)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(outdir.is_dir())

    def test_missing_mitos2_py_still_fails_closed(self):
        env = dict(os.environ)
        env.update({"MITOS2_PY": "", "MITOS2_REFDIR": str(self.tmp), "MITOS2_REFSEQVER": "x",
                    "MITOS2_EXTRA_PATH": "", "PROFILE": "test", "CONDA_ROOT": str(self.tmp)})
        result = subprocess.run(["bash", str(RUN_MITOS2), "-i", "/dev/null", "-o",
                                 str(self.tmp / "o"), "-c", "5"],
                                capture_output=True, text=True, env=env, timeout=180)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("MITOS2_PY", result.stdout + result.stderr)


class ToolCheckDocTests(unittest.TestCase):
    """F-02: the troubleshooting table must not blame the refdir for an outdir error."""

    def test_no_such_directory_row_mentions_outdir(self):
        text = TOOL_CHECK.read_text(encoding="utf-8")
        rows = [line for line in text.splitlines()
                if "no such directory" in line]
        self.assertTrue(rows, "tool_check.md 应保留该错误的信息行")
        joined = " ".join(rows)
        self.assertIn("outdir", joined.lower(),
                      "`no such directory` 的已知原因必须包含 --outdir（MITOS2 不创建它）")


# --------------------------------------------------------------------------- F-03
SYNTH_GENOME = (
    "ATG" + "AAA" * 5 + "TAA"                      # 1..21   CDS on + strand (7 aa)
    + "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"    # 22..59  spacer
    + "ACGTACGTACGT"                                # 60..71  spacer
    + "ACGTACGTACGTACGTACGTACGTACGTACGT"            # 72..107
)
CDS_START, CDS_END = 1, 21
GENOME_LEN = len(SYNTH_GENOME)
TRNA_START, TRNA_END = 30, 95          # synthetic tRNA span (names carry the evidence)
RRNA_START, RRNA_END = 100, GENOME_LEN


def _write_synthetic_mitos2(outdir, *, genome=None, gff_rows=None, fas_rows=None,
                            faa_rows=None, cds_start=None, cds_end=None,
                            protein_frame=0):
    """Build a minimal MITOS2-like output set (synthetic; no real sample data).

    ``protein_frame`` lets a test declare a CDS whose annotated span starts before
    the real start codon, so the correct ``/codon_start`` is 2 or 3 rather than 1.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    genome = genome if genome is not None else SYNTH_GENOME
    seqid = "synth"
    cds_start = CDS_START if cds_start is None else cds_start
    cds_end = CDS_END if cds_end is None else cds_end
    (outdir / "sequence.fas").write_text(">%s\n%s\n" % (seqid, genome), encoding="utf-8")

    from Bio.Seq import Seq
    cds_nt = genome[cds_start - 1:cds_end]
    cds_aa = str(Seq(cds_nt[protein_frame:]).translate(table=5))

    rows = gff_rows if gff_rows is not None else [
        ["seqid", "mitos", "region", "1", str(len(genome)), ".", "+", ".",
         "ID=%s:1..%d;Is_circular=False;Name=%s" % (seqid, len(genome), seqid)],
        [seqid, "mitos", "gene", str(cds_start), str(cds_end), ".", "+", ".",
         "ID=gene_nad2;Name=nad2;gene_id=nad2"],
        [seqid, "mitos", "exon", str(cds_start), str(cds_end), "1.0", "+", "0",
         "Parent=transcript_nad2;Name=nad2"],
        [seqid, "mitfi", "ncRNA_gene", str(TRNA_START), str(TRNA_END), ".", "+", ".",
         "ID=gene_trnW;Name=trnW;gene_id=trnW"],
        [seqid, "mitfi", "tRNA", str(TRNA_START), str(TRNA_END), ".", "+", ".",
         "ID=transcript_trnW(tca);Name=trnW(tca);Parent=gene_trnW(tca);gene_id=trnW(tca)"],
        [seqid, "mitos", "ncRNA_gene", str(RRNA_START), str(RRNA_END), ".", "-", ".",
         "ID=gene_rrnS;Name=rrnS;gene_id=rrnS"],
        [seqid, "mitos", "rRNA", str(RRNA_START), str(RRNA_END), ".", "-", ".",
         "ID=transcript_rrnS;Name=rrnS;Parent=gene_rrnS;gene_id=rrnS"],
    ]
    if gff_rows is None:
        (outdir / "result.gff").write_text(
            "\n".join("\t".join(r) for r in rows) + "\n", encoding="utf-8")
    else:
        (outdir / "result.gff").write_text(
            "\n".join("\t".join(r) for r in gff_rows) + "\n", encoding="utf-8")

    def extract(start, end, strand):
        from Bio.Seq import Seq
        nt = genome[start - 1:end]
        return str(Seq(nt).reverse_complement()) if strand == "-" else nt

    fas = fas_rows if fas_rows is not None else [
        ("nad2", cds_start, cds_end, "+", extract(cds_start, cds_end, "+")),
        ("trnW(tca)", TRNA_START, TRNA_END, "+", extract(TRNA_START, TRNA_END, "+")),
        ("rrnS", RRNA_START, RRNA_END, "-", extract(RRNA_START, RRNA_END, "-")),
    ]
    (outdir / "result.fas").write_text(
        "".join(">%s; %d-%d; %s; %s\n%s\n" % (seqid, s, e, st, name, nt)
                for name, s, e, st, nt in fas), encoding="utf-8")

    faa = faa_rows if faa_rows is not None else [("nad2", cds_start, cds_end, "+", cds_aa)]
    (outdir / "result.faa").write_text(
        "".join(">%s; %d-%d; %s; %s\n%s\n" % (seqid, s, e, st, name, aa)
                for name, s, e, st, aa in faa), encoding="utf-8")
    return outdir


@unittest.skipUnless(shutil.which("python3"), "python3 required")
class Mitos2ToGenbankTests(unittest.TestCase):
    """F-03: a MITOS2 -> GenBank conversion that adds no annotation of its own."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.tmp = Path(temp.name)
        self.out = _write_synthetic_mitos2(self.tmp / "mitos2")

    def _convert(self, *extra, outdir=None, expect=0, source=None):
        outdir = outdir or (self.tmp / "converted")
        outdir.mkdir(parents=True, exist_ok=True)
        src = source or self.out
        gb = outdir / "out.gb"
        result = subprocess.run(
            ["python3", str(CONVERTER),
             str(src / "result.gff"), str(src / "result.fas"),
             str(src / "result.faa"), str(gb)] + list(extra),
            capture_output=True, text=True, timeout=300)
        if expect == 0:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result, gb

    def test_script_exists_and_is_executable_by_python(self):
        self.assertTrue(CONVERTER.is_file(), "F-03 需要 scripts/mitos2_to_genbank.py")

    def test_features_and_names_are_preserved(self):
        _, gb = self._convert()
        from Bio import SeqIO
        record = SeqIO.read(gb, "genbank")
        by_type = {}
        for feature in record.features:
            by_type.setdefault(feature.type, []).append(feature)
        self.assertEqual(len(record.seq), len(SYNTH_GENOME))
        self.assertEqual(len(by_type.get("CDS", [])), 1)
        self.assertEqual(by_type["CDS"][0].qualifiers["gene"][0], "nad2")
        self.assertEqual([f.qualifiers["gene"][0] for f in by_type.get("tRNA", [])],
                         ["trnW(tca)"])
        self.assertEqual(by_type["CDS"][0].qualifiers["transl_table"][0], "5")
        # coordinates and strand come from MITOS2, 1-based inclusive in GFF
        self.assertEqual(int(by_type["CDS"][0].location.start) + 1, CDS_START)
        self.assertEqual(int(by_type["CDS"][0].location.end), CDS_END)
        self.assertEqual(by_type["rRNA"][0].location.strand, -1)

    def test_per_feature_sequence_matches_mitos2_extractions(self):
        _, gb = self._convert()
        from Bio import SeqIO
        from Bio.Seq import Seq
        record = SeqIO.read(gb, "genbank")
        for feature in record.features:
            if feature.type not in ("CDS", "tRNA", "rRNA"):
                continue
            extracted = str(feature.extract(record.seq))
            name = feature.qualifiers["gene"][0]
            mitos = None
            for line in (self.out / "result.fas").read_text(encoding="utf-8").splitlines():
                if line.startswith(">"):
                    current = line.split(";")[-1].strip()
                elif current == name:
                    mitos = (mitos or "") + line.strip()
            self.assertIsNotNone(mitos, "MITOS2 result.fas 应含 %s" % name)
            self.assertEqual(extracted.upper(), mitos.upper(),
                             "%s 的序列必须与 MITOS2 自身提取一致" % name)

    def test_codon_start_is_recovered_from_mitos2_protein(self):
        # The annotated span starts one base before the real ATG, so only frame 2
        # reproduces MITOS2's own protein.  /codon_start must therefore be 2 -- it
        # must never be silently assumed to be 1.
        shifted = "A" + SYNTH_GENOME
        out = _write_synthetic_mitos2(self.tmp / "shifted", genome=shifted,
                                      cds_start=1, cds_end=len(SYNTH_GENOME) + 1,
                                      protein_frame=1)
        _, gb = self._convert(outdir=self.tmp / "converted_shifted", source=out)
        from Bio import SeqIO
        record = SeqIO.read(gb, "genbank")
        cds = [f for f in record.features if f.type == "CDS"][0]
        self.assertEqual(cds.qualifiers["codon_start"][0], "2",
                         "必须用 MITOS2 自己的蛋白恢复读框，而不是默认 codon_start=1")

    def test_ambiguous_frame_is_reported_rather_than_guessed(self):
        # No protein result at all -> the frame cannot be established.  The record
        # may fall back to 1, but that uncertainty must be stated in the output.
        out = _write_synthetic_mitos2(self.tmp / "nofaa", faa_rows=[])
        result = subprocess.run(
            ["python3", str(CONVERTER), str(out / "result.gff"), str(out / "result.fas"),
             str(out / "result.faa"), str(self.tmp / "nofaa.gb")],
            capture_output=True, text=True, timeout=180)
        blob = result.stdout + result.stderr
        self.assertIn("codon_start", blob.lower())
        self.assertIn("不确定", blob)

    def test_genbank_round_trip_is_stable(self):
        _, gb = self._convert()
        from Bio import SeqIO
        first = SeqIO.read(gb, "genbank")
        again = self.tmp / "roundtrip.gb"
        SeqIO.write(first, again, "genbank")
        second = SeqIO.read(again, "genbank")
        self.assertEqual(str(first.seq), str(second.seq))
        self.assertEqual([(f.type, int(f.location.start), int(f.location.end),
                           f.location.strand) for f in first.features],
                         [(f.type, int(f.location.start), int(f.location.end),
                           f.location.strand) for f in second.features])

    def test_topology_is_never_inferred_from_the_mitos2_run_mode(self):
        result, gb = self._convert()
        from Bio import SeqIO
        record = SeqIO.read(gb, "genbank")
        self.assertEqual(record.annotations.get("topology", "linear"), "linear",
                         "MITOS2 的 circular 模式不等于物理环化证据；默认必须写 linear")
        self.assertIn("linear", (result.stdout + result.stderr).lower())

    def test_output_declares_it_is_qc_input_not_submission_quality(self):
        result, _ = self._convert()
        blob = (result.stdout + result.stderr)
        self.assertTrue("质检" in blob or "QC" in blob,
                        "输出必须声明它是注释质检输入，不是 NCBI 提交质量记录")
        self.assertIn("来源", blob)

    def test_missing_inputs_fail_closed(self):
        outdir = self.tmp / "fail"
        outdir.mkdir()
        result = subprocess.run(
            ["python3", str(CONVERTER), str(self.tmp / "missing.gff"),
             str(self.out / "result.fas"), str(self.out / "result.faa"),
             str(outdir / "o.gb")], capture_output=True, text=True, timeout=180)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((outdir / "o.gb").exists(),
                         "失败时不得留下半成品输出")

    def test_coordinates_outside_the_genome_fail_closed(self):
        outdir = _write_synthetic_mitos2(
            self.tmp / "badcoord",
            gff_rows=[
                ["synth", "mitos", "region", "1", str(GENOME_LEN), ".", "+", ".", "Name=synth"],
                ["synth", "mitos", "gene", "1", "99999", ".", "+", ".",
                 "ID=gene_nad2;Name=nad2;gene_id=nad2"],
            ])
        result = subprocess.run(
            ["python3", str(CONVERTER), str(outdir / "result.gff"),
             str(outdir / "result.fas"), str(outdir / "result.faa"),
             str(self.tmp / "bad.gb")], capture_output=True, text=True, timeout=180)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.tmp / "bad.gb").exists())

    def test_truncated_fas_sequence_fails_closed(self):
        outdir = _write_synthetic_mitos2(
            self.tmp / "badseq",
            fas_rows=[("nad2", CDS_START, CDS_END, "+", "ATG")])
        result = subprocess.run(
            ["python3", str(CONVERTER), str(outdir / "result.gff"),
             str(outdir / "result.fas"), str(outdir / "result.faa"),
             str(self.tmp / "badseq.gb")], capture_output=True, text=True, timeout=180)
        self.assertNotEqual(result.returncode, 0,
                            "result.fas 与坐标不一致必须失败，而不是静默产出错误的 GenBank")

    def test_unknown_feature_type_is_reported_not_silently_dropped(self):
        outdir = _write_synthetic_mitos2(
            self.tmp / "unknown",
            gff_rows=[
                ["synth", "mitos", "region", "1", str(GENOME_LEN), ".", "+", ".", "Name=synth"],
                ["synth", "mitos", "mystery_orf", "1", "21", ".", "+", ".",
                 "ID=gene_x;Name=x;gene_id=x"],
            ])
        result = subprocess.run(
            ["python3", str(CONVERTER), str(outdir / "result.gff"),
             str(outdir / "result.fas"), str(outdir / "result.faa"),
             str(self.tmp / "unknown.gb")], capture_output=True, text=True, timeout=180)
        blob = result.stdout + result.stderr
        self.assertIn("mystery_orf", blob,
                      "未识别的 feature 类型必须被点名报告，而不是静默丢弃")


if __name__ == "__main__":
    unittest.main()

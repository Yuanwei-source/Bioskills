import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class FastqQualityTests(unittest.TestCase):
    def test_stdout_reads_are_discarded_and_full_read_count_is_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            r1, r2 = root/'R1.fastq.gz', root/'R2.fastq.gz'
            r1.write_bytes(b'R1 raw')
            r2.write_bytes(b'R2 raw')
            manifest = root/'assembly_candidates.json'
            manifest.write_text(json.dumps({
                'format': 'mito-illumina-assembly-candidates-1', 'sample': 'T1', 'taxon': 'Test species',
                'read_count': 3,
                'inputs': [{'role': 'R1', 'path': str(r1), 'sha256': hashlib.sha256(r1.read_bytes()).hexdigest()},
                           {'role': 'R2', 'path': str(r2), 'sha256': hashlib.sha256(r2.read_bytes()).hexdigest()}],
            }))
            fastp = root/'fastp'
            fastp.write_text('#!/usr/bin/env python3\n'
                             'import json,sys\n'
                             'if "--version" in sys.argv: print("fastp 1.0.0"); raise SystemExit(0)\n'
                             'args=sys.argv\n'
                             'j=args[args.index("--json")+1]\n'
                             'h=args[args.index("--html")+1]\n'
                             'json.dump({"summary":{"before_filtering":{"total_reads":6}}, '
                             '"read1_before_filtering":{"total_reads":3,"q30_rate":0.99}, '
                             '"read2_before_filtering":{"total_reads":3,"q30_rate":0.98}},open(j,"w"))\n'
                             'open(h,"w").write("report")\n'
                             'print("FASTQ_OUTPUT_SHOULD_BE_DISCARDED")\n')
            fastp.chmod(0o755)
            outdir = root/'qc'
            result = subprocess.run([sys.executable, str(ROOT/'scripts/fastq_qc.py'),
                                     '--assembly-manifest', str(manifest), '--outdir', str(outdir),
                                     '--fastp', str(fastp)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn('FASTQ_OUTPUT_SHOULD_BE_DISCARDED', result.stdout)
            report = json.loads((outdir/'fastq_quality.json').read_text())
            self.assertEqual(report['status'], 'succeeded')
            self.assertEqual(report['fastp']['observed_read1_count'], 3)
            self.assertTrue(report['fastp']['read_count_matches_pair_validation'])

    def test_mismatched_fastp_count_is_not_reported_as_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            r1, r2 = root/'R1', root/'R2'
            r1.write_bytes(b'R1')
            r2.write_bytes(b'R2')
            manifest = root/'m.json'
            manifest.write_text(json.dumps({'format': 'mito-illumina-assembly-candidates-1', 'read_count': 10,
                'inputs': [{'role': 'R1', 'path': str(r1), 'sha256': hashlib.sha256(b'R1').hexdigest()},
                           {'role': 'R2', 'path': str(r2), 'sha256': hashlib.sha256(b'R2').hexdigest()}]}))
            fastp = root/'fastp'
            fastp.write_text('#!/usr/bin/env python3\nimport json,sys\n'
                             'if "--version" in sys.argv: print("1.0"); raise SystemExit(0)\n'
                             'j=sys.argv[sys.argv.index("--json")+1]; h=sys.argv[sys.argv.index("--html")+1]\n'
                             'json.dump({"read1_before_filtering":{"total_reads":2}},open(j,"w")); open(h,"w").close()\n')
            fastp.chmod(0o755)
            outdir = root/'qc'
            result = subprocess.run([sys.executable, str(ROOT/'scripts/fastq_qc.py'),
                                     '--assembly-manifest', str(manifest), '--outdir', str(outdir),
                                     '--fastp', str(fastp)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            report = json.loads((outdir/'fastq_quality.json').read_text())
            self.assertEqual(report['status'], 'failed_read_count_mismatch')


if __name__ == '__main__':
    unittest.main()

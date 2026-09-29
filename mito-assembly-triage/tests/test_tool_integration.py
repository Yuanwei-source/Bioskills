"""Real samtools/BLAST integration on generated public-free fixtures."""
import json
import gzip
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
import zlib

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from scripts.run_illumina_candidates import normalize_gzip_members, run_one

ROOT = Path(__file__).resolve().parents[1]


class ToolIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def bam(self, length=100, covered=100, count=10, chrom='mt', missing_quality=False, mapq=60):
        if not shutil.which('samtools'):
            self.skipTest('samtools not installed')
        sam = self.root/'input.sam'
        lines = ['@HD\tVN:1.6\tSO:coordinate', '@SQ\tSN:%s\tLN:%d' % (chrom, length)]
        for i in range(count):
            lines.append('\t'.join([str(i), str(16*(i%2)), chrom, '1', str(mapq), str(covered)+'M',
                                    '*', '0', '0', 'A'*covered, '*' if missing_quality else 'I'*covered]))
        sam.write_text('\n'.join(lines)+'\n')
        bam = self.root/'test.bam'
        subprocess.run(['samtools', 'sort', '-o', str(bam), str(sam)], check=True, capture_output=True)
        subprocess.run(['samtools', 'index', str(bam)], check=True, capture_output=True)
        return bam

    def run_depth(self, bam, length=100, chrom='mt', extra=()):
        fasta = self.root/'genome.fa'
        fasta.write_text('>'+chrom+'\n'+'A'*length+'\n')
        return subprocess.run([sys.executable, str(ROOT/'scripts/depth_analysis.py'), str(bam),
                               str(fasta), '--output-json', str(self.root/'depth.json'), *extra],
                              capture_output=True, text=True)

    def test_real_zero_coverage_returns_review(self):
        result = self.run_depth(self.bam(count=0))
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        report = json.loads((self.root/'depth.json').read_text())
        self.assertEqual(report['global_mean'], 0)
        self.assertEqual(report['zero_coverage_bases'], 100)

    def test_illumina_candidate_orchestrator_keeps_tool_outputs_and_does_not_choose(self):
        r1, r2 = self.root/'sample_R1.fastq', self.root/'sample_R2.fastq'
        record = '@read/1\nACGT\n+\nIIII\n'
        r1.write_text(record)
        r2.write_text(record.replace('/1', '/2'))
        outputs = {}
        for name, content in (
            ('getorganelle', '#!/bin/sh\necho "fake getorganelle 1.0"\n'
            'out=\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
            '[ ! -e "$out" ] || { echo "output directory already exists"; exit 9; }\n'
            'mkdir -p "$out"\nprintf ">mt_go\\nACGTACGT\\n" > "$out/sample.path_sequence.fasta"\n'),
            ('mitoflex', '#!/bin/sh\necho "fake mitoflex 1.0"\n'
             'mkdir -p run/sample/sample.temp/findmitoscaf\n'
             'printf ">mt_mf\\nACGTACGT\\n" > run/sample/sample.temp/findmitoscaf/sample.picked.fa\n')):
            exe = self.root/name
            exe.write_text(content)
            exe.chmod(0o755)
            outputs[name] = exe
        outdir = self.root/'sample_mitogenome'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
            '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
            '--sample', 'sample', '--taxon', 'Insect species',
            '--tools', 'getorganelle,mitoflex',
            '--getorganelle', str(outputs['getorganelle']), '--mitoflex', str(outputs['mitoflex']),
        ], capture_output=True, text=True, env=dict(os.environ, CONDA_ROOT=str(self.root/'empty-conda-root')))
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        manifest = json.loads((outdir/'intermediate/assemblies/assembly_candidates.json').read_text())
        self.assertEqual(manifest['candidate_count'], 2)
        self.assertEqual({item['producer'] for item in manifest['candidates']}, {'getorganelle', 'mitoflex'})
        self.assertEqual(manifest['assembly_decision'], 'not_selected')
        self.assertEqual(r1.read_text(), record)
        self.assertEqual(r2.read_text(), record.replace('/1', '/2'))
        flex_command = manifest['tools'][1]['command']
        self.assertIn('--disable-annotation', flex_command)
        self.assertNotIn('--disable-visualization', flex_command)
        self.assertIn('--threads', flex_command)
        self.assertFalse((outdir/'sample_assembly.fasta').exists())

    def test_illumina_candidate_orchestrator_protects_existing_output_dir(self):
        r1, r2 = self.root/'r1.fastq', self.root/'r2.fastq'
        r1.write_text('@r1\nA\n+\nI\n')
        r2.write_text('@r2\nA\n+\nI\n')
        outdir = self.root/'existing'
        outdir.mkdir()
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
            '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
            '--sample', 'sample', '--taxon', 'Insect species', '--tools', 'getorganelle',
        ], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn('输出目录必须不存在', result.stderr)
        self.assertEqual(list(outdir.iterdir()), [])

    def test_background_task_manager_runs_multitool_assembly_in_sample_intermediate(self):
        r1, r2 = self.root/'r1.fastq', self.root/'r2.fastq'
        r1.write_text('@read/1\nACGT\n+\nIIII\n')
        r2.write_text('@read/2\nACGT\n+\nIIII\n')
        getorganelle = self.root/'getorganelle'
        getorganelle.write_text(
            '#!/bin/sh\n'
            'if [ "$1" = "--version" ]; then echo "GetOrganelle 1.0"; exit 0; fi\n'
            'out=\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
            '[ ! -e "$out" ] || { echo "output directory already exists"; exit 9; }\n'
            'mkdir -p "$out"\nprintf ">mt_go\\nACGTACGT\\n" > "$out/sample.path_sequence.fasta"\n')
        getorganelle.chmod(0o755)
        mitoz = self.root/'mitoz'
        mitoz.write_text(
            '#!/bin/sh\n'
            'workdir=\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "--workdir" ]; then workdir="$2"; shift 2; else shift; fi; done\n'
            'mkdir -p "$workdir/backup"\nprintf ">mt_mz\\nACGTACGT\\n" > "$workdir/run.mitogenome.fa"\n'
            'cp "$workdir/run.mitogenome.fa" "$workdir/backup/copy.mitogenome.fa"\n')
        mitoz.chmod(0o755)
        outdir = self.root/'CMMC_mitogenome'
        launch = subprocess.run([
            sys.executable, str(ROOT/'scripts/task_manager.py'), 'start',
            '--task-root', str(outdir/'intermediate'/'tasks'), '--name', 'real-workflow',
            '--cwd', str(ROOT), '--', sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
            '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
            '--sample', 'sample', '--taxon', 'Test_species', '--tools', 'getorganelle,mitoz',
            '--table', '5', '--getorganelle', str(getorganelle), '--mitoz', str(mitoz), '--threads', '1',
        ], capture_output=True, text=True, check=True,
           env=dict(os.environ, CONDA_ROOT=str(self.root/'empty-conda-root')))
        task_dir = Path(json.loads(launch.stdout)['task_dir'])
        deadline = time.monotonic()+10
        while time.monotonic() < deadline:
            status = subprocess.run([
                sys.executable, str(ROOT/'scripts/task_manager.py'), 'status', '--task-dir', str(task_dir),
            ], capture_output=True, text=True, check=True)
            state = json.loads(status.stdout)
            if state['status'] not in ('queued', 'running'):
                break
            time.sleep(.05)
        self.assertEqual(state['status'], 'succeeded', state)
        self.assertEqual(state['exit_code'], 0)
        manifest = json.loads((outdir/'intermediate'/'assemblies'/'assembly_candidates.json').read_text())
        self.assertEqual(manifest['candidate_count'], 2)
        self.assertTrue((task_dir/'task.log').is_file())

    def test_illumina_assemblers_start_concurrently(self):
        r1, r2 = self.root/'r1.fastq', self.root/'r2.fastq'
        r1.write_text('@read/1\nACGT\n+\nIIII\n')
        r2.write_text('@read/2\nACGT\n+\nIIII\n')
        marker1, marker2 = self.root/'getorganelle.started', self.root/'mitoz.started'
        getorganelle = self.root/'getorganelle'
        getorganelle.write_text(
            '#!/bin/sh\n'
            'if [ "$1" = "--version" ]; then echo "GetOrganelle 1.0"; exit 0; fi\n'
            'out=\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
            'touch "%s"\n'
            'i=0; while [ ! -f "%s" ] && [ "$i" -lt 100 ]; do sleep 0.05; i=$((i+1)); done\n'
            '[ -f "%s" ] || { echo "MitoZ was not launched concurrently"; exit 9; }\n'
            'mkdir -p "$out"\nprintf \">mt_go\\nACGTACGT\\n" > "$out/sample.path_sequence.fasta"\n' %
            (marker1, marker2, marker2))
        getorganelle.chmod(0o755)
        mitoz = self.root/'mitoz'
        mitoz.write_text(
            '#!/bin/sh\n'
            'workdir=\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "--workdir" ]; then workdir="$2"; shift 2; else shift; fi; done\n'
            'touch "%s"\n'
            'i=0; while [ ! -f "%s" ] && [ "$i" -lt 100 ]; do sleep 0.05; i=$((i+1)); done\n'
            '[ -f "%s" ] || { echo "GetOrganelle was not launched concurrently"; exit 9; }\n'
            'mkdir -p "$workdir"\nprintf \">mt_mz\\nTTTTACGT\\n" > "$workdir/run.mitogenome.fa"\n' %
            (marker2, marker1, marker1))
        mitoz.chmod(0o755)
        outdir = self.root/'parallel_sample'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
            '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
            '--sample', 'sample', '--taxon', 'Test_species', '--tools', 'getorganelle,mitoz',
            '--getorganelle', str(getorganelle), '--mitoz', str(mitoz), '--table', '5',
            '--parallel-tools', '2', '--threads', '1',
        ], capture_output=True, text=True, env=dict(os.environ, CONDA_ROOT=str(self.root/'empty-conda-root')))
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        manifest = json.loads((outdir/'intermediate/assemblies/assembly_candidates.json').read_text())
        self.assertEqual([item['status'] for item in manifest['tools']], ['succeeded', 'succeeded'])
        self.assertEqual({item['producer'] for item in manifest['candidates']}, {'getorganelle', 'mitoz'})

    def test_mitoflex_python_site_supplement_is_inherited_by_child_process(self):
        package_root = self.root/'supplemental-site-packages'
        package = package_root/'compat_probe'
        package.mkdir(parents=True)
        (package/'__init__.py').write_text('VALUE = "supplemented"\n')
        workdir = self.root/'tool-work'
        workdir.mkdir()
        result = run_one(
            'mitoflex-probe', [sys.executable, '-c'],
            ['import compat_probe; print(compat_probe.VALUE)'], workdir,
            python_site_paths=[str(package_root)],
        )
        self.assertEqual(result['status'], 'succeeded')
        self.assertIn('supplemented', (workdir/'tool.log').read_text())

    def test_mitoflex_input_normalizer_flattens_concatenated_gzip_members(self):
        source = self.root/'lanes.fastq.gz'
        part1 = b'@read1/1\nACGT\n+\nIIII\n'
        part2 = b'@read2/1\nTGCA\n+\nIIII\n'
        source.write_bytes(gzip.compress(part1, mtime=0)+gzip.compress(part2, mtime=0))
        destination = self.root/'single.fastq.gz'
        result = normalize_gzip_members(source, destination)
        first_member = zlib.decompressobj(31)
        content = first_member.decompress(destination.read_bytes())
        self.assertTrue(first_member.eof)
        self.assertEqual(first_member.unused_data, b'')
        self.assertEqual(content, part1+part2)
        self.assertEqual(result['source'], str(source.resolve()))
        self.assertEqual(result['sha256'], __import__('hashlib').sha256(destination.read_bytes()).hexdigest())

    def test_illumina_assembly_resume_retries_failed_tool_without_overwriting_attempt(self):
        r1, r2 = self.root/'r1.fastq', self.root/'r2.fastq'
        r1.write_text('@read/1\nACGT\n+\nIIII\n')
        r2.write_text('@read/2\nACGT\n+\nIIII\n')
        marker = self.root/'failed-once'
        exe = self.root/'getorganelle'
        exe.write_text(
            '#!/bin/sh\n'
            'if [ "$1" = "--version" ]; then echo "GetOrganelle 1.0"; exit 0; fi\n'
            'if [ ! -f "%s" ]; then touch "%s"; echo "transient test failure"; exit 9; fi\n'
            'out=\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
            '[ ! -e "$out" ] || exit 8\nmkdir -p "$out"\nprintf ">mt\\nACGTACGT\\n" > "$out/mt.path_sequence.fasta"\n' %
            (marker, marker))
        exe.chmod(0o755)
        outdir = self.root/'resume_sample'
        common = [sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
                  '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
                  '--sample', 'sample', '--taxon', 'Test_species', '--tools', 'getorganelle',
                  '--getorganelle', str(exe), '--threads', '1']
        first = subprocess.run(common, capture_output=True, text=True)
        self.assertEqual(first.returncode, 2, first.stdout+first.stderr)
        root = outdir/'intermediate'/'assemblies'
        old_log = root/'getorganelle'/'tool.log'
        self.assertIn('transient test failure', old_log.read_text())
        second = subprocess.run(common+['--resume'], capture_output=True, text=True)
        self.assertEqual(second.returncode, 2, second.stdout+second.stderr)
        manifest = json.loads((root/'assembly_candidates.json').read_text())
        self.assertEqual([run['status'] for run in manifest['tools']], ['failed', 'succeeded'])
        self.assertTrue(old_log.is_file())
        self.assertEqual(manifest['candidate_count'], 1)
        third = subprocess.run(common+['--resume'], capture_output=True, text=True)
        self.assertEqual(third.returncode, 1)
        self.assertIn('不会覆盖已成功工具', third.stderr)

    def test_getorganelle_zero_exit_without_candidate_is_failed_and_resumable(self):
        r1, r2 = self.root/'r1.fastq', self.root/'r2.fastq'
        r1.write_text('@read/1\nACGT\n+\nIIII\n')
        r2.write_text('@read/2\nACGT\n+\nIIII\n')
        marker = self.root/'silent-error-once'
        exe = self.root/'getorganelle'
        exe.write_text(
            '#!/bin/sh\n'
            'if [ "$1" = "--version" ]; then echo "GetOrganelle 1.0"; exit 0; fi\n'
            'out=\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
            'if [ ! -f "%s" ]; then touch "%s"; echo "output directory existed"; exit 0; fi\n'
            'mkdir -p "$out"\nprintf \">mt\\nACGTACGT\\n" > "$out/mt.path_sequence.fasta"\n' % (marker, marker))
        exe.chmod(0o755)
        outdir = self.root/'silent_error_sample'
        common = [sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
                  '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
                  '--sample', 'sample', '--taxon', 'Test_species', '--tools', 'getorganelle',
                  '--getorganelle', str(exe), '--threads', '1']
        first = subprocess.run(common, capture_output=True, text=True)
        self.assertEqual(first.returncode, 2, first.stdout+first.stderr)
        root = outdir/'intermediate'/'assemblies'
        manifest = json.loads((root/'assembly_candidates.json').read_text())
        self.assertEqual(manifest['tools'][0]['status'], 'failed')
        self.assertIn('postcondition_error', manifest['tools'][0])
        second = subprocess.run(common+['--resume'], capture_output=True, text=True)
        self.assertEqual(second.returncode, 2, second.stdout+second.stderr)
        manifest = json.loads((root/'assembly_candidates.json').read_text())
        self.assertEqual([run['status'] for run in manifest['tools']], ['failed', 'succeeded'])
        self.assertEqual(manifest['candidate_count'], 1)

    def test_illumina_candidate_orchestrator_rejects_mispaird_or_truncated_fastq(self):
        r1, r2 = self.root/'r1.fastq', self.root/'r2.fastq'
        r1.write_text('@r1/1\nA\n+\nI\n@r2/1\nC\n+\nI\n')
        r2.write_text('@r1/2\nA\n+\nI\n@wrong/2\nC\n+\nI\n')
        outdir = self.root/'mismatch'
        cmd = [sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
               '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
               '--sample', 'sample', '--taxon', 'Insect species', '--tools', 'getorganelle']
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn('第 2 条 read ID 不配对', result.stderr)
        self.assertFalse(outdir.exists())
        r2.write_text('@r1/2\nA\n+\nI\n@r2/2\nC\n+\n')
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn('记录 2 不完整', result.stderr)
        self.assertFalse(outdir.exists())

    def test_illumina_candidate_orchestrator_runs_python_mitoflex_script_with_its_env(self):
        r1, r2 = self.root/'sample_R1.fastq', self.root/'sample_R2.fastq'
        record = '@read/1\nACGT\n+\nIIII\n'
        r1.write_text(record)
        r2.write_text(record.replace('/1', '/2'))
        script = self.root/'MitoFlex.py'
        script.write_text(
            'import pathlib, sys\n'
            'assert sys.argv[1] == "all"\n'
            'assert "--workname" in sys.argv and "--fastq1" in sys.argv and "--fastq2" in sys.argv\n'
            'assert "--disable-taxa" in sys.argv\n'
            'assert "--use-list" not in sys.argv\n'
            'base = pathlib.Path("sample/sample.temp/findmitoscaf")\n'
            'base.mkdir(parents=True, exist_ok=True)\n'
            '(base / "sample.picked.fa").write_text(">mt\\nACGTACGT\\n")\n'
            'print("MitoFlex v1.2.3")\n')
        interpreter = self.root/'mitoflex-python'
        interpreter.write_text(
            '#!/bin/sh\n'
            'if [ "$1" = "-B" ]; then exit 0; fi\n'
            'exec /usr/bin/python3 "$@"\n')
        interpreter.chmod(0o755)
        outdir = self.root/'sample_mitogenome'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
            '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
            '--sample', 'sample', '--taxon', 'Insect species', '--tools', 'mitoflex',
            '--mitoflex', str(script), '--mitoflex-python', str(interpreter),
        ], capture_output=True, text=True, env=dict(os.environ, CONDA_ROOT=str(self.root/'empty-conda-root')))
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        manifest = json.loads((outdir/'intermediate/assemblies/assembly_candidates.json').read_text())
        run = manifest['tools'][0]
        self.assertEqual(run['status'], 'succeeded', (outdir/'intermediate/assemblies/mitoflex/tool.log').read_text())
        self.assertEqual(run['command'][:2], [str(interpreter), str(script)])
        self.assertIn('--clade', run['command'])
        self.assertIn('--species-name', run['command'])
        self.assertIn('--disable-taxa', run['command'])
        self.assertNotIn('--use-list', run['command'])
        self.assertEqual(manifest['candidate_count'], 1, json.dumps(manifest, indent=2) +
                         '\nLOG=' + (outdir/'intermediate/assemblies/mitoflex/tool.log').read_text())

    def test_illumina_candidate_orchestrator_exposes_helpers_from_other_conda_env(self):
        r1, r2 = self.root/'r1.fastq', self.root/'r2.fastq'
        r1.write_text('@r/1\nACGT\n+\nIIII\n')
        r2.write_text('@r/2\nACGT\n+\nIIII\n')
        conda = self.root/'conda'
        executable_bin = conda/'envs/a/bin'
        helper_bin = conda/'envs/b/bin'
        executable_bin.mkdir(parents=True)
        helper_bin.mkdir(parents=True)
        helper = helper_bin/'mito-helper'
        helper.write_text('#!/bin/sh\necho helper-found\n')
        helper.chmod(0o755)
        exe = executable_bin/'get_organelle_from_reads.py'
        exe.write_text(
            '#!/bin/sh\n'
            'if [ "$1" = "--version" ]; then echo "GetOrganelle 1.0"; exit 0; fi\n'
            'command -v mito-helper >/dev/null || exit 7\n'
            'out=\nwhile [ "$#" -gt 0 ]; do if [ "$1" = "-o" ]; then out="$2"; shift 2; else shift; fi; done\n'
            'mkdir -p "$out"\nprintf ">mt\\nACGTACGT\\n" > "$out/x.path_sequence.fasta"\n')
        exe.chmod(0o755)
        outdir = self.root/'multi-env'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
            '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
            '--sample', 'sample', '--taxon', 'Insect species', '--tools', 'getorganelle',
        ], capture_output=True, text=True, env=dict(os.environ, CONDA_ROOT=str(conda)))
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        manifest = json.loads((outdir/'intermediate/assemblies/assembly_candidates.json').read_text())
        self.assertEqual(manifest['tools'][0]['status'], 'succeeded')
        self.assertEqual(manifest['candidate_count'], 1)

    def test_illumina_candidate_orchestrator_invokes_nonexecutable_novoplasty_with_perl(self):
        r1, r2 = self.root/'sample_R1.fastq', self.root/'sample_R2.fastq'
        record = '@read/1\nACGT\n+\nIIII\n'
        r1.write_text(record)
        r2.write_text(record.replace('/1', '/2'))
        seed = self.root/'seed.fa'
        seed.write_text('>seed\nACGT\n')
        script = self.root/'NOVOPlasty4.3.5.pl'
        script.write_text('placeholder; invoke with Perl\n')
        perl = self.root/'perl'
        perl.write_text(
            '#!/bin/sh\n'
            'config="$3"\n'
            'output_path=$(sed -n "s/^Output path[[:space:]]*=[[:space:]]*//p" "$config")\n'
            'case "$output_path" in */) ;; *) echo "missing output path separator" >&2; exit 23 ;; esac\n'
            'printf "NOVOPlasty: The Organelle Assembler\\nVersion 4.3.5\\n"\n'
            'printf ">mt_circular\\nACGTACGT\\n" > circular_contigs.fasta\n')
        perl.chmod(0o755)
        outdir = self.root/'sample_mitogenome'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_illumina_candidates.py'),
            '--r1', str(r1), '--r2', str(r2), '--outdir', str(outdir),
            '--sample', 'sample', '--taxon', 'Insect species', '--tools', 'novoplasty',
            '--novoplasty', str(script), '--novo-seed', str(seed), '--insert-size', '300',
            '--genome-range', '10000-60000',
        ], capture_output=True, text=True,
           env=dict(os.environ, CONDA_ROOT=str(self.root/'empty-conda-root'), PATH=str(self.root)+os.pathsep+os.environ.get('PATH', '')))
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        manifest = json.loads((outdir/'intermediate/assemblies/assembly_candidates.json').read_text())
        run = manifest['tools'][0]
        self.assertEqual(run['status'], 'succeeded', json.dumps(run, indent=2))
        self.assertEqual(run['version_from_log'], '4.3.5')
        self.assertEqual(run['command'][:2], [str(perl), str(script)])
        self.assertEqual(manifest['candidate_count'], 1)

    def test_candidate_annotation_runs_mitoz_without_selecting_a_winner(self):
        candidate = self.root/'candidate.fa'
        candidate.write_text('>mt\nACGTACGT\n')
        mito = self.root/'mitoz'
        mito.write_text('#!/bin/sh\nmkdir -p annotation_out\nprintf "##gff-version 3\\n" > annotation.gff3\n'
                        'printf "LOCUS       mt 8 bp DNA circular\\nORIGIN\\n        1 acgtacgt\\n//\\n" > annotation.gbf\n'
                        'echo "invoked $*"\n')
        mito.chmod(0o755)
        outdir = self.root/'annotations'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_candidate_annotations.py'),
            '--candidate-fasta', str(candidate), '--outdir', str(outdir),
            '--taxon', 'Insect species', '--table', '5', '--table-status', 'confirmed',
            '--tools', 'mitoz', '--mitoz', str(mito),
        ], capture_output=True, text=True, env=dict(os.environ, CONDA_ROOT=str(self.root/'empty-conda-root')))
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        manifest = json.loads((outdir/'annotation_candidates.json').read_text())
        self.assertEqual(manifest['annotation_decision'], 'not_selected')
        self.assertEqual(manifest['runs'][0]['status'], 'succeeded')
        self.assertTrue(any(path.endswith('.gbf') for path in manifest['runs'][0]['annotation_outputs']))
        self.assertIn('--genetic_code', manifest['runs'][0]['command'])
        self.assertEqual(candidate.read_text(), '>mt\nACGTACGT\n')

    def test_candidate_annotation_keeps_same_named_identical_inputs_separate(self):
        left = self.root/'assembler_a'/'candidate.fa'
        right = self.root/'assembler_b'/'candidate.fa'
        left.parent.mkdir(); right.parent.mkdir()
        left.write_text('>mt\nACGTACGT\n')
        right.write_text(left.read_text())
        mito = self.root/'mitoz'
        mito.write_text('#!/bin/sh\nprintf "##gff-version 3\\n" > annotation.gff3\necho "invoked $*"\n')
        mito.chmod(0o755)
        outdir = self.root/'annotations'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_candidate_annotations.py'),
            '--candidate-fasta', str(left), '--candidate-fasta', str(right),
            '--outdir', str(outdir), '--taxon', 'Insect species', '--table', '5',
            '--tools', 'mitoz', '--mitoz', str(mito),
        ], capture_output=True, text=True,
           env=dict(os.environ, CONDA_ROOT=str(self.root/'empty-conda-root')))
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        manifest = json.loads((outdir/'annotation_candidates.json').read_text())
        self.assertEqual(len(manifest['runs']), 2)
        self.assertEqual({run['status'] for run in manifest['runs']}, {'succeeded'})
        self.assertEqual(len({run['workdir'] for run in manifest['runs']}), 2)

    def test_candidate_annotation_aliases_long_mitoz_seqids_and_records_mapping(self):
        candidate = self.root/'candidate.fa'
        candidate.write_text('>a_sequence_id_longer_than_13\nACGT\n')
        mito = self.root/'mitoz'
        mito.write_text('#!/bin/sh\nprintf "##gff-version 3\\n" > annotation.gff3\necho "invoked $*"\n')
        mito.chmod(0o755)
        outdir = self.root/'annotations'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_candidate_annotations.py'),
            '--candidate-fasta', str(candidate), '--outdir', str(outdir),
            '--taxon', 'Insect species', '--table', '5', '--tools', 'mitoz', '--mitoz', str(mito),
        ], capture_output=True, text=True, env=dict(os.environ, CONDA_ROOT=str(self.root/'empty-conda-root')))
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        run = json.loads((outdir/'annotation_candidates.json').read_text())['runs'][0]
        self.assertEqual(run['status'], 'succeeded')
        self.assertEqual(run['sequence_id_alias']['original_header'], 'a_sequence_id_longer_than_13')
        self.assertEqual(run['sequence_id_alias']['mitoz_id'], 'mt01')

    def test_candidate_annotation_can_split_multirecord_assembly_with_provenance(self):
        candidate = self.root/'candidate.fa'
        content = '>contig_A description\nACGTACGT\n>contig_B\nTTTTCCCC\n'
        candidate.write_text(content)
        manifest_path = self.root/'assembly_candidates.json'
        manifest_path.write_text(json.dumps({
            'format': 'mito-illumina-assembly-candidates-1',
            'tools': [{'tool': 'getorganelle', 'status': 'succeeded'}],
            'candidates': [{'path': str(candidate), 'producer': 'getorganelle',
                            'sha256': __import__('hashlib').sha256(content.encode()).hexdigest()}],
        }))
        mito = self.root/'mitoz'
        mito.write_text('#!/bin/sh\nprintf "##gff-version 3\\n" > annotation.gff3\necho "invoked $*"\n')
        mito.chmod(0o755)
        outdir = self.root/'annotations'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_candidate_annotations.py'),
            '--assembly-manifest', str(manifest_path), '--outdir', str(outdir),
            '--taxon', 'Insect species', '--table', '5', '--tools', 'mitoz',
            '--mitoz', str(mito), '--split-multi-records',
        ], capture_output=True, text=True, env=dict(os.environ, CONDA_ROOT=str(self.root/'empty-conda-root')))
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        manifest = json.loads((outdir/'annotation_candidates.json').read_text())
        self.assertEqual(len(manifest['runs']), 2)
        self.assertEqual({run['candidate_provenance']['source_record_id'] for run in manifest['runs']},
                         {'contig_A', 'contig_B'})
        self.assertEqual({run['status'] for run in manifest['runs']}, {'succeeded'})
        self.assertTrue(all('不能证明完整' in run['candidate_provenance']['interpretation_limit']
                            for run in manifest['runs']))

    def test_candidate_annotation_rejects_manifest_fasta_hash_drift(self):
        candidate = self.root/'candidate.fa'
        original = '>mt\nACGTACGT\n'
        candidate.write_text(original)
        manifest_path = self.root/'assembly_candidates.json'
        manifest_path.write_text(json.dumps({
            'format': 'mito-illumina-assembly-candidates-1',
            'tools': [{'tool': 'getorganelle', 'status': 'succeeded'}],
            'candidates': [{'path': str(candidate), 'producer': 'getorganelle',
                            'sha256': __import__('hashlib').sha256(original.encode()).hexdigest()}],
        }))
        candidate.write_text('>mt\nTTTTTTTT\n')
        outdir = self.root/'annotations'
        result = subprocess.run([
            sys.executable, str(ROOT/'scripts/run_candidate_annotations.py'),
            '--assembly-manifest', str(manifest_path), '--outdir', str(outdir),
            '--taxon', 'Insect species', '--table', '5', '--tools', 'mitoz',
        ], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn('SHA256 不符', result.stderr)
        self.assertFalse(outdir.exists())

    def test_real_good_coverage(self):
        result = self.run_depth(self.bam())
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        report = json.loads((self.root/'depth.json').read_text())
        self.assertEqual(report['global_mean'], 10)

    def test_real_missing_quality_is_not_high_quality(self):
        result = self.run_depth(self.bam(missing_quality=True))
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        report = json.loads((self.root/'depth.json').read_text())
        self.assertEqual(report['global_mean'], 0)

    def test_real_unknown_mapq_excluded(self):
        result = self.run_depth(self.bam(mapq=255))
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        self.assertEqual(json.loads((self.root/'depth.json').read_text())['global_mean'], 0)

    def test_real_junction_does_not_count_gap_spans(self):
        if not shutil.which('samtools'):
            self.skipTest('samtools not installed')
        sys.path.insert(0, str(ROOT/'scripts'))
        from circularize import junction_evidence
        lines = ['@HD\tVN:1.6\tSO:coordinate', '@SQ\tSN:mt\tLN:300']
        for name, pos, cigar in [('deletion', 90, '10M100D10M'), ('skip', 90, '10M100N10M'),
                                 ('continuous', 95, '20M')]:
            lines.append('\t'.join([name, '0', 'mt', str(pos), '60', cigar,
                                     '*', '0', '0', 'A'*20, 'I'*20]))
        sam, bam = self.root/'junction.sam', self.root/'junction.bam'
        sam.write_text('\n'.join(lines)+'\n')
        subprocess.run(['samtools', 'sort', '-o', str(bam), str(sam)], check=True, capture_output=True)
        subprocess.run(['samtools', 'index', str(bam)], check=True, capture_output=True)
        evidence = junction_evidence(str(bam), 'mt:100-110')
        self.assertEqual(evidence['support'], 1)
        self.assertEqual(evidence['names'], ['continuous'])

    def test_real_reference_mismatch_is_error(self):
        result = self.run_depth(self.bam(chrom='other'))
        self.assertEqual(result.returncode, 1, result.stdout+result.stderr)
        self.assertFalse((self.root/'depth.json').exists())

    def test_real_last_window_weighting(self):
        result = self.run_depth(self.bam(length=501, covered=500), length=501)
        self.assertEqual(result.returncode, 2, result.stdout+result.stderr)
        report = json.loads((self.root/'depth.json').read_text())
        self.assertAlmostEqual(report['global_mean'], 5000/501)

    def test_real_missing_index_does_not_become_negative_evidence(self):
        bam = self.bam(count=0)
        Path(str(bam)+'.bai').unlink()
        result = self.run_depth(bam)
        self.assertEqual(result.returncode, 1, result.stdout+result.stderr)
        self.assertFalse((self.root/'depth.json').exists())

    def test_blast_two_copies_are_ambiguous_and_one_copy_works(self):
        if not all(shutil.which(tool) for tool in ('blastn', 'makeblastdb')):
            self.skipTest('BLAST+ not installed')
        rng = random.Random(713)
        sequence = ''.join(rng.choice('ACGT') for _ in range(300))
        record = SeqRecord(Seq(sequence), id='ref', annotations={'molecule_type': 'DNA'})
        record.features = [SeqFeature(FeatureLocation(0, len(sequence), strand=1),
                                      type='CDS', qualifiers={'gene': ['cox1']})]
        ref = self.root/'ref.gb'
        SeqIO.write(record, ref, 'genbank')
        reverse = str(Seq(sequence).reverse_complement())
        for name, sequences in [('forward', [sequence]), ('reverse', [reverse]),
                                ('two_forward', [sequence, sequence]), ('both', [sequence, reverse])]:
            with self.subTest(case=name):
                target = self.root/'target.fa'
                target.write_text('>mt\n'+('N'*100).join(sequences)+'\n')
                output = self.root/('order_%s.tsv' % name)
                result = subprocess.run([sys.executable, str(ROOT/'scripts/blast_genes.py'),
                    str(ref), str(target), '--out', str(output)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0 if len(sequences) == 1 else 1, result.stdout+result.stderr)
                self.assertEqual(output.exists(), len(sequences) == 1)
                if name == 'reverse':
                    self.assertEqual(output.read_text().split('\t')[4], '-')

    def test_blast_repeated_reference_names_use_coordinate_labels(self):
        if not all(shutil.which(tool) for tool in ('blastn', 'makeblastdb')):
            self.skipTest('BLAST+ not installed')
        rng = random.Random(9182)
        sequence = ''.join(rng.choice('ACGT') for _ in range(900))
        record = SeqRecord(Seq(sequence), id='ref', annotations={'molecule_type': 'DNA'})
        record.features = [
            SeqFeature(FeatureLocation(100, 250, strand=1), type='tRNA',
                       qualifiers={'product': ['tRNA-Leu']}),
            SeqFeature(FeatureLocation(500, 650, strand=-1), type='tRNA',
                       qualifiers={'product': ['tRNA-Leu']}),
        ]
        ref = self.root/'duplicate_reference.gb'
        SeqIO.write(record, ref, 'genbank')
        target = self.root/'target.fa'
        target.write_text('>mt\n'+sequence+'\n')
        output = self.root/'order.tsv'
        result = subprocess.run([sys.executable, str(ROOT/'scripts/blast_genes.py'),
            str(ref), str(target), '--out', str(output)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        rows = [line.split('\t') for line in output.read_text().splitlines()]
        self.assertTrue(all(len(row) == 8 for row in rows))
        self.assertTrue(all(row[7] == 'mt' for row in rows))
        self.assertEqual([row[0] for row in rows], ['tRNA-Leu@101..250:+', 'tRNA-Leu@501..650:-'])
        self.assertEqual([(row[2], row[3], row[4]) for row in rows],
                         [('101', '250', '+'), ('501', '650', '-')])


if __name__ == '__main__':
    unittest.main()

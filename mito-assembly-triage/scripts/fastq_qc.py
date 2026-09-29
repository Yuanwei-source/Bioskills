#!/usr/bin/env python3
"""Run read-only FASTQ quality reporting with fastp; never replace raw reads."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


FORMAT = 'mito-fastq-quality-report-1'


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description='FASTQ 双端质量报告；原始数据只读')
    parser.add_argument('--assembly-manifest', required=True,
                        help='由 run_illumina_candidates.py 生成，提供已配对验证与 reads SHA256')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--fastp', help='fastp 可执行文件；默认从 PATH 查找')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args(argv)
    try:
        if args.threads < 1:
            raise ValueError('--threads 必须为正整数')
        manifest_path = Path(args.assembly_manifest).resolve()
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        if manifest.get('format') != 'mito-illumina-assembly-candidates-1':
            raise ValueError('不支持的组装 manifest 格式')
        inputs = {entry.get('role'): entry for entry in manifest.get('inputs', [])}
        if not {'R1', 'R2'} <= set(inputs) or not manifest.get('read_count'):
            raise ValueError('manifest 缺少配对 reads 路径或已验证 read_count')
        reads = {role: Path(inputs[role]['path']).resolve() for role in ('R1', 'R2')}
        for role, path in reads.items():
            if not path.is_file():
                raise ValueError('%s 文件不存在: %s' % (role, path))
            if sha256_file(path) != inputs[role].get('sha256'):
                raise ValueError('%s FASTQ 与组装 manifest SHA256 不符: %s' % (role, path))
        fastp = args.fastp or shutil.which('fastp')
        if not fastp:
            raise ValueError('未找到 fastp；请指定 --fastp')
        fastp = str(Path(fastp).resolve())
        if not Path(fastp).is_file() or not os.access(fastp, os.X_OK):
            raise ValueError('fastp 不存在或不可执行: %s' % fastp)
        outdir = Path(args.outdir).resolve()
        if outdir.exists():
            raise ValueError('为保护已有结果，输出目录必须不存在: %s' % outdir)
        outdir.mkdir(parents=True)
        json_path, html_path = outdir / 'fastp.json', outdir / 'fastp.html'
        log_path = outdir / 'fastp.log'
        command = [fastp, '-i', str(reads['R1']), '-I', str(reads['R2']), '--stdout',
                   '--json', str(json_path), '--html', str(html_path),
                   '--report_title', '%s raw paired-end FASTQ' % manifest.get('sample', 'sample'),
                   '--thread', str(args.threads)]
        started = time.monotonic()
        with log_path.open('w', encoding='utf-8') as log:
            # fastp's --stdout sends the full transformed FASTQ stream to stdout;
            # discard that stream and retain only its diagnostics in the log.
            proc = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=log, check=False)
        result = {'format': FORMAT, 'created_at': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
                  'status': 'succeeded' if proc.returncode == 0 else 'failed',
                  'sample': manifest.get('sample'), 'taxon': manifest.get('taxon'),
                  'expected_pairs': manifest['read_count'],
                  'inputs': {role: {'path': str(path), 'sha256': inputs[role]['sha256']}
                             for role, path in reads.items()},
                  'fastp': {'path': fastp, 'version_output': subprocess.run(
                      [fastp, '--version'], capture_output=True, text=True, check=False).stdout.strip(),
                      'command': command, 'exit_code': proc.returncode,
                      'elapsed_seconds': round(time.monotonic() - started, 2), 'log': str(log_path),
                      'json': str(json_path), 'html': str(html_path)},
                  'interpretation_limit': '本报告只汇总 FASTQ 质量与 adapter 观测；未修剪数据，fastp 输出 reads 被丢弃。'}
        if json_path.is_file():
            qc = json.loads(json_path.read_text(encoding='utf-8'))
            summary = qc.get('summary', {})
            raw_r1 = (qc.get('read1_before_filtering') or {}).get('total_reads')
            raw_r2 = (qc.get('read2_before_filtering') or {}).get('total_reads')
            raw_summary = summary.get('before_filtering', {})
            result['fastp'].update({'read1_before_filtering': qc.get('read1_before_filtering'),
                                    'read2_before_filtering': qc.get('read2_before_filtering'),
                                    'summary': summary,
                                    'raw_quality_summary': raw_summary,
                                    'filtered_quality_summary': summary.get('after_filtering'),
                                    'observed_read1_count': raw_r1,
                                    'observed_read2_count': raw_r2,
                                    'read_count_matches_pair_validation':
                                    raw_r1 == manifest['read_count'] and raw_r2 == manifest['read_count']})
            if raw_r1 != manifest['read_count'] or raw_r2 != manifest['read_count']:
                result['status'] = 'failed_read_count_mismatch'
                result['error'] = ('fastp 读取的 R1 数量与完整 FASTQ 配对检查不同；可能未处理完所有 gzip members')
        elif proc.returncode == 0:
            result['status'] = 'failed_no_report'
            result['error'] = 'fastp 退出码为 0，但没有 JSON 质量报告'
        for role, path in reads.items():
            if sha256_file(path) != inputs[role]['sha256']:
                result['status'] = 'failed_input_changed'
                result['error'] = '%s 在 QC 期间发生变化' % role
        report = outdir / 'fastq_quality.json'
        report.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        print(json.dumps({'status': result['status'], 'report': str(report),
                          'expected_pairs': result['expected_pairs'],
                          'observed_read1_count': result.get('fastp', {}).get('observed_read1_count')},
                         ensure_ascii=False))
        return 0 if result['status'] == 'succeeded' else 2
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print('FASTQ 质量检查失败: %s' % exc, file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())

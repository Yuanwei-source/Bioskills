#!/usr/bin/env python3
"""Competitively map paired reads to assembly candidates and summarize evidence.

This is an evidence generator, not an assembly selector. All candidate FASTA
records are placed in one reference so reads are assigned competitively. The
report preserves candidate provenance and deliberately makes no winner call.
"""
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


FORMAT = 'mito-competitive-candidate-evidence-1'


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def parse_fasta(path):
    records, name, seq = [], None, []
    with open(path, encoding='ascii', errors='strict') as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    sequence = ''.join(seq).upper()
                    if not sequence:
                        raise ValueError('%s 第 %d 条记录为空' % (path, len(records) + 1))
                    records.append((name, sequence))
                name, seq = line[1:].split(maxsplit=1)[0], []
                if not name:
                    raise ValueError('%s:%d FASTA 标识为空' % (path, line_no))
            else:
                if name is None:
                    raise ValueError('%s:%d 在首个 FASTA 标识前出现序列' % (path, line_no))
                seq.append(''.join(line.split()))
    if name is not None:
        sequence = ''.join(seq).upper()
        if not sequence:
            raise ValueError('%s 最后一条记录为空' % path)
        records.append((name, sequence))
    if not records:
        raise ValueError('FASTA 没有序列记录: %s' % path)
    allowed = set('ACGTNRYKMSWBDHV.-')
    for rec_name, sequence in records:
        invalid = sorted(set(sequence) - allowed)
        if invalid:
            raise ValueError('%s 记录 %s 含非法字符: %s' % (path, rec_name, ''.join(invalid)))
    return records


def resolve_executable(explicit, name):
    value = explicit or shutil.which(name)
    if not value:
        raise ValueError('未找到 %s；请通过参数显式指定可执行文件' % name)
    path = Path(value).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError('%s 不存在或不可执行: %s' % (name, path))
    return str(path)


def run_logged(command, log_path, cwd=None):
    started = time.monotonic()
    with open(log_path, 'w', encoding='utf-8') as log:
        proc = subprocess.run(command, cwd=cwd, stdout=log, stderr=subprocess.STDOUT,
                              check=False, text=True)
    if proc.returncode:
        tail = Path(log_path).read_text(encoding='utf-8', errors='replace')[-3000:]
        raise RuntimeError('%s 失败，退出码 %d；日志末尾:\n%s' % (command[0], proc.returncode, tail))
    return {'command': list(map(str, command)), 'exit_code': proc.returncode,
            'elapsed_seconds': round(time.monotonic() - started, 2), 'log': str(log_path)}


def parse_coverage(text):
    result = {}
    lines = [line for line in text.splitlines() if line and not line.startswith('#')]
    for line in lines:
        fields = line.split('\t')
        if len(fields) < 9:
            continue
        result[fields[0]] = {
            'start': int(fields[1]), 'end': int(fields[2]),
            'mapped_reads': int(fields[3]), 'covered_bases': int(fields[4]),
            'breadth_percent': float(fields[5]), 'mean_depth': float(fields[6]),
            'mean_base_quality': float(fields[7]), 'mean_mapping_quality': float(fields[8]),
        }
    return result


def summarize_primary_pairs(bam, samtools, target_to_candidate, min_mapq, log_path):
    """Count fragments, not SAM alignment rows; require concordant primary mates."""
    command = [samtools, 'view', '-F', '2304', '-q', str(min_mapq), str(bam)]
    started = time.monotonic()
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    fragments = {}
    for line in proc.stdout:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 11:
            continue
        flag, qname, target = int(fields[1]), fields[0], fields[2]
        if flag & 0x4 or target == '*' or target not in target_to_candidate:
            continue
        candidate_key = target_to_candidate[target]
        slot = fragments.setdefault(qname, {})
        slot[1 if flag & 0x40 else 2 if flag & 0x80 else 0] = (candidate_key, bool(flag & 0x2))
    stderr = proc.stderr.read()
    status = proc.wait()
    Path(log_path).write_text('$ %s\n%s' % (' '.join(command), stderr), encoding='utf-8')
    if status:
        raise RuntimeError('samtools view 失败，退出码 %d；见 %s' % (status, log_path))
    counts = {key: {'primary_mapped_fragments': 0, 'proper_pair_fragments_same_candidate': 0,
                    'both_mates_mapq_pass_same_candidate': 0} for key in set(target_to_candidate.values())}
    for mates in fragments.values():
        assigned = {item[0] for item in mates.values()}
        if len(assigned) != 1:
            continue
        key = next(iter(assigned))
        counts[key]['primary_mapped_fragments'] += 1
        if len(mates) >= 2 and all(item[1] for item in mates.values()):
            counts[key]['proper_pair_fragments_same_candidate'] += 1
        if 1 in mates and 2 in mates:
            counts[key]['both_mates_mapq_pass_same_candidate'] += 1
    return counts, {'command': command, 'exit_code': 0,
                    'elapsed_seconds': round(time.monotonic() - started, 2),
                    'log': str(log_path)}


def main(argv=None):
    parser = argparse.ArgumentParser(description='将 reads 竞争比对到多组装候选并汇总证据；不自动选择组装')
    parser.add_argument('--assembly-manifest', required=True)
    parser.add_argument('--r1', required=True)
    parser.add_argument('--r2', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--bwa', help='bwa 可执行文件；默认从 PATH 查找')
    parser.add_argument('--samtools', help='samtools 可执行文件；默认从 PATH 查找')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--min-mapq', type=int, default=20)
    parser.add_argument('--min-baseq', type=int, default=20)
    args = parser.parse_args(argv)

    try:
        if args.threads < 1 or not 0 <= args.min_mapq <= 254 or not 0 <= args.min_baseq <= 93:
            raise ValueError('线程数、MAPQ 或碱基质量阈值非法')
        manifest_path = Path(args.assembly_manifest).resolve()
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        if manifest.get('format') != 'mito-illumina-assembly-candidates-1':
            raise ValueError('不支持的组装候选 manifest 格式')
        r1, r2 = Path(args.r1).resolve(), Path(args.r2).resolve()
        for reads in (r1, r2):
            if not reads.is_file():
                raise ValueError('reads 文件不存在: %s' % reads)
        read_hashes = {'r1': sha256_file(r1), 'r2': sha256_file(r2)}
        candidates = manifest.get('candidates', [])
        if not candidates:
            raise ValueError('manifest 中没有组装候选')
        outdir = Path(args.outdir).resolve()
        if outdir.exists():
            raise ValueError('为保护已有结果，输出目录必须不存在: %s' % outdir)
        outdir.mkdir(parents=True)
        bwa = resolve_executable(args.bwa, 'bwa')
        samtools = resolve_executable(args.samtools, 'samtools')

        targets, candidate_rows, target_to_candidate = [], [], {}
        for index, entry in enumerate(candidates, 1):
            fasta = Path(entry['path']).resolve()
            if not fasta.is_file():
                raise ValueError('候选 FASTA 不存在: %s' % fasta)
            observed_sha = sha256_file(fasta)
            if entry.get('sha256') and entry['sha256'] != observed_sha:
                raise ValueError('候选 FASTA 与 manifest SHA256 不符: %s' % fasta)
            producer = entry.get('producer', 'unknown')
            candidate_key = 'candidate_%03d' % index
            row = {'candidate_key': candidate_key, 'producer': producer,
                   'source_fasta': str(fasta), 'source_sha256': observed_sha,
                   'source_candidate': entry, 'records': []}
            for rec_index, (record_id, sequence) in enumerate(parse_fasta(fasta), 1):
                target = 't%04d_%03d' % (index, rec_index)
                targets.append((target, sequence))
                target_to_candidate[target] = candidate_key
                row['records'].append({'target_id': target, 'source_record_id': record_id,
                                       'length': len(sequence),
                                       'sequence_sha256': hashlib.sha256(sequence.encode('ascii')).hexdigest()})
            candidate_rows.append(row)

        target_fasta = outdir / 'competitive_candidates.fasta'
        with target_fasta.open('w', encoding='ascii') as handle:
            for name, sequence in targets:
                handle.write('>%s\n' % name)
                for start in range(0, len(sequence), 80):
                    handle.write(sequence[start:start + 80] + '\n')
        (outdir / 'target_map.json').write_text(json.dumps(candidate_rows, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        version_bwa = subprocess.run([bwa], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, check=False).stdout.strip().splitlines()[:4]
        version_samtools = subprocess.run([samtools, '--version'], stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT, text=True, check=False).stdout.splitlines()[:2]
        commands = []
        commands.append(run_logged([bwa, 'index', str(target_fasta)], outdir / 'bwa_index.log'))
        bam = outdir / 'competitive_candidates.sorted.bam'
        bwa_cmd = [bwa, 'mem', '-t', str(args.threads), str(target_fasta), str(r1), str(r2)]
        sam_cmd = [samtools, 'sort', '-@', str(args.threads), '-o', str(bam), '-']
        started = time.monotonic()
        with (outdir / 'bwa_mem.log').open('w', encoding='utf-8') as bwa_log, (outdir / 'samtools_sort.log').open('w', encoding='utf-8') as sam_log:
            bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=bwa_log)
            sort_proc = subprocess.Popen(sam_cmd, stdin=bwa_proc.stdout, stdout=sam_log, stderr=sam_log)
            bwa_proc.stdout.close()
            sort_status = sort_proc.wait()
            bwa_status = bwa_proc.wait()
        if bwa_status or sort_status:
            raise RuntimeError('竞争比对失败 bwa=%d samtools_sort=%d；检查日志' % (bwa_status, sort_status))
        commands.append({'command': bwa_cmd + ['|'] + sam_cmd, 'exit_code': 0,
                         'elapsed_seconds': round(time.monotonic() - started, 2),
                         'logs': [str(outdir / 'bwa_mem.log'), str(outdir / 'samtools_sort.log')]})
        commands.append(run_logged([samtools, 'index', str(bam)], outdir / 'samtools_index.log'))
        cov_cmd = [samtools, 'coverage', '-q', str(args.min_mapq), '-Q', str(args.min_baseq), str(bam)]
        cov_started = time.monotonic()
        cov_proc = subprocess.run(cov_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        (outdir / 'samtools_coverage.log').write_text(cov_proc.stderr, encoding='utf-8')
        if cov_proc.returncode:
            raise RuntimeError('samtools coverage 失败；见 samtools_coverage.log')
        commands.append({'command': cov_cmd, 'exit_code': cov_proc.returncode,
                         'elapsed_seconds': round(time.monotonic() - cov_started, 2),
                         'log': str(outdir / 'samtools_coverage.log'),
                         'output': str(outdir / 'coverage.tsv')})
        coverage = parse_coverage(cov_proc.stdout)
        (outdir / 'coverage.tsv').write_text(cov_proc.stdout, encoding='utf-8')
        pair_counts, pair_command = summarize_primary_pairs(bam, samtools, target_to_candidate, args.min_mapq,
                                                            outdir / 'samtools_view.log')
        commands.append(pair_command)
        for row in candidate_rows:
            rec_cov = [coverage.get(rec['target_id']) for rec in row['records']]
            rec_cov = [item for item in rec_cov if item]
            total_len = sum(rec['length'] for rec in row['records'])
            covered = sum(item['covered_bases'] for item in rec_cov)
            row['evidence'] = {
                'total_candidate_bases': total_len,
                'covered_bases_mapq_ge_%d_baseq_ge_%d' % (args.min_mapq, args.min_baseq): covered,
                'length_weighted_breadth_percent': round(100 * covered / total_len, 4) if total_len else 0,
                'length_weighted_mean_depth': round(sum(item['mean_depth'] * (item['end'] - item['start'] + 1)
                                                         for item in rec_cov) / total_len, 4) if total_len else 0,
                'records': [{'target_id': rec['target_id'], 'source_record_id': rec['source_record_id'],
                             'length': rec['length'], **coverage.get(rec['target_id'], {})}
                            for rec in row['records']],
                **pair_counts[row['candidate_key']],
            }
        for role, path in (('r1', r1), ('r2', r2)):
            if sha256_file(path) != read_hashes[role]:
                raise RuntimeError('%s 在竞争比对期间发生变化；本次结果不予登记为成功' % role.upper())
        report = {
            'format': FORMAT, 'created_at': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
            'inputs': {'assembly_manifest': str(manifest_path), 'assembly_manifest_sha256': sha256_file(manifest_path),
                       'r1': {'path': str(r1), 'sha256': read_hashes['r1']},
                       'r2': {'path': str(r2), 'sha256': read_hashes['r2']}},
            'tools': {'bwa': {'path': bwa, 'version_output': version_bwa},
                      'samtools': {'path': samtools, 'version_output': version_samtools}},
            'parameters': {'threads': args.threads, 'minimum_mapq': args.min_mapq,
                           'minimum_base_quality': args.min_baseq, 'mapping': 'BWA-MEM primary/best alignments'},
            'commands': commands, 'competitive_reference': str(target_fasta),
            'competitive_reference_sha256': sha256_file(target_fasta),
            'bam': str(bam), 'bam_sha256': sha256_file(bam), 'candidates': candidate_rows,
            'decision': 'not_selected',
            'interpretation_limits': [
                '覆盖与配对比对证据是候选相对支持度，不证明序列正确、完整、无 NUMT 或真实环化。',
                '当前统计按输入 FASTA 线性比对；未为圆形首尾边界补齐参考，因此不用于评价环化接缝支持。',
                '仅保留 BWA-MEM 主比对；高相似重复区域的竞争不确定性应结合 MAPQ 与原始图结构复核。',
                '来自同一软件或同一候选文件的 contig 不构成独立组装证据。',
                '不得单凭长度加权覆盖、reads 数或软件来源自动选定最终组装。',
            ],
        }
        (outdir / 'candidate_evidence.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        print(json.dumps({'status': 'succeeded', 'report': str(outdir / 'candidate_evidence.json'),
                          'bam': str(bam), 'candidate_count': len(candidate_rows)}, ensure_ascii=False))
        return 0
    except (OSError, ValueError, KeyError, RuntimeError, json.JSONDecodeError) as exc:
        print('候选竞争比对失败: %s' % exc, file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())

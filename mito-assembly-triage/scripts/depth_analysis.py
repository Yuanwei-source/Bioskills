#!/usr/bin/env python3
"""
reads 裁判: 覆盖度剖面 + mate 分布 + soft-clip 分析
用法: python3 depth_analysis.py <sample.bam> <genome.fasta> [--window 500]
输入: 已排序索引的 bam (bwa mem 产出), 参考 fasta
输出: 覆盖度剖面表 + 异常区报告
"""
import sys, subprocess, collections, re, os
import argparse
import json
import tempfile

# 前置门禁：本步骤所需依赖（唯一清单来源 config/dependencies.json）
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from _deps import require_stage  # noqa: E402
from _sam import cigar_operations, template_key
from _evidence_io import read_single_fasta, protect_output, atomic_write_text, InputArgumentParser


def reference_to_query_offset(cigar, reference_offset):
    reference_cursor = 0
    query_cursor = 0
    for length, operation in cigar_operations(cigar):
        if operation in 'M=X':
            if reference_cursor <= reference_offset < reference_cursor + length:
                return query_cursor + reference_offset - reference_cursor
            reference_cursor += length
            query_cursor += length
        elif operation in 'IS':
            query_cursor += length
        elif operation in 'DN':
            reference_cursor += length
    return None


def depth_quality_gate(global_mean, low):
    return global_mean > 0 and not low


def coverage_intervals(length, window):
    if length <= 0 or window <= 0:
        return []
    return [(start, min(start + window - 1, length)) for start in range(1, length + 1, window)]


def read_base_and_quality_at_reference(fields, reference_position):
    """Return (base, base_quality) in reference orientation, or (None, None)."""
    if len(fields) < 11:
        raise ValueError('SAM record has fewer than 11 fields')
    reference_offset = reference_position - int(fields[3])
    query_offset = reference_to_query_offset(fields[5], reference_offset)
    if query_offset is None:
        return None, None
    sequence = fields[9]
    qualities = fields[10] if len(fields) > 10 else ''
    # SAM already stores reverse-strand SEQ/QUAL in reference orientation.
    if sequence == '*' or query_offset < 0 or query_offset >= len(sequence):
        return None, None
    expected_length = sum(n for n, op in cigar_operations(fields[5]) if op in 'MIS=X')
    if expected_length != len(sequence):
        raise ValueError('SAM CIGAR/SEQ lengths differ')
    if qualities != '*' and len(qualities) != len(sequence):
        raise ValueError('SAM SEQ/QUAL lengths differ')
    base = sequence[query_offset].upper()
    quality = (ord(qualities[query_offset]) - 33) if qualities != '*' else None
    if quality is not None and not 0 <= quality <= 93:
        raise ValueError('invalid SAM base quality')
    return base, quality


def read_base_at_reference(fields, reference_position):
    return read_base_and_quality_at_reference(fields, reference_position)[0]


def base_support(records, position, min_mapq=20, min_baseq=20, min_depth=5,
                 max_strand_fraction=0.9):
    """Aggregate read-level evidence for one reference position.

    ``callable`` is True only when MAPQ, base quality, depth and strand balance all
    pass; otherwise the caller must not propose base replacement.  These four
    numbers are the documented evidence fields for a single-base repair.
    """
    bases = collections.Counter()
    strands = {'+': collections.Counter(), '-': collections.Counter()}
    excluded = dict.fromkeys(('duplicate', 'low_mapq', 'low_baseq', 'missing_baseq',
                              'no_base', 'secondary', 'qc_fail', 'mate_conflict'), 0)
    templates = collections.defaultdict(list)
    read_depth = 0
    for fields in records:
        if len(fields) < 11:
            raise ValueError('SAM record has fewer than 11 fields')
        flag = int(fields[1]); mapq = int(fields[4])
        if flag & 0x400:
            excluded['duplicate'] += 1; continue
        if flag & (0x4 | 0x100 | 0x800):
            excluded['secondary'] += 1; continue
        if flag & 0x200:
            excluded['qc_fail'] += 1; continue
        if mapq == 255 or mapq < min_mapq:
            excluded['low_mapq'] += 1; continue
        base, quality = read_base_and_quality_at_reference(fields, position)
        if base not in ('A', 'C', 'G', 'T'):
            excluded['no_base'] += 1; continue
        if quality is None:
            excluded['missing_baseq'] += 1; continue
        if quality < min_baseq:
            excluded['low_baseq'] += 1; continue
        templates[template_key(fields)].append((base, quality, flag, mapq))
        read_depth += 1
    template_mapqs, template_baseqs = [], []
    for calls in templates.values():
        if len({c[0] for c in calls}) != 1:
            excluded['mate_conflict'] += 1
            continue
        # Prefer R1 when both mates cover the locus, independent of SAM row order.
        # The chosen alignment supplies one strand observation per template.
        base, quality, flag, mapq = min(calls, key=lambda c: (not bool(c[2] & 0x40), -c[1], c[2]))
        bases[base] += 1
        strands['-' if flag & 0x10 else '+'][base] += 1
        template_mapqs.append(mapq)
        template_baseqs.append(quality)

    depth = sum(bases.values())
    dominant, dominant_count = bases.most_common(1)[0] if bases else (None, 0)
    support = dominant_count / depth if depth else 0.0
    reasons = []
    if depth == 0:
        reasons.append('深度为 0 (无可用 reads)')
    elif depth < min_depth:
        reasons.append('深度 %d < %d' % (depth, min_depth))
    if dominant_count:
        on_dominant = max(strands['+'][dominant], strands['-'][dominant])
        fraction = on_dominant / dominant_count
        if fraction > max_strand_fraction:
            reasons.append('链向偏倚: 优势碱基单链占比 %.2f > %.2f' % (fraction, max_strand_fraction))
    return {'position': position, 'base': dominant, 'depth': depth, 'support': support,
            'template_depth': depth, 'read_depth': read_depth, 'base_counts': dict(bases),
            'mapq': template_mapqs, 'base_quality': template_baseqs,
            'counting_unit': 'read_group+qname', 'duplicate_marking': 'not_verified',
            'strand_counts': {'+': dict(strands['+']), '-': dict(strands['-'])},
            'excluded': excluded, 'callable': not reasons, 'reasons': reasons}

def run_samtools(*args):
    result = subprocess.run(['samtools', *map(str, args)], capture_output=True, text=True,
                            check=True)
    return result.stdout


def reference_header(bam, chrom, length):
    header = run_samtools('view', '-H', bam)
    matches = []
    for line in header.splitlines():
        if line.startswith('@SQ\t'):
            tags = dict(field.split(':', 1) for field in line.split('\t')[1:] if ':' in field)
            if tags.get('SN') == chrom:
                matches.append(int(tags.get('LN', '0')))
    if matches != [length]:
        raise ValueError('BAM @SQ 与 FASTA 名称/长度不匹配: %s/%d' % (chrom, length))
    # Name and length do not prove sequence identity; report this limit explicitly.


def main():
    parser = InputArgumentParser(description=__doc__)
    parser.add_argument('bam')
    parser.add_argument('fasta')
    parser.add_argument('ambiguous_fasta', nargs='?', help='同坐标、仅模糊位点不同的 FASTA')
    parser.add_argument('--window', type=int, default=500)
    parser.add_argument('--min-depth', type=int, default=5, help='工程预警阈值，非生物学定律')
    parser.add_argument('--min-mapq', type=int, default=20)
    parser.add_argument('--min-baseq', type=int, default=20)
    parser.add_argument('--allow-base-replacement', action='store_true', help='仅输出候选建议，不改写序列')
    parser.add_argument('--output-json')
    args = parser.parse_args()
    if args.window <= 0 or args.min_depth <= 0 or not 0 <= args.min_mapq <= 254 or not 0 <= args.min_baseq <= 93:
        parser.error('window/min-depth 必须为正数，MAPQ 为 0–254，baseQ 为 0–93')
    require_stage('read_evidence', __file__)
    try:
        chrom, sequence = read_single_fasta(args.fasta)
        length = len(sequence)
        inputs = [args.bam, args.fasta]
        ambiguous = None
        if args.ambiguous_fasta:
            amb_name, ambiguous = read_single_fasta(args.ambiguous_fasta)
            inputs.append(args.ambiguous_fasta)
            if amb_name != chrom or len(ambiguous) != length or any(
                    a != b and a in 'ACGT' for a, b in zip(ambiguous, sequence)):
                raise ValueError('额外 FASTA 必须与 BAM 参考同名、同长度，且仅模糊碱基位点允许不同')
        if args.output_json:
            protect_output(args.output_json, inputs)
        reference_header(args.bam, chrom, length)
        # depth -s avoids counting overlapping mates twice. It is a coverage screen,
        # not the template-aware allele consensus implemented by base_support().
        with tempfile.TemporaryDirectory(prefix='mito-depth-') as directory:
            filtered = os.path.join(directory, 'quality.bam')
            # Missing BAM QUAL is encoded as 255, not a low Phred score. Filter
            # it explicitly before depth; raw qual in htslib expressions is 0-based.
            run_samtools('view', '-b', '-q', args.min_mapq, '-F', '0xF04',
                         '-e', 'mapq < 255 && length(qual) > 0 && max(qual) <= 93',
                         '-o', filtered, args.bam, chrom)
            output = run_samtools('depth', '-a', '-s', '-q', args.min_baseq, filtered)
        depths = {}
        for line in output.splitlines():
            fields = line.split('\t')
            if len(fields) != 3:
                raise ValueError('无法解析 samtools depth 输出')
            if fields[0] != chrom:
                continue
            pos, value = int(fields[1]), int(fields[2])
            if not 1 <= pos <= length or value < 0 or pos in depths:
                raise ValueError('非法或重复的 depth 坐标')
            depths[pos] = value
        global_mean = sum(depths.values()) / length
        intervals = coverage_intervals(length, args.window)
        windows, low = [], []
        for lo, hi in intervals:
            values = [depths.get(pos, 0) for pos in range(lo, hi + 1)]
            mean = sum(values) / len(values)
            item = dict(start=lo, end=hi, mean=mean, minimum=min(values), maximum=max(values))
            item['low'] = mean < args.min_depth or mean < 0.5 * global_mean
            windows.append(item)
            if item['low']:
                low.append((lo, hi))
        zero_bases = sum(depths.get(pos, 0) == 0 for pos in range(1, length+1))
        report = dict(format='mito-depth-1', reference=chrom, length=length, global_mean=global_mean,
                      zero_coverage_bases=zero_bases, windows=windows, local_observations=[],
                      base_support=[], reads_support='READS_CONSISTENT' if global_mean > 0 else 'NOT_ASSESSED',
                      parameters=dict(window=args.window, min_depth=args.min_depth,
                                      min_mapq=args.min_mapq, min_baseq=args.min_baseq),
                      limitations=['BAM 名称/长度已校验，序列身份需由比对 provenance 核验',
                                   '单候选比对不排除 NUMT；重复标记状态未自动验证',
                                   '覆盖筛查不是碱基共识；缺失质量与 MAPQ=255 已排除'])
        print('全局均值: %.2f×；零覆盖碱基: %d/%d' % (global_mean, zero_bases, length))
        print('位置\tmean\tmin\tmax\t状态')
        for item in windows:
            print('%d-%d\t%.2f\t%d\t%d\t%s' % (item['start'], item['end'], item['mean'],
                  item['minimum'], item['maximum'], '低覆盖待查' if item['low'] else '无窗口预警'))
        for lo, hi in low[:3]:
            records = [line.split('\t') for line in run_samtools(
                'view', args.bam, '%s:%d-%d' % (chrom, lo, hi)).splitlines() if line]
            mates = collections.Counter()
            total = clips = 0
            for fields in records:
                if len(fields) < 11:
                    raise ValueError('SAM record has fewer than 11 fields')
                flag, mapq = int(fields[1]), int(fields[4])
                if flag & 0xF04 or mapq == 255 or mapq < args.min_mapq:
                    continue
                total += 1
                clips += any(op == 'S' for _, op in cigar_operations(fields[5]))
                if not flag & 1 or flag & 8:
                    mates['无已比对mate'] += 1
                elif fields[6] not in ('=', chrom):
                    mates['其他参考'] += 1
                else:
                    mpos = int(fields[7])
                    mates['左侧' if mpos < lo else '右侧' if mpos > hi else '区域内'] += 1
            item = dict(start=lo, end=hi, reads=total, soft_clipped_reads=clips, mates=dict(mates))
            report['local_observations'].append(item)
            print('区域 %d-%d: reads=%d soft-clip=%d；仅为观测，需区分重复、表示方式、异质性与错接' %
                  (lo, hi, total, clips))
        if ambiguous is not None:
            for pos, base in enumerate(ambiguous, 1):
                if base in 'ACGT':
                    continue
                records = [line.split('\t') for line in run_samtools(
                    'view', args.bam, '%s:%d-%d' % (chrom, pos, pos)).splitlines() if line]
                evidence = base_support(records, pos, args.min_mapq, args.min_baseq, args.min_depth)
                evidence['original_base'] = base
                evidence['candidate_base'] = (evidence['base'] if args.allow_base_replacement
                    and evidence['callable'] and evidence['support'] > 0.95 else None)
                report['base_support'].append(evidence)
                print('位置 %d (%s): 模板=%d 有效读段=%d 优势碱基=%s 支持率=%.3f；%s' %
                      (pos, base, evidence['depth'], evidence['read_depth'], evidence['base'],
                       evidence['support'], '; '.join(evidence['reasons']) or '通过定点工程筛查'))
                if evidence['candidate_base']:
                    print('  候选碱基 %s；仍需核验竞争比对、混合等位与重复标记；未修改序列' % evidence['candidate_base'])
        review = bool(low or zero_bases)
        report['status'] = 'REVIEW' if review else 'NO_COVERAGE_WARNING'
        if args.output_json:
            atomic_write_text(args.output_json, json.dumps(report, ensure_ascii=False, indent=2) + '\n')
        print('判定: ' + ('覆盖不足或存在零覆盖位点，需核查' if review else '未发现覆盖预警；不等于序列或结构已验证'))
        return 2 if review else 0
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, 'stderr', None) or str(exc)
        print('ERROR: reads 检查未完成: %s' % detail, file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

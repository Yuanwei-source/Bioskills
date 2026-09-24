#!/usr/bin/env python3
"""
reads 裁判: 覆盖度剖面 + mate 分布 + soft-clip 分析
用法: python3 depth_analysis.py <sample.bam> <genome.fasta> [--window 500]
输入: 已排序索引的 bam (bwa mem 产出), 参考 fasta
输出: 覆盖度剖面表 + 异常区报告
"""
import sys, subprocess, collections, re, os

_COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def reference_to_query_offset(cigar, reference_offset):
    reference_cursor = 0
    query_cursor = 0
    for length_text, operation in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(length_text)
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


def read_base_at_reference(fields, reference_position):
    reference_offset = reference_position - int(fields[3])
    query_offset = reference_to_query_offset(fields[5], reference_offset)
    if query_offset is None:
        return None
    sequence = fields[9]
    if int(fields[1]) & 0x10:
        sequence = sequence.translate(_COMPLEMENT)[::-1]
    if query_offset < 0 or query_offset >= len(sequence):
        return None
    return sequence[query_offset].upper()

def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    bam, fasta = sys.argv[1], sys.argv[2]
    window = 500
    allow_replacement = '--allow-base-replacement' in sys.argv
    if '--window' in sys.argv:
        window = int(sys.argv[sys.argv.index('--window')+1])

    # 参考名
    chrom = None
    for line in open(fasta):
        if line.startswith('>'):
            chrom = line.strip()[1:].split()[0]
            break
    if not chrom:
        print('错误: 无法读取参考名'); sys.exit(1)

    # 序列长度
    L = sum(len(l.strip()) for l in open(fasta) if not l.startswith('>'))

    # 1. 深度剖面
    depth_call = subprocess.run(['samtools', 'depth', '-a', bam], capture_output=True, text=True)
    if depth_call.returncode != 0:
        print('samtools depth 失败:', depth_call.stderr.strip() or '未知错误', file=sys.stderr)
        sys.exit(1)
    out = depth_call.stdout
    depths = {}
    for line in out.split('\n'):
        if not line.strip(): continue
        p = line.split('\t')
        if len(p) >= 3 and p[0] == chrom:
            depths[int(p[1])] = int(p[2])

    print('=== 覆盖度剖面 (window=%d) ===' % window)
    print('%-8s %10s %8s %8s' % ('位置', 'mean', 'min', 'max'))
    print('-' * 40)
    means = []
    for lo, hi in coverage_intervals(L, window):
        seg = [depths.get(j, 0) for j in range(lo, hi + 1)]
        m, mn, mx = sum(seg)/len(seg), min(seg), max(seg)
        means.append(m)
        flag = ''
        if m < 0.5 * (sum(means)/len(means) if means else m): flag = '  <-- 低覆盖!'
        print('%6d-%6d %10.1f %8d %8d%s' % (lo, hi, m, mn, mx, flag))

    if not means:
        print('⚠ 无覆盖度数据 (检查 bam 染色体名: %s)' % chrom)
        sys.exit(1)
    global_mean = sum(means) / len(means)
    print('\n全局均值: %.1f×' % global_mean)

    # 2. 低覆盖区报告
    intervals = coverage_intervals(L, window)
    low = [(lo, hi) for (lo, hi), m in zip(intervals, means) if m < 0.5 * global_mean]
    if not depth_quality_gate(global_mean, low):
        print('\n覆盖度质量门失败: 全局均值必须大于 0 且不得存在低覆盖区', file=sys.stderr)
        sys.exit(2)
    print('\n低覆盖区 (mean < %.0f×): %s' % (0.5*global_mean,
          ', '.join('%d-%d' % interval for interval in low) if low else '无'))

    # 3. mate 分布 (检查连接区)
    print('\n=== 低覆盖区 mate 分布检查 ===')
    for lo, hi in low[:3]:
        out = subprocess.run(['samtools', 'view', bam, '%s:%d-%d' % (chrom, lo, hi)],
                             capture_output=True, text=True).stdout
        mates = collections.Counter()
        n = 0
        for line in out.split('\n'):
            if not line.strip(): continue
            f = line.split('\t')
            if int(f[1]) & 0x100: continue
            mpos = int(f[7]); n += 1
            if mpos < lo: mates['左侧(<%d)' % lo] += 1
            elif mpos > hi: mates['右侧(>%d)' % hi] += 1
            else: mates['区域内'] += 1
        if n:
            print('区域 %d-%d: %d reads, mate: %s' % (lo, hi, n, dict(mates)))

    # 4. soft-clip 统计 (编码区 clip = 缺失序列; CR 区 clip = 异质性, 正常)
    print('\n=== soft-clip 检查 (重点: CR 区域) ===')
    for lo, hi in low[:3]:
        out = subprocess.run(['samtools', 'view', bam, '%s:%d-%d' % (chrom, lo, hi)],
                             capture_output=True, text=True).stdout
        clip = 0; n = 0
        for line in out.split('\n'):
            if not line.strip(): continue
            f = line.split('\t')
            if int(f[1]) & 0x100: continue
            n += 1
            if 'S' in f[5]: clip += 1
        if n:
            pct = clip/n*100
            print('区域 %d-%d: %d reads, %d 有 soft-clip (%.0f%%)%s' % (
                lo, hi, n, clip, pct,
                '  → 若在控制区(AT富集): 正常长度异质性; 若在编码区: 缺失序列!' if pct > 20 else ''))

    # 5. 模糊碱基 reads 支持 (可选: 传入含模糊碱基的序列文件)
    amb = [l for l in sys.argv if l.endswith('.fasta') and l != fasta]
    if amb:
        print('\n=== 模糊碱基 reads 支持验证 ===')
        seq = ''.join(l.strip() for l in open(amb[0]) if not l.startswith('>')).upper()
        for i, ch in enumerate(seq):
            if ch not in 'ACGT':
                pos = i + 1
                out = subprocess.run(['samtools', 'view', bam, '%s:%d-%d' % (chrom, max(1,pos-30), pos+30)],
                                     capture_output=True, text=True).stdout
                bc = collections.Counter()
                for line in out.split('\n'):
                    if not line.strip(): continue
                    f = line.split('\t')
                    if int(f[1]) & 0x100: continue
                    base = read_base_at_reference(f, pos)
                    if base:
                        bc[base] += 1
                if bc:
                    dom, cnt = bc.most_common(1)[0]
                    support = cnt / sum(bc.values())
                    recommendation = '→ 伪影, 可替换为 %s' % dom if allow_replacement and support > 0.95 else '→ 需人工确认，未授权替换'
                    print('位置 %d (%s): reads 支持 %s=%d (%.0f%%) %s' % (
                        pos, ch, dom, cnt, support * 100, recommendation))

    if low:
        print('\n判定: 存在低覆盖区，未通过 reads 质量门')
        sys.exit(2)
    print('\n判定: 未发现低覆盖区；仍需结合 mate、soft-clip 和参考证据')

if __name__ == '__main__':
    main()

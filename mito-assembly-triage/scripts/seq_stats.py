#!/usr/bin/env python3
"""
序列体检: 长度/GC/AT/模糊碱基/header
用法: python3 seq_stats.py <assembly.fasta> [--window 500]
"""
import sys
from Bio import SeqIO
from collections import Counter

def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    fn = sys.argv[1]
    window = 500
    if '--window' in sys.argv:
        window = int(sys.argv[sys.argv.index('--window')+1])

    records = list(SeqIO.parse(fn, 'fasta'))
    if not records:
        print('错误: 空序列'); sys.exit(1)
    headers = ['>' + r.description for r in records]
    seqs = [(r.id, str(r.seq).upper()) for r in records]
    seq = ''.join(s for _, s in seqs)

    c = Counter(seq)
    amb = {k: v for k, v in c.items() if k not in 'ACGT'}
    acgt = sum(c[b] for b in 'ACGT')
    gc = (c['G'] + c['C']) / acgt * 100 if acgt else 0
    at = (c['A'] + c['T']) / acgt * 100 if acgt else 0

    print('文件: %s' % fn)
    print('contig 数: %d' % len(headers))
    for h in headers[:5]:
        print('  header: %s' % h)
    print('长度: %d bp' % len(seq))
    print('GC%%: %.2f   AT%%: %.2f (仅以 A/C/G/T 为分母)' % (gc, at))
    print('碱基组成:', dict(c))
    print('\n各 contig 统计:')
    for name, contig in seqs:
        cc = Counter(contig); ca = sum(cc[b] for b in 'ACGT')
        cg = (cc['G'] + cc['C']) / ca * 100 if ca else 0
        print('  %s: length=%d GC=%.2f%% N/IUPAC=%d' % (name, len(contig), cg, len(contig) - ca))
    if amb:
        print('⚠ 模糊碱基 (%d 个):' % sum(amb.values()))
        for a, n in sorted(amb.items()):
            print('  %s x%d' % (a, n))
        # 定位模糊碱基
        for a in amb:
            pos = 0
            found = 0
            while True:
                i = seq.find(a, pos)
                if i < 0 or found >= 10: break
                print('  %s 位于 %d (0-based)' % (a, i))
                pos = i + 1; found += 1
        print('  (reads 验证后可用支持碱基替换 — 见 depth_analysis.py)')
    else:
        print('✓ 无模糊碱基')

    # 滑窗 GC 均匀性 (粗筛异常区)
    print('\n滑窗 GC (window=%d):' % window)
    for i in range(0, len(seq) - window + 1, window):
        seg = seq[i:i+window]
        g = (seg.count('G') + seg.count('C')) / len(seg) * 100
        flag = '  <-- 低GC(可能AT富集区/CR)' if g < 18 else ''
        print('  %6d-%6d: GC=%.1f%%%s' % (i+1, i+window, g, flag))

if __name__ == '__main__':
    main()

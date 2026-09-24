#!/usr/bin/env python3
"""
外科手术拼接: 将 2 个 scaffold 拼接为正确的环状线粒体基因组, 并用参考验证
用法: python3 circularize.py <scaffold1.fasta> <scaffold2.fasta> <ref.gb> <outdir> --bam <bam> --junction-region <chr:start-end> --reads-validated
原理: 分析两个 scaffold 的基因组成和方向, 尝试 4 种组合(正序/反向互补),
      选择与参考基因顺序最一致的一种, 输出环状基因组并验证。
依赖: BioPython, blastn, minimap2
"""
import sys, os, subprocess, tempfile, itertools, re


def select_unique_candidate(candidates):
    if not candidates:
        return None
    top_score = max(score for score, *_ in candidates)
    top = [candidate for candidate in candidates if candidate[0] == top_score]
    return top[0] if len(top) == 1 else None


def circular_order_matches(query, reference):
    query = [str(gene).lower() for gene in query]
    reference = [str(gene).lower() for gene in reference]
    if len(query) != len(reference) or set(query) != set(reference):
        return False
    for oriented in (reference, list(reversed(reference))):
        for start in range(len(oriented)):
            if query == oriented[start:] + oriented[:start]:
                return True
    return False


def circular_annotation_matches(query, reference):
    query = [(str(name).lower(), strand) for name, strand in query]
    reference = [(str(name).lower(), strand) for name, strand in reference]
    if len(query) != len(reference) or {name for name, _ in query} != {name for name, _ in reference}:
        return False
    for reverse in (False, True):
        oriented = []
        for name, strand in (reversed(reference) if reverse else reference):
            oriented.append((name, ('-' if strand == '+' else '+') if reverse else strand))
        for start in range(len(oriented)):
            if query == oriented[start:] + oriented[:start]:
                return True
    return False


def sam_reference_span(fields):
    start = int(fields[3])
    end = start
    for length_text, operation in re.findall(r'(\d+)([MIDNSHP=X])', fields[5]):
        length = int(length_text)
        if operation in 'M=XDN':
            end += length
    return start, end


def junction_has_spanning_read(bam, region):
    try:
        chrom, coordinates = region.rsplit(':', 1)
        lo_text, hi_text = coordinates.split('-', 1)
        lo, hi = int(lo_text), int(hi_text)
    except (ValueError, AttributeError):
        raise ValueError('junction region must be chr:start-end')
    result = subprocess.run(['samtools', 'view', bam, '%s:%d-%d' % (chrom, lo, hi)], capture_output=True, text=True, check=True)
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split('\t')
        start, end = sam_reference_span(fields)
        if start <= lo and end >= hi:
            return True
    return False


def main():
    if len(sys.argv) < 5:
        print(__doc__); sys.exit(1)
    s1, s2, ref_gb, outdir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    reads_validated = '--reads-validated' in sys.argv
    bam = sys.argv[sys.argv.index('--bam') + 1] if '--bam' in sys.argv else None
    junction_region = sys.argv[sys.argv.index('--junction-region') + 1] if '--junction-region' in sys.argv else None
    os.makedirs(outdir, exist_ok=True)
    out_fa = os.path.join(outdir, 'genome_circular.fasta')
    if os.path.exists(out_fa):
        print('输出目录已有 genome_circular.fasta，请使用新目录避免误用旧结果', file=sys.stderr)
        sys.exit(1)

    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
    except ImportError:
        print('需要 BioPython'); sys.exit(1)

    def load_seq(fn):
        rec = next(SeqIO.parse(fn, 'fasta'))
        return rec.id, str(rec.seq)

    id1, seq1 = load_seq(s1)
    id2, seq2 = load_seq(s2)
    print('scaffold1: %s (%d bp)' % (id1, len(seq1)))
    print('scaffold2: %s (%d bp)' % (id2, len(seq2)))

    gb = SeqIO.read(ref_gb, 'genbank')
    ref_genes = []
    ref_annotation = []
    for f in sorted(gb.features, key=lambda feature: feature.location.start):
        if f.type in ('CDS', 'tRNA', 'rRNA') and f.location.strand is not None:
            name = f.qualifiers.get('gene', f.qualifiers.get('product', ['?']))[0].lower()
            strand = '+' if f.location.strand > 0 else '-'
            ref_genes.append(name)
            ref_annotation.append((name, strand))
    print('参考基因数:', len(ref_genes))

    # 生成 4 种候选组合
    variants = {
        's1+s2': seq1 + seq2,
        's1+rc(s2)': seq1 + str(Seq(seq2).reverse_complement()),
        'rc(s1)+s2': str(Seq(seq1).reverse_complement()) + seq2,
        'rc(s1)+rc(s2)': str(Seq(seq1).reverse_complement()) + str(Seq(seq2).reverse_complement()),
    }

    print('\n=== 候选组合验证 ===')
    ranked = []
    for name, seq in variants.items():
        tmp = os.path.join(outdir, 'cand_%s.fasta' % name.replace('(','').replace(')','').replace('+','_'))
        with open(tmp, 'w') as fh:
            fh.write('>candidate\n%s\n' % seq)
        # minimap2 对参考 (检查连续性)
        ref_fa = os.path.join(outdir, 'ref.fna')
        SeqIO.write(gb, ref_fa, 'fasta')
        paf_call = subprocess.run(['minimap2', '-x', 'asm5', ref_fa, tmp],
                                  capture_output=True, text=True, check=True)
        paf = paf_call.stdout
        blocks = []
        for line in paf.split('\n'):
            if not line.strip(): continue
            p = line.split('\t')
            blocks.append((int(p[2]), int(p[3]), p[4], int(p[7]), int(p[8])))
        total_cov = sum(b[1]-b[0] for b in blocks)
        fwd_blocks = sum(1 for b in blocks if b[2] == '+')
        score = total_cov * (1 if fwd_blocks == len(blocks) and fwd_blocks > 0 else 0.5)
        print('%-15s 比对块=%d 正向块=%d 覆盖=%dbp 分数=%.0f' % (name, len(blocks), fwd_blocks, total_cov, score))
        ranked.append((score, name, seq, blocks))

    best = select_unique_candidate(ranked)
    if best is None and ranked:
        print('\n⚠ 最高分候选不唯一，拒绝凭首个候选猜测方向')
        sys.exit(1)
    if best is None or best[0] == 0:
        print('\n⚠ 所有组合都不匹配参考 — 可能需要更多 scaffold 或重新组装')
        sys.exit(1)

    score, name, seq, blocks = best
    print('\n✓ 最佳组合: %s (分数 %.0f)' % (name, score))
    print('  比对块: %s' % blocks[:5])

    # 验证: 逐基因 blast
    print('\n=== 基因顺序验证 (blast_genes.py 逻辑) ===')
    validation_dir = os.path.join(outdir, 'validation')
    os.makedirs(validation_dir, exist_ok=True)
    validation_fa = os.path.join(validation_dir, 'candidate.fasta')
    with open(validation_fa, 'w') as fh:
        fh.write('>mitogenome_circular\n%s\n' % seq)
    gene_order_path = os.path.join(validation_dir, 'gene_order.txt')
    subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), 'blast_genes.py'),
                    ref_gb, validation_fa, '--out', gene_order_path], check=True)
    query_annotation = []
    with open(gene_order_path, encoding='utf-8') as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 5:
                query_annotation.append((fields[0], fields[4]))
    if not circular_annotation_matches(query_annotation, ref_annotation):
        print('候选基因顺序或方向与参考不一致，拒绝输出', file=sys.stderr)
        sys.exit(1)
    if not reads_validated or not bam or not junction_region:
        print('缺少 reads 接缝证据；需要 --bam、--junction-region 和 --reads-validated', file=sys.stderr)
        sys.exit(2)
    try:
        has_spanning_read = junction_has_spanning_read(bam, junction_region)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print('接缝证据验证失败:', exc, file=sys.stderr)
        sys.exit(1)
    if not has_spanning_read:
        print('BAM 中没有跨越指定接缝的 read，拒绝输出', file=sys.stderr)
        sys.exit(1)

    with open(out_fa, 'w') as fh:
        fh.write('>mitogenome_circular\n%s\n' % seq)
    print('输出: %s (%d bp)' % (out_fa, len(seq)))

    print('\n下一步: 用 reads 回贴验证覆盖度 (depth_analysis.py), 然后注释 (SKILL.md ⑦)')

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
两 scaffold 候选拼接：参考辅助排序，样本 reads 检查内部接缝；不证明物理闭环
用法: python3 circularize.py <scaffold1.fasta> <scaffold2.fasta> <ref.gb> <outdir> --bam <candidate.bam> --reads-validated [--min-junction-support 3] [--min-mapq 20]
原理: 分析两个 scaffold 的基因组成和方向, 尝试 4 种组合(正序/反向互补),
      选择参考辅助排序的候选，自动计算新增接缝并验证；输出仍是候选结构。
      候选基因顺序/方向与参考不一致时记 REVIEW(退出码 2)并保留候选,
      不作为自动失败条件, 也不允许据此接受为最终环化。
依赖: BioPython, blastn, minimap2
"""
import sys, os, subprocess, re
import argparse



# 前置门禁：本步骤所需依赖（唯一清单来源 config/dependencies.json）
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from _deps import require_stage  # noqa: E402
from _sam import cigar_operations, aligned_blocks, template_key
from _evidence_io import read_single_fasta, protect_output, InputArgumentParser


def join_anchor_region(left_length, overlap, candidate_length, flank):
    """Span the entire shared overlap and flank bases exclusive to each scaffold."""
    lo, hi = left_length - overlap - flank + 1, left_length + flank
    if flank < 1 or lo < 1 or hi > candidate_length:
        raise ValueError('候选没有足够的两侧锚点，不能用单区间验证该拼接')
    return 'mitogenome_candidate:%d-%d' % (lo, hi)

def select_unique_candidate(candidates):
    if not candidates:
        return None
    top_score = max(score for score, *_ in candidates)
    top = [candidate for candidate in candidates if candidate[0] == top_score]
    return top[0] if len(top) == 1 else None


def interval_union_length(intervals):
    """Return the union length of half-open intervals (never double-count HSPs)."""
    merged = []
    for start, end in sorted((min(a, b), max(a, b)) for a, b in intervals if b > a):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return sum(end - start for start, end in merged)


def suffix_prefix_overlap(left, right, minimum=20):
    for size in range(min(len(left), len(right)), minimum - 1, -1):
        if left[-size:] == right[:size]:
            return size
    return 0


def join_scaffolds(left, right, minimum_overlap=20):
    overlap = suffix_prefix_overlap(left, right, minimum_overlap)
    return left + right[overlap:], overlap


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
    for length, operation in cigar_operations(fields[5]):
        if operation in 'M=XDN':
            end += length
    return start, end


def parse_region(region):
    try:
        chrom, coordinates = region.rsplit(':', 1); lo_text, hi_text = coordinates.split('-', 1)
        lo, hi = int(lo_text), int(hi_text)
    except (ValueError, AttributeError): raise ValueError('junction region must be chr:start-end')
    if lo < 1 or hi < lo: raise ValueError('invalid junction coordinates')
    return chrom, lo, hi


def junction_evidence(bam, region, min_mapq=20):
    """Read-level evidence for one junction: molecules, MAPQ, strand, duplicates.

    Duplicate-flagged reads (0x400) are excluded, but be aware that this only
    removes reads already marked by ``samtools markdup``/``fixmate``; running
    markdup on the candidate BAM is a prerequisite for molecule-level counting.
    """
    chrom, lo, hi = parse_region(region)
    result = subprocess.run(['samtools', 'view', bam, '%s:%d-%d' % (chrom, lo, hi)],
                            capture_output=True, text=True, check=True)
    templates = {}
    strand_counts = {'+': 0, '-': 0}
    read_strand_counts = {'+': 0, '-': 0}
    mapqs = []
    duplicates = 0
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 11:
            raise ValueError('SAM record has fewer than 11 fields')
        flag = int(fields[1])
        if flag & 0x400:
            duplicates += 1
            continue
        mapq = int(fields[4])
        if flag & (0x4 | 0x100 | 0x200 | 0x800) or mapq == 255 or mapq < min_mapq:
            continue
        if fields[2] != chrom:
            raise ValueError('samtools returned a different reference')
        # hi is inclusive; block end is exclusive. A deletion, skip or insertion
        # inside the anchor interval cannot support this exact candidate adjacency.
        if any(start <= lo and end > hi for start, end in aligned_blocks(fields)):
            key = template_key(fields)
            read_strand_counts['-' if flag & 0x10 else '+'] += 1
            rank = (not bool(flag & 0x40), -mapq, flag)
            if key not in templates or rank < templates[key][0]:
                templates[key] = (rank, flag, mapq)
    for _, flag, mapq in templates.values():
        mapqs.append(mapq)
        strand_counts['-' if flag & 0x10 else '+'] += 1
    mapqs.sort()
    return {'region': region, 'support': len(templates), 'names': sorted({key[1] for key in templates}),
            'templates': [dict(read_group=rg, qname=q) for rg, q in sorted(templates)],
            'read_strand_counts': read_strand_counts, 'duplicate_marking': 'not_verified',
            'strand_counts': strand_counts, 'duplicates_excluded': duplicates,
            'mapq_min': mapqs[0] if mapqs else None,
            'mapq_median': mapqs[len(mapqs) // 2] if mapqs else None}


def junction_spanning_reads(bam, region, min_mapq=20):
    return set(junction_evidence(bam, region, min_mapq=min_mapq)['names'])


def main():
    parser = InputArgumentParser(description=__doc__)
    for name in ('scaffold1', 'scaffold2', 'reference', 'outdir'):
        parser.add_argument(name)
    parser.add_argument('--bam')
    parser.add_argument('--reads-validated', action='store_true')
    parser.add_argument('--accept-candidate', action='store_true')
    parser.add_argument('--junction-region')
    parser.add_argument('--min-junction-support', type=int, default=3)
    parser.add_argument('--min-mapq', type=int, default=20)
    parser.add_argument('--junction-flank', type=int, default=10)
    args = parser.parse_args()
    if args.min_junction_support < 1 or args.junction_flank < 1 or not 0 <= args.min_mapq <= 254:
        parser.error('support/flank 必须为正数，MAPQ 为 0–254')
    require_stage('circularize', __file__)
    from Bio import SeqIO
    from Bio.Seq import Seq
    s1, s2, ref_gb, outdir = args.scaffold1, args.scaffold2, args.reference, args.outdir
    reads_validated, bam = args.reads_validated, args.bam
    junction_region, min_support = args.junction_region, args.min_junction_support
    min_mapq, junction_flank = args.min_mapq, args.junction_flank
    os.makedirs(outdir, exist_ok=True)
    out_fa = os.path.join(outdir, 'genome_candidate.fasta')
    accept_candidate = args.accept_candidate
    inputs = [s1, s2, ref_gb] + ([bam] if bam else [])
    generated = ['genome_candidate.fasta', 'ref.fna', 'validation/candidate.fasta', 'validation/gene_order.txt',
                 'cand_s1_s2.fasta', 'cand_s1_rcs2.fasta', 'cand_rcs1_s2.fasta', 'cand_rcs1_rcs2.fasta']
    for filename in generated:
        protect_output(os.path.join(outdir, filename), inputs)
    id1, seq1 = read_single_fasta(s1)
    id2, seq2 = read_single_fasta(s2)
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
    variants = {}
    for name, left, right in (
        ('s1+s2', seq1, seq2),
        ('s1+rc(s2)', seq1, str(Seq(seq2).reverse_complement())),
        ('rc(s1)+s2', str(Seq(seq1).reverse_complement()), seq2),
        ('rc(s1)+rc(s2)', str(Seq(seq1).reverse_complement()), str(Seq(seq2).reverse_complement())),
    ):
        joined, overlap = join_scaffolds(left, right)
        variants[name] = (joined, overlap, len(left) - overlap)

    print('\n=== 候选组合验证 ===')
    ranked = []
    for name, (seq, overlap, junction) in variants.items():
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
        total_cov = interval_union_length([(b[0], b[1]) for b in blocks])
        fwd_blocks = sum(1 for b in blocks if b[2] == '+')
        score = total_cov * (1 if fwd_blocks == len(blocks) and fwd_blocks > 0 else 0.5)
        print('%-15s 端部重叠=%dbp 比对块=%d 正向块=%d 覆盖=%dbp 分数=%.0f' % (name, overlap, len(blocks), fwd_blocks, total_cov, score))
        ranked.append((score, name, seq, blocks, overlap, junction))

    best = select_unique_candidate(ranked)
    if best is None and ranked:
        print('\n⚠ 最高分候选不唯一，拒绝凭首个候选猜测方向')
        sys.exit(1)
    if best is None or best[0] == 0:
        print('\n⚠ 所有组合都不匹配参考 — 可能需要更多 scaffold 或重新组装')
        sys.exit(1)

    score, name, seq, blocks, _overlap, junction = best
    print('\n✓ 最佳组合: %s (分数 %.0f)' % (name, score))
    print('  比对块: %s' % blocks[:5])

    # 验证: 逐基因 blast
    print('\n=== 基因顺序验证 (blast_genes.py 逻辑) ===')
    validation_dir = os.path.join(outdir, 'validation')
    os.makedirs(validation_dir, exist_ok=True)
    validation_fa = os.path.join(validation_dir, 'candidate.fasta')
    with open(validation_fa, 'w') as fh:
        fh.write('>mitogenome_candidate\n%s\n' % seq)
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
        order_review = True
        print('候选基因顺序或方向与参考不一致: 记 REVIEW (保留候选与诊断记录)', file=sys.stderr)
        print('  顺序差异不是自动失败条件; 需人工核对参考亲缘度、反向块与接缝证据', file=sys.stderr)
    else:
        order_review = False
    # 候选先落盘: 顺序不一致属 REVIEW, 不能因此丢弃候选或当成拼接失败
    if os.path.exists(out_fa):
        existing = str(next(SeqIO.parse(out_fa, 'fasta')).seq)
        if existing != seq:
            print('已有候选与本次选出的结构不同，拒绝复用旧 BAM；请使用新输出目录', file=sys.stderr)
            sys.exit(2)
    else:
        with open(out_fa, 'w') as fh:
            fh.write('>mitogenome_candidate\n%s\n' % seq)
        print('输出候选结构: %s (%d bp)' % (out_fa, len(seq)))
    auto_region = join_anchor_region(len(seq1), _overlap, len(seq), junction_flank)
    if junction_region and parse_region(junction_region) != parse_region(auto_region):
        print('指定 --junction-region 与实际新增接缝不一致；期望 %s' % auto_region, file=sys.stderr)
        sys.exit(2)
    if not reads_validated or not bam:
        print('缺少 reads 接缝证据；需要候选序列回贴 BAM 和 --reads-validated；状态=UNRESOLVED', file=sys.stderr)
        sys.exit(2)
    try:
        evidence = junction_evidence(bam, auto_region, min_mapq=min_mapq)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print('接缝证据验证失败:', exc, file=sys.stderr)
        sys.exit(1)
    print('接缝证据 %s: 支持模板=%d  模板链向=%s  MAPQ(min/median)=%s/%s  重复标记剔除=%d'
          % (auto_region, evidence['support'], evidence['strand_counts'],
             evidence['mapq_min'], evidence['mapq_median'], evidence['duplicates_excluded']))
    if evidence['support'] < min_support:
        print('接缝 %s 仅有 %d 个模板达到 MAPQ>=%d，需要至少 %d 个'
              % (auto_region, evidence['support'], min_mapq, min_support), file=sys.stderr)
        print('  证据不足以宣称闭环; 请先对候选 BAM 跑 markdup, 并检查重复/低 MAPQ 读段'
              '与竞争结构', file=sys.stderr)
        sys.exit(1)

    if not accept_candidate:
        print('状态: PUTATIVE_CIRCULAR/REVIEW；当前单一区间证据不能宣称最终环化。')
        print('若已人工核对全部新增接缝、重复歧义和组装图，再显式使用 --accept-candidate。')
        sys.exit(2)
    if order_review:
        print('候选基因顺序/方向与参考不一致，不能接受为最终环化; 请人工核验后再决定', file=sys.stderr)
        sys.exit(2)
    print('状态: CANDIDATE_ACCEPTED（仍需独立注释与完整 provenance）')

    print('\n下一步: 核验所有连接、重复歧义和尾首闭合，并独立检查注释；本工具未验证完整闭环。')

if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print('ERROR: 候选检查未完成: %s' % (getattr(exc, 'stderr', None) or str(exc)), file=sys.stderr)
        sys.exit(1)

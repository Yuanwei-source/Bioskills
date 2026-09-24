#!/usr/bin/env python3
"""
逐基因定位: 从参考 GB 提取 37 基因, 逐个 blastn-short 到目标序列, 重建基因顺序
用法: python3 blast_genes.py <ref.gb> <target.fna> [--out out.txt] [--evalue 1e-5]
依赖: BioPython, blastn
"""
import sys, os, subprocess, tempfile


def select_gene_hits(hits, min_identity=80.0, min_coverage=0.8, expected_names=None):
    selected = {}
    errors = []
    names = set(expected_names or hits)
    for name in sorted(names):
        candidates = [candidate for candidate in hits.get(name, [])
                      if candidate['identity'] >= min_identity and candidate['aln_len'] / candidate['query_len'] >= min_coverage]
        if len(candidates) == 1:
            selected[name] = candidates[0]
        elif not candidates:
            errors.append('%s 未找到满足 identity>=%.1f%% 且 coverage>=%.0f%% 的唯一命中' % (name, min_identity, min_coverage * 100))
        else:
            errors.append('%s 存在多个满足阈值的命中，不能判定唯一定位' % name)
    return selected, errors

def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    ref_gb, target = sys.argv[1], sys.argv[2]
    out = 'gene_order.txt'
    evalue = '1e-5'
    min_identity = 80.0
    min_coverage = 0.8
    if '--out' in sys.argv:
        out = sys.argv[sys.argv.index('--out')+1]
    if '--evalue' in sys.argv:
        evalue = sys.argv[sys.argv.index('--evalue')+1]
    if '--min-identity' in sys.argv:
        min_identity = float(sys.argv[sys.argv.index('--min-identity')+1])
    if '--min-coverage' in sys.argv:
        min_coverage = float(sys.argv[sys.argv.index('--min-coverage')+1])

    try:
        from Bio import SeqIO
    except ImportError:
        print('需要 BioPython: pip install biopython'); sys.exit(1)

    # 1. 提取参考基因
    gb = SeqIO.read(ref_gb, 'genbank')
    genes = []
    for f in gb.features:
        if f.type in ('CDS', 'tRNA', 'rRNA') and f.location.strand is not None:
            name = f.qualifiers.get('gene', f.qualifiers.get('product', ['?']))[0]
            genes.append((name, f.type, f.extract(gb.seq)))

    # 2. 建目标库
    tmpdir = tempfile.mkdtemp(prefix='blastgenes_')
    try:
        db = os.path.join(tmpdir, 'target')
        subprocess.run(['makeblastdb', '-in', target, '-dbtype', 'nucl', '-out', db, '-parse_seqids'],
                       check=True, capture_output=True)
        # 3. 逐个 blast (短基因用 blastn-short)
        hits = {}
        for name, ftype, seq in genes:
            q = os.path.join(tmpdir, 'q.fna')
            with open(q, 'w') as fh:
                fh.write('>query\n%s\n' % seq)
            r = subprocess.run(['blastn', '-query', q, '-db', db, '-task', 'blastn-short',
                                '-outfmt', '6 qseqid pident length sstart send evalue',
                                '-evalue', evalue, '-word_size', '7', '-gapopen', '5', '-gapextend', '2',
                                '-max_target_seqs', '2'], capture_output=True, text=True, check=True)
            candidates = hits.setdefault(name, [])
            for line in r.stdout.split('\n'):
                if not line.strip(): continue
                p = line.split('\t')
                ss, se = int(p[3]), int(p[4])
                candidates.append({
                    'evalue': float(p[5]), 'identity': float(p[1]),
                    'aln_len': int(p[2]), 'query_len': len(seq),
                    'start': min(ss, se), 'end': max(ss, se),
                    'strand': '+' if ss < se else '-', 'type': ftype,
                })
    finally:
        subprocess.run(['rm', '-rf', tmpdir])

    selected, errors = select_gene_hits(hits, min_identity, min_coverage, [gene[0] for gene in genes])
    if errors:
        for error in errors:
            print('ERROR: %s' % error)
        sys.exit(1)

    print('%-14s %-12s %8s %6s %6s %s' % ('基因', '类型', '位置', 'strand', 'pid%', 'cov%'))
    print('-' * 60)
    ordered = sorted(selected.items(), key=lambda item: item[1]['start'])
    for name, hit in ordered:
        coverage = hit['aln_len'] / hit['query_len'] * 100
        print('%-14s %-12s %6d-%6d  %s   %5.1f  %5.1f%%' % (name[:14], hit['type'], hit['start'], hit['end'], hit['strand'], hit['identity'], coverage))

    with open(out, 'w') as fh:
        for name, hit in ordered:
            coverage = hit['aln_len'] / hit['query_len'] * 100
            fh.write('%s\t%s\t%d\t%d\t%s\t%.1f\t%.1f\n' % (name, hit['type'], hit['start'], hit['end'], hit['strand'], hit['identity'], coverage))
    print('\n结果写入: %s (%d/%d 基因命中)' % (out, len(selected), len(genes)))

    # 提示: 与标准顺序对照
    print('\n→ 对照 references/standard_gene_order.md 检查基因顺序与方向')

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
逐基因定位: 从参考 GB 提取 37 基因, 逐个 blastn-short 到目标序列, 重建基因顺序
用法: python3 blast_genes.py <ref.gb> <target.fna> [--out out.txt] [--evalue 1e-5]
依赖: BioPython, blastn
"""
import sys, os, subprocess, tempfile
import argparse
import math
import shutil
from collections import Counter



# 前置门禁：本步骤所需依赖（唯一清单来源 config/dependencies.json）
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from _deps import require_stage  # noqa: E402
from _evidence_io import protect_output, atomic_write_text, InputArgumentParser

def hsps_follow(previous, current):
    """Conservative local chain; gap limits are engineering screens, not gene laws."""
    if previous['strand'] != current['strand'] or previous.get('target') != current.get('target'):
        return False
    qgap = current['qstart'] - previous['qend'] - 1
    sgap = (current['start'] - previous['end'] - 1 if current['strand'] == '+'
            else previous['start'] - current['end'] - 1)
    max_gap = max(20, int(previous['query_len'] * 0.2))
    return 0 <= qgap <= max_gap and 0 <= sgap <= max_gap and abs(qgap-sgap) <= max_gap


def merge_hsps(hsps):
    """Merge one unambiguous chain; never reuse query or target bases."""
    if not hsps:
        return None
    ordered = sorted(hsps, key=lambda h: (h['qstart'], h['qend']))
    if any(not hsps_follow(a, b) for a, b in zip(ordered, ordered[1:])):
        return None
    covered = sum(h['qend'] - h['qstart'] + 1 for h in ordered)
    aligned = sum(h['aln_len'] for h in ordered)
    identity = sum(h['identity'] * h['aln_len'] for h in ordered) / aligned
    return dict(ordered[0], identity=identity, aln_len=covered, alignment_columns=aligned,
                evalue=max(h['evalue'] for h in ordered),
                qstart=ordered[0]['qstart'], qend=ordered[-1]['qend'],
                start=min(h['start'] for h in ordered), end=max(h['end'] for h in ordered),
                hsps=ordered)


def locus_candidates(hsps):
    """Chain only mutually unique neighbours; retain alternative loci/branches."""
    unique = {}
    for hit in hsps:
        key = (hit['qstart'], hit['qend'], hit['start'], hit['end'], hit['strand'], hit.get('target'))
        old = unique.get(key)
        if old is None or hit['evalue'] < old['evalue']:
            unique[key] = hit
    ordered = sorted(unique.values(), key=lambda h: (h['qstart'], h['start']))
    successors = {i: [j for j in range(len(ordered)) if hsps_follow(hit, ordered[j])]
                  for i, hit in enumerate(ordered)}
    predecessors = {j: [i for i in successors if j in successors[i]] for j in successors}
    if any(len(v) > 1 for v in list(successors.values()) + list(predecessors.values())):
        return [dict(hit, ambiguous_chain=True) for hit in ordered]
    result = []
    for start in predecessors:
        if predecessors[start]:
            continue
        chain, index = [], start
        while True:
            chain.append(ordered[index])
            if not successors[index]:
                break
            index = successors[index][0]
        result.append(merge_hsps(chain))
    return result


def select_gene_hits(hits, min_identity=80.0, min_coverage=0.8, expected_names=None):
    selected = {}
    errors = []
    names = set(expected_names or hits)
    for name in sorted(names):
        if any(hit.get('ambiguous_chain') for hit in hits.get(name, [])):
            errors.append('%s HSP 链存在分支，无法唯一汇总' % name)
            continue
        candidates = [candidate for candidate in hits.get(name, [])
                      if candidate['identity'] >= min_identity and candidate['aln_len'] / candidate['query_len'] >= min_coverage]
        if len(candidates) == 1:
            selected[name] = candidates[0]
        elif not candidates:
            errors.append('%s 未找到满足 identity>=%.1f%% 且 coverage>=%.0f%% 的唯一命中' % (name, min_identity, min_coverage * 100))
        else:
            errors.append('%s 存在多个满足阈值的命中，不能判定唯一定位' % name)
    return selected, errors


def reference_gene_labels(features):
    """Keep repeated feature names distinct without assigning biological copy IDs."""
    names = [feature.qualifiers.get('gene', feature.qualifiers.get('product', ['?']))[0]
             for feature in features]
    counts = Counter(names)
    labels = []
    for name, feature in zip(names, features):
        if counts[name] == 1:
            labels.append(name)
            continue
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand = '+' if feature.location.strand == 1 else '-'
        labels.append('%s@%d..%d:%s' % (name, start, end, strand))
    if len(labels) != len(set(labels)):
        raise ValueError('参考中重复 feature 的名称与坐标均相同，无法生成唯一追踪标签')
    return labels

def main():
    parser = InputArgumentParser(description=__doc__)
    parser.add_argument('reference')
    parser.add_argument('target')
    parser.add_argument('--out', default='gene_order.txt')
    parser.add_argument('--evalue', type=float, default=1e-5)
    parser.add_argument('--min-identity', type=float, default=80)
    parser.add_argument('--min-coverage', type=float, default=0.8)
    args = parser.parse_args()
    if not 0 <= args.min_identity <= 100 or not 0 < args.min_coverage <= 1 or not math.isfinite(args.evalue) or args.evalue <= 0:
        parser.error('identity 必须为 0–100，coverage 为 (0,1]，evalue 为有限正数')
    require_stage('gene_locating', __file__)
    from Bio import SeqIO
    ref_gb, target, out = args.reference, args.target, args.out
    evalue, min_identity, min_coverage = str(args.evalue), args.min_identity, args.min_coverage
    protect_output(out, [ref_gb, target])

    # 1. 提取参考基因
    gb = SeqIO.read(ref_gb, 'genbank')
    features = []
    for f in gb.features:
        if f.type in ('CDS', 'tRNA', 'rRNA') and f.location.strand is not None:
            features.append(f)
    names = reference_gene_labels(features)
    genes = [(name, feature.type, feature.extract(gb.seq))
             for name, feature in zip(names, features)]
    if not genes:
        raise ValueError('参考中没有可定位的 CDS/tRNA/rRNA')

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
                                '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send evalue',
                                '-evalue', evalue, '-word_size', '7', '-gapopen', '5', '-gapextend', '2',
                                '-max_target_seqs', '2'], capture_output=True, text=True, check=True)
            candidates = hits.setdefault(name, [])
            raw = {}
            for line in r.stdout.split('\n'):
                if not line.strip(): continue
                p = line.split('\t')
                ss, se = int(p[6]), int(p[7])
                key = (p[1], '+' if ss < se else '-')
                raw.setdefault(key, []).append({
                    'evalue': float(p[8]), 'identity': float(p[2]),
                    'aln_len': int(p[3]), 'query_len': len(seq),
                    'qstart': int(p[4]), 'qend': int(p[5]),
                    'start': min(ss, se), 'end': max(ss, se),
                    'strand': '+' if ss < se else '-', 'type': ftype, 'target': p[1],
                })
            for group in raw.values():
                candidates.extend(locus_candidates(group))
    finally:
        shutil.rmtree(tmpdir)

    selected, errors = select_gene_hits(hits, min_identity, min_coverage, [gene[0] for gene in genes])
    if errors:
        for error in errors:
            print('ERROR: %s' % error)
        sys.exit(1)

    print('%-30s %-12s %8s %6s %6s %s' % ('参考特征', '类型', '位置', 'strand', 'pid%', 'cov%'))
    print('-' * 78)
    ordered = sorted(selected.items(), key=lambda item: item[1]['start'])
    for name, hit in ordered:
        coverage = hit['aln_len'] / hit['query_len'] * 100
        print('%-30s %-12s %6d-%6d  %s   %5.1f  %5.1f%%' % (name[:30], hit['type'], hit['start'], hit['end'], hit['strand'], hit['identity'], coverage))

    lines = []
    for name, hit in ordered:
        coverage = hit['aln_len'] / hit['query_len'] * 100
        lines.append('%s\t%s\t%d\t%d\t%s\t%.1f\t%.1f\t%s\n' %
                     (name, hit['type'], hit['start'], hit['end'], hit['strand'], hit['identity'], coverage, hit['target']))
    atomic_write_text(out, ''.join(lines))
    print('\n结果写入: %s (%d/%d 基因命中)' % (out, len(selected), len(genes)))

    # 提示: 与标准顺序对照
    print('\n→ 对照 references/standard_gene_order.md 检查基因顺序与方向')

if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print('ERROR: 基因定位未完成: %s' % (getattr(exc, 'stderr', None) or str(exc)), file=sys.stderr)
        sys.exit(1)

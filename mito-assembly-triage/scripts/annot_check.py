#!/usr/bin/env python3
"""
注释质检脚本: 对注释好的线粒体基因组 GB 文件做系统性检查

检查项:
  1. 基因完整性: 13 CDS + 22 tRNA + 2 rRNA = 37 基因
  2. CDS 翻译: 无内部终止密码子 (遗传密码表 5, 默认; 可 --table 指定)
  3. CDS 起始终止密码子报告 (含不完整终止 T/TA 标注)
  4. tRNA: 数量 22, 长度 60-75bp (异常短 <50bp 报错), 反密码子, 链
  5. rRNA: 数量 2, 长度范围 (rrnL 1100-1500, rrnS 600-850)
  6. 重叠检查: CDS-tRNA 重叠 >8bp 报错 (NCBI: 只能重叠几个bp)
              其他类型重叠 >8bp 警告 (ND4-ND4L 等短重叠正常)
  7. 链分布: 正链 CDS 数量 (昆虫标准 9+/4-)
  8. 基因顺序 (可选 --ref ref.gb 对比近缘参考)
  9. CDS 长度与参考对比 (可选 --ref, ±20% 警告)

用法:
  python3 annot_check.py <annotation.gb>
  python3 annot_check.py <annotation.gb> --ref <reference.gb> --table 5
  python3 annot_check.py <annotation.gb> --tolerate-overlap "ND5,trnH"   # 人工确认的真实重叠

--tolerate-overlap: 人工确认的 CDS-tRNA 重叠对 (多工具一致 + 证据充分), 从 ERROR 降为 ACCEPTED。
  可重复使用, 例如 --tolerate-overlap "ND5,trnH" --tolerate-overlap "ND4,trnH2"

退出码: 0=全部通过 1=有错误 2=有警告但无错误
"""
import sys
import re
from collections import Counter

EXPECTED_CDS = {'atp6', 'atp8', 'cox1', 'cox2', 'cox3', 'cytb', 'nd1', 'nd2', 'nd3', 'nd4', 'nd4l', 'nd5', 'nd6'}
EXPECTED_TRNA = {'trna', 'trnc', 'trnd', 'trne', 'trnf', 'trng', 'trnh', 'trni', 'trnk', 'trnl1', 'trnl2', 'trnm', 'trnn', 'trnp', 'trnq', 'trnr', 'trns1', 'trns2', 'trnt', 'trnv', 'trnw', 'trny'}
EXPECTED_RRNA = {'rrnl', 'rrns'}


def feature_gene(feature):
    values = feature.qualifiers.get('gene') or feature.qualifiers.get('product') or ['?']
    return str(values[0])


def canonical_gene(feature):
    key = re.sub(r'[^a-z0-9]', '', feature_gene(feature).lower())
    if key in ('16s', 'srrna'):
        return 'rrns'
    if key in ('12s', 'lrrna', 'lrna'):
        return 'rrnl'
    if key == 'trnl':
        return 'trnl1'
    if key == 'trns':
        return 'trns1'
    return key


def gene_identity_errors(cds, trnas, rnas):
    errors = []
    for label, features, expected in (('CDS', cds, EXPECTED_CDS), ('tRNA', trnas, EXPECTED_TRNA), ('rRNA', rnas, EXPECTED_RRNA)):
        names = [canonical_gene(feature) for feature in features]
        counts = Counter(names)
        missing = sorted(expected - set(names))
        extra = sorted(set(names) - expected)
        duplicates = sorted(name for name, count in counts.items() if count > 1)
        if missing:
            errors.append('%s 缺少基因身份: %s' % (label, ', '.join(missing)))
        if extra:
            errors.append('%s 出现非标准/未知基因身份: %s' % (label, ', '.join(extra)))
        if duplicates:
            errors.append('%s 基因身份重复: %s' % (label, ', '.join(duplicates)))
    return errors

def load_gb(fn):
    try:
        from Bio import SeqIO
        return SeqIO.read(fn, 'genbank')
    except Exception as e:
        print('ERROR: 无法读取 GB 文件 %s: %s' % (fn, e))
        sys.exit(1)

def main():
    args = sys.argv[1:]
    if not args or args[0] in ('-h', '--help'):
        print(__doc__); sys.exit(0)
    fn = args[0]
    table = 5
    ref_fn = None
    require_circular = '--require-circular' in args
    if '--table' in args:
        table = int(args[args.index('--table') + 1])
    if '--ref' in args:
        ref_fn = args[args.index('--ref') + 1]
    tolerate = []
    for i, a in enumerate(args):
        if a == '--tolerate-overlap':
            pair = tuple(x.strip() for x in args[i + 1].split(','))
            tolerate.append(pair)

    gb = load_gb(fn)
    L = len(gb.seq)
    print('=' * 68)
    print('线粒体基因组注释质检: %s' % fn)
    topology = gb.annotations.get('topology', '')
    print('基因组: %d bp  |  遗传密码表: %d  |  拓扑: %s' % (L, table, topology or '?'))
    print('=' * 68)

    errors, warnings = [], []
    if require_circular and topology.lower() != 'circular':
        errors.append('要求环状拓扑，但文件 topology=%s' % (topology or '未声明'))
    from Bio.Data import CodonTable
    tbl = CodonTable.unambiguous_dna_by_id[table]

    # ---------- 1. 基因完整性 ----------
    cds = [f for f in gb.features if f.type == 'CDS']
    trnas = [f for f in gb.features if f.type == 'tRNA']
    rnas = [f for f in gb.features if f.type == 'rRNA']
    print('\n[1] 基因完整性: CDS=%d tRNA=%d rRNA=%d (标准 13/22/2)' % (len(cds), len(trnas), len(rnas)))
    if len(cds) != 13: errors.append('CDS 数量 %d != 13' % len(cds))
    if len(trnas) != 22: errors.append('tRNA 数量 %d != 22 (NCBI 要求 22)' % len(trnas))
    if len(rnas) != 2: errors.append('rRNA 数量 %d != 2' % len(rnas))
    print('\n[1b] 基因身份:')
    identity_errors = gene_identity_errors(cds, trnas, rnas)
    if identity_errors:
        for error in identity_errors:
            errors.append(error)
            print('  ERROR: %s' % error)
    else:
        print('  ✓ 13 CDS + 22 tRNA + 2 rRNA 的基因身份完整且无重复')

    # ---------- 2. CDS 翻译 ----------
    print('\n[2] CDS 翻译验证 (密码表 %d):' % table)
    for f in sorted(cds, key=lambda x: x.location.start):
        g = f.qualifiers.get('gene', ['?'])[0]
        lo, hi = f.location.start + 1, f.location.end
        strand = '+' if f.location.strand > 0 else '-'
        cs = f.extract(gb.seq)
        aa = str(cs.translate(table=tbl))
        internal = aa[:-1].count('*')
        # 起始终止密码子
        start_codon = str(cs[:3]).upper()
        last3 = str(cs[-3:]).upper()
        valid_starts = set(tbl.start_codons)
        valid_stops = set(tbl.stop_codons)
        incomplete = last3 not in valid_stops
        if start_codon not in valid_starts:
            errors.append('%s [%d..%d] 起始密码子 %s 不符合遗传密码表' % (g, lo, hi, start_codon))
        if incomplete and not f.qualifiers.get('note'):
            errors.append('%s [%d..%d] 缺少完整终止密码子 %s' % (g, lo, hi, last3))
        elif incomplete:
            warnings.append('%s [%d..%d] 终止密码子不完整，已按 note 记录' % (g, lo, hi))
        status = 'OK' if internal == 0 else 'ERROR'
        if internal != 0: errors.append('%s [%d..%d] 有 %d 个内部终止密码子' % (g, lo, hi, internal))
        print('  %-6s [%5d..%5d] %s %3daa  起始=%s 终止=%s%s %s' % (
            g, lo, hi, strand, len(aa), start_codon, last3,
            ' (不完整)' if incomplete else '', status))

    # ---------- 3. tRNA ----------
    print('\n[3] tRNA 检查 (22 个, 长度 60-75bp):')
    short_trnas = []
    for f in sorted(trnas, key=lambda x: x.location.start):
        g = f.qualifiers.get('gene', ['?'])[0]
        lo, hi = f.location.start + 1, f.location.end
        ln = hi - lo + 1
        flag = ''
        if ln < 50:
            # 有人工确认 note 的短 tRNA (如 blast 保守区锚定的缩短型) -> WARN; 否则 ERROR
            has_note = 'note' in f.qualifiers
            if has_note:
                flag = ' WARN-过短(人工确认)'
                warnings.append('tRNA %s 长度 %dbp <50 (note 说明人工确认: %s)' % (g, ln, str(f.qualifiers['note'][0])[:60]))
            else:
                flag = ' ERROR-过短'
                errors.append('tRNA %s 长度 %dbp 异常短 (<50), 无人工确认 note' % (g, ln))
                short_trnas.append(g)
        elif ln < 60 or ln > 75:
            flag = ' WARN-长度偏离'
            warnings.append('tRNA %s 长度 %dbp 偏离 60-75' % (g, ln))
        # 反密码子 (有 anticodon qualifier 时)
        anti = f.qualifiers.get('anticodon', ['?'])[0]
        astr = (' 反密码子=%s' % anti) if anti != '?' else ''
        if anti == '?':
            warnings.append('tRNA %s 缺少反密码子 qualifiers' % g)
        print('  %-6s [%5d..%5d] %3dbp%s%s' % (g, lo, hi, ln, astr, flag))
    if short_trnas:
        print('  ! 过短 tRNA: %s (检查坐标来源, 应从 arwen/MITOS2 直接提取)' % ', '.join(short_trnas))

    # ---------- 4. rRNA ----------
    print('\n[4] rRNA 检查:')
    for f in rnas:
        g = f.qualifiers.get('gene', ['?'])[0]
        lo, hi = f.location.start + 1, f.location.end
        ln = hi - lo + 1
        rng = '1100-1500' if g.lower() in ('rrnl', '16s', 'lrna') else '600-850'
        if not (int(rng.split('-')[0]) <= ln <= int(rng.split('-')[1])):
            warnings.append('rRNA %s 长度 %d 超出范围 %s' % (g, ln, rng))
            print('  %-6s [%5d..%5d] %4dbp (期望 %s) WARN' % (g, lo, hi, ln, rng))
        else:
            print('  %-6s [%5d..%5d] %4dbp OK' % (g, lo, hi, ln))

    # ---------- 5. 重叠检查 ----------
    print('\n[5] 重叠检查 (>8bp):')
    feats = sorted([f for f in cds + trnas + rnas], key=lambda f: f.location.start)
    n_overlap = 0
    for i, f1 in enumerate(feats):
        for f2 in feats[i + 1:]:
            if f2.location.start >= f1.location.end:
                break
            ov = f1.location.end - f2.location.start
            if ov <= 8:
                continue
            g1 = f1.qualifiers.get('gene', ['?'])[0]
            g2 = f2.qualifiers.get('gene', ['?'])[0]
            t1, t2 = f1.type, f2.type
            n_overlap += 1
            if {t1, t2} == {'CDS', 'tRNA'}:
                if (g1, g2) in tolerate or (g2, g1) in tolerate:
                    print('  ACCEPT: %-6s[%s] <-> %-6s[%s] 重叠 %3dbp (人工确认的真实特征)' % (g1, t1, g2, t2, ov))
                else:
                    errors.append('%s(%s) 与 %s(%s) 重叠 %dbp (CDS-tRNA 只允许几个bp; 若为真实特征用 --tolerate-overlap "%s,%s" 确认)' % (g1, t1, g2, t2, ov, g1, g2))
                    print('  ERROR: %-6s[%s] <-> %-6s[%s] 重叠 %3dbp (CDS-tRNA)' % (g1, t1, g2, t2, ov))
            else:
                warnings.append('%s(%s) 与 %s(%s) 重叠 %dbp (需人工确认是否真实特征)' % (g1, t1, g2, t2, ov))
                print('  WARN : %-6s[%s] <-> %-6s[%s] 重叠 %3dbp (人工确认)' % (g1, t1, g2, t2, ov))
    if n_overlap == 0:
        print('  (无 >8bp 重叠)')

    # ---------- 6. 链分布 ----------
    plus = sum(1 for f in cds if f.location.strand > 0)
    minus = len(cds) - plus
    print('\n[6] 链分布: 正链 CDS=%d 负链 CDS=%d (昆虫标准 9+/4-)' % (plus, minus))
    if plus != 9:
        warnings.append('正链 CDS 数量 %d != 9 (昆虫标准, 需确认)' % plus)

    # ---------- 7. 基因顺序与 CDS 长度 (参考对比) ----------
    if ref_fn:
        ref = load_gb(ref_fn)
        print('\n[7] 参考对比 (%s):' % ref_fn)
        ref_cds = {f.qualifiers.get('gene', ['?'])[0]: f for f in ref.features if f.type == 'CDS'}
        # 基因顺序
        order_q = [f.qualifiers.get('gene', ['?'])[0] for f in cds]
        order_r = [f.qualifiers.get('gene', ['?'])[0] for f in sorted(ref_cds.values(), key=lambda x: x.location.start)]
        q_set = [g for g in order_q if g in order_r]
        print('  基因顺序匹配(CDS): %s' % ('一致' if q_set == order_r else '不同! 查询=%s' % '>'.join(q_set)))
        if q_set != order_r:
            warnings.append('CDS 基因顺序与参考不同: %s vs %s' % ('>'.join(q_set), '>'.join(order_r)))
        # 长度对比
        for f in sorted(cds, key=lambda x: x.location.start):
            g = f.qualifiers.get('gene', ['?'])[0]
            ln = f.location.end - f.location.start + 1
            rf = ref_cds.get(g)
            if rf:
                rln = rf.location.end - rf.location.start + 1
                d = abs(ln - rln) / rln * 100
                flag = ' WARN' if d > 20 else ' OK'
                if d > 20: warnings.append('%s 长度 %d vs 参考 %d (差 %.0f%%)' % (g, ln, rln, d))
                print('  %-6s %4dbp vs 参考 %4dbp (差 %4.0f%%)%s' % (g, ln, rln, d, flag))

    # ---------- 汇总 ----------
    print('\n' + '=' * 68)
    if errors:
        print('结果: 错误 %d 个 (必须修复), 警告 %d 个' % (len(errors), len(warnings)))
        for e in errors:
            print('  [ERROR] %s' % e)
        for w in warnings:
            print('  [WARN ] %s' % w)
        print('=' * 68)
        sys.exit(1)
    elif warnings:
        print('结果: 通过 (无错误), 警告 %d 个 (需人工确认)' % len(warnings))
        for w in warnings:
            print('  [WARN ] %s' % w)
        print('=' * 68)
        sys.exit(2)
    else:
        print('结果: 全部通过 (37 基因完整, CDS 干净, 无异常重叠)')
        print('=' * 68)
        sys.exit(0)

if __name__ == '__main__':
    main()

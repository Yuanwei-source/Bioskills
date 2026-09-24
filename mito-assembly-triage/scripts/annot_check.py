#!/usr/bin/env python3
"""
注释质检脚本: 对注释好的线粒体基因组 GB 文件做系统性检查

判据来源分三层, 不要混为一谈:
  (a) 一般生物学预期      —— 典型后生动物 13 CDS + 22 tRNA + 2 rRNA; 允许有证据的类群例外
  (b) NCBI 提交审查要求   —— 要求提供基因/CDS 等注释, 差异需说明; 未规定统一的 bp 阈值
  (c) 本工具的工程预警阈值 —— 如 >8bp 重叠、tRNA 60-75bp、rrnL/rrnS 长度区间、9+/4- 链分布
      (c) 是本工具的默认值, 不是领域公理; 用 --overlap-severity / --allow-atypical 等显式放宽。

检查项:
  1. 基因完整性: 13 CDS + 22 tRNA + 2 rRNA (典型后生动物预期, 见 (a))
     (CDS 同义命名归一: nad* == nd*, cob == cytb, coi/coii/coiii == cox1/cox2/cox3)
     裸名 trnL / trnS 不自动归一到 trnL1/trnS1, 记为 UNDETERMINED_TRNA
  2. CDS 翻译: 内部终止 / 起始密码子 / 终止密码子 (含真实 partial CDS 的 /codon_start)
     非典型起始密码子报 NONCANONICAL_START_REVIEW, 由 --tolerate-start 按基因+密码子显式确认
  3. tRNA: 数量, 长度 (60-75bp 预警; <50bp 且无 note 报错), 反密码子, 类型判定
  4. rRNA: 数量, 长度区间 (工程预警; 按归一后身份 rrns/rrnl 选择区间)
  5. 重叠: 逐对按分段坐标计算并全部记录; >8bp 默认 REVIEW (可 --overlap-severity error)
  6. 链分布: 正链 CDS 数量 (昆虫背景预警, 非跨类群标准); 整条序列反向互补时不误报
  7. 基因顺序 (可选 --ref, 按归一后基因身份比较)
  8. CDS 长度与参考对比 (可选 --ref, ±20% 预警, 按分段长度)

用法:
  python3 annot_check.py <annotation.gb>
  python3 annot_check.py <annotation.gb> --ref <reference.gb> --table 5
  python3 annot_check.py <annotation.gb> --tolerate-overlap "nad5,trnH"      # 已人工审核并保留该注释
  python3 annot_check.py <annotation.gb> --tolerate-start "cox1:CGA"         # 类群+基因+密码子限定的非典型起始
  python3 annot_check.py <annotation.gb> --overlap-severity error            # 把 >8bp 重叠升级为错误
  python3 annot_check.py <annotation.gb> --allow-atypical "<理由>"           # 基因集差异记待核查(需给理由)

--tolerate-overlap GENE1,GENE2
  表示该基因对的重叠"已人工审核并保留该注释", 仍会记录, 但不再声明其功能真实性。
  名称按归一后身份比较, 大小写不敏感; 未匹配到任何重叠对的条目会报警(避免静默失效)。
--tolerate-start GENE:CODON
  表示该基因在该密码子上的非典型起始"已有类群/文献/同源证据", 记为 NONCANONICAL_START_REVIEW。
  必须逐条给出, 不能全局放行; 未匹配到任何 CDS 的条目会报警。
--overlap-severity error|warn
  >8bp 重叠的默认等级, 默认 warn(=待核查 REVIEW)。任何重叠都不构成自动裁剪序列的理由。
--allow-atypical "<理由>"
  仅把「基因集数量与身份」差异从 ERROR 降为待核查, 必须给出理由并逐项确认;
  不放松起始密码子、重叠、tRNA/rRNA 长度与链分布检查。

退出码: 0=无任何发现 1=有错误(含参数错误) 2=仅有待核查/警告
注意: 退出码 0 只表示"注释内部一致", 不构成"序列被 reads 支持"的证据。
"""
import sys
import re
from collections import Counter

EXPECTED_CDS = {'atp6', 'atp8', 'cox1', 'cox2', 'cox3', 'cytb', 'nd1', 'nd2', 'nd3', 'nd4', 'nd4l', 'nd5', 'nd6'}
EXPECTED_TRNA = {'trna', 'trnc', 'trnd', 'trne', 'trnf', 'trng', 'trnh', 'trni', 'trnk', 'trnl1', 'trnl2', 'trnm', 'trnn', 'trnp', 'trnq', 'trnr', 'trns1', 'trns2', 'trnt', 'trnv', 'trnw', 'trny'}
EXPECTED_RRNA = {'rrnl', 'rrns'}

#: 裸名 -> 该裸名可能对应的多个身份。裸名不足以判定类型, 必须等反密码子/结构/同源证据。
AMBIGUOUS_TRNA = {'trnl': ('trnl1', 'trnl2'), 'trns': ('trns1', 'trns2')}
#: 反密码子 (DNA 写法) -> 归一身份。动物线粒体标准: trnL1(CUN)=tag, trnL2(UUR)=taa,
#: trnS1(AGN)=gct, trnS2(UCN)=tga
ANTICODON_IDENTITY = {'tag': 'trnl1', 'taa': 'trnl2', 'gct': 'trns1', 'tga': 'trns2'}

TRNA_LENGTH_RANGE = (60, 75)
TRNA_HARD_MIN = 50
RRNA_LENGTH_RANGE = {'rrnl': (1100, 1500), 'rrns': (600, 850)}
OVERLAP_REVIEW_THRESHOLD = 8
CDS_LENGTH_TOLERANCE_PCT = 20
EXPECTED_PLUS_STRAND_CDS = 9


# --------------------------------------------------------------------------- names
def feature_gene(feature):
    values = feature.qualifiers.get('gene') or feature.qualifiers.get('product') or ['?']
    return str(values[0])


def _canonical_key(text):
    """Normalise a free-text gene name to a comparison key (never guesses trnL/trnS)."""
    key = re.sub(r'\(.*?\)', '', str(text)).lower()      # drop "(UUR)" / "(AGN)"
    key = re.sub(r'[^a-z0-9]', '', key)
    if key in ('16s', '16srrna', 'lrrna', 'lrna', 'rrnl'):
        return 'rrnl'
    if key in ('12s', '12srrna', 'srrna', 'rrns'):
        return 'rrns'
    if re.fullmatch(r'nad[1-6]l?', key):                 # nad* (新写法) == nd*
        return 'nd' + key[3:]
    if key in ('cob', 'cytb', 'cytb2'):
        return 'cytb'
    if key in ('coi', 'cox1', 'cox1p'):
        return 'cox1'
    if key in ('coii', 'cox2'):
        return 'cox2'
    if key in ('coiii', 'cox3'):
        return 'cox3'
    return key                                            # trnl / trns stay bare


def canonical_gene(feature):
    return _canonical_key(feature_gene(feature))


def parse_anticodon(feature):
    raw = feature.qualifiers.get('anticodon')
    if not raw:
        return None
    match = re.search(r'seq\s*:\s*([acgtu]+)', str(raw[0]), re.I)
    if not match:
        return None
    return match.group(1).lower().replace('u', 't')


def resolve_trna_identity(feature):
    """Return a concrete tRNA identity, or None when the type is not determined."""
    key = canonical_gene(feature)
    if key in EXPECTED_TRNA:
        return key
    if key in AMBIGUOUS_TRNA:
        anticodon = parse_anticodon(feature)
        candidate = ANTICODON_IDENTITY.get(anticodon or '')
        if candidate in AMBIGUOUS_TRNA[key]:
            return candidate
    return None


# ------------------------------------------------------------------- coordinates
def feature_segments(feature):
    """Real base intervals (0-based half-open), one entry per location part.

    Compound (``join``) locations keep their parts, so a gene spanning the origin
    reports its true length instead of the min/max span across the whole genome.
    """
    return [(int(part.start), int(part.end)) for part in feature.location.parts]


def feature_length(feature):
    return sum(end - start for start, end in feature_segments(feature))


def segment_overlap(left, right):
    total = 0
    for start_a, end_a in left:
        for start_b, end_b in right:
            total += max(0, min(end_a, end_b) - max(start_a, start_b))
    return total


def feature_span(feature):
    segments = feature_segments(feature)
    return min(start for start, _ in segments) + 1, max(end for _, end in segments)


# ---------------------------------------------------------------------- findings
class Findings:
    """Three-level finding collector.  Errors block, reviews need accounting."""

    def __init__(self):
        self.errors = []
        self.reviews = []
        self.infos = []

    def error(self, message):
        self.errors.append(message)

    def review(self, message):
        self.reviews.append(message)

    def info(self, message):
        self.infos.append(message)


def identity_findings(findings, cds, trnas, rnas, atypical_reason):
    """Gene-set completeness and identity, with bare trnL/trnS left undetermined."""
    sink = findings.review if atypical_reason else findings.error
    tag = ('待核查(非典型类群: %s)' % atypical_reason) if atypical_reason else 'ERROR'

    for label, features, expected in (('CDS', cds, EXPECTED_CDS),
                                      ('rRNA', rnas, EXPECTED_RRNA)):
        names = [canonical_gene(feature) for feature in features]
        counts = Counter(names)
        missing = sorted(expected - set(names))
        extra = sorted(set(names) - expected)
        duplicates = sorted(name for name, count in counts.items() if count > 1)
        if missing:
            sink('%s %s 缺少基因身份: %s' % (tag, label, ', '.join(missing)))
        if extra:
            sink('%s %s 出现非标准/未知基因身份: %s' % (tag, label, ', '.join(extra)))
        if duplicates:
            sink('%s %s 基因身份重复: %s' % (tag, label, ', '.join(duplicates)))

    resolved, pending = [], []
    for feature in trnas:
        identity = resolve_trna_identity(feature)
        (resolved if identity else pending).append(identity or canonical_gene(feature))
    missing = set(EXPECTED_TRNA) - set(resolved)

    # A bare trnL/trnS cannot be assigned to L1/L2 or S1/S2 without anticodon,
    # structure or homology evidence, so the group's completeness is undetermined.
    for bare, group in AMBIGUOUS_TRNA.items():
        count = pending.count(bare)
        if not count:
            continue
        for member in group:
            missing.discard(member)
        findings.review(
            'UNDETERMINED_TRNA: %d 个 %s 未给出类型(裸名), 不得自动归一为 %s; '
            '无法判定 %s 是否缺失 —— 需反密码子 / 结构 / 同源证据确认'
            % (count, bare, '/'.join(group), '/'.join(group)))

    missing = sorted(missing)
    extra = sorted(set(resolved) - EXPECTED_TRNA)
    counts = Counter(resolved)
    duplicates = sorted(name for name, count in counts.items() if count > 1)
    if missing:
        sink('%s tRNA 缺少基因身份: %s' % (tag, ', '.join(missing)))
    if extra:
        sink('%s tRNA 出现非标准/未知基因身份: %s' % (tag, ', '.join(extra)))
    if duplicates:
        sink('%s tRNA 基因身份重复: %s' % (tag, ', '.join(duplicates)))
    if not (missing or extra or duplicates or pending):
        findings.info('tRNA 身份完整且无重复 (22 个)')


def cds_findings(findings, gb, cds, tbl, start_exceptions, used_start_exceptions):
    valid_starts = set(tbl.start_codons)
    valid_stops = set(tbl.stop_codons)
    print('\n[2] CDS 翻译验证 (密码表 %d):' % tbl.id)
    for feature in sorted(cds, key=lambda item: feature_segments(item)[0][0]):
        display = feature_gene(feature)
        canonical = canonical_gene(feature)
        segments = feature_segments(feature)
        low, high = feature_span(feature)
        strand = '+' if (feature.location.strand or 0) > 0 else '-'
        length = feature_length(feature)
        sequence = feature.extract(gb.seq)

        try:
            codon_start = int(str(feature.qualifiers.get('codon_start', ['1'])[0]))
        except ValueError:
            codon_start = 1
        if codon_start not in (1, 2, 3):
            findings.error('INVALID_CDS: %s 的 /codon_start=%d 非法 (只能是 1/2/3)' % (display, codon_start))
            codon_start = 1
        partial_5p = codon_start != 1
        coding = sequence[codon_start - 1:]
        remainder = len(coding) % 3
        complete = coding[:len(coding) - remainder]
        protein = str(complete.translate(table=tbl)) if len(complete) else ''
        start_codon = str(coding[:3]).upper()
        last3 = str(coding[-3:]).upper() if len(coding) >= 3 else ''

        if remainder == 0 and last3 in valid_stops:
            terminal = ('complete', last3)
        elif remainder == 1 and str(coding[-1:]).upper() == 'T':
            terminal = ('incomplete', str(coding[-1:]).upper())
        elif remainder == 2 and str(coding[-2:]).upper() == 'TA':
            terminal = ('incomplete', str(coding[-2:]).upper())
        else:
            terminal = ('nonstop', last3 or str(coding).upper())

        internal = protein.count('*') - (1 if terminal[0] == 'complete' else 0)

        if partial_5p:
            findings.info("PARTIAL_CDS_5P: %s 5' 端不完整 (/codon_start=%d); 不检查起始密码子"
                          % (display, codon_start))
        elif start_codon not in valid_starts:
            key = (canonical, start_codon)
            if key in start_exceptions:
                used_start_exceptions.add(key)
                findings.review(
                    'NONCANONICAL_START_REVIEW: %s 起始密码子 %s 不在密码表 %d 的合法起始集合内, '
                    '已按 --tolerate-start 声明为类群已知例外; 仍需保留证据等级与不确定性'
                    % (display, start_codon, tbl.id))
            else:
                findings.error(
                    'INVALID_CDS: %s 起始密码子 %s 不在密码表 %d 的合法起始集合内 '
                    '(NONCANONICAL_START_REVIEW); 若为类群已知例外, 请用 --tolerate-start "%s:%s" '
                    '逐条确认并附文献/同源证据'
                    % (display, start_codon, tbl.id, canonical, start_codon))

        if terminal[0] == 'incomplete':
            findings.review(
                '%s 终止密码子不完整 (T/TA 前缀: %s); 与转录后多聚腺苷酸化相容, '
                '但需转录本或近缘全长同源证据; note 属自述, 须交叉核验' % (display, terminal[1]))
        elif terminal[0] == 'nonstop':
            findings.error(
                'INVALID_CDS: %s 末端 %s 既非合法终止密码子, 也不是 T/TA 前缀; '
                'note 不能豁免, 需重核读码框与边界' % (display, terminal[1]))

        if internal > 0:
            if 'transl_except' in feature.qualifiers:
                findings.review(
                    '%s 有 %d 个内部终止密码子, 但声明了 /transl_except; 记 REVIEW, 需独立证据'
                    % (display, internal))
            else:
                findings.error('%s 有 %d 个内部终止密码子' % (display, internal))

        flag = 'OK' if not (internal or terminal[0] == 'nonstop') else 'CHECK'
        location = 'join(%s)' % ','.join('%d-%d' % (s + 1, e) for s, e in segments) \
            if len(segments) > 1 else '%d-%d' % (low, high)
        print('  %-8s %-12s %s %4dbp %4daa  起始=%s 终止=%s%s %s' % (
            display[:8], location, strand, length, len(protein), start_codon or '-', terminal[1],
            ' (不完整)' if terminal[0] == 'incomplete' else
            (" (5' partial)" if partial_5p else ''), flag))


def trna_findings(findings, trnas):
    print('\n[3] tRNA 检查 (%d 个, 长度 %d-%dbp):' % (len(trnas), *TRNA_LENGTH_RANGE))
    short = []
    for feature in sorted(trnas, key=lambda item: feature_segments(item)[0][0]):
        display = feature_gene(feature)
        identity = resolve_trna_identity(feature)
        length = feature_length(feature)
        flag = ''
        if length < TRNA_HARD_MIN:
            if 'note' in feature.qualifiers:
                flag = ' WARN-过短(有 note, 需核坐标)'
                findings.review('tRNA %s 长度 %dbp < %d (有 note: %s); 坐标来源需核验'
                                % (display, length, TRNA_HARD_MIN, str(feature.qualifiers['note'][0])[:60]))
            else:
                flag = ' ERROR-过短'
                findings.error('tRNA %s 长度 %dbp 异常短 (<%d), 无 note 说明; 核坐标来源'
                               % (display, length, TRNA_HARD_MIN))
                short.append(display)
        elif not (TRNA_LENGTH_RANGE[0] <= length <= TRNA_LENGTH_RANGE[1]):
            flag = ' WARN-长度偏离'
            findings.review('tRNA %s 长度 %dbp 偏离 %d-%dbp (长度只触发检查, 不用于判定身份)'
                            % (display, length, *TRNA_LENGTH_RANGE))
        anticodon = feature.qualifiers.get('anticodon', ['?'])[0]
        if anticodon == '?':
            findings.review('tRNA %s 缺少 anticodon qualifiers' % display)
        if identity is None:
            findings.review('UNDETERMINED_TRNA: %s 类型未确定 (需反密码子/结构/同源证据)' % display)
        print('  %-10s %4dbp  %-12s %s%s' % (
            display[:10], length, identity or '待确定', str(anticodon)[:28], flag))
    if short:
        print('  ! 过短 tRNA: %s (应从 arwen/MITOS2 直接提取坐标)' % ', '.join(short))


def rrna_findings(findings, rnas):
    print('\n[4] rRNA 检查 (长度区间为工程预警):')
    for feature in rnas:
        display = feature_gene(feature)
        identity = canonical_gene(feature)
        low, high = feature_span(feature)
        length = feature_length(feature)
        low_bound, high_bound = RRNA_LENGTH_RANGE.get(identity, RRNA_LENGTH_RANGE['rrns'])
        if not (low_bound <= length <= high_bound):
            findings.review('rRNA %s 长度 %d 超出范围 %d-%d (身份=%s; 只触发检查, 不把参考边界当真值)'
                            % (display, length, low_bound, high_bound, identity))
            print('  %-10s [%5d..%5d] %4dbp (预警 %d-%d) WARN' % (display[:10], low, high, length, low_bound, high_bound))
        else:
            print('  %-10s [%5d..%5d] %4dbp OK' % (display[:10], low, high, length))


def overlap_findings(findings, features, tolerate_pairs, used_tolerate, severity):
    """Record every overlapping pair; >8bp is REVIEW by default (never auto-trim)."""
    print('\n[5] 重叠检查 (逐对记录; >%dbp 触发详细审查):' % OVERLAP_REVIEW_THRESHOLD)
    entries = [(feature, feature_gene(feature), canonical_gene(feature), feature.type,
                feature_segments(feature)) for feature in features]
    long_overlaps = 0
    for index, (feature_a, display_a, canon_a, type_a, segments_a) in enumerate(entries):
        for feature_b, display_b, canon_b, type_b, segments_b in entries[index + 1:]:
            overlap = segment_overlap(segments_a, segments_b)
            if overlap <= 0:
                continue
            pair = frozenset((canon_a, canon_b))
            types = '-'.join(sorted((type_a, type_b)))
            if overlap <= OVERLAP_REVIEW_THRESHOLD:
                findings.info('OVERLAP_SHORT: %s(%s) <-> %s(%s) 重叠 %dbp (<=%dbp, 已记录; '
                              'atp8-atp6 / nad4-nad4l 等属真实特征)'
                              % (display_a, type_a, display_b, type_b, overlap,
                                 OVERLAP_REVIEW_THRESHOLD))
                print('  INFO : %-8s <-> %-8s 重叠 %3dbp (<=%dbp, 记录)' % (
                    display_a, display_b, overlap, OVERLAP_REVIEW_THRESHOLD))
                continue
            long_overlaps += 1
            if pair in tolerate_pairs:
                used_tolerate.add(pair)
                findings.review('OVERLAP_ACCEPTED: %s(%s) <-> %s(%s) 重叠 %dbp (%s); '
                                '已人工审核并保留该注释 —— 不代表已证明该重叠具有功能真实性'
                                % (display_a, type_a, display_b, type_b, overlap, types))
                print('  ACCEPT: %-8s <-> %-8s 重叠 %3dbp (已审核保留, 未证明功能真实性)' % (
                    display_a, display_b, overlap))
            elif severity == 'error':
                findings.error('OVERLAP_LONG: %s(%s) 与 %s(%s) 重叠 %dbp (%s); 已按要求升级为错误'
                               % (display_a, type_a, display_b, type_b, overlap, types))
                print('  ERROR: %-8s <-> %-8s 重叠 %3dbp (%s)' % (display_a, display_b, overlap, types))
            else:
                findings.review('OVERLAP_LONG: %s(%s) 与 %s(%s) 重叠 %dbp (%s); 待核查 —— '
                                '重叠本身不构成裁剪序列的理由' % (display_a, type_a, display_b, type_b, overlap, types))
                print('  WARN : %-8s <-> %-8s 重叠 %3dbp (%s, 待核查)' % (
                    display_a, display_b, overlap, types))
    if not long_overlaps:
        print('  (无 >%dbp 重叠)' % OVERLAP_REVIEW_THRESHOLD)


def strand_findings(findings, cds, expected_plus=EXPECTED_PLUS_STRAND_CDS):
    plus = sum(1 for feature in cds if (feature.location.strand or 0) > 0)
    minus = len(cds) - plus
    print('\n[6] 链分布: 正链 CDS=%d 负链 CDS=%d (昆虫背景预警 %d+/%d-, 非跨类群标准)'
          % (plus, minus, expected_plus, len(cds) - expected_plus))
    if plus == expected_plus:
        return
    if minus == expected_plus:
        findings.info('ORIENTATION: 序列整体方向与参照相反 (原样 %d+/%d-, 反向互补后 %d+/%d-); '
                      '属表示差异, 不记异常' % (plus, minus, minus, plus))
        print('  INFO : 整体为反向互补方向, 归一后满足 %d+/%d-, 不记异常' % (minus, plus))
        return
    findings.review('正链 CDS 数量 %d != %d (昆虫背景预警, 需确认类群与整条序列方向)' % (plus, expected_plus))


def reference_findings(findings, cds, ref, tolerance=CDS_LENGTH_TOLERANCE_PCT):
    print('\n[7] 参考对比 (按归一后基因身份):')
    ref_cds = {}
    for feature in sorted(ref.features, key=lambda item: feature_segments(item)[0][0]):
        if feature.type == 'CDS':
            ref_cds.setdefault(canonical_gene(feature), feature)
    order_query = [canonical_gene(feature) for feature in cds]
    order_ref = [canonical_gene(feature) for feature in
                 sorted((f for f in ref.features if f.type == 'CDS'),
                        key=lambda item: feature_segments(item)[0][0])]
    shared = [name for name in order_query if name in order_ref]
    if shared == order_ref:
        print('  基因顺序匹配(CDS): 一致')
    else:
        print('  基因顺序匹配(CDS): 不同! 查询=%s' % '>'.join(shared))
        findings.review('CDS 基因顺序与参考不同: %s vs %s (顺序差异为 REVIEW, 不自动失败)'
                        % ('>'.join(shared), '>'.join(order_ref)))
    compared = 0
    for feature in sorted(cds, key=lambda item: feature_segments(item)[0][0]):
        identity = canonical_gene(feature)
        reference = ref_cds.get(identity)
        if reference is None:
            continue
        compared += 1
        length = feature_length(feature)
        ref_length = feature_length(reference)
        delta = abs(length - ref_length) / ref_length * 100 if ref_length else 0
        if delta > tolerance:
            findings.review('%s 长度 %d vs 参考 %d (差 %.0f%%)' % (identity, length, ref_length, delta))
        print('  %-8s %5dbp vs 参考 %5dbp (差 %4.0f%%)%s' % (
            identity, length, ref_length, delta, ' WARN' if delta > tolerance else ' OK'))
    if compared == 0:
        findings.info('参考与查询没有可比的共享 CDS 身份 (命名差异已归一, 仍无可比项)')


# -------------------------------------------------------------------------- main
def parse_arguments(argv):
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        sys.exit(0)
    options = {'fn': argv[0], 'table': 5, 'ref': None, 'require_circular': False,
               'tolerate_overlap': [], 'tolerate_start': [], 'overlap_severity': 'warn',
               'allow_atypical': None}
    index = 1

    def need(flag):
        nonlocal index
        if index + 1 >= len(argv) or argv[index + 1].startswith('--'):
            print('ERROR: %s 需要一个参数' % flag)
            sys.exit(1)
        index += 1
        return argv[index]

    while index < len(argv):
        token = argv[index]
        if token == '--table':
            try:
                options['table'] = int(need(token))
            except ValueError:
                print('ERROR: --table 需要整数'); sys.exit(1)
        elif token == '--ref':
            options['ref'] = need(token)
        elif token == '--require-circular':
            options['require_circular'] = True
        elif token == '--tolerate-overlap':
            options['tolerate_overlap'].append(need(token))
        elif token == '--tolerate-start':
            options['tolerate_start'].append(need(token))
        elif token == '--overlap-severity':
            value = need(token)
            if value not in ('error', 'warn'):
                print('ERROR: --overlap-severity 只能是 error 或 warn'); sys.exit(1)
            options['overlap_severity'] = value
        elif token == '--allow-atypical':
            options['allow_atypical'] = need(token)
        else:
            print('ERROR: 未知参数 %s' % token)
            sys.exit(1)
        index += 1
    return options


def load_gb(fn):
    try:
        from Bio import SeqIO
        return SeqIO.read(fn, 'genbank')
    except Exception as exc:
        print('ERROR: 无法读取 GB 文件 %s: %s' % (fn, exc))
        sys.exit(1)


def main():
    options = parse_arguments(sys.argv[1:])
    from Bio.Data import CodonTable
    try:
        tbl = CodonTable.unambiguous_dna_by_id[options['table']]
    except KeyError:
        print('ERROR: 未知遗传密码表 %s' % options['table']); sys.exit(1)

    genbank = load_gb(options['fn'])
    findings = Findings()

    print('=' * 72)
    print('线粒体基因组注释质检: %s' % options['fn'])
    topology = genbank.annotations.get('topology', '')
    print('基因组: %d bp  |  遗传密码表: %d (%s)  |  拓扑: %s'
          % (len(genbank.seq), tbl.id, tbl.names[0], topology or '?'))
    print('=' * 72)

    cds = [f for f in genbank.features if f.type == 'CDS']
    trnas = [f for f in genbank.features if f.type == 'tRNA']
    rnas = [f for f in genbank.features if f.type == 'rRNA']

    if options['require_circular'] and topology.lower() != 'circular':
        findings.error('INVALID_TOPOLOGY: 要求环状拓扑, 但文件 topology=%s' % (topology or '未声明'))

    print('\n[1] 基因完整性: CDS=%d tRNA=%d rRNA=%d (典型后生动物预期 13/22/2)'
          % (len(cds), len(trnas), len(rnas)))
    if options['allow_atypical']:
        print('  ! --allow-atypical 理由: %s (仅基因集数量/身份降为待核查, 不放松其他检查)'
              % options['allow_atypical'])
    else:
        if len(cds) != 13: findings.error('CDS 数量 %d != 13 (典型后生动物预期)' % len(cds))
        if len(trnas) != 22: findings.error('tRNA 数量 %d != 22 (典型后生动物预期)' % len(trnas))
        if len(rnas) != 2: findings.error('rRNA 数量 %d != 2 (典型后生动物预期)' % len(rnas))
    if options['allow_atypical']:
        for label, actual, expected in (('CDS', len(cds), 13), ('tRNA', len(trnas), 22),
                                        ('rRNA', len(rnas), 2)):
            if actual != expected:
                findings.review('待核查(非典型类群: %s) %s 数量 %d != %d'
                                % (options['allow_atypical'], label, actual, expected))
    print('\n[1b] 基因身份:')
    identity_findings(findings, cds, trnas, rnas, options['allow_atypical'])

    start_exceptions = set()
    for entry in options['tolerate_start']:
        gene, _, codon = entry.partition(':')
        if not codon:
            print('ERROR: --tolerate-start 需要 "基因:密码子" 形式, 收到 %s' % entry); sys.exit(1)
        start_exceptions.add((_canonical_key(gene), codon.strip().upper()))
    used_start_exceptions = set()

    cds_findings(findings, genbank, cds, tbl, start_exceptions, used_start_exceptions)
    for gene, codon in sorted(start_exceptions - used_start_exceptions):
        findings.review('未匹配任何起始密码子: --tolerate-start "%s:%s" (配置未生效)' % (gene, codon))

    trna_findings(findings, trnas)
    rrna_findings(findings, rnas)

    tolerate_pairs = set()
    for entry in options['tolerate_overlap']:
        parts = [part.strip() for part in entry.split(',')]
        if len(parts) != 2:
            print('ERROR: --tolerate-overlap 需要 "基因1,基因2" 形式, 收到 %s' % entry); sys.exit(1)
        tolerate_pairs.add(frozenset(_canonical_key(part) for part in parts))
    used_tolerate = set()
    overlap_findings(findings, cds + trnas + rnas, tolerate_pairs, used_tolerate,
                     options['overlap_severity'])
    for pair in sorted(tolerate_pairs - used_tolerate, key=lambda item: sorted(item)):
        findings.review('未匹配任何重叠对: --tolerate-overlap "%s" (配置未生效)' % ','.join(sorted(pair)))

    strand_findings(findings, cds)

    if options['ref']:
        reference_findings(findings, cds, load_gb(options['ref']))

    print('\n' + '=' * 72)
    if findings.errors:
        print('结果: 错误 %d 个 (必须处理), 待核查 %d 个, 记录 %d 条'
              % (len(findings.errors), len(findings.reviews), len(findings.infos)))
        for message in findings.errors:
            print('  [ERROR] %s' % message)
        for message in findings.reviews:
            print('  [WARN ] %s' % message)
        for message in findings.infos:
            print('  [INFO ] %s' % message)
        print('=' * 72)
        print('注意: 本结果为注释内部一致性的判断, 不构成样本 reads 支持该序列的证据; '
              'reads 支持字段为 NOT_ASSESSED。')
        sys.exit(1)
    if findings.reviews:
        print('结果: 无错误, 待核查 %d 个, 记录 %d 条 (退出码 2 不等于通过质量门: '
              '每条待核查项须在案例中逐条入账)' % (len(findings.reviews), len(findings.infos)))
        for message in findings.reviews:
            print('  [WARN ] %s' % message)
        for message in findings.infos:
            print('  [INFO ] %s' % message)
        print('=' * 72)
        print('注意: 本结果为注释内部一致性的判断, 不构成样本 reads 支持该序列的证据; '
              'reads 支持字段为 NOT_ASSESSED。')
        sys.exit(2)
    print('结果: 无错误、无待核查项 (记录 %d 条)' % len(findings.infos))
    for message in findings.infos:
        print('  [INFO ] %s' % message)
    print('=' * 72)
    print('注意: 本结果为注释内部一致性的判断, 不构成样本 reads 支持该序列的证据; '
          'reads 支持字段为 NOT_ASSESSED。')
    sys.exit(0)


if __name__ == '__main__':
    main()

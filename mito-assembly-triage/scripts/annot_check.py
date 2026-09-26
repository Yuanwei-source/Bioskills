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
  2. CDS 翻译: 内部终止 / 起始密码子 / 终止密码子
     partial 由 GenBank location 的 < / > (Biopython 位置对象) 判定, 不由 /codon_start 决定;
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
import json
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
#: 身份未确定的 rRNA 名称: 不能归为 rrnS/rrnL, 也不套用任何长度区间
UNDETERMINED_RRNA_KEYS = {'rrna'}
OVERLAP_REVIEW_THRESHOLD = 8
CDS_LENGTH_TOLERANCE_PCT = 20
EXPECTED_PLUS_STRAND_CDS = 9

#: /product 里的三字母氨基酸 -> 单字母 (用于 tRNA-Leu 这类写法)
AMINO_ACID_3_TO_1 = {
    'ala': 'a', 'arg': 'r', 'asn': 'n', 'asp': 'd', 'cys': 'c', 'gln': 'q',
    'glu': 'e', 'gly': 'g', 'his': 'h', 'ile': 'i', 'leu': 'l', 'lys': 'k',
    'met': 'm', 'phe': 'f', 'pro': 'p', 'ser': 's', 'thr': 't', 'trp': 'w',
    'tyr': 'y', 'val': 'v',
}
#: GenBank location 里的分段与 fuzzy 标记
#: Biopython 打印为 "[<0:100](+)" 而原始文件是 "<1..100", 两种写法都要认
LOCATION_PART_RE = re.compile(r'([<>]?)(\d+)\s*(?:\.\.|:)\s*([<>]?)(\d+)')


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
    if key in ('cob', 'cytb', 'cytb2', 'cytochromeb'):
        return 'cytb'
    if key in ('coi', 'cox1', 'cox1p'):
        return 'cox1'
    if key in ('coii', 'cox2'):
        return 'cox2'
    if key in ('coiii', 'cox3'):
        return 'cox3'
    # ---- /product 自然语言名 (NCBI 常见写法) ----
    match = re.fullmatch(r'nadh?dehydrogenasesubunit(\d+)(l?)', key)
    if match:
        return 'nd' + match.group(1) + match.group(2)
    match = re.fullmatch(r'cytochromecoxidasesubunit(i{1,3}|\d)', key)
    if match:
        return 'cox' + {'i': '1', 'ii': '2', 'iii': '3'}.get(match.group(1), match.group(1))
    match = re.fullmatch(r'atp(?:ase)?synthase[a-z0-9]*?subunit(\d+)', key)
    if match:
        return 'atp' + match.group(1)
    match = re.fullmatch(r'(1[26])sribosomalrna', key)
    if match:
        return 'rrnl' if match.group(1) == '16' else 'rrns'
    match = re.fullmatch(r'trna([a-z]{3})', key)
    if match and match.group(1) in AMINO_ACID_3_TO_1:
        # tRNA-Leu -> trnl (裸名): 类型 (CUN/UUR) 仍需反密码子/结构证据
        return 'trn' + AMINO_ACID_3_TO_1[match.group(1)]
    if key in ('rrna', 'ribosomalrna'):
        return 'rrna'
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


def _partial_from_location_string(feature):
    """Fallback: read partiality from the printed location string.

    Only used when no location part carries a fuzzy position object (e.g. a
    location built by hand).  Note the printed form uses ``:``
    (``[<0:100](+)``) while the GenBank file uses ``..`` (``<1..100``);
    accepting both is what makes this fallback trustworthy.
    """
    parts = LOCATION_PART_RE.findall(str(feature.location))
    if not parts:
        return (False, False)
    first, last = parts[0], parts[-1]
    if (feature.location.strand or 0) >= 0:
        return (first[0] == '<', last[2] == '>')
    return (first[2] == '>', last[0] == '<')


def _object_partial(parts, strand):
    """Partiality from Biopython position objects, or None when they carry none."""
    try:
        from Bio.SeqFeature import AfterPosition, BeforePosition
    except ImportError:
        return None
    if not any(isinstance(part.start, (BeforePosition, AfterPosition))
               or isinstance(part.end, (BeforePosition, AfterPosition)) for part in parts):
        return None                      # no fuzzy information in the objects at all
    first, last = parts[0], parts[-1]
    if strand >= 0:
        return (isinstance(first.start, BeforePosition),
                isinstance(last.end, AfterPosition))
    return (isinstance(first.end, AfterPosition),
            isinstance(last.start, BeforePosition))


def feature_partial_detail(feature):
    """Partiality with its provenance.

    The position objects are **authoritative**.  The location string is only
    consulted when the parts carry no fuzzy position information at all
    (hand-built locations).  When both are available and disagree, the object
    result wins *and* the conflict is reported -- it is never silently resolved
    in favour of the string, because that would silently change a "complete" CDS
    into a "partial" one (or the reverse).

    Returns {'five', 'three', 'source' in ('position_objects', 'location_string',
    'none'), 'conflict', 'string_result'}.
    """
    try:
        parts = list(feature.location.parts)
    except (AttributeError, TypeError):
        return {'five': False, 'three': False, 'source': 'none', 'conflict': False,
                'string_result': (False, False)}
    string_result = _partial_from_location_string(feature)
    if not parts:
        return {'five': string_result[0], 'three': string_result[1], 'source': 'location_string',
                'conflict': False, 'string_result': string_result}
    object_result = _object_partial(parts, feature.location.strand or 0)
    if object_result is None:
        return {'five': string_result[0], 'three': string_result[1],
                'source': 'location_string', 'conflict': False, 'string_result': string_result}
    return {'five': object_result[0], 'three': object_result[1],
            'source': 'position_objects',
            'conflict': string_result != object_result,
            'string_result': string_result}


def feature_partial(feature):
    """(five_prime_partial, three_prime_partial); position objects win.

    ``<``/``>`` are positional (lower/higher coordinate); the strand decides which
    biological end they refer to, and a minus-strand gene's 5' end is the *high*
    coordinate.  Location parts are in transcription order, so for a compound
    (cross-origin ``join``) feature the first part carries the 5' end and the last
    part the 3' end:
      plus  : 5' <= BeforePosition on the first part's start,
              3' <= AfterPosition  on the last part's end
      minus : 5' <= AfterPosition  on the first part's end,
              3' <= BeforePosition on the last part's start
    """
    detail = feature_partial_detail(feature)
    return (detail['five'], detail['three'])


# /transl_except token whitelist: standard three-letter codes plus the INSDC
# special tokens.  TERM is a legal TOKEN but it means "translation stops here",
# so it can never explain an internal stop as a readable codon.
TRANSL_EXCEPT_AA_TOKENS = frozenset(name.lower() for name in (
    'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu',
    'Lys', 'Met', 'Phe', 'Pro', 'Pyl', 'Sec', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val',
    'OTHER', 'fMet', 'TERM'))
TRANSL_EXCEPT_TERMINATION_TOKENS = frozenset(('term',))


def _split_location_list(text):
    """Split on commas that are not nested inside parentheses."""
    pieces, depth, current = [], 0, []
    for char in text:
        if char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
        elif char == ',' and depth == 0:
            pieces.append(''.join(current))
            current = []
            continue
        current.append(char)
    pieces.append(''.join(current))
    return [piece.strip() for piece in pieces]


def _location_positions(text):
    """(1-based genomic positions, is_complement, mixed) for a location, else None.

    Grammar kept deliberately small: ``a``, ``a..b``, ``complement(...)`` and
    ``join(...)``/``order(...)``.  Anything else -- fuzzy markers, ``^``, bare
    text -- returns None so the caller reports the qualifier as unparsed rather
    than silently accepting one substring of it.

    ``mixed`` is True when the leaves disagree about ``complement()``.  Collapsing
    them into a single OR-ed boolean is what previously let
    ``join(complement(10..11),12)`` be treated as an entirely minus-strand
    position, hiding the plus-strand leaf.
    """
    text = text.strip()
    if not text:
        return None
    if text.startswith('complement(') and text.endswith(')'):
        inner = _location_positions(text[len('complement('):-1])
        return None if inner is None else (inner[0], not inner[1], inner[2])
    for keyword in ('join(', 'order('):
        if text.startswith(keyword) and text.endswith(')'):
            positions, flags = [], set()
            for piece in _split_location_list(text[len(keyword):-1]):
                parsed = _location_positions(piece)
                if parsed is None:
                    return None
                positions.extend(parsed[0])
                flags.add(parsed[1])
            if not positions:
                return None
            return (positions, flags == {True}, len(flags) > 1)
    if re.fullmatch(r'\d+', text):
        return ([int(text)], False, False)
    match = re.fullmatch(r'(\d+)\s*\.\.\s*(\d+)', text)
    if match:
        start, end = int(match.group(1)), int(match.group(2))
        return (list(range(start, end + 1)), False, False) if start <= end else None
    return None


def _parse_transl_except_entry(text):
    """(positions, complemented, mixed, amino_acid) from one '(pos:...,aa:...)' chunk."""
    text = text.strip()
    if not (text.startswith('(') and text.endswith(')')):
        return None
    pieces = _split_location_list(text[1:-1])
    if len(pieces) != 2 or not pieces[0].startswith('pos:'):
        return None
    match = re.fullmatch(r'aa:\s*([A-Za-z]+)', pieces[1])
    if match is None:
        return None
    parsed = _location_positions(pieces[0][len('pos:'):])
    if parsed is None:
        return None
    return (parsed[0], parsed[1], parsed[2], match.group(1))


def parse_transl_except(feature):
    """Parsed /transl_except entries plus a count of unreadable ones.

    Each qualifier must be consumed **completely**: one or more ``(pos:...,aa:...)``
    chunks separated by commas.  A qualifier that does not match this grammar in
    full -- trailing text, a missing ``pos:`` token, an unsupported location form --
    is counted in ``unparsed``, never silently dropped and never partially
    accepted.  ``entry['explains_stop']`` separates "legal aa token" (TERM is one)
    from "can stand in for a stop codon" (TERM cannot).
    """
    entries = []
    unparsed = 0
    for value in feature.qualifiers.get('transl_except', []) or []:
        pending = str(value).strip()
        if not pending:
            unparsed += 1
            continue
        local, broken = [], 0
        while pending:
            if not pending.startswith('('):
                broken += 1
                break
            depth, close = 0, -1
            for index, char in enumerate(pending):
                if char == '(':
                    depth += 1
                elif char == ')':
                    depth -= 1
                    if depth == 0:
                        close = index
                        break
            if close < 0:
                broken += 1
                break
            chunk, rest = pending[:close + 1], pending[close + 1:].lstrip()
            parsed = _parse_transl_except_entry(chunk)
            if parsed is None:
                broken += 1
            else:
                positions, complemented, mixed, amino_acid = parsed
                token = amino_acid.strip().lower()
                local.append({
                    'text': chunk.strip(), 'positions': positions,
                    'complemented': complemented, 'mixed': mixed,
                    'amino_acid': amino_acid,
                    'valid_token': token in TRANSL_EXCEPT_AA_TOKENS,
                    'explains_stop': (token in TRANSL_EXCEPT_AA_TOKENS
                                      and token not in TRANSL_EXCEPT_TERMINATION_TOKENS),
                })
            if not rest:
                break
            if not rest.startswith(','):
                broken += 1
                break
            pending = rest[1:].lstrip()
            if not pending:
                # a dangling separator means an entry is missing: the qualifier was
                # not consumed completely and must not be trusted
                broken += 1
                break
        if broken:
            # A partially readable qualifier is NOT partially trusted: salvaging the
            # well-formed prefix would let trailing corruption carry an exception
            # into the analysis.  The whole qualifier is rejected instead.
            unparsed += 1
            continue
        entries.extend(local)
    return entries, unparsed


def coding_position_list(feature):
    """1-based genomic positions of the CDS in TRANSCRIPTION order.

    Verified against ``feature.extract(record)``: for a compound location the
    parts are already in transcription order, and a minus-strand part runs from
    its high coordinate down to its low coordinate.
    """
    positions = []
    for part in feature.location.parts:
        low, high = int(part.start), int(part.end)
        if (part.strand or 0) >= 0:
            positions.extend(range(low + 1, high + 1))
        else:
            positions.extend(range(high, low, -1))
    return positions


def cds_codon_positions(feature, codon_start=1):
    """[[genomic positions of codon 0], ...] in transcription order.

    Matching is by exact position set: a /transl_except whose declared range is
    not *exactly* the three genomic positions of one codon explains nothing.  A
    wrong-frame 3 nt window (e.g. 11..13 for a stop at 10..12) must not be
    snapped onto the nearest codon by an integer division.
    """
    positions = coding_position_list(feature)[max(0, codon_start - 1):]
    return [positions[index:index + 3] for index in range(0, len(positions) - 2, 3)]


def terminal_is_complete(last3, remainder, valid_stops):
    return remainder == 0 and last3 in valid_stops


def oriented_gene_cycle(features):
    """[(canonical_name, strand)] in coordinate order."""
    ordered = sorted(features, key=lambda item: feature_segments(item)[0][0])
    return [(canonical_gene(feature), '+' if (feature.location.strand or 0) > 0 else '-')
            for feature in ordered]


def _reverse_cycle(cycle):
    return [(name, '-' if strand == '+' else '+') for name, strand in reversed(cycle)]


def _induced_pairs(cycle, shared):
    """Cyclic adjacency pairs over the genes present on both sides."""
    nodes = [node for node in cycle if node[0] in shared]
    if len(nodes) < 2:
        return set()
    return {(nodes[index], nodes[(index + 1) % len(nodes)]) for index in range(len(nodes))}


def cycles_equivalent(cycle_a, cycle_b):
    """Equal up to cyclic rotation and full reverse complement (representation changes)."""
    names_a = {name for name, _ in cycle_a}
    names_b = {name for name, _ in cycle_b}
    if len(cycle_a) != len(cycle_b) or names_a != names_b:
        return False
    pairs_a = _induced_pairs(cycle_a, names_a)
    return (pairs_a == _induced_pairs(cycle_b, names_a)
            or pairs_a == _induced_pairs(_reverse_cycle(cycle_b), names_a))


def adjacency_diff(cycle_a, cycle_b):
    """Adjacency pairs differing between two arrangements (rotation/RC tolerant).

    An empty result means the shared genes keep the same cyclic arrangement;
    gene-set loss or gain is reported separately and does not show up here.
    """
    shared = {name for name, _ in cycle_a} & {name for name, _ in cycle_b}
    reference = _induced_pairs(cycle_a, shared)
    forward = _induced_pairs(cycle_b, shared)
    if reference == forward or reference == _induced_pairs(_reverse_cycle(cycle_b), shared):
        return []
    return sorted(reference.symmetric_difference(forward))


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

    for label, features, expected, undetermined_key in (('CDS', cds, EXPECTED_CDS, None),
                                                        ('rRNA', rnas, EXPECTED_RRNA, 'rrna')):
        names = [canonical_gene(feature) for feature in features]
        counts = Counter(names)
        missing = set(expected - set(names))
        extra = set(names) - expected
        duplicates = sorted(name for name, count in counts.items() if count > 1)
        if undetermined_key and counts.get(undetermined_key):
            # 名称未确定 -> 无法判定该类是否缺失, 也不套用具体区间
            for member in expected:
                missing.discard(member)
            extra.discard(undetermined_key)
            findings.review(
                'UNDETERMINED_RRNA: %d 个 rRNA 名称未确定 (?%s?), 无法判定 %s 是否缺失, '
                '也不套用任何长度区间 —— 需同源/结构/邻域证据'
                % (counts[undetermined_key], undetermined_key, '/'.join(sorted(expected))))
        missing = sorted(missing)
        extra = sorted(extra)
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


def _registry_site(item, path, index):
    """('index', n) | ('positions', frozenset) | ('gene_wide', None) for one record.

    A /transl_except record must say WHICH stop it covers.  Without that binding a
    single record silently becomes a gene-wide reassignment rule and clears every
    matching stop in the CDS.  Reusing a record across sites is only allowed when
    the broader scope is declared explicitly and audited.
    """
    scope = item.get('scope')
    has_index = item.get('codon_index') is not None
    has_pos = item.get('pos') is not None
    if scope is not None:
        if scope.strip() != 'gene_wide' or has_index or has_pos:
            print('ERROR: 例外记录 %s 的 exceptions[%d].scope 只支持 "gene_wide", '
                  '且不得与 codon_index/pos 同时出现' % (path, index))
            sys.exit(1)
        return ('gene_wide', None)
    if has_index and has_pos:
        print('ERROR: 例外记录 %s 的 exceptions[%d] 只能给 codon_index 或 pos 之一'
              % (path, index))
        sys.exit(1)
    if has_index:
        return ('index', int(item['codon_index']))
    if has_pos:
        parsed = _location_positions(item['pos'])
        if parsed is None or parsed[2]:
            print('ERROR: 例外记录 %s 的 exceptions[%d].pos 无法解析或链方向混杂: %r'
                  % (path, index, item['pos']))
            sys.exit(1)
        return ('positions', frozenset(parsed[0]))
    print('ERROR: 例外记录 %s 的 exceptions[%d] 缺少位点绑定: 必须给 codon_index 或 pos, '
          '或显式 scope="gene_wide"; 不绑定位点的记录不得用来验证任意 stop'
          % (path, index))
    sys.exit(1)


def _transl_except_evidence(records, table, expected_taxon, codon_index, positions):
    """(record, scope_text) for a declared stop-codon exception; record None = rejected.

    Layered on purpose: matching a position and reading frame is a SYNTAX fact while
    accepting an exception is a BIOLOGICAL claim.  Codon meaning depends on the
    taxon and the genetic code, so no global codon->amino-acid table is consulted.
    Acceptance is fail-closed on three axes:

    * sample binding -- an explicit ``--taxon`` must be given and must match the
      record; a record whose taxon cannot be tied to this sample proves nothing;
    * site binding   -- the record must cover THIS stop (codon index or exact
      genomic positions).  Only an explicit ``scope="gene_wide"`` record may serve
      more than one site;
    * completeness   -- taxon/source/rationale/transl_table present and consistent.
    """
    if not records:
        return None, ' (\u672a\u767b\u8bb0\u8be5 gene/codon/aa \u7ec4\u5408)'
    if not expected_taxon:
        return None, (' (\u672a\u63d0\u4f9b --taxon: registry \u7684 taxon \u65e0\u6cd5\u7ed1\u5b9a\u5230\u672c\u6837\u672c, '
                      '\u65e0\u7c7b\u7fa4\u524d\u63d0\u7684\u8bb0\u5f55\u4e0d\u5f97\u4f5c\u4e3a\u5df2\u9a8c\u8bc1\u4f8b\u5916)')
    wanted = frozenset(positions)
    notes = []
    for record in records:
        kind, value = record['site']
        if kind == 'index' and value != codon_index:
            notes.append('\u4f4d\u70b9 codon_index=%s \u4e0d\u5339\u914d (\u672c\u5bc6\u7801\u5b50\u4e3a %s)'
                         % (value, codon_index))
            continue
        if kind == 'positions' and value != wanted:
            notes.append('\u4f4d\u70b9 pos=%s \u4e0d\u5339\u914d' % record.get('pos'))
            continue
        scope_text = 'gene_wide' if kind == 'gene_wide' else str(
            record.get('pos') or ('codon_index=%s' % value))
        incomplete = [name for name in ('taxon', 'source', 'rationale')
                      if not str(record.get(name) or '').strip()]
        if record.get('transl_table') is None:
            incomplete.append('transl_table')
        if incomplete:
            notes.append('\u8bb0\u5f55\u7f3a\u5c11 %s' % '/'.join(incomplete))
            continue
        if int(record['transl_table']) != int(table.id):
            notes.append('\u8bb0\u5f55 transl_table=%s \u4e0e\u672c\u6b21 --table %s \u4e0d\u4e00\u81f4'
                         % (record['transl_table'], table.id))
            continue
        if str(record['taxon']).strip().lower() != str(expected_taxon).strip().lower():
            notes.append('\u8bb0\u5f55 taxon=%s \u4e0e\u672c\u6b21 --taxon %s \u4e0d\u4e00\u81f4'
                         % (record['taxon'], expected_taxon))
            continue
        return record, scope_text
    return None, ' (%s)' % '; '.join(notes)


def cds_findings(findings, gb, cds, tbl, start_exceptions, used_start_exceptions,
                 registry=None, expected_taxon=None):
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

        partial_detail = feature_partial_detail(feature)
        five_p, three_p = partial_detail['five'], partial_detail['three']
        if partial_detail['conflict']:
            findings.review(
                'PARTIAL_SOURCE_CONFLICT: %s 的位置对象与 location 字符串不一致 '
                '(位置对象=%s, 字符串=%s); 以位置对象为准, 请核对注释来源'
                % (display, (five_p, three_p), partial_detail['string_result']))
        declared_table = str(feature.qualifiers.get('transl_table', [''])[0])
        if declared_table and declared_table.strip() and declared_table.strip() != str(tbl.id):
            findings.review('TABLE_CONFLICT: %s 声明 /transl_table=%s 与 --table %d 不一致; '
                            '本检查按 --table %d 翻译'
                            % (display, declared_table, tbl.id, tbl.id))
        try:
            codon_start = int(str(feature.qualifiers.get('codon_start', ['1'])[0]))
        except ValueError:
            codon_start = 1
        if codon_start not in (1, 2, 3):
            findings.error('INVALID_CDS: %s 的 /codon_start=%d 非法 (只能是 1/2/3)' % (display, codon_start))
            codon_start = 1
        if codon_start != 1 and not five_p:
            # 完整标注的 CDS 与 /codon_start 矛盾: 是真 5' 缺失, 还是注释不一致?
            findings.review("CODON_START_CONFLICT: %s 的 location 标注为完整, 却设 /codon_start=%d; "
                            "注释自相矛盾, 需先确认 5' 端是否真的缺失" % (display, codon_start))
        coding = sequence[codon_start - 1:]
        remainder = len(coding) % 3
        complete = coding[:len(coding) - remainder]
        protein = str(complete.translate(table=tbl)) if len(complete) else ''
        start_codon = str(coding[:3]).upper()
        last3 = str(coding[-3:]).upper() if len(coding) >= 3 else ''

        # /transl_except 必须解释“具体哪个内部终止”, 而不是只要存在就全免.
        # 位置/读框一致只是语法事实; “是否接受为生物学例外”必须有审计证据.
        exceptions, unparsed = parse_transl_except(feature)
        stop_codons = [index for index, amino_acid in enumerate(protein) if amino_acid == '*']
        if terminal_is_complete(last3, remainder, valid_stops) and protein.endswith('*'):
            stop_codons = stop_codons[:-1]
        codon_positions = cds_codon_positions(feature, codon_start)
        minus_strand = (feature.location.strand or 0) < 0
        transl_registry = (registry or {}).get('transl_except') or {}
        explained = set()
        for entry in exceptions:
            amino_acid = entry['amino_acid']
            position_text = entry['text']
            if not entry['valid_token']:
                index, reason = None, ' (aa:%s 不是合法的例外氨基酸)' % amino_acid
            elif not entry['explains_stop']:
                index, reason = None, (' (aa:TERM 表示终止, 不能被当作可继续翻译的氨基酸'
                                       '来解释内部终止)')
            elif entry['mixed']:
                index, reason = None, (' (pos 由不同链方向的片段组成, 不能作为链一致的单密码'
                                       '子位置)')
            elif entry['complemented'] != minus_strand:
                index, reason = None, (' (pos 的链方向与 CDS 不一致: %s链 CDS 需要 pos:%s)'
                                       % ('负' if minus_strand else '正',
                                          'complement(a..b)' if minus_strand else 'a..b'))
            elif len(entry['positions']) != 3:
                index, reason = None, \
                    ' (pos 范围 %d nt, 不是单个密码子)' % len(entry['positions'])
            else:
                wanted = set(entry['positions'])
                index = next((number for number, codon in enumerate(codon_positions)
                              if codon and set(codon) == wanted), None)
                reason = ('' if index is not None else
                          ' (pos:%s 与任何一个密码子的基因组位置都不完全相同)' % position_text)
            if index is None or index not in stop_codons:
                findings.review('TRANSL_EXCEPT_UNEXPLAINED: %s 声明的 /transl_except '
                                '%s 未对应任何内部终止密码子%s'
                                % (display, position_text, reason))
                continue
            actual_codon = str(coding[3 * index:3 * index + 3]).upper()
            findings.info('TRANSL_EXCEPT_MATCHED: %s 位置 %s 精确对应密码子 %d 的内部终止 %s '
                          '(位置/读框事实, 单独的 MATCHED 不等于已接受)'
                          % (display, position_text, index + 1, actual_codon))
            if not canonical:
                # defence in depth: an unidentified CDS must never have a
                # gene-keyed record applied to it (see the loader's empty-selector
                # rejection)
                record, note = None, (' (该 CDS 缺少可识别的 gene 身份 (/gene 与 /product 归一化后'
                                      '仍为空), registry 证据无法绑定到基因)')
            else:
                records = transl_registry.get((canonical, actual_codon,
                                               amino_acid.strip().lower()))
                record, note = _transl_except_evidence(records, tbl, expected_taxon,
                                                       index + 1, entry['positions'])
            if record is None:
                findings.review(
                    'TRANSL_EXCEPT_DECLARED_UNVERIFIED: %s 的 /transl_except %s 声明密码子 %d 的'
                    ' %s 由 %s 替代, 但缺少适用于本样本与本位点的已审计证据%s; '
                    '该内部终止仍按 ERROR 处理 (请用 --exception-registry 记录 gene/codon/'
                    'amino_acid + 位点绑定 (codon_index 或 pos, 或显式 scope="gene_wide") + '
                    'transl_table/taxon/source/rationale, 并用 --taxon 给出本样本类群)'
                    % (display, position_text, index + 1, actual_codon, amino_acid, note))
            else:
                explained.add(index)
                findings.info(
                    'TRANSL_EXCEPT_VALIDATED: %s 密码子 %d 的内部终止 %s->%s 已按审计记录接受 '
                    '(scope=%s taxon=%s transl_table=%s source=%s rationale=%s)'
                    % (display, index + 1, actual_codon, amino_acid, note, record.get('taxon'),
                       record.get('transl_table'), record.get('source'), record.get('rationale')))
        if unparsed:
            findings.review('TRANSL_EXCEPT_UNPARSED: %s 的 /transl_except 有 %d 条无法完整解析 '
                            '(要求 (pos:a..b|pos:complement(a..b)|pos:join(...), aa:三字母代码)); '
                            '未能采纳的例外不构成豁免' % (display, unparsed))
        if exceptions and not stop_codons:
            findings.review('TRANSL_EXCEPT_UNEXPLAINED: %s 声明了 /transl_except, '
                            '但该 CDS 没有内部终止密码子, 声明无对应异常' % display)
        unexplained_stops = [index for index in stop_codons if index not in explained]

        if remainder == 0 and last3 in valid_stops:
            terminal = ('complete', last3)
        elif remainder == 1 and str(coding[-1:]).upper() == 'T':
            terminal = ('incomplete', str(coding[-1:]).upper())
        elif remainder == 2 and str(coding[-2:]).upper() == 'TA':
            terminal = ('incomplete', str(coding[-2:]).upper())
        else:
            terminal = ('nonstop', last3 or str(coding).upper())

        internal = protein.count('*') - (1 if terminal[0] == 'complete' else 0)

        if five_p:
            findings.info("PARTIAL_CDS_5P: %s 5' 端不完整 (location 标记); 不检查起始密码子" % display)
        elif start_codon not in valid_starts:
            key = (canonical, start_codon)
            if key in start_exceptions:
                used_start_exceptions.add(key)
                record = ((registry or {}).get('start') or {}).get(key)
                if record:
                    findings.review(
                        'NONCANONICAL_START_REVIEW: %s 起始密码子 %s 按已审计例外记录接受; '
                        'taxon=%s source=%s rationale=%s'
                        % (display, start_codon, record.get('taxon'), record.get('source'),
                           record.get('rationale')))
                    incomplete = [name for name in ('taxon', 'source', 'rationale')
                                  if not str(record.get(name) or '').strip()]
                    if incomplete:
                        findings.review(
                            'EXCEPTION_RECORD_INCOMPLETE: %s:%s 的已审计记录缺少 %s; '
                            '该例外只能维持 REVIEW' % (canonical, start_codon, '/'.join(incomplete)))
                    if expected_taxon and record.get('taxon') and \
                            str(record['taxon']).strip().lower() != str(expected_taxon).strip().lower():
                        findings.review(
                            'EXCEPTION_TAXON_MISMATCH: %s:%s 的记录 taxon=%s 与本次 --taxon %s '
                            '不一致; 该例外只能维持 REVIEW'
                            % (canonical, start_codon, record.get('taxon'), expected_taxon))
                    if not expected_taxon:
                        # same fail-closed rule as /transl_except: a taxon that cannot be
                        # tied to this sample must not be presented as audited evidence
                        findings.review(
                            'EXCEPTION_TAXON_UNVERIFIED: %s:%s 引用了已审计记录, 但未提供 '
                            '--taxon, 无法确认该记录的类群适用于本样本'
                            % (canonical, start_codon))
                else:
                    findings.review(
                        'NONCANONICAL_START_REVIEW: %s 起始密码子 %s 不在密码表 %d 的合法起始集合内, '
                        '已按 --tolerate-start 接受' % (display, start_codon, tbl.id))
                    findings.review(
                        'EXCEPTION_NOT_REGISTERED: %s:%s 未引用已审计例外记录 '
                        '(--exception-registry), 不得作为已验证结论; 请在案例中记录 '
                        'taxon/source/rationale' % (canonical, start_codon))
            else:
                findings.error(
                    'INVALID_CDS: %s 起始密码子 %s 不在密码表 %d 的合法起始集合内 '
                    '(NONCANONICAL_START_REVIEW); 若为类群已知例外, 请用 --tolerate-start "%s:%s" '
                    '逐条确认并附文献/同源证据'
                    % (display, start_codon, tbl.id, canonical, start_codon))

        if three_p:
            findings.info("PARTIAL_CDS_3P: %s 3' 端不完整 (location 标记); 不要求终止密码子" % display)
        elif terminal[0] == 'incomplete':
            findings.review(
                '%s 终止密码子不完整 (T/TA 前缀: %s); 与转录后多聚腺苷酸化相容, '
                '但需转录本或近缘全长同源证据; note 属自述, 须交叉核验' % (display, terminal[1]))
        elif terminal[0] == 'nonstop':
            findings.error(
                'INVALID_CDS: %s 末端 %s 既非合法终止密码子, 也不是 T/TA 前缀; '
                'note 不能豁免, 需重核读码框与边界' % (display, terminal[1]))

        if unexplained_stops:
            findings.error('%s 有 %d 个内部终止密码子 (未被 /transl_except 解释)'
                           % (display, len(unexplained_stops)))

        flag = 'OK' if not (unexplained_stops or terminal[0] == 'nonstop') else 'CHECK'
        location = 'join(%s)' % ','.join('%d-%d' % (s + 1, e) for s, e in segments) \
            if len(segments) > 1 else '%d-%d' % (low, high)
        print('  %-8s %-12s %s %4dbp %4daa  起始=%s 终止=%s%s %s' % (
            display[:8], location, strand, length, len(protein), start_codon or '-', terminal[1],
            ' (不完整)' if terminal[0] == 'incomplete' else
            (" (5' partial)" if five_p else (" (3' partial)" if three_p else '')), flag))


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
        if identity not in RRNA_LENGTH_RANGE:
            # 身份未确定: 不能套用 rrnS/rrnL 的区间
            findings.info('UNDETERMINED_RRNA: %s 身份未确定 (%s), 跳过长度区间判定'
                          % (display, identity))
            print('  %-10s [%5d..%5d] %4dbp (身份未确定, 跳过区间) INFO'
                  % (display[:10], low, high, length))
            continue
        low_bound, high_bound = RRNA_LENGTH_RANGE[identity]
        if not (low_bound <= length <= high_bound):
            findings.review('rRNA %s 长度 %d 超出范围 %d-%d (身份=%s; 只触发检查, 不把参考边界当真值)'
                            % (display, length, low_bound, high_bound, identity))
            print('  %-10s [%5d..%5d] %4dbp (预警 %d-%d) WARN'
                  % (display[:10], low, high, length, low_bound, high_bound))
        else:
            print('  %-10s [%5d..%5d] %4dbp OK' % (display[:10], low, high, length))


def overlap_findings(findings, features, tolerate_pairs, used_tolerate, severity):
    """Record every overlapping pair; >8bp is REVIEW by default (never auto-trim)."""
    print('\n[5] 重叠检查 (逐对记录; >%dbp 触发详细审查):' % OVERLAP_REVIEW_THRESHOLD)
    entries = [(feature, feature_gene(feature), canonical_gene(feature), feature.type,
                feature_segments(feature)) for feature in features]
    name_counts = Counter(entry[2] for entry in entries)
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
                if any(name_counts.get(name, 0) > 1 for name in pair):
                    findings.review(
                        'TOLERATE_OVERLAP_AMBIGUOUS: --tolerate-overlap "%s" 涉及的基因身份在记录中出现多次, '
                        '无法定位到具体 feature/坐标对; 请改用坐标形式并逐条审核'
                        % ','.join(sorted(pair)))
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
    print('\n[7] 参考对比 (按归一后基因身份; 顺序按环状邻接比较):')
    ref_cds = {}
    for feature in sorted(ref.features, key=lambda item: feature_segments(item)[0][0]):
        if feature.type == 'CDS':
            ref_cds.setdefault(canonical_gene(feature), feature)
    query_cycle = oriented_gene_cycle(cds)
    ref_cycle = oriented_gene_cycle([f for f in ref.features if f.type == 'CDS'])
    query_names = {name for name, _ in query_cycle}
    ref_names = {name for name, _ in ref_cycle}
    missing = sorted(ref_names - query_names)
    extra = sorted(query_names - ref_names)
    if missing or extra:
        findings.review('GENE_SET_DIFF: 与参考相比 缺失=%s 多余=%s '
                        '(基因集差异与顺序差异分开报告)'
                        % (missing or '无', extra or '无'))
    differences = adjacency_diff(ref_cycle, query_cycle)
    if differences:
        findings.review('ARRANGEMENT_DIFF: 共有基因的环状邻接关系与参考不同 (%d 处): %s'
                        % (len(differences), differences[:6]))
        print('  邻接关系: 不同 (差异 %d 处)' % len(differences))
    else:
        print('  邻接关系一致 (已归一环状旋转、起点与整链反向互补; 共享基因 %d 个)'
              % len(query_names & ref_names))
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
def load_exception_registry(path):
    """Audited exceptions, in two sections.

    ``start``         : (gene, codon) -> taxon/source/rationale  (non-canonical starts)
    ``transl_except`` : (gene, codon, amino_acid) -> taxon/transl_table/source/rationale

    An entry carrying ``amino_acid`` is a /transl_except record, anything else is a
    start-codon record.  The CLI flag only *selects* an entry; the biological
    justification lives in the registry so that accepting an exception is never the
    same as proving it.
    """
    try:
        with open(path, encoding='utf-8') as handle:
            data = json.load(handle)
    except (OSError, ValueError) as exc:
        print('ERROR: 无法读取例外记录 %s: %s' % (path, exc))
        sys.exit(1)
    if not isinstance(data, dict):
        print('ERROR: 例外记录 %s 的顶层必须是 JSON 对象 ({"exceptions": [...]}), 实际为 %s'
              % (path, type(data).__name__))
        sys.exit(1)
    raw = data.get('exceptions', [])
    if not isinstance(raw, list):
        print('ERROR: 例外记录 %s 的 "exceptions" 必须是数组, 实际为 %s'
              % (path, type(raw).__name__))
        sys.exit(1)
    registry = {'start': {}, 'transl_except': {}}
    for index, item in enumerate(raw):
        if not isinstance(item, dict):
            print('ERROR: 例外记录 %s 的 exceptions[%d] 必须是对象, 实际为 %s'
                  % (path, index, type(item).__name__))
            sys.exit(1)
        # Evidence fields describe taxon/source/rationale: an array, object or number
        # there would be stringified into a non-empty "audited" value, so the type is
        # validated at LOAD time and a malformed registry never applies partially.
        for key in ('gene', 'codon', 'amino_acid', 'taxon', 'source', 'rationale',
                    'pos', 'scope'):
            value = item.get(key)
            if value is not None and not isinstance(value, str):
                print('ERROR: 例外记录 %s 的 exceptions[%d].%s 必须是字符串, 实际为 %s'
                      % (path, index, key, type(value).__name__))
                sys.exit(1)
        for key in ('transl_table', 'codon_index'):
            value = item.get(key)
            if value is not None and (isinstance(value, bool) or not isinstance(value, int)
                                      or value < 1):
                print('ERROR: 例外记录 %s 的 exceptions[%d].%s 必须是正整数, 实际为 %r'
                      % (path, index, key, value))
                sys.exit(1)
        # Selectors must not be null: the PRESENCE of the key is what makes a record
        # a transl_except record, so "amino_acid": null is a malformed transl_except
        # record, not a start-codon record.
        for key in ('gene', 'codon', 'amino_acid'):
            if key in item and item[key] is None:
                print('ERROR: 例外记录 %s 的 exceptions[%d].%s 键存在但为 null; selector 必须给出'
                      '有效值 (键存在即表示这是 transl_except 记录, 不得降格为 start 记录)'
                      % (path, index, key))
                sys.exit(1)
        gene = _canonical_key(item.get('gene', ''))
        codon = str(item.get('codon', '')).upper()
        # Classify by KEY PRESENCE, not truthiness: an explicit amino_acid of "" or
        # null is a malformed transl_except record, not a start-codon record.
        has_amino_acid = 'amino_acid' in item
        if has_amino_acid:
            amino_acid = item.get('amino_acid')
            # An empty selector matches an empty selector: a CDS without /gene or
            # /product canonicalises to '' too, so an empty gene name would bind
            # "evidence" to no gene identity at all and clear the stop ERROR.
            if not gene:
                print('ERROR: 例外记录 %s 的 exceptions[%d].gene 归一化后为空 (%r); '
                      'transl_except 证据必须绑定一个可识别的基因身份'
                      % (path, index, item.get('gene')))
                sys.exit(1)
            if not re.fullmatch(r'[ACGTURYKMSWBDHVN]{3}', codon):
                print('ERROR: 例外记录 %s 的 exceptions[%d].codon 必须是三个 IUPAC 碱基, '
                      '实际为 %r' % (path, index, item.get('codon')))
                sys.exit(1)
            amino_key = str(amino_acid).strip().lower()
            if amino_key not in TRANSL_EXCEPT_AA_TOKENS:
                print('ERROR: 例外记录 %s 的 exceptions[%d].amino_acid 不是合法的例外氨基酸 '
                      'token, 实际为 %r' % (path, index, amino_acid))
                sys.exit(1)
            record = {'taxon': item.get('taxon'), 'source': item.get('source'),
                      'rationale': item.get('rationale'),
                      'transl_table': item.get('transl_table'),
                      'site': _registry_site(item, path, index), 'pos': item.get('pos')}
            key = (gene, codon, amino_key)
            bucket = registry['transl_except'].setdefault(key, [])
            # The duplicate identity must include EVERY axis the runtime uses to
            # SELECT a record, otherwise a shared registry cannot hold the same
            # homologous site for several taxa or genetic codes.
            identity = (record['site'], str(record.get('taxon') or '').strip().lower(),
                        record.get('transl_table'))
            if any(existing['identity'] == identity for existing in bucket):
                print('ERROR: 例外记录 %s 的 exceptions[%d] 与已有记录重复 '
                      '(同一 gene/codon/amino_acid/位点/taxon/transl_table): %s/%s/%s %s '
                      '-> 重复 key 会被静默覆盖, 因此按格式错误处理'
                      % (path, index, gene, codon, amino_acid, identity))
                sys.exit(1)
            record['identity'] = identity
            bucket.append(record)
        else:
            registry['start'][(gene, codon)] = {
                'taxon': item.get('taxon'), 'source': item.get('source'),
                'rationale': item.get('rationale')}
    return registry


def parse_arguments(argv):
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        sys.exit(0)
    options = {'fn': argv[0], 'table': 5, 'ref': None, 'require_circular': False,
               'tolerate_overlap': [], 'tolerate_start': [], 'overlap_severity': 'warn',
               'allow_atypical': None, 'exception_registry': None, 'taxon': None}
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
        elif token == '--exception-registry':
            options['exception_registry'] = need(token)
        elif token == '--taxon':
            options['taxon'] = need(token)
        else:
            print('ERROR: 未知参数 %s' % token)
            sys.exit(1)
        index += 1
    if options['require_circular']:
        print('CIRCULAR_DECLARATION_CHECK: --require-circular 只校验 GenBank 的 topology 声明, '
              '不等于物理环化证据; 环化证据必须来自接缝 reads / 组装图')
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
    registry = load_exception_registry(options['exception_registry']) \
        if options['exception_registry'] else {}

    cds_findings(findings, genbank, cds, tbl, start_exceptions, used_start_exceptions,
                 registry=registry, expected_taxon=options['taxon'])
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

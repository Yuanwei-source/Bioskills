#!/usr/bin/env python3

"""
COX1 物种鉴定: 提取 COX1 (或给定区域), 提交 NCBI blastn, 解析**结构化**结果

必须显式给出 --coords (本脚本不自动识别 COX1) 与 --allow-public-upload。

输出与判读:
  - 报告每个候选的 identity 与 **query coverage** (多 HSP 的 query 区间并集),
    而不是单个 HSP 的 Identities 分母;
  - 结果只给候选归属 (provisional_candidate / ambiguous / insufficient),
    不给物种级确定结论;
  - 阈值 (85% identity / 80% coverage / top 与次优相差 1%) 是工程启发式,
    无类群特异性依据。

退出码 (任务状态与判读状态分开):
  0 = 得到判读; 1 = insufficient / no_match (分析结论, 不是程序故障);
  2 = 拒绝执行 (缺 --allow-public-upload 或坐标非法); 3 = 网络或结果格式故障。

用法:
  python3 cox1_id.py <genome.fasta> --allow-public-upload --coords 1353,2891
                    [--max-results 5] [--output-json results.json]
依赖: 网络 (NCBI BLAST URL API), 可选 BioPython

--output-json 会把每个候选的 accession、逐 HSP 原始坐标/identity/bitscore 与
blocker 原因持久化; 无论最终判读是什么都会写出。

坐标与数值边界（#5）：越界/非有限/方向矛盾 = 格式故障（退出码 3）；
长度字段缺失只记 unchecked；identity==0 记为 blocker（zero_identity）而不是格式故障。
"""
import json
import math
import sys
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ElementTree

BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
MIN_IDENTITY = 85.0
MIN_QUERY_COVERAGE = 0.8
MIN_QUERY_LENGTH = 400
MAX_QUERY_LENGTH = 5000
POLL_ATTEMPTS = 40
POLL_INTERVAL = 12


# ------------------------------------------------------------------------- input
def read_single_fasta(fn):
    records = []
    current = []
    for line in open(fn, encoding='utf-8'):
        line = line.strip()
        if line.startswith('>'):
            if current:
                records.append(''.join(current).upper())
                current = []
        else:
            current.append(line)
    if current:
        records.append(''.join(current).upper())
    if len(records) != 1:
        raise ValueError('输入必须包含且仅包含一条 FASTA 序列')
    return records[0]


# ---------------------------------------------------------- BLAST XML dialects
# This tool requests XML2 (FORMAT_TYPE=XML2, i.e. ``blastn -outfmt 16``) but must
# also keep reading the legacy dialect (-outfmt 5 / FORMAT_TYPE=XML):
#
#   XML2   root <BlastXML2 xmlns="http://www.ncbi.nlm.nih.gov">
#          query-len | description/HitDescr/id + title | len |
#          hsps/Hsp/bit-score + identity + query-from + query-to + hit-from + hit-to + align-len
#   legacy root <BlastOutput>
#          BlastOutput_query-len | Hit_accession + Hit_def | Hit_len |
#          Hit_hsps/Hsp/Hsp_bit-score + Hsp_identity + Hsp_query-from + ...
#
# Lookup is by namespace-stripped local name so that a namespace change cannot
# silently turn real hits into an empty result set.
XML2_ROOT = 'BlastXML2'
LEGACY_XML_ROOT = 'BlastOutput'
_XML_MARKERS = frozenset((
    'query-len', 'BlastOutput_query-len', 'Iteration_query-len',
    'Hit', 'HitDescr', 'hsps', 'Hsp', 'BlastOutput_iterations', 'Iteration_hits'))


# ------------------------------------------------------------- structured parse
def _text(element):
    return element.text.strip() if element is not None and element.text else ''


def _int_text(element):
    value = _text(element)
    try:
        return int(value)
    except ValueError:
        return None


def _float_text(element):
    value = _text(element)
    try:
        return float(value)
    except ValueError:
        return None


def _local_name(tag):
    """Tag name with any ``{namespace}`` prefix removed."""
    return tag.rsplit('}', 1)[-1] if isinstance(tag, str) else ''


def _child(element, *names):
    """First DIRECT child whose local name is in ``names`` (else None).

    Never rely on ``Element`` truthiness: an element with no children (e.g.
    ``<len>1200</len>``) is falsy and ``a or b`` would silently skip it.
    """
    if element is None:
        return None
    for child in element:
        if _local_name(child.tag) in names:
            return child
    return None


def _value(element, *names):
    """Text of the first matching direct child ('' when absent)."""
    return _text(_child(element, *names))


def _int_value(element, *names):
    return _int_text(_child(element, *names))


def _find_all(element, name):
    """All descendants (self excluded) whose local name is exactly ``name``."""
    return [item for item in element.iter()
            if item is not element and _local_name(item.tag) == name]


def _int_descendant(element, *names):
    """First descendant int matching any local name, in the order of ``names``.

    Unlike ``_int_value`` this searches the whole subtree: ``query-len`` sits
    below ``Search`` in XML2 and ``Iteration_query-len`` below ``Iteration`` in
    the legacy dialect, so a direct-child lookup silently returns None.
    """
    for name in names:
        for item in element.iter():
            if _local_name(item.tag) == name:
                return _int_text(item)
    return None


def union_length(intervals):
    """Union length of half-open-free 1-based inclusive intervals (never double counts)."""
    merged = []
    for start, end in sorted((min(a, b), max(a, b)) for a, b in intervals if b > a or a == b):
        if start > end:
            continue
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return sum(end - start + 1 for start, end in merged)


def parse_search_info(text):
    """(status, rid) from an NCBI SearchInfo payload.

    NCBI answers ``FORMAT_OBJECT=SearchInfo`` with the **QBlastInfo plain text**
    block (``RID = ...`` on a line of its own, ``Status=WAITING`` / ``Status=READY``),
    not with XML -- Biopython's ``NCBIWWW`` parses exactly those textual fields.
    An XML ``BlastSearchInfo`` payload is still accepted.  A ``CMD=Put`` response
    carries only the RID, so a missing status is not an error; an unrecognised
    payload returns ``(None, None)`` so the caller can fail with exit code 3
    instead of polling a status it will never see.
    """
    text = text or ''
    tag_free = re.sub(r'<[^>]+>', ' ', text)

    def _pick(*patterns):
        for pattern in patterns:
            for candidate in (text, tag_free):
                match = re.search(pattern, candidate)
                if match:
                    return match.group(1).strip()
        return None

    status = _pick(r'<Status>\s*([A-Za-z_]+)\s*</Status>',
                   r'\bStatus\s*=\s*([A-Za-z_]+)')
    rid = _pick(r'<RID>\s*([^<\s]+)\s*</RID>',
                r'\bRID\s*=\s*([A-Za-z0-9_.-]+)')
    return (status.upper() if status else None, rid)


def subject_overlap_bases(hsps):
    """Target bases reused by more than one HSP (direction-normalised).

    Query-side union coverage is not enough: two HSPs can map to overlapping
    stretches of the SAME subject record while covering disjoint query regions
    (a repeat-collapsed reference, or duplicated query sequence).  Aggregating
    those into one precise identity would present reused target evidence as an
    independent second alignment.
    """
    intervals = [(min(hsp['hit_from'], hsp['hit_to']), max(hsp['hit_from'], hsp['hit_to']))
                 for hsp in hsps
                 if hsp.get('hit_from') is not None and hsp.get('hit_to') is not None]
    spanned = sum(end - start + 1 for start, end in intervals)
    return max(0, spanned - union_length(intervals))


def hsp_collinearity(hsps, max_span_ratio=3.0, hit_len=None, origin_margin=200):
    """Check whether a hit's HSPs describe ONE continuous, collinear alignment.

    BLAST coordinate facts, verified with real local blastn runs (do not
    re-derive): ``Hsp_query-from`` pairs with ``Hsp_hit-from``; and for a hit to
    the subject's MINUS strand ``hit_from > hit_to`` while ``query_from <
    query_to``.  Measured examples, subject 1..4000, a contiguous minus-strand
    region split into two HSPs by an insert in the query::

        contiguous minus : q 1..400 -> s 1800..1401 ; q 501..900 -> s 1400..1001
        scrambled minus  : q 1..400 -> s 1400..1001 ; q 500..900 -> s 1801..1401
        contiguous plus  : q 1..400 -> s 1001..1400 ; q 501..900 -> s 1401..1800

    So the test is NOT "lower target coordinate increases".  It is: order the
    HSPs by ``query_from`` and require ``hit_from`` to advance **in the target
    strand's own direction** -- non-decreasing on the plus strand, NON-INCREASING
    on the minus strand.  The scrambled case above is the one that must fail.

    Two separate checks, deliberately not merged:
      1. orientation -- every HSP must share one (query, target) signature; a
         genuine mix of subject strands is contradictory evidence;
      2. collinearity -- the direction-aware monotonic test above, plus the
         ``max_span_ratio`` engineering guard (an extra warning, never a
         substitute for the direction check).

    A jump that fails (2) while touching both ends of the subject is additionally
    flagged ``cross_origin_candidate``: on a circular reference that may be a
    legitimate origin-spanning arrangement, but it must NOT be accepted as a
    continuous linear alignment without explicit coordinates and structural
    evidence.
    """
    usable = [hsp for hsp in hsps
              if None not in (hsp.get('query_from'), hsp.get('query_to'),
                              hsp.get('hit_from'), hsp.get('hit_to'))]
    if len(usable) < 2:
        # A single HSP cannot be checked for collinearity, but its DERIVED strand
        # and paired coordinate must still reflect the raw coordinates: reporting
        # '+' for a subject-minus hit makes --output-json contradict the very HSP
        # it summarises.
        subject_strand = '+' if not usable else \
            ('+' if usable[0]['hit_to'] >= usable[0]['hit_from'] else '-')
        return {'consistent_direction': True, 'monotonic': True, 'span_ratio': 1.0,
                'collinear': True, 'max_span_ratio': max_span_ratio,
                'subject_strand': subject_strand, 'cross_origin_candidate': False,
                'paired_target_coordinates': [hsp['hit_from'] for hsp in usable],
                'pairing': 'query_from<->hit_from'}
    signatures = set()
    for hsp in usable:
        query_direction = '+' if hsp['query_to'] >= hsp['query_from'] else '-'
        target_direction = '+' if hsp['hit_to'] >= hsp['hit_from'] else '-'
        signatures.add((query_direction, target_direction))
    if len(signatures) != 1:
        # e.g. one HSP on the subject's plus strand and another on its minus strand
        return {'consistent_direction': False, 'monotonic': False, 'span_ratio': 0.0,
                'collinear': False, 'max_span_ratio': max_span_ratio,
                'subject_strand': 'mixed', 'cross_origin_candidate': False,
                'paired_target_coordinates': [], 'pairing': 'query_from<->hit_from'}
    _query_direction, target_direction = signatures.pop()
    ordered = sorted(usable, key=lambda hsp: hsp['query_from'])
    paired = [hsp['hit_from'] for hsp in ordered]
    if target_direction == '+':
        monotonic = all(later >= earlier for earlier, later in zip(paired, paired[1:]))
    else:
        monotonic = all(later <= earlier for earlier, later in zip(paired, paired[1:]))
    aligned = sum(abs(hsp['query_to'] - hsp['query_from']) + 1 for hsp in usable)
    coordinates = [hsp['hit_from'] for hsp in usable] + [hsp['hit_to'] for hsp in usable]
    span = max(coordinates) - min(coordinates) + 1
    ratio = (span / aligned) if aligned else 0.0
    cross_origin = False
    if not monotonic and hit_len:
        cross_origin = min(coordinates) <= origin_margin and max(coordinates) >= hit_len - origin_margin
    return {'consistent_direction': True, 'monotonic': monotonic, 'span_ratio': ratio,
            'collinear': monotonic and ratio <= max_span_ratio,
            'max_span_ratio': max_span_ratio, 'subject_strand': target_direction,
            'cross_origin_candidate': cross_origin,
            'paired_target_coordinates': paired, 'pairing': 'query_from<->hit_from'}


_HSP_FIELDS = (
    ('query-from', ('Hsp_query-from', 'query-from'), int),
    ('query-to', ('Hsp_query-to', 'query-to'), int),
    ('hit-from', ('Hsp_hit-from', 'hit-from'), int),
    ('hit-to', ('Hsp_hit-to', 'hit-to'), int),
    ('align-len', ('Hsp_align-len', 'align-len'), int),
    ('identity', ('Hsp_identity', 'identity'), int),
    ('bit-score', ('Hsp_bit-score', 'bit-score'), float),
)


def _hsp_numbers(hsp):
    """Validated numeric fields of one ``<Hsp>``; ValueError on any gap.

    All of these are required to decide whether the HSPs describe one alignment.
    Without the subject coordinates the direction, collinearity and target-reuse
    checks cannot run at all, so a half-readable hit must NOT fall through to the
    "single unverifiable HSP" fast path and become an auto-selectable candidate.
    """
    numbers, missing, invalid = {}, [], []
    for label, names, caster in _HSP_FIELDS:
        text = _value(hsp, *names)
        if text == '':
            missing.append(label)
            continue
        try:
            numbers[label] = caster(text)
        except ValueError:
            invalid.append('%s=%r' % (label, text))
    if missing or invalid:
        raise ValueError('BLAST XML 的 <Hsp> 字段不完整: %s%s%s'
                         % ('缺少 ' + ', '.join(missing) if missing else '',
                            '; ' if missing and invalid else '',
                            '非法 ' + ', '.join(invalid) if invalid else ''))
    if min(numbers['query-from'], numbers['query-to'],
           numbers['hit-from'], numbers['hit-to']) < 1:
        raise ValueError('BLAST XML 的 <Hsp> 坐标必须 >= 1: %r' % numbers)
    if numbers['align-len'] <= 0:
        raise ValueError('BLAST XML 的 <Hsp> align-len 必须 > 0: %r' % numbers)
    if not 0 <= numbers['identity'] <= numbers['align-len']:
        raise ValueError('BLAST XML 的 <Hsp> identity 必须在 0..align-len 之间: %r' % numbers)
    if not math.isfinite(numbers['bit-score']):
        raise ValueError('BLAST XML 的 <Hsp> bit-score 必须是有限数值（nan/inf 不可用）: %r'
                         % numbers)
    if numbers['bit-score'] < 0:
        raise ValueError('BLAST XML 的 <Hsp> bit-score 不能为负: %r' % numbers)
    span = abs(numbers['query-to'] - numbers['query-from']) + 1
    if span > numbers['align-len']:
        # 带 gap 时 align-len >= query 跨度；反过来说明同一份记录内部矛盾。
        raise ValueError('BLAST XML 的 <Hsp> query 跨度(%d) 大于 align-len(%d)，记录自相矛盾: %r'
                         % (span, numbers['align-len'], numbers))
    return numbers


_STRAND_FIELDS = (('query', ('Hsp_query-strand', 'query-strand'), ('Hsp_query-frame', 'query-frame')),
                  ('hit', ('Hsp_hit-strand', 'hit-strand'), ('Hsp_hit-frame', 'hit-frame')))


def _hsp_orientation(hsp, numbers):
    """坐标方向与文本方向（*-strand / *-frame）必须自洽。

    文本字段在旧实现里完全没被使用；既然回复里给了它，就让它与坐标方向对账——
    同一事实的两种表述互相矛盾属于格式故障（退出码 3）。字段缺失不算错误，
    但会记入 unchecked，使"未校验"可见而不是默认正确。
    """
    unchecked = []
    seen = False
    for side, strand_names, frame_names in _STRAND_FIELDS:
        if side == 'query':
            increasing = numbers['query-to'] >= numbers['query-from']
        else:
            increasing = numbers['hit-to'] >= numbers['hit-from']
        text = _value(hsp, *strand_names)
        if text:
            seen = True
            expected = 'plus' if increasing else 'minus'
            if text.strip().lower() != expected:
                raise ValueError('BLAST XML 的 <%s-strand>=%r 与坐标方向不符（坐标是 %s）'
                                 % (side, text, expected))
            continue
        frame = _value(hsp, *frame_names)
        if frame:
            seen = True
            try:
                sign = int(frame)
            except ValueError:
                raise ValueError('BLAST XML 的 <%s-frame> 无法解析: %r' % (side, frame))
            expected_sign = 1 if increasing else -1
            if sign != expected_sign:
                raise ValueError('BLAST XML 的 <%s-frame>=%r 与坐标方向不符（坐标方向应为 %+d）'
                                 % (side, frame, expected_sign))
    if not seen:
        unchecked.append('strand_absent')
    return unchecked


def _hsp_ranges(numbers, query_len, hit_len):
    """坐标是否落在声明的长度内。长度缺失时不报错，但记入 unchecked。"""
    unchecked = []
    if query_len:
        if max(numbers['query-from'], numbers['query-to']) > query_len:
            raise ValueError('BLAST XML 的 query 坐标(%d..%d) 超出 query-len=%d'
                             % (numbers['query-from'], numbers['query-to'], query_len))
    else:
        unchecked.append('query_len_absent')
    if hit_len is not None and hit_len <= 0:
        raise ValueError('BLAST XML 的 Hit_len 必须 > 0（声明为 %d 属格式故障）' % hit_len)
    if hit_len:
        if max(numbers['hit-from'], numbers['hit-to']) > hit_len:
            raise ValueError('BLAST XML 的 hit 坐标(%d..%d) 超出 Hit_len=%d'
                             % (numbers['hit-from'], numbers['hit-to'], hit_len))
    else:
        unchecked.append('hit_len_absent')
    return unchecked


def parse_blast_xml(text):
    """Parse an NCBI BLAST XML response into per-hit union-coverage records.

    Reads both dialects this project can meet (see the constants above): the XML2
    dialect this tool requests, and the legacy ``-outfmt 5`` dialect used by the
    older local fixtures.  Returns {'query_len': int|None, 'hits': [...]} sorted
    by strength; per hit: accession, description, identity (None when the HSPs
    cannot be reliably aggregated), aligned/spanned/redundant/subject-overlap
    bases, coverage (query-interval union), bitscore, hsp_count, blockers and the
    raw HSP list.

    A document that is XML but is not a recognisable BLAST report -- or that
    contains ``<Hit>`` elements this parser cannot read -- raises ``ValueError``:
    the caller must report a format failure (exit code 3), not ``no_match``
    (exit code 1).  A protocol mismatch is not evidence about the database.
    """
    try:
        root = ElementTree.fromstring(text)
    except ElementTree.ParseError as exc:
        raise ValueError('BLAST 结果不是有效 XML: %s' % exc)
    root_name = _local_name(root.tag)
    if root_name not in (XML2_ROOT, LEGACY_XML_ROOT):
        raise ValueError('无法识别的 BLAST XML 格式: 根元素 <%s> (支持 %s=XML2 / %s=legacy XML)'
                         % (root_name, XML2_ROOT, LEGACY_XML_ROOT))
    if not any(_local_name(item.tag) in _XML_MARKERS for item in root.iter()):
        raise ValueError('BLAST XML 结构无法识别: 根元素 <%s> 中没有任何已知字段' % root_name)

    query_len = None
    for name in ('BlastOutput_query-len', 'Iteration_query-len', 'query-len'):
        query_len = _int_descendant(root, name)
        if query_len is not None:
            break
    if query_len is not None and query_len <= 0:
        raise ValueError('BLAST XML 的 query-len 必须 > 0（声明为 %d 属格式故障）' % query_len)

    hit_elements = _find_all(root, 'Hit')
    hits = []
    for hit in hit_elements:
        accession = _value(hit, 'Hit_accession', 'Hit_id')
        description = ' '.join(_value(hit, 'Hit_def').split())
        if not accession or not description:
            # XML2 keeps them under <description><HitDescr><id>/<title>
            descr = _child(hit, 'description')
            item = _child(descr, 'HitDescr')
            if item is None:
                item = descr
            accession = accession or _value(item, 'id', 'accession')
            description = description or ' '.join(_value(item, 'title', 'def').split())
        hit_len = _int_value(hit, 'Hit_len', 'len')
        intervals, hsps, identity_total, align_total, best_bits, hsp_count = [], [], 0, 0, 0.0, 0
        unchecked = []
        for hsp in _find_all(hit, 'Hsp'):
            numbers = _hsp_numbers(hsp)
            unchecked.extend(_hsp_orientation(hsp, numbers))
            unchecked.extend(_hsp_ranges(numbers, query_len, hit_len))
            hsp_count += 1
            intervals.append((numbers['query-from'], numbers['query-to']))
            hsps.append({'query_from': numbers['query-from'], 'query_to': numbers['query-to'],
                         'hit_from': numbers['hit-from'], 'hit_to': numbers['hit-to'],
                         'identity': numbers['identity'], 'align_len': numbers['align-len'],
                         'identity_pct': numbers['identity'] / numbers['align-len'] * 100,
                         'bitscore': numbers['bit-score']})
            identity_total += numbers['identity']
            align_total += numbers['align-len']
            best_bits = max(best_bits, numbers['bit-score'])
        # Non-redundant query coverage (overlapping HSPs are counted once).  A
        # composite identity is only aggregated when the HSPs neither overlap
        # (query OR target side) nor conflict: otherwise the same bases would be
        # weighted twice, or reused target positions would be presented as two
        # independent alignments.
        aligned_bases = union_length(intervals)
        spanned_bases = sum(abs(query_to - query_from) + 1 for query_from, query_to in intervals)
        redundant_bases = max(0, spanned_bases - aligned_bases)
        target_reuse = subject_overlap_bases(hsps)
        collinearity = hsp_collinearity(hsps, hit_len=hit_len)
        blockers = []
        if any(hsp['identity'] == 0 for hsp in hsps):
            # 有长度却没有一个匹配碱基：这是"弱到无用的命中"，不是格式故障；
            # 它不得被汇总成 identity，也不得参与自动择优。
            blockers.append('zero_identity')
        if redundant_bases > 0:
            blockers.append('overlapping_hsps')
        if target_reuse > 0:
            blockers.append('subject_overlap')
        if not collinearity['collinear']:
            blockers.append('non_collinear_hsps')
        aggregate = not blockers
        per_hsp = [item['identity_pct'] for item in hsps]
        hits.append({
            'accession': accession,
            'description': description,
            'identity': ((identity_total / align_total * 100) if align_total else 0.0)
                        if aggregate else None,
            'identity_method': 'aligned-length-weighted over non-overlapping collinear HSPs'
                               if aggregate
                               else 'not_aggregated: %s -> 不汇总为单一 identity'
                                    % ', '.join(blockers),
            'identity_range': (min(per_hsp), max(per_hsp)) if per_hsp else None,
            'aligned_bases': aligned_bases,
            'spanned_bases': spanned_bases,
            'redundant_bases': redundant_bases,
            'subject_overlap_bases': target_reuse,
            'subject_overlap': target_reuse > 0,
            # metrics could not be reliably aggregated -- NOT a statement that the
            # hit is an invalid candidate
            'ambiguous_alignment': bool(blockers),
            'conflicting_alignment': not collinearity['collinear'],
            'ambiguity_reasons': blockers,
            'collinearity': collinearity,
            'hit_len': hit_len,
            'cross_origin_candidate': collinearity['cross_origin_candidate'],
            'coverage': (aligned_bases / query_len) if query_len else 0.0,
            'bitscore': best_bits,
            'hsp_count': hsp_count,
            'ranges_checked': bool(query_len) and bool(hit_len),
            'unchecked': sorted(set(unchecked)),
            'hsps': hsps,
        })
    if hit_elements and not any(item['hsp_count'] for item in hits):
        raise ValueError('BLAST XML 含 %d 个 <Hit> 但未能解析出任何 HSP; 结果结构不被支持'
                         % len(hit_elements))
    hits.sort(key=lambda item: (-(item['bitscore'] or 0.0), -((item['identity'] or 0.0))))
    return {'query_len': query_len, 'hits': hits}


def distinct_candidates(hits):
    """One candidate per accession (duplicate records are not separate species)."""
    best = {}
    for hit in hits:
        key = hit.get('accession') or hit.get('description')
        score = (hit.get('identity') or 0.0, hit.get('bitscore') or 0.0)
        current = best.get(key)
        current_score = ((current.get('identity') or 0.0), (current.get('bitscore') or 0.0)) \
            if current else None
        if current is None or score > current_score:
            best[key] = hit
    return sorted(best.values(),
                  key=lambda item: (-(item.get('bitscore') or 0.0),
                                    -(item.get('identity') or 0.0)))


# ------------------------------------------------------------------------ verdict
def _hit_parts(hit):
    """Accept dict records and the legacy (identity, description) /
    (identity, align_len, description) tuples."""
    if isinstance(hit, dict):
        return float(hit.get('identity') or 0.0), hit.get('aligned_bases'), hit.get('description')
    if len(hit) >= 3:
        return float(hit[0]), hit[1], hit[2]
    return float(hit[0]), None, hit[1]


def select_supported_hit(hits, min_identity=MIN_IDENTITY, min_coverage=None, query_len=None):
    """Return the unique best hit, or None when the evidence is not discriminating.

    Hits below ``min_identity`` -- or, when ``min_coverage`` and ``query_len`` are
    given, hits covering less than that fraction of the query -- are dropped before
    deciding.  A top hit not separated by at least 1% identity from the next is
    treated as ambiguous.
    """
    supported = []
    for hit in hits:
        identity, align_len, _ = _hit_parts(hit)
        if identity < min_identity:
            continue
        if min_coverage and query_len:
            if not align_len or align_len / query_len < min_coverage:
                continue
        supported.append(hit)
    if not supported:
        return None
    top_identity = max(_hit_parts(hit)[0] for hit in supported)
    lower = [_hit_parts(hit)[0] for hit in supported if _hit_parts(hit)[0] < top_identity]
    if lower and top_identity - max(lower) < 1.0:
        return None
    top = [hit for hit in supported if _hit_parts(hit)[0] == top_identity]
    return top[0] if len(top) == 1 else None


def interpret(identity, coverage, min_identity=MIN_IDENTITY,
              min_coverage=MIN_QUERY_COVERAGE):
    """Provisional, non-species verdict for one barcode query."""
    if identity is None or coverage is None or coverage < min_coverage:
        return {'status': 'insufficient',
                'message': 'query coverage %.2f < %.2f: 覆盖率不足, 不能得出任何条码结论'
                           % (coverage or 0.0, min_coverage)}
    if identity < min_identity:
        return {'status': 'ambiguous',
                'message': 'identity %.1f%% < %.0f%%: 参考物种可能不合适'
                           % (identity, min_identity)}
    return {'status': 'provisional_candidate',
            'message': 'identity %.1f%% / query coverage %.0f%%: 仅候选归属, 需形态、多位点与'
                       '文献证据; identity 不等于物种鉴定'
                       % (identity, coverage * 100)}


# ------------------------------------------------------------------------- remote
def _request(params, timeout=120):
    data = urllib.parse.urlencode(params).encode()
    request = urllib.request.Request(BLAST_URL, data=data,
                                     headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode('utf-8', 'replace')


def request_search(query, max_results=10):
    """CMD=Put with FORMAT_OBJECT=SearchInfo -> RID (no results yet)."""
    return _request({'CMD': 'Put', 'PROGRAM': 'blastn', 'DATABASE': 'nt',
                     'QUERY': '>query\n' + query, 'FORMAT_OBJECT': 'SearchInfo',
                     'HITLIST_SIZE': str(max_results), 'EXPECT': '1e-30',
                     'WORD_SIZE': '11', 'MEGABLAST': 'on'})


def poll_search_info(rid):
    """CMD=Get with FORMAT_OBJECT=SearchInfo -> status only (kept separate from results)."""
    return _request({'CMD': 'Get', 'FORMAT_OBJECT': 'SearchInfo', 'RID': rid}, timeout=60)


def fetch_results_xml(rid, max_results=10):
    """CMD=Get with FORMAT_TYPE=XML2 -> structured alignment blocks."""
    return _request({'CMD': 'Get', 'FORMAT_TYPE': 'XML2', 'RID': rid,
                     'ALIGNMENTS': str(max_results), 'DESCRIPTIONS': str(max_results)},
                    timeout=180)


# ---------------------------------------------------------------------------- main
def main():
    argv = sys.argv[1:]
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__)
        sys.exit(0)
    fn = argv[0]
    if '--allow-public-upload' not in argv:
        print('拒绝上传序列到公共 NCBI 服务；确认后请显式添加 --allow-public-upload',
              file=sys.stderr)
        sys.exit(2)
    if '--coords' not in argv:
        print('必须用 --coords start,end 给出 COX1 的本地定位坐标；'
              '固定 DNA 模式不能证明基因身份', file=sys.stderr)
        sys.exit(2)
    max_results = int(argv[argv.index('--max-results') + 1]) if '--max-results' in argv else 5
    out_json = argv[argv.index('--output-json') + 1] if '--output-json' in argv else None

    try:
        sequence = read_single_fasta(fn)
    except (OSError, ValueError) as exc:
        print('输入错误:', exc, file=sys.stderr)
        sys.exit(2)
    try:
        start_text, end_text = argv[argv.index('--coords') + 1].split(',')
        start, end = int(start_text), int(end_text)
    except (ValueError, IndexError):
        print('--coords 必须是 start,end 两个整数', file=sys.stderr)
        sys.exit(2)
    if start < 1 or end < start or end > len(sequence):
        print('--coords 超出输入序列范围', file=sys.stderr)
        sys.exit(2)
    query = sequence[start - 1:end]
    if len(query) < MIN_QUERY_LENGTH:
        print('查询片段至少需要 %d bp' % MIN_QUERY_LENGTH, file=sys.stderr)
        sys.exit(2)
    if len(query) > MAX_QUERY_LENGTH:
        print('拒绝上传超过 %d bp 的查询；请用 --coords 提供较小的 COX1 区域' % MAX_QUERY_LENGTH,
              file=sys.stderr)
        sys.exit(2)
    print('使用指定区域 [%d..%d] (%d bp)' % (start, end, len(query)))

    try:
        info = request_search(query, max_results)
    except (urllib.error.URLError, OSError, ValueError) as exc:
        print('NCBI 提交失败(网络/服务故障): %s' % exc, file=sys.stderr)
        sys.exit(3)
    status, rid = parse_search_info(info)
    if not rid:
        print('NCBI 未返回 RID(格式故障); 响应前 300 字符: %s' % info[:300], file=sys.stderr)
        sys.exit(3)
    print('RID: %s (轮询 SearchInfo...)' % rid)

    for _ in range(POLL_ATTEMPTS):
        time.sleep(POLL_INTERVAL)
        try:
            status, _ = parse_search_info(poll_search_info(rid))
        except (urllib.error.URLError, OSError, ValueError) as exc:
            print('轮询失败(网络故障): %s' % exc, file=sys.stderr)
            sys.exit(3)
        if status == 'READY':
            break
        if status in ('FAILED', 'ERROR', 'UNKNOWN'):
            print('NCBI 检索失败: %s' % status, file=sys.stderr)
            sys.exit(3)
    else:
        print('NCBI 检索超时(最后状态=%s)' % status, file=sys.stderr)
        sys.exit(3)

    try:
        xml_text = fetch_results_xml(rid, max_results)
    except (urllib.error.URLError, OSError, ValueError) as exc:
        print('结果下载失败(网络故障): %s' % exc, file=sys.stderr)
        sys.exit(3)
    try:
        parsed = parse_blast_xml(xml_text)
    except ValueError as exc:
        print('结果解析失败(格式故障): %s' % exc, file=sys.stderr)
        sys.exit(3)

    query_len = parsed['query_len'] or len(query)
    hits = distinct_candidates(parsed['hits'])
    results = {
        'query_length': query_len, 'blast_url_api': BLAST_URL, 'rid': rid,
        'candidates': [{
            'accession': hit['accession'], 'description': hit['description'],
            'identity': hit['identity'], 'identity_range': hit['identity_range'],
            'coverage': hit['coverage'], 'hsp_count': hit['hsp_count'],
            'bitscore': hit['bitscore'], 'hit_len': hit['hit_len'],
            'aligned_bases': hit['aligned_bases'], 'spanned_bases': hit['spanned_bases'],
            'redundant_bases': hit['redundant_bases'],
            'subject_overlap_bases': hit['subject_overlap_bases'],
            'ambiguity_reasons': hit['ambiguity_reasons'],
            'collinearity': hit['collinearity'], 'hsps': hit['hsps'],
        } for hit in hits],
        'selection': None, 'verdict': None,
        'notes': ['阈值 (%.0f%% identity / %.0f%% coverage / top-次优 1%%) 是工程启发式, '
                  '无类群特异性依据' % (MIN_IDENTITY, MIN_QUERY_COVERAGE * 100),
                  'identity 不等于物种鉴定; 候选必须人工核对原始 HSP 与坐标'],
    }

    def finish(verdict, selection=None, code=0):
        """Record the verdict, persist --output-json (if any), then exit."""
        results['verdict'] = verdict
        results['selection'] = selection
        if out_json:
            try:
                with open(out_json, 'w', encoding='utf-8') as handle:
                    json.dump(results, handle, ensure_ascii=False, indent=2)
            except OSError as exc:
                print('结果文件写入失败 %s: %s' % (out_json, exc), file=sys.stderr)
                sys.exit(2)
            print('完整候选与逐 HSP 原始记录已写入: %s' % out_json)
        sys.exit(code)

    print('\n=== 物种鉴定结果 (COX1 blastn vs nt, XML2) ===')
    print('查询长度: %d bp | 数据库: nt (NCBI BLAST URL API, 该接口不暴露库版本号)' % query_len)
    print('%-40s %8s %8s %6s' % ('命中描述', 'Ident%', 'cov%', 'HSP'))
    for hit in hits[:max_results]:
        identity = hit['identity']
        identity_text = ('%.1f%%' % identity) if identity is not None else '未汇总'
        print('%-40s %8s %7.0f%% %6d'
              % (hit['description'][:40], identity_text, hit['coverage'] * 100, hit['hsp_count']))
    ambiguous_hits = [hit for hit in hits if hit.get('ambiguity_reasons')]
    for hit in ambiguous_hits[:3]:
        print('AMBIGUOUS_ALIGNMENT: %s 的多 HSP 指标无法可靠汇总 (原因: %s; q-overlap=%dbp, '
              'subject-overlap=%dbp, target span/aligned=%.1f) —— 这不等于该 hit 不是有效候选; '
              '未汇总指标不参与自动择优; 原始 HSP 见下方明细与 --output-json'
              % (hit['accession'], ','.join(hit['ambiguity_reasons']), hit['redundant_bases'],
                 hit['subject_overlap_bases'], hit['collinearity']['span_ratio']), file=sys.stderr)
        for hsp in hit['hsps'][:20]:
            print('    HSP q%s..%s -> s%s..%s  identity=%s/%s  bits=%.1f'
                  % (hsp['query_from'], hsp['query_to'], hsp['hit_from'], hsp['hit_to'],
                     hsp['identity'], hsp['align_len'], hsp['bitscore']), file=sys.stderr)
        if len(hit['hsps']) > 20:
            print('    ... 其余 %d 个 HSP 见 --output-json' % (len(hit['hsps']) - 20),
                  file=sys.stderr)
        if hit.get('conflicting_alignment'):
            print('CONFLICTING_ALIGNMENT: %s 的 HSP 在 query/hit 方向或目标共线性上互相冲突; '
                  '不能视为一条连续可靠的对齐' % hit['accession'], file=sys.stderr)
        if hit.get('cross_origin_candidate'):
            print('CROSS_ORIGIN_CANDIDATE: %s 的 HSP 跳变同时触及目标两端; 若参考为环状, '
                  '这可能是跨原点排列, 但必须用明确坐标与结构证据确认, '
                  '不得直接按连续线性比对接受' % hit['accession'], file=sys.stderr)
    if not hits:
        print('no_match: 当前数据库与检索条件下未发现可比命中 '
              '(这只说明本次检索无候选, 不等于无近缘物种证据)', file=sys.stderr)
        finish('no_match', code=1)
    print('候选 accession 数: %d (最优 accession 不等于已完成物种鉴定; 跨 accession 候选全部保留)'
          % len(hits))

    best = select_supported_hit(hits, min_identity=MIN_IDENTITY,
                                min_coverage=MIN_QUERY_COVERAGE, query_len=query_len)
    if best is None:
        if ambiguous_hits:
            print('METRICS_NOT_AGGREGATED: 存在 %d 个候选 accession, 但其多 HSP 指标无法可靠汇总; '
                  '这**不等于**这些 hit 不是有效候选 —— 请人工核对原始 HSP 后判断'
                  % len({hit['accession'] for hit in ambiguous_hits}), file=sys.stderr)
        print('insufficient: 没有唯一、达到 %.0f%% identity 且 query coverage >= %.0f%% 的可自动汇总命中'
              % (MIN_IDENTITY, MIN_QUERY_COVERAGE * 100), file=sys.stderr)
        finish('insufficient', code=1)
    verdict = interpret(best['identity'], best['coverage'])
    print('\n最佳候选: %s (accession=%s)' % (best['description'][:60], best['accession']))
    print('  identity=%.1f%%  query coverage=%.0f%%  HSP=%d'
          % (best['identity'], best['coverage'] * 100, best['hsp_count']))
    print('  [%s] %s' % (verdict['status'], verdict['message']))
    print('\n注意: 本节不给出物种级确定结论; 阈值分档为工程启发式, 无类群特异性依据; '
          '阈值失败只过滤候选, 不构成分类学结论。')
    finish(verdict['status'], best['accession'], 0)


if __name__ == '__main__':
    main()

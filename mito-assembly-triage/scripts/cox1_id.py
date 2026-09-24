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
  python3 cox1_id.py <genome.fasta> --allow-public-upload --coords 1353,2891 [--max-results 5]
依赖: 网络 (NCBI BLAST URL API), 可选 BioPython
"""
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
    """(status, rid) from an NCBI SearchInfo XML payload ('not xml' -> (None, None))."""
    status = re.search(r'<Status>\s*([A-Za-z]+)\s*</Status>', text or '')
    rid = re.search(r'<RID>\s*([^<\s]+)\s*</RID>', text or '')
    return (status.group(1).upper() if status else None,
            rid.group(1) if rid else None)


def parse_blast_xml(text):
    """Parse an NCBI BLAST XML2 response into per-hit union-coverage records.

    Returns {'query_len': int|None, 'hits': [{accession, description, identity,
    aligned_bases, coverage, bitscore, hsp_count}]} sorted by strength.

    identity is the alignment-length-weighted identity over all HSPs of the hit;
    coverage is the query-interval union (overlapping HSPs are not double counted).
    """
    try:
        root = ElementTree.fromstring(text)
    except ElementTree.ParseError as exc:
        raise ValueError('BLAST 结果不是有效 XML: %s' % exc)
    query_len = (_int_text(root.find('.//BlastOutput_query-len'))
                 or _int_text(root.find('.//Iteration_query-len')))
    hits = []
    for hit in root.iter('Hit'):
        accession = _text(hit.find('Hit_accession')) or _text(hit.find('Hit_id'))
        description = ' '.join((_text(hit.find('Hit_def'))).split())
        intervals, identity_total, align_total, best_bits, hsp_count = [], 0, 0, 0.0, 0
        for hsp in hit.findall('./Hit_hsps/Hsp'):
            query_from = _int_text(hsp.find('Hsp_query-from'))
            query_to = _int_text(hsp.find('Hsp_query-to'))
            align_len = _int_text(hsp.find('Hsp_align-len')) or 0
            identity = _int_text(hsp.find('Hsp_identity')) or 0
            bits = _float_text(hsp.find('Hsp_bit-score')) or 0.0
            if query_from is None or query_to is None:
                continue
            hsp_count += 1
            intervals.append((query_from, query_to))
            identity_total += identity
            align_total += align_len
            best_bits = max(best_bits, bits)
        aligned_bases = union_length(intervals)
        hits.append({
            'accession': accession,
            'description': description,
            'identity': (identity_total / align_total * 100) if align_total else 0.0,
            'aligned_bases': aligned_bases,
            'coverage': (aligned_bases / query_len) if query_len else 0.0,
            'bitscore': best_bits,
            'hsp_count': hsp_count,
        })
    hits.sort(key=lambda item: (-item['bitscore'], -item['identity']))
    return {'query_len': query_len, 'hits': hits}


def distinct_candidates(hits):
    """One candidate per accession (duplicate records are not separate species)."""
    best = {}
    for hit in hits:
        key = hit.get('accession') or hit.get('description')
        current = best.get(key)
        if current is None or (hit['identity'], hit['bitscore']) > (current['identity'],
                                                                   current['bitscore']):
            best[key] = hit
    return sorted(best.values(), key=lambda item: (-item['bitscore'], -item['identity']))


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
    print('\n=== 物种鉴定结果 (COX1 blastn vs nt, XML2) ===')
    print('查询长度: %d bp | 数据库: nt (NCBI BLAST URL API, 该接口不暴露库版本号)' % query_len)
    print('%-40s %8s %8s %6s' % ('命中描述', 'Ident%', 'cov%', 'HSP'))
    for hit in hits[:max_results]:
        print('%-40s %7.1f%% %7.0f%% %6d'
              % (hit['description'][:40], hit['identity'], hit['coverage'] * 100, hit['hsp_count']))
    if not hits:
        print('no_match: 当前数据库与检索条件下未发现可比命中 '
              '(这只说明本次检索无候选, 不等于无近缘物种证据)', file=sys.stderr)
        sys.exit(1)

    best = select_supported_hit(hits, min_identity=MIN_IDENTITY,
                                min_coverage=MIN_QUERY_COVERAGE, query_len=query_len)
    if best is None:
        print('insufficient: 没有唯一、达到 %.0f%% identity 且 query coverage >= %.0f%% 的命中'
              % (MIN_IDENTITY, MIN_QUERY_COVERAGE * 100), file=sys.stderr)
        sys.exit(1)
    verdict = interpret(best['identity'], best['coverage'])
    print('\n最佳候选: %s (accession=%s)' % (best['description'][:60], best['accession']))
    print('  identity=%.1f%%  query coverage=%.0f%%  HSP=%d'
          % (best['identity'], best['coverage'] * 100, best['hsp_count']))
    print('  [%s] %s' % (verdict['status'], verdict['message']))
    print('\n注意: 本节不给出物种级确定结论; 阈值分档为工程启发式, 无类群特异性依据; '
          '阈值失败只过滤候选, 不构成分类学结论。')


if __name__ == '__main__':
    main()

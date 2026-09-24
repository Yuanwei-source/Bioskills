#!/usr/bin/env python3
"""
COX1 物种鉴定: 提取 COX1 (或给定区域), 提交 NCBI blastn, RID 轮询, 解析最佳命中

必须显式给出 --coords (本脚本不自动识别 COX1) 与 --allow-public-upload。
结果报告 identity、**query coverage** 与多个候选; 不做物种级确定结论:
identity/coverage 只能给出候选归属, 阈值分档属工程启发式, 无类群特异性依据。

用法: python3 cox1_id.py <genome.fasta> --allow-public-upload --coords 1353,2891 [--max-results 5]
依赖: 网络 (NCBI eutils/blast), 可选 BioPython
"""
import sys, re, time, urllib.request, urllib.parse


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

MIN_IDENTITY = 85.0
MIN_QUERY_COVERAGE = 0.8


def _hit_parts(hit):
    """Accept (identity, description) and (identity, align_len, description)."""
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


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    fn = sys.argv[1]
    max_res = 5
    if '--max-results' in sys.argv:
        max_res = int(sys.argv[sys.argv.index('--max-results')+1])
    if '--allow-public-upload' not in sys.argv:
        print('拒绝上传序列到公共 NCBI 服务；确认后请显式添加 --allow-public-upload', file=sys.stderr)
        sys.exit(2)

    try:
        seq = read_single_fasta(fn)
    except (OSError, ValueError) as exc:
        print('输入错误:', exc, file=sys.stderr)
        sys.exit(2)
    coords = None
    if '--coords' in sys.argv:
        try:
            a, b = sys.argv[sys.argv.index('--coords')+1].split(',')
            start, end = int(a), int(b)
        except (ValueError, IndexError):
            print('--coords 必须是 start,end 两个整数', file=sys.stderr)
            sys.exit(2)
        if start < 1 or end < start or end > len(seq):
            print('--coords 超出输入序列范围', file=sys.stderr)
            sys.exit(2)
        coords = (start, end)
        query = seq[start-1:end]
        print('使用指定区域 [%d..%d] (%d bp)' % (start, end, len(query)))
    else:
        print('未提供 COX1 的可靠定位。固定 DNA 模式不能证明基因身份；请先用本地注释/tblastn 定位后以 --coords 提取。', file=sys.stderr)
        sys.exit(2)
    if len(query) < 400:
        print('查询片段至少需要 400 bp', file=sys.stderr)
        sys.exit(2)
    if not query or len(query) > 5000:
        print('拒绝上传超过 5000 bp 的查询；请用 --coords 提供较小的 COX1 区域', file=sys.stderr)
        sys.exit(2)

    # NCBI blast
    url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
    params = {
        'CMD': 'Put', 'PROGRAM': 'blastn', 'DATABASE': 'nt',
        'QUERY': '>query\n' + query,
        'FORMAT_TYPE': 'Text', 'HITLIST_SIZE': str(max_res),
        'EXPECT': '1e-30', 'WORD_SIZE': '11', 'MEGABLAST': 'on',
    }
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(url, data=data, headers={'User-Agent': 'Mozilla/5.0'})
    resp = urllib.request.urlopen(req, timeout=120).read().decode()
    m = re.search(r'RID = (\S+)', resp)
    if not m:
        print('提交失败:', resp[:300]); sys.exit(1)
    rid = m.group(1)
    print('RID: %s (轮询中...)' % rid)

    status = None
    for i in range(40):
        time.sleep(12)
        q = urllib.parse.urlencode({'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'Text',
                                    'ALIGNMENTS': str(max_res), 'DESCRIPTIONS': str(max_res)})
        r = urllib.request.Request(url + '?' + q, headers={'User-Agent': 'Mozilla/5.0'})
        out = urllib.request.urlopen(r, timeout=60).read().decode()
        status_match = re.search(r'Status=([A-Z]+)', out)
        status = status_match.group(1) if status_match else None
        if status == 'READY':
            break
        if status in ('FAILED', 'ERROR'):
            print('NCBI 检索失败:', status, file=sys.stderr)
            sys.exit(1)
    if status != 'READY':
        print('NCBI 检索超时或未返回 READY', file=sys.stderr)
        sys.exit(1)

    # 解析
    print('\n=== 物种鉴定结果 (COX1 blastn vs nt) ===')
    print('查询长度: %d bp | 数据库: nt (NCBI BLAST URL API 提交, 该接口不暴露库版本号)' % len(query))
    hits = re.findall(r'>(\S+?) (.+?)\nLength=\d+\n\n Score = (\d+) bits.*?Identities = (\d+)/(\d+)', out, re.S)
    query_len = len(query)
    parsed_hits = [(int(idn) / int(tot) * 100, int(tot), desc.strip())
                   for _, desc, _, idn, tot in hits]
    print('%-45s %8s %8s' % ('命中描述', 'Ident%', 'cov%'))
    for identity, align_len, description in parsed_hits[:max_res]:
        print('%-45s %7.1f%% %7.0f%%' % (description[:45], identity, align_len / query_len * 100))
    best_hit = select_supported_hit(parsed_hits, min_identity=MIN_IDENTITY,
                                    min_coverage=MIN_QUERY_COVERAGE, query_len=query_len)
    if best_hit is None:
        print('没有唯一、达到 %.0f%% identity 且 query coverage >= %.0f%% 的 COI 命中; 状态: insufficient'
              % (MIN_IDENTITY, MIN_QUERY_COVERAGE * 100), file=sys.stderr)
        sys.exit(1)
    identity, align_len, description = best_hit
    coverage = align_len / query_len
    verdict = interpret(identity, coverage)
    print('\n最佳候选: %s' % description[:60])
    print('  identity=%.1f%%  query coverage=%.0f%%' % (identity, coverage * 100))
    print('  [%s] %s' % (verdict['status'], verdict['message']))
    print('\n注意: 本节不给出物种级确定结论; 阈值分档为工程启发式, 无类群特异性依据。')

if __name__ == '__main__':
    main()

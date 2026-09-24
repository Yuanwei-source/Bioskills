#!/usr/bin/env python3
"""
COX1 物种鉴定: 提取 COX1 (或给定区域), 提交 NCBI blastn, RID 轮询, 解析最佳命中
用法: python3 cox1_id.py <genome.fasta> --allow-public-upload [--coords 1353,2891] [--max-results 5]
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

def select_supported_hit(hits, min_identity=85.0):
    supported = [(identity, description) for identity, description in hits if identity >= min_identity]
    if not supported:
        return None
    top_identity = max(identity for identity, _ in supported)
    lower = [identity for identity, _ in supported if identity < top_identity]
    if lower and top_identity - max(lower) < 1.0:
        return None
    top = [hit for hit in supported if hit[0] == top_identity]
    return top[0] if len(top) == 1 else None


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
        marker = r'GGAGG[ACGT]{4}GGAACTG'
        m = re.search(marker, seq)
        if m:
            start = m.start() + 1
            query = seq[start-1:start-1+1400]
            print('自动定位 COX1: 位置 %d, 长度 %d bp' % (start, len(query)))
        else:
            print('未能自动定位 COX1, 请用 --coords 指定', file=sys.stderr)
            sys.exit(1)
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
    print('%-45s %8s %6s' % ('物种', 'Score', 'Ident%'))
    print('-' * 65)
    hits = re.findall(r'>(\S+?) (.+?)\nLength=\d+\n\n Score = (\d+) bits.*?Identities = (\d+)/(\d+)', out, re.S)
    parsed_hits = [(int(idn) / int(tot) * 100, desc.strip()) for _, desc, _, idn, tot in hits if int(tot) >= 400]
    best_hit = select_supported_hit(parsed_hits)
    if best_hit is None:
        print('没有唯一且达到 85% identity 的 COI 命中', file=sys.stderr)
        sys.exit(1)
    identity, description = best_hit
    print('%-45s %5.1f%%' % (description[:45], identity))

    print('\n判读标准:')
    print('  >97%  同种 (intraspecific)')
    print('  85-97% 同属不同种')
    print('  <85%  不同属/科 → 参考物种需更换')

if __name__ == '__main__':
    main()

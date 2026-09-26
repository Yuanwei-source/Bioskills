#!/usr/bin/env python3
"""MITOS2 result.gff + result.fas + result.faa  ->  GenBank, for annotation QC.

WHY THIS EXISTS
  MITOS2's ``runmitos.py`` (driven by ``scripts/run_mitos2.sh``) writes
  result.gff / result.fas / result.faa / result.mitos / result.bed but **no
  GenBank**, while ``scripts/annot_check.py`` accepts only GenBank, and the skill
  ships no converter.  MitoFinder can emit GenBank but requires a reference
  GenBank as input.  Without this bridge the documented annotation -> QC loop
  cannot be closed at all (real-data validation finding F-03).

WHAT IT IS, AND WHAT IT IS NOT
  * It is a **format bridge**: every coordinate, strand, phase and gene name comes
    from MITOS2's own output files.  It adds no gene, no boundary and no base.
  * It is an **annotation-QC input**, with provenance recorded on every record.
    Its output is **not** an NCBI-submission-quality annotation, and passing
    ``annot_check.py`` on it would still not make it one.
  * It never infers sample bases from a reference and never completes a missing
    gene: what MITOS2 did not annotate stays absent.
  * It never infers physical circularity.  MITOS2's circular *mode* is an input
    handling choice, not evidence of a closed molecule; ``--topology`` is an
    explicit declaration (default ``linear``) and is reported as such.
  * Fragmented genes keep their original names and separate coordinates
    (``nad6_0`` / ``nad6_1`` / ``nad6_2`` are NOT merged): merging would hide the
    fragmentation that the QC step exists to surface.
  * ``/codon_start`` has an explicit source.  It is recovered by translating frames
    1/2/3 of MITOS2's own extracted CDS and matching MITOS2's own protein; when no
    protein is available (or the best frame is not unique) the value is reported as
    **uncertain** instead of being silently assumed to be 1.
  * Partial CDS markers are **not** invented: MITOS2 does not export ``<``/``>``
    partiality, so every location is written as a complete span and that limitation
    is stated in the output.

Cross-validation performed before writing (all must agree, else fail closed):
  GFF coordinates/strand  ==  result.fas header  ==  result.fas sequence
  and  result.fas sequence == the genome slice it claims to be.

Usage:
  mitos2_to_genbank.py <result.gff> <result.fas> <result.faa> <out.gb>
                       [--genome <input.fasta>] [--topology linear|circular]
                       [--locus <name>]
"""
import re
import sys
from pathlib import Path


# 前置门禁：本步骤所需依赖（唯一清单来源 config/dependencies.json）
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from _deps import require_stage  # noqa: E402
require_stage('mitos2_bridge', __file__)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

CDS_PATTERN = re.compile(r'^(nad\d+l?|cox[123]|cob|cytb|atp[68])(_\d+)?$')
TRNA_PATTERN = re.compile(r'^(trn[A-Za-z]\d?)\((?:anticodon:)?([acgtu]{3})\)$')
RRNA_NAMES = ('rrnl', 'rrns')
KNOWN_NON_GENE_TYPES = ('region', 'gene', 'exon', 'ncRNA_gene', 'tRNA', 'rRNA',
                        'CDS', 'origin_of_replication', 'rep_origin', 'D-loop')
PRODUCTS = {
    'nad1': 'NADH dehydrogenase subunit 1', 'nad2': 'NADH dehydrogenase subunit 2',
    'nad3': 'NADH dehydrogenase subunit 3', 'nad4': 'NADH dehydrogenase subunit 4',
    'nad4l': 'NADH dehydrogenase subunit 4L', 'nad5': 'NADH dehydrogenase subunit 5',
    'nad6': 'NADH dehydrogenase subunit 6', 'cox1': 'cytochrome c oxidase subunit 1',
    'cox2': 'cytochrome c oxidase subunit 2', 'cox3': 'cytochrome c oxidase subunit 3',
    'cob': 'cytochrome b', 'cytb': 'cytochrome b',
    'atp6': 'ATP synthase F0 subunit 6', 'atp8': 'ATP synthase F0 subunit 8',
    'rrnl': '16S ribosomal RNA', 'rrns': '12S ribosomal RNA',
}
PROVENANCE = 'MITOS2-derived; converted for annotation QC by mitos2_to_genbank.py'


def fail(message, code=2):
    print('错误: %s' % message, file=sys.stderr)
    sys.exit(code)


def parse_gff(path):
    """Minimal GFF3 reader -> list of dicts.  Only the fields MITOS2 uses."""
    rows, region = [], None
    try:
        text = Path(path).read_text(encoding='utf-8')
    except OSError as exc:
        fail('无法读取 GFF %s: %s' % (path, exc))
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 9:
            fail('GFF 行少于 9 列: %r' % line[:120])
        try:
            start, end = int(parts[3]), int(parts[4])
        except ValueError:
            fail('GFF 坐标不是整数: %r' % line[:120])
        attrs = {}
        for item in parts[8].split(';'):
            if '=' in item:
                key, _, value = item.partition('=')
                attrs[key.strip()] = value.strip()
        row = {'seqid': parts[0], 'source': parts[1], 'type': parts[2],
               'start': start, 'end': end, 'strand': parts[6], 'phase': parts[7],
               'attrs': attrs, 'name': attrs.get('Name') or attrs.get('gene_id') or ''}
        if parts[2] == 'region':
            region = row
        else:
            rows.append(row)
    return rows, region


def parse_mitos_fasta(path):
    """MITOS2 result.fas / result.faa -> {name: (start, end, strand, seq)}."""
    entries = {}
    try:
        lines = Path(path).read_text(encoding='utf-8').splitlines()
    except OSError as exc:
        fail('无法读取 %s: %s' % (path, exc))
    name, seq, meta = None, [], None
    def flush():
        if name is not None:
            entries[name] = (meta[0], meta[1], meta[2], ''.join(seq))
    for line in lines:
        if line.startswith('>'):
            flush()
            fields = [p.strip() for p in line[1:].split(';')]
            if len(fields) >= 4:
                coord, strand, name = fields[1], fields[2], fields[3]
            elif len(fields) == 3:
                coord, strand, name = fields
            else:
                fail('MITOS2 FASTA-like header 无法解析: %r' % line[:120])
            try:
                start, end = (int(x) for x in coord.split('-'))
            except ValueError:
                fail('MITOS2 header 坐标无法解析: %r' % line[:120])
            if strand not in ('+', '-'):
                fail('MITOS2 header 链向非法: %r' % line[:120])
            meta, seq = (start, end, strand), []
        else:
            seq.append(line.strip())
    flush()
    return entries


def best_protein_frame(nt, protein):
    """(offset, matched) best reproducing MITOS2's own protein; offset None if unsure."""
    if not protein:
        return None, 0
    scores = []
    for offset in range(3):
        trimmed = nt[offset:len(nt) - ((len(nt) - offset) % 3)]   # avoid partial-codon warnings
        translated = str(Seq(trimmed).translate(table=5)) if trimmed else ''
        scores.append(sum(1 for a, b in zip(translated, protein) if a == b and a != '*'))
    best = max(scores)
    if best <= 0 or scores.count(best) > 1:
        return None, best
    return scores.index(best), best


def main():
    args = [a for a in sys.argv[1:]]
    if len(args) < 4 or args[0] in ('-h', '--help'):
        print(__doc__)
        sys.exit(0)
    gff_path, fas_path, faa_path, out_path = (Path(args[0]), Path(args[1]),
                                              Path(args[2]), Path(args[3]))
    genome_path, topology, locus = None, 'linear', None
    rest = args[4:]
    while rest:
        flag = rest.pop(0)
        if flag == '--genome' and rest:
            genome_path = Path(rest.pop(0))
        elif flag == '--topology' and rest:
            topology = rest.pop(0).lower()
        elif flag == '--locus' and rest:
            locus = rest.pop(0)
        else:
            fail('未知或缺少取值的参数: %s' % flag)

    if topology not in ('linear', 'circular'):
        fail('--topology 只能是 linear 或 circular，收到 %r' % topology)

    if not gff_path.is_file():
        fail('result.gff 不存在: %s' % gff_path)
    if not fas_path.is_file():
        fail('result.fas 不存在: %s' % fas_path)
    if not faa_path.is_file():
        fail('result.faa 不存在: %s' % faa_path)

    # sequence source: the annotated genome itself (MITOS2 deletes its intermediate
    # sequence.fas when it finishes, so the caller's assembly is authoritative)
    if genome_path is None:
        for candidate in (gff_path.parent / 'sequence.fas', gff_path.parent / 'sequence.fas-0'):
            if candidate.is_file():
                genome_path = candidate
                break
    if genome_path is None or not genome_path.is_file():
        fail('未找到基因组序列: 请用 --genome <input.fasta> 指定（MITOS2 会删除中间文件 sequence.fas）')
    genome = ''.join(line.strip() for line in genome_path.read_text(encoding='utf-8').splitlines()
                     if not line.startswith('>'))
    genome = re.sub(r'[^A-Za-z]', '', genome).upper()
    if not genome:
        fail('基因组序列为空: %s' % genome_path)

    rows, region = parse_gff(gff_path)
    fas = parse_mitos_fasta(fas_path)
    faa = parse_mitos_fasta(faa_path)
    if region is not None and region['end'] > len(genome):
        fail('GFF region 声明 %d bp，但序列只有 %d bp' % (region['end'], len(genome)))

    features, unknown_types, warnings, notes = [], [], [], []
    for row in rows:
        name = row['name']
        if row['type'] not in KNOWN_NON_GENE_TYPES:
            unknown_types.append('%s(%s)' % (row['type'], name or '?'))
            continue
        if row['type'] in ('gene', 'CDS'):
            if not CDS_PATTERN.match(name.lower()):
                # `gene` rows that are not CDS genes (e.g. OH_*) are reported, not used
                if name and not name.lower().startswith(('oh_', 'oh')):
                    unknown_types.append('%s(%s)' % (row['type'], name))
                continue
            kind = 'CDS'
        elif row['type'] == 'tRNA':
            if not TRNA_PATTERN.match(name):
                warnings.append('tRNA 名称格式未识别: %r（保留原始名称）' % name)
            kind = 'tRNA'
        elif row['type'] == 'rRNA':
            kind = 'rRNA'
        else:
            continue

        start, end, strand = row['start'], row['end'], row['strand']
        if start < 1 or end < start or end > len(genome):
            fail('%s %s 的坐标 %d-%d 超出序列范围 1-%d'
                 % (kind, name or '?', start, end, len(genome)))
        if strand not in ('+', '-'):
            fail('%s %s 的链向非法: %r' % (kind, name, strand))

        sliced = genome[start - 1:end]
        if strand == '-':
            sliced = str(Seq(sliced).reverse_complement())
        if name in fas:
            f_start, f_end, f_strand, f_seq = fas[name]
            if (f_start, f_end, f_strand) != (start, end, strand):
                fail('%s %s 在 result.gff 与 result.fas 中的坐标/链不一致: %s vs %s'
                     % (kind, name, (start, end, strand), (f_start, f_end, f_strand)))
            if not f_seq:
                fail('%s %s 在 result.fas 中序列为空' % (kind, name))
            if f_seq.upper().replace('U', 'T') != sliced.replace('U', 'T'):
                fail('%s %s 的 result.fas 序列与基因组对应区间不一致（长度 %d vs %d）'
                     % (kind, name, len(f_seq), len(sliced)))
        else:
            notes.append('%s %s 在 result.fas 中缺失，序列取自基因组区间' % (kind, name))

        qualifiers = {'gene': [name]}
        if kind == 'CDS':
            base = CDS_PATTERN.match(name.lower()).group(1)
            offset, matched = best_protein_frame(sliced, faa.get(name, (0, 0, '+', ''))[3])
            if offset is None:
                codon_start = 1
                warnings.append('%s %s 的 codon_start 不确定（result.faa 缺失或读框不唯一，'
                                '暂按 1 写出并在 note 中标注）' % (kind, name))
                qualifiers['note'] = ['codon_start uncertain: %s' % PROVENANCE]
            else:
                codon_start = offset + 1
            qualifiers['product'] = [PRODUCTS.get(base, name)]
            qualifiers['transl_table'] = ['5']
            qualifiers['codon_start'] = [str(codon_start)]
        elif kind == 'tRNA':
            match = TRNA_PATTERN.match(name)
            if match:
                anticodon = match.group(2).upper().replace('U', 'T')
                qualifiers['product'] = ['tRNA-%s' % name[3:4].upper()]
                qualifiers['anticodon'] = ['(pos:%s,aa:%s)' % (anticodon, name[3:4].upper())]
            else:
                qualifiers['product'] = ['tRNA']
        else:
            qualifiers['product'] = [PRODUCTS.get(name.lower(), 'ribosomal RNA')]
        features.append(SeqFeature(
            FeatureLocation(start - 1, end, strand=1 if strand == '+' else -1),
            type=kind, qualifiers=qualifiers))

    if unknown_types:
        fail('存在未识别的 feature 类型/名称，拒绝静默丢弃: %s' % ', '.join(sorted(set(unknown_types))))
    if not features:
        fail('GFF 中没有任何可转换的 CDS/tRNA/rRNA feature')

    record = SeqRecord(Seq(genome), id=(locus or 'MITOS2_annotation'),
                       name=(locus or 'MITOS2_annotation'),
                       description='MITOS2 annotation, converted for annotation QC',
                       annotations={'molecule_type': 'DNA', 'topology': topology,
                                    'data_file_division': 'INV'})
    record.features.append(SeqFeature(
        FeatureLocation(0, len(genome), strand=1), type='source',
        qualifiers={'organism': ['not determined by this conversion'],
                    'mol_type': ['genomic DNA'],
                    'note': [PROVENANCE] + notes + warnings}))
    record.features.extend(sorted(features, key=lambda f: int(f.location.start)))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + '.tmp')
    SeqIO.write(record, tmp, 'genbank')
    tmp.replace(out_path)

    print('来源: %s + %s + %s; 序列: %s (%d bp)'
          % (gff_path.name, fas_path.name, faa_path.name, genome_path.name, len(genome)))
    print('转换: %d 个 feature（未新增任何注释；碎片化基因保持原名称与独立坐标）'
          % (len(features)))
    print('topology 声明: %s（MITOS2 的 circular 模式不是物理环化证据；如需声明 circular '
          '必须显式 --topology circular，且 --require-circular 也只校验声明）' % topology)
    print('partial: MITOS2 不导出 < / > 部分标记，因此所有 location 均按完整区间写出（保留该不确定性）')
    for message in warnings:
        print('warn: %s' % message)
    print('注意: 本输出是**注释质检输入**，带来源可追溯；annot_check.py 通过它也不代表'
          '达到 NCBI 提交质量，也不代表 reads 支持该序列。')
    print('已写出: %s' % out_path)


if __name__ == '__main__':
    main()

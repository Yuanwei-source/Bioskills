#!/usr/bin/env python3
"""MITOS2 result.gff + result.fas + result.faa  ->  GenBank, for annotation QC.

WHY THIS EXISTS
  MITOS2's ``runmitos.py`` (driven by ``scripts/run_mitos2.sh``) writes
  result.gff / result.fas / result.faa / result.mitos / result.bed but **no
  GenBank**, while ``scripts/annot_check.py`` accepts only GenBank, so this converter bridges their file formats.  MitoFinder can emit GenBank but requires a reference
  GenBank as input.  Without this bridge the documented annotation -> QC loop
  cannot be closed at all (real-data validation finding F-03).

WHAT IT IS, AND WHAT IT IS NOT
  * It is a **format bridge**: every coordinate, strand and gene name comes
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
    protein fully matches (or the full-match frame is not unique) and no CDS phase is supplied the value is reported as
    **uncertain** instead of being silently assumed to be 1.
  * Partial CDS markers are **not** invented: MITOS2 does not export ``<``/``>``
    partiality, so every location is written as a complete span and that limitation
    is stated in the output.

Cross-validation performed before writing when the feature is present in FAS (otherwise noted):
  GFF coordinates/strand  ==  result.fas header  ==  result.fas sequence
  and  result.fas sequence == the genome slice it claims to be.

Usage:
  mitos2_to_genbank.py <result.gff> <result.fas> <result.faa> <out.gb>
                       [--genome <input.fasta>] [--topology linear|circular]
                       --table <confirmed genetic code> [--locus <name>]
"""
import re
import sys
import argparse
import io
from pathlib import Path


# 前置门禁：本步骤所需依赖（唯一清单来源 config/dependencies.json）
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from _deps import require_stage  # noqa: E402
from _evidence_io import read_single_fasta, protect_output, atomic_write_text, InputArgumentParser

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
            if name in entries:
                fail('MITOS2 FASTA 中重复的 feature 名称: %s' % name)
            entries[name] = (meta[0], meta[1], meta[2], ''.join(seq))
    for line in lines:
        if line.startswith('>'):
            flush()
            fields = [p.strip() for p in line[1:].split(';')]
            # MITOS2 headers may include a semicolon-delimited sequence ID
            # (e.g. ``scaffold;len=...;topology=...;start-end;+;gene``).
            # The coordinate/strand/name triplet is the stable suffix.
            if len(fields) < 3:
                fail('MITOS2 FASTA-like header 无法解析: %r' % line[:120])
            coord, strand, name = fields[-3:]
            try:
                start, end = (int(x) for x in coord.split('-'))
            except ValueError:
                fail('MITOS2 header 坐标无法解析: %r' % line[:120])
            if strand not in ('+', '-'):
                fail('MITOS2 header 链向非法: %r' % line[:120])
            meta, seq = (start, end, strand), []
        else:
            if line.strip() and name is None:
                fail('MITOS2 FASTA 序列出现在 header 之前')
            seq.append(line.strip())
    flush()
    return entries


def best_protein_frame(nt, protein, table):
    """Only a unique full protein match establishes the frame; never a best guess.

    One terminal stop may be omitted. Initiator M is permitted only for a start
    codon in the explicit code table. Unknown residues/internal stops fail closed.
    """
    from Bio.Seq import Seq
    from Bio.Data import CodonTable
    if not protein:
        return None, 0
    protein = protein.upper()
    if protein.endswith('*'):
        protein = protein[:-1]
    if not protein or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', protein):
        return None, 0
    matches = []
    starts = CodonTable.unambiguous_dna_by_id[table].start_codons
    for offset in range(3):
        trimmed = nt[offset:len(nt) - ((len(nt) - offset) % 3)]   # avoid partial-codon warnings
        translated = str(Seq(trimmed).translate(table=table)) if trimmed else ''
        if translated.endswith('*'):
            translated = translated[:-1]
        initiated = 'M' + translated[1:] if translated and trimmed[:3].upper() in starts else translated
        if protein in (translated, initiated):
            matches.append(offset)
    return (matches[0], len(protein)) if len(matches) == 1 else (None, 0)


def main():
    parser = InputArgumentParser(description=__doc__)
    for key in ('gff', 'fas', 'faa', 'output'):
        parser.add_argument(key, type=Path)
    parser.add_argument('--genome', type=Path)
    parser.add_argument('--topology', choices=('linear', 'circular'), default='linear')
    parser.add_argument('--locus')
    parser.add_argument('--organism', help='调用者确认的样本学名；不提供时明确写 not determined')
    parser.add_argument('--taxid', type=int, help='调用者确认的 NCBI TaxID；不会自动推断')
    parser.add_argument('--table', type=int, required=True, help='与 MITOS2 注释使用的密码表一致')
    args = parser.parse_args()
    require_stage('mitos2_bridge', __file__)
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Data import CodonTable
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    if args.table not in CodonTable.unambiguous_dna_by_id:
        fail('未知遗传密码表: %s' % args.table)
    gff_path, fas_path, faa_path, out_path = args.gff, args.fas, args.faa, args.output
    genome_path, topology, locus = args.genome, args.topology, args.locus
    if args.taxid is not None and args.taxid <= 0:
        fail('--taxid 必须是正整数')

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
    genome_id, genome = read_single_fasta(genome_path)
    protect_output(out_path, [genome_path, gff_path, fas_path, faa_path])

    rows, region = parse_gff(gff_path)
    fas = parse_mitos_fasta(fas_path)
    faa = parse_mitos_fasta(faa_path)
    if any(row['seqid'] != genome_id for row in rows + ([region] if region else [])):
        fail('GFF seqid 与输入 FASTA 标识不一致')
    if region is not None and region['end'] > len(genome):
        fail('GFF region 声明 %d bp，但序列只有 %d bp' % (region['end'], len(genome)))

    features, unknown_types, warnings, notes = [], [], [], []
    seen = set()
    # Some GFF writers emit a parent gene and its CDS. Preserve one CDS, giving
    # the explicit CDS row (including phase) priority over the parent gene row.
    cds_rows = {(r['name'], r['start'], r['end'], r['strand']) for r in rows if r['type'] == 'CDS'}
    for row in rows:
        name = row['name']
        key = (name, row['start'], row['end'], row['strand'])
        if row['type'] == 'gene' and key in cds_rows:
            continue
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
        if (kind, key) in seen:
            fail('重复 GFF feature: %s %s' % (kind, name))
        seen.add((kind, key))

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
            protein = faa.get(name)
            if protein and protein[:3] != (start, end, strand):
                fail('%s 的 result.faa 坐标/链与 GFF 不一致' % name)
            offset, matched = best_protein_frame(sliced, protein[3] if protein else '', args.table)
            phase = row['phase'] if row['type'] == 'CDS' else '.'
            if phase not in ('.', '0', '1', '2'):
                fail('%s 的 CDS phase 非法: %s' % (name, phase))
            if phase != '.':
                if protein and (offset is None or offset != int(phase)):
                    fail('%s 的 CDS phase 与蛋白重现读框冲突' % name)
                # Explicit phase is source data, not inferred from protein identity.
                offset = int(phase)
                qualifiers['note'] = ['codon_start source: GFF CDS phase']
            if offset is None:
                warnings.append('%s %s 的 codon_start 不确定（无唯一完整蛋白匹配且无 CDS phase）；未写出 codon_start'
                                % (kind, name))
                qualifiers['note'] = ['codon_start uncertain: %s' % PROVENANCE]
            else:
                qualifiers['codon_start'] = [str(offset + 1)]
                qualifiers.setdefault('note', ['codon_start source: unique full MITOS2 protein match'])
            qualifiers['product'] = [PRODUCTS.get(base, name)]
            qualifiers['transl_table'] = [str(args.table)]
        elif kind == 'tRNA':
            match = TRNA_PATTERN.match(name)
            if match:
                anticodon = match.group(2).upper().replace('U', 'T')
                qualifiers['product'] = ['tRNA-%s' % name[3:4].upper()]
                # MITOS name exports sequence, not the anticodon's genomic position.
                # Do not manufacture INSDC pos; QC can use this exact provenance note.
                qualifiers['note'] = ['MITOS2 anticodon sequence: %s; position not exported' % anticodon]
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
    source_qualifiers = {'organism': [args.organism or 'not determined by this conversion'],
                         'mol_type': ['genomic DNA'],
                         'note': [PROVENANCE] + notes + warnings}
    if args.organism:
        source_qualifiers['note'].append('organism supplied by caller: %s' % args.organism)
    if args.taxid is not None:
        source_qualifiers['db_xref'] = ['taxon:%d' % args.taxid]
        source_qualifiers['note'].append('TaxID supplied by caller: %d' % args.taxid)
    record.features.append(SeqFeature(
        FeatureLocation(0, len(genome), strand=1), type='source',
        qualifiers=source_qualifiers))
    record.features.extend(sorted(features, key=lambda f: int(f.location.start)))

    buffer = io.StringIO()
    SeqIO.write(record, buffer, 'genbank')
    atomic_write_text(out_path, buffer.getvalue())

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
    try:
        main()
    except (OSError, ValueError) as exc:
        fail(str(exc))

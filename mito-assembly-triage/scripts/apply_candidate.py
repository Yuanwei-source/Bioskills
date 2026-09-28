#!/usr/bin/env python3
"""Apply explicit, human-selected sequence substitutions or simple CDS boundary edits.

This utility does not infer edits or establish biological validity. Every edit must
state the exact expected source value, rationale and a pre-existing evidence reference.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import re
import sys
from pathlib import Path

from _deps import require_stage
from _evidence_io import InputArgumentParser, atomic_write_text, protect_output


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_edits(path):
    def reject_constant(value):
        raise ValueError('非法 JSON 数值: %s' % value)
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('重复 JSON 字段: %s' % key)
            result[key] = value
        return result
    with open(path, encoding='utf-8') as handle:
        data = json.load(handle, object_pairs_hook=pairs)
    if not isinstance(data, dict) or set(data) != {'changes'} or not isinstance(data['changes'], list) or not data['changes']:
        raise ValueError('编辑文件须为包含非空 changes 数组的 JSON object')
    return data['changes']


def text_record(record_id, sequence):
    return '>%s\n%s\n' % (record_id, '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60)))


def base_candidate(source, edits):
    from Bio import SeqIO
    from Bio.Seq import Seq
    records = list(SeqIO.parse(str(source), 'fasta'))
    if not records:
        raise ValueError('FASTA 至少需要一个序列记录')
    by_id = {}
    sequences = {}
    for record in records:
        if not record.id or record.id in by_id:
            raise ValueError('FASTA 记录 ID 必须非空且唯一')
        sequence = str(record.seq).upper()
        if not sequence or any(base not in 'ACGTRYSWKMBDHVN' for base in sequence):
            raise ValueError('FASTA 序列含空记录或非法 IUPAC DNA 碱基: %s' % record.id)
        by_id[record.id] = record
        sequences[record.id] = list(sequence)
    seen = set()
    applied = []
    for change in edits:
        if not isinstance(change, dict) or set(change) != {'sequence_id', 'position', 'old', 'new', 'evidence', 'rationale'}:
            raise ValueError('base change 字段须为 sequence_id/position/old/new/evidence/rationale')
        sequence_id = change['sequence_id']
        if not isinstance(sequence_id, str) or not isinstance(change['old'], str) or not isinstance(change['new'], str):
            raise ValueError('sequence_id/old/new 必须为字符串')
        pos, old, new = change['position'], change['old'].upper(), change['new'].upper()
        if sequence_id not in sequences:
            raise ValueError('未知 FASTA sequence_id: %s' % sequence_id)
        chars = sequences[sequence_id]
        if not isinstance(pos, int) or isinstance(pos, bool) or not 1 <= pos <= len(chars):
            raise ValueError('碱基坐标须为该序列内有效的 1-based 位置')
        key = (sequence_id, pos)
        if key in seen:
            raise ValueError('同一候选不能重复编辑 %s:%d' % key)
        if len(old) != 1 or len(new) != 1 or old not in 'ACGTRYSWKMBDHVN' or new not in 'ACGT' or old == new:
            raise ValueError('旧碱基须为 IUPAC DNA 字符，新碱基须为不同的明确 A/C/G/T')
        if chars[pos-1] != old:
            raise ValueError('%s:%d 的原碱基为 %s，与声明的 %s 不符' % (sequence_id, pos, chars[pos-1], old))
        if not isinstance(change['evidence'], str) or not change['evidence'].strip():
            raise ValueError('每项碱基编辑必须引用已有证据')
        if not isinstance(change['rationale'], str) or not change['rationale'].strip():
            raise ValueError('每项碱基编辑必须说明理由')
        chars[pos-1] = new
        seen.add(key)
        applied.append(dict(sequence_id=sequence_id, position=pos, old=old, new=new,
                            evidence=change['evidence'], rationale=change['rationale']))
    output = io.StringIO()
    for record in records:
        record.seq = Seq(''.join(sequences[record.id]))
    SeqIO.write(records, output, 'fasta')
    return output.getvalue(), applied


def annotation_candidate(source, edits, table):
    from Bio import SeqIO
    from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition
    try:
        record = SeqIO.read(str(source), 'genbank')
    except (ValueError, OSError) as exc:
        raise ValueError('需要且只接受一个 GenBank 记录: %s' % exc)
    applied = []
    seen = set()
    for change in edits:
        required = {'selector', 'old', 'new', 'old_codon_start', 'new_codon_start', 'evidence', 'rationale'}
        if not isinstance(change, dict) or set(change) != required:
            raise ValueError('CDS 边界变更须含 selector/old/new/old_codon_start/new_codon_start/evidence/rationale')
        selector = change['selector']
        if not isinstance(selector, dict) or set(selector) != {'qualifier', 'value'} or selector['qualifier'] not in ('locus_tag', 'gene'):
            raise ValueError('selector 仅支持 locus_tag 或 gene 精确匹配')
        matches = []
        for feature in record.features:
            values = feature.qualifiers.get(selector['qualifier'], [])
            if feature.type == 'CDS' and selector['value'] in values:
                matches.append(feature)
        if len(matches) != 1:
            raise ValueError('selector 必须唯一命中一个 CDS，实际命中 %d 个' % len(matches))
        feature = matches[0]
        if isinstance(feature.location, CompoundLocation) or not isinstance(feature.location, FeatureLocation):
            raise ValueError('不自动编辑复合或不确定边界 CDS')
        old, new = change['old'], change['new']
        if not isinstance(old, dict) or set(old) != {'start', 'end', 'strand'} or not isinstance(new, dict) or set(new) != {'start', 'end', 'strand'}:
            raise ValueError('old/new 须含 1-based inclusive start/end 与 strand')
        if not isinstance(feature.location.start, ExactPosition) or not isinstance(feature.location.end, ExactPosition):
            raise ValueError('不自动编辑已有模糊 CDS 边界')
        observed = {'start': int(feature.location.start)+1, 'end': int(feature.location.end), 'strand': feature.location.strand}
        if observed != old:
            raise ValueError('CDS 当前边界 %s 与声明的 old %s 不符' % (observed, old))
        start, end, strand = new['start'], new['end'], new['strand']
        if any(not isinstance(x, int) or isinstance(x, bool) for x in (start, end, strand)) or strand not in (-1, 1):
            raise ValueError('新边界须为整数且 strand 为 +1 或 -1')
        if start < 1 or end < start or end > len(record.seq):
            raise ValueError('新 CDS 边界超出序列范围')
        if not isinstance(change['evidence'], str) or not change['evidence'].strip() or not isinstance(change['rationale'], str) or not change['rationale'].strip():
            raise ValueError('每项边界编辑必须引用证据并说明理由')
        key = (selector['qualifier'], selector['value'])
        if key in seen:
            raise ValueError('同一 CDS 不能在一个候选中重复编辑')
        old_codon_start = change['old_codon_start']
        new_codon_start = change['new_codon_start']
        if any(not isinstance(value, int) or isinstance(value, bool) or value not in (1, 2, 3)
               for value in (old_codon_start, new_codon_start)):
            raise ValueError('old_codon_start/new_codon_start 必须显式为 1、2 或 3')
        qualifier_values = feature.qualifiers.get('codon_start', ['1'])
        if len(qualifier_values) != 1 or not qualifier_values[0].isdigit() or int(qualifier_values[0]) != old_codon_start:
            raise ValueError('old_codon_start 与 GenBank 当前 codon_start 不符')
        before_translation = str(feature.translate(record.seq, table=table, cds=False))
        feature.location = FeatureLocation(start-1, end, strand=strand)
        feature.qualifiers['codon_start'] = [str(new_codon_start)]
        after_translation = str(feature.translate(record.seq, table=table, cds=False))
        seen.add(key)
        applied.append(dict(selector=selector, old=old, new=new, old_codon_start=old_codon_start,
                            new_codon_start=new_codon_start, table=table,
                            translation_before=before_translation, translation_after=after_translation,
                            evidence=change['evidence'], rationale=change['rationale']))
    output = io.StringIO()
    SeqIO.write(record, output, 'genbank')
    return output.getvalue(), applied


def qualifier_candidate(source, edits):
    """Change exact GenBank qualifier values without touching sequence or coordinates."""
    from Bio import SeqIO
    try:
        record = SeqIO.read(str(source), 'genbank')
    except (ValueError, OSError) as exc:
        raise ValueError('需要且只接受一个 GenBank 记录: %s' % exc) from exc
    applied = []
    seen = set()
    for change in edits:
        required = {'selector', 'qualifier', 'old', 'new', 'evidence', 'rationale'}
        if not isinstance(change, dict) or set(change) != required:
            raise ValueError('注释字段编辑须含 selector/qualifier/old/new/evidence/rationale')
        selector = change['selector']
        if not isinstance(selector, dict) or selector not in ({'type': 'source'},) and (
                set(selector) != {'type', 'qualifier', 'value'} or
                selector.get('type') not in ('CDS', 'tRNA', 'rRNA', 'gene') or
                selector.get('qualifier') not in ('gene', 'locus_tag', 'product') or
                not isinstance(selector.get('value'), str) or not selector['value']):
            raise ValueError('selector 须为唯一 source，或用 type + gene/locus_tag/product 精确定位特征')
        key = change['qualifier']
        if not isinstance(key, str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', key):
            raise ValueError('qualifier 名称无效')
        old, new = change['old'], change['new']
        if any(not isinstance(values, list) or any(not isinstance(v, str) or not v for v in values)
               for values in (old, new)) or old == new:
            raise ValueError('old/new 须为不同的非空字符串数组；空数组表示不存在该字段')
        if any(not isinstance(change[name], str) or not change[name].strip()
               for name in ('evidence', 'rationale')):
            raise ValueError('每项注释字段编辑必须引用证据并说明理由')
        matches = [feature for feature in record.features if feature.type == selector['type'] and
                   (selector['type'] == 'source' or selector['value'] in feature.qualifiers.get(selector['qualifier'], []))]
        if len(matches) != 1:
            raise ValueError('selector 必须唯一命中一个特征，实际命中 %d 个' % len(matches))
        feature = matches[0]
        feature_key = (id(feature), key)
        if feature_key in seen:
            raise ValueError('同一特征的同一 qualifier 不能重复编辑')
        observed = feature.qualifiers.get(key, [])
        if observed != old:
            raise ValueError('%s 当前 %s=%s 与声明的 old=%s 不符' % (selector, key, observed, old))
        if new:
            feature.qualifiers[key] = list(new)
        else:
            feature.qualifiers.pop(key, None)
        if selector['type'] == 'source' and key == 'organism':
            if len(new) != 1:
                raise ValueError('source /organism 必须恰有一个新值')
            record.annotations['organism'] = new[0]
            record.annotations['source'] = new[0]
        seen.add(feature_key)
        applied.append(dict(selector=selector, qualifier=key, old=old, new=new,
                            evidence=change['evidence'], rationale=change['rationale']))
    output = io.StringIO()
    SeqIO.write(record, output, 'genbank')
    return output.getvalue(), applied


def write_candidate(args):
    require_stage('apply_candidate', __file__)
    source, edits_path, output, manifest = map(Path, (args.input, args.edits, args.output, args.manifest))
    protect_output(edits_path, (source,))
    source_hash_before, edits_hash_before = digest(source), digest(edits_path)
    changes = read_edits(edits_path)
    protect_output(output, (source, edits_path, manifest))
    protect_output(manifest, (source, edits_path, output))
    if output.exists() or manifest.exists():
        raise ValueError('候选或 manifest 输出已存在，拒绝覆盖')
    if args.kind == 'base':
        candidate, applied = base_candidate(source, changes)
    elif args.kind == 'qualifier':
        candidate, applied = qualifier_candidate(source, changes)
    else:
        candidate, applied = annotation_candidate(source, changes, args.table)
    source_hash_after, edits_hash_after = digest(source), digest(edits_path)
    if source_hash_before != source_hash_after or edits_hash_before != edits_hash_after:
        raise ValueError('输入或编辑清单在生成候选期间发生变化')
    source_hash = source_hash_after
    atomic_write_text(output, candidate)
    manifest_data = {
        'schema_version': '1.0', 'kind': args.kind,
        'created_at': datetime.now(timezone.utc).isoformat(),
        'source': {'path': str(source.resolve()), 'sha256': source_hash},
        'edit_list': {'path': str(edits_path.resolve()), 'sha256': edits_hash_after},
        'candidate': {'path': str(output.resolve()), 'sha256': hashlib.sha256(candidate.encode('utf-8')).hexdigest()},
        'changes': applied,
        'scientific_status': 'candidate_only_not_validated'
    }
    atomic_write_text(manifest, json.dumps(manifest_data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print('已生成候选（尚未验证）: %s' % output)
    print('变更清单: %s' % manifest)


def main(argv=None):
    parser = InputArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=('base', 'cds-boundary', 'qualifier'))
    parser.add_argument('input')
    parser.add_argument('edits')
    parser.add_argument('output')
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--table', type=int)
    args = parser.parse_args(argv)
    if args.kind == 'cds-boundary' and args.table is None:
        parser.error('cds-boundary 必须显式指定确认过的 --table')
    try:
        write_candidate(args)
        return 0
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        print('ERROR: %s' % exc, file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

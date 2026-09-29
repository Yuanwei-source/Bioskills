#!/usr/bin/env python3
"""Run independent animal-mitogenome annotation candidates on assembly FASTAs."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from env_check import probe_command, probe_python_module  # noqa: E402


def has_only_background_tasks(outdir):
    intermediate = outdir / 'intermediate'
    tasks = intermediate / 'tasks'
    try:
        return (outdir.is_dir() and set(p.name for p in outdir.iterdir()) == {'intermediate'}
                and intermediate.is_dir() and set(p.name for p in intermediate.iterdir()) == {'tasks'}
                and tasks.is_dir())
    except OSError:
        return False


def fasta_info(path):
    ids, records = [], 0
    try:
        with open(path, encoding='ascii', errors='replace') as handle:
            for line in handle:
                if line.startswith('>'):
                    header = line[1:].strip()
                    if not header:
                        raise ValueError('FASTA 有空序列标识: %s' % path)
                    ids.append(header.split()[0])
                    records += 1
    except OSError as exc:
        raise ValueError('无法读取 FASTA: %s (%s)' % (path, exc)) from exc
    if not records:
        raise ValueError('FASTA 没有序列记录: %s' % path)
    return ids


def fasta_records(path):
    records, header, sequence = [], None, []
    with open(path, encoding='ascii', errors='strict') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(sequence).upper()))
                header, sequence = line[1:].strip(), []
                if not header:
                    raise ValueError('FASTA 有空序列标识: %s' % path)
            else:
                if header is None:
                    raise ValueError('FASTA 序列出现在首个标识之前: %s' % path)
                sequence.append(''.join(line.split()))
    if header is not None:
        records.append((header, ''.join(sequence).upper()))
    if any(not sequence for _, sequence in records):
        raise ValueError('FASTA 含空序列记录: %s' % path)
    return records


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def mitoz_safe_fasta(source, destination, topology):
    """Change only the FASTA identifier for MitoZ's 13-character limit."""
    chunks, header = [], None
    with open(source, encoding='ascii', errors='replace') as handle:
        for line in handle:
            if line.startswith('>'):
                header = line[1:].strip()
            elif line.strip():
                chunks.append(''.join(line.split()).upper())
    if not header or not chunks:
        raise ValueError('无法为 MitoZ 准备单记录 FASTA: %s' % source)
    sequence = ''.join(chunks)
    alias = 'mt01'
    with open(destination, 'w', encoding='ascii') as handle:
        handle.write('>%s topology=%s\n' % (alias, topology))
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start:start + 80] + '\n')
    return {'original_header': header, 'mitoz_id': alias,
            'sequence_sha256': hashlib.sha256(sequence.encode('ascii')).hexdigest()}


def run(command, cwd, env, log):
    start = time.monotonic()
    try:
        with open(log, 'w', encoding='utf-8') as handle:
            proc = subprocess.run(command, cwd=cwd, env=env, stdout=handle,
                                  stderr=subprocess.STDOUT, check=False)
        return {'status': 'succeeded' if proc.returncode == 0 else 'failed',
                'exit_code': proc.returncode, 'command': command, 'log': str(log),
                'elapsed_seconds': round(time.monotonic() - start, 2)}
    except OSError as exc:
        Path(log).write_text('无法启动工具：%s\n' % exc, encoding='utf-8')
        return {'status': 'failed_to_start', 'exit_code': None,
                'command': command, 'log': str(log), 'error': str(exc),
                'elapsed_seconds': round(time.monotonic() - start, 2)}


def main(argv=None):
    parser = argparse.ArgumentParser(description='对组装候选运行 MITOS2 与 MitoZ 注释，保留各自原始输出供裁决')
    candidate_source = parser.add_mutually_exclusive_group(required=True)
    candidate_source.add_argument('--candidate-fasta', action='append',
                                  help='一个组装候选 FASTA；可重复指定，必须是单记录 FASTA')
    candidate_source.add_argument('--assembly-manifest',
                                  help='run_illumina_candidates.py 生成的 assembly_candidates.json')
    parser.add_argument('--outdir', required=True, help='新建输出目录')
    parser.add_argument('--taxon', required=True, help='用户确认的学名')
    parser.add_argument('--table', required=True, type=int, help='已核实/暂定的 NCBI 遗传密码表编号')
    parser.add_argument('--table-status', choices=('provisional', 'confirmed'), default='provisional')
    parser.add_argument('--clade', default='Arthropoda',
                        choices=('Chordata', 'Arthropoda', 'Echinodermata', 'Annelida-segmented-worms',
                                 'Bryozoa', 'Mollusca', 'Nematoda', 'Nemertea-ribbon-worms', 'Porifera-sponges'))
    parser.add_argument('--topology', choices=('linear', 'circular'), default='linear',
                        help='注释器输入声明；circular 仅是记录声明，不证明物理闭合')
    parser.add_argument('--split-multi-records', action='store_true',
                        help='逐条注释多记录 FASTA；片段注释不等于可独立采用的完整组装')
    parser.add_argument('--tools', default='mitos2,mitoz', help='逗号分隔：mitos2,mitoz')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--mitoz', help='MitoZ 可执行文件路径')
    parser.add_argument('--mitos2-python', help='可 import mitos 的 Python；默认按 CONDA_ROOT 自动发现')
    args = parser.parse_args(argv)

    try:
        if args.table < 1 or args.threads < 1:
            raise ValueError('--table 与 --threads 必须为正整数')
        candidate_metadata = {}
        skipped_candidates = []
        if args.assembly_manifest:
            assembly_manifest = Path(args.assembly_manifest).resolve()
            try:
                assembly_data = json.loads(assembly_manifest.read_text(encoding='utf-8'))
            except (OSError, ValueError) as exc:
                raise ValueError('组装候选 manifest 无法读取: %s (%s)' % (assembly_manifest, exc)) from exc
            if assembly_data.get('format') != 'mito-illumina-assembly-candidates-1':
                raise ValueError('不支持的组装候选 manifest 格式')
            succeeded_tools = {item.get('tool') for item in assembly_data.get('tools', [])
                               if item.get('status') == 'succeeded'}
            candidate_values = []
            for entry in assembly_data.get('candidates', []):
                raw_path = entry.get('path')
                if not raw_path or entry.get('producer') not in succeeded_tools:
                    skipped_candidates.append({'path': raw_path, 'reason': '候选缺少成功工具来源记录'})
                    continue
                path = Path(raw_path).resolve()
                if not path.is_file():
                    raise ValueError('manifest 候选 FASTA 不存在: %s' % path)
                observed = sha256_file(path)
                if not entry.get('sha256') or observed != entry['sha256']:
                    raise ValueError('候选 FASTA 与 assembly manifest 的 SHA256 不符: %s' % path)
                candidate_values.append(str(path))
                candidate_metadata[str(path)] = dict(entry)
        else:
            candidate_values = args.candidate_fasta or []
        candidates = list(dict.fromkeys(Path(item).resolve() for item in candidate_values))
        if not candidates:
            raise ValueError('组装候选 manifest 中没有 FASTA 候选')
        eligible = []
        for path in candidates:
            if not path.is_file():
                raise ValueError('候选 FASTA 不存在: %s' % path)
            ids = fasta_info(path)
            if len(ids) != 1 and not args.split_multi_records:
                skipped_candidates.append({'path': str(path), 'reason':
                                           '当前注释器适配器逐序列运行；需要显式传 --split-multi-records 才会逐条生成片段注释'})
                continue
            eligible.append(path)
        candidates = eligible
        if not candidates:
            raise ValueError('没有适合注释的候选 FASTA；多记录输入可显式使用 --split-multi-records 逐条注释')
        if len({str(item) for item in candidates}) != len(candidates):
            raise ValueError('候选 FASTA 路径重复')
        outdir = Path(args.outdir).resolve()
        precreated_for_background = outdir.exists() and has_only_background_tasks(outdir)
        if outdir.exists() and not precreated_for_background:
            raise ValueError('为保护旧结果，输出目录必须不存在: %s' % outdir)
        requested = [item.strip().lower() for item in args.tools.split(',') if item.strip()]
        if not requested or len(set(requested)) != len(requested) or set(requested) - {'mitos2', 'mitoz'}:
            raise ValueError('--tools 只能是 mitos2,mitoz，且不能重复')
        if not precreated_for_background:
            outdir.mkdir(parents=True)
        candidate_provenance = {}
        expanded_candidates = []
        for source in candidates:
            records = fasta_records(source)
            if len(records) == 1:
                expanded_candidates.append(source)
                if str(source) in candidate_metadata:
                    candidate_provenance[str(source)] = {
                        'source_candidate': candidate_metadata[str(source)],
                        'source_record_id': records[0][0].split()[0],
                        'source_fasta': str(source),
                    }
                continue
            if not args.split_multi_records:
                skipped_candidates.append({'path': str(source), 'reason':
                                           '多记录候选逐条注释需要显式启用 --split-multi-records；片段注释不代表可独立采用的完整组装'})
                continue
            split_root = outdir / 'candidate_inputs'
            split_root.mkdir(parents=True, exist_ok=True)
            source_hash = sha256_file(source)
            for record_index, (header, sequence) in enumerate(records, 1):
                split_path = split_root / ('%s-%s-record_%04d.fasta' % (source.stem, source_hash[:10], record_index))
                with split_path.open('w', encoding='ascii') as handle:
                    handle.write('>%s\n' % header.split()[0])
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start:start + 80] + '\n')
                expanded_candidates.append(split_path)
                candidate_provenance[str(split_path)] = {
                    'source_candidate': candidate_metadata.get(str(source)),
                    'source_fasta': str(source), 'source_fasta_sha256': source_hash,
                    'source_record_id': header.split()[0], 'record_index': record_index,
                    'source_record_count': len(records),
                    'interpretation_limit': '逐条 contig 注释，仅作基因线索；不能证明完整/可独立采用的组装',
                }
        candidates = expanded_candidates
        if not candidates:
            raise ValueError('没有可注释的单记录序列；若输入是多记录 FASTA，请传 --split-multi-records')
        runs, produced = [], []
        mitoz_exe = None
        if 'mitoz' in requested:
            mitoz_exe = args.mitoz or os.environ.get('MITOZ') or os.environ.get('MITOZ_EXE')
            if not (mitoz_exe and Path(mitoz_exe).is_file() and os.access(mitoz_exe, os.X_OK)):
                mitoz_exe = probe_command('mitoz')[0]
        mitos_python = args.mitos2_python or os.environ.get('MITOS2_PY')
        mitos_version = None
        if 'mitos2' in requested and not mitos_python:
            mitos_python, mitos_version = probe_python_module('mitos', distribution='mitos', search_conda=True)
        if mitos_python and mitos_version is None:
            _, mitos_version = probe_python_module('mitos', interpreter=mitos_python, distribution='mitos')

        for candidate in candidates:
            candidate_hash = hashlib.sha256(candidate.read_bytes()).hexdigest()
            candidate_key = (re.sub(r'[^A-Za-z0-9_.-]+', '_', candidate.stem).strip('._') or 'candidate')
            path_hash = hashlib.sha256(str(candidate).encode('utf-8')).hexdigest()[:8]
            candidate_key += '-' + candidate_hash[:10] + '-' + path_hash
            for tool in requested:
                workdir = outdir / candidate_key / tool
                workdir.mkdir(parents=True)
                if tool == 'mitos2':
                    if not mitos_python:
                        runs.append({'candidate': str(candidate), 'tool': tool, 'status': 'unavailable',
                                     'reason': 'CONDA_ROOT 中未找到可 import mitos 的 Python'})
                        continue
                    if not os.environ.get('MITOS2_REFDIR'):
                        runs.append({'candidate': str(candidate), 'tool': tool, 'status': 'not_run',
                                     'reason': 'MITOS2_REFDIR 未配置；MITOS2 数据库不能从 Python 包推断'})
                        continue
                    env = os.environ.copy()
                    env['MITOS2_PY'] = str(Path(mitos_python).resolve())
                    env['PATH'] = str(Path(mitos_python).resolve().parent) + os.pathsep + env.get('PATH', '')
                    for var in ('MITOS2_REFDIR', 'MITOS2_REFSEQVER', 'MITOS2_EXTRA_PATH',
                                'MITOS2_RSCRIPT', 'CONDA_ROOT'):
                        if os.environ.get(var):
                            env[var] = os.environ[var]
                    extra_bins = [Path(item).resolve().parent for name in ('cmsearch', 'RNAplot')
                                  if (item := probe_command(name)[0])]
                    if extra_bins:
                        env['MITOS2_EXTRA_PATH'] = os.pathsep.join(
                            [str(path) for path in extra_bins] + ([env['MITOS2_EXTRA_PATH']] if env.get('MITOS2_EXTRA_PATH') else []))
                    command = ['bash', str(ROOT / 'scripts/run_fasta_annotation.sh'),
                               '--fasta', str(candidate), '--outdir', str(workdir / 'route'),
                               '--table', str(args.table), '--table-status', args.table_status,
                               '--organism', args.taxon, '--topology', 'linear']
                    result = run(command, ROOT, env, workdir / 'tool.log')
                    result['version'] = mitos_version
                else:
                    if not mitoz_exe:
                        runs.append({'candidate': str(candidate), 'tool': tool, 'status': 'unavailable',
                                     'reason': '当前 PATH、CONDA_ROOT/envs/*/bin 和显式路径均未发现 MitoZ'})
                        continue
                    mitoz_input = workdir / 'mitoz_input.fasta'
                    alias_map = mitoz_safe_fasta(candidate, mitoz_input, args.topology)
                    env = os.environ.copy()
                    env['PATH'] = str(Path(mitoz_exe).resolve().parent) + os.pathsep + env.get('PATH', '')
                    command = [mitoz_exe, 'annotate', '--workdir', str(workdir), '--outprefix', 'annotation',
                               '--thread_number', str(args.threads), '--fastafiles', str(mitoz_input),
                               '--species_name', args.taxon, '--genetic_code', str(args.table), '--clade', args.clade]
                    result = run(command, workdir, env, workdir / 'tool.log')
                    log = Path(result['log']).read_text(encoding='utf-8', errors='replace')
                    version_match = re.search(r'(?i)MitoZ[^\n]{0,60}?([0-9]+(?:\.[0-9]+){1,2})', log)
                    result['version'] = version_match.group(1) if version_match else None
                result.update({'candidate': str(candidate), 'candidate_sha256': candidate_hash,
                               'tool': tool, 'workdir': str(workdir)})
                if str(candidate) in candidate_provenance:
                    result['candidate_provenance'] = candidate_provenance[str(candidate)]
                if tool == 'mitoz':
                    result['sequence_id_alias'] = alias_map
                annotation_files = sorted(str(p.resolve()) for p in workdir.rglob('*')
                                          if p.is_file() and p.suffix.lower() in
                                          ('.gb', '.gbk', '.gbf', '.genbank', '.gff', '.gff3'))
                result['annotation_outputs'] = annotation_files
                if result.get('status') == 'succeeded' and not annotation_files:
                    result['status'] = 'succeeded_no_annotation_files'
                    result['reason'] = '工具退出码为 0，但未发现 GenBank/GFF 注释输出'
                produced.extend({'candidate': str(candidate), 'tool': tool, 'path': item} for item in annotation_files)
                runs.append(result)
        manifest = {
            'format': 'mito-annotation-candidates-1',
            'created_at': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
            'taxon': args.taxon, 'genetic_code_table': args.table,
            'genetic_code_status': args.table_status, 'clade': args.clade,
            'assembly_manifest': str(Path(args.assembly_manifest).resolve()) if args.assembly_manifest else None,
            'candidates': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                            'provenance': candidate_provenance.get(str(path))}
                           for path in candidates],
            'skipped_candidates': skipped_candidates,
            'runs': runs, 'outputs': produced, 'annotation_decision': 'not_selected',
            'limitations': [
                '各注释器输出是预测候选；不能以工具多数票决定基因身份或边界。',
                '必须统一坐标/链方向/密码表后逐基因比较，并结合参考同源、结构与 reads 证据。',
                '注释工具通过自身质检不证明组装真实性、样本碱基正确或序列闭环。',
                'MITOS2 数据库路径须显式配置；MitoZ 输入仅将 FASTA ID 改为短别名，序列不变且原 ID 映射保存在 manifest。',
                '线性/环状拓扑是注释输入声明，不构成物理闭环证据。'],
        }
        path = outdir / 'annotation_candidates.json'
        path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        succeeded = sum(item.get('status') == 'succeeded' for item in runs)
        print('注释候选清单：%s' % path)
        print('成功工具运行：%d/%d；注释输出文件：%d；自动选择：无' % (succeeded, len(runs), len(produced)))
        return 0 if succeeded else 2
    except (OSError, ValueError) as exc:
        print('错误：%s' % exc, file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())

#!/usr/bin/env python3
"""Generate and register animal mitogenome candidates from paired Illumina reads.

The script is an adapter/orchestrator, not an assembler selector. Each tool
runs in its own output directory and the manifest preserves command, version,
exit status, log, and discovered FASTA candidates. It never edits input reads
or chooses a final assembly.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from env_check import probe_command  # noqa: E402

TOOL_IDS = ('getorganelle', 'mitoflex', 'novoplasty', 'mitoz', 'mitofinder')
FASTA_SUFFIXES = ('.fa', '.fasta', '.fna')


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def has_only_background_tasks(outdir):
    """Allow task_manager.py's task records to pre-create the sample root."""
    intermediate = outdir / 'intermediate'
    tasks = intermediate / 'tasks'
    try:
        return (outdir.is_dir() and set(p.name for p in outdir.iterdir()) == {'intermediate'}
                and intermediate.is_dir() and set(p.name for p in intermediate.iterdir()) == {'tasks'}
                and tasks.is_dir())
    except OSError:
        return False


def iter_fastq(path):
    """Yield validated records from standard four-line FASTQ, streaming input."""
    opener = gzip.open if str(path).endswith(('.gz', '.bgz')) else open
    try:
        with opener(path, 'rt', encoding='ascii', errors='strict') as handle:
            number = 0
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                number += 1
                if not sequence or not plus or not quality:
                    raise ValueError('FASTQ 记录 %d 不完整: %s' % (number, path))
                header, sequence = header.rstrip('\r\n'), sequence.rstrip('\r\n')
                plus, quality = plus.rstrip('\r\n'), quality.rstrip('\r\n')
                if not header.startswith('@') or not plus.startswith('+'):
                    raise ValueError('FASTQ 第 %d 条记录格式错误: %s' % (number, path))
                if len(sequence) != len(quality):
                    raise ValueError('FASTQ 第 %d 条序列与质量长度不等: %s' % (number, path))
                yield header[1:].split()[0], len(sequence)
    except (OSError, EOFError, UnicodeError) as exc:
        raise ValueError('FASTQ 无法完整读取或压缩流损坏 %s: %s' % (path, exc)) from exc


def normalized_read_id(header_id):
    return re.sub(r'/[12]$', '', header_id)


def validate_fastq_pair(r1, r2):
    count, first_length = 0, None
    left, right = iter_fastq(r1), iter_fastq(r2)
    while True:
        try:
            record1 = next(left)
        except StopIteration:
            record1 = None
        try:
            record2 = next(right)
        except StopIteration:
            record2 = None
        if record1 is None or record2 is None:
            if record1 != record2:
                raise ValueError('R1/R2 FASTQ 记录数不同；输入可能不完整或不是成对文件')
            break
        count += 1
        if normalized_read_id(record1[0]) != normalized_read_id(record2[0]):
            raise ValueError('R1/R2 第 %d 条 read ID 不配对: %s vs %s' % (count, record1[0], record2[0]))
        if first_length is None:
            first_length = record1[1]
        if record1[1] != record2[1]:
            raise ValueError('R1/R2 第 %d 条 read 长度不同: %d vs %d' % (count, record1[1], record2[1]))
    if not count:
        raise ValueError('FASTQ 文件为空')
    return {'read_count': count, 'read_length_first_record': first_length}


def parse_fasta(path):
    records, name, length = [], None, 0
    try:
        with open(path, encoding='ascii', errors='replace') as handle:
            for line in handle:
                if line.startswith('>'):
                    if name is not None:
                        records.append({'id': name, 'length': length})
                    header = line[1:].strip()
                    name, length = (header.split(maxsplit=1)[0], 0) if header else (None, 0)
                else:
                    length += sum(not char.isspace() for char in line)
        if name is not None:
            records.append({'id': name, 'length': length})
    except OSError:
        return []
    return records


def resolve(name, explicit=None, env_names=()):
    for value in (explicit, *(os.environ.get(env) for env in env_names)):
        if value and Path(value).is_file() and os.access(value, os.X_OK):
            return str(Path(value).resolve())
    found, _ = probe_command(name)
    return found


def resolve_script(name, explicit=None, env_names=()):
    for value in (explicit, *(os.environ.get(env) for env in env_names)):
        if value and Path(value).is_file():
            return str(Path(value).resolve())
    conda_root = os.environ.get('CONDA_ROOT')
    if conda_root:
        for env_bin in sorted(Path(conda_root).glob('envs/*/bin')):
            candidate = env_bin / name
            if candidate.is_file():
                return str(candidate.resolve())
    return None


def mitoflex_python(explicit=None):
    choices = [explicit or os.environ.get('MITOFLEX_PYTHON')]
    conda_root = os.environ.get('CONDA_ROOT')
    if conda_root:
        for env_bin in sorted(Path(conda_root).glob('envs/*/bin')):
            choices.extend([str(env_bin / 'python'), str(env_bin / 'python3')])
    modules = ('numpy', 'pandas', 'ete3', 'Bio', 'psutil')
    roots = sorted({str(path.parents[1]) for path in Path(conda_root).glob('envs/*/lib/python*/site-packages/ete3/__init__.py')}) if conda_root else []
    for candidate in dict.fromkeys(item for item in choices if item):
        if not os.path.isfile(candidate) or not os.access(candidate, os.X_OK):
            continue
        try:
            direct = subprocess.run(
                [candidate, '-B', '-c', 'import numpy,pandas,ete3,Bio,psutil'],
                cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=15)
            if direct.returncode == 0:
                return str(Path(candidate).resolve()), []
            probe = ('import importlib.util,json; print(json.dumps([m for m in %r '
                     'if importlib.util.find_spec(m) is None]))') % (modules,)
            result = subprocess.run([candidate, '-B', '-c', probe], cwd=ROOT,
                                    capture_output=True, text=True, timeout=15)
            missing = json.loads(result.stdout.strip()) if result.returncode == 0 else list(modules)
            # Only supplement ete3 when this interpreter does not have it at
            # all. A broken, incompatible ete3 already on sys.path must not be
            # hidden by putting another package copy in front of native libs.
            if any(module != 'ete3' for module in missing):
                continue
            for root in roots:
                trial = ('import sys; sys.path.append(%r); '
                         'import numpy,pandas,ete3,Bio,psutil') % root
                tested = subprocess.run([candidate, '-B', '-c', trial], cwd=ROOT,
                                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                         timeout=15)
                if tested.returncode == 0:
                    return str(Path(candidate).resolve()), [root]
        except (OSError, subprocess.SubprocessError):
            continue
    return None, []


def normalize_gzip_members(source, destination):
    """Rewrite concatenated gzip members as one standards-compliant member.

    MitoFlex's Rust filter uses a single-member GzDecoder, so plain `cat` of
    multiple .fastq.gz files otherwise makes it silently ignore later lanes.
    Python's gzip reader consumes every member; GzipFile writes one member.
    """
    source, destination = Path(source), Path(destination)
    with gzip.open(str(source), 'rb') as reader, destination.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', compresslevel=1, mtime=0) as writer:
            shutil.copyfileobj(reader, writer, length=8 * 1024 * 1024)
    return {'source': str(source.resolve()), 'path': str(destination.resolve()),
            'bytes': destination.stat().st_size, 'sha256': sha256(destination),
            'normalization': 'all concatenated gzip members rewritten into one gzip member'}


def with_tool_path(executable, extra_paths=()):
    env = os.environ.copy()
    paths = list(extra_paths)
    launch = executable if isinstance(executable, (list, tuple)) else [executable]
    paths.extend(str(Path(item).resolve().parent) for item in launch
                 if item and os.path.isfile(item))
    # Composite tools may be installed in separate conda environments. Keep
    # the launcher environment first, then expose other env bins for helpers.
    # The selected executable itself is always called by its absolute path.
    conda_root = env.get('CONDA_ROOT')
    if conda_root:
        paths.extend(str(item.resolve()) for item in sorted(Path(conda_root).glob('envs/*/bin'))
                     if item.is_dir())
    env['PATH'] = os.pathsep.join(paths + [env.get('PATH', '')])
    return env


def version_from_log(tool, log_text):
    if tool == 'novoplasty':
        match = re.search(r'(?im)^\s*Version\s+v?([0-9]+(?:\.[0-9]+){1,3})\b', log_text[:8192])
        if match:
            return match.group(1)
    match = re.search(r'(?i)(?:getorganelle|mitoz|mitoflex|mitofinder|novoplasty)[^\r\n]{0,80}?v?([0-9]+(?:\.[0-9]+){1,3})', log_text)
    return match.group(1) if match else None


def run_one(tool, exe, argv, cwd, extra_paths=(), python_site_paths=()):
    started = time.monotonic()
    log_path = cwd / 'tool.log'
    command = (list(exe) if isinstance(exe, (list, tuple)) else [exe]) + list(argv)
    try:
        env = with_tool_path(exe, extra_paths)
        if python_site_paths:
            shim = cwd / '.pythonpath_compat'
            shim.mkdir(exist_ok=True)
            (shim / 'sitecustomize.py').write_text(
                'import os, sys\n'
                'for p in os.environ.get("MITO_EXTRA_SITE_PACKAGES", "").split(os.pathsep):\n'
                '    if p and p not in sys.path:\n'
                '        sys.path.append(p)\n', encoding='utf-8')
            env['MITO_EXTRA_SITE_PACKAGES'] = os.pathsep.join(python_site_paths)
            env['PYTHONPATH'] = os.pathsep.join([str(shim.resolve()), env.get('PYTHONPATH', '')])
        with open(log_path, 'w', encoding='utf-8') as log:
            proc = subprocess.run(command, cwd=cwd, env=env,
                                  stdout=log, stderr=subprocess.STDOUT, check=False)
        log_text = log_path.read_text(encoding='utf-8', errors='replace')
        result = {'tool': tool, 'status': 'succeeded' if proc.returncode == 0 else 'failed',
                'exit_code': proc.returncode, 'command': command,
                'log': str(log_path), 'version_from_log': version_from_log(tool, log_text),
                'elapsed_seconds': round(time.monotonic() - started, 2)}
        if python_site_paths:
            result['supplemental_python_site_packages'] = list(python_site_paths)
        return result
    except OSError as exc:
        log_path.write_text('无法启动工具：%s\n' % exc, encoding='utf-8')
        return {'tool': tool, 'status': 'failed_to_start', 'exit_code': None,
                'command': command, 'log': str(log_path), 'error': str(exc),
                'elapsed_seconds': round(time.monotonic() - started, 2)}


def candidate_files(tool, workdir):
    files = sorted(p for p in workdir.rglob('*') if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES)
    if tool == 'getorganelle':
        files = [p for p in files if 'path_sequence' in p.name.lower()]
    elif tool == 'novoplasty':
        files = [p for p in files if 'contig' in p.name.lower() or 'circular' in p.name.lower()]
    elif tool == 'mitofinder':
        files = [p for p in files if 'contig' in p.name.lower() or 'mito' in p.name.lower()]
    elif tool == 'mitoflex':
        # MitoFlex's findmitoscaf module writes the selected mitochondrial
        # sequence set as <workname>.picked.fa; raw MEGAHIT contigs stay as
        # intermediate artifacts and must not be promoted to candidates.
        files = [p for p in files if p.name.endswith('.picked.fa')]
    elif tool == 'mitoz':
        # MitoZ's final selected assembly is emitted by circle_check as
        # <outprefix>.mitogenome.fa. Exclude raw contigs and annotation FASTAs.
        files = [p for p in files if p.name.endswith('.mitogenome.fa')]
    output = []
    seen_hashes = set()
    for path in files:
        records = parse_fasta(path)
        digest = sha256(path) if records and any(record['length'] > 0 for record in records) else None
        if digest and digest not in seen_hashes:
            seen_hashes.add(digest)
            output.append({'path': str(path.resolve()), 'records': records,
                           'sha256': digest, 'bytes': path.stat().st_size})
    return output


def novoplasty_config(args, outdir, read_length):
    if not args.novo_seed:
        raise ValueError('NOVOPlasty 需要 --novo-seed；不会自动从另一组装器的结果造种子')
    if not args.insert_size:
        raise ValueError('NOVOPlasty 需要 --insert-size（可由文库信息或 insert-size 估计工具获得）')
    if not args.genome_range:
        raise ValueError('NOVOPlasty 需要 --genome-range；不自动假定所有动物线粒体长度相同')
    seed = Path(args.novo_seed).resolve()
    if not seed.is_file():
        raise ValueError('NOVOPlasty 种子不存在: %s' % seed)
    seed_copy = outdir / 'input_seed.fasta'
    shutil.copyfile(seed, seed_copy)
    config = outdir / 'config.txt'
    values = [
        'Project:', '-----------------------', 'Project name          = %s' % args.sample,
        'Type                  = mito', 'Genome Range          = %s' % args.genome_range,
        'K-mer                 = 33', 'Max memory            =', 'Extended log          = 1',
        'Save assembled reads  = no', 'Seed Input            = %s' % seed_copy,
        'Extend seed directly  = no', 'Reference sequence    =', 'Variance detection    =',
        'Chloroplast sequence  =', '', 'Dataset 1:', '-----------------------',
        'Read Length           = %s' % read_length, 'Insert size            = %s' % args.insert_size,
        'Platform               = illumina', 'Single/Paired          = PE',
        'Combined reads         =', 'Forward reads          = %s' % Path(args.r1).resolve(),
        'Reverse reads          = %s' % Path(args.r2).resolve(), 'Store Hash             =', '',
        'Heteroplasmy:', '-----------------------', 'MAF                    =',
        'HP exclude list        =', 'PCR-free               =', '', 'Optional:',
        '-----------------------', 'Insert size auto       = no',
        'Use Quality Scores     = yes', 'Reduce ambigious N\'s  = no',
        # NOVOPlasty parses Output path as a directory and rejects the
        # otherwise-valid path unless it ends in a path separator.
        'Output path            = %s' % (str(outdir.resolve()) + os.sep),
    ]
    config.write_text('\n'.join(values) + '\n', encoding='utf-8')
    return config


def build_tools(args, outdir, read_length):
    candidates = []
    results = []
    requested = [item.strip().lower() for item in args.tools.split(',') if item.strip()]
    if len(set(requested)) != len(requested) or any(item not in TOOL_IDS for item in requested):
        raise ValueError('--tools 只能使用且不能重复: ' + ','.join(TOOL_IDS))
    if args.parallel_tools < 1:
        raise ValueError('--parallel-tools 必须为正整数')
    executor = None
    jobs_by_future = {}
    candidates_by_tool = {}
    for tool in requested:
        tool_root = outdir / tool
        if args.resume:
            attempts = [p for p in tool_root.glob('attempt_*') if p.is_dir()] if tool_root.is_dir() else []
            workdir = tool_root / ('attempt_%02d' % (len(attempts) + 1))
        else:
            workdir = tool_root
        workdir.mkdir(parents=True, exist_ok=True)
        launcher, extra_path, python_site_paths = None, [], []
        input_normalization = None
        if tool == 'getorganelle':
            exe = resolve('get_organelle_from_reads.py', args.getorganelle,
                          ('GETORGANELLE', 'GETORGANELLE_EXE'))
            # GetOrganelle refuses an existing -o path unless resuming or
            # overwriting. Keep logs in workdir but give it a fresh child path.
            organelle_outdir = workdir / 'getorganelle_output'
            argv = ['-1', str(Path(args.r1).resolve()), '-2', str(Path(args.r2).resolve()),
                    '-F', 'animal_mt', '-R', str(args.rounds), '-k', args.kmers,
                    '-t', str(args.threads), '-o', str(organelle_outdir)]
        elif tool == 'mitoflex':
            flex_root = os.environ.get('MITOFLEX_ROOT')
            root_exe = Path(flex_root) / 'MitoFlex.py' if flex_root else None
            root_exe = str(root_exe) if root_exe and root_exe.is_file() else None
            exe = resolve_script('MitoFlex.py', args.mitoflex or root_exe,
                                 ('MITOFLEX', 'MITOFLEX_SCRIPT'))
            if exe and Path(exe).suffix == '.py':
                flex_python, python_site_paths = mitoflex_python(args.mitoflex_python)
                if flex_python:
                    launcher = [flex_python, exe]
                    extra_path = [str(Path(flex_python).resolve().parent)]
            elif exe and os.access(exe, os.X_OK):
                launcher = exe
            flex_r1, flex_r2 = Path(args.r1), Path(args.r2)
            if launcher and flex_r1.suffix.lower() == '.gz' and flex_r2.suffix.lower() == '.gz':
                normalized_dir = workdir / 'normalized_reads'
                normalized_dir.mkdir(parents=True, exist_ok=True)
                normalized_r1 = normalized_dir / 'R1.single_member.fastq.gz'
                normalized_r2 = normalized_dir / 'R2.single_member.fastq.gz'
                input_normalization = [normalize_gzip_members(flex_r1, normalized_r1),
                                       normalize_gzip_members(flex_r2, normalized_r2)]
                flex_r1, flex_r2 = normalized_r1, normalized_r2
            argv = ['all', '--workname', args.sample, '--basedir', str(workdir.resolve()),
                    '--clade', args.mitoflex_clade, '--species-name', args.taxon,
                    '--fastq1', str(flex_r1.resolve()), '--fastq2', str(flex_r2.resolve()),
                    '--threads', str(args.threads), '--disable-visualization']
            if not args.mitoflex_use_taxonomy_filter:
                argv.append('--disable-taxa')
            if args.table:
                argv.extend(['--genetic-code', str(args.table)])
        elif tool == 'novoplasty':
            script = resolve_script('NOVOPlasty4.3.5.pl', args.novoplasty,
                                    ('NOVOPLASTY', 'NOVOPLASTY_SCRIPT'))
            perl = probe_command('perl', want_version=False)[0] or shutil.which('perl')
            exe = script
            launcher = [perl, script] if script and perl else None
            extra_path = [str(Path(perl).resolve().parent)] if perl else []
            needed = [flag for flag, value in (('--novo-seed', args.novo_seed),
                      ('--insert-size', args.insert_size), ('--genome-range', args.genome_range)) if not value]
            if script and not launcher:
                results.append({'tool': tool, 'status': 'unavailable', 'reason': '找到 NOVOPlasty 脚本但未找到 Perl'})
                workdir.rmdir()
                continue
            if script and needed:
                results.append({'tool': tool, 'status': 'not_run',
                                'reason': 'NOVOPlasty 需显式提供 ' + ', '.join(needed)})
                workdir.rmdir()
                continue
            config = novoplasty_config(args, workdir, read_length) if script else None
            argv = ['-c', str(config)] if config else []
        elif tool == 'mitoz':
            exe = resolve('mitoz', args.mitoz, ('MITOZ', 'MITOZ_EXE'))
            argv = ['all', '--outprefix', 'run', '--thread_number', str(args.threads),
                    '--workdir', str(workdir), '--clade', args.mitoz_clade,
                    '--genetic_code', str(args.table), '--species_name', args.taxon,
                    '--fq1', str(Path(args.r1).resolve()), '--fq2', str(Path(args.r2).resolve()),
                    '--requiring_taxa', args.mitoz_requiring_taxa or args.mitoz_clade,
                    '--requiring_relax', '6', '--min_abundance', '0']
        else:  # MitoFinder needs a reference GenBank and records its assembler.
            exe = resolve('mitofinder', args.mitofinder, ('MITOFINDER', 'MITOFINDER_EXE'))
            if not args.reference_genbank:
                results.append({'tool': tool, 'status': 'not_run',
                                'reason': 'MitoFinder reads route requires --reference-genbank'})
                workdir.rmdir()
                continue
            argv = ['--megahit', '-j', args.sample, '-1', str(Path(args.r1).resolve()),
                    '-2', str(Path(args.r2).resolve()), '-r', str(Path(args.reference_genbank).resolve()),
                    '-o', str(args.table or 5), '-p', str(args.threads), '-m', str(args.memory_gb)]
        if not exe:
            results.append({'tool': tool, 'status': 'unavailable',
                            'reason': '当前 PATH、CONDA_ROOT/envs/*/bin 和显式路径均未发现可执行程序'})
            continue
        if tool == 'mitoflex' and launcher is None:
            results.append({'tool': tool, 'status': 'not_run',
                            'reason': '找到 MitoFlex.py，但 CONDA_ROOT 中未找到可导入 numpy/pandas/ete3/Biopython/psutil 的 Python；可用 --mitoflex-python 显式指定'})
            workdir.rmdir()
            continue
        if tool == 'mitoz' and not args.table:
            results.append({'tool': tool, 'status': 'not_run',
                            'reason': 'MitoZ 注释需要显式 --table，不使用默认密码表'})
            continue
        if tool not in ('mitoflex', 'novoplasty'):
            launcher = exe
        if executor is None:
            executor = ThreadPoolExecutor(max_workers=min(args.parallel_tools, len(requested)))
        result_index = len(results)
        results.append(None)
        print('开始组装器 %s（工作目录：%s）' % (tool, workdir.resolve()), flush=True)
        future = executor.submit(run_one, tool, launcher, argv, workdir, extra_path, python_site_paths)
        jobs_by_future[future] = (result_index, tool, workdir, python_site_paths, input_normalization)
    for future in as_completed(jobs_by_future):
        result_index, tool, workdir, python_site_paths, input_normalization = jobs_by_future[future]
        result = future.result()
        result['workdir'] = str(workdir.resolve())
        if input_normalization:
            result['input_normalization'] = input_normalization
        outputs = candidate_files(tool, workdir)
        # Some wrappers (notably GetOrganelle) can print a fatal error and still
        # exit 0. An assembler is successful only if it produced a usable FASTA.
        if result.get('status') == 'succeeded' and not outputs:
            result['status'] = 'failed'
            result['postcondition_error'] = '工具退出码为 0，但未发现可用的候选 FASTA；请检查 tool.log'
        result['candidates'] = outputs
        candidates_by_tool[tool] = outputs
        results[result_index] = result
        print('组装器 %s 结束：状态=%s，候选 FASTA=%d，退出码=%s' %
              (tool, result['status'], len(outputs), result.get('exit_code')), flush=True)
    if executor is not None:
        executor.shutdown(wait=True)
    for tool in requested:
        candidates.extend(dict(candidate, producer=tool) for candidate in candidates_by_tool.get(tool, []))
    return results, candidates


def main(argv=None):
    parser = argparse.ArgumentParser(description='动物线粒体 Illumina 双端多组装候选编排（不自动挑选最终组装）')
    parser.add_argument('--r1', required=True, help='R1 FASTQ/FASTQ.GZ')
    parser.add_argument('--r2', required=True, help='R2 FASTQ/FASTQ.GZ')
    parser.add_argument('--outdir', required=True, help='新建样本交付目录；工具结果写入其中的 intermediate/assemblies')
    parser.add_argument('--sample', required=True)
    parser.add_argument('--taxon', required=True, help='用户确认的物种名，传递给支持该字段的工具')
    parser.add_argument('--tools', default='getorganelle,mitoflex,novoplasty,mitoz',
                        help='逗号分隔；可选 getorganelle,mitoflex,novoplasty,mitoz,mitofinder')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--parallel-tools', type=int, default=3,
                        help='同时运行的组装软件数量；每个软件仍各用 --threads 指定线程数')
    parser.add_argument('--rounds', type=int, default=10, help='GetOrganelle 招募轮数')
    parser.add_argument('--kmers', default='21,45,65,85,105', help='GetOrganelle k-mer 列表')
    parser.add_argument('--table', type=int, help='必须由用户/类群证据确认；MitoZ 路线必需')
    parser.add_argument('--mitoz-clade', default='Arthropoda')
    parser.add_argument('--mitoz-requiring-taxa', help='MitoZ 的目标类群过滤标签；默认等于 --mitoz-clade')
    parser.add_argument('--reference-genbank', help='MitoFinder 专用近缘参考 GenBank')
    parser.add_argument('--memory-gb', type=float, default=16)
    parser.add_argument('--novo-seed', help='NOVOPlasty 种子序列；不会自动从其它组装候选生成')
    parser.add_argument('--insert-size', type=int, help='NOVOPlasty 文库插入片段大小')
    parser.add_argument('--genome-range', help='NOVOPlasty 预期线粒体长度范围，如 10000-60000；需按类群提供')
    parser.add_argument('--getorganelle', help='GetOrganelle 可执行文件路径；优先于 PATH/conda 环境发现')
    parser.add_argument('--mitoflex', help='MitoFlex.py 可执行文件路径')
    parser.add_argument('--mitoflex-python', help='运行 MitoFlex.py 的兼容 Python 环境解释器')
    parser.add_argument('--mitoflex-clade', default='Arthropoda',
                        choices=('Platyhelminthes', 'Mollusca', 'Echinodermata', 'Nematoda',
                                 'Chordata', 'Annelida', 'Nemertea', 'Bryozoa', 'Porifera', 'Arthropoda'),
                        help='MitoFlex 内部检索用类群，需与其支持列表一致')
    parser.add_argument('--mitoflex-use-taxonomy-filter', action='store_true',
                        help='启用 MitoFlex 近缘类群过滤；默认跳过，避免数据库缺类群时误删真实候选')
    parser.add_argument('--novoplasty', help='NOVOPlasty Perl 脚本路径')
    parser.add_argument('--mitoz', help='MitoZ 可执行文件路径')
    parser.add_argument('--mitofinder', help='MitoFinder 可执行文件路径')
    parser.add_argument('--resume', action='store_true',
                        help='在相同输入与既有 assembly_candidates.json 上追加重跑失败/未执行的工具；不重跑已成功工具')
    args = parser.parse_args(argv)

    try:
        r1, r2 = Path(args.r1).resolve(), Path(args.r2).resolve()
        if not r1.is_file() or not r2.is_file():
            raise ValueError('R1/R2 必须是现存文件')
        if r1 == r2:
            raise ValueError('R1 与 R2 不能是同一个文件')
        if args.threads < 1 or args.rounds < 1 or args.memory_gb <= 0:
            raise ValueError('threads、rounds、memory-gb 必须大于 0')
        args.r1, args.r2 = str(r1), str(r2)
        outdir = Path(args.outdir).resolve()
        previous_manifest_path = outdir / 'intermediate' / 'assemblies' / 'assembly_candidates.json'
        prior = None
        if args.resume:
            if not outdir.is_dir() or not previous_manifest_path.is_file():
                raise ValueError('--resume 要求已有样本目录和 assembly_candidates.json')
            try:
                prior = json.loads(previous_manifest_path.read_text(encoding='utf-8'))
            except (OSError, ValueError) as exc:
                raise ValueError('旧 assembly manifest 无法读取: %s' % exc) from exc
            if prior.get('format') != 'mito-illumina-assembly-candidates-1':
                raise ValueError('--resume 不支持此 assembly manifest 格式')
            if prior.get('sample') != args.sample or prior.get('taxon') != args.taxon:
                raise ValueError('--resume 的样本名/物种与既有 assembly manifest 不一致')
            if any(item.get('role') == 'R1' and item.get('path') != str(r1) for item in prior.get('inputs', [])) or \
               any(item.get('role') == 'R2' and item.get('path') != str(r2) for item in prior.get('inputs', [])):
                raise ValueError('--resume 的 R1/R2 路径与既有 assembly manifest 不一致')
        else:
            precreated_for_background = outdir.exists() and has_only_background_tasks(outdir)
            if outdir.exists() and not precreated_for_background:
                raise ValueError('为保护既有结果，输出目录必须不存在: %s' % outdir)
        if outdir == r1 or outdir == r2:
            raise ValueError('输出目录不能覆盖输入 reads')
        if prior:
            print('续跑前校验 R1/R2 文件哈希…', flush=True)
            input_records = [{'role': role, 'path': str(path), 'bytes': path.stat().st_size,
                              'sha256': sha256(path)} for role, path in (('R1', r1), ('R2', r2))]
            if prior.get('inputs') != input_records:
                raise ValueError('--resume 的 FASTQ 内容或大小与既有 manifest 不一致')
            fastq_summary = {'read_count': prior.get('read_count'),
                             'read_length_first_record': prior.get('read_length_first_record')}
        else:
            print('完整校验配对 FASTQ 结构、read ID、质量长度和压缩流…', flush=True)
            fastq_summary = validate_fastq_pair(r1, r2)
            input_records = [{'role': role, 'path': str(path), 'bytes': path.stat().st_size,
                              'sha256': sha256(path)} for role, path in (('R1', r1), ('R2', r2))]
        print('FASTQ 校验完成：%s reads；开始调度组装器。' % fastq_summary['read_count'], flush=True)
        if args.table is not None and args.table <= 0:
            raise ValueError('--table 必须为正整数')
        if not 0 < args.memory_gb <= 4096:
            raise ValueError('--memory-gb 必须在 (0,4096] 范围内')
        if args.genome_range and not __import__('re').fullmatch(r'\d+-\d+', args.genome_range):
            raise ValueError('--genome-range 格式应为正整数范围，如 10000-60000')
        if not args.resume and not precreated_for_background:
            outdir.mkdir(parents=True)
        intermediates = outdir / 'intermediate' / 'assemblies'
        intermediates.mkdir(parents=True, exist_ok=True)
        requested = [item.strip().lower() for item in args.tools.split(',') if item.strip()]
        previous_tools = prior.get('tools', []) if prior else []
        previous_success = {
            item.get('tool') for item in previous_tools
            if item.get('status') == 'succeeded' and item.get('candidates')
        }
        if args.resume and previous_success.intersection(requested):
            raise ValueError('--resume 不会覆盖已成功工具: ' + ','.join(sorted(previous_success.intersection(requested))) +
                             '；请只指定失败或未执行的工具')
        previous_candidates = prior.get('candidates', []) if prior else []
        for entry in previous_candidates:
            candidate_path = Path(entry.get('path', '')).resolve()
            if not candidate_path.is_file() or sha256(candidate_path) != entry.get('sha256'):
                raise ValueError('旧候选缺失或 SHA256 漂移，拒绝 resume: %s' % candidate_path)
        tool_results, candidates = build_tools(args, intermediates, fastq_summary['read_length_first_record'])
        all_tool_results = previous_tools + tool_results
        all_candidates = []
        seen_candidate_hashes = set()
        for candidate in previous_candidates + candidates:
            key = (candidate.get('producer'), candidate.get('sha256'))
            if not key[0] or not key[1] or key in seen_candidate_hashes:
                continue
            seen_candidate_hashes.add(key)
            all_candidates.append(candidate)
        successes = sum(item.get('status') == 'succeeded' for item in all_tool_results)
        independent = len({candidate['producer'] for candidate in all_candidates})
        manifest = {
            'format': 'mito-illumina-assembly-candidates-1',
            'created_at': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
            'sample': args.sample, 'taxon': args.taxon,
            'platform': 'Illumina paired-end short reads',
            'inputs': input_records,
            'read_count': fastq_summary['read_count'],
            'read_length_first_record': fastq_summary['read_length_first_record'],
            'tools': all_tool_results, 'candidate_count': len(all_candidates), 'candidates': all_candidates,
            'resume_count': (prior.get('resume_count', 0) + 1) if prior else 0,
            'assembly_decision': 'not_selected',
            'interpretation': ('多工具候选已生成，须进入 reads 回贴、组装图/接缝和候选差异裁决。'
                               if independent >= 2 else
                               '独立组装器不足两个或没有候选；不能把本次结果称为多工具比较完成。'),
            'limitations': [
                '成功退出不证明组装正确、唯一、完整或闭环。',
                '候选不会按软件投票；序列、图结构、read 支持和竞争解释须联合判断。',
                'NOVOPlasty 仅在用户提供种子、插入长度和类群适用的长度范围时运行。',
                'MitoFinder 仅在提供近缘参考 GenBank 时运行；参考不用于强制改写样本碱基。',
                '原始 reads 未修改；本步骤不选择或覆盖最终交付 FASTA。'],
        }
        manifest_path = intermediates / 'assembly_candidates.json'
        temp_manifest = manifest_path.with_suffix('.json.tmp')
        temp_manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        temp_manifest.replace(manifest_path)
        print('组装候选清单：%s' % manifest_path)
        print('成功运行记录：%d；候选 FASTA：%d；独立候选来源：%d' %
              (successes, len(all_candidates), independent))
        if independent < 2:
            print('状态：证据不足以完成多组装比较；详见各工具状态和日志。')
            return 2
        return 0
    except (OSError, ValueError) as exc:
        print('错误：%s' % exc, file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())

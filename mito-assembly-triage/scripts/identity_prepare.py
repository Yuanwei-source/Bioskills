#!/usr/bin/env python3
"""Prepare sample-derived identity markers without uploading sample data.

Exit 0: usable marker queries; 1: no reliable query (not a non-insect verdict);
2: input/configuration refusal; 3: tool/network/format failure.
"""
import argparse
import difflib
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import urllib.request

from _identity_markers import (check_annotation, fasta_records, fastq_records,
                               obtain_protein_references, paired_records, run_tool, sample_pairs,
                               save_json, supported_segments, write_fasta)
from run_illumina_candidates import resolve, sha256

ROOT = Path(__file__).resolve().parents[1]


def discover_mge(args):
    if args.mge:
        path = Path(args.mge).resolve()
        if not path.is_file():
            raise ValueError('MitoGeneExtractor executable does not exist')
        return str(path)
    for name in ('MitoGeneExtractor', 'MitoGeneExtractor-v1.9.6', 'mitogeneextractor'):
        path = resolve(name, env_names=('MITOGENEEXTRACTOR',))
        if path:
            return path
    raise ValueError('Install MitoGeneExtractor >=1.9.6 or provide --mge')


def recover_ssu(args, outdir, reads, tools):
    """Optional Rfam model recruitment + SPAdes; experimental short-read route."""
    cmsearch = resolve('cmsearch', args.cmsearch, ('CMSEARCH',))
    spades = resolve('spades.py', args.spades, ('SPADES',))
    if not cmsearch or not spades:
        raise ValueError('--recover-18s requires cmsearch and spades.py')
    directory = outdir / 'ssu'; directory.mkdir()
    model = directory / 'RF01960.cm'
    if args.ssu_model:
        model.write_bytes(Path(args.ssu_model).read_bytes())
        source = str(Path(args.ssu_model).resolve())
    else:
        source = 'https://rfam.org/family/RF01960/cm'
        with urllib.request.urlopen(source, timeout=120) as response:
            model.write_bytes(response.read())
    if 'RF01960' not in model.read_text() or not model.read_text().startswith('INFERNAL'):
        raise ValueError('Expected Rfam RF01960 Infernal covariance model')
    nucleotide = directory / 'sample.fasta'
    write_fasta(nucleotide, [(name, seq) for path in reads for name, seq, _ in fastq_records(path)])
    table = directory / 'recruitment.tbl'
    tools.append(run_tool([cmsearch, '--anytrunc', '--cpu', str(args.threads), '--noali',
                           '-E', '0.0001', '--tblout', str(table), str(model), str(nucleotide)],
                          directory, 'cmsearch_reads'))
    accepted = set()
    for line in table.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            accepted.add(line.split()[0].rsplit('/', 1)[0])
    if not accepted:
        return [], {'status': 'no_recruited_ssu', 'source': source, 'sha256': sha256(model)}
    selected = []
    for mate, path in enumerate(reads, 1):
        target = directory / ('ssu_R%d.fastq' % mate)
        with target.open('w') as handle:
            for name, seq, qual in fastq_records(path):
                if name.rsplit('/', 1)[0] in accepted:
                    handle.write('@%s\n%s\n+\n%s\n' % (name, seq, qual))
        selected.append(target)
    assembly = directory / 'assembly'
    tools.append(run_tool([spades, '--only-assembler', '-1', str(selected[0]), '-2', str(selected[1]),
                           '-o', str(assembly), '-t', str(args.threads), '-m', '4'], directory, 'spades_ssu'))
    contigs = assembly / 'contigs.fasta'
    table = directory / 'contigs.tbl'
    tools.append(run_tool([cmsearch, '--anytrunc', '--cpu', str(args.threads), '--noali',
                           '-E', '0.0001', '--tblout', str(table), str(model), str(contigs)],
                          directory, 'cmsearch_contigs'))
    sequences = dict(fasta_records(contigs))
    queries = []
    from Bio.Seq import Seq
    for line in table.read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        row = line.split(); start, end = int(row[7]), int(row[8])
        seq = sequences[row[0]][min(start, end)-1:max(start, end)]
        if row[9] == '-':
            seq = str(Seq(seq).reverse_complement())
        if len(seq) >= args.min_length and set(seq) <= set('ACGT'):
            queries.append({'sequence': seq, 'marker': '18S', 'origin': row[0],
                            'evidence': 'experimental_rfam_recruitment_spades',
                            'identity_verified': False,
                            'limitations': ['Not independently validated for mixed samples',
                                            'Contig support requires follow-up; not a host certificate']})
    return queries, {'status': 'experimental', 'model_source': source, 'sha256': sha256(model),
                     'recruited_templates': len(accepted)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument('--fasta', help='One or more assembled FASTA records')
    sources.add_argument('--r1', action='append', help='R1 FASTQ; repeat for each lane')
    parser.add_argument('--r2', action='append', help='Matching R2; repeat in same order')
    parser.add_argument('--annotation', help='Optional sequence-matched GenBank for independent comparison')
    parser.add_argument('--outdir', required=True, help='New intermediate/identity_screen directory')
    parser.add_argument('--expected-taxon', default='Insecta', help='User expectation, not identification')
    parser.add_argument('--sample-pairs', type=int, default=1000000)
    parser.add_argument('--seed', type=int, default=19)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--min-length', type=int, default=300)
    parser.add_argument('--min-templates', type=int, default=5)
    parser.add_argument('--tables', default='5,2', help='Discovery codes; never auto-assign annotation code')
    parser.add_argument('--protein-ref', help='Small broad-class protein FASTA; default pinned discovery panel')
    parser.add_argument('--mge'); parser.add_argument('--exonerate'); parser.add_argument('--fastp')
    parser.add_argument('--recover-18s', action='store_true', help='Experimental optional nuclear marker route')
    parser.add_argument('--ssu-model'); parser.add_argument('--cmsearch'); parser.add_argument('--spades')
    args = parser.parse_args(argv)
    outdir = Path(args.outdir).resolve()
    created = False
    manifest = {'format': 'mito-identity-markers-1', 'expected_taxon': args.expected_taxon,
                'queries': [], 'inputs': [], 'tools': [], 'observations': [], 'status': 'preparing'}
    try:
        tables = list(dict.fromkeys(int(value) for value in args.tables.split(',')))
        from Bio.Data import CodonTable
        if (not tables or any(t not in CodonTable.unambiguous_dna_by_id for t in tables)
                or args.sample_pairs < 1 or args.threads < 1 or args.min_length < 100
                or args.min_templates < 2 or not args.expected_taxon.strip()):
            raise ValueError('Invalid table, sampling, marker length or expected taxon')
        manifest['discovery_parameters'] = {'tables': tables, 'min_length': args.min_length,
                                            'min_templates': args.min_templates, 'min_fraction': .9}
        if args.r1 and (not args.r2 or len(args.r1) != len(args.r2)):
            raise ValueError('Each --r1 must have one --r2')
        if args.fasta and args.r2 or args.r1 and args.annotation:
            raise ValueError('Annotation needs --fasta; --r2 needs --r1')
        if args.recover_18s and not args.r1:
            raise ValueError('--recover-18s currently supports paired reads only')
        if outdir.exists():
            raise ValueError('Output directory already exists; use a new directory')
        inputs = ([args.fasta] if args.fasta else args.r1 + args.r2)
        if args.r1:
            for index, left in enumerate(inputs):
                for right in inputs[index + 1:]:
                    if os.path.samefile(left, right):
                        raise ValueError('FASTQ files must be distinct, including across lanes')
        if args.annotation:
            inputs += [args.annotation]
        if args.protein_ref:
            inputs += [args.protein_ref]
        if args.ssu_model:
            inputs += [args.ssu_model]
        for name in inputs:
            path = Path(name).resolve()
            if not path.is_file():
                raise ValueError('Missing input: %s' % path)
            if outdir == path or path in outdir.parents:
                raise ValueError('Output overlaps input')
            manifest['inputs'].append({'path': str(path), 'sha256': sha256(path)})
        if args.r1:
            raw_hashes = [entry['sha256'] for entry in manifest['inputs'][:len(args.r1) + len(args.r2)]]
            if len(set(raw_hashes)) != len(raw_hashes):
                raise ValueError('FASTQ inputs include duplicate file content; copied lanes are not independent')
        mge = discover_mge(args)
        version_run = subprocess.run([mge, '--version'], capture_output=True, text=True, check=True)
        version_output = version_run.stdout + version_run.stderr
        version = re.search(r'(\d+)\.(\d+)\.(\d+)', version_output)
        if not version or tuple(map(int, version.groups())) < (1, 9, 6) or 'beta' in version_output.lower():
            raise ValueError('MitoGeneExtractor >=1.9.6 stable required')
        help_run = subprocess.run([mge, '--help'], capture_output=True, text=True, check=True)
        mge_help = help_run.stdout + help_run.stderr
        if '--vulgar_directory' not in mge_help and '--vulgar_file' not in mge_help:
            raise ValueError('Unsupported MitoGeneExtractor CLI; inspect --help')
        if '--treat-references-as-individual' not in mge_help:
            raise ValueError('MGE must support independent reference discovery')
        manifest['mge_cli'] = {'path': mge, 'help_sha256': hashlib.sha256(mge_help.encode()).hexdigest(),
                               'version': version_output.strip(),
                               'vulgar_directory': '--vulgar_directory' in mge_help}
        exonerate = resolve('exonerate', args.exonerate, ('EXONERATE',))
        if not exonerate:
            raise ValueError('Missing exonerate')
        fastp = resolve('fastp', args.fastp, ('FASTP',)) if args.r1 else None
        if args.r1 and not fastp:
            raise ValueError('FASTQ identity preparation requires fastp')
        genomes = fasta_records(args.fasta, 'ACGTRYSWKMBDHVN') if args.fasta else None
        annotated = check_annotation(args.annotation, genomes) if args.annotation else []
        outdir.mkdir(parents=True); created = True
        reference_dir = outdir / 'references'; reference_dir.mkdir()
        if args.protein_ref:
            protein = reference_dir / 'custom.faa'
            protein.write_bytes(Path(args.protein_ref).read_bytes())
            fasta_records(protein, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*')
            manifest['references'] = [{'path': str(protein), 'sha256': sha256(protein),
                                       'source': 'user supplied discovery reference; coverage not verified'}]
        else:
            protein, manifest['references'] = obtain_protein_references(
                reference_dir, ROOT / 'config/identity_reference_sources.json')
        datasets = []
        if args.r1:
            sample_dir = outdir / 'sample'; sample_dir.mkdir()
            raw, manifest['sampling'] = sample_pairs(
                [(Path(l).resolve(), Path(r).resolve()) for l, r in zip(args.r1, args.r2)],
                args.sample_pairs, args.seed, sample_dir)
            filtered = [sample_dir / ('filtered_R%d.fastq' % i) for i in (1, 2)]
            manifest['tools'].append(run_tool([fastp, '-i', str(raw[0]), '-I', str(raw[1]),
                '-o', str(filtered[0]), '-O', str(filtered[1]), '-q', '20', '-u', '0',
                '-n', '0', '--thread', str(args.threads), '--json', str(sample_dir / 'fastp.json'),
                '--html', str(sample_dir / 'fastp.html')], sample_dir, 'fastp'))
            qc = json.loads((sample_dir / 'fastp.json').read_text())
            if qc['summary']['before_filtering']['total_reads'] != 2 * manifest['sampling']['sampled_pairs']:
                raise ValueError('fastp did not process all sampled mates')
            manifest['sampling']['filtered_reads'] = qc['summary']['after_filtering']['total_reads']
            if 2 * sum(1 for _ in paired_records(*filtered)) != manifest['sampling']['filtered_reads']:
                raise ValueError('Filtered mate counts disagree with fastp report')
            datasets = [('sampled_reads', ['-q', str(filtered[0]), '-q', str(filtered[1])], 'reads')] if manifest['sampling']['filtered_reads'] else []
        else:
            for index, (name, seq) in enumerate(genomes):
                path = outdir / ('source_%d.fasta' % index)
                write_fasta(path, [('source_%d' % index, seq)])
                datasets.append((name, ['-d', str(path)], 'fasta'))
        for index, (origin, source_args, mode) in enumerate(datasets):
            for table in tables:
                directory = outdir / ('mge_%d_table%d' % (index, table)); directory.mkdir()
                vulgar = directory if '--vulgar_directory' in mge_help else directory / 'hits.vulgar'
                manifest['tools'].append(run_tool([mge, *source_args, '-p', str(protein),
                    '--treat-references-as-individual',
                    '-e', exonerate, '-o', str(directory / 'alignment_'),
                    '-c', str(directory / 'raw_consensus_'), '-V', str(vulgar),
                    '-C', str(table), '-n', '0', '--temporaryDirectory', str(directory)], directory, 'mge'))
                for alignment in sorted(directory.glob('alignment_*.fas')):
                    if not alignment.read_text().strip():
                        manifest['observations'].append({'origin': origin, 'table': table,
                            'alignment': str(alignment), 'status': 'no_alignment', 'usable_segments': 0})
                        continue
                    segments, qc = supported_segments(alignment, mode, args.min_length, args.min_templates)
                    save_json(alignment.with_suffix('.qc.json'), qc)
                    manifest['observations'].append({'origin': origin, 'table': table,
                        'alignment': str(alignment), 'mixed_columns': qc['mixed_columns'],
                        'usable_segments': len(segments)})
                    for sequence in segments:
                        matches = []
                        for annotation in annotated:
                            if annotation['record'] == origin:
                                matched_bp = difflib.SequenceMatcher(None, sequence, annotation['sequence'],
                                                                    autojunk=False).find_longest_match(
                                                                    0, len(sequence), 0, len(annotation['sequence'])).size
                                if matched_bp >= args.min_length:
                                    matches.append(dict(annotation, independent_overlap_bp=matched_bp,
                                                        precise_boundary_verified=False))
                        manifest['queries'].append({'sequence': sequence, 'marker': 'COX1',
                            'origin': origin, 'table_used_for_discovery': table,
                            'evidence': 'template_supported' if mode == 'reads' else 'protein_homology',
                            'annotation_matches': [{k: v for k, v in a.items() if k != 'sequence'} for a in matches],
                            'identity_verified': False, 'qc': str(alignment.with_suffix('.qc.json'))})
        if args.recover_18s:
            queries, manifest['ssu'] = recover_ssu(args, outdir, filtered, manifest['tools'])
            manifest['queries'].extend(queries)
        # Deduplicate representations, retaining every source/branch provenance.
        from Bio.Seq import Seq
        unique = {}
        for query in manifest['queries']:
            seq = query.pop('sequence')
            key = min(seq, str(Seq(seq).reverse_complement()))
            unique.setdefault((query['marker'], key), {'sequence': seq, 'provenance': []})['provenance'].append(query)
        manifest['queries'] = []
        for index, ((marker, _), value) in enumerate(unique.items(), 1):
            path = outdir / ('query_%03d.fasta' % index)
            write_fasta(path, [('query_%03d' % index, value['sequence'])])
            manifest['queries'].append({'id': 'query_%03d' % index, 'marker': marker,
                'path': str(path), 'sha256': sha256(path), 'length': len(value['sequence']),
                'provenance': value['provenance']})
        manifest['mode'] = 'reads' if args.r1 else ('assembly_annotation' if args.annotation else 'assembly')
        manifest['scope'] = ('Marker evidence in sampled reads; not whole-library composition' if args.r1
                             else 'Assembly candidate origin only; not original library or host verification')
        manifest['status'] = 'prepared' if manifest['queries'] else 'insufficient'
        for entry in manifest['inputs']:
            if sha256(entry['path']) != entry['sha256']:
                raise ValueError('Source input changed during marker preparation')
        save_json(outdir / 'markers.json', manifest)
        print(json.dumps({'manifest': str(outdir / 'markers.json'), 'status': manifest['status'],
                          'queries': len(manifest['queries'])}, ensure_ascii=False))
        return 0 if manifest['queries'] else 1
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, IndexError, ImportError,
            subprocess.CalledProcessError) as exc:
        manifest.update(status='preparation_failed', error=str(exc))
        if created:
            save_json(outdir / 'markers.json', manifest)
        print(str(exc), file=sys.stderr)
        return 3 if created else 2


if __name__ == '__main__':
    raise SystemExit(main())

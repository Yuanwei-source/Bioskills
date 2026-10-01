"""Local marker discovery; discovery references never establish sample identity."""
import collections
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import random
import re
import subprocess
import urllib.request

from _evidence_io import atomic_write_text
from run_illumina_candidates import normalized_read_id, sha256, with_tool_path


def save_json(path, value):
    atomic_write_text(path, json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def fasta_records(path, alphabet='ACGTRYSWKMBDHVN-~', unique=True):
    records, identifier, pieces = [], None, []
    for line in Path(path).read_text(encoding='ascii').splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if not line[1:].split():
                raise ValueError('Missing FASTA identifier: %s' % path)
            if identifier is not None:
                records.append((identifier, ''.join(pieces).upper()))
            identifier, pieces = line[1:].split()[0], []
        else:
            if identifier is None or any(c.upper() not in alphabet for c in line):
                raise ValueError('Invalid FASTA sequence: %s' % path)
            pieces.append(line)
    if identifier is not None:
        records.append((identifier, ''.join(pieces).upper()))
    if not records or any(not seq for _, seq in records):
        raise ValueError('Empty FASTA record: %s' % path)
    if unique and len({name for name, _ in records}) != len(records):
        raise ValueError('Duplicate FASTA identifiers: %s' % path)
    return records


def write_fasta(path, records):
    atomic_write_text(path, ''.join('>%s\n%s\n' % row for row in records))


def fastq_records(path):
    opener = gzip.open if str(path).endswith(('.gz', '.bgz')) else open
    with opener(path, 'rt', encoding='ascii') as stream:
        while True:
            header = stream.readline()
            if not header:
                return
            seq, plus, qual = (stream.readline().rstrip('\r\n') for _ in range(3))
            if (not header.startswith('@') or not header[1:].split() or not plus.startswith('+') or not seq
                    or len(seq) != len(qual) or not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+', seq)
                    or not re.fullmatch('[!-~]+', qual)):
                raise ValueError('Invalid or truncated FASTQ: %s' % path)
            yield header[1:].split()[0], seq.upper(), qual


def paired_records(r1, r2):
    for left, right in itertools.zip_longest(fastq_records(r1), fastq_records(r2)):
        if left is None or right is None or normalized_read_id(left[0]) != normalized_read_id(right[0]):
            raise ValueError('FASTQ mate counts/identifiers do not match: %s, %s' % (r1, r2))
        yield left, right


def sample_pairs(lanes, budget, seed, directory):
    """Two passes bound memory; validate every member, then sample each lane."""
    if budget < len(lanes):
        raise ValueError('Sample budget must allow at least one pair per lane')
    counts = [sum(1 for _ in paired_records(*lane)) for lane in lanes]
    print('FASTQ validation complete: lane pairs=%s; sampling budget=%s' % (counts, budget), flush=True)
    if any(count == 0 for count in counts):
        raise ValueError('Empty FASTQ lane')
    allocations = [min(count, budget // len(lanes)) for count in counts]
    remaining = min(budget, sum(counts)) - sum(allocations)
    for i, count in enumerate(counts):
        extra = min(remaining, count - allocations[i])
        allocations[i] += extra
        remaining -= extra
    paths = [directory / 'sample_R1.fastq', directory / 'sample_R2.fastq']
    rng = random.Random(seed)
    with paths[0].open('w') as one, paths[1].open('w') as two:
        for lane_idx, (lane, count, amount) in enumerate(zip(lanes, counts, allocations)):
            chosen = set(rng.sample(range(count), amount))
            for index, pair in enumerate(paired_records(*lane)):
                if index in chosen:
                    for mate, (record, handle) in enumerate(zip(pair, (one, two)), 1):
                        # Bind to original QNAME, not row number: duplicate SAM/
                        # FASTQ records must not become independent templates.
                        template = int(hashlib.sha256(normalized_read_id(record[0]).encode()).hexdigest(), 16)
                        name = 'L%d_P%d/%d' % (lane_idx, template, mate)
                        handle.write('@%s\n%s\n+\n%s\n' % (name, record[1], record[2]))
            print('Sampled lane %d: %d paired reads' % (lane_idx + 1, amount), flush=True)
    return paths, {'seed': seed, 'budget_pairs': budget, 'lane_counts': counts,
                   'sampled_per_lane': allocations, 'sampled_pairs': sum(allocations)}


def run_tool(command, directory, name):
    log = directory / (name + '.log')
    print('Running %s; log=%s' % (name, log), flush=True)
    with log.open('w') as handle:
        result = subprocess.run(command, cwd=directory, env=with_tool_path(command),
                                stdin=subprocess.DEVNULL, stdout=handle, stderr=subprocess.STDOUT)
    evidence = {'command': command, 'exit_code': result.returncode, 'log': str(log)}
    if result.returncode:
        raise RuntimeError('%s failed with exit %s; see %s' % (name, result.returncode, log))
    return evidence


def obtain_protein_references(directory, config):
    """Pinned tiny discovery panel; not an identification database."""
    data = json.loads(Path(config).read_text())
    records, sources = [], []
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*'
    for reference in data['references']:
        with urllib.request.urlopen(reference['url'], timeout=60) as response:
            payload = response.read()
        path = directory / (reference['id'] + '.faa')
        path.write_bytes(payload)
        if sha256(path) != reference['sha256']:
            raise ValueError('Discovery reference checksum mismatch: %s' % reference['id'])
        entries = fasta_records(path, alphabet)
        if len(entries) != 1:
            raise ValueError('Expected one discovery protein per reference')
        records.append((reference['id'], entries[0][1]))
        sources.append(dict(reference, path=str(path)))
    target = directory / 'discovery_proteins.faa'
    write_fasta(target, records)
    return target, sources


def supported_segments(alignment, mode, min_length=300, min_templates=5, fraction=.9):
    """Call only template-supported columns; refuse conflicting mixture consensus.

    MGE pads extracted regions against a protein reference; never concatenate
    across unsupported columns. Different records in assemblies are processed
    separately upstream. Mate conflicts contribute neither depth nor a base.
    """
    rows = fasta_records(alignment, unique=False)
    lengths = {len(seq) for _, seq in rows}
    if len(lengths) != 1:
        raise ValueError('MGE alignment rows have unequal lengths')
    if mode == 'fasta':
        if len(rows) != 1:
            raise ValueError('Assembly extraction must contain one source sequence')
        # Keep coordinate gaps as unknown so disjoint HSPs are not glued together.
        sequence = re.sub('[^ACGT]', 'N', rows[0][1])
        return [m.group() for m in re.finditer('[ACGT]{%d,}' % min_length, sequence)], {
            'mode': mode, 'mixed_columns': 0, 'quality_limit': 'No reads support assessed'}
    groups = collections.defaultdict(list)
    for name, seq in rows:
        # Prepared read identifiers are controlled; MGE can append alignment info.
        match = re.search(r'L\d+_P\d+', name)
        if not match:
            raise ValueError('Cannot bind MGE alignment row to a sampled template: %s' % name)
        groups[match.group()].append(seq)
    sequence, depths, conflicts = [], [], []
    for index in range(next(iter(lengths))):
        counts = collections.Counter()
        for reads in groups.values():
            bases = {read[index] for read in reads if read[index] in 'ACGT'}
            if len(bases) == 1:
                counts[next(iter(bases))] += 1
        depth = sum(counts.values())
        depths.append(depth)
        ordered = counts.most_common()
        if len(ordered) > 1 and ordered[1][1] >= 3 and ordered[1][1] / depth >= .2:
            conflicts.append(index + 1)
        sequence.append(ordered[0][0] if depth >= min_templates
                        and ordered[0][1] / depth >= fraction else 'N')
    qc = {'mode': mode, 'templates': len(groups), 'column_depth': depths,
          'min_templates': min_templates, 'min_fraction': fraction,
          'mixed_columns': len(conflicts), 'mixed_column_positions': conflicts,
          'limit': 'Template support is not evidence excluding NUMTs or mixed haplotypes'}
    if conflicts:
        return [], qc
    return [m.group() for m in re.finditer('[ACGT]{%d,}' % min_length, ''.join(sequence))], qc


def check_annotation(path, genomes):
    from Bio import SeqIO
    records = list(SeqIO.parse(str(path), 'genbank'))
    if not records or len({r.id for r in records}) != len(records):
        raise ValueError('Empty or duplicate annotation records')
    genome_map = dict(genomes)
    if set(genome_map) != {r.id for r in records}:
        raise ValueError('Annotation and FASTA record IDs must match')
    candidates = []
    for record in records:
        if str(record.seq).upper() != genome_map[record.id]:
            raise ValueError('Annotation sequence does not match FASTA: %s' % record.id)
        for feature in record.features:
            gene = ' '.join(feature.qualifiers.get('gene', [])).lower()
            product = ' '.join(feature.qualifiers.get('product', [])).lower()
            is_cox1 = gene in ('cox1', 'co1', 'coi', 'cox i') or (
                'cytochrome' in product and re.search(r'(?:subunit|oxidase)\s+(?:i|1)\b', product))
            if feature.type == 'CDS' and is_cox1:
                candidates.append({'record': record.id, 'sequence': str(feature.extract(record.seq)).upper(),
                                   'location': str(feature.location), 'gene': gene})
    return candidates

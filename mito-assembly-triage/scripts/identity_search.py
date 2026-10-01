#!/usr/bin/env python3
"""Resumable NCBI marker BLAST and conservative taxonomic consistency screen.

Upload needs --allow-public-upload. Calls are rate-limited and durable RID state
is saved before polling. No reference-name or first-hit species assignment.
Exit 0: screening verdict; 1: insufficient; 2: refusal; 3: service/format failure.
"""
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import cox1_id
from _identity_markers import fasta_records, save_json
from run_illumina_candidates import sha256

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
ALIASES = {'昆虫': 'Insecta', 'insect': 'Insecta', 'insects': 'Insecta', '节肢动物': 'Arthropoda'}


def local(tag):
    return tag.rsplit('}', 1)[-1]


def first(node, name):
    return next((child for child in node.iter() if local(child.tag) == name), None)


def content(node, name):
    item = first(node, name)
    return item.text.strip() if item is not None and item.text else ''


def parse_batch_xml(payload, queries):
    root = ET.fromstring(payload)
    searches = [node for node in root.iter() if local(node.tag) == 'Search']
    if not searches:
        raise ValueError('Expected XML2 Search records')
    results = {}
    for search in searches:
        name = content(search, 'query-title').split()[0] if content(search, 'query-title') else ''
        if name not in queries or name in results:
            raise ValueError('BLAST query ID missing/duplicate/unexpected: %s' % name)
        length = int(content(search, 'query-len'))
        if length != queries[name]['length']:
            raise ValueError('BLAST query length differs from uploaded sequence')
        wrapper = ET.Element('BlastXML2'); wrapper.append(search)
        parsed = cox1_id.parse_blast_xml(ET.tostring(wrapper, encoding='unicode'))
        xml_hits = [node for node in search.iter() if local(node.tag) == 'Hit']
        # Existing parser sorts hits, so match metadata by its original accession.
        metadata = {}
        for hit in xml_hits:
            descriptions = [node for node in hit.iter() if local(node.tag) == 'HitDescr']
            if descriptions:
                key = content(descriptions[0], 'id') or content(descriptions[0], 'accession')
                records = [{'accession': content(d, 'accession') or content(d, 'id'),
                            'taxid': int(content(d, 'taxid')) if content(d, 'taxid').isdigit() else None,
                            'title': content(d, 'title')} for d in descriptions]
                metadata[key] = {'accession': records[0]['accession'], 'descriptions': records,
                                 'taxids': list(dict.fromkeys(record['taxid'] for record in records))}
        for hit in parsed['hits']:
            info = metadata.get(hit['accession'], {})
            hit.update(info)
            hit.setdefault('taxids', [None])
            hit['ranking_bitscore'] = sum(hsp['bitscore'] for hsp in hit['hsps'])
        results[name] = {'query_length': length, 'hits': parsed['hits'],
                         'database_actual': content(root, 'db') or None}
    if set(results) != set(queries):
        raise ValueError('BLAST returned only part of the uploaded marker batch')
    return results


def parse_taxonomy(payload):
    root = ET.fromstring(payload)
    if root.tag != 'TaxaSet':
        raise ValueError('Invalid NCBI taxonomy response')
    taxa = {}
    for node in root.findall('Taxon'):
        taxid = int(node.findtext('TaxId'))
        lineage = [{'taxid': int(item.findtext('TaxId')),
                    'name': item.findtext('ScientificName'), 'rank': item.findtext('Rank')}
                   for item in node.findall('LineageEx/Taxon')]
        lineage.append({'taxid': taxid, 'name': node.findtext('ScientificName'),
                        'rank': node.findtext('Rank')})
        taxa[taxid] = {'taxid': taxid, 'name': node.findtext('ScientificName'),
                       'rank': node.findtext('Rank'), 'lineage': lineage}
    return taxa


def classify_marker(result, taxa, expected, min_identity=75, min_coverage=.8, near_best=.95):
    """Near-best taxonomic LCA; thresholds are explicit engineering heuristics."""
    hits = [h for h in result['hits'] if not h.get('blockers') and h.get('identity') is not None
            and h['identity'] >= min_identity and h['coverage'] >= min_coverage]
    if not hits:
        return {'status': 'insufficient', 'reason': 'No quality-qualified homology hit'}
    best = max(h['ranking_bitscore'] for h in hits)
    ambiguous = [h for h in result['hits'] if h.get('blockers')
                 and h['ranking_bitscore'] >= best * near_best]
    if ambiguous:
        return {'status': 'insufficient', 'reason': 'Strong competing hits cannot be reliably summarized',
                'excluded_ambiguous_hits': ambiguous}
    close = [h for h in hits if h['ranking_bitscore'] >= best * near_best]
    if len(close) >= 30:
        return {'status': 'insufficient', 'reason': 'Near-best hits saturate the requested hit-list limit',
                'near_best_hits': close, 'saturated_near_best': True}
    taxids = {taxid for h in close for taxid in h.get('taxids', [None])}
    if None in taxids or any(taxid not in taxa for taxid in taxids):
        return {'status': 'insufficient', 'reason': 'Near-best hit taxonomy incomplete',
                'near_best_hits': close}
    paths = [[node['taxid'] for node in taxa[taxid]['lineage']] for taxid in taxids]
    common = set(paths[0]).intersection(*map(set, paths[1:]))
    if not common:
        return {'status': 'insufficient', 'reason': 'No common taxonomy lineage'}
    lca_id = next(t for t in reversed(paths[0]) if t in common)
    exemplar = taxa[next(iter(taxids))]
    lca_node = next(node for node in exemplar['lineage'] if node['taxid'] == lca_id)
    lca_path = exemplar['lineage'][:exemplar['lineage'].index(lca_node) + 1]
    expected_ids = {node['taxid'] for node in expected['lineage']}
    if expected['taxid'] in {node['taxid'] for node in lca_path}:
        status = 'compatible'
    elif lca_id in expected_ids:
        status = 'insufficient'
    else:
        status = 'discordant'
    return {'status': status, 'provisional_lca': lca_node, 'near_best_hits': close,
            'reason': 'Near-best qualified hits share this taxonomic lineage; not species identification'}


def summarize(verdicts, mode, failures=()):
    statuses = {v['status'] for v in verdicts.values()}
    if {'compatible', 'discordant'} <= statuses:
        status = 'suspected_mixed'
    elif failures or not verdicts or 'insufficient' in statuses:
        status = 'insufficient'
    elif statuses == {'discordant'}:
        status = 'suspected_mismatch'
    elif statuses == {'compatible'}:
        status = 'compatible'
    else:
        status = 'insufficient'
    return {'status': status, 'mode': mode,
            'scope': 'Marker-origin consistency only; not specimen authentication or read abundance',
            'species_identified': False, 'company_error_established': False,
            'distinct_marker_types': sorted({v['marker'] for v in verdicts.values() if v.get('marker')}),
            'limitations': ['Multiple discovery branches or windows of one marker are not independent markers',
                            'Discordant queries flag candidate-origin conflicts; mixed specimens are not established',
                            'Database coverage and alternative recovery can affect marker classification'],
            'nuclear_marker_present': any(v.get('marker') == '18S' for v in verdicts.values()),
            'failures': list(failures), 'queries': verdicts}


class RemoteClient:
    """One durable request clock; also spaces taxonomy calls conservatively."""
    def __init__(self, directory, state, email):
        self.directory, self.state, self.email = directory, state, email

    def persist(self):
        save_json(self.directory / 'search_state.json', self.state)

    def request(self, url, params):
        wait = 10 - (time.time() - self.state.get('last_request', 0))
        if wait > 0:
            time.sleep(wait)
        self.state['last_request'] = time.time(); self.persist()
        params = dict(params, tool='mito-assembly-triage-identity')
        if self.email:
            params['email'] = self.email
        request = urllib.request.Request(url, data=urllib.parse.urlencode(params).encode(),
                                         headers={'User-Agent': 'mito-assembly-triage/identity'})
        with urllib.request.urlopen(request, timeout=180) as response:
            return response.read().decode('utf-8', 'strict')

    def taxonomy(self, ids):
        taxa = {}
        ids = sorted(set(ids))
        for offset in range(0, len(ids), 50):
            batch = ids[offset:offset + 50]
            raw = self.request(EUTILS + 'efetch.fcgi', {'db': 'taxonomy',
                               'id': ','.join(map(str, batch)), 'retmode': 'xml'})
            path = self.directory / ('taxonomy_%s.xml' % hashlib.sha256(','.join(map(str, batch)).encode()).hexdigest()[:12])
            path.write_text(raw)
            parsed = parse_taxonomy(raw)
            if set(parsed) != set(batch):
                raise ValueError('Incomplete taxonomy records')
            taxa.update(parsed)
        return taxa


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--markers', required=True)
    parser.add_argument('--outdir', required=True, help='New search directory; --resume to reuse')
    parser.add_argument('--allow-public-upload', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--email', default=os.environ.get('NCBI_EMAIL', ''), help='Optional NCBI contact email; also NCBI_EMAIL')
    parser.add_argument('--max-polls', type=int, default=20)
    parser.add_argument('--max-queries', type=int, default=20)
    args = parser.parse_args(argv)
    if not args.allow_public_upload:
        print('NCBI sequence upload needs --allow-public-upload', file=sys.stderr); return 2
    directory, marker_path = Path(args.outdir).resolve(), Path(args.markers).resolve()
    state = None
    lock = None
    writable = False
    try:
        if args.max_polls < 1 or not 1 <= args.max_queries <= 20 or (args.email and '@' not in args.email):
            raise ValueError('Invalid polling/query limit or contact email')
        manifest = json.loads(marker_path.read_text())
        if manifest.get('format') != 'mito-identity-markers-1' or manifest.get('status') not in ('prepared', 'insufficient'):
            raise ValueError('Marker preparation has not succeeded')
        if marker_path == directory or marker_path in directory.parents:
            raise ValueError('Output overlaps marker input')
        for entry in manifest['inputs']:
            if sha256(entry['path']) != entry['sha256']:
                raise ValueError('Source input changed since marker preparation')
        queries = {}
        for query in manifest['queries']:
            path = Path(query['path'])
            if sha256(path) != query['sha256']:
                raise ValueError('Marker query changed after preparation')
            records = fasta_records(path, 'ACGT')
            if len(records) != 1 or records[0][0] != query['id'] or len(records[0][1]) != query['length']:
                raise ValueError('Query identifier or length does not match manifest')
            if query['id'] in queries or not 100 <= query['length'] <= 5000:
                raise ValueError('Duplicate or oversized query')
            queries[query['id']] = dict(query, sequence=records[0][1])
        if len(queries) > args.max_queries:
            raise ValueError('Too many queries; review preparation before upload, no silent truncation')
        fingerprint = sha256(marker_path)
        if directory.exists():
            if not args.resume:
                raise ValueError('Search directory exists; use --resume')
        else:
            if args.resume:
                raise ValueError('Cannot resume missing search directory')
            directory.mkdir(parents=True)
            state = {'format': 'mito-identity-search-1', 'markers_sha256': fingerprint,
                     'expected_taxon': manifest['expected_taxon'], 'status': 'new',
                     'database_requested': 'core_nt',
                     'created_at': dt.datetime.now(dt.timezone.utc).isoformat(),
                     'thresholds': {'min_identity': 75, 'min_coverage': .8, 'near_best_bitscore': .95,
                                    'ranking': 'sum of non-overlapping collinear HSP bit scores',
                                    'kind': 'engineering heuristics, not species identification cutoffs'}}
        client = RemoteClient(directory, state, args.email)
        lock = (directory / '.search.lock').open('a')
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            state = None
            raise ValueError('Another process is running this search; do not resubmit')
        if (directory / 'search_state.json').exists():
            state = json.loads((directory / 'search_state.json').read_text())
            if state.get('markers_sha256') != fingerprint:
                raise ValueError('Resume marker manifest changed')
            client.state = state
        if state is None:
            raise ValueError('Cannot resume a directory without search_state.json')
        if state.get('status') == 'completed':
            result_path = directory / 'identity_result.json'
            if sha256(result_path) != state.get('result_sha256'):
                raise ValueError('Completed result was changed')
            result = json.loads(result_path.read_text())
            print(json.dumps({'status': result['status'], 'result': str(result_path)}))
            return 1 if result['status'] == 'insufficient' else 0
        result_file = directory / 'blast.xml'
        if result_file.exists() and sha256(result_file) != state.get('blast_xml_sha256'):
            raise ValueError('Stored BLAST response changed or lacks a bound checksum')
        writable = True
        client.persist()
        conflicts = ['Preparation has unresolved mixed columns: %s' % item['alignment']
                     for item in manifest.get('observations', []) if item.get('mixed_columns', 0) > 0]
        if not queries:
            result_path = directory / 'identity_result.json'
            result = summarize({}, manifest['mode'], conflicts)
            result['marker_manifest_sha256'] = fingerprint
            save_json(result_path, result)
            state.update(status='completed', result_sha256=sha256(result_path))
            client.persist()
            return 1
        if state.get('status') == 'submitting':
            raise RuntimeError('Previous submission outcome unknown; inspect raw submission before retrying')
        if not state.get('rid'):
            state['status'] = 'submitting'; client.persist()
            raw = client.request(cox1_id.BLAST_URL, {'CMD': 'Put', 'PROGRAM': 'blastn',
                'DATABASE': 'core_nt', 'QUERY': ''.join('>%s\n%s\n' % (name, query['sequence']) for name, query in queries.items()),
                'HITLIST_SIZE': '30', 'EXPECT': '1e-10', 'WORD_SIZE': '11'})
            (directory / 'submission.txt').write_text(raw)
            _, rid = cox1_id.parse_search_info(raw)
            if not rid:
                raise ValueError('NCBI did not return a RID')
            state.update(rid=rid, status='waiting', last_poll=time.time()); client.persist()
        if not result_file.exists():
            ready = False
            for _ in range(args.max_polls):
                delay = 60 - (time.time() - state.get('last_poll', 0))
                if delay > 0:
                    time.sleep(delay)
                state['last_poll'] = time.time(); client.persist()
                raw = client.request(cox1_id.BLAST_URL, {'CMD': 'Get', 'FORMAT_OBJECT': 'SearchInfo', 'RID': state['rid']})
                status, _ = cox1_id.parse_search_info(raw)
                (directory / 'poll.txt').write_text(raw)
                if status not in ('READY', 'WAITING', 'FAILED', 'UNKNOWN'):
                    raise ValueError('NCBI polling response has no recognized search status')
                if status == 'READY':
                    ready = True; break
                if status in ('FAILED', 'UNKNOWN'):
                    raise RuntimeError('NCBI RID %s is %s; not a taxonomic verdict' % (state['rid'], status))
            if not ready:
                state['status'] = 'waiting'; client.persist()
                print('NCBI search still pending; reuse --resume', file=sys.stderr); return 3
            raw = client.request(cox1_id.BLAST_URL, {'CMD': 'Get', 'FORMAT_TYPE': 'XML2_S', 'RID': state['rid'],
                                                   'ALIGNMENTS': '30', 'DESCRIPTIONS': '30'})
            # Save even malformed results for diagnosis; never turn them into no-match.
            result_file.write_text(raw)
            state['blast_xml_sha256'] = sha256(result_file); client.persist()
        parsed = parse_batch_xml(result_file.read_text(), queries)
        save_json(directory / 'parsed_hits.json', parsed)
        expected_name = ALIASES.get(manifest['expected_taxon'].lower(), manifest['expected_taxon'])
        raw = client.request(EUTILS + 'esearch.fcgi', {'db': 'taxonomy',
                               'term': expected_name + '[Scientific Name]', 'retmode': 'json'})
        ids = json.loads(raw)['esearchresult']['idlist']
        if len(ids) != 1:
            raise ValueError('Expected taxon cannot be uniquely resolved: %s' % expected_name)
        taxids = {int(ids[0])}
        for result in parsed.values():
            taxids.update(t for h in result['hits'] for t in h.get('taxids', []) if t)
        taxa = client.taxonomy(taxids)
        expected = taxa[int(ids[0])]
        verdicts = {name: dict(classify_marker(result, taxa, expected), marker=queries[name]['marker'],
                              provenance=queries[name]['provenance'], database_actual=result['database_actual'])
                    for name, result in parsed.items()}
        result = summarize(verdicts, manifest['mode'], conflicts)
        result.update(expected_taxon=expected, rid=state['rid'], thresholds=state['thresholds'],
                      marker_manifest_sha256=fingerprint)
        save_json(directory / 'identity_result.json', result)
        state['status'] = 'completed'; state['result_sha256'] = sha256(directory / 'identity_result.json'); client.persist()
        print(json.dumps({'status': result['status'], 'result': str(directory / 'identity_result.json')}, ensure_ascii=False))
        return 1 if result['status'] == 'insufficient' else 0
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, IndexError, ET.ParseError) as exc:
        if writable:
            state['error'] = str(exc); save_json(directory / 'search_state.json', state)
            save_json(directory / 'identity_result.json', {'status': 'service_or_format_failure',
                       'scientific_verdict': None, 'error': str(exc)})
        print(str(exc), file=sys.stderr)
        return 3 if writable else 2
    finally:
        if lock is not None:
            lock.close()


if __name__ == '__main__':
    raise SystemExit(main())

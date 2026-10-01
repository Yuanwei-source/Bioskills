"""Identity-screen scientific counterexamples and input/network contracts."""
import copy
import importlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
markers = importlib.import_module('_identity_markers')
search = importlib.import_module('identity_search')
prepare = importlib.import_module('identity_prepare')


def taxonomy(taxid, name, nodes):
    return {'taxid': taxid, 'name': name, 'rank': 'species',
            'lineage': [{'taxid': n, 'name': str(n), 'rank': 'no rank'} for n in nodes]}


EXPECTED = taxonomy(50557, 'Insecta', [1, 2759, 33208, 6656, 50557])
TAXA = {7: taxonomy(7, 'insect A', [1, 2759, 33208, 6656, 50557, 7]),
        8: taxonomy(8, 'insect B', [1, 2759, 33208, 6656, 50557, 8]),
        9606: taxonomy(9606, 'Homo sapiens', [1, 2759, 33208, 7711, 9606])}


def hit(taxids, score=500, coverage=1, identity=99):
    return {'taxids': taxids, 'ranking_bitscore': score, 'coverage': coverage,
            'identity': identity, 'blockers': []}


def xml(description_taxids=(7,), title='query_001', hsps=None):
    if hsps is None:
        hsps = [(1, 300, 1, 300, 500)]
    descriptions = ''.join('<HitDescr><id>A%d</id><accession>A%d.1</accession><taxid>%d</taxid>'
                           '<title>record</title></HitDescr>' % (i, i, t)
                           for i, t in enumerate(description_taxids))
    hsp = ''.join('<Hsp><bit-score>%s</bit-score><identity>%d</identity><align-len>%d</align-len>'
                  '<query-from>%d</query-from><query-to>%d</query-to><hit-from>%d</hit-from>'
                  '<hit-to>%d</hit-to></Hsp>' % (score, qe-qs+1, qe-qs+1, qs, qe, hs, he)
                  for qs, qe, hs, he, score in hsps)
    return ('<BlastXML2><Search><query-title>%s</query-title><query-len>300</query-len>'
            '<hits><Hit><description>%s</description><len>1000</len><hsps>%s</hsps>'
            '</Hit></hits></Search></BlastXML2>') % (title, descriptions, hsp)


class TaxonomicScreenTests(unittest.TestCase):
    def marker_fixture(self, root):
        query = root / 'q.fasta'; markers.write_fasta(query, [('query_001', 'A'*300)])
        manifest = root / 'markers.json'
        markers.save_json(manifest, {'format': 'mito-identity-markers-1', 'status': 'prepared',
            'expected_taxon': 'Insecta', 'mode': 'assembly', 'inputs': [],
            'queries': [{'id': 'query_001', 'length': 300, 'marker': 'COX1', 'path': str(query),
                         'sha256': markers.sha256(query), 'provenance': []}]})
        return manifest

    def test_remote_roundtrip_resume_uses_bound_result_without_resubmit(self):
        def request(client, url, params):
            if params.get('CMD') == 'Put':
                return 'RID = TEST123\nRTOE = 1\n'
            if params.get('FORMAT_OBJECT') == 'SearchInfo':
                return 'Status=READY\n'
            if params.get('FORMAT_TYPE') == 'XML2_S':
                return xml()
            if url.endswith('esearch.fcgi'):
                return '{"esearchresult":{"idlist":["50557"]}}'
            if url.endswith('efetch.fcgi'):
                nodes = []
                for taxid in params['id'].split(','):
                    item = EXPECTED if taxid == '50557' else TAXA[int(taxid)]
                    lin = ''.join('<Taxon><TaxId>%d</TaxId><ScientificName>%s</ScientificName>'
                                  '<Rank>no rank</Rank></Taxon>' % (n['taxid'], n['name'])
                                  for n in item['lineage'][:-1])
                    nodes.append('<Taxon><TaxId>%d</TaxId><ScientificName>%s</ScientificName>'
                                 '<Rank>species</Rank><LineageEx>%s</LineageEx></Taxon>'
                                 % (item['taxid'], item['name'], lin))
                return '<TaxaSet>'+''.join(nodes)+'</TaxaSet>'
            raise AssertionError(params)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); manifest = self.marker_fixture(root); out = root/'online'
            argv = ['--markers', str(manifest), '--outdir', str(out), '--allow-public-upload']
            with patch.object(search.RemoteClient, 'request', request), patch.object(search.time, 'sleep'):
                self.assertEqual(search.main(argv), 0)
            result = json.loads((out/'identity_result.json').read_text())
            self.assertEqual(result['status'], 'compatible')
            state = json.loads((out/'search_state.json').read_text())
            self.assertEqual(state['blast_xml_sha256'], markers.sha256(out/'blast.xml'))
            with patch.object(search.RemoteClient, 'request') as remote:
                self.assertEqual(search.main(argv+['--resume']), 0)
                remote.assert_not_called()

    def test_uncertain_submission_is_not_automatically_uploaded_twice(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); manifest = self.marker_fixture(root); out = root/'online'
            argv = ['--markers', str(manifest), '--outdir', str(out), '--allow-public-upload']
            with patch.object(search.RemoteClient, 'request', side_effect=OSError('connection lost')):
                self.assertEqual(search.main(argv), 3)
            with patch.object(search.RemoteClient, 'request') as remote:
                self.assertEqual(search.main(argv+['--resume']), 3)
                remote.assert_not_called()
            self.assertEqual(json.loads((out/'identity_result.json').read_text())['scientific_verdict'], None)

    def test_query_drift_prevents_upload(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); manifest = self.marker_fixture(root)
            (root/'q.fasta').write_text('>query_001\n'+'G'*300+'\n')
            with patch.object(search.RemoteClient, 'request') as remote:
                code = search.main(['--markers', str(manifest), '--outdir', str(root/'online'),
                                    '--allow-public-upload'])
            self.assertEqual(code, 2)
            remote.assert_not_called()
            self.assertFalse((root/'online').exists())

    def test_near_best_insects_support_class_not_species(self):
        result = search.classify_marker({'hits': [hit([7]), hit([8], 490)]}, TAXA, EXPECTED)
        self.assertEqual(result['status'], 'compatible')
        self.assertEqual(result['provisional_lca']['taxid'], 50557)

    def test_strong_ambiguous_hit_prevents_weaker_directional_verdict(self):
        ambiguous = hit([7], 700); ambiguous['blockers'] = ['overlapping_hsps']
        verdict = search.classify_marker({'hits': [ambiguous, hit([9606], 500)]}, TAXA, EXPECTED)
        self.assertEqual(verdict['status'], 'insufficient')
        self.assertEqual(verdict['excluded_ambiguous_hits'], [ambiguous])

    def test_saturated_near_best_does_not_certify_class(self):
        verdict = search.classify_marker({'hits': [hit([7])] * 30}, TAXA, EXPECTED)
        self.assertEqual(verdict['status'], 'insufficient')

    def test_invalid_resume_preserves_existing_outputs_even_with_active_lock(self):
        import fcntl
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); manifest = self.marker_fixture(root); out = root/'online'; out.mkdir()
            markers.save_json(out/'search_state.json', {'markers_sha256': 'another manifest', 'status': 'completed'})
            (out/'identity_result.json').write_text('{"status":"compatible"}')
            originals = {p: p.read_bytes() for p in (out/'search_state.json', out/'identity_result.json')}
            argv = ['--markers', str(manifest), '--outdir', str(out), '--allow-public-upload', '--resume']
            with patch.object(search.RemoteClient, 'request') as remote:
                self.assertEqual(search.main(argv), 2)
                with (out/'.search.lock').open('a') as lock:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    self.assertEqual(search.main(argv), 2)
                remote.assert_not_called()
            for path, content in originals.items():
                self.assertEqual(path.read_bytes(), content)

    def test_resume_directory_without_state_is_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); manifest = self.marker_fixture(root); out = root/'online'; out.mkdir()
            keep = out/'keep.txt'; keep.write_text('original')
            self.assertEqual(search.main(['--markers', str(manifest), '--outdir', str(out),
                                         '--allow-public-upload', '--resume']), 2)
            self.assertEqual(keep.read_text(), 'original')
            self.assertFalse((out/'identity_result.json').exists())

    def test_no_queries_persists_insufficiency_and_conflict_without_network(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); manifest = self.marker_fixture(root); out = root/'online'
            data = json.loads(manifest.read_text()); data.update(status='insufficient', queries=[],
                observations=[{'alignment': 'mixed_alignment.fas', 'mixed_columns': 5}])
            markers.save_json(manifest, data)
            argv = ['--markers', str(manifest), '--outdir', str(out), '--allow-public-upload']
            with patch.object(search.RemoteClient, 'request') as remote:
                self.assertEqual(search.main(argv), 1)
                before = (out/'identity_result.json').read_bytes()
                self.assertEqual(search.main(argv + ['--resume']), 1)
                self.assertEqual((out/'identity_result.json').read_bytes(), before)
                remote.assert_not_called()
            self.assertEqual(len(json.loads(before)['failures']), 1)

    def test_grouped_identical_sequences_keep_all_taxa(self):
        result = search.parse_batch_xml(xml((7, 9606)), {'query_001': {'length': 300}})
        hits = result['query_001']['hits']
        self.assertEqual(set(hits[0]['taxids']), {7, 9606})
        verdict = search.classify_marker(result['query_001'], TAXA, EXPECTED)
        self.assertEqual(verdict['status'], 'insufficient')
        self.assertEqual(verdict['provisional_lca']['taxid'], 33208)

    def test_split_hsp_scores_are_summed(self):
        result = search.parse_batch_xml(xml(hsps=[(1, 150, 1, 150, 250),
                                                  (151, 300, 501, 650, 250)]),
                                        {'query_001': {'length': 300}})
        self.assertEqual(result['query_001']['hits'][0]['ranking_bitscore'], 500)

    def test_missing_taxonomy_and_low_coverage_are_inconclusive(self):
        for candidate in (hit([None]), hit([9999]), hit([7], coverage=.2)):
            with self.subTest(candidate=candidate):
                self.assertEqual(search.classify_marker({'hits': [candidate]}, TAXA, EXPECTED)['status'],
                                 'insufficient')

    def test_mammalian_marker_is_mismatch_not_company_error(self):
        verdict = search.classify_marker({'hits': [hit([9606])]}, TAXA, EXPECTED)
        summary = search.summarize({'q': verdict}, 'reads')
        self.assertEqual(summary['status'], 'suspected_mismatch')
        self.assertFalse(summary['company_error_established'])
        self.assertFalse(summary['species_identified'])

    def test_mixed_signals_are_not_erased_by_an_ambiguous_query(self):
        summary = search.summarize({'a': {'status': 'compatible'},
                                    'b': {'status': 'discordant'}, 'c': {'status': 'insufficient'}}, 'reads')
        self.assertEqual(summary['status'], 'suspected_mixed')

    def test_preparation_mixture_blocks_clean_compatibility(self):
        summary = search.summarize({'a': {'status': 'compatible'}}, 'reads', ['mixed columns'])
        self.assertEqual(summary['status'], 'insufficient')

    def test_missing_extra_duplicate_query_or_bad_xml_is_format_failure(self):
        cases = [('<html>Error</html>', {'query_001': {'length': 300}}),
                 (xml(), {'query_001': {'length': 301}}),
                 (xml(), {'query_001': {'length': 300}, 'query_002': {'length': 300}}),
                 (xml(title='wrong'), {'query_001': {'length': 300}})]
        for payload, queries in cases:
            with self.subTest(payload=payload):
                with self.assertRaises(ValueError):
                    search.parse_batch_xml(payload, queries)

    def test_upload_without_authorization_makes_no_network_call(self):
        with patch.object(search.urllib.request, 'urlopen') as remote:
            code = search.main(['--markers', '/missing', '--outdir', '/missing-output'])
        self.assertEqual(code, 2)
        remote.assert_not_called()


class MarkerSupportTests(unittest.TestCase):
    def alignment(self, rows):
        directory = tempfile.TemporaryDirectory(); self.addCleanup(directory.cleanup)
        path = Path(directory.name) / 'alignment.fas'
        markers.write_fasta(path, rows)
        return path

    def test_overlapping_mates_are_not_two_templates(self):
        path = self.alignment([('L0_P1/1', 'A'*300), ('L0_P1/2', 'A'*300)])
        sequences, qc = markers.supported_segments(path, 'reads')
        self.assertEqual(sequences, [])
        self.assertEqual(max(qc['column_depth']), 1)

    def test_supported_template_consensus_and_mixture(self):
        path = self.alignment([('L0_P%d/1' % i, 'A'*300) for i in range(5)])
        self.assertEqual(markers.supported_segments(path, 'reads')[0], ['A'*300])
        path = self.alignment([('L0_P%d/1' % i, ('A' if i < 5 else 'G')*300) for i in range(10)])
        sequences, qc = markers.supported_segments(path, 'reads')
        self.assertEqual(sequences, [])
        self.assertEqual(qc['mixed_columns'], 300)

    def test_gaps_do_not_join_disjoint_regions(self):
        path = self.alignment([('source_0', 'A'*200+'~'*3+'T'*200)])
        self.assertEqual(markers.supported_segments(path, 'fasta')[0], [])

    def test_mate_disagreement_does_not_support_either_base(self):
        path = self.alignment([('L0_P1/1', 'A'*300), ('L0_P1/2', 'G'*300)])
        self.assertEqual(max(markers.supported_segments(path, 'reads')[1]['column_depth']), 0)

    def test_sampling_preserves_lanes_and_validates_entire_inputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); lanes = []
            for lane in range(2):
                paths = [root / ('lane%d_%d.fastq' % (lane, mate)) for mate in (1, 2)]
                for mate, path in enumerate(paths, 1):
                    path.write_text(''.join('@p%d/%d\nACGT\n+\nIIII\n' % (i, mate) for i in range(20)))
                lanes.append(paths)
            paths, stats = markers.sample_pairs(lanes, 6, 19, root)
            self.assertEqual(stats['sampled_per_lane'], [3, 3])
            self.assertEqual(len(list(markers.paired_records(*paths))), 6)
            lanes[0][1].write_text(lanes[0][1].read_text() + '@bad/2\nACGT\n+\nI\n')
            with self.assertRaises(ValueError):
                markers.sample_pairs(lanes, 6, 19, root)

    def test_annotation_sequence_mismatch_is_refused(self):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'wrong.gb'
            record = SeqRecord(Seq('AAAA'), id='sample', annotations={'molecule_type': 'DNA'})
            SeqIO.write(record, path, 'genbank')
            with self.assertRaises(ValueError):
                markers.check_annotation(path, [('sample', 'TTTT')])

    def test_duplicate_lane_paths_are_refused_before_output(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); left, right = root/'r1.fastq', root/'r2.fastq'
            left.write_text('@p/1\nACGT\n+\nIIII\n'); right.write_text('@p/2\nACGT\n+\nIIII\n')
            out = root / 'out'
            code = prepare.main(['--r1', str(left), '--r2', str(right), '--r1', str(left),
                                  '--r2', str(right), '--outdir', str(out)])
            self.assertEqual(code, 2)
            self.assertFalse(out.exists())

    def test_copied_lanes_do_not_inflate_template_support(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); argv = []
            for lane in range(2):
                for mate in (1, 2):
                    path = root/('lane%d_R%d.fastq' % (lane, mate))
                    path.write_text('@p/%d\nACGT\n+\nIIII\n' % mate)
                    argv += ['--r%d' % mate, str(path)]
            out = root/'out'
            self.assertEqual(prepare.main(argv + ['--outdir', str(out)]), 2)
            self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()

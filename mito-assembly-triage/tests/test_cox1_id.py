#!/usr/bin/env python3
"""Regression coverage for the NCBI remote COX1 BLAST response contract."""
import importlib.util
import pathlib
import unittest
from unittest.mock import patch


ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('cox1_id_under_test', ROOT / 'scripts' / 'cox1_id.py')
COX1 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COX1)


XML2_SINGLE = '''<?xml version="1.0"?>
<BlastXML2 xmlns="http://www.ncbi.nlm.nih.gov">
  <BlastOutput2><report><Report><results><Results><search><Search>
    <query-len>800</query-len>
    <hits><Hit>
      <description><HitDescr><id>ACC000001.1</id><title>Test insect mitochondrion</title></HitDescr></description>
      <len>1200</len>
      <hsps><Hsp>
        <bit-score>1478.44</bit-score><identity>800</identity>
        <query-from>1</query-from><query-to>800</query-to>
        <query-strand>Plus</query-strand><hit-from>201</hit-from><hit-to>1000</hit-to>
        <hit-strand>Plus</hit-strand><align-len>800</align-len>
      </Hsp></hsps>
    </Hit></hits>
  </Search></search></Results></results></Report></report></BlastOutput2>
</BlastXML2>'''


class RemoteResponseContract(unittest.TestCase):
    def test_fetch_requests_single_file_xml2(self):
        observed = {}

        def fake_request(params, timeout=120):
            observed.update(params)
            observed['timeout'] = timeout
            return XML2_SINGLE

        with patch.object(COX1, '_request', side_effect=fake_request):
            self.assertEqual(COX1.fetch_results_xml('RID123'), XML2_SINGLE)

        self.assertEqual(observed['FORMAT_TYPE'], 'XML2_S')
        self.assertEqual(observed['CMD'], 'Get')
        self.assertEqual(observed['RID'], 'RID123')

    def test_single_file_xml2_response_is_parsed_as_a_candidate(self):
        result = COX1.parse_blast_xml(XML2_SINGLE)
        self.assertEqual(result['query_len'], 800)
        self.assertEqual(len(result['hits']), 1)
        hit = result['hits'][0]
        self.assertEqual(hit['accession'], 'ACC000001.1')
        self.assertEqual(hit['identity'], 100.0)
        self.assertEqual(hit['coverage'], 1.0)

    def test_remote_poll_interval_obeys_ncbi_one_minute_guidance(self):
        self.assertGreaterEqual(COX1.POLL_INTERVAL, 60)


if __name__ == '__main__':
    unittest.main()

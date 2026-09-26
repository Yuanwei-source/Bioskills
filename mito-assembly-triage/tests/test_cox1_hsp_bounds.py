#!/usr/bin/env python3
"""#5：BLAST HSP 的坐标与数值边界校验（mutation 式）。

策略（本文件即验收标准）：

| 情况 | 处置 | 理由 |
|---|---|---|
| 坐标越界（`query-from/to` > query-len、`hit-from/to` > Hit_len） | `ValueError` → **退出码 3** | 同一份回复内部自相矛盾 → 不能信任其中任何 HSP；与"半可读 HSP 即格式故障"的既有立场一致 |
| 长度字段缺失（`query-len` / `Hit_len` 不存在） | **不报错**，但记录"未校验" | 无法检查不等于数据错误；必须让"未校验"可见，而不是默认它正确 |
| 非有限 `bit-score`（nan/inf） | `ValueError` → 3 | 数值本身不可用 |
| `align-len` 小于 query 跨度 | `ValueError` → 3 | 带 gap 时 align-len ≥ 跨度，反之为内部矛盾 |
| `identity == 0` 且 `align-len > 0` | **blocker**（`zero_identity`），不退出 3 | 这是"弱到无用的命中"，不是格式问题；它不得被汇总成 candidate |
| `*-strand` / `*-frame` 文本与坐标方向不符 | `ValueError` → 3 | 同一事实的两种表述互相矛盾 |

**退出码边界必须保持**：格式故障 = 3，`no_match`/`insufficient` = 1。
"""
import contextlib
import importlib.util
import io
import json
import pathlib
import sys
import tempfile
import unittest
from types import SimpleNamespace

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIXTURES = ROOT / 'tests' / 'fixtures'
XML2 = FIXTURES / 'blast_xml2_outfmt16.xml'
LEGACY = FIXTURES / 'blast_xml_outfmt5.xml'


def load_module(name='cox1_bounds'):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / 'cox1_id.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DEFAULT_HSP = {
    'bit-score': '1478.44', 'identity': '800', 'align-len': '800',
    'query-from': '1', 'query-to': '800', 'hit-from': '201', 'hit-to': '1000',
    'query-strand': 'Plus', 'hit-strand': 'Plus',
}


def xml2_payload(hsp=(), query_len=1200, hit_len=1200, strand_style='strand'):
    """构造最小 XML2 载荷；hsp 覆盖默认字段，值设为 None 表示删掉该字段。"""
    fields = dict(DEFAULT_HSP)
    fields.update(dict(hsp))
    body = []
    for key, value in fields.items():
        if value is None:
            continue
        body.append('<%s>%s</%s>' % (key, value, key))
    if strand_style == 'frame':
        body = [item for item in body if 'strand' not in item]
        body.append('<hit-frame>%s</hit-frame>' % ('1' if fields.get('hit-strand') == 'Plus' else '-1'))
    hit_len_tag = '' if hit_len is None else '<len>%d</len>' % hit_len
    query_len_tag = '' if query_len is None else '<query-len>%d</query-len>' % query_len
    return ('<BlastXML2 xmlns="http://www.ncbi.nlm.nih.gov">'
            '%s<hits><Hit><num>1</num>'
            '<description><HitDescr><id>gnl|BL_ORD_ID|0</id>'
            '<title>synthetic subject</title></HitDescr></description>'
            '%s<hsps><Hsp>%s</Hsp></hsps></Hit></hits></BlastXML2>'
            % (query_len_tag, hit_len_tag, ''.join(body)))


class _Cli(unittest.TestCase):
    """CLI 级：格式故障必须保持退出码 3（不得伪装成 no_match=1）。"""

    def setUp(self):
        self.module = load_module('cox1_cli_bounds')
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.dir = pathlib.Path(temp.name)
        self.fasta = self.dir / 'genome.fa'
        self.fasta.write_text('>g\n%s\n' % ('ACGT' * 300), encoding='utf-8')

    def run_cli(self, payload, extra=()):
        module = self.module
        saved = (module.request_search, module.poll_search_info,
                 module.fetch_results_xml, module.time)
        module.request_search = lambda query, max_results=10: 'RID = RID1\nRTOE = 0'
        module.poll_search_info = lambda rid: 'Status=READY'
        module.fetch_results_xml = lambda rid, max_results=10: payload
        module.time = SimpleNamespace(sleep=lambda seconds: None)
        argv = sys.argv
        sys.argv = ['cox1_id.py', str(self.fasta), '--allow-public-upload',
                    '--coords', '1,1200'] + list(extra)
        try:
            with contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    module.main()
            return caught.exception.code
        finally:
            (module.request_search, module.poll_search_info,
             module.fetch_results_xml, module.time) = saved
            sys.argv = argv


class ParserBoundsTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module('cox1_bounds_parse')

    def parse(self, payload):
        return self.module.parse_blast_xml(payload)

    # --- 方言基线：真实 blastn 输出仍须可解析 -----------------------------
    def test_real_fixtures_still_parse(self):
        for fixture in (XML2, LEGACY):
            with self.subTest(fixture=fixture.name):
                parsed = self.parse(fixture.read_text(encoding='utf-8'))
                self.assertTrue(parsed['hits'])
                self.assertTrue(parsed['hits'][0]['ranges_checked'],
                                '真实 fixture 带长度字段，范围检查应当真的执行过')

    # --- 越界即格式故障 ---------------------------------------------------
    def test_query_coordinate_beyond_query_len_is_a_format_failure(self):
        payload = xml2_payload({'query-to': '1600'}, query_len=1200)
        with self.assertRaises(ValueError) as caught:
            self.parse(payload)
        self.assertIn('query', str(caught.exception))

    def test_hit_coordinate_beyond_hit_len_is_a_format_failure(self):
        payload = xml2_payload({'hit-to': '5000'}, hit_len=1200)
        with self.assertRaises(ValueError) as caught:
            self.parse(payload)
        self.assertIn('hit', str(caught.exception))

    def test_zero_length_fields_are_format_failures(self):
        with self.assertRaises(ValueError):
            self.parse(xml2_payload(query_len=0))

    # --- 长度字段缺失：不报错，但必须标记"未校验" ------------------------
    def test_missing_length_fields_skip_the_range_check_but_record_it(self):
        # 长度字段缺失时范围检查无法执行（但 HSP 本身自洽，避免触发跨度规则）
        parsed = self.parse(xml2_payload(query_len=None, hit_len=None))
        hit = parsed['hits'][0]
        self.assertFalse(hit['ranges_checked'])
        self.assertIn('query_len_absent', hit['unchecked'])
        self.assertIn('hit_len_absent', hit['unchecked'])
        self.assertIsNone(parsed['query_len'])

    # --- 数值边界 ---------------------------------------------------------
    def test_non_finite_bit_score_is_a_format_failure(self):
        for value in ('nan', 'inf', '-inf'):
            with self.subTest(bitscore=value):
                with self.assertRaises(ValueError):
                    self.parse(xml2_payload({'bit-score': value}))

    def test_align_len_smaller_than_the_query_span_is_a_format_failure(self):
        payload = xml2_payload({'query-from': '1', 'query-to': '800', 'align-len': '500'})
        with self.assertRaises(ValueError) as caught:
            self.parse(payload)
        self.assertIn('align-len', str(caught.exception))

    def test_align_len_larger_than_the_span_is_allowed_with_gaps(self):
        parsed = self.parse(xml2_payload({'align-len': '900', 'gaps': '100'}))
        self.assertEqual(parsed['hits'][0]['hsp_count'], 1)

    # --- 弱命中是 blocker，不是格式故障 ----------------------------------
    def test_zero_identity_is_a_blocker_not_a_format_failure(self):
        parsed = self.parse(xml2_payload({'identity': '0'}))
        hit = parsed['hits'][0]
        self.assertIn('zero_identity', hit['ambiguity_reasons'])
        self.assertIsNone(hit['identity'], '无信息的命中不得被汇总成一个 identity')

    # --- 方向文本与坐标方向必须自洽 --------------------------------------
    def test_strand_text_contradicting_the_coordinates_is_a_format_failure(self):
        payload = xml2_payload({'hit-strand': 'Minus'})   # hit-from < hit-to 却是 Minus
        with self.assertRaises(ValueError) as caught:
            self.parse(payload)
        self.assertIn('strand', str(caught.exception).lower())

    def test_matching_strands_are_accepted(self):
        plus = self.parse(xml2_payload())
        self.assertEqual(plus['hits'][0]['hsp_count'], 1)
        minus = self.parse(xml2_payload({'hit-from': '1000', 'hit-to': '201',
                                         'hit-strand': 'Minus'}))
        self.assertEqual(minus['hits'][0]['hsp_count'], 1)

    def test_absent_strand_fields_are_recorded_as_unchecked(self):
        parsed = self.parse(xml2_payload({'query-strand': None, 'hit-strand': None}))
        self.assertIn('strand_absent', parsed['hits'][0]['unchecked'])

    def test_frame_fields_are_accepted_for_the_legacy_dialect(self):
        payload = xml2_payload({'hit-strand': 'Plus'}, strand_style='frame')
        self.assertEqual(self.parse(payload)['hits'][0]['hsp_count'], 1)
        bad = xml2_payload({'hit-strand': 'Minus'}, strand_style='frame')  # frame=-1 但坐标递增
        with self.assertRaises(ValueError):
            self.parse(bad)


class ExitCodeBoundaryTests(_Cli):
    def test_bounds_violation_is_exit_3_not_no_match(self):
        payload = xml2_payload({'hit-to': '99999'}, hit_len=1200)
        self.assertEqual(self.run_cli(payload), 3)

    def test_real_fixture_still_reaches_a_verdict(self):
        self.assertEqual(self.run_cli(XML2.read_text(encoding='utf-8')), 0)

    def test_output_json_keeps_the_raw_hsp_evidence_on_a_blocker(self):
        out = self.dir / 'raw.json'
        code = self.run_cli(xml2_payload({'identity': '0'}), extra=['--output-json', str(out)])
        self.assertIn(code, (1, 2))
        data = json.loads(out.read_text(encoding='utf-8'))
        blob = json.dumps(data, ensure_ascii=False)
        self.assertIn('zero_identity', blob, '被阻断的命中仍须保留原始 HSP 与阻断原因')
        self.assertIn('bitscore', blob)
        self.assertIn('"hsps"', blob, '原始 HSP 必须留在 --output-json 里供人工核对')


if __name__ == '__main__':
    unittest.main()

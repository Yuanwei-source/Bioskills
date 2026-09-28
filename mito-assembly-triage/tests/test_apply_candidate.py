import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/apply_candidate.py'


class CandidateApplication(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def run_cli(self, *args, code=0):
        result = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True)
        self.assertEqual(result.returncode, code, result.stdout + result.stderr)
        return result

    def json_file(self, name, value):
        path = self.root / name
        path.write_text(json.dumps(value))
        return path

    def test_base_edit_records_exact_change_and_never_overwrites_source(self):
        source = self.root / 'source.fa'; source.write_text('>mt demo\nACGTN\n')
        edits = self.json_file('edits.json', {'changes': [
            {'sequence_id': 'mt', 'position': 4, 'old': 'T', 'new': 'C', 'evidence': 'E1/reads', 'rationale': 'read support'}]})
        output, manifest = self.root / 'candidate.fa', self.root / 'manifest.json'
        self.run_cli('base', source, edits, output, '--manifest', manifest)
        self.assertEqual(source.read_text(), '>mt demo\nACGTN\n')
        self.assertEqual(output.read_text(), '>mt demo\nACGCN\n')
        record = json.loads(manifest.read_text())
        self.assertEqual(record['changes'][0]['position'], 4)
        self.assertEqual(record['scientific_status'], 'candidate_only_not_validated')

    def test_base_edit_targets_one_record_in_multi_fasta(self):
        source = self.root/'multi.fa'; source.write_text('>mt1 first\nACGT\n>mt2 second\nTTAA\n')
        edits = self.json_file('multi-edits.json', {'changes': [
            {'sequence_id':'mt2','position':2,'old':'T','new':'C','evidence':'E2/reads','rationale':'unique read support'}]})
        output = self.root/'multi-candidate.fa'; manifest = self.root/'multi-manifest.json'
        self.run_cli('base', source, edits, output, '--manifest', manifest)
        with output.open() as handle:
            records = list(SeqIO.parse(handle, 'fasta'))
        self.assertEqual([(r.id, str(r.seq)) for r in records], [('mt1','ACGT'),('mt2','TCAA')])

    def test_base_edit_can_resolve_ambiguous_source_with_explicit_evidence(self):
        source = self.root/'ambiguous.fa'; source.write_text('>mt\nACNT\n')
        edits = self.json_file('ambiguous-edits.json', {'changes': [
            {'sequence_id':'mt','position':3,'old':'N','new':'G','evidence':'E1/reads','rationale':'independent high-quality reads support G'}]})
        output = self.root/'resolved.fa'; manifest = self.root/'resolved.json'
        self.run_cli('base', source, edits, output, '--manifest', manifest)
        self.assertEqual(output.read_text(), '>mt\nACGT\n')
        self.assertEqual(json.loads(manifest.read_text())['changes'][0]['old'], 'N')

    def test_base_edit_refuses_wrong_reference_base(self):
        source = self.root / 'source.fa'; source.write_text('>mt\nACGT\n')
        edits = self.json_file('edits.json', {'changes': [
            {'sequence_id': 'mt', 'position': 2, 'old': 'T', 'new': 'A', 'evidence': 'E1/x', 'rationale': 'why'}]})
        output, manifest = self.root / 'candidate.fa', self.root / 'manifest.json'
        self.run_cli('base', source, edits, output, '--manifest', manifest, code=1)
        self.assertFalse(output.exists()); self.assertFalse(manifest.exists())

    def test_base_edit_rejects_duplicate_coordinate(self):
        source = self.root / 'source.fa'; source.write_text('>mt\nACGT\n')
        item = {'sequence_id': 'mt', 'position': 1, 'old': 'A', 'new': 'G', 'evidence': 'E1/x', 'rationale': 'why'}
        edits = self.json_file('edits.json', {'changes': [item, item]})
        self.run_cli('base', source, edits, self.root/'candidate.fa', '--manifest', self.root/'m.json', code=1)

    def test_base_edit_rejects_non_string_allele_without_traceback(self):
        source = self.root/'source.fa'; source.write_text('>mt\nACGT\n')
        edits = self.json_file('edits.json', {'changes': [
            {'sequence_id':'mt','position':1,'old':1,'new':'G','evidence':'E1/x','rationale':'why'}]})
        result = self.run_cli('base', source, edits, self.root/'candidate.fa', '--manifest', self.root/'m.json', code=1)
        self.assertNotIn('Traceback', result.stderr)

    def test_candidate_outputs_are_never_silently_overwritten(self):
        source = self.root/'source.fa'; source.write_text('>mt\nACGT\n')
        edits = self.json_file('edits.json', {'changes': [
            {'sequence_id':'mt','position':1,'old':'A','new':'G','evidence':'E1/x','rationale':'why'}]})
        output, manifest = self.root/'candidate.fa', self.root/'m.json'
        output.write_text('keep this candidate\n'); manifest.write_text('keep this manifest\n')
        self.run_cli('base', source, edits, output, '--manifest', manifest, code=1)
        self.assertEqual(output.read_text(), 'keep this candidate\n')
        self.assertEqual(manifest.read_text(), 'keep this manifest\n')

    def test_base_edit_refuses_output_alias(self):
        source = self.root / 'source.fa'; source.write_text('>mt\nACGT\n')
        edits = self.json_file('edits.json', {'changes': [
            {'sequence_id': 'mt', 'position': 1, 'old': 'A', 'new': 'G', 'evidence': 'E1/x', 'rationale': 'why'}]})
        self.run_cli('base', source, edits, source, '--manifest', self.root/'m.json', code=1)
        self.assertEqual(source.read_text(), '>mt\nACGT\n')

    def test_cds_start_shift_requires_explicit_current_codon_start(self):
        record = SeqRecord(Seq('AATGAAATAA'), id='mt', name='mt', description='synthetic')
        record.annotations['molecule_type'] = 'DNA'
        record.features = [SeqFeature(FeatureLocation(0, 10, strand=1), type='CDS',
                                      qualifiers={'locus_tag': ['cox1'], 'codon_start': ['2']})]
        source = self.root/'shift.gb'; SeqIO.write(record, source, 'genbank')
        edit = {'selector': {'qualifier':'locus_tag','value':'cox1'},
                'old': {'start':1,'end':10,'strand':1}, 'new': {'start':2,'end':10,'strand':1},
                'old_codon_start':1,'new_codon_start':1,'evidence':'E1/ref','rationale':'test'}
        edits = self.json_file('shift.json', {'changes':[edit]})
        self.run_cli('cds-boundary', source, edits, self.root/'bad.gb', '--manifest', self.root/'bad.json', '--table', 1, code=1)
        edit['old_codon_start'] = 2; edit['new_codon_start'] = 1
        edits = self.json_file('shift.json', {'changes':[edit]})
        output=self.root/'shift-candidate.gb'; manifest=self.root/'shift-manifest.json'
        self.run_cli('cds-boundary', source, edits, output, '--manifest', manifest, '--table', 1)
        result=SeqIO.read(output,'genbank')
        cds=next(f for f in result.features if f.type=='CDS')
        self.assertEqual(cds.qualifiers['codon_start'], ['1'])
        self.assertEqual(json.loads(manifest.read_text())['changes'][0]['old_codon_start'], 2)

    def test_cds_boundary_change_requires_table_and_records_translation(self):
        record = SeqRecord(Seq('ATGAAATAA'), id='mt', name='mt', description='synthetic')
        record.annotations['molecule_type'] = 'DNA'
        record.features = [SeqFeature(FeatureLocation(0, 9, strand=1), type='CDS',
                                      qualifiers={'locus_tag': ['cox1'], 'codon_start': ['1']})]
        source = self.root / 'source.gb'; SeqIO.write(record, source, 'genbank')
        edits = self.json_file('edits.json', {'changes': [{
            'selector': {'qualifier': 'locus_tag', 'value': 'cox1'},
            'old': {'start': 1, 'end': 9, 'strand': 1},
            'new': {'start': 1, 'end': 6, 'strand': 1},
            'old_codon_start': 1, 'new_codon_start': 1,
            'evidence': 'E1/ref', 'rationale': 'synthetic boundary example'}]})
        self.run_cli('cds-boundary', source, edits, self.root/'candidate.gb', '--manifest', self.root/'m.json', code=1)
        self.run_cli('cds-boundary', source, edits, self.root/'candidate.gb', '--manifest', self.root/'m.json', '--table', 1)
        manifest = json.loads((self.root/'m.json').read_text())
        self.assertEqual(manifest['changes'][0]['translation_before'], 'MK*')
        self.assertEqual(manifest['changes'][0]['translation_after'], 'MK')

    def test_qualifier_candidate_corrects_source_and_removes_stale_cds_fields(self):
        record = SeqRecord(Seq('ATGAAATAA'), id='mt', name='mt', description='synthetic')
        record.annotations.update(molecule_type='DNA', organism='.', topology='linear')
        record.features = [
            SeqFeature(FeatureLocation(0, 9), type='source', qualifiers={
                'organism': ['Donor species'], 'db_xref': ['taxon:1'], 'collection_date': ['old date']}),
            SeqFeature(FeatureLocation(0, 9, strand=1), type='CDS', qualifiers={
                'gene': ['cox1'], 'protein_id': ['DONOR.1'], 'transl_except': ['(pos:9,aa:TERM)']}),
        ]
        source = self.root/'source.gb'; SeqIO.write(record, source, 'genbank')
        edits = self.json_file('qualifier-edits.json', {'changes': [
            {'selector': {'type': 'source'}, 'qualifier': 'organism', 'old': ['Donor species'],
             'new': ['Sample species'], 'evidence': 'user sample declaration', 'rationale': 'correct source'},
            {'selector': {'type': 'source'}, 'qualifier': 'db_xref', 'old': ['taxon:1'],
             'new': ['taxon:2'], 'evidence': 'taxon registry', 'rationale': 'match sample'},
            {'selector': {'type': 'source'}, 'qualifier': 'collection_date', 'old': ['old date'],
             'new': [], 'evidence': 'donor record', 'rationale': 'remove donor metadata'},
            {'selector': {'type': 'CDS', 'qualifier': 'gene', 'value': 'cox1'},
             'qualifier': 'protein_id', 'old': ['DONOR.1'], 'new': [],
             'evidence': 'donor record', 'rationale': 'remove donor accession'},
            {'selector': {'type': 'CDS', 'qualifier': 'gene', 'value': 'cox1'},
             'qualifier': 'transl_except', 'old': ['(pos:9,aa:TERM)'], 'new': [],
             'evidence': 'translation check', 'rationale': 'no internal stop'},
        ]})
        output=self.root/'candidate.gb'; manifest=self.root/'manifest.json'
        self.run_cli('qualifier', source, edits, output, '--manifest', manifest)
        candidate=SeqIO.read(output, 'genbank')
        self.assertEqual(str(candidate.seq), str(record.seq))
        self.assertEqual(candidate.annotations['organism'], 'Sample species')
        self.assertEqual(candidate.annotations['source'], 'Sample species')
        self.assertEqual(candidate.features[0].qualifiers['organism'], ['Sample species'])
        self.assertEqual(candidate.features[0].qualifiers['db_xref'], ['taxon:2'])
        self.assertNotIn('collection_date', candidate.features[0].qualifiers)
        self.assertNotIn('protein_id', candidate.features[1].qualifiers)
        self.assertNotIn('transl_except', candidate.features[1].qualifiers)
        self.assertEqual(len(json.loads(manifest.read_text())['changes']), 5)
        self.assertEqual(SeqIO.read(source, 'genbank').features[0].qualifiers['organism'], ['Donor species'])

    def test_qualifier_candidate_rejects_stale_value_and_nonunique_selector(self):
        record = SeqRecord(Seq('ATGAAATAA'), id='mt', name='mt', description='synthetic')
        record.annotations['molecule_type'] = 'DNA'
        record.features = [SeqFeature(FeatureLocation(0, 9), type='source',
                                      qualifiers={'organism': ['Donor species']})]
        source=self.root/'source.gb'; SeqIO.write(record, source, 'genbank')
        edit={'selector': {'type':'source'}, 'qualifier':'organism', 'old':['Wrong species'],
              'new':['Sample species'], 'evidence':'E1', 'rationale':'correct source'}
        edits=self.json_file('edits.json', {'changes':[edit]})
        output=self.root/'candidate.gb'; manifest=self.root/'manifest.json'
        self.run_cli('qualifier', source, edits, output, '--manifest', manifest, code=1)
        self.assertFalse(output.exists()); self.assertFalse(manifest.exists())
        record.features.append(SeqFeature(FeatureLocation(0, 9), type='source',
                                          qualifiers={'organism': ['Donor species']}))
        SeqIO.write(record, source, 'genbank')
        edit['old']=['Donor species']
        self.json_file('edits.json', {'changes':[edit]})
        self.run_cli('qualifier', source, edits, output, '--manifest', manifest, code=1)
        self.assertFalse(output.exists()); self.assertFalse(manifest.exists())


if __name__ == '__main__':
    unittest.main()

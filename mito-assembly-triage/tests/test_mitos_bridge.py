import importlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import contextlib
import io
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts'
sys.path.insert(0, str(SCRIPTS))
with patch.dict(os.environ, {'MITO_ENV_GATE': 'off'}):
    bridge = importlib.import_module('mitos2_to_genbank')
    annot = importlib.import_module('annot_check')


class BridgeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def fixture(self, seq='ATGTGATAA', protein='MW', name='cox1', kind='gene', phase='.'):
        (self.root/'genome.fa').write_text('>mt\n'+seq+'\n')
        (self.root/'result.gff').write_text('mt\tMITOS\t%s\t1\t%d\t.\t+\t%s\tName=%s\n' %
                                           (kind, len(seq), phase, name))
        header = '>mt;1-%d;+;%s\n' % (len(seq), name)
        (self.root/'result.fas').write_text(header+seq+'\n')
        (self.root/'result.faa').write_text(header+protein+'\n' if protein else '')

    def run_bridge(self, table='2', out='out.gb', extra=()):
        args = [str(SCRIPTS/'mitos2_to_genbank.py'), str(self.root/'result.gff'),
                str(self.root/'result.fas'), str(self.root/'result.faa'), str(self.root/out),
                '--genome', str(self.root/'genome.fa')]
        if table is not None:
            args += ['--table', table]
        return subprocess.run([sys.executable, *args, *extra], capture_output=True, text=True,
                              env=dict(os.environ, MITO_ENV_GATE='off'))

    def test_explicit_table_preserved(self):
        self.fixture()
        result = self.run_bridge()
        self.assertEqual(result.returncode, 0, result.stderr)
        cds = next(f for f in SeqIO.read(self.root/'out.gb', 'genbank').features if f.type == 'CDS')
        self.assertEqual(cds.qualifiers['transl_table'], ['2'])
        self.assertEqual(cds.qualifiers['codon_start'], ['1'])

    def test_caller_supplied_organism_and_taxid_are_explicit(self):
        self.fixture()
        result = self.run_bridge(extra=('--organism', 'Cemus macaoensis', '--taxid', '871496'))
        self.assertEqual(result.returncode, 0, result.stderr)
        source = next(f for f in SeqIO.read(self.root/'out.gb', 'genbank').features if f.type == 'source')
        self.assertEqual(source.qualifiers['organism'], ['Cemus macaoensis'])
        self.assertEqual(source.qualifiers['db_xref'], ['taxon:871496'])

    def test_table_required(self):
        self.fixture()
        self.assertNotEqual(self.run_bridge(table=None).returncode, 0)
        self.assertFalse((self.root/'out.gb').exists())

    def test_multiple_genome_records_rejected(self):
        self.fixture()
        with (self.root/'genome.fa').open('a') as handle:
            handle.write('>second\nACGT\n')
        self.assertNotEqual(self.run_bridge().returncode, 0)
        self.assertFalse((self.root/'out.gb').exists())

    def test_no_input_overwrite(self):
        self.fixture()
        original = (self.root/'genome.fa').read_bytes()
        self.assertNotEqual(self.run_bridge(out='genome.fa').returncode, 0)
        self.assertEqual((self.root/'genome.fa').read_bytes(), original)

    def test_symlink_and_hardlink_cannot_alias_input(self):
        self.fixture()
        original = (self.root/'genome.fa').read_bytes()
        (self.root/'alias.gb').symlink_to(self.root/'genome.fa')
        self.assertNotEqual(self.run_bridge(out='alias.gb').returncode, 0)
        os.link(self.root/'genome.fa', self.root/'hardlink.gb')
        self.assertNotEqual(self.run_bridge(out='hardlink.gb').returncode, 0)
        self.assertEqual((self.root/'genome.fa').read_bytes(), original)

    def test_table_changes_frame_evidence(self):
        self.assertEqual(bridge.best_protein_frame('ATGAGATAA', 'M', 2)[0], None)
        # AGA is internal stop in table 2, Ser in table 5.
        self.assertEqual(bridge.best_protein_frame('ATGAGATAA', 'MS', 5)[0], 0)
        self.assertIsNone(bridge.best_protein_frame('ATGAGATAA', 'MS', 2)[0])

    def test_initiator_methionine_and_terminal_stop(self):
        self.assertEqual(bridge.best_protein_frame('ATTAAATAA', 'MK*', 5)[0], 0)

    def test_weak_frame_is_uncertain(self):
        self.fixture(protein='MQQQQQQ')
        result = self.run_bridge()
        self.assertEqual(result.returncode, 0, result.stderr)
        cds = next(f for f in SeqIO.read(self.root/'out.gb', 'genbank').features if f.type == 'CDS')
        self.assertNotIn('codon_start', cds.qualifiers)
        self.assertIn('codon_start uncertain', ' '.join(cds.qualifiers['note']))
        findings = annot.Findings()
        from Bio.Data import CodonTable
        record = SeqIO.read(self.root/'out.gb', 'genbank')
        with contextlib.redirect_stdout(io.StringIO()):
            annot.cds_findings(findings, record, [cds], CodonTable.unambiguous_dna_by_id[2], set(), set())
        self.assertTrue(any('CODON_START_UNCERTAIN' in item for item in findings.reviews))
        self.assertFalse(findings.errors)

    def test_anticodon_preserved_without_fake_position(self):
        self.fixture(seq='A'*60, protein='', name='trnL(tag)', kind='tRNA')
        result = self.run_bridge()
        self.assertEqual(result.returncode, 0, result.stderr)
        trna = next(f for f in SeqIO.read(self.root/'out.gb', 'genbank').features if f.type == 'tRNA')
        self.assertNotIn('anticodon', trna.qualifiers)
        self.assertEqual(annot.parse_anticodon(trna), 'tag')
        self.assertEqual(annot.resolve_trna_identity(trna), 'trnl1')

    def test_duplicate_fasta_names_rejected(self):
        self.fixture()
        path = self.root/'result.fas'
        path.write_text(path.read_text()*2)
        self.assertNotEqual(self.run_bridge().returncode, 0)

    def test_mitos_headers_with_semicolon_sequence_metadata(self):
        self.fixture()
        header = '>scaffold313;len=15030;topology=linear; 1-9; +; cox1\n'
        (self.root/'result.fas').write_text(header+'ATGTGATAA\n')
        (self.root/'result.faa').write_text(header+'MW\n')
        (self.root/'genome.fa').write_text('>scaffold313;len=15030;topology=linear\nATGTGATAA\n')
        (self.root/'result.gff').write_text(
            'scaffold313;len=15030;topology=linear\tMITOS\tgene\t1\t9\t.\t+\t.\tName=cox1\n')
        result = self.run_bridge()
        self.assertEqual(result.returncode, 0, result.stderr)
        record = SeqIO.read(self.root/'out.gb', 'genbank')
        cds = next(feature for feature in record.features if feature.type == 'CDS')
        self.assertEqual(str(record.seq), 'ATGTGATAA')
        self.assertEqual(cds.qualifiers['gene'], ['cox1'])

    def test_foreign_seqid_rejected(self):
        self.fixture()
        path = self.root/'result.gff'
        path.write_text(path.read_text().replace('mt\t', 'other\t'))
        self.assertNotEqual(self.run_bridge().returncode, 0)

    def test_gene_and_cds_parent_do_not_duplicate(self):
        self.fixture()
        path = self.root/'result.gff'
        path.write_text(path.read_text()+path.read_text().replace('\tgene\t', '\tCDS\t'))
        result = self.run_bridge()
        self.assertEqual(result.returncode, 0, result.stderr)
        features = SeqIO.read(self.root/'out.gb', 'genbank').features
        self.assertEqual(sum(f.type == 'CDS' for f in features), 1)

    def test_phase_and_protein_conflict_rejected(self):
        self.fixture(kind='CDS', phase='1')
        self.assertNotEqual(self.run_bridge().returncode, 0)

    def test_phase_does_not_hide_unmatched_protein(self):
        self.fixture(protein='MQQQQQQ', kind='CDS', phase='0')
        self.assertNotEqual(self.run_bridge().returncode, 0)

    def test_explicit_phase_without_protein_is_source_data(self):
        self.fixture(protein='', kind='CDS', phase='0')
        result = self.run_bridge()
        self.assertEqual(result.returncode, 0, result.stderr)
        cds = next(f for f in SeqIO.read(self.root/'out.gb', 'genbank').features if f.type == 'CDS')
        self.assertEqual(cds.qualifiers['codon_start'], ['1'])
        self.assertIn('GFF CDS phase', ' '.join(cds.qualifiers['note']))


if __name__ == '__main__':
    unittest.main()

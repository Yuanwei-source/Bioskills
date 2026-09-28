"""Task ledger CLI: complete workflows, evidence conflicts and immutable history."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / 'scripts/case.py'


class CaseFixtures:
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.case = self.root/'task'
        self.sequence = self.root/'input.fa'
        self.sequence.write_text('>mt\nACGT\n')
        self.log = self.root/'check.log'
        self.log.write_text('Observation: boundary candidate; reads unavailable\n')

    def cli(self, command, *args, expected=0):
        result = subprocess.run([sys.executable, str(CLI), command, str(self.case), *map(str, args)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout+result.stderr)
        return result

    def init(self, case_type='abnormal_case'):
        self.cli('init', '--case-id', 'sample-1', '--issue', 'boundary review',
                 '--input', 'assembly', self.sequence, '--taxon', 'Insecta',
                 '--case-type', case_type)

    def event(self, event_id='E1', status='completed', code=2):
        record = dict(id=event_id, action='annotation check', command=['annot_check.py', 'input.gb', '--table', '5'],
                      working_directory=str(self.root), tool_version='test-1', exit_code=code,
                      execution_status=status, input_ids=['I1'], summary='boundary candidate',
                      artifacts=[dict(id='log', role='log', path=str(self.log))])
        path = self.root/(event_id+'.json')
        path.write_text(json.dumps(record))
        self.cli('event', '--record', path)

    def conclusion(self, **overrides):
        record = dict(id='C1', claim='Boundary remains uncertain', conclusion_type='annotation_quality',
                      status='UNRESOLVED', confidence='moderate', reads_support='NOT_ASSESSED',
                      evidence_for=[dict(event_id='E1', artifact_id='log', observation='Homology supports one boundary')],
                      evidence_against=[], not_tested=['No transcript evidence'], scope='cox1 boundary',
                      limitations=['No physical sequence validation'], rationale='Alternatives remain',
                      next_steps=['Compare conserved protein boundary'], decided_by='test-agent',
                      review=dict(level='ASSIST', status='pending', reviewer=None, note='Needs boundary review'),
                      competitive_reference_ids=[])
        record.update(overrides)
        return record

    def write_conclusion(self, record=None, expected=0, update=False):
        path = self.root/'conclusion.json'
        path.write_text(json.dumps(record or self.conclusion()))
        self.cli('conclusion', '--record', path, *(['--update'] if update else []), expected=expected)

    def state(self):
        return json.loads((self.case/'case.json').read_text())


class CaseWorkflow(CaseFixtures, unittest.TestCase):
    def test_partial_support_survives_report(self):
        self.init(); self.event(); self.write_conclusion()
        self.cli('validate', '--verify-files')
        self.cli('report')
        text = (self.case/'report.md').read_text()
        self.assertIn('Homology supports one boundary', text)
        self.assertIn('No transcript evidence', text)
        self.assertIn('UNRESOLVED', text)
        self.assertNotIn('无证据支持', text)
        self.assertIn('不证明科学结论', text)

    def test_conflicting_evidence_both_rendered(self):
        self.init(); self.event(); self.event('E2')
        self.write_conclusion(self.conclusion(evidence_against=[dict(event_id='E2', artifact_id='log',
                                                                   observation='Alternative boundary remains plausible')]))
        self.cli('report')
        self.assertIn('Alternative boundary remains plausible', (self.case/'report.md').read_text())

    def test_no_change_without_evidence_not_pass(self):
        self.init()
        self.write_conclusion(self.conclusion(status='NO_CHANGE', confidence='not_assessable',
                                              evidence_for=[], rationale='Data unavailable'))
        self.cli('report')
        self.assertIn('Data unavailable', (self.case/'report.md').read_text())
        self.assertNotIn('质量通过', (self.case/'report.md').read_text())

    def test_dangling_event_rejected_without_mutation(self):
        self.init(); before = (self.case/'case.json').read_bytes()
        self.write_conclusion(expected=1)
        self.assertEqual((self.case/'case.json').read_bytes(), before)

    def test_failed_event_cannot_support_scientific_claim(self):
        self.init(); self.event(status='failed', code=1)
        self.write_conclusion(expected=1)

    def test_failed_event_can_support_tool_failure_claim(self):
        self.init('tool_failure_case'); self.event(status='failed', code=1)
        self.write_conclusion(self.conclusion(conclusion_type='tool_execution'))

    def test_duplicate_event_and_conclusion_rejected(self):
        self.init(); self.event(); self.write_conclusion()
        before = (self.case/'case.json').read_bytes()
        self.cli('event', '--record', self.root/'E1.json', expected=1)
        self.write_conclusion(expected=1)
        self.assertEqual((self.case/'case.json').read_bytes(), before)

    def test_update_preserves_old_claim(self):
        self.init(); self.event(); self.write_conclusion()
        self.write_conclusion(self.conclusion(claim='New claim'), update=True)
        state = self.state()
        self.assertEqual(state['conclusions'][0]['claim'], 'New claim')
        self.assertEqual(state['conclusion_history'][0]['conclusion']['claim'], 'Boundary remains uncertain')

    def test_input_drift_blocks_event_and_report(self):
        self.init(); self.event(); self.write_conclusion()
        self.sequence.write_text('>changed\nTTTT\n')
        self.cli('validate', '--verify-files', expected=1)
        self.cli('report', expected=1)
        self.assertFalse((self.case/'report.md').exists())

    def test_artifact_drift_blocks_report(self):
        self.init(); self.event(); self.write_conclusion()
        self.log.write_text('different output')
        self.cli('report', expected=1)

    def test_report_cannot_overwrite_inputs_or_case(self):
        self.init(); self.event(); self.write_conclusion()
        for path in (self.sequence, self.log, self.case/'case.json', self.case/'.case.lock'):
            with self.subTest(path=path):
                before = path.read_bytes()
                self.cli('report', '--output', path, expected=1)
                self.assertEqual(path.read_bytes(), before)

    def test_resolved_needs_support_but_not_reads_for_annotation(self):
        self.init(); self.event()
        self.write_conclusion(self.conclusion(status='RESOLVED', confidence='high', evidence_for=[]), expected=1)
        self.write_conclusion(self.conclusion(status='RESOLVED', confidence='high'))

    def test_discriminating_needs_registered_competitor(self):
        self.init(); self.event()
        self.write_conclusion(self.conclusion(reads_support='READS_DISCRIMINATING'), expected=1)

    def test_normal_case_rejects_unresolved_conclusion(self):
        self.init('normal_validation_case'); self.event()
        self.write_conclusion(expected=1)

    def test_init_never_overwrites_existing_task(self):
        self.init(); before = (self.case/'case.json').read_bytes()
        self.cli('init', '--case-id', 'other', '--issue', 'other', '--input', 'assembly', self.sequence, expected=1)
        self.assertEqual((self.case/'case.json').read_bytes(), before)

    def test_candidate_registration_and_verification_trace_to_candidate_file(self):
        self.init(); self.event(event_id='E1', code=0)
        reads = self.root/'reads.fq'; reads.write_text('@r\nACGG\n+\nIIII\n')
        self.cli('input','--role','reads','--path',reads)
        candidate_file = self.root/'candidate.fa'; candidate_file.write_text('>mt\nACGG\n')
        producer = dict(id='E2', action='apply_candidate', command=['apply_candidate.py'],
                        working_directory=str(self.root), tool_version='test-1', exit_code=0,
                        execution_status='completed', input_ids=['I1'], summary='one explicit edit',
                        artifacts=[dict(id='candidate', role='candidate', path=str(candidate_file))])
        path = self.root/'producer.json'; path.write_text(json.dumps(producer))
        self.cli('event', '--record', path)
        proposal = dict(id='R1', change_type='sequence', base_input_id='I1', producer_event_id='E2',
                        artifact_id='candidate', rationale='synthetic', expected_effect='test ledger', risks=['synthetic'])
        path = self.root/'candidate.json'; path.write_text(json.dumps(proposal))
        self.cli('candidate', '--record', path)
        verifier = dict(id='E3', action='validate_candidate', command=['validator'],
                        working_directory=str(self.root), tool_version='test-1', exit_code=0,
                        execution_status='completed', input_ids=['I2','I3'], summary='validated candidate',
                        artifacts=[dict(id='validation', role='log', path=str(self.log))])
        path = self.root/'verifier.json'; path.write_text(json.dumps(verifier))
        self.cli('event', '--record', path)
        check = dict(candidate_id='R1', event_id='E3', result='pass', criteria='synthetic criteria', notes='synthetic result')
        path = self.root/'check.json'; path.write_text(json.dumps(check))
        self.cli('candidate-check', '--record', path)
        self.cli('validate', '--verify-files'); self.cli('report')
        state = self.state()
        self.assertEqual(state['candidates'][0]['status'], 'verified')
        self.assertEqual(state['candidates'][0]['candidate_input_id'], 'I3')
        self.assertNotIn('candidate_id', state['candidates'][0]['checks'][0])
        self.assertIn('验证事件 E3', (self.case/'report.md').read_text())

    def test_nonzero_candidate_check_cannot_pass(self):
        self.init(); self.event(event_id='E1', code=0)
        reads = self.root/'reads.fq'; reads.write_text('@r\nACGG\n+\nIIII\n')
        self.cli('input','--role','reads','--path',reads)
        candidate_file = self.root/'candidate.fa'; candidate_file.write_text('>mt\nACGG\n')
        producer = dict(id='E2', action='apply_candidate', command=['apply_candidate.py'],
                        working_directory=str(self.root), tool_version='test-1', exit_code=0,
                        execution_status='completed', input_ids=['I1'], summary='one explicit edit',
                        artifacts=[dict(id='candidate', role='candidate', path=str(candidate_file))])
        path = self.root/'producer.json'; path.write_text(json.dumps(producer)); self.cli('event','--record',path)
        proposal = dict(id='R1', change_type='sequence', base_input_id='I1', producer_event_id='E2',
                        artifact_id='candidate', rationale='synthetic', expected_effect='test ledger', risks=[])
        path = self.root/'candidate.json'; path.write_text(json.dumps(proposal)); self.cli('candidate','--record',path)
        verifier = dict(id='E3', action='validate_candidate', command=['validator'],
                        working_directory=str(self.root), tool_version='test-1', exit_code=2,
                        execution_status='completed', input_ids=['I2','I3'], summary='review required',
                        artifacts=[dict(id='validation', role='log', path=str(self.log))])
        path = self.root/'verifier.json'; path.write_text(json.dumps(verifier)); self.cli('event','--record',path)
        check = dict(candidate_id='R1', event_id='E3', result='pass', criteria='synthetic criteria', notes='must reject')
        path = self.root/'check.json'; path.write_text(json.dumps(check))
        self.cli('candidate-check','--record',path,expected=1)
        verifier['id'] = 'E4'; verifier['exit_code'] = 0; verifier['input_ids'] = ['I3']
        path = self.root/'verifier-no-reads.json'; path.write_text(json.dumps(verifier)); self.cli('event','--record',path)
        check.update(event_id='E4', notes='must reject without reads')
        path = self.root/'check-no-reads.json'; path.write_text(json.dumps(check))
        self.cli('candidate-check','--record',path,expected=1)
        self.assertEqual(self.state()['candidates'][0]['status'], 'proposed')


class CaseIntegrity(CaseFixtures, unittest.TestCase):
    def test_candidate_status_requires_every_check_to_pass(self):
        sys.path.insert(0, str(ROOT/'scripts'))
        from _case_store import candidate_status
        self.assertEqual(candidate_status([]), 'proposed')
        self.assertEqual(candidate_status(['pass']), 'verified')
        self.assertEqual(candidate_status(['pass', 'inconclusive']), 'inconclusive')
        self.assertEqual(candidate_status(['pass', 'fail']), 'inconclusive')
        self.assertEqual(candidate_status(['fail', 'inconclusive']), 'inconclusive')
        self.assertEqual(candidate_status(['fail', 'fail']), 'rejected')

    def test_missing_dependency_leaves_no_task(self):
        result = subprocess.run([sys.executable, '-S', str(CLI), 'init', str(self.case),
            '--case-id', 'x', '--issue', 'x', '--input', 'assembly', str(self.sequence)],
            capture_output=True, text=True, env=dict(os.environ, MITO_ENV_GATE='on'))
        self.assertEqual(result.returncode, 3, result.stderr)
        self.assertFalse(self.case.exists())

    def test_unknown_artifact_rejected(self):
        self.init(); self.event()
        record = self.conclusion()
        record['evidence_for'][0]['artifact_id'] = 'absent'
        self.write_conclusion(record, expected=1)

    def test_stale_revision_cannot_overwrite_new_conclusion(self):
        self.init(); self.event(); self.write_conclusion()
        before = (self.case/'case.json').read_bytes()
        self.cli('conclusion', '--record', self.root/'conclusion.json', '--update',
                 '--expected-revision', '1', expected=1)
        self.assertEqual((self.case/'case.json').read_bytes(), before)

    def test_simultaneous_events_not_lost(self):
        self.init(); self.event()
        record = json.loads((self.root/'E1.json').read_text())
        processes = []
        for identifier in ('E2', 'E3'):
            record['id'] = identifier
            path = self.root/(identifier+'.json')
            path.write_text(json.dumps(record))
            processes.append(subprocess.Popen([sys.executable, str(CLI), 'event', str(self.case),
                '--record', str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
        for process in processes:
            out, err = process.communicate(timeout=30)
            self.assertEqual(process.returncode, 0, out+err)
        self.assertEqual({e['id'] for e in self.state()['events']}, {'E1','E2','E3'})

    def test_reads_support_requires_reads_input(self):
        self.init(); self.event()
        self.write_conclusion(self.conclusion(reads_support='READS_CONSISTENT'), expected=1)

    def test_discriminating_with_reads_and_competitor(self):
        self.init()
        reads, competitor = self.root/'reads.fq', self.root/'nuclear.fa'
        reads.write_text('@r\nACGT\n+\nIIII\n'); competitor.write_text('>nuclear\nTGCA\n')
        self.cli('input', '--role', 'reads', '--path', reads)
        self.cli('input', '--role', 'competitive_reference', '--path', competitor)
        self.event()
        record = json.loads((self.root/'E1.json').read_text())
        record['id'] = 'E2'; record['input_ids'] = ['I1','I2','I3']
        path = self.root/'E2.json'; path.write_text(json.dumps(record))
        self.cli('event', '--record', path)
        conclusion = self.conclusion(reads_support='READS_DISCRIMINATING', competitive_reference_ids=['I3'])
        conclusion['evidence_for'][0]['event_id'] = 'E2'
        self.write_conclusion(conclusion)
        self.cli('report')

    def test_schema_v2_not_silently_migrated(self):
        self.init(); path = self.case/'case.json'
        state = self.state(); state['schema_version'] = '2.0'; path.write_text(json.dumps(state))
        before = path.read_bytes()
        self.cli('validate', expected=1)
        self.assertEqual(path.read_bytes(), before)

    def test_missing_file_blocks_report_without_touching_case(self):
        self.init(); self.event(); self.write_conclusion()
        before = (self.case/'case.json').read_bytes()
        self.log.unlink()
        self.cli('report', expected=1)
        self.assertEqual((self.case/'case.json').read_bytes(), before)

    def test_report_cannot_alias_input(self):
        self.init(); self.event(); self.write_conclusion()
        for name, symlink in [('alias.md',True),('hardlink.md',False)]:
            target = self.root/name
            if symlink: target.symlink_to(self.sequence)
            else: os.link(self.sequence,target)
            self.cli('report', '--output', target, expected=1)
        self.assertEqual(self.sequence.read_text(), '>mt\nACGT\n')

    def test_broken_json_and_unreviewed_identity_fail_cleanly(self):
        self.init(); self.event()
        path = self.root/'bad.json'; path.write_text('{"id":"C1","id":"C2"}')
        result = self.cli('conclusion','--record',path,expected=1)
        self.assertNotIn('Traceback', result.stderr)
        record = self.conclusion(review=dict(level='EXPERT',status='reviewed',reviewer=None,note='claimed'))
        self.write_conclusion(record, expected=1)

    def test_atomic_write_failure_preserves_case(self):
        sys.path.insert(0, str(ROOT/'scripts'))
        from _case_store import CaseStore
        from unittest.mock import patch
        self.init(); before = (self.case/'case.json').read_bytes()
        with patch('_evidence_io.os.replace', side_effect=OSError('simulated full disk')):
            with self.assertRaises(OSError):
                CaseStore(self.case).mutate('input', dict(role='other', path=str(self.log)))
        self.assertEqual((self.case/'case.json').read_bytes(), before)
        self.assertEqual(sorted(p.name for p in self.case.iterdir()), ['.case.lock','case.json'])


if __name__ == '__main__':
    unittest.main()

"""Local lesson memory: authorization, evidence links, review, search and revisions."""
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / 'scripts/memory.py'
CASE_CLI = ROOT / 'scripts/case.py'


class MemoryFixtures:
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.library = self.root / 'knowledge'
        self.sequence = self.root / 'input.fa'
        self.sequence.write_text('>mt\nACGT\n')
        self.log = self.root / 'check.log'
        self.log.write_text('synthetic annotation observation\n')

    def cli(self, action, *args, expected=0, library=None, python=sys.executable, env=None):
        result = subprocess.run([python, str(CLI), '--library', str(library or self.library),
                                 action, *map(str, args)], capture_output=True, text=True, env=env)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def case_cli(self, case, action, *args, expected=0):
        result = subprocess.run([sys.executable, str(CASE_CLI), action, str(case), *map(str, args)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def make_case(self, name='sample-1', status='RESOLVED', case_type='abnormal_case', directory=None):
        case = Path(directory) if directory is not None else self.root / name
        self.case_cli(case, 'init', '--case-id', name, '--issue', 'memory fixture',
                      '--input', 'assembly', self.sequence, '--taxon', 'Insecta',
                      '--case-type', case_type)
        failed = case_type == 'tool_failure_case'
        event = dict(id='E1', action='annotation check', command=['annot_check.py', '--table', '5'],
                     working_directory=str(self.root), tool_version='test-1',
                     exit_code=1 if failed else (0 if status == 'RESOLVED' else 2),
                     execution_status='failed' if failed else 'completed', input_ids=['I1'],
                     summary='synthetic observation',
                     artifacts=[dict(id='log', role='log', path=str(self.log))])
        event_path = self.root / (name + '-event.json')
        event_path.write_text(json.dumps(event))
        self.case_cli(case, 'event', '--record', event_path)
        conclusion = dict(id='C1', claim='Synthetic claim',
            conclusion_type='tool_execution' if failed else 'annotation_quality', status=status,
            confidence='high' if status == 'RESOLVED' else 'moderate', reads_support='NOT_ASSESSED',
            evidence_for=[dict(event_id='E1', artifact_id='log', observation='Synthetic evidence')],
            evidence_against=[], not_tested=['Synthetic untested item'], scope='nad6',
            limitations=['Synthetic limitation'], rationale='Synthetic rationale', next_steps=['Recheck'],
            decided_by='test-agent',
            review=dict(level='ASSIST', status='pending', reviewer=None, note='Synthetic review'),
            competitive_reference_ids=[])
        conclusion_path = self.root / (name + '-conclusion.json')
        conclusion_path.write_text(json.dumps(conclusion))
        self.case_cli(case, 'conclusion', '--record', conclusion_path)
        return case

    def draft(self, **overrides):
        record = dict(title='nad6 断裂后的竞争回贴',
                      claim='单一注释断裂不能证明样本碱基错误',
                      lesson_domain='annotation',
                      diagnostic_clues=['nad6 断裂', '局部覆盖下降'],
                      suggested_next_test='检查组装图并对竞争候选回贴',
                      applicable_when=dict(taxa=['Insecta'], data_types=['illumina_paired_end'],
                                           conditions=['记录组装器和版本']),
                      not_applicable_when=['没有当前样本证据时直接修改碱基'],
                      counterexample_case_ids=[], conflicts=[])
        record.update(overrides)
        return record

    def write_draft(self, record=None, name='lesson.json'):
        path = self.root / name
        path.write_text(json.dumps(record or self.draft(), ensure_ascii=False))
        return path

    def propose(self, case, record=None, expected=0, authorize=True, allow=False, library=None):
        path = self.write_draft(record)
        args = ['--case', case, '--record', path]
        if authorize:
            args.append('--authorize-local')
        if allow:
            args.append('--allow-unresolved')
        return self.cli('propose', *args, expected=expected, library=library)

    def verify(self, lesson_id='L000001', expected=0, library=None):
        return self.cli('review', lesson_id, '--status', 'verified', '--reviewer', 'human',
                        '--reason', 'checked synthetic evidence', '--authorize-local',
                        expected=expected, library=library)


class MemoryWorkflow(MemoryFixtures, unittest.TestCase):
    def test_propose_requires_authorization_and_preserves_case(self):
        case = self.make_case()
        before = (case / 'case.json').read_bytes()
        self.propose(case, authorize=False, expected=1)
        self.assertFalse(self.library.exists())
        self.propose(case)
        self.assertEqual((case / 'case.json').read_bytes(), before)
        document = json.loads((self.library / 'lessons/L000001.json').read_text())
        revision = document['revisions'][0]
        self.assertEqual(revision['validation_status'], 'candidate')
        self.assertEqual(revision['supporting_case_ids'], ['sample-1'])
        self.assertEqual(revision['sources'][0]['path'], str(case / 'case.json'))
        self.assertEqual(revision['generalization_scope'], 'single_case')
        self.assertEqual(revision['transferability'], 'none')
        self.assertTrue((self.library / 'index.sqlite3').is_file())
        for path in (self.library, self.library / 'lessons'):
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o700)
        for path in (self.library / 'lessons/L000001.json', self.library / 'index.sqlite3'):
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o600)

    def test_unresolved_and_duplicate_are_fail_closed(self):
        case = self.make_case(status='UNRESOLVED')
        self.propose(case, expected=1)
        self.assertFalse(self.library.exists())
        self.propose(case, allow=True)
        self.propose(case, allow=True, expected=1)
        self.assertEqual(len(list((self.library / 'lessons').glob('*.json'))), 1)
        self.verify(expected=1)
        self.cli('review', 'L000001', '--status', 'rejected', '--reviewer', 'human',
                 '--reason', 'unresolved source', '--authorize-local')

    def test_review_search_show_and_exact_scope_filters(self):
        case = self.make_case()
        self.propose(case)
        self.assertNotIn('L000001', self.cli('search', '--query', 'nad6').stdout)
        self.assertIn('L000001', self.cli('search', '--status', 'candidate', '--query', 'nad6').stdout)
        self.verify()
        self.assertIn('L000001', self.cli('search').stdout)
        self.assertNotIn('L000001', self.cli('search', '--query', '!!!').stdout)
        result = self.cli('search', '--query', 'nad6 断裂', '--taxon', 'insecta',
                          '--data-type', 'illumina_paired_end', '--domain', 'annotation')
        self.assertIn('L000001', result.stdout)
        self.assertIn('不替代当前样本检查', result.stdout)
        self.assertNotIn('L000001', self.cli('search', '--query', 'nad6', '--taxon', 'Diptera').stdout)
        self.assertNotIn('L000001', self.cli('search', '--query', 'nad6', '--data-type', 'assembly_fasta').stdout)
        self.cli('show', 'L000001', '--revision', '0', expected=1)
        shown = json.loads(self.cli('show', 'L000001').stdout)
        self.assertEqual(shown['validation_status'], 'verified')
        self.assertEqual(shown['review_history'][0]['reviewer'], 'human')

    def test_tool_failure_cannot_become_biology_lesson(self):
        case = self.make_case(case_type='tool_failure_case', status='UNRESOLVED')
        self.propose(case, self.draft(lesson_domain='biology'), allow=True, expected=1)
        self.assertFalse(self.library.exists())
        self.propose(case, self.draft(lesson_domain='tool'), allow=True)
        shown = json.loads(self.cli('show', 'L000001').stdout)
        self.assertEqual(shown['source_case_type'], 'tool_failure_case')
        self.assertEqual(shown['lesson_domain'], 'tool')

    def test_supersede_preserves_old_revision_and_resets_review(self):
        case = self.make_case()
        self.propose(case)
        self.verify()
        changed = self.draft(claim='修订后的经验仍须重新审核')
        path = self.write_draft(changed)
        self.cli('supersede', 'L000001', '--case', case, '--record', path, '--authorize-local')
        old = json.loads(self.cli('show', 'L000001', '--revision', '1').stdout)
        current = json.loads(self.cli('show', 'L000001').stdout)
        self.assertEqual(old['validation_status'], 'superseded')
        self.assertEqual(current['revision'], 2)
        self.assertEqual(current['validation_status'], 'candidate')
        self.assertEqual(current['supersedes'], {'lesson_id': 'L000001', 'revision': 1})
        self.assertNotIn('L000001', self.cli('search', '--query', 'nad6').stdout)
        self.cli('supersede', 'L000001', '--case', case, '--record', path,
                 '--authorize-local', expected=1)

    def test_source_case_drift_blocks_verification(self):
        case = self.make_case()
        draft_path = self.write_draft()
        self.cli('propose', '--case', case, '--record', draft_path, '--authorize-local')
        extra = self.root / 'extra.txt'
        extra.write_text('additional registered input\n')
        self.case_cli(case, 'input', '--role', 'other', '--path', extra)
        self.verify(expected=1)
        shown = json.loads(self.cli('show', 'L000001').stdout)
        self.assertEqual(shown['validation_status'], 'candidate')
        self.cli('supersede', 'L000001', '--case', case, '--record', draft_path, '--authorize-local')
        refreshed = json.loads(self.cli('show', 'L000001').stdout)
        self.assertEqual(refreshed['revision'], 2)
        self.assertNotEqual(refreshed['content_hash'], shown['content_hash'])
        self.verify()

    def test_stale_or_missing_index_requires_authorized_rebuild(self):
        case = self.make_case()
        self.propose(case)
        self.verify()
        source = self.library / 'lessons/L000001.json'
        source.write_text(source.read_text() + '\n')
        self.cli('search', '--query', 'nad6', expected=1)
        self.cli('rebuild', expected=1)
        self.cli('rebuild', '--authorize-local')
        self.assertIn('L000001', self.cli('search', '--query', 'nad6').stdout)
        (self.library / 'index.sqlite3').unlink()
        self.cli('search', '--query', 'nad6', expected=1)
        self.assertFalse((self.library / 'index.sqlite3').exists())
        self.cli('rebuild', '--authorize-local')
        self.assertIn('L000001', self.cli('search', '--query', 'nad6').stdout)

    def test_private_or_unstructured_drafts_are_rejected(self):
        case = self.make_case()
        bad_values = [
            self.draft(title='见 /home/dell/private/sample.fa'),
            self.draft(claim='联系 researcher@example.org'),
            self.draft(suggested_next_test='设置 api_key=secret-value'),
            self.draft(diagnostic_clues=['>mt\nACGTACGTACGTACGTACGT']),
            self.draft(claim='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'),
            self.draft(validation_status='verified'),
        ]
        for index, record in enumerate(bad_values):
            with self.subTest(index=index):
                library = self.root / ('private-%d' % index)
                self.propose(case, record, expected=1, library=library)
                self.assertFalse(library.exists())

    def test_conflicting_lessons_remain_separate(self):
        first = self.make_case('sample-1')
        second = self.make_case('sample-2')
        self.propose(first)
        self.propose(second, self.draft(title='nad6 断裂的反例', claim='当前候选不能排除 NUMT',
                                        conflicts=['L000001']))
        self.verify('L000001')
        self.verify('L000002')
        result = self.cli('search', '--query', 'nad6', '--limit', '10').stdout
        self.assertIn('L000001', result)
        self.assertIn('L000002', result)
        shown = json.loads(self.cli('show', 'L000002').stdout)
        self.assertEqual(shown['conflicts'], ['L000001'])

    def test_same_case_id_in_different_directories_cannot_share_revision_chain(self):
        first = self.make_case('shared-id')
        second = self.make_case('shared-id', directory=self.root / 'independent-case')
        self.propose(first)
        self.propose(second)
        self.assertTrue((self.library / 'lessons/L000002.json').is_file())
        changed = self.write_draft(self.draft(title='不同案例的修订不能串接'))
        self.cli('supersede', 'L000001', '--case', second, '--record', changed,
                 '--authorize-local', expected=1)
        self.assertEqual(json.loads(self.cli('show', 'L000001').stdout)['revision'], 1)

    def test_review_waits_for_case_lock_and_rechecks_hash(self):
        sys.path.insert(0, str(ROOT / 'scripts'))
        from _case_store import CaseStore, now, snapshot
        case = self.make_case()
        self.propose(case)
        store = CaseStore(case)
        command = [sys.executable, str(CLI), '--library', str(self.library), 'review', 'L000001',
                   '--status', 'verified', '--reviewer', 'human', '--reason', 'concurrent check',
                   '--authorize-local']
        with store.locked():
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            self.addCleanup(lambda: process.poll() is None and (process.kill(), process.communicate()))
            try:
                early = process.communicate(timeout=0.2)
            except subprocess.TimeoutExpired:
                early = None
            extra = self.root / 'concurrent.txt'
            extra.write_text('registered while review waits\n')
            state = store.load()
            state['inputs'].append(snapshot(extra, 'I2', 'other'))
            state['revision'] += 1
            state['updated_at'] = now()
            store.save(state)
        out, err = early if early is not None else process.communicate(timeout=30)
        self.assertEqual(process.returncode, 1, out + err)
        self.assertIn('来源案例自提案后已变化', err)
        self.assertEqual(json.loads(self.cli('show', 'L000001').stdout)['validation_status'], 'candidate')

    def test_normal_validation_case_is_searchable_without_becoming_universal(self):
        case = self.make_case(case_type='normal_validation_case')
        self.propose(case, self.draft(title='正常验证对照', lesson_domain='biology'))
        self.verify()
        shown = json.loads(self.cli('show', 'L000001').stdout)
        self.assertEqual(shown['source_case_type'], 'normal_validation_case')
        self.assertEqual(shown['generalization_scope'], 'single_case')
        self.assertEqual(shown['transferability'], 'none')

    def test_missing_dependency_creates_no_library(self):
        case = self.make_case()
        path = self.write_draft()
        env = dict(os.environ, MITO_ENV_GATE='on')
        result = subprocess.run([sys.executable, '-S', str(CLI), '--library', str(self.library),
                                 'propose', '--case', str(case), '--record', str(path), '--authorize-local'],
                                capture_output=True, text=True, env=env)
        self.assertEqual(result.returncode, 3, result.stdout + result.stderr)
        self.assertFalse(self.library.exists())


if __name__ == '__main__':
    unittest.main()

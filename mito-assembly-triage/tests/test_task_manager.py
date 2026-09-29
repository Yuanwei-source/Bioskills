"""Background task launcher lifecycle tests."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
MANAGER = ROOT/'scripts/task_manager.py'


class BackgroundTaskTests(unittest.TestCase):
    def test_list_reports_running_and_finished_tasks(self):
        with tempfile.TemporaryDirectory() as temp:
            task_root = Path(temp)/'tasks'
            started = subprocess.run([
                sys.executable, str(MANAGER), 'start', '--task-root', str(task_root),
                '--name', 'listed task', '--cwd', str(ROOT), '--', sys.executable, '-c',
                'import time; time.sleep(.15)',
            ], capture_output=True, text=True, check=True)
            task_dir = Path(json.loads(started.stdout)['task_dir'])
            deadline = time.monotonic()+10
            while time.monotonic() < deadline:
                result = subprocess.run([
                    sys.executable, str(MANAGER), 'list', '--task-root', str(task_root),
                ], capture_output=True, text=True, check=True)
                tasks = json.loads(result.stdout)
                if tasks and tasks[0]['status'] == 'succeeded':
                    break
                time.sleep(.05)
            self.assertEqual(tasks[0]['task_dir'], str(task_dir))
            self.assertEqual(tasks[0]['status'], 'succeeded')
            self.assertEqual(tasks[0]['exit_code'], 0)

    def test_list_marks_dead_supervisor_interrupted(self):
        with tempfile.TemporaryDirectory() as temp:
            task_root = Path(temp)/'tasks'; task_root.mkdir()
            task_dir = task_root/'orphan'; task_dir.mkdir()
            (task_dir/'task.json').write_text(json.dumps({
                'format': 'mito-background-task-1', 'name': 'orphan', 'status': 'running',
                'supervisor_pid': 99999999, 'task_dir': str(task_dir),
            }))
            result = subprocess.run([
                sys.executable, str(MANAGER), 'list', '--task-root', str(task_root),
            ], capture_output=True, text=True, check=True)
            tasks = json.loads(result.stdout)
            self.assertEqual(tasks[0]['status'], 'interrupted')
            self.assertIn('退出码', tasks[0]['error'])

    def test_log_tail_handles_carriage_return_progress_and_bounds_large_lines(self):
        with tempfile.TemporaryDirectory() as temp:
            task_root = Path(temp)/'tasks'
            started = subprocess.run([
                sys.executable, str(MANAGER), 'start', '--task-root', str(task_root),
                '--name', 'cr-progress', '--cwd', str(ROOT), '--', sys.executable, '-c',
                'print("old progress\\rnew progress\\rfinal marker")',
            ], capture_output=True, text=True, check=True)
            task_dir = Path(json.loads(started.stdout)['task_dir'])
            deadline = time.monotonic()+10
            while time.monotonic() < deadline:
                state = json.loads(subprocess.check_output([
                    sys.executable, str(MANAGER), 'status', '--task-dir', str(task_dir),
                ], text=True))
                if state['status'] not in ('queued', 'running'):
                    break
                time.sleep(.05)
            log = subprocess.run([
                sys.executable, str(MANAGER), 'log', '--task-dir', str(task_dir), '--lines', '3',
            ], capture_output=True, text=True, check=True)
            self.assertIn('old progress', log.stdout)
            self.assertIn('new progress', log.stdout)
            self.assertIn('final marker', log.stdout)
            with (task_dir/'task.log').open('a', encoding='utf-8') as handle:
                handle.write('x'*2_500_000+'\rend-marker\n')
            log = subprocess.run([
                sys.executable, str(MANAGER), 'log', '--task-dir', str(task_dir), '--lines', '3',
            ], capture_output=True, text=True, check=True)
            self.assertLess(len(log.stdout), 2_100_000)
            self.assertTrue('end-marker' in log.stdout)

    def test_detached_command_records_log_and_final_exit_status(self):
        with tempfile.TemporaryDirectory() as temp:
            task_root = Path(temp)/'intermediate'/'tasks'
            started = subprocess.run([
                sys.executable, str(MANAGER), 'start', '--task-root', str(task_root),
                '--name', 'smoke run', '--cwd', str(ROOT), '--', sys.executable, '-c',
                'import time; print("background-ok", flush=True); time.sleep(.2)',
            ], capture_output=True, text=True, check=True)
            launch = json.loads(started.stdout)
            task_dir = Path(launch['task_dir'])
            deadline = time.monotonic()+10
            while time.monotonic() < deadline:
                status = subprocess.run([
                    sys.executable, str(MANAGER), 'status', '--task-dir', str(task_dir),
                ], capture_output=True, text=True, check=True)
                state = json.loads(status.stdout)
                if state['status'] in ('succeeded', 'failed', 'failed_to_start', 'interrupted'):
                    break
                time.sleep(.05)
            self.assertEqual(state['status'], 'succeeded', state)
            self.assertEqual(state['exit_code'], 0)
            log = subprocess.run([
                sys.executable, str(MANAGER), 'log', '--task-dir', str(task_dir),
            ], capture_output=True, text=True, check=True)
            self.assertIn('background-ok', log.stdout)
            self.assertEqual(state['command'][0], sys.executable)

    def test_nonzero_exit_is_not_reported_as_success(self):
        with tempfile.TemporaryDirectory() as temp:
            task_root = Path(temp)/'tasks'
            started = subprocess.run([
                sys.executable, str(MANAGER), 'start', '--task-root', str(task_root),
                '--name', 'expected failure', '--cwd', str(ROOT), '--', sys.executable,
                '-c', 'raise SystemExit(7)',
            ], capture_output=True, text=True, check=True)
            task_dir = Path(json.loads(started.stdout)['task_dir'])
            deadline = time.monotonic()+10
            while time.monotonic() < deadline:
                output = subprocess.run([
                    sys.executable, str(MANAGER), 'status', '--task-dir', str(task_dir),
                ], capture_output=True, text=True, check=True)
                state = json.loads(output.stdout)
                if state['status'] not in ('queued', 'running'):
                    break
                time.sleep(.05)
            self.assertEqual(state['status'], 'failed', state)
            self.assertEqual(state['exit_code'], 7)


if __name__ == '__main__':
    unittest.main()

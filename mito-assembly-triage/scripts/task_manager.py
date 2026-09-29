#!/usr/bin/env python3
"""Launch and track long-running skill commands independently of the session.

Task state and logs live under a caller-provided intermediate/ directory. The
launcher uses argument vectors (never a shell command string), starts a
detached supervisor, and records the final exit status when the child finishes.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import uuid


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')


def read_state(task_dir):
    return json.loads((task_dir / 'task.json').read_text(encoding='utf-8'))


def write_state(task_dir, state):
    temp = task_dir / 'task.json.tmp'
    temp.write_text(json.dumps(state, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    temp.replace(task_dir / 'task.json')


def pid_alive(pid):
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
        return True
    except (ProcessLookupError, ValueError):
        return False
    except PermissionError:
        return True


def process_group_alive(pgid):
    if not pgid:
        return False
    try:
        os.killpg(int(pgid), 0)
        return True
    except (ProcessLookupError, ValueError):
        return False
    except PermissionError:
        return True


def supervise(task_dir):
    state = read_state(task_dir)
    state.update({'status': 'running', 'started_at': now(), 'supervisor_pid': os.getpid(),
                  'process_group_id': os.getpgrp()})
    write_state(task_dir, state)
    log_path = task_dir / 'task.log'
    try:
        with log_path.open('ab', buffering=0) as log:
            proc = subprocess.Popen(state['command'], cwd=state['cwd'], env=os.environ.copy(),
                                    stdin=subprocess.DEVNULL, stdout=log,
                                    stderr=subprocess.STDOUT, close_fds=True)
            state.update({'child_pid': proc.pid, 'process_group_id': os.getpgid(proc.pid)})
            write_state(task_dir, state)
            exit_code = proc.wait()
        state.update({'status': 'succeeded' if exit_code == 0 else 'failed',
                      'exit_code': exit_code, 'finished_at': now()})
    except OSError as exc:
        with log_path.open('a', encoding='utf-8') as log:
            log.write('启动命令失败: %s\n' % exc)
        state.update({'status': 'failed_to_start', 'exit_code': None,
                      'error': str(exc), 'finished_at': now()})
    write_state(task_dir, state)
    return 0


def tail(path, lines):
    if not path.is_file():
        return ''
    # Some bioinformatics tools refresh one progress line with carriage returns
    # instead of newlines. Bound the read and treat both as line separators.
    with path.open('rb') as handle:
        handle.seek(0, os.SEEK_END)
        size = handle.tell()
        handle.seek(max(0, size - 2 * 1024 * 1024))
        data = handle.read()
    content = data.decode('utf-8', errors='replace').replace('\r\n', '\n').replace('\r', '\n')
    return '\n'.join(content.splitlines()[-lines:]) + ('\n' if content else '')


def refresh_orphaned_state(task_dir, state):
    supervisor_pid = state.get('supervisor_pid')
    if not supervisor_pid:
        try:
            supervisor_pid = int((task_dir / 'supervisor.pid').read_text(encoding='ascii').strip())
        except (OSError, ValueError):
            supervisor_pid = None
    if state.get('status') in ('queued', 'running') and supervisor_pid and not pid_alive(supervisor_pid):
        pgid = state.get('process_group_id')
        if process_group_alive(pgid):
            state.update({'status': 'orphaned_running', 'observed_at': now(),
                          'error': '后台管理进程已退出，但任务进程组仍活动；请勿重启，先检查进程与日志'})
        else:
            state.update({'status': 'supervisor_lost', 'observed_at': now(),
                          'error': '后台管理进程已退出且未观察到任务进程组；最终退出码未知'})
        write_state(task_dir, state)
    return state


def main(argv=None):
    parser = argparse.ArgumentParser(description='后台启动并跟踪长时间运行的 skill 任务')
    commands = parser.add_subparsers(dest='action', required=True)
    start = commands.add_parser('start', help='在新任务目录中启动命令')
    start.add_argument('--task-root', required=True, help='任务存储目录，通常是样本 intermediate/tasks')
    start.add_argument('--name', required=True, help='任务简短名称')
    start.add_argument('--cwd', default=str(Path.cwd()), help='命令工作目录')
    start.add_argument('command', nargs=argparse.REMAINDER, help='命令及参数，放在 -- 后')
    status = commands.add_parser('status', help='查看任务状态')
    status.add_argument('--task-dir', required=True)
    listing = commands.add_parser('list', help='列出任务根目录下的任务及当前状态')
    listing.add_argument('--task-root', required=True, help='包含任务子目录的目录')
    log = commands.add_parser('log', help='显示任务日志末尾')
    log.add_argument('--task-dir', required=True)
    log.add_argument('--lines', type=int, default=60)
    internal = commands.add_parser('_worker', help=argparse.SUPPRESS)
    internal.add_argument('--task-dir', required=True)
    args = parser.parse_args(argv)

    if args.action == '_worker':
        return supervise(Path(args.task_dir).resolve())
    if args.action == 'start':
        command = args.command[1:] if args.command[:1] == ['--'] else args.command
        if not command:
            parser.error('start 必须在 -- 后提供要运行的命令')
        task_root = Path(args.task_root).resolve()
        cwd = Path(args.cwd).resolve()
        task_root.mkdir(parents=True, exist_ok=True)
        safe_name = re.sub(r'[^A-Za-z0-9_.-]+', '_', args.name).strip('._') or 'task'
        task_dir = task_root / ('%s-%s' % (safe_name, uuid.uuid4().hex[:8]))
        task_dir.mkdir()
        state = {'format': 'mito-background-task-1', 'name': args.name,
                 'created_at': now(), 'status': 'queued', 'exit_code': None,
                 'command': command, 'cwd': str(cwd), 'task_dir': str(task_dir),
                 'log': str(task_dir / 'task.log')}
        write_state(task_dir, state)
        proc = subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), '_worker',
             '--task-dir', str(task_dir)], cwd=task_dir,
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, close_fds=True, start_new_session=True,
            env=os.environ.copy())
        (task_dir / 'supervisor.pid').write_text(str(proc.pid) + '\n', encoding='ascii')
        print(json.dumps({'task_dir': str(task_dir), 'task_file': str(task_dir / 'task.json'),
                          'log': str(task_dir / 'task.log'), 'supervisor_pid': proc.pid,
                          'status': 'queued'}, ensure_ascii=False, indent=2))
        return 0
    if args.action == 'list':
        task_root = Path(args.task_root).resolve()
        tasks = []
        if task_root.is_dir():
            for task_dir in sorted((p for p in task_root.iterdir() if p.is_dir()), key=lambda p: p.name):
                try:
                    state = refresh_orphaned_state(task_dir, read_state(task_dir))
                except (OSError, ValueError):
                    continue
                tasks.append({key: state.get(key) for key in
                              ('name', 'status', 'exit_code', 'created_at', 'started_at',
                               'finished_at', 'task_dir', 'log', 'supervisor_pid', 'child_pid',
                               'process_group_id', 'error')
                              if key in state})
        print(json.dumps(tasks, ensure_ascii=False, indent=2))
        return 0
    task_dir = Path(args.task_dir).resolve()
    try:
        state = read_state(task_dir)
    except (OSError, ValueError) as exc:
        print('任务状态无法读取: %s' % exc, file=sys.stderr)
        return 1
    if args.action == 'status':
        state = refresh_orphaned_state(task_dir, state)
        print(json.dumps(state, ensure_ascii=False, indent=2))
        return 0
    if args.lines < 1:
        parser.error('--lines 必须为正整数')
    print(tail(task_dir / 'task.log', args.lines), end='')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

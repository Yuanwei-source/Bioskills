#!/usr/bin/env python3
"""
进化引擎: 管理 skill 的知识层 (knowledge/)

命令:
  search --query "关键词"           # 检索历史案例 (按信号/物种/内容)
  add-case --sample X --species Y   # 追加案例 (写入 MITO_KNOWLEDGE_DIR/cases/X.md)
  update-signals --signal S --judgment J --ref 案例   # 追加信号
  update-pitfalls --pitfall P --fix F                # 追加坑
  suggest-promotions                # 模式提炼: 统计重复信号/动作, 建议固化
  stats                             # 显示知识库统计
"""
import sys, os, re, argparse, datetime

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DATA_ROOT = os.environ.get('XDG_DATA_HOME') or os.path.join(os.path.expanduser('~'), '.local', 'share')
KNOWLEDGE = os.environ.get('MITO_KNOWLEDGE_DIR') or os.path.join(DEFAULT_DATA_ROOT, 'mito-assembly-triage')
CASES = os.path.join(KNOWLEDGE, 'cases')
SIGNALS = os.path.join(KNOWLEDGE, 'signals.md')
PITFALLS = os.path.join(KNOWLEDGE, 'pitfalls.md')
STATS = os.path.join(KNOWLEDGE, 'stats.md')


def case_files():
    if not os.path.isdir(CASES):
        return []
    return sorted(fn for fn in os.listdir(CASES) if fn.endswith('.md'))


def ensure_cases_dir():
    os.makedirs(CASES, exist_ok=True)


def ensure_knowledge_files():
    ensure_cases_dir()
    defaults = [
        (SIGNALS, '# 信号映射表\n\n| 信号 | 判定 | 处理 | 来源 |\n|---|---|---|---|\n'),
        (PITFALLS, '# 常见工具坑\n\n| 坑 | 现象 | 解决 | 来源 |\n|---|---|---|---|\n'),
        (STATS, '# 样品处理统计\n\n暂无样品记录。\n'),
    ]
    for path, content in defaults:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if not os.path.exists(path):
            with open(path, 'x', encoding='utf-8') as fh:
                fh.write(content)

def search(query):
    print('=== 检索案例 (关键词: %s) ===' % query)
    kw = [k.lower() for k in query.split()]
    found = 0
    for fn in case_files():
        content = open(os.path.join(CASES, fn)).read()
        if all(k in content.lower() for k in kw):
            found += 1
            print('\n▶ %s' % fn)
            # 打印 frontmatter 摘要
            m = re.search(r'^---\n(.*?)\n---', content, re.S)
            if m:
                for line in m.group(1).split('\n'):
                    if ':' in line and not line.startswith(' '):
                        print('   %s' % line.strip())
            # 打印信号/判定/新知识小节
            for sec in ['信号', '判定', '新知识']:
                mm = re.search(r'## %s.*?\n(.*?)(?=\n## |\Z)' % sec, content, re.S)
                if mm:
                    lines = [l.strip() for l in mm.group(1).split('\n') if l.strip()][:3]
                    print('   [%s] %s' % (sec, ' | '.join(lines)))
    if not found:
        print('未找到匹配案例 — 这是新的问题类型! 处理完请 add-case 沉淀')
    return found

def case_path(sample):
    if not sample or not isinstance(sample, str) or '/' in sample or '\\' in sample or '\x00' in sample or os.path.basename(sample) != sample or sample in ('.', '..'):
        raise ValueError('sample must be a single file name')
    path = os.path.join(CASES, '%s.md' % sample)
    if os.path.islink(path):
        raise ValueError('case path must not be a symlink')
    return path

def add_case(args):
    ensure_knowledge_files()
    fn = case_path(args.sample)
    force = getattr(args, 'force', False)
    if os.path.exists(fn) and not force:
        raise FileExistsError(fn)
    content = """---
sample: %(sample)s
species: %(species)s
date: %(date)s
tools: %(tools)s
status: 进行中
signals: [%(signals)s]
---

# %(sample)s 线粒体组装诊断案例

## 信号 (Signal)
%(signals)s

## 判定 (Diagnosis)
%(diagnosis)s

## 动作 (Actions)
%(actions)s

## 新知识 (Lessons)
%(lessons)s
""" % {
        'sample': args.sample, 'species': args.species or '待鉴定',
        'date': datetime.date.today().isoformat(), 'tools': args.tools or '待记录',
        'signals': args.signals or '待补充', 'diagnosis': args.diagnosis or '待补充',
        'actions': args.actions or '待补充', 'lessons': args.lessons or '待补充',
    }
    with open(fn, 'w' if force else 'x') as fh:
        fh.write(content)
    print('✓ 案例写入: %s' % fn)
    print('  记得同时更新 signals.md / pitfalls.md (有新知识时) 和 stats.md')


def update_signals(args):
    ensure_knowledge_files()
    if not os.path.exists(SIGNALS):
        print('✗ 找不到 signals.md'); return
    line = '| %s | %s | %s | %s |' % (args.signal, args.judgment, '待处理', args.ref or '手动')
    with open(SIGNALS, 'a') as fh:
        fh.write(line + '\n')
    print('✓ 已追加信号: %s' % args.signal)

def update_pitfalls(args):
    ensure_knowledge_files()
    if not os.path.exists(PITFALLS):
        print('✗ 找不到 pitfalls.md'); return
    with open(PITFALLS, 'a') as fh:
        fh.write('| %s | %s | %s | %s |\n' % (args.pitfall, '见案例', args.fix or '待解决', datetime.date.today().isoformat()))
    print('✓ 已追加坑: %s' % args.pitfall)

def suggest_promotions():
    print('=== 模式提炼 (统计 signals.md 出现频次) ===')
    # 统计案例数
    n_cases = len(case_files())
    print('案例数: %d' % n_cases)
    # 从案例的 frontmatter 统计 signals 出现频次
    sig_count = {}
    for fn in case_files():
        content = open(os.path.join(CASES, fn)).read()
        m = re.search(r'^signals: \[(.*?)\]', content)
        if m:
            for s in m.group(1).split(','):
                s = s.strip()
                sig_count[s] = sig_count.get(s, 0) + 1
    if sig_count:
        print('\n信号频次:')
        for s, c in sorted(sig_count.items(), key=lambda x: -x[1]):
            mark = '  ★ 出现%d次 → 建议提升到 SKILL.md 主流程!' % c if c >= 3 else ''
            print('  %s: %d%s' % (s, c, mark))
    # 脚本固化建议
    print('\n建议: 当某手工操作重复 3 次以上, 固化为 scripts/ 下的脚本;')
    print('当某信号重复 3 次以上, 提升为 SKILL.md 的固定检查项。')

def stats():
    files = case_files()
    print('知识库位置: %s' % KNOWLEDGE)
    print('案例数: %d' % len(files))
    for f in files:
        c = open(os.path.join(CASES, f)).read()
        m = re.search(r'^status: (.*)', c, re.M)
        print('  - %s [%s]' % (f.replace('.md',''), (m.group(1) if m else '?')))
    print('\n文件:')
    for d in [CASES, SIGNALS, PITFALLS, STATS]:
        print('  %s' % d)

def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(0)
    cmd = sys.argv[1]
    if cmd == 'search':
        if '--query' not in sys.argv:
            print('需要 --query'); sys.exit(1)
        search(sys.argv[sys.argv.index('--query')+1])
    elif cmd == 'add-case':
        ap = argparse.ArgumentParser()
        for a in ['sample','species','tools','signals','diagnosis','actions','lessons']:
            ap.add_argument('--%s' % a)
        ap.add_argument('--force', action='store_true')
        try:
            add_case(ap.parse_args(sys.argv[2:]))
        except (ValueError, FileExistsError) as exc:
            print('✗ %s' % exc, file=sys.stderr)
            sys.exit(1)
    elif cmd == 'update-signals':
        ap = argparse.ArgumentParser()
        ap.add_argument('--signal', required=True); ap.add_argument('--judgment', required=True)
        ap.add_argument('--ref')
        update_signals(ap.parse_args(sys.argv[2:]))
    elif cmd == 'update-pitfalls':
        ap = argparse.ArgumentParser()
        ap.add_argument('--pitfall', required=True); ap.add_argument('--fix')
        update_pitfalls(ap.parse_args(sys.argv[2:]))
    elif cmd == 'suggest-promotions':
        suggest_promotions()
    elif cmd == 'stats':
        stats()
    else:
        print('未知命令: %s' % cmd); print(__doc__)

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""代码 ↔ 文档一致性审计（常驻回归护栏）。

为什么需要它：这个 skill 的规则同时活在两个地方——代码（真正强制执行的常量、退出码、
fail-closed 判定）和文档（SKILL.md 与 references/*.md，是 agent 与人的唯一计划依据）。
两边漂移过一次：一次文档压缩把 6 个 CLI flag 的**唯一**文档位置删掉了，其中 4 个是
`required=True`，照文档构造命令必然报错。这类丢失不会被单元测试发现，因为代码没变。

本工具做三件事，任何一件失败都以非 0 退出：

1. **可执行 token 覆盖**：从代码抽取 CLI flag（`--x`）、判定/错误码（`UPPER_SNAKE`）、
   枚举值（`CASE_TYPES` 等常量），要求每个 token 要么**有文档**，要么在
   `tools/doc_contract_baseline.json` 的 `allow_undocumented` 里**显式**豁免。
   新出现的 token 既没文档也没豁免 → 失败（强制"写文档，或明确豁免"这一次决定）。
   已经进过 `documented` 的 token 消失 → 失败（这就是上面那次漂移）。
2. **规则表**：一组精选的"代码强制行为"必须仍能在文档里检索到；
   每条规则自带**代码凭据**，凭据本身若不存在也失败（防止规则表腐化）。
3. **基线同步**：`--check-baseline` 时当前状态必须与基线文件一致。

用法：
    python3 tools/audit_code_docs.py                 # 审计（CI 用）
    python3 tools/audit_code_docs.py --json          # 机器可读
    python3 tools/audit_code_docs.py --update-baseline   # 明确豁免/收编后更新基线

更新基线是一个**需要理由的动作**：新增豁免意味着"这个 token 故意不写进文档"。
"""
import argparse
import json
import pathlib
import re
import sys

SKILL_ROOT = pathlib.Path(__file__).resolve().parent.parent
BASELINE_NAME = 'doc_contract_baseline.json'
BASELINE_FORMAT = 'mito-doc-contract-baseline-1'
DOC_FILES = ('SKILL.md', 'references')

# 代码里会被抽取 token 的文件（相对 skill 根）
CODE_FILES = (
    'tools/experience.py',
    'scripts/reference_registry.py',
    'scripts/annot_check.py',
    'scripts/cox1_id.py',
    'scripts/blast_genes.py',
    'scripts/depth_analysis.py',
    'scripts/circularize.py',
    'scripts/mitos2_to_genbank.py',
    'scripts/run_mitos2.sh',
    'scripts/run_bg.sh',
    'scripts/check_env.sh',
)

# 枚举来源：名字 -> 定义这些值的常量名
ENUM_CONSTANTS = (
    'CASE_TYPES', 'LESSON_DOMAINS', 'CASE_TYPE_LESSON_DOMAINS', 'GENERALIZATION_SCOPES',
    'TRANSFERABILITY_LEVELS', 'REFERENCE_KINDS', 'REFERENCE_PURPOSES', 'REFERENCE_LEVELS',
    'DECISION_STATUSES', 'CONFIDENCE_LEVELS', 'READS_SUPPORT_LEVELS', 'LESSON_STATUSES',
)

FLAG_RE = re.compile(r"['\"]?(--[a-z][a-z0-9-]{2,})")
CODE_RE = re.compile(r"['\"]([A-Z][A-Z0-9_]{3,})['\"]")
CONST_RE = re.compile(r"^([A-Z][A-Z0-9_]{3,})\s*=", re.M)
ENUM_VALUE_RE = re.compile(r"['\"]([A-Za-z0-9_]+)['\"]")

# (规则名, 代码文件, 代码凭据, [文档正则…])
RULES = (
    ("annot_check 退出码 2 ≠ 通过", 'scripts/annot_check.py', 'sys.exit(2)', (r'退出码?\s*2', r'exit\s*2')),
    ("拓扑声明≠物理闭环", 'scripts/annot_check.py', 'CIRCULAR_DECLARATION_CHECK', (r'CIRCULAR_DECLARATION_CHECK',)),
    ("transl_except 匹配≠接受", 'scripts/annot_check.py', 'TRANSL_EXCEPT_MATCHED', (r'TRANSL_EXCEPT_MATCHED',)),
    ("transl_except 未验证时内部终止保持 ERROR", 'scripts/annot_check.py', 'TRANSL_EXCEPT_DECLARED_UNVERIFIED',
     (r'DECLARED_UNVERIFIED',)),
    ("裸名 trnL/trnS 不自动归一", 'scripts/annot_check.py', 'UNDETERMINED_TRNA', (r'UNDETERMINED_TRNA', r'裸名')),
    (">8bp 重叠默认 REVIEW", 'scripts/annot_check.py', 'OVERLAP_REVIEW_THRESHOLD', (r'8\s*bp',)),
    ("重叠豁免=人工审核，非功能真实", 'scripts/annot_check.py', 'OVERLAP_ACCEPTED', (r'OVERLAP_ACCEPTED',)),
    ("不允许把重叠当裁剪理由", 'scripts/annot_check.py', 'OVERLAP_LONG', (r'裁剪',)),
    ("tRNA 长度是工程阈值", 'scripts/annot_check.py', 'TRNA_LENGTH_RANGE', (r'60\s*[–\-~]\s*75',)),
    ("rRNA 长度是工程阈值", 'scripts/annot_check.py', 'RRNA_LENGTH_RANGE', (r'1100', r'长度区间')),
    ("9+/4− 是昆虫预警", 'scripts/annot_check.py', 'EXPECTED_PLUS_STRAND_CDS', (r'9\s*\+', r'ORIENTATION')),
    ("例外未登记只维持 REVIEW", 'scripts/annot_check.py', 'EXCEPTION_NOT_REGISTERED', (r'EXCEPTION_NOT_REGISTERED',)),
    ("例外缺字段/taxon 不符", 'scripts/annot_check.py', 'EXCEPTION_RECORD_INCOMPLETE',
     (r'EXCEPTION_RECORD_INCOMPLETE', r'EXCEPTION_TAXON_MISMATCH')),
    ("未给 --taxon 无法升级", 'scripts/annot_check.py', 'EXCEPTION_TAXON_UNVERIFIED', (r'EXCEPTION_TAXON_UNVERIFIED',)),
    ("codon_start 与 location 矛盾", 'scripts/annot_check.py', 'CODON_START_CONFLICT', (r'CODON_START_CONFLICT', r'codon_start')),
    ("transl_table 冲突", 'scripts/annot_check.py', 'TABLE_CONFLICT', (r'TABLE_CONFLICT',)),
    ("基因集差异与顺序差异分开", 'scripts/annot_check.py', 'GENE_SET_DIFF', (r'GENE_SET_DIFF', r'ARRANGEMENT_DIFF')),
    ("顺序对照需显式 --ref", 'scripts/annot_check.py', '--ref', (r'--ref',)),
    ("异常入口只降基因集/身份", 'scripts/annot_check.py', '--allow-atypical', (r'allow-atypical',)),
    ("cox1_id 退出码 3 = 格式/网络故障", 'scripts/cox1_id.py', 'FORMAT', (r'退出码?\s*3', r'exit\s*3')),
    ("circularize 不匹配不得接受候选", 'scripts/circularize.py', 'accept-candidate', (r'accept-candidate', r'REVIEW')),
    ("案例类型三值", 'tools/experience.py', 'CASE_TYPES', (r'case_type', r'normal_validation_case')),
    ("逐异常判定写入 anomalies[]", 'tools/experience.py', 'case_anomaly', (r'anomalies', r'case-anomaly')),
    ("未登记参考不得引用", 'tools/experience.py', 'case_reference', (r'case-reference',)),
    ("未声明用途不得扩大使用", 'tools/experience.py', 'purposes', (r'未声明用途', r'等级限制用途')),
    ("lesson 领域由 case_type 限定", 'tools/experience.py', 'CASE_TYPE_LESSON_DOMAINS', (r'lesson_domain',)),
    ("推广上限默认 single_case/none", 'tools/experience.py', 'TRANSFERABILITY_MIN_CASES', (r'single_case', r'transferability')),
    ("公共 lesson 缺失范围字段保守补齐", 'tools/experience.py', 'validate_public_lesson', (r'公共知识', r'validate_public_lesson')),
    ("参考必须带版本", 'scripts/reference_registry.py', 'ACCESSION_VERSION', (r'ACCESSION\.VERSION', r'版本')),
    ("等级限制用途矩阵", 'scripts/reference_registry.py', 'LEVEL_PURPOSES', (r'LEVEL_PURPOSES', r'L5')),
    ("版本漂移即错误", 'scripts/reference_registry.py', 'hash', (r'漂移',)),
    ("数据库记录需版本", 'scripts/reference_registry.py', 'directory_sha256', (r'record-database', r'数据库')),
    ("verify 检测缺失/变更", 'scripts/reference_registry.py', 'verify_registry', (r'verify',)),
    ("注册表不进 skill 安装目录", 'scripts/reference_registry.py', 'MITO_KNOWLEDGE_DIR', (r'MITO_KNOWLEDGE_DIR',)),
    ("后台任务不算完成", 'scripts/run_bg.sh', 'status', (r'run_bg', r'完成')),
)


def read_text(path):
    try:
        return path.read_text(encoding='utf-8', errors='replace')
    except OSError:
        return ''


def doc_text(root):
    parts = []
    skill_md = root / 'SKILL.md'
    if skill_md.is_file():
        parts.append(read_text(skill_md))
    references = root / 'references'
    if references.is_dir():
        for item in sorted(references.rglob('*.md')):
            parts.append(read_text(item))
    return '\n'.join(parts)


def code_text(root):
    out = {}
    for rel in CODE_FILES:
        path = root / rel
        if path.is_file():
            out[rel] = read_text(path)
    return out


def extract_tokens(code):
    flags, codes, enums = set(), set(), set()
    for rel, text in code.items():
        for match in FLAG_RE.finditer(text):
            flags.add(match.group(1))
        for match in CODE_RE.finditer(text):
            codes.add(match.group(1))
        for match in CONST_RE.finditer(text):
            codes.add(match.group(1))
        for name in ENUM_CONSTANTS:
            found = re.search(re.escape(name) + r"\s*=\s*[\(\{]([^\)\}]*)", text)
            if found:
                for value in ENUM_VALUE_RE.findall(found.group(1)):
                    enums.add(value)
    # 过于通用、不承载接口语义的 token 不进审计
    trivial = {'--help', '--version'}
    tokens = {}
    for tok in sorted(flags - trivial):
        tokens[tok] = 'flag'
    for tok in sorted(codes):
        if len(tok) >= 4:
            tokens[tok] = 'code'
    for tok in sorted(enums):
        if len(tok) >= 3:
            tokens.setdefault(tok, 'enum')
    return tokens


def audit(root, baseline):
    code = code_text(root)
    docs = doc_text(root)
    tokens = extract_tokens(code)
    documented_now = sorted(t for t in tokens if t in docs)
    undocumented_now = sorted(t for t in tokens if t not in docs)
    allowed = set(baseline.get('allow_undocumented') or [])
    known_documented = set(baseline.get('documented') or [])

    problems = []
    # 1) 曾经有文档的 token 不得变成无文档（那次漂移就是这个）
    for token in sorted(known_documented):
        if token not in tokens:
            problems.append(('token-removed-from-code', token,
                             '基线认为它有文档，但代码里已不存在该 token'))
        elif token not in docs:
            problems.append(('documentation-lost', token,
                             '该 token 曾有文档，现在文档里找不到（丢失或改写到无法检索）'))
    # 2) 新 token 必须"写文档 或 明确豁免"
    for token in undocumented_now:
        if token in allowed or token in known_documented:
            continue
        problems.append(('new-token-undocumented', token,
                         '代码里新出现但文档中不存在；写进文档，或在基线 allow_undocumented 中显式豁免'))
    # 3) 规则表必须仍可检索，且其代码凭据必须真实存在
    for name, rel, evidence, patterns in RULES:
        text = code.get(rel, '')
        if evidence not in text:
            problems.append(('rule-evidence-missing', name,
                             '规则表引用的代码凭据在 %s 中找不到: %s' % (rel, evidence)))
        if not any(re.search(p, docs) for p in patterns):
            problems.append(('rule-undocumented', name,
                             '该代码强制行为在文档中检索不到（模式: %s）' % ' | '.join(patterns)))
    return {
        'tokens_total': len(tokens),
        'documented': documented_now,
        'undocumented': undocumented_now,
        'problems': problems,
        'baseline': {'documented': known_documented, 'allow_undocumented': allowed},
    }


def baseline_path(root):
    return root / 'tools' / BASELINE_NAME


def load_baseline(root):
    path = baseline_path(root)
    if not path.is_file():
        return {'format': BASELINE_FORMAT, 'documented': [], 'allow_undocumented': []}
    data = json.loads(read_text(path))
    if data.get('format') != BASELINE_FORMAT:
        raise ValueError('基线格式不兼容（期望 %s）: %s' % (BASELINE_FORMAT, path))
    return data


def save_baseline(root, result):
    data = {
        'format': BASELINE_FORMAT,
        'comment': ('documented: token 曾出现在文档里，后续改动不得让它消失；'
                    'allow_undocumented: 明确决定不写进文档的标识；'
                    '新增 token 必须二选一，否则审计失败。'),
        'exemption_policy': ('逐工具的调参 flag（如 --min-mapq/--evalue/--junction-region）'
                             '不逐条写进 references：SKILL.md 已声明"参数以 --help 为准"，'
                             'tool-catalog.md 只描述各工具能回答什么、不能证明什么。'
                             '豁免的是"参数细节"，不是"接口语义"：会改变结论语义的 flag'
                             '（--authorize/--allow-public-upload/--input/--next-test 等）'
                             '必须保持有文档。'),
        'documented': result['documented'],
        'allow_undocumented': result['undocumented'],
    }
    baseline_path(root).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n',
                                   encoding='utf-8')


def main(argv=None):
    parser = argparse.ArgumentParser(description='代码 ↔ 文档一致性审计')
    parser.add_argument('--root', default=str(SKILL_ROOT), help='skill 根目录')
    parser.add_argument('--json', action='store_true')
    parser.add_argument('--update-baseline', action='store_true')
    parser.add_argument('--check-baseline', action='store_true',
                        help='额外要求当前状态与基线完全一致')
    args = parser.parse_args(argv)
    root = pathlib.Path(args.root).resolve()
    baseline = load_baseline(root)
    result = audit(root, baseline)

    if args.update_baseline:
        save_baseline(root, result)
        print('基线已更新: %s（documented=%d, allow_undocumented=%d）'
              % (baseline_path(root), len(result['documented']), len(result['undocumented'])))
        return 0

    if args.json:
        print(json.dumps({'tokens_total': result['tokens_total'],
                          'documented': len(result['documented']),
                          'undocumented': len(result['undocumented']),
                          'problems': [{'kind': k, 'subject': s, 'detail': d}
                                       for k, s, d in result['problems']]},
                         ensure_ascii=False, indent=1))
    else:
        print('token 总数 %d：有文档 %d，无文档 %d'
              % (result['tokens_total'], len(result['documented']), len(result['undocumented'])))
        if result['problems']:
            print('\n发现 %d 个问题：' % len(result['problems']))
            for kind, subject, detail in result['problems']:
                print('  [%s] %s —— %s' % (kind, subject, detail))
        else:
            print('无漂移：所有 token 都有文档或已显式豁免，规则表可检索。')

    if args.check_baseline:
        current = {'documented': result['documented'], 'allow_undocumented': result['undocumented']}
        saved = {'documented': sorted(baseline.get('documented') or []),
                 'allow_undocumented': sorted(baseline.get('allow_undocumented') or [])}
        if current['documented'] != saved['documented'] or \
                current['allow_undocumented'] != saved['allow_undocumented']:
            print('基线文件与当前状态不一致；请用 --update-baseline 明确更新', file=sys.stderr)
            return 1

    return 1 if result['problems'] else 0


if __name__ == '__main__':
    sys.exit(main())

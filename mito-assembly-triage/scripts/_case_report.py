"""Render evidence without inferring its strength from decision status."""
import html
import json
import re


def escape(value):
    text = html.escape(str(value), quote=False)
    text = re.sub(r'([\\`*_{}\[\]()#+.!|>-])', r'\\\1', text)
    return text.replace('\r', '').replace('\n', '<br>')


def render(case, case_hash):
    lines = ['# 线粒体诊断任务报告', '',
             '本报告整理已登记的证据。格式、关联和文件哈希校验通过不证明科学结论；'
             '执行元数据及判定/复核身份为登记内容，本引擎没有代执行命令或认证身份。', '',
             '- 案例：' + escape(case['case_id']),
             '- 类型：' + escape(case['case_type']),
             '- 类群：' + escape(case['taxon']),
             '- 任务：' + escape(case['issue']),
             '- 修订号：' + str(case['revision']),
             '- case.json SHA-256：' + case_hash, '',
             '## 检查范围与逐条结论', '']
    if not case['conclusions']:
        lines += ['尚未登记结论；不能据此判断样本正常。', '']
    else:
        lines += ['| ID | 命题 | 类型 | 状态 | 置信度 | reads 支持 |', '|---|---|---|---|---|---|']
        for conclusion in case['conclusions']:
            lines.append('| ' + ' | '.join(escape(conclusion[key]) for key in
                ('id', 'claim', 'conclusion_type', 'status', 'confidence', 'reads_support')) + ' |')
        lines += ['', 'NO_CHANGE 表示本次不修改；UNRESOLVED 可具有部分支持或相互冲突的证据。', '']
    events = {event['id']: event for event in case['events']}
    for conclusion in case['conclusions']:
        lines += ['### ' + escape(conclusion['id']) + '：' + escape(conclusion['claim']), '',
                  '**适用范围：** ' + escape(conclusion['scope']), '',
                  '**判定理由：** ' + escape(conclusion['rationale']), '']
        for field, title in [('evidence_for', '支持证据'), ('evidence_against', '反对证据')]:
            lines += ['**' + title + '：**', '']
            for item in conclusion[field]:
                event = events[item['event_id']]
                artifact = next(a for a in event['artifacts'] if a['id'] == item['artifact_id'])
                lines.append('- ' + escape(item['observation']) + '；事件 ' + escape(event['id']) +
                             ' / 产物 ' + escape(artifact['id']) + '；文件 ' + escape(artifact['path']) +
                             '；SHA-256 ' + artifact['sha256'] +
                             ('；坐标 ' + escape(item['coordinates']) if item.get('coordinates') else ''))
            if not conclusion[field]:
                lines.append('- 未登记；不表示已排除相反解释。')
            lines.append('')
        for field, title in [('not_tested', '未检查'), ('limitations', '局限'), ('next_steps', '下一步')]:
            lines += ['**' + title + '：**', '']
            lines += ['- ' + escape(value) for value in conclusion[field]] or ['- 未登记。']
            lines.append('')
        review = conclusion['review']
        lines += ['**判定人（登记）：** ' + escape(conclusion['decided_by']), '',
                  '**复核记录：** ' + escape(review['level']) + ' / ' + escape(review['status']) +
                  '；复核人 ' + escape(review['reviewer'] or '未登记') + '；' + escape(review['note']), '',
                  '**竞争参考输入 ID：** ' + escape(', '.join(conclusion['competitive_reference_ids']) or '未登记'), '']
    lines += ['## 修复候选与验证', '']
    if not case.get('candidates'):
        lines.append('尚未登记候选；本报告不表示无需修复。')
    for candidate in case.get('candidates', []):
        item = next(i for i in case['inputs'] if i['id'] == candidate['candidate_input_id'])
        lines += ['### ' + escape(candidate['id']) + '：' + escape(candidate['status']) +
                  ('（仅表示登记的验证标准结果）' if candidate['status'] == 'verified' else ''), '',
                  '- 修改类型：' + escape(candidate['change_type']),
                  '- 基础输入：' + escape(candidate['base_input_id']) + '；候选输入：' + escape(item['id']),
                  '- 生成事件/产物：' + escape(candidate['producer_event_id']) + ' / ' + escape(candidate['artifact_id']),
                  '- 候选文件：' + escape(item['path']) + '；SHA-256 ' + item['sha256'],
                  '- 修改理由：' + escape(candidate['rationale']),
                  '- 预期影响：' + escape(candidate['expected_effect']),
                  '- 风险：' + escape('; '.join(candidate['risks']) or '未登记'), '']
        if not candidate['checks']:
            lines.append('尚无修复后验证事件；候选仍未验证。')
        for check in candidate['checks']:
            lines.append('- 验证事件 ' + escape(check['event_id']) + '：' + escape(check['result']) +
                         '；标准：' + escape(check['criteria']) + '；记录：' + escape(check['notes']))
        lines.append('')
    lines += ['## 输入清单', '', '| ID | 角色 | 路径 | 大小 | SHA-256 |', '|---|---|---|---|---|']
    for item in case['inputs']:
        lines.append('| ' + ' | '.join(escape(item[k]) for k in ('id', 'role', 'path', 'size', 'sha256')) + ' |')
    lines += ['', '## 执行事件（登记记录）', '']
    for event in case['events']:
        lines += ['### ' + escape(event['id']) + '：' + escape(event['action']), '',
                  '- 登记时间：' + escape(event['recorded_at']),
                  '- 执行状态：' + escape(event['execution_status']) + '；实际退出码：' + escape(event['exit_code']),
                  '- 命令参数数组：' + escape(json.dumps(event['command'], ensure_ascii=False)),
                  '- 工作目录：' + escape(event['working_directory']),
                  '- 工具版本：' + escape(event['tool_version']),
                  '- 输入 ID：' + escape(', '.join(event['input_ids'])),
                  '- 观察记录：' + escape(event['summary'])]
        for artifact in event['artifacts']:
            lines.append('- 产物 ' + escape(artifact['id']) + ' (' + escape(artifact['role']) + ')：' +
                         escape(artifact['path']) + '；SHA-256 ' + artifact['sha256'])
        lines.append('')
    if case['conclusion_history']:
        lines += ['## 结论修订历史', '', '完整旧证据和理由保存在 case.json 的 conclusion_history。', '']
        for entry in case['conclusion_history']:
            old = entry['conclusion']
            lines.append('- 修订 ' + str(entry['revision']) + '，' + escape(old['id']) + '：' +
                         escape(old['claim']) + ' (' + escape(old['status']) + ')；替换登记时间 ' + escape(entry['replaced_at']))
    return '\n'.join(lines) + '\n'

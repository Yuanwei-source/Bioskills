# 经验学习与共享政策

> 当前实现是本地 MVP：`scripts/memory.py` 从版本 3 案例提出经验，经审核后供后续样本检索。
> 跨任务写入必须显式授权；公共导出、上传和同步尚未实现。经验不得自动修改 `SKILL.md`、脚本或科学规则。

## 1. 事实源、引用与索引

| 层 | 位置 | 用途与边界 |
|---|---|---|
| 任务案例 | 用户任务目录的 v3 `case.json` 与报告 | 保存单次样本证据；`memory.py` 只读，不复制或修改案例 |
| 本地经验 | `<经验库>/lessons/<lesson-id>.json` | 每个文件保存同一经验的全部修订；这是事实源 |
| 派生索引 | `<经验库>/index.sqlite3` | FTS5 检索索引，可删除后用 `rebuild` 重建，不是独立证据 |
| 原始测序数据 | 用户自己的存储 | 不进入经验库；经验只保留泛化描述及来源案例路径、ID 和哈希 |

经验库必须由 `--library` 或 `MITO_KNOWLEDGE_DIR` 显式指定，不写入 skill 安装目录。
所有新增、审核、修订和索引重建都要求 `--authorize-local`。`search` 与 `show` 只读。
本地存储授权不等于共享授权；当前没有公共导出或上传命令。

## 2. 每个样本如何形成经验

先完成 `case.py` 的结论和必要验证，再撰写一个经验草稿 JSON。最小字段示例：

```json
{
  "title": "nad6 断裂后的竞争回贴",
  "claim": "单一注释断裂不能证明样本碱基错误",
  "lesson_domain": "annotation",
  "diagnostic_clues": ["nad6 断裂", "局部覆盖下降"],
  "suggested_next_test": "检查组装图并对竞争候选回贴",
  "applicable_when": {
    "taxa": ["Insecta"],
    "data_types": ["illumina_paired_end"],
    "conditions": ["记录实际组装器与版本"]
  },
  "not_applicable_when": ["没有当前样本证据时直接修改碱基"],
  "counterexample_case_ids": [],
  "conflicts": []
}
```

不要写工具自动生成的状态、案例路径、哈希、审核历史或推广级别；程序会补全。
MVP 每条修订只关联一个来源案例，并固定为 `generalization_scope=single_case`、
`transferability=none`。不同样本分别提出经验，不通过频次自动合并或推广。

案例类型保留在经验中：

- `abnormal_case`：支持异常诊断线索，但不证明适用于其他类群或数据类型；
- `normal_validation_case`：作为正常验证对照，标签本身不证明真实阴性；
- `tool_failure_case`：只能提出 `lesson_domain=tool` 的经验，不得升级为生物学结论。

案例含 `UNRESOLVED` 结论时，只有显式 `--allow-unresolved` 才能生成候选；
含未解决结论的来源不能审核为 `verified`。

## 3. 当前命令

```bash
# 提出候选；案例和草稿先校验，失败不创建经验库
python3 scripts/memory.py --library /path/to/knowledge propose \
  --case /path/to/task-case --record lesson-draft.json --authorize-local
# 含 UNRESOLVED 时另加 --allow-unresolved

# 审核当前修订；reviewer 是自述字段，不构成独立验证证据
python3 scripts/memory.py --library /path/to/knowledge review L000001 \
  --status verified --reviewer '姓名或审核标识' --reason '实际核查内容' --authorize-local

# 默认只检索 verified；类群和数据类型为大小写不敏感的精确范围匹配
python3 scripts/memory.py --library /path/to/knowledge search \
  --query 'nad6 断裂' --taxon Insecta --data-type illumina_paired_end --domain annotation
# 查候选用 --status candidate；查当前全部状态用 --status all

# 读取完整经验、证据引用和审核历史
python3 scripts/memory.py --library /path/to/knowledge show L000001

# 同一来源案例的新修订；新修订重新从 candidate 开始，旧修订保留
python3 scripts/memory.py --library /path/to/knowledge supersede L000001 \
  --case /path/to/task-case --record revised-lesson.json --authorize-local
python3 scripts/memory.py --library /path/to/knowledge show L000001 --revision 1

# 删除索引或源文件变化后，从 JSON 事实源重建
python3 scripts/memory.py --library /path/to/knowledge rebuild --authorize-local
```

`--help` 为当前参数的唯一权威。旧 `experience.py`、`propose-lesson`、`review-lesson`、
`case-*` 和公共 `sync` 命令仍不可执行，不能按历史文档调用。

## 4. 生命周期与证据边界

```text
candidate ──review──► verified ──review──► deprecated / withdrawn
    │                      │
    └──review──► rejected  └──supersede──► 新 candidate 修订
```

- `supersede` 可用于当前任意状态；旧修订标为 `superseded`，内容仍可 `show --revision` 追溯。
- 不允许的状态转换会失败，不会静默覆盖经验。
- 来源案例在提案后发生合法修订也会改变哈希；审核为 `verified` 前必须重新核对并生成新修订。
  科学内容不变也可用 `supersede` 刷新来源哈希；新修订仍从 `candidate` 开始。
- 相同来源与科学内容的当前经验会去重；内容哈希包含来源路径和案例哈希，只用于身份与去重，不证明科学正确。
  不同目录中的相同 `case_id` 不会被视为同一案例，也不能串接修订链。
- `verified` 仅表示记录了审核状态，不表示普遍适用或已经独立验证。
- `conflicts` 保留相反经验 ID；检索不会按支持次数投票或删除任一方。
- 反例和未解决项必须保留；不能通过改写标题把它们包装成定稿。

召回的经验只是下一次诊断的候选线索。当前样本仍需执行相应检查，经验不能替代
reads、组装图、注释翻译或其他必要证据，也不能自动改变遗传密码表或修复阈值。

## 5. 检索与索引约束

检索先应用状态、领域、类群和数据类型过滤，再按 FTS5 关键词排序，默认返回最多 5 张摘要卡片。
`show` 才读取完整来源与审核记录。MVP 不推断分类层级：`Insecta` 不会自动匹配只标记
`Diptera` 的经验；需要更宽适用范围时，必须在草稿中明确写出并经审核。

索引记录源文件指纹。JSON 被修改或索引缺失时，搜索拒绝静默使用旧结果，要求 `rebuild`。
重建会重新校验所有经验文档。索引不是科学证据，也不能防止有文件写权限者篡改源文件。

## 6. 隐私与共享

经验草稿只写泛化后的诊断模式。写入端会拒绝常见绝对路径、邮箱、凭据模式、FASTA 片段及
整段 DNA 字符串，但这是辅助检查，不替代人工隐私审查。来源案例的绝对路径、案例 ID 和哈希
只用于本地回溯，也可能关联未公开项目，不得直接公开。

不要把原始 FASTQ、FASTA、BAM、GFA、未公开序列或采集元数据放入经验草稿。
新建经验目录权限为 `0700`，JSON 与 SQLite 索引权限为 `0600`；已有目录的父级权限仍需用户自行核查。
公共导出、去标识化、贡献审核和同步仍属后续阶段；在其实现前，不得把整个本地经验库当作已脱敏内容分享。

## 7. 晋升与后续阶段

稳定经验晋升到 `SKILL.md`、固定检查或脚本，必须另经人工审核、合适反例和回归测试；
同一 AI 的重复判断不是独立证据。工具版本变化只触发重新核查，不自动删除生物学经验。

先用真实样本评估召回是否有用、范围过滤是否过严以及是否干扰当前诊断，再决定是否增加
多案例合并、语义向量重排或公共共享。自动修改 Skill、按频次晋升和 LLM 自动裁决冲突不属于后续默认路线。

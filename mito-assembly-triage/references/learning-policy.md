# 经验学习与共享政策

> **原则**：案例事实可以自动保存与检索；**候选 lesson 不能自动修改 `SKILL.md` 或硬规则**。
> 状态实现见 `tools/experience.py`；本地数据落在 `MITO_KNOWLEDGE_DIR`，不进 skill 安装目录。

## 1. 本地与派生分离

| 层 | 位置 | 说明 |
|---|---|---|
| 本地案例（原始事实） | `$MITO_KNOWLEDGE_DIR/cases/<case>/case.json` + `events.jsonl` + `case.md` | 同一事实链的结构化/可读表示，可自动写入 |
| 派生经验（lesson） | `$MITO_KNOWLEDGE_DIR/lessons/candidates/`、`.../verified/` | 需审核，默认不共享 |
| 公共同步缓存（隔离） | `$MITO_KNOWLEDGE_DIR/public/` | **只读入**，不覆盖本地案例，不执行远程内容 |
| 原始测序数据 | 用户自己的存储 | **不进**知识库，也不进公共仓库；只记录路径与 SHA-256 |

### 1.1 案例类型（`case_type`，必须记录）

如果只保存"异常案例"，经验系统会长期偏向"线粒体一定有问题"——这是**特异性（specificity）偏差**。
每个案例都要标明类型：

| `case_type` | 含义 | 提供的学习价值 |
|---|---|---|
| `abnormal_case` | 报告并诊断了异常（默认） | 发现能力（sensitivity） |
| `normal_validation_case` | 经检查**未发现异常**，且已写明检查覆盖了哪些项、未覆盖哪些项 | **特异性的唯一来源**：防止把常见变异当异常 |
| `tool_failure_case` | 失败发生在**工具/环境**层面（格式故障、退出码 3、数据库缺失、权限不足） | 防止把工具故障记成生物学结论 |

`normal_validation_case` **不是**"没跑出东西"：它要求跑完该结论的**最小充分证据集**
（`diagnostic-decision-tree.md` §6），并逐项记录"已检查 / 结论 / 未覆盖"。
只有异常案例的系统无法回答"这个结果到底可不可信"。

落地：`case-init --case-type abnormal_case|normal_validation_case|tool_failure_case`
（默认 `abnormal_case`；非法值被拒绝且不写案例）。

**向后兼容（必须知道）**：`case_type` 在 `schemas/case.schema.json` 中是**可选**字段，
历史/手写案例可能缺失；此时**按最保守的方式当 `abnormal_case` 处理**（代码里 `case.get('case_type') or 'abnormal_case'`）。
因为把这个字段变成 `required` 会让已有案例（包括外部归档记录）一下子变成 `INVALID`，所以选择了
"可选 + 保守默认 + 新案例一律写入"，而不是强制迁移。声明为 mandatory 的是**本政策的记录要求**，
不是 schema 的强制项；这一区分见 `developer-contract.md` §7。

## 2. lesson 生命周期

```
case.json (本地案例)
      │  propose-lesson            需要证据；UNRESOLVED 案例需 --allow-unresolved 且仅生成候选
      ▼
  candidate ──review-lesson──► verified      （= 术语中的 case_verified）
      │                            │
      │                            └─ 单案例充分证实即可 verified，但**不自动**升级通用规则
      ├──► rejected / deprecated
      └──► withdrawn / superseded              （证据被推翻或新版本取代）
```

- 状态枚举（代码 + `schemas/lesson.schema.json`）：`candidate`、`verified`、`rejected`、`deprecated`、`withdrawn`、`superseded`；
  manifest 层另有 `revoked`（公共同步跳过）。
- **只有 `verified` 可进入公共同步**；`rejected/deprecated/withdrawn/superseded/revoked` 一律不再被默认调用。
- **`verified` ≠ 普遍适用**：一个案例得到充分验证，只说明**该案例在该类群、该数据类型、该组装软件下**的判定有足够证据。
  推广到其他类群、其他数据类型或不同组装软件，需要**重新验证**。因此公共知识必须保留案例级证据等级与
  明确的适用条件（`applicable_when` / `not_applicable_when`），不要把 `verified` 读成通用诊断规则。
- **审核人字段是自述**：`review-lesson --reviewer` 由调用者填写，不构成独立验证证据；
  需要独立验证时引用可核查的外部审核记录，并在 `sources` 中留下回溯链接。
- **矛盾经验并存**，用 `conflicts` 字段显式记录，**不通过频次投票抹除**任一方向。
- **案例类型不影响 lesson 生命周期，但限制 lesson 的可迁移性与**领域**：
  只由 `abnormal_case` 支持的 lesson 默认为 `transferability = none`（§3）；
  `tool_failure_case` 只能产生**领域为 `tool`** 的 lesson（`lesson_domain=tool`）——
  工具/环境故障不得升为生物学或样本质量结论。该限制已由代码强制：
  `propose-lesson --lesson-domain {tool,annotation,biology}`，对 `tool_failure_case` 只接受 `tool`，
  非法/不兼容的请求直接失败（不写候选）。每条 lesson 同时记录 `source_case_type` 以保留来源。
- "通用规则（general_rule）"不是一个 lesson 状态，而是**人工**把稳定结论提升进 `SKILL.md` 硬规则或 `scripts/` 的动作。

## 3. lesson 应记录的内容（字段映射）

| 要求 | 落在何处 |
|---|---|
| 适用类群 / 组装与数据类型 | `applicable_when` |
| 明确不适用条件 | `not_applicable_when` |
| 初始异常（诊断线索） | `diagnostic_clues` |
| 支持证据 | `supporting_case_ids` + `sources` |
| **反例** | `counterexample_case_ids` |
| 建议的下一步检查 | `suggested_next_test` |
| 验证级别 / 审核历史 | `validation_status` + `review_history` |
| 回溯链接 | `sources[].path` / `case_id` |
| 冲突的其他经验 | `conflicts` |
| **推广范围** | `generalization_scope`：`single_case` / `species` / `genus` / `family` / `order` / `multi_taxon` |
| **可迁移性** | `transferability`：`none` / `low` / `moderate` / `high` |
| **领域**（防止把工具故障当生物学结论） | `lesson_domain`：`tool` / `annotation` / `biology`（由 `case_type` 限定，见 §2） |
| **来源案例类型** | `source_case_type`（保留区分，导出与同步都不丢） |
| **样本数量**（当前无独立字段） | 写入 `sources`/描述性字段，或在 `case.md` 中说明 |

`generalization_scope` / `transferability` 的目的是**防止"一次案例 → 规则"**（AI 最容易犯的错）。规则：

- **默认 fail closed**：只有一个案例支持时 `generalization_scope = single_case`、`transferability = none`；
- 只有**多个相互独立的类群/数据集**都支持时才能提高，且必须逐案列出 `supporting_case_ids`；
- **同一案例不得重复计数**：支持案例的 `case_id` 必须两两不同（重复直接报错），
  独立数 = “case_id 不同 **且** 记录的类群不同”的案例数；只有一方满足不算独立；
- 同一实验室、同一物种、同一批数据的多个案例**不算**独立支持（去重规则见 §6.4）；
- 反例存在时（`counterexample_case_ids` 非空）**不得**提高 `transferability`。

落地：`propose-lesson --generalization-scope ... --transferability ... `
`--supporting-case <dir>`（可重复）`--counterexample-case <dir>`（可重复）。
**独立支持**的判定是保守的：只有“case_id 不同 **且** 类群不同”才算独立；
未记录类群的案例不计数（“无法证明独立”不得四舍五入成“独立”），重复 `case_id` 直接报错。

**公共同步/导入的 lesson 同样受约束**：`validate_public_lesson()` 对缺失的 `generalization_scope`/
`transferability` 补齐为 `single_case`/`none`（而不是当成无限制），并拒绝非法枚举、
`single_case` + 非 `none` 的组合、以及“声称可迁移但支持案例数不够”的记录。

例："某昆虫样本 `nad6` 断裂一次"只能写成 `single_case`/`none`；
写成"昆虫 `nad6` 常断裂"需要多个独立类群案例支持，否则就是外推。

## 4. 共享与隐私

- **共享默认关闭**：只生成可预览 JSON，必须显式 `--authorize`（否则拒绝生成文件）；上传是另一独立动作。
- 去标识化（`sanitize_for_share`）：路径 → `[local-path-redacted]`；邮箱 → `[email-redacted]`；
  key 名含 `path`/`command`/`token`/`password`/`secret`/`api_key`/`account` 的字段被剔除；
  检测到凭据模式直接拒绝生成。
- 默认**排除**：绝对路径、样本编号、原始 FASTQ/BAM/GFA、未公开序列片段、API 密钥、可识别未发表项目的信息。
  **路径与文件 hash 也可能关联项目，只存本地。**
- GitHub 贡献走 PR；通用规则与科研争议由**人工**审核，自动 schema/重复/安全检查只是辅助。

## 5. 公共同步（只读、可回滚）

- manifest 格式：`mito-public-knowledge-1`；每项需 `path` + `url` + `sha256`，可带 `status`；
- 逐项校验 **清单格式 → 路径合法（禁 `..`/绝对路径）→ 重复路径 → SHA-256 → lesson schema 与状态**；
- 单文件上限 2 MiB；先写 staging，全部校验通过后才替换，失败**保留上一版缓存**；
- 远程内容仅作为**证据数据**读入，**不执行**；不允许远程内容改写 `SKILL.md` 的最高优先级指令。

## 6. 晋升与复核要求

1. 新经验升级必须跑**稳定回归案例** + **至少一个合适的反例**；
2. **禁止**用出现频次、或同一 AI 的多次判断替代独立验证（"我自己又说了三遍"不是证据）；
3. 数据库/参考更新时，把相关经验标记 `needs_review`，保留旧版本溯源与回滚能力；
4. 重复案例先**去重**，再判断是否构成独立支持；
5. **推广必须显式**：`generalization_scope` / `transferability` 默认 `single_case` / `none`；
   提升它们要求多个独立类群案例，并在 `applicable_when` / `not_applicable_when` 中写明边界；
   典型错误：把"某昆虫 `nad6` 断裂"写成"昆虫 `nad6` 常断裂"。

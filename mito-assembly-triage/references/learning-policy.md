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
- **矛盾经验并存**，用 `conflicts` 字段显式记录，**不通过频次投票抹除**任一方向。
- "通用规则（general_rule）"不是一个 lesson 状态，而是**人工**把稳定结论提升进 `SKILL.md` 硬规则/`scripts/` 的动作。

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
| **样本数量**（当前无独立字段） | 写入 `sources`/描述性字段，或在 `case.md` 中说明 |

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
4. 重复案例先**去重**，再判断是否构成独立支持。

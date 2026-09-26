# 结论层与报告（你改变了什么）

> 异常清单（`anomalies[]`）是**过程记录**，不是结论。专家/合作者/审稿人真正需要的是：
> **"这个结果改变了什么？"** 本文件定义结论层、证据矩阵与固定输出模板。
> 结论的证据门槛见 `evidence-standard.md`；分流见 `diagnostic-decision-tree.md`。

## 1. 为什么要有结论层

一个案例里可以同时有"零覆盖缺口""缺 trnI""nad2 读框崩坏""末端低覆盖"四条异常，
但它们对读者的意义是分层级的：

| 异常（过程） | 结论（结果） |
|---|---|
| 64 bp 零覆盖缺口、末端 0.1× + soft-clip | **组装不完整**（`assembly_quality`） |
| 缺 tRNA、基因碎片化、内部 stop | **注释不可靠**（`annotation_quality`） |
| 高深度局部区块 | **来源未定**（`biological_interpretation`，`UNRESOLVED`） |
| BLAST 高 identity | **身份线索**，不是物种鉴定（`gene_identity`） |

结论必须能被单独引用、单独降级；**不允许**用案例级 `decision` 代替逐条结论。

## 2. 结论类型（`conclusion_type`）

| 类型 | 结论在说什么 | 典型证据 | 天然上限 |
|---|---|---|---|
| `assembly_quality` | 组装是否完整/是否可信（缺口、接缝、末端、重复） | 覆盖剖面、跨接缝 reads、端部唯一性 | 无 reads → `not_assessable` |
| `annotation_quality` | 注释是否自洽、基因集是否可靠 | `annot_check.py`、重译、同源、RNA 结构 | 注释自洽 ≠ 序列正确 |
| `gene_identity` | 某个基因/片段的身份与边界 | 同源覆盖、反密码子、邻域、参考等级 | 参考距离限制（决策树 §5） |
| `sequence_accuracy` | 具体碱基/片段是否正确 | reads pileup 四项（深度/MAPQ/碱基质量/链向） | 无 reads → 不得声称已校正 |
| `biological_interpretation` | 生物学解释（丢失、重排、退化、来源） | 先排除技术与注释原因 + 类群文献 | 技术原因未排除 → 不得作为结论 |

一条结论**只能属于一个类型**；跨类型时拆成两条（例如"该区无 reads 支持"(assembly) 与
"该基因真实缺失"(biological) 是两条，证据不同）。

## 3. 结论对象的必备字段

每条结论至少包含（与 `SKILL.md` 的"每项结论至少包含"一致）：

| 字段 | 含义 | 常见错误 |
|---|---|---|
| `conclusion_type` | 上表的类型 | 把 `annotation_quality` 写成 `biological_interpretation` |
| `claim` | 一句话结论（可被反驳的陈述） | 写成"需要进一步分析"这类不可反驳句 |
| `status` | `RESOLVED` / `NO_CHANGE` / `UNRESOLVED` | 用 `RESOLVED` 覆盖仍存疑的部分 |
| `confidence` | `high`/`moderate`/`low`/`not_assessable`（门槛见 `evidence-standard.md` §2） | 一个全局 `high` 覆盖全部结论 |
| `evidence_for` / `evidence_against` / `not_tested` | 支持/反对/未测 | 只写支持，把"没测"写成"没有" |
| `scope` / `limitations` | 适用类群/数据类型 + 局限 | 把单样本结论写成类群规律 |
| `decided_by` | `AUTO` / `ASSIST` / `EXPERT`（决策树 §4） | `EXPERT` 级事项无人确认却标 `RESOLVED` |
| 证据指针 | 引用**真实存在**的事件 action / 产物路径（`case-anomaly --event-action`） | 结论悬空，无法回溯到工具输出 |

`case-validate` 只校验**记录格式**；`decided_by=EXPERT` 的结论还需要人工确认记录，格式合法不构成科学验收。

## 4. 证据矩阵（报告首页）

把散落在事件与假设里的证据压成一张表——这是专家读报告时最先看的部分：

| 问题 | 关键证据 | 结论 | 限制 |
|---|---|---|---|
| 是否完整 mt genome | 覆盖剖面（零覆盖缺口、末端低覆盖） | 否 | 末端覆盖过低，阴性判别力弱 |
| 是否污染 | 同源检索 | 弱支持"无外源" | 无核参考，不能排除 NUMT |
| 是否真实缺 tRNA | 多工具 + 结构 + 同源 | 未检出（等级 B） | 工具灵敏度依赖 |
| 是否环化 | 跨接缝 reads | 无支持证据 | 低覆盖使阴性功效低，不足以证明"不相接" |

**规则**：矩阵中每一行的"结论"都必须能追溯到 §3 的结论对象；**限制**列不得为空——
填不出限制，通常说明该结论还没到可报告的程度。

## 5. 强制输出模板（AI 回答与报告统一结构）

任何最终回答/报告固定四段，顺序不变：

```text
1. 观察事实（Observed facts）
   - 只写工具/数据直接产出的事实：数值、坐标、命令、退出码、文件 hash
   - 每条事实标明来源（工具 + 命令 + 产物路径/事件 action）

2. 有证据支持的结论（Supported conclusions）
   - 每条：conclusion_type + claim + status + confidence + evidence_for/against + limitations
   - 每条必须引用第 1 段中的事实（没有对应事实的结论不得出现在这一段）

3. 无证据支持的声明（Unsupported claims）
   - 明确列出被排除/被拒绝的推断（"不能排除 NUMT""不足以证明不相接"）
   - 也列出本来想说但没有证据的假设 —— 这一段的目的是防止过度解释

4. 下一步最小实验（Next minimum experiment）
   - 最多 1–3 项，每项写明"能改变哪条结论/置信档位"，以及缺什么输入
   - 若已到停止条件（决策树 §3），写"建议停止"并给出理由
```

**禁止**：把第 3 段省略；用频率、"通常/可能属于"、单工具结果把第 3 段的内容挪进第 2 段。

## 6. 报告顺序与案例记录的映射

- 报告顺序：**事实 → 证据 → 结论 → 限制 → 下一步**（与 §5 一致）。
- 结论层是**派生物**，不是第二套数据格式：结论落在 `case.json` 的 `anomalies[]`
  （`id`/`claim`/`status`/`confidence`/可选 `reads_support`）+ `events.jsonl` 的事件链上，
  由 `tools/experience.py` 写入，不新建平行文件。
- **`case-report` 已实现本契约**：`python3 tools/experience.py case-report <dir>` 生成的 `case.md` 依次包含
  ① 观察事实（事件） ② **证据矩阵**（逐异常：claim/状态/置信/reads 支持/证据事件）
  ③ 有证据支持的结论 ④ **无证据支持的声明** ⑤ 下一步最小实验 ⑥ 案例级汇总
  ⑦ 假设与未测项 ⑧ 证据边界；没有逐异常记录时明写"尚无逐异常判定"，而不是静默省略。
  用例：`tests/test_review_round_lesson_domain_and_report.py::CaseReportLayerTests`。
- 任何结论在报告里都必须能反查到 case 记录（`--event-action` 指向真实事件）。
- 正常样本（`case_type = normal_validation_case`）也要按同一结构写结论（"确认无异常"本身是一条结论，
  且必须写明检查覆盖了哪些项、未覆盖哪些项）—— 见 `learning-policy.md` §2。

## 7. 与其他文件的关系

- 证据门槛 / `confidence` / `reads_support` / 阴性证据强度：`evidence-standard.md`
- 分流、优先级、停止条件、参考等级、最小证据集：`diagnostic-decision-tree.md`
- 注释判据：`annotation_quality.md`；坐标与顺序：`standard_gene_order.md`
- 工具参数与退出码：`tool-catalog.md`；实现约束：`developer-contract.md`

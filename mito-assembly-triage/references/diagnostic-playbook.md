# 动态诊断手册

诊断循环：`INTAKE → HYPOTHESIZE → CHOOSE_TEST → EXECUTE → UPDATE → DECIDE → VERIFY → LEARN`

> **适用边界**：针对用户报告的**单个异常**按需诊断，默认不重组装、不跑全部工具。
> 结构化事实链由 `tools/experience.py` 的 `case-*` 子命令承载（`case.json` / `events.jsonl` / `case.md`）。

## INTAKE（先收语境，再动手）

必须记录：

- 用户报告的**具体异常**（现象 + 位置 + 已知前提）；
- 类群（目/科，尽量到属）与遗传密码表；
- 输入类型与实际可用数据：`FASTA` / `GFF` / `GB` / `FASTQ` / `BAM` / `GFA`；
- 工具与版本、是否部分组装、是否多 contig；
- 用户允许的资源开销、联网与对外上传权限（默认**不允许**向公共库上传）。

同时**先检查已有产出**（`logs/`、既有 BAM/GB/候选目录），避免重跑与重复消耗。

## HYPOTHESIZE（至少留一个替代解释）

每个异常都要保留竞争假设与可能的反证。最低覆盖：

| 异常 | 至少考虑的替代解释 |
|---|---|
| 基因缺失 | 漏注释、组装断裂、错接、真实丢失、NUMT/污染 |
| 内部 stop / 移码 | 密码表、边界、碱基错误、真实生物例外 |
| tRNA 异常 | 命名退化、真实结构退化、工具灵敏度、错接 |
| rRNA 边界异常 | 参考边界本身有误、类群长度差异、组装嵌合 |
| 接缝/控制区异常 | 重复错接、真实长度异质性、未闭合 |
| 覆盖异常 | NUMT、污染、多倍型/异质性、组装错误、比对歧义 |
| 重排/顺序不同 | 参照物远缘、方向/旋转表示差异、真实重排 |

写入 `case.json` → `hypotheses[]`，每条含 `id` / `explanation` / `support` / `against` / `unknown`（schema 强约束）。

## CHOOSE_TEST（选能区分假设、且成本低的检查）

- 优先选**能改变假设排序**、成本低、可复现的检查；
- 即使某项检查不能改变排序，**若它是安全关键修复的验收条件，也不得省略**；
- 没有可用数据时不循环空跑；无法区分时直接进入 DECIDE 并写 `UNRESOLVED`；
- 先查工具目录（`tool-catalog.md`）确认输入要求与"不能证明什么"。

## EXECUTE

- 只调用现有脚本与经典工具；记录**命令、参数、版本、输入 hash、输出路径**（写入 `events.jsonl`）；
- **禁止直接覆盖原始序列**；候选另存新目录；
- 对外查询/上传样本序列必须先获得显式许可（如 `cox1_id.py --allow-public-upload`）。

## UPDATE（按假设分别记账）

对每条假设写出：支持证据、反证、证据来源、依赖关系、尚缺数据。

**独立性按产生机制判断，不按软件个数**：共用同一参考、或共用同一份初始组装的多个软件结果，
**不按独立证据重复计数**（见 `evidence-standard.md`）。

## DECIDE（分级下结论，不用一个状态掩盖局部）

三层状态，别混用：

| 层级 | 取值 | 落点 |
|---|---|---|
| 假设级 | `SUPPORTED` / `REFUTED` / `UNRESOLVED` / `NOT_TESTED` | `hypotheses[]` 的 `support`/`against`/`unknown`，以及事件 `impact`（如 `H1:against`） |
| 操作级 | `NO_CHANGE` / `ANNOTATION_CORRECTED` / `SEQUENCE_CORRECTED` / `STRUCTURE_CORRECTED` / `UNRESOLVED` | `modifications[]` 与事件记录 |
| 案例级 | `RESOLVED` / `NO_CHANGE` / `UNRESOLVED`（**schema 强约束**） | `case.json` → `decision.status` |

案例级还要给 `decision.confidence` ∈ `high` / `moderate` / `low` / `not_assessable`（评分标准见 `evidence-standard.md` §2）。

**停止条件**：无法区分时必须停止，清晰指出"最少还需什么证据"，不要为了收尾而猜测。

## VERIFY

- 修复必须通过**修复前定义的**验收检查（否则等于事后找理由）；
- 区分三类修改所需证据：注释校正 / 单碱基修改 / 结构连接；
- 修改后的 FASTA 与注释**分别保存**，原件不动；`annot_check.py --require-circular` 作为注释侧验收（见 SKILL.md G9）。

## LEARN

- 写入成功、失败与未解决案例（`case-report`），说明新案例是**支持还是挑战**既有 lesson；
- **AI 自己重复提出同一解释不构成可信度提升**（频次不是证据）；
- 政策边界见 `learning-policy.md`。

## 低成本检查速查（按需，不是固定流水线）

| 异常 | 先做（低成本） | 升级检查（需更多数据） |
|---|---|---|
| 基因缺失 | 名称规范化、全长同源检索、邻域与跨环位置复核 | 候选 contig、GFA、独立注释；必要时 reads |
| 内部 stop / 移码 | 按类群密码表重译原始坐标、同源蛋白比对 | 若怀疑碱基错误 → 含**竞争参考**的 reads pileup |
| tRNA 异常 | 反密码子、专用结构模型、近缘同源与邻域 | 分歧时复核序列、局部 reads（若可用）与文献 |
| rRNA 边界异常 | 同源保守区、相邻基因、参考注释准确性 | 必要时结构/转录证据；**不得只按长度改边界** |
| 接缝/控制区/重排 | 多条类群适当参考、端部唯一重叠、结构图 | 唯一锚定的接缝 reads、pair/long reads、竞争结构检查 |
| 覆盖异常 / 疑似 NUMT | MAPQ 分布、soft-clip、异常 mate、局部重复 | 竞争性比对到核与线粒体候选；**核组装缺失时须注明不能完全排除 NUMT** |

## 相关文件

- 证据等级与 `raw-read-supported` 要求：`evidence-standard.md`
- 注释侧硬阈值：`annotation_quality.md`
- 顺序与坐标约定：`standard_gene_order.md`
- 工具输入/产出与局限：`tool-catalog.md`；环境：`tool_check.md`
- 经验晋升与共享边界：`learning-policy.md`

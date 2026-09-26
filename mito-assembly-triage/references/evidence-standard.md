# 证据标准

> 核心原则：**按结论类型定证据门槛，而不是统一要求"两个来源"**。
> 配套：注释判据与工程预警见 `annotation_quality.md`；诊断记账见 `diagnostic-playbook.md`；
> 分流/优先级/停止条件/参考等级见 `diagnostic-decision-tree.md`；结论层与报告模板见 `conclusion-report.md`；
> schema、解析、退出码等**实现约束**见 `developer-contract.md`（不要把实现约束当生物学规则）。

## 1. 结论类型 → 最低证据

| 结论 | 最低记录/检查 | 不得越界宣称 |
|---|---|---|
| FASTA / 注释内部一致 | 原始文件 hash、坐标/链/密码表、重译、同源或结构证据 | ≠ 原始 reads 证明该碱基或结构正确 |
| 基因存在/身份/边界 | 对应的同源覆盖、ORF 或 tRNA 结构、邻域与替代定位；并审查参考可靠性 | 两个依赖**同一数据库**的预测不算两份独立真值 |
| 单碱基序列校正 | 原始 reads pileup：碱基质量、链向、重复、替代等位与比对歧义；尽量外部或不同技术独立复核 | 仅"深度高"或**只比对到 mt 候选**，不能排除 NUMT / 混样 |
| 新增接缝 / 环化 | 每个新增与闭合连接的坐标、唯一锚定、支持分子数、竞争路径、重叠裁剪、数据分辨力 | short reads 跨不过重复区时，不能声称唯一闭环 |
| 真实重排 / 丢失 | 排除漏注释与方向/坐标误差；多参考比较；并有可用的结构或原始数据支持 | 参考与样本不一致**不能单独**证明错误，也不能单独证明演化事件 |

## 2. 置信等级（逐条结论的 `confidence`）

`confidence` 必须**针对具体结论**，不能属于整个样本；禁止用一个 global `high` 覆盖案例中所有问题。
同一案例里可以同时存在 `high` 的注释身份结论与 `not_assessable` 的结构结论。

| 等级 | 判据 |
|---|---|
| `high` | 有**直接证据**，或**多种相互独立的机制证据**支持；**关键竞争解释**已被当前数据排除到该结论所需的程度，并留证 |
| `moderate` | 有多条相互独立的旁证（同源 + 结构 + 邻域等）一致，但缺少直接 reads 或存在未测替代解释 |
| `low` | 只有单类证据（或同机制的多软件重复），替代解释仍在 |
| `not_assessable` | 该结论所需的关键输入缺失（如无 FASTQ/BAM、无可用参考），无法评估 |

**`high` 不要求排除所有理论可能性**，只要求排除**会改变当前结论的主要替代解释**——
即"若这个解释成立，结论就要改写"的那几个候选；纯逻辑上可能存在但没有样本内证据支持的猜测不必逐一排除。
否则 cox1/rrnL 身份、tRNA 类型这类结构性结论会永远停在 `moderate`，`high` 形同虚设。
"排除了哪些替代解释、依据是什么"必须写进 `evidence_against` / 事件的 `impact`，"未被排除"的写进 `not_tested`。
这一条与下面 §2.1 一致：无核基因组限制的是**与 NUMT 相关**的结论，不限制与 NUMT 无关的注释身份结论。

### 2.1 上限按结论分别生效，不设全局上限

不要设置"无 reads 则全部结论 ≤ moderate"或"无核基因组则全部结论 ≤ moderate"这类**全局**上限。
它们应分别作用于相关结论：

| 缺失的输入 | 限制什么 | **不**限制什么 |
|---|---|---|
| 无 FASTQ/BAM | 样本碱基真实性与物理接缝/闭环的 reads 验证无法评估：`UNRESOLVED` / `not_assessable`，`reads_support = NOT_ASSESSED` | 基因身份、边界、tRNA/rRNA 结构等可由序列/同源/RNA 结构充分支持的注释结论 |
| 无核基因组 | NUMT 排除能力（声明受限，不得称已排除） | 与 NUMT 无关的注释结论 |

并且始终保持两个区分（见 §3）："有 reads 与候选一致" **不等于** "reads 足以排除其他候选结构"。

“候选与同源/图结构等间接证据相容”是另一条命题，可按已有证据给出相应置信度；
不能将它改写为“样本结构已验证”。文件中可直接观察的 N、断点或注释矛盾也仍可报告。
缺少该命题必要输入才用 `not_assessable`；已有可评估但不足或冲突的证据时用 `low/moderate`
并保留 `UNRESOLVED`，不要用一个状态覆盖所有结构相关陈述。

落地位置：每个异常各自把 `status` / `confidence` / `reads_support` 写入 `case.json` 的 `anomalies[]`
（schema 校验这三个枚举）；案例级 `decision` 只作汇总结论，不能代替逐异常判定。

## 3. `raw-read-supported` 的门槛与写法

声称 `raw-read-supported` 必须同时给出：

- 输入 hash（FASTQ/BAM）、参考候选（含其来源与版本）；
- 工具 + 参数、唯一比对策略（是否多映射过滤）；
- MAPQ 与碱基质量分布（不只报均值深度）；
- 链向偏好、证据坐标、**重复标记状态**（是否已 `markdup`/去重）。

这些字段现在由脚本实际产出，而不是只写在文档里：

| 结论类型 | 产出点 | 关键字段 |
|---|---|---|
| 单碱基校正 | `scripts/depth_analysis.py` 的 `base_support()` | `depth`、`support`、`strand_counts`、`mapq`、`excluded{duplicate,low_mapq,low_baseq,secondary}`；深度/MAPQ/碱基质量/链向四项全部通过才 `callable` |
| 接缝 / 环化 | `scripts/circularize.py` 的 `junction_evidence()` | `support`（独立 QNAME 数）、`strand_counts`、`mapq_min`/`mapq_median`、`duplicates_excluded` |

`reads_support` 只取三个值：

- `NOT_ASSESSED`：无 reads 或未做该检查（缺 reads **不妨碍**有充分序列/同源证据的注释纠错）；
- `READS_CONSISTENT`：有 reads 与候选一致；
- `READS_DISCRIMINATING`：reads 足以排除其他候选结构。

两句话必须分开写："**有 reads 与候选序列一致**"（`READS_CONSISTENT`）与
"**reads 足以排除其他候选结构**"（`READS_DISCRIMINATING`）。只有前者时，不得宣称后者。
还要标记 PCR / 光学重复与共享片段偏差。

**没有核基因组时，必须写明 NUMT 排除能力的限制。**
**"无内部终止"不能证明序列来自线粒体**：numt 可以不携带 in-frame 终止密码子（移码、整块缺失同样常见）。

### 3.1 `READS_DISCRIMINATING` 的前置条件：竞争参考

`READS_DISCRIMINATING`（reads 足以**排除其他候选结构**）必须声明竞争参考，因为"排除了谁"取决于"比过谁"：

| `competitive_reference` | 含义 | 允许的最高 `reads_support` |
|---|---|---|
| `none` | 只有单一线粒体候选参考（无核基因组，也无其他同源位点） | **最多 `READS_CONSISTENT`** |
| `nuclear` | 有核基因组/核组装参与竞争比对 | 可到 `READS_DISCRIMINATING`（须记录比对结果与参数） |
| `organelle+related` | 有同源细胞器/相关序列（其他位点、重复拷贝、近缘 paralog）参与竞争比对 | 可到 `READS_DISCRIMINATING`（须记录比对结果与参数） |

- 证据记录里要写明 `competitive_reference` 取值与**具体参考文件 + hash**；**未声明时按 `none` 处理**。
- **没有竞争参考就不得写 `READS_DISCRIMINATING`**：单 mt 参考下的高 MAPQ 不能排除 NUMT——
  真实案例中局部 4–5× 深度超额在去重、高 MAPQ、低错配条件下**仍然无法判定来源**，正是缺竞争参考所致。
- 有竞争参考、但比对结果模棱两可（两处都无显著差异）时，同样**不得**升级，只保留超额/一致性观测。

### 3.2 阴性证据的强度（"没有发现"不等于"不存在"）

阴性结果（无 reads 支持 / 无同源命中 / 无跨接缝 reads）能支持什么，取决于该检验的**功效**：

| 阴性观测 | 可否定当前候选结构？ | 依据 |
|---|---|---|
| 可唯一比对且按文库/读长预期有足够跨接支持的区域内无支持 | **较强**（可以据此降档或否定） | 明确该检验在候选成立时的期望支持；平均深度本身不能保证阴性功效 |
| 低覆盖区（如 <5×，局部末端 <1×）内无支持 | **弱** | 期望命中数本就很少，"没有"在两种假设下都常见 |
| 重复/低复杂度区内无支持 | **弱** | 多映射被过滤，命中率与环境无关地下降 |
| 短读长无法跨越的区间（如远长于读长的重复、接缝闭合） | **不可**作为否定依据 | 现有数据原理上不具备分辨力 |
| 单一工具未发现（如只跑了 MITOS2 没找到某 tRNA） | **弱**（只记"未检出"） | 工具灵敏度与命名退化（见 `annotation_quality.md` §3） |

写法：阴性结论必须同时给出**覆盖/深度、读长、多映射处理、工具版本**这些功效条件，
并写成"**没有证据支持 X**"，而不是"**证明 Y**"。真实案例中"两端无跨接缝 read pair"曾被当成
"两端不相接"的证据，实际末端只有 0.09–0.15× 覆盖，该检验**功效极低**。

适用范围：本节适用于**所有工具产生的阴性结果**（reads、同源检索、注释工具、接缝检索），不只是 reads。

## 4. "独立性"的定义

- 按**产生机制**判断独立性，不按软件个数：
  - BLAST 与另一个从**同一参考库**训练/推断的注释 → 不是两份独立真值；
  - 两个比对软件跑**同一份 BAM 输入**的同一区域 → 不构成独立复核。
- 真正的独立来源示例：样本 reads、独立实验技术、结构/转录证据、不同数据库且各自人工审校的注释。

## 5. 状态与失败语义

- **只有当某条具体结论缺少它自己所需的证据时**，才把该结论标为 `UNRESOLVED` / `not_assessable`。
  不要因为"没有 reads"或"没有核基因组"就把**整个案例**判成 `UNRESOLVED`：
  无 reads 不妨碍注释身份/边界类结论在充分同源、翻译或 RNA 结构证据下定为 `RESOLVED`，
  只是该结论的 `reads_support = NOT_ASSESSED`。
- 证据冲突、重复无法唯一解析或缺必要输入时，该命题保持 `UNRESOLVED`；
  置信度按 §2 的证据档位给出，无法评估时才用 `not_assessable`；
- 修改与原件分离保存；候选文件**不得覆盖**原始文件；
- 参考相似度、基因顺序、单软件结果只能作**定位线索**。

## 6. 来源与引用管理

判据分三层，措辞必须对上：**(a) 一般生物学预期**（需来源 + 覆盖类群）、
**(b) NCBI 提交审查要求**（需引用官方页面，不得把工具默认值写成"NCBI 要求…"）、
**(c) 本工具的工程阈值**（需记录选择依据与被测范围，并标为启发式）。

硬规则（有生物学含义的判据）应附：来源 URL/DOI、覆盖类群、更新时间、已知例外。
**工具参数的默认值是工程启发式，不得伪装成领域公理**（例如 `>8 bp` 重叠阈值、tRNA 长度区间、`9+/4−` 链分布）。

核心权威资料（知识来源，不是硬阈值）：

| 主题 | 来源 |
|---|---|
| NCBI 细胞器基因组提交与注释要求 | https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/ |
| NCBI 遗传密码表（按类群选表） | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/ |
| MITOS / Bernt et al. 2013（历史注释的系统性误差） | https://pubmed.ncbi.nlm.nih.gov/22982435/ |
| 基因边界不确定性与不完整终止 | Donath et al. 2019, *NAR* 47(20):10543–10552, DOI 10.1093/nar/gkz833（= PMC6847864） |
| 昆虫线粒体测序/注释、基因排列与 tRNA 结构变异 | Cameron 2014, *Syst. Entomol.* 39:400–411, DOI 10.1111/syen.12071 |
| 非典型线粒体 tRNA 结构与注释（Jühling et al. 2012） | https://pmc.ncbi.nlm.nih.gov/articles/PMC3326299/ |
| metazoa 线粒体 tRNA 数量/结构变异 | https://pmc.ncbi.nlm.nih.gov/articles/PMC11571959/ |
| GetOrganelle 方法与组装图相关限制 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7488116/ |
| MitoHiFi 方法（长读长语境，概念参考，非必需依赖） | DOI 10.1186/s12859-023-05385-y（= PMC10354987） |

引用时写明**该来源覆盖的类群**；跨类群外推时必须标注为外推。

## 7. 实现与格式限制

案例 schema、双路径验证和报告生成器的限制见
[developer-contract.md](developer-contract.md) §1/§7；这些实现要求不构成生物学证据。

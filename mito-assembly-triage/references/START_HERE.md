# 从这里开始（入口层）

> 这是一个**问题导向**的入口：先用「你有什么数据 / 你被什么现象叫来」定位路线，再按需读其他 references。
> 不要一上来就通读全部文件，也不要默认跑全套工具。科学判据在 `evidence-standard.md` /
> `annotation_quality.md`；**在"组装问题"和"注释问题"之间怎么分流**见 `diagnostic-decision-tree.md`；
> 报告怎么组织见 `conclusion-report.md`。

## 1. 我手里有什么数据？→ 能问到哪一层

| 你有的数据 | 入口动作 | **能**回答 | **不能**回答（证据上限） |
|---|---|---|---|
| 只有 FASTA | `scripts/seq_stats.py <assembly.fasta>` 看长度/GC/模糊碱基 | 序列基本体检：长度、contig 数、GC、`N`/模糊碱基位置 | 不能判断碱基是否正确、不能判断是否完整、不能判断是否环化 |
| FASTA + GB/GFF | 先 `mitos2_to_genbank.py`（MITOS2 输出）或直接用 GB → `annot_check.py` | 注释内部一致性：基因数量/身份、CDS 翻译与起始终止、tRNA/rRNA、重叠、链分布、（可选）顺序对照 | 注释自洽 **≠** 序列正确；**≠** reads 支持；`--require-circular` 只校验声明 |
| FASTA + FASTQ（或已有 BAM） | 比对 → 排序去重 → `scripts/depth_analysis.py <sample.bam> <genome.fasta>` | 覆盖剖面、零覆盖缺口、末端低覆盖、soft-clip、mate 分布、逐碱基支持（深度/MAPQ/碱基质量/链向） | 单参考比对**不能排除 NUMT**（除非有竞争参考，见 `evidence-standard.md` §3）；低覆盖 **≠** 必然错接 |
| FASTA + FASTQ + 核基因组/核组装 | 同上 + **竞争性比对**（mt 候选 vs 核候选） | 可以把 `READS_DISCRIMINATING` 用在 NUMT 相关结论上 | 没有这一步时，NUMT 相关结论只能是 `UNRESOLVED` |
| 需要一个参考 GenBank | 按参考等级选一条**已注释**公共记录（§5.1）；登记 accession/版本/日期/类群/hash | 顺序、边界、基因身份的对照 | 参考不是真值；**上传样本**才是需要单独授权的动作 |
| 多个 contig（含 GFA） | `circularize.py`（两 scaffold 场景）+ 端部唯一性检查 | 4 种拼接组合中的唯一最高分方向、新增接缝的 `junction_evidence()` | **单一**接缝证据不能宣称最终闭环；单 contig 的闭环需 terminal overlap + 跨接缝 reads（见 `standard_gene_order.md` §3） |
| 需要一个近缘参考 | 按 `diagnostic-decision-tree.md` §5 的**参考等级**选；本地库优先；需新取时按 §5.1（**下载公共参考不需要上传授权**，但必须登记 accession/版本/hash） | 顺序/边界/身份的对照线索 | 远缘参考不能用于边界判断；参考注释自身也可疑 |

**先查已有产出**：`logs/`、既有 BAM/GB/候选目录。重跑一遍不增加证据，只增加噪声。

## 2. 我是被哪个现象叫来的？→ 问题导向入口

| 用户的说法 | 第一步（低成本） | 必读 | 必须先想到的竞争解释 |
|---|---|---|---|
| "少了一个基因 / 少 tRNA" | 名称规范化 → 全长同源检索 → 邻域与跨原点位置复核；再看 reads 覆盖 | `annotation_quality.md` §2/§3、`diagnostic-decision-tree.md` §1 | 漏注释、组装断裂、错接、真实丢失、NUMT/污染；tRNA 还要过 §3 的"搜索失败等级" |
| "内部 stop / 移码" | 按类群密码表重译原始坐标 → 同源蛋白比对 → 查 `/transl_except` | `annotation_quality.md` §2、`SKILL.md` 的 `transl_except` 段 | 密码表、边界与 `codon_start`、碱基错误、真实生物例外（须 registry 证据） |
| "环不起来 / 两端接不上" | 端部是否唯一（重复？）→ 有无跨接缝 reads → 是否只是拓扑声明 | `standard_gene_order.md` §3、`tool-catalog.md` 的 `circularize` 行 | 端部重复导致无法唯一锚定、数据分辨力不足（低覆盖阴性）、真实未闭合 |
| "深度异常 / 疑似 NUMT" | MAPQ 分布、soft-clip、mate 异常 → 局部重复检查 → **竞争参考**比对 | `evidence-standard.md` §3、`diagnostic-decision-tree.md` §1 | NUMT、污染、异质性/多倍型、组装错误、比对歧义；**来源归属**往往留 `UNRESOLVED` |
| "基因顺序不一样" | 归一旋转与整链反向互补 → 换类群适当的参考 → 顺序差异与基因集差异分开报告 | `standard_gene_order.md` §1/§3/§6 | 参考远缘、方向/旋转表示差异、真实重排 |
| "一片多余的 N / 模糊碱基" | `seq_stats.py` 定位 → reads 覆盖与 pileup → 该区是否低复杂度/重复 | `evidence-standard.md` §1、`SKILL.md` 的单碱基校正段 | 组装空缺、覆盖不足、比对歧义、真实非 A/C/G/T |
| "cox1 是不是这个物种" | `cox1_id.py --allow-public-upload --coords ...`（需显式同意） | `tool-catalog.md` 的 COX1 行 | identity/coverage **不等于**物种鉴定；多 HSP 不可汇总时必须放弃单值 |

现象不在表里时：先写清"现象 + 位置 + 已知前提 + 可用数据"，再进 §3 的首次诊断。

## 3. 首次诊断（约 30 分钟，不必跑全套）

```bash
# Step 0 环境：先确认工具/数据库就位，避免把"环境缺失"当成"数据问题"
bash scripts/check_env.sh

# Step 1 序列体检（只读）
python3 scripts/seq_stats.py <assembly.fasta>

# Step 2 注释质检（有 GB 时）
python3 scripts/annot_check.py <ann.gb> --table 5 --taxon '<目/科>'   # 有近缘参考再加 --ref <ref.gb>
#   退出码 0 = 无发现；2 = 有待核查项（逐条入账，不等于通过）；1 = 有错误

# Step 3 reads 裁判（有 FASTQ/BAM 时；BAM 需已排序+建索引）
python3 scripts/depth_analysis.py <sample.bam> <genome.fasta>
#   ⚠ 单参考比对不能排除 NUMT；要区分来源必须有竞争参考

# Step 4 记账：把现象、假设与逐异常判定写进结构化案例
python3 tools/experience.py case-init work/case-001 --issue <类型> --observation '<现象>' \
  --hypothesis 'H1: <解释一>' --hypothesis 'H2: <解释二>'   # 已核实无异常则加 --case-type normal_validation_case
python3 tools/experience.py case-event work/case-001 --action annot_check \
  --result '<结论>' --impact H1:against
python3 tools/experience.py case-anomaly work/case-001 --id A1 --claim '<claim>' \
  --status UNRESOLVED --confidence low --event-action annot_check

# Step 5 决定下一步：走 decision tree，而不是继续堆检查
```

Step 5 之后只需回答一个问题：**当前异常更像组装问题还是注释问题**（`diagnostic-decision-tree.md` §1）。
分流错了，后面所有检查都不会改变结论。若两个方向都说得通 → 保持 `UNRESOLVED` 并写清"最少还需什么证据"。

## 4. 六个最容易踩的坑（每一条都在真实案例中被踩过）

1. `annot_check.py` **退出码 2 ≠ 通过**：它是"无 ERROR 但有待核查"，每条必须逐条入账。
2. `READS_CONSISTENT` **≠** `READS_DISCRIMINATING`：前者只是"reads 与候选一致"，后者才是"reads 排除了其他结构"。
3. `case-validate` **只校验记录格式**（schema/字段/枚举），**不证明科学结论**。
4. MITOS2 的 **circular 模式与 `topology=circular` 头都不是物理环化证据**。
5. **没有核基因组/竞争参考时**，NUMT 相关结论只能是受限的 `moderate`/`UNRESOLVED`，不得写成"已排除"。
6. **单个案例不构成规则**：一个样本的观察不能升级为类群规律（见 `learning-policy.md` §2/§3）。

## 5. 与其他文件的关系

- 分流与优先级、停止条件、参考等级、最小证据集：`diagnostic-decision-tree.md`
- 证据门槛、`confidence`、`reads_support`、阴性证据强度：`evidence-standard.md`
- 注释判据（CDS/tRNA/rRNA/重叠/链/例外）：`annotation_quality.md`
- 报告与结论层、输出模板、证据矩阵：`conclusion-report.md`
- 工具输入/产出/退出码：`tool-catalog.md`；环境：`tool_check.md`
- 案例/经验/共享边界：`learning-policy.md`
- 实现与格式约束（schema、selector 归一化、原子写入、测试对应）：`developer-contract.md`

---
name: mito-assembly-triage
description: 线粒体基因组组装质检、诊断与修复流水线。用于处理 NOVOPlasty/GetOrganelle/SPAdes 组装结果的基因方向错误、注释失败、模糊碱基、控制区异常等问题。执行 reads 验证、基因顺序重建、环化修复、多工具融合注释、论文级环形基因图绘制，并在每次任务后沉淀案例到外部 `MITO_KNOWLEDGE_DIR` 实现 skill 自我进化。当用户提到"线粒体组装有问题/注释不出来/方向不对/基因顺序异常"或需要从 reads 重新组装线粒体时使用。
---

# 线粒体组装诊断与修复流水线（进化式）

## V2 动态诊断入口

本 Skill 按需处理用户报告的单个异常，不默认重组装或运行全部工具。诊断循环为：
`INTAKE → HYPOTHESIZE → CHOOSE_TEST → EXECUTE → UPDATE → DECIDE → VERIFY → LEARN`。
先读取 `references/diagnostic-playbook.md` 和 `references/evidence-standard.md`，再按异常选择 V1 脚本；工具目录见 `references/tool-catalog.md`，经验边界见 `references/learning-policy.md`。

每个异常只能判为 `RESOLVED`、`NO_CHANGE` 或 `UNRESOLVED`，另记录 `high/moderate/low/not_assessable` 置信等级。没有 reads 时不得报告 raw-read-supported；参考和单软件结果不能代替样本证据。

结构化案例可用：
```bash
python3 scripts/v2_case.py init work/case-001 --issue internal_stop \
  --observation 'nad5 出现内部 stop' --input assembly_fasta assembly.fasta
python3 scripts/v2_case.py event work/case-001 --action annot_check \
  --result 'table 5 下仍有内部 stop' --impact H1:against
python3 scripts/v2_case.py validate work/case-001
python3 scripts/v2_case.py report work/case-001
```
案例目录应位于工作目录或外部知识目录，不覆盖输入；`case.json`、`events.jsonl` 和 `case.md` 是同一事实链的结构化/可读表示。

## 设计哲学

- **稳定层与知识层分离**：SKILL.md 只含稳定的主流程；经验数据写入 `MITO_KNOWLEDGE_DIR`（默认位于用户 XDG 数据目录），不写入 skill 安装目录
- **每次任务强制"结束仪式"**：处理完一个样品且用户确认持久化后，写结构化案例 -> skill 越用越聪明
- **先检索后执行**：新任务开始前先查外部案例库，有相似案例先读历史做法
- **准确度优先于速度（最高原则）**：
  - 耗时不是问题，产出必须经得起推敲。任何"图快"的捷径（如用参考序列拼接缺失片段、用宽松参数跳过验证）都违背本 skill 宗旨
  - 每个修复动作必须能回答：**这个序列是 reads 支持的，还是我"借用"的？** 参考只能定位候选差异，不能凭覆盖度或相似度把参考序列引入样本；每个位点记录 pileup、碱基质量、MAPQ、链偏好和竞争候选。
  - 长任务（GetOrganelle 数小时、全量比对 1-2h、metaSPAdes 组装）用后台化策略等待，不因耗时偷工减料
  - 凡是要下结论的节点，必须过 Quality Gates
- **经验不足时搜索文献（辅证原则）**：
  - 当 `experience.py search` 无匹配案例、或现有案例无法解释新信号时，**必须搜索文献/资料**（文献检索、资料检索、工具官方文档）作为判断依据
  - 文献结论与经验冲突时，以文献为准并记录到案例（经验是样本性的，文献是普遍性的）
  - 不要仅凭直觉或"上次这样做的"就下结论
- **不要重复造轮子（工具优先原则）**：
  - **优先使用专业生物信息工具**，不手写脚本替代现成工具：
    - 注释：MITOS2（金标准）、MitoFinder、arwen（tRNA）、tRNAscan-SE、NCBI 工具
    - 组装/修复：GetOrganelle、metaSPAdes（优于 SPAdes）、NOVOPlasty、MitoZ
    - 比对/验证：blast、minimap2、bwa、samtools
  - 反面教训（历史案例）：手写"蛋白读码框法"精修 CDS 边界导致 ND2 被切短一半（160 vs 正确 333aa），而 MITOS2 直接给出正确结果
  - **何时才写自定义代码**：仅当现成工具无法解决时，针对特定问题写个性化高级代码作为补充——且必须验证其输出（如 ND5 homopolymer 区 reads 组装、蛋白交叉验证，是 MITOS2 盲区的补充）
  - 自定义代码的产出要能与专业工具交叉验证，不能替代专业工具成为主流程

## 任务开始前

### Step 0: 检索历史经验
```bash
python3 tools/experience.py search --query "<样品信号关键词>"
# 例: python3 tools/experience.py search --query "模糊碱基"
# 例: python3 tools/experience.py search --query "反向块 覆盖度"
```
有相似案例 -> 先读 `$MITO_KNOWLEDGE_DIR/cases/*.md`，复用其判定和动作。
**无匹配案例或信号陌生 -> 搜索文献/资料作为辅证**（文献检索 / 资料检索），
将文献结论写入案例的"新知识"栏——文献是普遍性证据，案例是样本性证据。

### Step 0b: 环境检查
```bash
bash scripts/check_env.sh
```
检查 blastn/minimap2/bwa/samtools/GetOrganelle/MitoFinder/MITOS2 及数据库路径。
**机器路径通过环境变量或 `config/env.sh` 配置（不要把路径写回 skill）**；工具细节与已知坑见 `references/tool_check.md`。

### Step 0c: 假设登记（借鉴 omics-skills，探索前先列候选解释）

诊断开始前，显式列出**至少 3 个候选解释**，后续每步证据用于支持/排除；记录到案例中：

| 候选假设 | 预期信号 | 排除/确认证据 |
|----------|----------|---------------|
| A. 组装 misassembly（重复区错接） | 反向块 + 连接区覆盖度塌陷 | reads mate 分布/覆盖度剖面 |
| B. 参考物种错误 | 与参考相似度 <90% | COX1 鉴定 |
| C. 后处理伪影（模糊碱基等） | CDS 内 R/Y/M/W/N | reads 支持单一碱基 |
| D. 真实生物学变异（基因重排/CR 异质性） | 顺序与近缘种不同但 reads 支持 | 多物种对比 + reads 跨边界 |
| E. 组装不完整（缺 CR/片段化） | 长度偏短 + 末端低覆盖 | 比对确认缺失区域 |

> 规则：**假设被证据排除即划掉；被确认即定论**。不要带着未排除的假设进入修复阶段。

## 输入要求 (Input Requirements)

- 必需：组装 fasta（NOVOPlasty/GetOrganelle/SPAdes 输出）
- 强烈建议：原始 reads（R1/R2 fastq.gz）——无 reads 则跳过 Step⑤，修复后验证降级为仅参考比对
- 可选：参考/近缘物种 GenBank 文件（只作定位和比较，不作为样本序列来源）
- 工具：blastn, minimap2, bwa, samtools, GetOrganelle, MitoFinder, MITOS2（见 check_env.sh）

## Quick Reference（快速定位）

| 任务 | 动作 |
|------|------|
| 快速体检序列 | `python3 scripts/seq_stats.py <assembly.fasta>` |
| 重建基因顺序 | `python3 scripts/blast_genes.py <ref.gb> <target.fna>` |
| 鉴定物种 | `python3 scripts/cox1_id.py <genome.fasta>` |
| reads 验证组装 | bwa mem -> `python3 scripts/depth_analysis.py <bam> <fasta>` |
| 生成环化候选 | `python3 scripts/circularize.py <s1> <s2> <ref.gb> <outdir> --bam <bam> --junction-region <chr:start-end> --reads-validated`（默认只输出 candidate；人工核对全部接缝后才可加 `--accept-candidate`） |
| 后台跑长任务 | `bash scripts/run_bg.sh <名> -- <命令> [参数...]` + `bash scripts/check_bg.sh <名>` |
| 检索历史案例 | `python3 tools/experience.py search --query "..."` |
| 沉淀本次经验 | `python3 tools/experience.py add-case --sample ...` |

## 长任务后台化（防 AI 前端超时的运行策略）

生物信息步骤耗时差异大，**超过 5 分钟的任务必须后台运行**，禁止在前端阻塞等待：

```bash
# 启动后台任务 (自动写 logs/<名称>.log/.pid/.status)
bash scripts/run_bg.sh <任务名> -- <命令> [参数...]
bash scripts/run_bg.sh <任务名> --trusted-shell "<已确认可信的完整命令>"
# 轮询状态 (运行中/完成/失败 + 日志尾部)
bash scripts/check_bg.sh <任务名> [--tail 20]
```

默认模式按参数逐个执行，不解析 shell 语法；只有确认命令和其中的路径均可信时，才使用 `--trusted-shell` 执行管道或重定向。

步骤耗时分级与建议运行方式：

| 步骤 | 典型耗时 | 运行方式 |
|------|----------|----------|
| ① seq_stats | <1s | 直接跑 |
| ② minimap2/blastn 比对 | 秒~分钟 | 数据小直接跑；大参考建议后台 |
| ③ blast_genes 逐基因 | 秒~分钟 | 直接跑 |
| ④ cox1_id (NCBI 轮询) | 1-5 分钟 | 自带等待逻辑，直接跑 |
| ⑤ bwa 全 reads 比对 | **5-30 分钟** | ! **必须后台** |
| ⑥ GetOrganelle 重组装 | **10-60 分钟** | ! **必须后台** |
| ⑥ circularize 外科拼接 | 秒~分钟 | 直接跑 |
| ⑦ MitoFinder/MITOS2 注释 | 2-15 分钟 | 建议后台 |
| ⑧ 环形图 | <1 分钟 | 直接跑 |

**轮询策略（避免空转）**：后台任务启动后先返回用户/做其他步骤；检查状态时用 `check_bg.sh`，状态为 running 时**不要连续空转等待**——先继续可并行的工作（如参考下载、文献检索），隔几分钟再查。

**宿主工具增强**：若支持后台监控，可用日志监控或完成通知；run_bg.sh 方案作为跨环境默认。

## 主流程（8 步决策流水线）

### ① 统计体检
```bash
python3 scripts/seq_stats.py <assembly.fasta>
```
检查：长度、GC%、AT%、模糊碱基（R/Y/M/W/N 等）及位置、header 格式。
- **模糊碱基在 CDS 内 -> 记录为待修复伪影**（NOVOPlasty 正常只出 ACGT，出现即是异常信号）

### ② 参考比对（建立参照系）
```bash
# 下载近缘物种注释 GB（物种未知时先跳过，用 Step④ 鉴定后补）
# 全基因组比对读方向:
minimap2 -x asm5 <ref.fna> <contig.fna> > mm.paf   # 看链方向和断点
blastn -query <contig.fna> -db <refdb> -outfmt "6 qseqid pident length qstart qend sstart send"
# sstart > send 表示反向; 多个方向矛盾块 = 结构异常
```
- **相似度 <90% -> 提示参考可能不合适**，转 Step④ 鉴定物种；参考差异本身不是样本错误。

### ③ 逐基因定位（降维重建基因顺序）
```bash
# 从参考 GB 提取 37 基因 -> 逐个 blastn-short 到 contig -> 按坐标排序重建顺序
python3 scripts/blast_genes.py <ref.gb> <contig.fna> [--out out.txt]
```
对照 `references/standard_gene_order.md` 检查：
- 反向块（连续多个基因方向颠倒）-> misassembly 信号
- 缺失/多余基因 -> 记录

### ④ 物种鉴定（COX1 条形码）
```bash
python3 scripts/cox1_id.py <contig.fna> --allow-public-upload [--coords 1,2]   # 默认自动找 COX1
```
该命令会把选定片段上传到公共 NCBI；仅在确认数据可公开后使用。该命令会将选定片段上传到公共 NCBI；仅在确认数据可公开后执行。结果解读：
- COX1 结果必须同时报告 identity、alignment coverage、多个近似命中和数据库版本；固定序列模式不能定位 COX1，也不自动下确定物种结论。
- 拿到近缘种后回到 Step② 换参考重新比对

### ⑤ reads 裁判（组装真伪验证，reads 可用时必做）
```bash
# ! 长任务: 全 reads 比对建议后台 (参考上方后台化策略)
bash scripts/run_bg.sh bwa_map --trusted-shell "bwa mem <contig.fna> R1.fastq.gz R2.fastq.gz | samtools sort -o c.bam && samtools index c.bam"
bash scripts/check_bg.sh bwa_map   # 完成后:
python3 scripts/depth_analysis.py c.bam <contig.fna> [--allow-base-replacement]
```
三个证据层次：
1. **覆盖度剖面**：<0.5×均值且持续 >500bp 的区域 = 可疑连接点
2. **mate 分布**：连接区 reads 的 mate 位置散乱 = 错接
3. **soft-clip**：CR 区大量 clip = 长度异质性（**正常**，非错误）；编码区 clip = 缺失序列
4. **模糊碱基 reads 支持**：清晰支持单一碱基 -> 伪影，可修复

### ⑥ 修复（两条路）
```bash
# 路线A: GetOrganelle 重组装（reads 可用时）— ! 长任务, 必须后台
bash scripts/run_bg.sh getorganelle -- get_organelle_from_reads.py -1 R1.fq -2 R2.fq -s <seed.fasta> -o GO_out -F animal_mt -R 15
bash scripts/check_bg.sh getorganelle

# 路线B: 外科手术拼接（重组装卡在重复区时）— 快, 直接跑
python3 scripts/circularize.py <scaffold1.fasta> <scaffold2.fasta> <ref.gb> <outdir> --bam <candidate.bam> --junction-region <chr:start-end> --reads-validated
```

未完成 reads 回贴和接缝验证时不要添加 `--reads-validated`；脚本会拒绝输出最终环状序列。

**准确度要求（ND5 类基因的教训）**：
- 任何"参考拼接"得到的基因片段（如参考[:372]+reads片段），**必须验证该片段的 reads 覆盖**：
  ```bash
  # 全量 reads 比对 (后台) -> 检查该片段覆盖度
  bash scripts/run_bg.sh bwa_all --trusted-shell "bwa mem -t 16 genome.fasta R1.fastq.gz R2.fastq.gz | samtools sort -o all.bam && samtools index all.bam"
  samtools depth -r <chr>:<start>-<end> all.bam | awk '{s+=$3;n++; if($3<100) low++} END {print "mean:", s/n, "低覆盖位置:", low}'
  ```
  低覆盖、低 MAPQ、链偏好或存在竞争候选 -> **不得接受参考拼接**，改用：
  1. 全量 reads 比对 -> 提取基因区域（含侧翼）reads
  2. **metaSPAdes**（优于 SPAdes, Allio et al. 2020）组装
  3. 6 读码框 × blastp 定位干净 ORF -> 精修起始终止 -> 纯 reads 基因
- 提取 reads 区域要**够宽**（覆盖完整基因 + 侧翼 ±200bp），否则缺失端段
**双重验证（必做）**：
1. 与近缘种参考：用于发现候选差异；基因顺序不同应标记 REVIEW，真实重排不得自动失败
2. reads 回贴：编码区覆盖度均匀（CR 区异质性除外）

### ⑦ 注释（工具优先，勿造轮子）
```bash
# MITOS2 金标准 — 统一入口 (python/refdir/refseqver/cmsearch PATH 全从 config/env.sh 读取)
bash scripts/run_mitos2.sh -i <genome.fasta> -o mitos2_out -c 5
# 等价于: $MITOS2_PY -m mitos.scripts.runmitos -i ... --refdir $MITOS2_REFDIR/ --refseqver $MITOS2_REFSEQVER
# 缺 cmsearch/RNAplot 等时检查 config/env.sh 的 MITOS2_EXTRA_PATH; 重装修复见 references/tool_check.md
```
**工具优先流程**：
1. 首选 **MITOS2** 完整注释（CDS/tRNA/rRNA/起始终止密码子）
2. 交叉验证：MitoFinder、arwen（tRNA）、参考锚定比对
3. **仅当 MITOS2 对个别基因失败**（如 homopolymer 区移码导致内部终止）时，
   才用 reads 区域组装 + 蛋白验证补充（见 ⑥ 的 ND5 方法）——这是补充，不是替代
4. 所有 CDS 必须**翻译验证**（无内部终止）+ 蛋白 blastp 交叉验证
- 注意：MITOS2 重装后需修复 drawmitos wrapper 和依赖（见 references/tool_check.md）

### ⑧ 绘图
```bash
# 论文级环形基因图 (正负链双环 + 呼吸链复合体配色 + GC环)
# 绘图 python 从 config/env.sh 的 PLOT_PY 读取 (需 BioPython+matplotlib)
bash scripts/run_circular_map.sh <final.gb> <out.png> --title "<物种> mitochondrial genome"
```
输出 300dpi PNG + SVG 矢量。配色按呼吸链复合体（蓝=ND、红=COX、橙=ATP、绿=CytB、金=rRNA、灰=tRNA）。
（注：如报 "需要 BioPython + matplotlib"，改 config/env.sh 的 PLOT_PY 指向含这两个库的 python）

## Quality Gates（质量门——借鉴 omics-skills，每个关键节点强制检查）

通过所有 gate 才可进入下一阶段；失败必须记录原因和处置。

| # | Gate | 通过标准 | 失败处置 |
|---|------|----------|----------|
| G1 | 输入完整 | fasta 可读、长度合理（线粒体 14-20kb）、有 header | 修复文件/重新导出 |
| G2 | 统计体检 | 模糊碱基已定位并归类（伪影 vs 真实） | 记录待修复 |
| G3 | 参考匹配 | 相似度 ≥90% 或已确认参考为近缘种 | COX1 鉴定换参考 |
| G4 | 基因完整性 | 37 基因定位结果、重复候选和缺失均可解释；顺序差异为 REVIEW | 记录缺失/反向/可能重排 |
| G5 | reads 支持（如有 reads） | 编码区覆盖度均匀，无 unexplained 低覆盖区 | 标记可疑连接点 |
| G6 | 修复验证 | 样本 reads 支持每个修改位点；参考比对仅作辅助 | 换修复路线 |
| G7 | 注释翻译 | 13 CDS 无内部终止，氨基酸长度接近参考 | 精修边界 |
| G8 | 环化 | 每个新增接缝均有独立 reads/pair/组装图证据，端部重叠已去冗余 | 输出 PUTATIVE/UNRESOLVED，不强制环化 |
| G9 | 注释质检 | `annot_check.py <gb> --require-circular` 全检通过（基因身份/翻译/起止/tRNA/rRNA/重叠/链分布）；真实例外需显式确认 | 修复坐标/确认重叠 |
| G10 | NCBI 预检 | 可选：table2asn 本地验证无内部终止等（NCBI organelle 提交前必做） | 按报错修复 |

> 全部通过后执行结束仪式；任何 gate 失败都应在案例的 actions 中记录处置方式（这正是知识进化的原料）。

## 结束仪式（强制！skill 进化的核心）

任务结束**在用户确认持久化且数据可存储时**执行：

```bash
export MITO_KNOWLEDGE_DIR="${MITO_KNOWLEDGE_DIR:-$HOME/.local/share/mito-assembly-triage}"
# 1. 写入外部知识库
python3 tools/experience.py add-case \
  --sample <样品名> --species "<鉴定结果>" --tools "<工具链>" \
  --signals "<信号列表, 逗号分隔>" --diagnosis "<判定>" \
  --actions "<采取的动作>" --lessons "<新学到的东西>"

# 2. 新信号/新坑增量更新（如有）
python3 tools/experience.py update-signals --signal "<新信号>" --judgment "<判定>" --ref "<案例>"
python3 tools/experience.py update-pitfalls --pitfall "<新坑>" --fix "<解决办法>"

# 3. 每积累 5-10 个案例，运行模式提炼
python3 tools/experience.py suggest-promotions
# -> 重复≥3次的手工操作固化为 scripts/; 重复信号提升为主流程检查项
```

## 内置判断标准（从历史案例沉淀）

| 信号 | 判定 | 处理 |
|------|------|------|
| NOVOPlasty 输出含 R/Y/M/W/N | 后处理伪影 | reads 验证后替换 |
| 与参考相似度 75~85% | 参考物种不对 | COX1 鉴定换参考 |
| 全基因组 blast 反向块 | 组装 misassembly | 重组装/外科修复 |
| 覆盖度持续骤降区 | 可疑连接点 | mate 分布确认 |
| CR 区低覆盖+soft-clip | **正常**（长度异质性） | 统计拷贝数取优势型 |
| CDS 翻译内部终止 | 边界错误 | 参考锚定+搜索起始终止 |
| CDS 3' 端不完整 (T/TA) | **正常** | 标注说明 |
| trnP 在 ND6 前/后 | **都可**（种间可变） | 以近缘种为准 |

## 输出 (Output)

- 修复后的环状基因组 fasta（trnI 起点，15-17kb）
- 融合注释 GB（13 CDS + 22 tRNA + 2 rRNA，翻译验证通过）
- 论文级环形基因图（PNG 300dpi + SVG）
- 诊断报告（triage 记录：信号、假设排除、判定、动作、gates 状态）
- 外部 `MITO_KNOWLEDGE_DIR` 案例 + signals/pitfalls/stats 更新（进化沉淀）

## 常见工具坑速查（详见 references/tool_check.md）

- MITOS2 本地报 `no such directory` -> 需 `--refdir <目录>/ --refseqver refseq89m`
- MITOS2 报 `cmsearch`/`plotprot.R`/`RNAplot`/`drawmitos` 找不到 -> 按 `references/tool_check.md` 配 PATH 或 wrapper
- MitoFinder 报 `install.sh.ok`/`Mitofinder.config` 缺失 -> 手动创建/复制，路径加尾斜杠

# 现有工具箱（按需调用）

> 运行前先做两件事：**① 查是否已有可复用的产出**（`logs/`、既有 BAM/GB/候选目录），避免重算；
> ② 记录**工具版本、完整命令、参数、输入 hash、输出路径与实际退出状态**。
> 本文件只描述现有脚本，不引入平行路径。参数以脚本 `--help` 为准；环境见 `tool_check.md`。

| 任务 | 工具 | 输入要求 | 可回答的问题 | **不能**证明的结论 | 产出 / 退出码 |
|---|---|---|---|---|---|
| 序列体检 | `scripts/seq_stats.py <assembly.fasta> [--window 500]` | 任意 FASTA | 长度、contig 数、GC/AT（仅以 ACGT 为分母）、模糊碱基数量与**逐 contig 位置（0-based）**、滑窗 GC | 不证明 reads 来源，不证明注释正确 | 只读；退出码 0 |
| 注释质检 | `scripts/annot_check.py <ann.gb> [--ref <ref.gb>] [--table 5] [--taxon <name>] [--require-circular] [--tolerate-overlap "G1,G2"] [--tolerate-start "GENE:CODON"] [--exception-registry <json>] [--overlap-severity error\|warn] [--allow-atypical "<理由>"]` | 带注释的 GenBank | 基因数量/身份/重复、tRNA 类型是否可判定、CDS 翻译与起始终止、location 的 5'/3' partial、`/codon_start` 矛盾、`/transl_table` 冲突、`/transl_except` 是否真解释了该内部终止、tRNA/rRNA 长度、**逐对**重叠（分段坐标）、链分布与方向、（可选）环状邻接与长度对比 | 注释内部一致 ≠ 序列被 reads 支持（脚本结尾固定打印该免责声明）；`--require-circular` 只校验 topology 声明 | 只读；`0` 无发现 / `1` 有错误（含参数错误）/ `2` 仅待核查 |
| 基因定位 | `scripts/blast_genes.py <ref.gb> <target.fna> [--out gene_order.txt] [--evalue 1e-5] [--min-identity 80] [--min-coverage 0.8]` | 近缘参考 GB + 目标 FASTA；需 `blastn` | 每个参考基因在目标上的**唯一**定位、链与顺序（identity ≥ 80% **且**覆盖 ≥ 80%） | 短/局部命中被刻意排除，故**不能**单独判定真实重排；顺序差异须人工 REVIEW | 写 `gene_order.txt`（1-based）；任一基因无唯一命中或命中歧义即退出 `1` 且**不写部分结果** |
| reads 裁判 | `scripts/depth_analysis.py <sample.bam> <genome.fasta> [--window 500] [--allow-base-replacement] [<ambiguous.fasta>]` | **已排序并建索引**的 BAM + 参考 FASTA；需 `samtools` | 覆盖度剖面与低覆盖区、mate 分布、soft-clip 比例；对每个模糊碱基给出 `base_support()`：(深度, MAPQ, 碱基质量, 链向, 支持率, 各类排除数) | 低覆盖 ≠ 必然错接；排除 NUMT 需核参考；四项证据任一项不过就**不得**替换碱基 | 只读；输入/BAM 错误 → `1`；存在低覆盖区 → `2` |
| 环化候选 | `scripts/circularize.py <s1.fa> <s2.fa> <ref.gb> <outdir> --bam <candidate.bam> --reads-validated [--min-junction-support 3] [--min-mapq 20] [--junction-region chr:s-e] [--junction-flank 10] [--accept-candidate]` | 两个 scaffold、参考 GB、候选回贴 BAM；需 BioPython + `blastn` + `minimap2` + **`samtools`** | 4 种拼接组合中**唯一最高分**方向、自动定位的新增接缝、该接缝的 `junction_evidence()`（独立分子数、链向、MAPQ min/median、剔除的重复标记数） | 默认输出 `PUTATIVE_CIRCULAR/REVIEW`：单一区间证据**不能**宣称最终环化；候选顺序/方向与参考不一致时记 **REVIEW** 且**不允许** `--accept-candidate`；`--min-junction-support` 是工程下限 | 写 `outdir/genome_candidate.fasta` + 验证目录；`2` = REVIEW/需人工，`1` = 接缝证据不足 |
| COX1 查询 | `scripts/cox1_id.py <genome.fasta> --allow-public-upload --coords start,end [--max-results 5] [--output-json results.json]` | **必须**显式给 `--coords`（不自动识别 COX1）；查询 400–5000 bp；**会向公共 NCBI 上传该片段** | 候选取样片的 identity、**query coverage**（多 HSP query 区间并集）、多个候选与 `provisional_candidate`/`ambiguous`/`insufficient` 判读；请求 `FORMAT_TYPE=XML2` 并**同时**支持 XML2（`-outfmt 16`）与旧版 XML（`-outfmt 5`）；`SearchInfo` 按 NCBI 实际返回的 QBlast 文本（`RID =`/`Status=`）解析；原始 HSP 打印在 CLI 上并用 `--output-json` 持久化 | identity/coverage **不等于**物种鉴定；**多 HSP 指标不可靠时不汇总单一 identity**：query 区间重叠 → `overlapping_hsps`；**target 区间重叠**（两条 query 打到目标同一段）→ `subject_overlap`；同一目标链但坐标**逆链方向**推进（负链要求 `hit_from` 递减）→ `non_collinear_hsps` + `CONFLICTING_ALIGNMENT`；正负链混合 → 矛盾证据；跳变触及两端 → `CROSS_ORIGIN_CANDIDATE`（环状参考需额外结构证据）。以上均**不参与自动择优**，但**不等于该 hit 不是有效候选**；覆盖率 <80% 或 top 与次优相差 <1% 不给结论 | 上传到 NCBI；退出码 `0` 得判读 / `1` insufficient 或 no_match / `2` 拒绝执行 / `3` 网络或结果格式故障 |
| 经验与案例 | `tools/experience.py search \| case-init/case-event/case-validate/case-report \| propose-lesson \| review-lesson \| export-contribution \| sync-public` | 工作目录可写；`MITO_KNOWLEDGE_DIR` 可写 | 检索历史案例/lesson、结构化案例记账（`case-init` **必须**给 `--hypothesis`）、`case-validate` 按 `schemas/case.schema.json` 校验、候选 lesson 与审核、可预览贡献、公共同步 | `case-validate` 只证明**记录格式**合格（schema/字段/枚举），不证明科学结论；公共同步只读入隔离缓存，**不执行**远程内容 | 写 `MITO_KNOWLEDGE_DIR`；校验失败 → `1` |
| 环境自检 | `bash scripts/check_env.sh` | `config/env.sh` 可读 | 核心工具、组装/注释工具、MITOS2 模块与数据库是否就位 | 不验证数据正确性 | 只读；`2` = 核心依赖缺失 |
| 后台化 | `bash scripts/run_bg.sh <名> -- <命令...>` / `--trusted-shell "<命令>"`；`bash scripts/check_bg.sh <名> [--tail N]` | 可信命令/参数 | 长任务（>5 分钟）非阻塞运行与状态轮询 | 默认模式不解析 shell 语法；只有确认可信时才用 `--trusted-shell`；不得把 running/timeout 描述为完成 | 写 `logs/<名>.log/.pid/.status` |
| 注释（统一入口） | `bash scripts/run_mitos2.sh -i <genome.fasta> -o mitos2_out -c 5` | `config/env.sh` 的 `MITOS2_PY`/`MITOS2_REFDIR`/`MITOS2_REFSEQVER` | MITOS2 完整注释（CDS/tRNA/rRNA/起始终止） | 单工具结果不是真值；个别基因失败时需 reads + 蛋白补充 | 写注释输出目录 |
| 环形图 | `bash scripts/run_circular_map.sh <final.gb> <out.png> --title "..."` | GB **必须声明 `topology=circular`**，否则脚本直接退出；需 `PLOT_PY`（BioPython + matplotlib） | 论文级正负链双环基因图 | 不影响序列/注释正确性 | 写 300dpi PNG + SVG |
| NCBI 提交预检 | `table2asn`（外部工具） | 需自行提供 template；**不在本 skill 的环境自检范围内** | 提交前的本地格式/内部终止等预检 | 通过预检 ≠ 成果已被 NCBI 接受 | 人工执行，不进 `check_env.sh` |

## 退出码接口约定（合并后不要误读）

| 现象 | **不是** | 正确读法 |
|---|---|---|
| `annot_check.py` 退出码 `2` | ❌ 注释完全通过 | ✅ 无 ERROR，但**存在待核查项**；每条须逐条入账；只有 `0` 才是"无发现" |
| `cox1_id.py` 退出码 `3` | ❌ 没有找到相似序列 | ✅ 网络或结果格式**故障**；必须与 `1`（insufficient / no_match，分析结论）严格区分 |
| 后台任务 "子进程正常退出" | ❌ 结果已通过科学验收 | ✅ 只表示命令执行成功；科学验收由相应的证据标准另行判定 |

任何脚本、CI 或 AI 指令都不得把上表左列当成右列。`baseline_report.md` 是 v1 历史审计快照，
描述的是当时行为，不作为现行接口依据。

## 使用原则

1. **工具优先**：优先专业工具（MITOS2 / MitoFinder / GetOrganelle / metaSPAdes / bwa / minimap2 / samtools），
   只在现成工具无法解决时才写补充代码，且补充代码必须能交叉验证。
2. **按诊断价值调用**：新增检查必须能改变假设排序或构成验收条件；否则不跑。
3. **阈值分层**：`>8bp` 重叠、tRNA `60–75bp`、rrnL/rrnS 区间、`9+/4−` 都是本工具的**工程预警**，
   不是 NCBI 硬性标准；放宽要显式使用对应参数并记录理由。
4. **不覆盖原始数据**：所有候选写入新目录，原件保持只读。
5. **不把私有路径写回 skill**：机器路径只存在于 `config/env.sh` / 环境变量。

## 相关文件

- 判据与阈值：`annotation_quality.md`；顺序与坐标：`standard_gene_order.md`
- 证据门槛与 `confidence`：`evidence-standard.md`
- 环境变量与依赖自检：`tool_check.md`
- 已知坑：`MITO_KNOWLEDGE_DIR` 的 `pitfalls.md`

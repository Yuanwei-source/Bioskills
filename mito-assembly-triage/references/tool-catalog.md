# 现有工具箱（按需调用）

> 运行前先做两件事：**① 查是否已有可复用的产出**（`logs/`、既有 BAM/GB/候选目录），避免重算；
> ② 记录**工具版本、完整命令、参数、输入 hash、输出路径**。
> 本文件只描述现有脚本；不引入平行路径。环境变量与依赖自检见 `tool_check.md`。

| 任务 | 工具 | 输入要求 | 可回答的问题 | **不能**证明的结论 | 产出 / 是否写文件 |
|---|---|---|---|---|---|
| 序列体检 | `scripts/seq_stats.py <assembly.fasta> [--window 500]` | 任意 FASTA | 长度、contig 数、GC/AT（仅以 ACGT 为分母）、模糊碱基数量与**逐 contig** 位置、滑窗 GC 均匀性 | 不证明 reads 来源，不证明坐标正确性 | 只读，stdout |
| 注释质检 | `scripts/annot_check.py <ann.gb> [--ref <ref.gb>] [--table 5] [--require-circular] [--tolerate-overlap "G1,G2"] [--overlap-severity error\|warn] [--allow-atypical]` | 带注释的 GenBank | 基因数量/身份/重复、CDS 翻译与起始终止、tRNA/rRNA 长度、重叠分级、链分布、（可选）顺序与长度对比 | 注释内部一致 ≠ 序列被 reads 支持；典型动物基因集不是绝对真值 | 只读，stdout + 退出码 0/1/2 |
| 基因定位 | `scripts/blast_genes.py <ref.gb> <target.fna> [--out gene_order.txt] [--evalue 1e-5] [--min-identity 80] [--min-coverage 0.8]` | 近缘参考 GB + 目标 FASTA；需生物信息 `blastn` | 每个参考基因在目标上的**唯一**全长定位、链与顺序；命中歧义直接报错 | 短/局部同源命中被刻意排除，故**不能**单独判定真实重排；顺序差异须人工 REVIEW | 写 `gene_order.txt`（1-based closed） |
| reads 裁判 | `scripts/depth_analysis.py <sample.bam> <genome.fasta> [--window 500] [--allow-base-replacement] [<ambiguous.fasta>]` | **已排序并建索引**的 BAM（bwa mem 产出）+ 参考 FASTA；需 `samtools` | 覆盖度剖面与低覆盖区、低覆盖区 mate 分布、soft-clip 比例、模糊碱基的 reads 支持 | 低覆盖≠必然错接；排除 NUMT 需核参考；模糊碱基替换需支持率 >95% **且**显式 `--allow-base-replacement` | 只读，stdout + 退出码（有低覆盖区 → 2） |
| 环化候选 | `scripts/circularize.py <s1.fasta> <s2.fasta> <ref.gb> <outdir> --bam <candidate.bam> --reads-validated [--min-junction-support 3] [--min-mapq 20] [--junction-region chr:s-e] [--junction-flank 10] [--accept-candidate]` | 两个 scaffold、参考 GB、候选序列回贴 BAM；需 BioPython + `blastn` + `minimap2` | 4 种拼接组合中**唯一最高分**的候选方向、自动定位的新增接缝、该接缝的独立跨接 read 数 | 默认输出 `PUTATIVE_CIRCULAR/REVIEW`：单一区间证据**不能**宣称最终环化；需人工核对全部新增接缝后才能 `--accept-candidate` | 写 `outdir/genome_candidate.fasta` 与验证目录 |
| COX1 查询 | `scripts/cox1_id.py <genome.fasta> --allow-public-upload --coords start,end [--max-results 5]` | **必须**显式给 `--coords`（脚本不自动找 COX1）；查询 400–5000 bp；**会向公共 NCBI 上传该片段** | 候选取样片的 top 命中 identity 与近似物种 | identity ≠ 确定物种；无唯一且 ≥85% 的命中会直接失败；倾斜命中（top 与次优 <1%）判为不可判 | 上传到 NCBI，stdout |
| 经验与案例 | `tools/experience.py search \| case-init/case-event/case-validate/case-report \| propose-lesson \| review-lesson \| export-contribution \| sync-public` | 工作目录可写；`MITO_KNOWLEDGE_DIR` 可写 | 检索历史案例/lesson、结构化案例记账、候选 lesson 与审核、可预览贡献、公共同步 | 不代替诊断证据；公共同步只读入隔离缓存，**不执行**远程内容 | 写 `MITO_KNOWLEDGE_DIR`（本地）/ 显式 `--output` |
| 环境自检 | `bash scripts/check_env.sh` | `config/env.sh` 可读 | 核心工具、组装/注释工具、MITOS2 模块与数据库是否就位 | 不验证数据正确性 | 只读，退出码 2 = 核心依赖缺失 |
| 后台化 | `bash scripts/run_bg.sh <名> -- <命令...>` / `--trusted-shell "<命令>"`；`bash scripts/check_bg.sh <名> [--tail N]` | 可信命令/参数 | 长任务（>5 分钟）非阻塞运行与轮询 | 默认模式不解析 shell 语法；只有确认可信时才用 `--trusted-shell` | 写 `logs/<名>.log/.pid/.status` |
| 注释（统一入口） | `bash scripts/run_mitos2.sh -i <genome.fasta> -o mitos2_out -c 5` | `config/env.sh` 的 `MITOS2_PY`/`MITOS2_REFDIR`/`MITOS2_REFSEQVER` | MITOS2 完整注释（CDS/tRNA/rRNA/起始终止） | 单工具结果不是真值；个别基因失败时需 reads+蛋白补充 | 写注释输出目录 |
| 环形图 | `bash scripts/run_circular_map.sh <final.gb> <out.png> --title "..."` | GB + `PLOT_PY`（BioPython + matplotlib） | 论文级正负链双环基因图 | 不影响序列/注释正确性 | 写 300dpi PNG + SVG |

## 使用原则

1. **工具优先**：优先专业工具（MITOS2 / MitoFinder / GetOrganelle / metaSPAdes / bwa / minimap2 / samtools），
   只在现成工具无法解决时才写补充代码，且补充代码必须能交叉验证。
2. **按诊断价值调用**：新增检查必须能改变假设排序或构成验收条件；否则不跑。
3. **不覆盖原始数据**：所有候选写入新目录，原件保持只读。
4. **不把私有路径写回 skill**：机器路径只存在于 `config/env.sh` / 环境变量。

## 相关文件

- 硬阈值与判定等级：`annotation_quality.md`
- 证据门槛：`evidence-standard.md`
- 环境变量与依赖自检：`tool_check.md`
- 已知坑：SKILL.md「常见工具坑速查」+ `MITO_KNOWLEDGE_DIR` 的 `pitfalls.md`

# 从这里开始

用户可提供原始 reads、已有组装 FASTA，或组装与注释结果，无须先知道哪里错。先盘点文件并核对样本身份与文件对应关系，
再选择路线；读取相关资料以判断每项异常是否能修。详细循环见
[diagnostic-playbook.md](diagnostic-playbook.md)，不要默认运行全部工具。

## 1. 按输入选路线

| 已有数据 | 入口动作 | 证据边界 |
|---|---|---|
| FASTQ/read only | 先识别平台、读长、文库和样本信息；选兼容的线粒体组装工具生成候选，再进入 FASTA-only 路线 | 组装成功不等于图唯一、序列完整或已闭环；不能支持的技术类型须停止并说明 |
| FASTA only | 见下方“只有组装 FASTA”路线 | 可生成注释候选并联合复核；无 reads 时不能验证样本碱基/连接支持 |
| FASTA（仅做统计） | `python3 scripts/seq_stats.py <assembly.fasta>` | 统计工具；多记录时不应将汇总 GC 当作一条分子属性 |
| FASTA + GenBank | 按类群确认密码表，运行 `annot_check.py` | 注释一致性不等于样本序列正确 |
| MITOS2 输出 | 用 `mitos2_to_genbank.py --table <已确认密码表>` 转换后质检 | 此转换器针对 MITOS2 文件组合，不是通用 GFF 转换器；拓扑仍只是声明 |
| FASTA + FASTQ/BAM | 按读长与文库比对、排序、建索引，记录重复标记/处理；定点检查 reads | 先确认 BAM 与候选参考匹配；单参考不能排除 NUMT |
| reads + 核/其他竞争参考 | 比较候选、核位点或竞争路径 | 存在竞争参考不自动获得判别力，须检查实际结果 |
| 多 contig / GFA | 先看图与连接假设；两 scaffold 场景可用 `circularize.py` | 脚本只自动验证一个内部接缝，闭合连接另需验证 |
| 需要公共参考 | 先看本地资源；按 [reference-policy.md](reference-policy.md) §5 选等级和用途，获取后登记 | 参考是对照；序列下载与样本上传分开授权 |

先核对已有产物的输入、版本和用途；兼容的产物可复用，条件变化时才重跑。
FASTA 和 GenBank 同时存在时，先逐碱基核对 GenBank 序列与 FASTA，再审查 source、taxon、
采集字段与用户声明的样本是否一致。能确认的注释字段错误可用 `apply_candidate.py qualifier`
按旧值精确修改；输出候选后复跑注释检查，并记录仍需人工判断的基因身份与边界。

## 2. 只有动物／昆虫 Illumina 双端 reads：多软件候选组装，再诊断和注释

首版 reads-only 自动组装范围限定为动物／昆虫 Illumina paired-end 短读长。先收集样本名/学名、读长、插入片段信息、数据文件及 lane/批次关系；核验 FASTQ 配对、完整性和质量。多个 lane 属于同一文库时按 read 方向合并，不能把 R1 与 R2 互相拼接。原始 reads 不改写。

不同生物信息软件可分别安装在多个 conda 环境。设置 `CONDA_ROOT` 后，依赖检查和调度脚本会在 `CONDA_ROOT/envs/*/bin` 定位可执行文件并以绝对路径调用，不需要把软件装进同一个环境或激活一个“大环境”。MITOS2 的 Python 与数据库路径依旧分别使用 `MITOS2_PY` 和 `MITOS2_REFDIR` 配置。

若 MitoFlex/NOVOPlasty 是源码脚本而不在 `envs/*/bin`，分别设置 `MITOFLEX_ROOT` +（可选）`MITOFLEX_PYTHON` 和 `NOVOPLASTY`。调度器会检查 MitoFlex Python 的 `numpy`、`pandas`、`ete3`、Biopython、`psutil`；若只缺 `ete3`，会验证能否从其它 conda 环境补充其 site-packages。NOVOPlasty 脚本由 Perl 启动，不要求可执行位。机器路径放在本地环境变量或未跟踪的机器配置中，不提交到 skill。

MitoFlex 的 Rust FASTQ 过滤器只读取一个 gzip member；若原始 `*.fastq.gz` 是用 `cat` 合并的多个 gzip members，调度器会在该工具的中间目录重写为单-member gzip，并在 manifest 记录源/规范化文件的哈希。原始 reads 保持不变。该转换避免 MitoFlex 静默只读第一个 lane。

先运行 `CONDA_ROOT=<conda根目录> python3 tools/env_check.py --stage assembly_compare_illumina --json` 检查基础多工具路线。候选组装器默认尝试 GetOrganelle、MitoFlex、NOVOPlasty 和 MitoZ；可运行的组装器默认最多同时运行 3 个，每个工具分别使用 `--threads` 指定的线程数。MitoZ 需显式提供遗传密码表。若工具或 NOVOPlasty 所需的种子/文库参数不可用，清单会注明未执行；有近缘参考时可显式加入 MitoFinder。GetOrganelle 动物线粒体模式为 `animal_mt`，还需初始化对应数据库。

长时间流程应通过任务管理器后台运行，任务目录放在样本目录的 `intermediate/tasks/`。组装器会在调度器创建的任务目录旁写入 `intermediate/assemblies/`，不会覆盖旧结果：

```bash
python3 scripts/task_manager.py start \
  --task-root <样本名>_mitogenome/intermediate/tasks \
  --name illumina_assembly --cwd "$PWD" -- \
  env CONDA_ROOT=<conda根目录> <python解释器> scripts/run_illumina_candidates.py \
  --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> \
  --sample <样本名> --taxon '<确认的学名>' \
  --outdir <样本名>_mitogenome \
  --tools getorganelle,mitoflex,novoplasty,mitoz --table <遗传密码表编号>
```

NOVOPlasty 只有在额外提供 `--novo-seed`、`--insert-size` 和经类群判断的 `--genome-range` 后才运行；它不会暗中拿 GetOrganelle 结果当独立种子。MitoZ 路线需要显式 `--table`，MitoFinder 需要 `--reference-genbank`。MitoFlex 默认加 `--disable-taxa` 跳过内部近缘分类过滤，避免未知/缺失类群被误删；并用上游支持的 `--disable-annotation`，将后续注释交给独立比较流程；若用户有与其内部数据库匹配的分类配置，可显式传 `--mitoflex-use-taxonomy-filter`。调度器返回任务目录后，用它查询状态和日志；若显示 `orphaned_running` 或 `supervisor_lost`，先确认实际进程状态，不要直接重复提交；状态文件记录最终退出码：

```bash
python3 scripts/task_manager.py list --task-root <样本名>_mitogenome/intermediate/tasks
python3 scripts/task_manager.py status --task-dir <返回的任务目录>
python3 scripts/task_manager.py log --task-dir <返回的任务目录> --lines 100
```

生成 `assembly_candidates.json` 后、判读候选前，完成原始 reads 的基础质量报告：

```bash
python3 scripts/task_manager.py start \\
  --task-root <样本名>_mitogenome/intermediate/tasks \\
  --name raw_fastq_qc --cwd "$PWD" -- \\
  python3 scripts/fastq_qc.py \\
  --assembly-manifest <样本名>_mitogenome/intermediate/assemblies/assembly_candidates.json \\
  --outdir <样本名>_mitogenome/intermediate/quality/raw_fastp --threads 4
```

质量 JSON 会校验 fastp 处理的 read count 是否等于完整 FASTQ 配对检查值；若不相等（例如工具只读了一个 gzip member），本次质量结果标为失败，不能拿部分数据代表全量 reads。

组装任务终止后，先做同一竞争参考下的候选 reads 回贴；该任务也应后台运行：

```bash
python3 scripts/task_manager.py start \\
  --task-root <样本名>_mitogenome/intermediate/tasks \\
  --name candidate_competitive_mapping --cwd "$PWD" -- \\
  python3 scripts/compare_assembly_candidates.py \\
  --assembly-manifest <样本名>_mitogenome/intermediate/assemblies/assembly_candidates.json \\
  --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> \\
  --outdir <样本名>_mitogenome/intermediate/validation/candidate_mapping \\
  --threads 4 --min-mapq 20 --min-baseq 20
```

结构化结果在 `candidate_evidence.json`，同时保留竞争参考、target 映射表、排序 BAM、索引、覆盖表、输入哈希、版本、参数与命令日志。它只比较线性候选上的主比对证据，不计算圆形首尾接缝、不排除 NUMT、不自动选出最终组装；接缝和拓扑必须另行验证。

某个组装器失败时，先看 `intermediate/assemblies/<工具名>/tool.log`。修正环境或参数后，对同一 FASTQ 使用 `--resume --tools <失败工具>` 重试；resume 会验证输入路径与 SHA256、保留旧日志和失败记录，并拒绝覆盖已经成功的工具结果。

每个可信组装候选进入独立注释比较：

```bash
CONDA_ROOT=<conda根目录> python3 scripts/run_candidate_annotations.py \
  --candidate-fasta <候选1.fasta> --candidate-fasta <候选2.fasta> \
  --outdir <样本名>_mitogenome/intermediate/annotation_candidates \
  --taxon '<用户确认的学名>' --table <确认或暂定的NCBI密码表> \
  --table-status <confirmed|provisional> --clade Arthropoda
```

MITOS2 和 MitoZ 的每个输出都是候选；最终注释需统一坐标、方向、密码表后逐基因比较，不能按工具多数票决定。随后按 [FASTQ/BAM 路线](evidence-standard.md) 将原始 reads 回贴到竞争候选，比较碱基、覆盖和每个连接；检查组装图、重复和闭环声明。报告应保留工具失败和未执行原因。GetOrganelle 上游 reads 路线面向 Illumina 单端/双端 FASTQ；长读长和混合数据不属于首版支持范围。

注释比较可直接使用 `--assembly-manifest <assembly_candidates.json>`，脚本会校验候选 FASTA 的 SHA256。多记录候选默认安全跳过；加 `--split-multi-records` 可让 MITOS2 和 MitoZ 分别注释每条 contig，manifest 保留它来自哪一个组装及原始 record ID。逐条 contig 注释只用于发现局部基因线索，不等于把碎片恢复为完整组装，也不应直接作为最终交付。

每个可解释的组装候选都进入“只有组装 FASTA”路线：初检、按物种与密码表生成注释候选、联合检查、证据支持的修复、复验及结果交付。reads 组装路线的图、日志、候选对比和比对文件收在最终 `<样本名>_mitogenome/intermediate/` 中；根目录仍只提供 FASTA、GFF3、GenBank 和报告。

## 3. 按异常找资料

| 现象 | 最先检查 | 读取 |
|---|---|---|
| 少基因 / 少 tRNA | 名称、同源、邻域与跨原点位置；区分未检出与丢失 | [annotation_quality.md](annotation_quality.md) §2–§3 |
| 内部 stop / 移码 | 密码表、坐标、读框与例外证据 | [annotation_quality.md](annotation_quality.md) §2 |
| 环不起来 / 两端接不上 | 重复、唯一锚定、跨接证据与数据分辨力 | [standard_gene_order.md](standard_gene_order.md) §3 |
| 深度异常 / NUMT | 比对质量、重复、mate 与竞争参考 | [evidence-standard.md](evidence-standard.md) §3 |
| 顺序不同 | 旋转、整链反向互补与类群适当参考 | [standard_gene_order.md](standard_gene_order.md) |
| N / 模糊碱基 | 定位、局部 reads、混合信号与比对歧义 | [evidence-standard.md](evidence-standard.md) §1 |
| COX1 物种线索 | 本地确定坐标；远程上传须授权 | [tool-catalog.md](tool-catalog.md) 的 COX1 行 |

陌生现象先写“现象 + 位置 + 前提 + 可用数据”，再列竞争解释。

## 4. 只有组装 FASTA：从组装初检到注释联合复核

这条路线面向“已有人工组装，但尚无 GenBank 注释”的情况。它不从 reads 从头组装；
用户有 FASTQ/BAM 时把它作为后续独立证据。先建立任务目录并登记 FASTA：

```bash
python3 scripts/case.py init <workdir> --case-id <sample-id> \
  --issue '线粒体组装与注释联合审核' --taxon '<已确认学名>' \
  --input assembly <assembly.fasta>
python3 scripts/assembly_intake.py <assembly.fasta> --json <workdir>/intake.json
```

初检会按 FASTA 记录分别报告长度、GC（只以 A/C/G/T 为分母）、IUPAC 模糊碱基区间和输入哈希；
它不串接 contig、不判定序列正确性或环状闭合。多记录时先厘清记录关系和连接；
不能只因线粒体通常为环状而拼接。模糊碱基只能由样本 reads 等证据支持替换；缺 reads 时保留并标注，
仍可注释其余可判读区域，但跨模糊区的基因/边界应视为受限。

确定物种适用的遗传密码表及其依据后，用一条 FASTA 运行候选生成路线：

```bash
bash scripts/run_fasta_annotation.sh --fasta <assembly.fasta> \
  --outdir <sample>_mitogenome/intermediate/annotation_route --table <NCBI-table> \
  --table-status <confirmed|provisional> --organism '<sample scientific name>' \
  --topology linear
```

遗传密码表来源不足时标记 `provisional`；路线不会替用户确认它。`--organism` 按用户声明写入 source；
只有 TaxID 有可靠依据时才传 `--taxid`。`--topology` 是输出记录中的声明，不是闭环证据。目录必须是新目录；路线保留 MITOS2 原始文件、
初检、GenBank 转换日志、`annot_check.py` 日志与机器可读的 `route.json`。多记录输入会安全停止，
因为当前桥接器要求单记录；先单独诊断 contig/连接，不能简单拼成一条序列。若 MITOS2/数据库/Infernal
未配置，查看日志及 `bash scripts/check_env.sh --stage annot_independent`，补齐项目 `config/env.sh`
所指环境后再运行；也可在单次命令中传入对应环境变量而不改项目配置。失败不代表 FASTA 或注释通过。
已有 MITOS2 `result.gff/fas/faa` 时可传 `--mitos-dir` 重用，只重新执行转换与质检。

随后不能把 `annotation.gb` 直接交付为专家注释。必须逐项联合审查：

1. 核对 GB 记录序列与输入 FASTA 完全一致、密码表/物种/来源字段正确；核对 MITOS2 的基因身份、
   CDS 读框与边界、partial/碎片、tRNA 反密码子和类群已知例外。`annot_check` 只筛内部不一致，
   其 0 退出码不是专家验收或真实性证明；昆虫默认阈值不能直接推广到其他类群。
2. 把初检发现与注释问题一起形成问题清单，逐条提出竞争解释、所需证据和处理决定。
   FASTA 暴露的碱基异常、端部重复或接缝疑点需要 reads、装配图或独立证据；不能让参考序列替用户样本定碱基。
3. 证据充分时才用候选工具分别修改碱基、注释或结构；候选不得覆盖输入。修改后重新运行对应检查，
   并核对坐标/相邻基因/序列哈希变化。证据不足时保留 `UNRESOLVED`，而不是为了得到完整文件强改。
4. 按任务记录流程登记 MITOS2、转换、检查及修复验证事件；将日志、清单和机器可读记录放入 `intermediate/`。
   最终目录根部交付 `<sample>_assembly.fasta`、`<sample>_annotation.gff3`、`<sample>_annotation.gb` 和
   `<sample>_report.md`。报告按 [report-template.md](report-template.md) 填写，说明修改差异、证据、用途建议和未解决事项。reads 未提供时明确标注
   “序列/注释一致性已审查，样本碱基与连接支持未评估”。

这个自动入口只实现“FASTA 初检 → MITOS2 候选 → GenBank 桥接 → 初步质检”。所谓专家级结果必须来自之后
的类群化人工/证据复核与修复后验证；没有一个工具能仅凭单条 FASTA 证明组装的样本真实性或每个注释边界正确。

## 5. 条件式首次检查

1. 明确本次目标与类群，查现有结果；只检查将要调用的工具依赖。
   第一次/换机/升级后：`bash scripts/check_env.sh --setup`（全量盘点；essential 齐全时记录 lock，
   缺失项给出用途与安装方式）；以后日常用 `--daily`（读 lock 做轻量检查，环境变了会报错，
   不会静默继续）；单步可用 `--stage annot_check`（缺依赖则以退出码 3 停下并注明该步骤未执行）。
   缺少与当前问题无关的工具不阻断当前检查。
2. 需要序列属性时运行 `seq_stats.py`；需要注释检查且有 GB 时执行下例。
   已明确是注释问题可先检查注释，注明结构尚未验证。
3. 怀疑碱基/连接问题且有 reads 时，再做定点 reads 检查；缺输入时只限制相应命题。
4. 记录事实与竞争解释；进入 [diagnostic-decision-tree.md](diagnostic-decision-tree.md)
   决定继续、改变路线还是停止。

```bash
# 占位符须换成已确认的值；不能默认将表 5 用于所有动物
python3 scripts/annot_check.py <ann.gb> --table <已确认的密码表编号> --taxon '<类群>'
# 返回 2 表示有待核查项，返回 1 表示错误；具体含义见工具目录
```

任务 case/日志可写入工作目录。首次记账的完整示例见诊断手册；
授权与跨任务经验边界见 [learning-policy.md](learning-policy.md)。
输出前阅读 [conclusion-report.md](conclusion-report.md)，用 case.py report 生成证据报告，
再按实际证据核对科学解释。

## 任务记录与交付

完整诊断从 `scripts/case.py init` 开始登记输入；实际检查后登记事件和逐条结论，
再生成报告。复制可用的命令及 JSON 模板见 [task-records.md](task-records.md)。

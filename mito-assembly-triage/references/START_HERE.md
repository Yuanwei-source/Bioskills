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

## 2. 只有原始 reads：候选组装，再诊断和注释

先收集样本名/学名、测序平台、单双端与文库信息、读长、数据文件及其 lane/批次关系；核验 FASTQ 配对、完整性和质量。多个 lane 属于同一文库时按 read 方向合并，不能把 R1 与 R2 互相拼接。原始 reads 不改写。

先运行 `python3 tools/env_check.py --stage assembly_from_illumina_reads` 确认依赖。GetOrganelle 上游当前说明其 reads 路线面向 Illumina 单端/双端 FASTQ，动物线粒体用 `animal_mt`；需安装并初始化对应数据库。示例命令如下，参数仍须按当前版本 `--help`、read 长度/质量和文库特征评估：

```bash
get_organelle_from_reads.py -1 <R1.fastq.gz> -2 <R2.fastq.gz> \
  -R 10 -k 21,45,65,85,105 -F animal_mt \
  -o <sample>_mitogenome/intermediate/GetOrganelle
```

保留并检查版本、数据库、运行参数、组装图和所有候选连接。GetOrganelle 的 reads 接口不直接支持 ONT/PacBio 原始 reads；长读长或混合数据须选用适配平台的组装工具，项目目前没有统一的跨平台 reads-only 自动入口。缺少匹配工具、数据质量不足或图结构有多解时，停止在证据支持的候选层级并列明需要补充什么。

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

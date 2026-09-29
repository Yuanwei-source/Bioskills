# 工具目录（按需调用）

先核对已有产物是否可复用，再选择当前检查所需工具。记录版本、完整命令、参数、输入
hash、产物路径和实际退出状态。参数以各工具 `--help` 为准，依赖见
[tool_check.md](tool_check.md)。以下路径相对 skill 根目录；实际输出放在用户任务目录。

## 诊断与候选工具

| 任务 / 工具 | 必要输入与用途 | 关键限制 | 输出与退出码 |
|---|---|---|---|
| 明确编辑清单候选：`apply_candidate.py` | FASTA 碱基替换、唯一 CDS 边界调整、GenBank 注释字段精确编辑；需 edits JSON 和证据引用 | 只写新候选和变更 manifest，不推断编辑、不做科学验收；边界需显式密码表 | 0 生成候选 / 1 输入、旧值或输出保护失败 / 3 缺依赖 |
| 原始 FASTQ 质量：`fastq_qc.py` | 已完成配对/文件哈希登记的 assembly manifest；fastp | 质量报告前后 counts 必须与全量 FASTQ 检查一致；stdout reads 被丢弃，原始 FASTQ 只读 | `fastq_quality.json` + fastp JSON/HTML/log；COUNT 漂移或工具故障非 0 |
| Illumina reads-only 组装：GetOrganelle `get_organelle_from_reads.py` | 配对或单端 Illumina FASTQ；`-F animal_mt` 数据库 | 只作为候选组装；检查 graph/path、重复、分支和所有连接，工具输出 circular 不能单独证明拓扑；不适用于直接输入 ONT/PacBio 原始 reads | 0/非 0 依上游版本；须保留完整日志、graph 和所有候选序列 |
| Illumina 组装：MitoFlex `MitoFlex.py all` | 配对 FASTQ、MitoFlex clade 与遗传密码表 | skill 默认关闭其近缘分类过滤以避免数据库缺类群时误删候选；仍保留 HMM/深度等上游筛选，最终必须独立比对和裁决 | 原始 MitoFlex 工作目录与 FASTA；命令记录 `--disable-taxa`，可用 `--mitoflex-use-taxonomy-filter` 显式启用 taxonomy 过滤 |
| 多组装候选竞争回贴：`compare_assembly_candidates.py` | 组装 manifest、配对 FASTQ；BWA-MEM + samtools；对所有候选 FASTA 记录建立一个竞争参考 | 仅保存主比对；覆盖与 proper-pair 是相对支持证据，不排除 NUMT、不判闭环、不自动排序或选择；当前不评价圆形原点接缝 | `intermediate/validation/` 下 BAM、coverage TSV、target map 与结构化证据 JSON；成功仍需人工/下游逐候选解释 |
| FASTA-only 初检：`assembly_intake.py` | 单/多记录组装 FASTA；逐条长度、GC、模糊碱基区间、SHA-256 | 只报告 FASTA 可见事实；不判断分子身份、碱基真实性、方向或闭环 | 0 报告 / 1 FASTA 格式或读写错误；可用 `--json` 保存结构化报告 |
| 任务记录/修复追溯：`case.py` | 输入、事件、结论及已有工具生成的候选；[用法](task-records.md) | 不生成或采纳修复；verified 仅表示登记标准结果，不证明科学结论 | 0 完成 / 1 非法输入或文件漂移 / 3 缺依赖 |
| 序列体检：`seq_stats.py` | FASTA；长度、GC/AT、模糊碱基和滑窗 | 模糊位置是 0-based；不验证样本真实性 | 统计输出；正常 0 |
| 注释质检：`annot_check.py` | GenBank；显式传类群确认的 `--table`，必要时 `--taxon`、`--ref` | 无 reads 验证；`--taxon` 不自动切换昆虫阈值，拓扑检查只校声明 | 0 无发现 / 1 错误（含参数错误）/ 2 待核查 |
| 注释格式桥：`mitos2_to_genbank.py` | MITOS2 GFF/FAS/FAA、单记录组装 FASTA、必需 `--table`；生成可质检 GB | 非通用 GFF 转换；不补基因/碱基，不导出未知 partial，不是提交质量记录 | 写 GB；输入/坐标异常非 0，无半成品 |
| FASTA-only 注释候选路线：`run_fasta_annotation.sh` | 单记录组装 FASTA、显式密码表及证据状态、可用 MITOS2 环境；依次初检、注释、桥接、质检 | 不修组装、不证明样本真实性或闭环；MITOS2 注释仍需类群化独立复核；可 `--mitos-dir` 重用已有 MITOS2 结果 | 0 无初检发现 / 2 有待核查或多记录安全停止 / 1 错误 / 3 MITOS2 环境门禁未满足 |
| 多候选多工具注释：`run_candidate_annotations.py` | 单记录候选或 assembly manifest；MITOS2、MitoZ；可选 `--split-multi-records` 将 contig 分别注释 | contig 级注释只提供局部基因线索，不合并 contig、不代表完整组装；注释结果不自动择优 | 原始工具输出 + `annotation_candidates.json`；没有合格序列或输出时非 0/明确记录 |
| 基因定位：`blast_genes.py` | 已注释参考 GB + 目标 FASTA；需 blastn | 默认 identity ≥80% 且 coverage ≥80%；独立定位分开，HSP 仅无歧义共线时合并，identity 按比对列加权；不独自判重排或丢失 | 1-based 定位（TSV 第 8 列为 target）；任一无唯一命中返回 1，不写部分结果 |
| reads 检查：`depth_analysis.py` | 排序并建索引 BAM + 匹配 FASTA；覆盖、soft-clip、碱基支持 | 单 mt 参考不排除 NUMT；callable 的工程检查不替代全部科学验收 | 0 无覆盖预警 / 1 输入或 BAM 错误 / 2 低覆盖或零覆盖位点；可加 `--output-json` |
| 两 scaffold 候选：`circularize.py` | 两 FASTA、参考 GB；首次可无 BAM 生成候选，回贴后验证需 BAM；需 blastn/minimap2/samtools/Biopython | 只自动验证一个内部接缝，不验证尾首闭合；方向/顺序与参考冲突时拒绝 accept | 0 CANDIDATE_ACCEPTED / 2 REVIEW（含待回贴）/ 1 失败或接缝不足；0 也非物理闭环证明 |
| COX1 线索：`cox1_id.py` | 本地确定 `--coords start,end`，片段 400–5000 bp；须有上传授权再传 `--allow-public-upload` | 向 NCBI 上传片段；不自动定位 COX1，identity/coverage 不等于物种鉴定 | 0 得判读 / 1 insufficient 或 no_match / 2 拒绝 / 3 网络或格式故障 |
**HSP 坐标与数值边界（格式故障 = 退出码 3，不是 `no_match`）**：坐标**越界**（`query-from/to` 超出
`query-len`、`hit-from/to` 超出 `Hit_len`）、`query-len`/`Hit_len` 声明为 `<= 0`、非有限 `bit-score`
（`nan`/`inf`）、`align-len` 小于 query 跨度、以及 `*-strand`/`*-frame` 与坐标方向矛盾，都判为
**格式故障**（退出码 3）——同一份回复内部自相矛盾时，其中任何 HSP 都不可信，不得降级成 `no_match`。
长度字段**缺失**不报错，但会把 `ranges_checked=false` 与 `unchecked` 原因写进 `--output-json`，
使"未校验"可见而不是默认正确。`identity == 0` 且 `align-len > 0` 属**弱到无用的命中**：记为 blocker
（`zero_identity`）、不汇总 identity、不参与自动择优，但**不**退出 3，原始 HSP 始终保留在
`--output-json` 中供人工核对。


所有 Python 工具位于 `scripts/`。注释选项与例外接受条件见
[annotation_quality.md](annotation_quality.md)；解析细节与测试见
[developer-contract.md](developer-contract.md) §2/§8/§9。

`circularize.py` 两阶段使用：首次生成候选后回贴 reads，按实际候选构建 BAM；
再次运行给 `--bam <candidate.bam> --reads-validated`。标记不替代证据。
`--accept-candidate` 仅在已核对全部新增接缝、重复歧义与组装图且满足采纳授权后使用；
当前工具提示要求人工核对，不因退出 0 而宣称全部连接已自动验证。

`depth_analysis.py BAM FASTA [AMBIGUOUS_FASTA] --output-json FILE` 保存覆盖与位点支持。
额外 FASTA 与参考必须同名、同长度，明确碱基不能改变；其他扩展名也可用。
`--min-depth 5 --min-mapq 20 --min-baseq 20` 为工程默认值。重叠 mate 仅计一个模板，
冲突、未知质量、MAPQ=255 不计定点支持；覆盖也先显式排除未知质量和 MAPQ。
全零覆盖不能通过。samtools 需支持 `view -e` 及质量表达式；不支持时明确失败。
`--allow-base-replacement` 只影响候选建议，不写修复 FASTA，也不授权采纳；soft-clip 仅报观测。

COX1 多 HSP 出现 query/target 重叠、非共线或混链时不自动汇总择优；
这不等于命中已被证明无效。查看逐 HSP 证据，使用 `--output-json <文件>` 保存详情。
环状参考跨原点候选另需坐标与结构检查。

## 案例、参考与运行支持

| 工具 | 用途 | 输出与限制 |
|---|---|---|
| 历史 `case-*` / lesson / sync 命令 | 当前分发缺少 `experience.py`，均不可执行 | 任务内记录改用 case.py；旧命令不映射为新命令，详见 task-records.md |
| `scripts/reference_registry.py` | acquire/register/record-database/list/verify | 固定版本、hash 与用途；不判参考是否可靠，见 [reference-policy.md](reference-policy.md) |
| 依赖盘点/门禁：`python3 tools/env_check.py --setup|--daily|--stage NAME [--json --strict --network --path DIR --lock PATH]` | skill 目录 + `config/dependencies.json`（无需网络，`--network` 才探测联网依赖） | 三级依赖（essential/extended/optional）是否就绪、某个 stage 能否执行、环境相对上次是否变化 | **不验证数据正确性**；lock 是缓存不是信任凭证（日常仍真实探测） | `0` 就绪 / `2` essential 缺失 / `3` 该 stage 缺依赖；写 `$MITO_KNOWLEDGE_DIR/environment.lock.json` |
| `bash scripts/check_env.sh` | 全环境盘点，可选 | 0 核心依赖就绪 / 2 核心依赖缺失；只阻断实际依赖缺失工具的步骤 |
| `python3 scripts/task_manager.py start --task-root <样本目录>/intermediate/tasks --name <名> -- <命令...>` | 脱离会话启动长任务 | 任务目录内保存命令、日志、PID 与最终退出码；进程以参数数组启动 |
| `python3 scripts/task_manager.py list --task-root <样本目录>/intermediate/tasks` | 列出任务及当前状态 | 若管理进程丢失但任务进程组仍活动，标记 `orphaned_running`；否则标记 `supervisor_lost`，退出码未知；两者都不能当作科学任务成功 |
| `python3 scripts/task_manager.py status --task-dir <任务目录>` / `log --task-dir <任务目录>` | 查询状态与日志 | 运行中不等于完成；成功退出也不等于科学验收 |
| `bash scripts/run_mitos2.sh -i <FASTA> -o <目录> -c <已确认密码表>` | 注释；环境指定 Python 与数据库 | 自动建输出目录，不生成 GB；需要时用格式桥转换 |
| `bash scripts/run_circular_map.sh <GB> <PNG> --title <标题>` | 正负链基因图 | 需要 GB 声明 circular、Biopython 与 matplotlib；声明/绘图不证明环化 |
| 外部 `table2asn` | NCBI 提交预检，需 template | 不在 check_env 范围；通过不等于已被 NCBI 接受 |

## 执行原则

- 只做能改变判断或属于验收必要条件的检查；无关工具缺失不要求安装。
- 原始数据只读，候选与产物另存；工具路径来自环境配置，不写回公共技能文件。
- 阈值是工程预警，不能包装成 NCBI 硬要求或类群通则。
- 退出码 2 的含义因工具不同，按本表解释；失败/未执行不能写成阴性科学结果。
- **退出码 3 = 依赖/环境故障**（缺必需依赖、网络或结果格式故障）：该步骤**未执行**，不是科学结论。`scripts/*.py` 的前置门禁、`env_check.py --stage`、以及缺 `jsonschema` 的案例命令都用 3。
- `baseline_report.md` 是历史审计快照，不作为现行接口依据。

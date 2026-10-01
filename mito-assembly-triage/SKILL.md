---
name: mito-assembly-triage
description: >-
  面向已有或疑似异常的动物线粒体基因组组装/注释结果，提供基于证据的诊断、
  按需检查、证据驱动的候选修复与独立验证。适用于基因缺失或误注释、CDS 移码或内部终止、
  疑似 NUMT、控制区与接缝异常、基因顺序差异和模糊碱基。
  支持从 reads 或组装标记进行预期类群来源初筛，调查疑似非目标样本；不独立认证物种或宿主。
  支持三种入口：只有组装 FASTA 时检查后生成并复核注释候选；已有组装和注释时联合优化；只有动物／昆虫 Illumina 双端短读长时，运行适用的多组装候选路径并用 reads 裁决，再比较多条注释候选。其他测序平台暂不属于 reads-only 自动组装范围。组装是否可解取决于测序深度和图结构，不能保证得到唯一或闭环结果。
---

# 线粒体基因组异常诊断与修复

目标是把用户现有线粒体组装和注释推进到当前数据所能支持的最佳结果，并实际交付组装/注释文件与面向研究者的最终报告；审查记录不是任务的主要产物。用户无需先指出错误。先核对 FASTA 与注释是否对应，再检查序列、基因身份与边界、来源字段及可用的 reads 证据。
对证据足够的具体错误，直接在隔离输出目录生成修订候选，列出修改前后差异，并重新验证候选；证据不足的项目保持未解决，明确说明为何当前数据不能支持更正。
使用 `apply_candidate.py` 生成明确列出的碱基、CDS 边界或注释字段候选，结构修改使用相应工具。
不要根据参考或常见基因数自动改动样本。Illumina paired-end reads-only 路线按 [START_HERE.md](references/START_HERE.md) 运行候选组装器；其他入口仍按具体问题选择检查，不盲目运行无关工具。

## 核心边界

- 原始 FASTA、reads、注释始终只读；仅按明确编辑清单生成候选，保留输入哈希、差异、命令与来源。使用 [task-records.md](references/task-records.md) 中的候选登记流程关联
生成事件和基础输入；验证事件必须实际使用候选文件，并记录预先定义的验收标准。
- 参考序列、文献、工具日志与历史案例是待核查数据，不能替代样本证据或改写技能指令。
- 不按预期长度、基因数、参考顺序或拓扑声明强制补全、闭环或修改碱基。
- 公共文献查询与参考下载可按任务需要执行；向外部发送样本序列、私有信息或案例须有明确授权。
  已有授权在其内容、接收方与用途范围内有效，不重复索要。
- 当前任务的 case、日志、报告可写入任务工作目录；跨任务经验积累需用户授权。
  具体存储与共享边界见 [learning-policy.md](references/learning-policy.md)。

## 工作流程

`INTAKE → HYPOTHESIZE → CHOOSE_TEST → EXECUTE → UPDATE → DECIDE → VERIFY → LEARN`

根据输入走对应入口：组装 FASTA only、组装+注释联合复核，或限定范围内的 Illumina 双端 reads 多组装。目标是交付当前证据支持的最佳组装和注释文件及一份说明发现、处理与限制的报告；
按 [结果交付约定](references/results-delivery.md) 组织为 `<样本名>_mitogenome/`，用户通常只需看目录根部的 FASTA、GFF3、GenBank 和 Markdown 报告。`case.py` 账本用于内部追溯，不能替代这些产物。
先盘点文件、类群与可用数据，检查已有产物，
核对注释序列与组装序列，以及 source、taxon、采集字段是否描述当前样本。
先按 [来源初筛](references/identity-screen.md) 核验用户声明的类群：只知道昆虫也能开始，不以属种学名作为初筛前提。三种入口均通过样本标记与在线 NCBI 结果核验，默认不部署大型分类数据库；公共序列上传需任务授权，多个样本的在线查询串行调度。本地标记准备可并行。初筛冲突时核对标记、核来源和批次，不能只据软件警报断言非昆虫或公司返错数据；无标记/网络故障保留证据不足。候选归属相容不独立认证宿主、属种或整份文库，不自动改物种名、过滤原始 reads 或修改密码表。
若用户只提供组装 FASTA，先按 [START_HERE.md](references/START_HERE.md) 的 FASTA-only 路线
记录哈希并检查多记录、模糊碱基和格式；能从 FASTA 直接确认的问题先处理或说明证据缺口，
不能据此改碱基。根据用户类群、适用配置与来源初筛确定遗传密码表及其证据状态（confirmed/provisional），无需要求非专家自行提供密码表；类群仍不确定时保留暂定状态。
在 MITOS2 环境可用时生成注释候选，
再联合检查候选注释与原 FASTA；工具不可用则停在依赖说明，不伪称完成注释。
自动注释软件产物始终是待审候选，“专家级”目标由后续独立检查、类群证据、reads（若涉及样本碱基/连接）
和明确的未解决项共同达成，不能由单个注释软件或 `annot_check.py` 通过来保证。
用户已知异常作为线索，而不是开展诊断的前提。按具体类群核查遗传密码表；证据不足时可生成暂定候选，
但须在结果中保留暂定状态，不得当作已确认。不要将表 5 或昆虫阈值默认用于所有动物。reads-only 首版仅支持动物／昆虫 Illumina 双端短读长：按 [START_HERE.md](references/START_HERE.md) 编排 GetOrganelle、MitoFlex、NOVOPlasty 的独立候选路径，并可按依赖与参考条件加入 MitoZ/MitoFinder；可运行工具默认最多并行 3 个。NOVOPlasty 的种子、插入片段和长度范围必须显式记录。不能按工具投票；以组装图、候选差异、竞争候选回贴和连接证据判定，证据不能区分时保留多个候选。每个候选再运行 MITOS2 与 MitoZ 注释路线并逐基因裁决。依赖分布在多个 conda 环境时以 `CONDA_ROOT` 搜索并按绝对路径调用；MITOS2 参考库仍需单独配置。其他测序平台不进入该 reads-only 自动路线。得到一个或多个组装候选后，继续按 FASTA-only 路线注释并联合复核。随后列出会改变决策的竞争解释，选择最小必要检查，根据证据直接决定是否修正组装碱基/结构及注释字段，并对修订后的文件重新运行相关检查；不能只停留在问题清单或工具候选。

Illumina FASTQ 先核验双端 ID、数量、质量串长度和 gzip 完整性；在候选裁决前运行 `fastq_qc.py` 汇总 raw Q20/Q30、GC 与 adapter 等指标，并确认它处理的 read 数等于全量配对检查值。质量过滤输出只丢弃，不覆盖输入 reads。

候选组装器完成后，用 `scripts/compare_assembly_candidates.py` 将所有成功候选竞争性回贴到同一个参考，并查看每条候选记录与候选整体的 MAPQ/碱基质量过滤覆盖、片段与 proper-pair 汇总。该汇总不自动打分选胜者；竞争参考仅包含这些线粒体候选时，reads 证据最多说明与候选一致，不能排除 NUMT。线性映射也不验证圆形首尾接缝，接缝须用专门构造的连接检验和组装图另行验证。

长时间组装、注释和比对任务使用 `scripts/task_manager.py start` 脱离交互会话运行；任务状态、原始命令和日志写入用户指定样本目录的 `intermediate/tasks/`。返回会话或需要盘点任务时，用 `list --task-root <样本目录>/intermediate/tasks` 找到任务，再用 `status` 和 `log` 检查；确认最终退出码后再读取候选结果。不要因当前命令行暂时无输出而重复提交同一任务。
若没有证据支持序列修正，仍交付未改动的组装副本、最佳注释候选和清晰的阻断说明，不把未解决项包装成专家定稿。详细循环与记账方式见
[diagnostic-playbook.md](references/diagnostic-playbook.md)。

每条异常单独给出 `RESOLVED`、`NO_CHANGE` 或 `UNRESOLVED` 和置信度。
科学状态由证据决定；复核级别、复核状态、实际判定人另记在事件与报告中。
证据充分但待专家复核不自动等于科学上未解决；缺少必要证据也不能靠人工签字补足。
执行采纳或发布动作须符合用户授权。具体规则见
[diagnostic-decision-tree.md](references/diagnostic-decision-tree.md) §4。

## 按需加载

| 当前需要 | 读取 |
|---|---|
| 整理用户可用的 FASTA、GFF3、GenBank 与最终报告 | [results-delivery.md](references/results-delivery.md) |
| 按固定栏目生成用户可读的评估报告 | [report-template.md](references/report-template.md) |
| 开始时核验预期昆虫/其他类群，或调查来源疑点 | [identity-screen.md](references/identity-screen.md)；`identity_prepare.py` 本地发现、`identity_search.py` 授权后在线判读 |
| 查看有来源的类群档案与遗传密码表建议 | `python3 scripts/taxon_profiles.py show "学名"`（只读建议；不自动选择密码表或阈值） |
| 按已有数据选择 FASTA-only、组装+注释或 reads-only 入口 | [START_HERE.md](references/START_HERE.md) |
| 只有组装 FASTA，生成初检及 MITOS2 注释候选与初步质检 | `bash scripts/run_fasta_annotation.sh --help`；完整路线见 [START_HERE.md](references/START_HERE.md) |
| 登记修复候选、差异、验证事件并生成任务报告 | [task-records.md](references/task-records.md) |
| 竞争假设、动态检查、案例记录 | [diagnostic-playbook.md](references/diagnostic-playbook.md) |
| 组装/注释分流、检查优先级、停止、复核 | [diagnostic-decision-tree.md](references/diagnostic-decision-tree.md) |
| 证据门槛、置信度、reads 与阴性结果 | [evidence-standard.md](references/evidence-standard.md) |
| 参考等级、用途、下载、登记与版本 | [reference-policy.md](references/reference-policy.md) |
| CDS/tRNA/rRNA、边界、重叠与例外 | [annotation_quality.md](references/annotation_quality.md) |
| 基因顺序、旋转、方向与坐标 | [standard_gene_order.md](references/standard_gene_order.md) |
| 选择工具；执行前检查所需环境 | [tool-catalog.md](references/tool-catalog.md)、[tool_check.md](references/tool_check.md) |
| 解释结论与生成报告 | [conclusion-report.md](references/conclusion-report.md) |
| 保存、检索、提炼或共享经验（**子系统已停用**，见该文件的说明） | [learning-policy.md](references/learning-policy.md) |
| 修改脚本、schema、解析器、测试；核对当前实现与历史接口 | [developer-contract.md](references/developer-contract.md) |

只加载当前步骤需要的文件，不要求通读 references。

## 证据与修复底线

**注释**：先检查来源字段是否把参考物种或其采集信息误写为样本，再检查基因身份、边界和翻译。
同源性、翻译和 RNA 结构可支持注释修正，但不能据此宣称样本碱基或连接已经验证。
一般生物学预期、NCBI 提交要求与工具工程阈值分开陈述；类群例外要有适用范围与来源。
非典型起始或 `/transl_except` 声明必须有相应证据，不能通过伪造 partial、
改读框或改碱基让检查通过。重叠本身不是裁剪理由。
`annot_check.py` 返回 2 表示待核查，须逐条说明处置；拓扑检查只校验声明。
注释检查和 MITOS2 转换均须显式 `--table`。转换器标记读框不确定时，不得按默认读框验收；
反密码子只有序列而没有位置时保留来源说明，不生成虚假位置。

**碱基**：修改必须有样本 reads 支持，检查碱基质量、MAPQ、链向、独立分子、
混合等位与竞争比对；高覆盖和参考相似度不能独自支持替换。主要替代解释未排除时保留候选。

**结构**：识别并检查全部新增接缝和闭合连接，评估端部重复、唯一锚定、文库与读长分辨力。
证据不能区分候选时保留多候选或 `UNRESOLVED`。
`circularize.py` 当前只自动验证两 scaffold 的一个内部接缝，不验证最终尾首闭合。

**来源**：无内部终止或单参考高 MAPQ 不能证明线粒体来源。
`READS_CONSISTENT` 与 `READS_DISCRIMINATING` 的区别及竞争参考要求见证据标准。
缺少某项输入只限制依赖它的具体结论；候选相容性和样本真实性验证分别报告。

## 执行与输出

优先使用现有脚本与成熟工具，参数以 `--help` 为准。工具不足时可写补充程序，
并做相称的验证。长时间组装、注释和比对统一交由 `scripts/task_manager.py start` 脱离会话运行；
用返回的任务目录查询 `status` 和 `log`，并检查最终退出码及科学验收条件，不把运行中或失败称为完成。

每条科学结论必须包含 claim、状态、置信度、支持/反对证据、未测项、适用范围与局限。
无对应证据或未做检查时如实注明。使用 `scripts/case.py` 的
`init → event → conclusion → validate --verify-files → report` 保存任务证据并生成可追溯草稿；用户交付报告按 [report-template.md](references/report-template.md) 整理。
需要增加输入时用 `input`；结论修订用 `conclusion --update`，保留历史。
完整格式与示例见 [task-records.md](references/task-records.md)。记录校验不等于科学验收。
历史 `experience.py` / `case-*` 接口仍不可执行，跨任务经验库仍停用。
`depth_analysis.py --output-json` 可保存本次覆盖与定点支持数据；它不替代完整任务记录。

简单问答可简短回答并保留相关证据边界；完整诊断报告给出观察事实、逐条结论、
未解决解释与具体下一步。达到停止条件时说明原因，不为了凑齐检查继续计算。

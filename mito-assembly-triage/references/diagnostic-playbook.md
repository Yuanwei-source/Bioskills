# 动态诊断手册

`INTAKE → HYPOTHESIZE → CHOOSE_TEST → EXECUTE → UPDATE → DECIDE → VERIFY → LEARN`

按用户目标处理一个或多个异常，不默认重组装。入口见
[START_HERE.md](START_HERE.md)，工具选择见 [tool-catalog.md](tool-catalog.md)。

## INTAKE

记录具体异常、坐标与观察来源，区分事实和用户猜测。确认类群、遗传密码表、
组装状态、可用 FASTA/GB/GFF/FASTQ/BAM/GFA、软件/数据库版本与输入哈希。
检查已有产物能否复用，了解资源限制与上传授权。若需参考先记用途，取得时再登记。

当前任务的结构化案例可写入任务工作目录；长期知识库另按
[learning-policy.md](learning-policy.md) 授权。不要把输出写进 skill 安装目录。
案例类型使用 `abnormal_case`、`normal_validation_case` 或 `tool_failure_case`；
正常案例必须明确检查范围，工具失败不得当作生物学异常。

## HYPOTHESIZE

列出当前证据相容且会影响决策的解释，以及各自的预期观测与反证；不按固定数量凑假设。
基因未检出需考虑命名/漏注释/检索灵敏度/断裂/真实丢失；内部 stop 需考虑密码表、
边界/读框、碱基错误和有证据的生物例外；结构异常需考虑表示方式、重复与错接。
将竞争解释写入事件 `summary` 与结论 `rationale`；支持、反对与未知分别记录。

## CHOOSE_TEST

按 [diagnostic-decision-tree.md](diagnostic-decision-tree.md) 选择当前最有判别力的检查。
优先低成本且能改变结论的检查，修复验收必需项不得省略。参考等级与用途按
[reference-policy.md](reference-policy.md) §5；只检查所用工具的依赖。
可检索已授权知识库的适用经验，但经验频次不能替代当前证据。

## EXECUTE

优先现有脚本与成熟工具；补充代码需相称验证。
记录完整命令、参数、版本、输入 hash、输出路径、真实退出状态和日志。
原件只读，候选另存。公共资料查询/下载可按任务执行，发送样本或私有信息须明确授权。

## UPDATE

对每项检查记录支持、反对、无法区分或工具故障，更新竞争解释；新证据可改变原分流。
按产生机制判断证据独立性，不按软件数量计数，见
[evidence-standard.md](evidence-standard.md) §4。

## DECIDE

每条异常单独记录 claim、状态、置信度与 reads 支持，存入 `conclusions[]`，不生成案例级总置信度。
`NO_CHANGE` 仅表示处置不修改；注明理由，不能据此推出样本质量通过。
置信度与缺失输入的规则统一见证据标准 §2；复核与判定人的记录见决策树 §4。

`annot_check.py` 的待核查项逐条说明证据与处置；INFO 项保留日志即可，
除非影响本次结论。结论需要的支持/反证、未测项、范围与局限即使无专用 CLI 字段，
也必须写入事件及报告。四类停止条件只在决策树 §3 定义，不为收尾强行给肯定结论。

### 任务内记录（当前可用方式）

使用 [task-records.md](task-records.md) 中的 `scripts/case.py` 建立任务，登记输入、
实际执行事件和逐条结论，再校验文件并生成报告。事件与结论同存于版本 3 的 `case.json`。
`depth_analysis.py --output-json` 的产物可作为事件 artifact 登记。
遗传密码表及确认依据写入事件 summary；修改差异和修改后验证分别登记产物与事件。
碱基和 CDS 边界改动可按证据清单用 `apply_candidate.py` 生成候选；结构改动使用对应结构工具。
候选另存并用 `case.py candidate` 关联基础输入和生成事件；验证事件必须实际使用候选，
写明预先定义的验收标准。此工具不自动生成、采纳修复，也不维护跨任务经验库。

## VERIFY

用修复前定义的验收标准分别检查注释、碱基或结构修改，不能用候选来源证明候选正确。
注释侧返回 0，或返回 2 且待核查项已按证据逐条处理，均不自动验证 reads/结构。
每个新增接缝和尾首闭合另需证据，拓扑声明不替代物理验证。

## LEARN

任务报告可记录成功、失败和未解决事项。获准后才把案例纳入跨任务经验库、提出 lesson，
并说明它支持还是挑战既有经验。重复判断不提高可信度，推广和共享见经验政策。

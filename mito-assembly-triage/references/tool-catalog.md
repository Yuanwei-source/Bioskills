# 工具目录（按需调用）

先核对已有产物是否可复用，再选择当前检查所需工具。记录版本、完整命令、参数、输入
hash、产物路径和实际退出状态。参数以各工具 `--help` 为准，依赖见
[tool_check.md](tool_check.md)。以下路径相对 skill 根目录；实际输出放在用户任务目录。

## 诊断与候选工具

| 任务 / 工具 | 必要输入与用途 | 关键限制 | 输出与退出码 |
|---|---|---|---|
| 序列体检：`seq_stats.py` | FASTA；长度、GC/AT、模糊碱基和滑窗 | 模糊位置是 0-based；不验证样本真实性 | 统计输出；正常 0 |
| 注释质检：`annot_check.py` | GenBank；显式传类群确认的 `--table`，必要时 `--taxon`、`--ref` | 无 reads 验证；`--taxon` 不自动切换昆虫阈值，拓扑检查只校声明 | 0 无发现 / 1 错误（含参数错误）/ 2 待核查 |
| 注释格式桥：`mitos2_to_genbank.py` | MITOS2 GFF/FAS/FAA 与组装 FASTA；生成可质检 GB | 非通用 GFF 转换；不补基因/碱基，不导出未知 partial，不是提交质量记录 | 写 GB；输入/坐标异常非 0，无半成品 |
| 基因定位：`blast_genes.py` | 已注释参考 GB + 目标 FASTA；需 blastn | 默认 identity ≥80% 且 coverage ≥80%；短/局部命中被排除，不独自判重排或丢失 | 1-based 定位；任一无唯一命中返回 1，不写部分结果 |
| reads 检查：`depth_analysis.py` | 排序并建索引 BAM + 匹配 FASTA；覆盖、soft-clip、碱基支持 | 单 mt 参考不排除 NUMT；callable 的工程检查不替代全部科学验收 | 0 正常 / 1 输入或 BAM 错误 / 2 低覆盖 |
| 两 scaffold 候选：`circularize.py` | 两 FASTA、参考 GB；首次可无 BAM 生成候选，回贴后验证需 BAM；需 blastn/minimap2/samtools/Biopython | 只自动验证一个内部接缝，不验证尾首闭合；方向/顺序与参考冲突时拒绝 accept | 0 CANDIDATE_ACCEPTED / 2 REVIEW（含待回贴）/ 1 失败或接缝不足；0 也非物理闭环证明 |
| COX1 线索：`cox1_id.py` | 本地确定 `--coords start,end`，片段 400–5000 bp；须有上传授权再传 `--allow-public-upload` | 向 NCBI 上传片段；不自动定位 COX1，identity/coverage 不等于物种鉴定 | 0 得判读 / 1 insufficient 或 no_match / 2 拒绝 / 3 网络或格式故障 |

所有 Python 工具位于 `scripts/`。注释选项与例外接受条件见
[annotation_quality.md](annotation_quality.md)；解析细节与测试见
[developer-contract.md](developer-contract.md) §2/§8/§9。

`circularize.py` 两阶段使用：首次生成候选后回贴 reads，按实际候选构建 BAM；
再次运行给 `--bam <candidate.bam> --reads-validated`。标记不替代证据。
`--accept-candidate` 仅在已核对全部新增接缝、重复歧义与组装图且满足采纳授权后使用；
当前工具提示要求人工核对，不因退出 0 而宣称全部连接已自动验证。

`depth_analysis.py` 的可选模糊碱基支持检查需要额外的 `.fasta` 输入（与 BAM 参考坐标一致）；
`--allow-base-replacement` 影响替换建议，不自动写入修复 FASTA。其控制区 soft-clip 提示
仍有“正常长度异质性”等过强措辞，只作为待查解释，不能原样当成结论。

COX1 多 HSP 出现 query/target 重叠、非共线或混链时不自动汇总择优；
这不等于命中已被证明无效。查看逐 HSP 证据，使用 `--output-json <文件>` 保存详情。
环状参考跨原点候选另需坐标与结构检查。

## 案例、参考与运行支持

| 工具 | 用途 | 输出与限制 |
|---|---|---|
| `tools/experience.py case-init/case-event/case-anomaly/case-reference/case-validate/case-report` | 案例、事件、逐异常记录、参考关联、记录校验（`case-validate` 只校记录格式，不证明科学结论）；case-init 必须有 hypothesis；缺 `jsonschema` 时校验/写入口以**退出码 3** 失败（提示 pip 安装，不降级） | 写指定任务目录；例子见 [diagnostic-playbook.md](diagnostic-playbook.md)、确切参数见 [developer-contract.md](developer-contract.md) §3.1 |
| `case-validate / case-report` | 格式校验 / 从已有记录生成 case.md | 校验失败 1；报告是草稿，按状态分类有缺陷，见 [conclusion-report.md](conclusion-report.md) §5 |
| `search-structured / propose-lesson / review-lesson` | 检索、提炼、审核经验 | 写跨任务知识库前需授权；verified 不是通用规则 |
| `export-contribution / sync-public` | 预览脱敏贡献 / 同步隔离缓存 | 生成预览与实际上传分开；授权与校验见 [learning-policy.md](learning-policy.md) |
| `scripts/reference_registry.py` | acquire/register/record-database/list/verify | 固定版本、hash 与用途；不判参考是否可靠，见 [reference-policy.md](reference-policy.md) |
| 代码↔文档审计：`python3 tools/audit_code_docs.py [--json] [--check-baseline] [--update-baseline]` | skill 目录本身（无需网络） | 代码里的 flag/错误码/枚举值是否仍有文档、规则表是否仍可检索、基线是否同步 | **不判断文档写得对不对**（语义仍须人工评审）；`--update-baseline` 只写豁免清单 | 非 0 = 发现漂移；写 `tools/doc_contract_baseline.json` |
| `bash scripts/check_env.sh` | 全环境盘点，可选 | 0 核心依赖就绪 / 2 核心依赖缺失；只阻断实际依赖缺失工具的步骤 |
| `bash scripts/run_bg.sh <名> -- <命令...>` | 长任务后台运行 | logs 中记录 log/pid/status；可信 shell 才用 `--trusted-shell` |
| `bash scripts/check_bg.sh <名> [--tail N]` | 查询实际状态与日志 | running/timeout 不等于完成；正常退出不等于科学验收 |
| `bash scripts/run_mitos2.sh -i <FASTA> -o <目录> -c <已确认密码表>` | 注释；环境指定 Python 与数据库 | 自动建输出目录，不生成 GB；需要时用格式桥转换 |
| `bash scripts/run_circular_map.sh <GB> <PNG> --title <标题>` | 正负链基因图 | 需要 GB 声明 circular、Biopython 与 matplotlib；声明/绘图不证明环化 |
| 外部 `table2asn` | NCBI 提交预检，需 template | 不在 check_env 范围；通过不等于已被 NCBI 接受 |

## 执行原则

- 只做能改变判断或属于验收必要条件的检查；无关工具缺失不要求安装。
- 原始数据只读，候选与产物另存；工具路径来自环境配置，不写回公共技能文件。
- 阈值是工程预警，不能包装成 NCBI 硬要求或类群通则。
- 退出码 2 的含义因工具不同，按本表解释；失败/未执行不能写成阴性科学结果。
- `baseline_report.md` 是历史审计快照，不作为现行接口依据。

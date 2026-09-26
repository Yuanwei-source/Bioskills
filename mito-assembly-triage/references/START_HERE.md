# 从这里开始

按当前数据和异常选入口，只读取相关资料。详细循环见
[diagnostic-playbook.md](diagnostic-playbook.md)，不要默认运行全部工具。

## 1. 按输入选路线

| 已有数据 | 入口动作 | 证据边界 |
|---|---|---|
| FASTA | `python3 scripts/seq_stats.py <assembly.fasta>` | 可统计序列属性，不能验证样本碱基与物理闭环 |
| FASTA + GenBank | 按类群确认密码表，运行 `annot_check.py` | 注释一致性不等于样本序列正确 |
| MITOS2 输出 | 用 `mitos2_to_genbank.py` 转换后质检 | 此转换器针对 MITOS2 文件组合，不是通用 GFF 转换器；拓扑仍只是声明 |
| FASTA + FASTQ/BAM | 按读长与文库比对、排序、建索引，记录重复标记/处理；定点检查 reads | 先确认 BAM 与候选参考匹配；单参考不能排除 NUMT |
| reads + 核/其他竞争参考 | 比较候选、核位点或竞争路径 | 存在竞争参考不自动获得判别力，须检查实际结果 |
| 多 contig / GFA | 先看图与连接假设；两 scaffold 场景可用 `circularize.py` | 脚本只自动验证一个内部接缝，闭合连接另需验证 |
| 需要公共参考 | 先看本地资源；按 [reference-policy.md](reference-policy.md) §5 选等级和用途，获取后登记 | 参考是对照；序列下载与样本上传分开授权 |

先核对已有产物的输入、版本和用途；兼容的产物可复用，条件变化时才重跑。

## 2. 按异常找资料

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

## 3. 条件式首次检查

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
输出前阅读 [conclusion-report.md](conclusion-report.md)：自动报告须复核后使用。

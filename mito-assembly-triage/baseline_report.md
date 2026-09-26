# mito-assembly-triage V1 冻结审计

日期：2026-09-24
基线提交：`749ce2d`；冻结标签：`v1-fixed`
后续演进：在同一 Skill 与同一 Git 历史中增量实现

## 当前 V1 已确认的修复

- `depth_analysis.py` 在低覆盖时仍会完成 mate、soft-clip 和模糊碱基诊断；reverse FLAG 只做一次方向转换。
- `seq_stats.py` 按 FASTA contig 统计，GC/AT 使用 A/C/G/T 分母。
- `blast_genes.py` 保留 target ID、链方向，并对同一 target/strand 的 query HSP 做 union。
- `circularize.py` 做端部精确重叠裁剪和 PAF 区间 union；默认输出候选而非宣称最终环化。
- `annot_check.py` 修正 16S/12S 到 rrnL/rrnS 的映射和嵌套 overlap 计算。
- `cox1_id.py` 不再用固定 DNA 模式宣称 COX1 定位；无可靠坐标时停止并要求 REVIEW。
- 后台任务使用锁和原子状态文件；环境检查在核心依赖缺失时返回非零。
- `experience.py` 新案例默认标记“待审核”，信号统计使用多行匹配。

## V1 测试与环境

使用 `<MITOS2 环境>/bin/python`（Biopython 1.87）运行（原始记录里的本地绝对路径已改为占位符）：

```text
Ran 26 tests
OK
```

核心可用性检查发现：`blastn`、`makeblastdb`、`samtools`、`minimap2` 可用；当前 generic 配置未找到 `bwa`、`seqkit`、GetOrganelle、MitoFinder、MITOS2 依赖和 MITOS2 Python 配置。该结果只影响依赖这些工具的路径，不影响 FASTA/注释纯 Python 检查。

## V1 基线与后续增量的差距

1. `SKILL.md` 仍偏向静态流程和质量门，缺少 INTAKE → HYPOTHESIZE → CHOOSE_TEST → EXECUTE → UPDATE → DECIDE → VERIFY → LEARN 的最小诊断循环。
2. 没有结构化 `case.json`、`events.jsonl` 的契约和校验/写出工具。
3. 经验模块只有旧版 Markdown 案例、关键词检索和频次建议；缺少异常字段检索、候选 lesson、验证状态和版本化发布/撤销。
4. 没有统一报告模板来区分观察、证据来源、推断、冲突、未知、修改和验证状态。
5. 旧案例没有自动破坏性迁移需求；V2 应采用可逆映射，并将真实运行数据放在工作目录而非公共仓库。

## V1 遗留限制（不在 Phase 0 混入修复）

- 环化仍需真实 reads、全部接缝和组装图证据；当前脚本只是候选结构生成器。
- 缺少 reads 时不能报告 raw-read-supported，也不能排除 NUMT 或唯一结构。
- 依赖工具是否可用由运行环境决定；本分支不替换环境、不新增组装器。

## Phase 0 结论

V1 既有功能路径通过现有回归测试，已知修复项未发现回归。后续功能可在同一 Skill 中增量加入结构化案例/事件记录和规范化诊断报告；不需要重写现有分析脚本。

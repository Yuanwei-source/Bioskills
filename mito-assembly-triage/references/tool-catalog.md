# V1 工具目录（按需调用）

| 任务 | 工具 | 修改输入 | 局限 |
|---|---|---:|---|
| 序列体检 | `scripts/seq_stats.py` | 否 | 不证明 reads 来源 |
| 注释质检 | `scripts/annot_check.py` | 否 | 典型动物基因集合不是绝对真值 |
| 基因定位 | `scripts/blast_genes.py` + BLAST | 否 | 局部同源命中需检查重复/覆盖 |
| reads 诊断 | `scripts/depth_analysis.py` + samtools | 否 | 依赖正确比对和参考 |
| 环化候选 | `scripts/circularize.py` + minimap2 | 写候选文件 | 不是最终环化证明 |
| COX1 查询 | `scripts/cox1_id.py` | 向公共服务上传需显式同意 | identity 不是确定物种 |

运行前记录工具版本、命令、参数和输入 hash。

---
name: Gene Enrichment Analysis (ORA / GSEA)
description: >
  对基因列表或差异表达结果执行 KEGG 和 GO 富集分析。
  自动识别输入格式，根据是否提供排序指标（logFC 等）动态切换 ORA 或 GSEA 模式。
  支持多物种、自动 ID 类型探测（Symbol / Ensembl / Entrez），
  并输出散点图、山脊联合图、GSEA 曲线图、Barcode 综合图等高质量可视化结果。
---

# Skill: 智能基因富集分析

## 前提条件

底层 R 脚本路径：`z:/Project/skills/Enrichment/Enrichment.R`

必须已安装的 R 包：
- `clusterProfiler`, `enrichplot`, `ggplot2`, `ggridges`
- `dplyr`, `forcats`, `stringr`, `data.table`, `optparse`
- 对应物种的 OrgDb（见物种参考表）

---

## AI 执行 SOP（三步走，必须严格遵循）

> **禁止跳过第 1 步直接执行！列名猜错会导致整个分析崩溃。**

### Step 1 — 数据嗅探（Data Exploration）

使用文件查看工具读取用户提供输入文件的**前 10 行**，记录以下信息：

| 要找的信息       | 描述                                                                 |
| ------------ | -------------------------------------------------------------------- |
| **基因列**   | 哪一列储存基因名？（如 `Gene_Symbol`、`GeneID`、`SYMBOL`、`gene_name` 等）   |
| **排序指标列** | 是否有代表差异倍数的列？（如 `log2FoldChange`、`logFC`、`FC`）优先选 log 转换后的值 |
| **p 值列**   | 是否有 `padj`、`adj.P.Val`、`FDR` 列？ORA 模式下可用于过滤                     |
| **ID 类型** | 基因名是 Symbol（TP53）、Ensembl（ENSG...）还是纯数字（Entrez）？脚本会自动探测，但你需判断是否合理 |

**若文件只有 1 列（纯基因列表）**：直接进入 ORA 模式，没有 metric / pval。

### Step 2 — 分析模式路由（Mode Selection）

根据 Step 1 的结果决定分析模式：

```
┌─ 用户提供了含所有基因的完整差异结果表（包含 logFC）
│     → GSEA 模式：提取 metric_col，传入 --metric_col
│
└─ 用户只提供了显著基因名单（或纯基因列表 txt）
      → ORA 模式：不传 --metric_col
          ├─ 有 padj 列 → 传入 --pval_col 和 --pval_cutoff
          └─ 无 padj 列 → 不传，默认全部输入基因视为显著
```

### Step 3 — 精准调用底层脚本（Execution）

基于以上信息，组装并执行命令：

**ORA 示例（有筛选）：**
```bash
Rscript "z:/Project/skills/Enrichment/Enrichment.R" \
  --input    "/path/to/deseq2_results.csv" \
  --gene_col "Gene_Symbol" \
  --pval_col "padj" \
  --pval_cutoff 0.05 \
  --organism "ssc" \
  --orgdb    "org.Ss.eg.db" \
  --outdir   "/path/to/output/" \
  --task_id  "Group_A_vs_B"
```

**GSEA 示例（全量差异表）：**
```bash
Rscript "z:/Project/skills/Enrichment/Enrichment.R" \
  --input      "/path/to/all_genes_DE_result.csv" \
  --gene_col   "Gene_Symbol" \
  --metric_col "log2FoldChange" \
  --organism   "hsa" \
  --orgdb      "org.Hs.eg.db" \
  --outdir     "/path/to/output/" \
  --task_id    "Treatment_vs_Control"
```

**纯基因列表 ORA（单列 txt）：**
```bash
Rscript "z:/Project/skills/Enrichment/Enrichment.R" \
  --input   "/path/to/gene_list.txt" \
  --outdir  "/path/to/output/"
```

---

## 输出文件说明

### ORA 模式输出
| 文件名 | 说明 |
|---|---|
| `00.enrichment_summary.csv` | 全分析统计摘要 |
| `01.kegg_results.rds/.csv` | KEGG ORA 完整结果 |
| `02.go_results.rds/.csv` | GO ORA 完整结果（BP/CC/MF） |
| `03.kegg_scatter.pdf/.png` | KEGG 气泡散点图 |
| `04.kegg_ridge.pdf/.png` | KEGG 山脊联合图（需 logFC）|
| `05.go_scatter.pdf/.png` | GO 气泡散点图（分 Ontology 分面）|
| `06.go_ridge.pdf/.png` | GO 山脊联合图（需 logFC）|
| `07.enrichment_barcode.pdf/.png` | KEGG+GO 综合 Barcode 图 |

### GSEA 模式输出
| 文件名 | 说明 |
|---|---|
| `00.enrichment_summary.csv` | 全分析统计摘要 |
| `01.gsea_kegg_results.rds/.csv` | KEGG GSEA 完整结果 |
| `02.gsea_go_results.rds/.csv` | GO GSEA 完整结果 |
| `03.gsea_kegg_dotplot.pdf/.png` | KEGG GSEA Dotplot（Activated/Suppressed 分面）|
| `04.gsea_kegg_ridge.pdf/.png` | KEGG GSEA Ridge 分布图 |
| `05.gsea_kegg_score.pdf/.png` | KEGG Top5 通路 Enrichment Score 曲线 |
| `06.gsea_go_dotplot.pdf/.png` | GO GSEA Dotplot |
| `07.gsea_go_ridge.pdf/.png` | GO GSEA Ridge 分布图 |
| `08.gsea_go_score.pdf/.png` | GO Top5 通路 Enrichment Score 曲线 |

---

## 完整参数参考表

| 参数 | 类型 | 是否必填 | 说明 |
|---|---|---|---|
| `--input` | string | 必填 | 输入文件路径（csv/tsv/txt 均可，自动识别分隔符）|
| `--gene_col` | string | 可选 | 基因名列名，未指定则自动探测 |
| `-m`, `--metric_col` | string | 可选 | 排序指标列 (如 log2FoldChange)。提供后默认触发 GSEA，或供 ORA 画山脊图。 |
| `--mode` | string | 可选 | auto, ORA, GSEA。默认 auto (有指标则 GSEA)。选择 ORA 可强制按基因集分析。 |
| `-p`, `--pval_col` | string | 可选 | ORA 模式下用于过滤显著基因的 p 值列 (如 padj)。设定后会先过滤再分析。 |
| `--pval_cutoff` | float | 可选 | p 值阈值，默认 `0.05` |
| `--background` | string | 可选 | 背景基因文件路径（单列 txt），ORA 严谨分析必备 |
| `--organism` | string | 可选 | KEGG 物种代码，默认 `ssc` |
| `--orgdb` | string | 可选 | OrgDb R 包名，默认 `org.Ss.eg.db` |
| `--outdir` | string | 可选 | 输出目录，默认 `./enrichment_out` |
| `--task_id` | string | 可选 | 图表标题标识符，默认 `Enrichment` |
| `--kegg_p_cutoff` | float | 可选 | KEGG 富集 p 值阈值，默认 `0.05` |
| `--go_p_cutoff` | float | 可选 | GO 富集 p 值阈值，默认 `0.05` |
| `--go_q_cutoff` | float | 可选 | GO 富集 q 值阈值，默认 `0.2` |
| `--show_n` | int | 可选 | 各图表展示前 N 个 term，默认 `15` |

---

## 常用物种参考表

| 物种                  | `--organism` | `--orgdb`            |
|-----------------------|--------------|----------------------|
| 猪 (Sus scrofa)       | `ssc`        | `org.Ss.eg.db`       |
| 人 (Homo sapiens)     | `hsa`        | `org.Hs.eg.db`       |
| 小鼠 (Mus musculus)   | `mmu`        | `org.Mm.eg.db`       |
| 大鼠 (Rattus norvegicus) | `rno`     | `org.Rn.eg.db`       |
| 牛 (Bos taurus)       | `bta`        | `org.Bt.eg.db`       |
| 拟南芥 (Arabidopsis)  | `ath`        | `org.At.tair.db`     |
| 酵母 (S. cerevisiae)  | `sce`        | `org.Sc.sgd.db`      |
| 羊 (Ovis aries)       | `oas`        | `org.Oar.eg.db`      |
| 鸡 (Gallus gallus)    | `gga`        | `org.Gg.eg.db`       |
| 果蝇 (Drosophila melanogaster) | `dme` | `org.Dm.eg.db`       |
| 疟蚊 (Anopheles gambiae) | `aga`     | `org.Ag.eg.db`       |

---

## 常见问题与排坑指南

**Q: 脚本报错 "Entrez IDs found = 0"？**  
A: 最可能原因是 ID 类型不匹配。检查基因名是否为 ENSEMBL（`ENSSSCG0000...`）而非 Symbol，
或该物种的基因在 OrgDb 中覆盖率低。可尝试用 `bitr()` 手动测试几个基因。

**Q: KEGG 报错"Connection timeout"？**  
A: 脚本已配置国内镜像 `https://rest.kegg.cn/`。若仍失败，检查服务器网络访问权限。

**Q: GSEA 结果为空 / 没有显著通路？**  
A: GSEA 对基因数量要求高（通常需要 ≥ 1000 个输入基因）。确认传入的是**全体基因**（非经过过滤的显著基因）的完整差异结果表。

**Q: 山脊图（Ridge Plot）没有生成？**  
A: ORA 模式下的山脊图需要提供 `--metric_col`（logFC 列），无此列时自动跳过，这是正常的优雅降级行为。

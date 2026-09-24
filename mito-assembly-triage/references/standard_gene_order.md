# 基因顺序与坐标规范

> **适用边界**：本文件给出**比较规则**和一份**典型昆虫排列参照**；参照不是全类群通则。
> 目标：把"表示差异"与"真实结构差异"分开，避免假缺失、假重复、假重排。

## 1. 比较前先规范化

| 维度 | 处理 |
|---|---|
| 拓扑 | 明确单条/多条、完整/部分组装、是否候选环状 |
| 环状旋转 | 起点不同属于**表示差异**，不是重排 |
| 整条反向互补 | 属于**表示差异**；比较前先统一方向 |
| 真正倒位/转位/复制 | 必须**另外**用结构或 reads 证据验证，不能只凭顺序变化 |

## 2. 坐标约定（代码实际行为，改代码前先改这里）

| 场景 | 约定 |
|---|---|
| GenBank / BLAST outfmt / samtools 区间 / `blast_genes.py` 输出 | **1-based closed**（`location.start + 1` .. `location.end`） |
| 区间合并与重叠数学（`circularize.py` 的 `interval_union_length`、接缝计算） | **0-based half-open** |
| 对外导出 GFF/GenBank 前 | 必须显式转换，转换点写进注释或事件记录 |
| 跨原点基因 | 必须按环状语义处理（`join`/分段坐标），否则产生**假缺失**与**假重叠** |

保留原始坐标与规范化坐标两套映射；不要把规范化结果写回原始文件。

## 3. 检索与判定

1. 使用**不止一条**亲缘匹配且注释已核查的参考，比较**同源身份 + 链 + 邻接关系**；
2. **不能用短 BLAST 命中拼接出假重排**（`blast_genes.py` 要求 identity ≥ 80% 且覆盖 ≥ 80% 的**唯一**命中，
   多个阈值内命中直接报错而不是任选一个）；
3. 顺序差异一律先记为 **REVIEW**，不得自动失败或自动"修复"；
4. 数据无法跨越重复区、或没有唯一锚定新增连接时 → 保持候选结构并标记 `UNRESOLVED`（见 `circularize.py` 的
   `PUTATIVE_CIRCULAR/REVIEW` 输出）。

## 4. 典型昆虫排列（参照锚点，非通则）

以果蝇型（*Drosophila* 型）为代表的泛甲壳类/昆虫近祖排列。**仅作对照起点**；任何偏离都必须查类群文献 + 结构证据。

| # | 基因 | 类型 | 链 | # | 基因 | 类型 | 链 |
|---|---|---|---|---|---|---|---|
| 1 | `trnI` | tRNA | + | 20 | `trnN` | tRNA | + |
| 2 | `trnQ` | tRNA | − | 21 | `trnS1`(AGN) | tRNA | + |
| 3 | `trnM` | tRNA | + | 22 | `trnE` | tRNA | + |
| 4 | `nad2` | CDS | + | 23 | `trnF` | tRNA | − |
| 5 | `trnW` | tRNA | + | 24 | `nad5` | CDS | − |
| 6 | `trnC` | tRNA | − | 25 | `trnH` | tRNA | − |
| 7 | `trnY` | tRNA | − | 26 | `nad4` | CDS | − |
| 8 | `cox1` | CDS | + | 27 | `nad4l` | CDS | − |
| 9 | `trnL2`(UUR) | tRNA | + | 28 | `trnT` | tRNA | + |
| 10 | `cox2` | CDS | + | 29 | `trnP` | tRNA | − |
| 11 | `trnK` | tRNA | + | 30 | `nad6` | CDS | + |
| 12 | `trnD` | tRNA | + | 31 | `cob` | CDS | + |
| 13 | `atp8` | CDS | + | 32 | `trnS2`(UCN) | tRNA | + |
| 14 | `atp6` | CDS | + | 33 | `nad1` | CDS | − |
| 15 | `cox3` | CDS | + | 34 | `trnL1`(CUN) | tRNA | − |
| 16 | `trnG` | tRNA | + | 35 | `rrnL`(16S) | rRNA | − |
| 17 | `nad3` | CDS | + | 36 | `trnV` | tRNA | − |
| 18 | `trnA` | tRNA | + | 37 | `rrnS`(12S) | rRNA | − |
| 19 | `trnR` | tRNA | + | — | 控制区 CR | — | — |

自检：CDS = 13（正链 9 / 负链 4）；tRNA = 22；rRNA = 2。
`annot_check.py` 的"链分布 9+/4−"WARN 与"基因顺序"WARN 均以本类排列为默认对照。

## 5. 可变位置与常见陷阱

- `trnP` / `trnT` 在**部分**昆虫类群是重排热点（位置因类群而异）；**不是所有类群都变**，
  按类群资料启用，不要把它变成全类群硬规则。
- `trnL1`/`trnL2`、`trnS1`/`trnS2` 必须保留**反密码子 + 命名映射**；否则会造成假重复/假缺失。
- CDS 命名两套写法并存：本参照用 `nad*`/`cob`（MITOS2/NCBI 新写法），历史注释常用 `nd*`/`cytb`；
  `annot_check.py` 已归一两者，比较时不要因为它们不同就报错。
- 控制区（CR，A+T 富集）通常位于 `rrnS` 与 `trnI` 之间，但长度与重复结构高度可变，
  低覆盖 + 大量 soft-clip 往往是**正常长度异质性**，不是错误。
- `atp8-atp6`、`nad4-nad4l` 的短重叠是真实特征（见 `annotation_quality.md` §4）。

## 6. 来源

昆虫排列与重排热点结论需按类群引用文献；通用来源索引见 `evidence-standard.md` §6。

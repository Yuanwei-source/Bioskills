# 基因顺序与坐标规范

> **适用边界**：本文件给出**比较规则**、一份**典型昆虫排列参照**和一份**鳞翅目参照**；
> 参照不是全类群通则。目标：把"表示差异"与"真实结构差异"分开，避免假缺失、假重复、假重排。
> 注释侧阈值见 `annotation_quality.md`。

## 1. 比较前先规范化

| 维度 | 处理 |
|---|---|
| 拓扑 | 明确单条/多条、完整/部分、是否候选环状 |
| 环状旋转 | 起点不同属于**表示差异**，不是重排 |
| 整条反向互补 | 属于**表示差异**；比较前先统一方向（链分布检查已按此处理） |
| 真正倒位/转位/复制 | 必须**另外**用结构或 reads 证据验证，不能只凭顺序变化 |

## 2. 坐标约定（代码实际行为，改代码前先改这里）

| 场景 | 约定 |
|---|---|
| GenBank / BLAST outfmt / samtools 区间 / `blast_genes.py` 输出 / `cox1_id.py --coords` | **1-based closed** |
| `annot_check.py` 打印的坐标 | **1-based closed**；跨原点 `join()` 基因按各段分别打印 |
| 内部区间数学（`annot_check.py` 的重叠/长度、`circularize.py` 的 `interval_union_length`） | **0-based half-open** |
| `seq_stats.py` 的模糊碱基位置 | **0-based**（输出里已标注） |
| `depth_analysis.py` 的窗口与位置 | **1-based**（samtools 语义） |
| `circularize.py` 的接缝索引 `junction` | 0-based 偏移，输出区间时按 1-based closed 打印 |
| 对外导出 GFF/GenBank 前 | 必须显式转换，转换点写进注释或事件记录 |
| 跨原点基因 | 必须按**分段位置**处理，否则产生**假缺失**与**假重叠** |

## 3. 检索与判定

1. 使用**不止一条**亲缘匹配且注释已核查的参考，比较**同源身份 + 链 + 邻接关系**；
2. **不能用短 BLAST 命中拼接出假重排**（`blast_genes.py` 要求 identity ≥ 80% 且覆盖 ≥ 80% 的**唯一**命中；
   多个阈值内命中或任一基因无唯一命中时直接退出 1，**不写部分结果**）；
3. 顺序差异一律先记为 **REVIEW**，不得自动失败或自动"修复"；
4. `annot_check.py` 的顺序与长度对比**仅在提供 `--ref` 时进行**，比较对象是该参考，
   **不是**内建的排列锚点；两者不要混说。
5. 顺序比较是**有方向的环状邻接比较**，已归一环状旋转、起点与整链反向互补：
   - 共有基因的邻接关系相同 → 输出"邻接关系一致"（不因起点不同或整体 RC 误报重排）；
   - 邻接关系不同 → `ARRANGEMENT_DIFF`（REVIEW）；
   - **基因缺失/多余单独报告为 `GENE_SET_DIFF`**，不与重排混为一谈（非共有基因不进入邻接比较）。
5. 数据无法跨越重复区、或没有唯一锚定新增连接时 → 保持候选结构并标记 `UNRESOLVED`。
   `circularize.py` 在候选顺序/方向与参考不一致时记 **REVIEW**、保留候选与诊断记录并退出 2，
   且**不允许**据此 `--accept-candidate`。
   ⚠️ 该脚本目前只自动验证**一个**内部接缝（两 scaffold 场景），**不**验证最终尾部→首部的闭合连接；
   多接缝/闭环场景必须人工逐接缝取证，不能声称"全部接缝已自动验证"。

## 4. 典型昆虫排列参照（锚点，仅六足为主）

以下为果蝇型（*Drosophila* 型）排列，是**六足/泛甲壳类的近祖排列假说**，
常用作对照起点，**不是**跨类群定律。只作对照，任何偏离都必须查类群文献 + 结构证据。

线性顺序（起点取 `trnI`，旋转与整链反向互补属表示差异）：

`trnI → trnQ → trnM → nad2 → trnW → trnC → trnY → cox1 → trnL2(UUR) → cox2 → trnK → trnD → atp8 → atp6 → cox3 → trnG → nad3 → trnA → trnR → trnN → trnS1(AGN) → trnE → trnF → nad5 → trnH → nad4 → nad4l → trnT → trnP → nad6 → cob → trnS2(UCN) → nad1 → trnL1(CUN) → rrnL → trnV → rrnS`

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
参照来源与讨论见 §7；`tests/gb_fixtures.py` 的 `ORDER` 与本表逐项一致（已由回归测试锚定）。

## 5. 鳞翅目参照（与上表并存，不能互相替代）

许多鳞翅目物种的 `trnM` 位置**与昆虫近祖排列不同**，常见为 `trnM–trnI–trnQ`
（近祖为 `trnI–trnQ–trnM`），鳞翅目内部亦存在例外。因此：

- 不要用果蝇型锚点直接否掉鳞翅目的顺序；
- 鳞翅目样本的顺序差异应优先与**鳞翅目内部参考**比较（可用 `annot_check.py --ref <鳞翅目参考.gb>`）；
- 这一层背景也解释了为什么该目的 `cox1` 常见非典型起始密码子，
  参见 `annotation_quality.md` §9 与 §2.4。

## 6. 可变位置与常见陷阱

- `trnP` / `trnT` 及 tRNA 簇区（如 `trnA-trnR-trnN-trnS1-trnE-trnF`、`trnI-trnQ-trnM`）在部分类群是重排热点；
  **不是所有类群都变**，按类群资料启用，不要把它变成全类群硬规则。
- `trnL1`/`trnL2`、`trnS1`/`trnS2` 必须保留**反密码子 + 命名映射**；裸名不得自动归类（见 `annotation_quality.md` §3）。
- CDS 命名两套写法并存：本参照用 `nad*`/`cob`（MITOS2/NCBI 新写法），历史注释常用 `nd*`/`cytb`；
  `annot_check.py` 已归一两者，比较时不要因为它们不同就报错。
- 控制区（CR，A+T 富集）通常位于 `rrnS` 与 `trnI` 之间，长度与重复结构高度可变；
  低覆盖 + 大量 soft-clip 往往是**正常长度异质性**，不是错误。
- `atp8-atp6`、`nad4-nad4l` 的短重叠是真实特征；`atp8-atp6` 之间也可能是间隔区。
  两者都只触发检查（见 `annotation_quality.md` §5）。

## 7. 来源

- Cameron 2014, *Systematic Entomology* 39:400–411, DOI 10.1111/syen.12071 —
  昆虫线粒体测序/注释方法、基因排列与重排、"tRNA 数量与结构变异"。
- Boore/Bernt 等关于泛甲壳类/昆虫近祖排列与保守基因块的工作（§4 锚点作为**假说**使用）。
- 鳞翅目 `trnM` 位置差异：见 *Hyphantria cunea* 等鳞翅目 mitogenome 论文（§5）。
- 完整来源索引与引用管理要求：`evidence-standard.md` §6。

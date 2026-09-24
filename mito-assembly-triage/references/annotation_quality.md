# 注释质量标准

> **执行者**：`scripts/annot_check.py`。本文标注"代码"的条目由它执行，标注"人工"的必须人工判断。
> **证据要求**见 `evidence-standard.md`；顺序与坐标约定见 `standard_gene_order.md`。

## 0. 三类判据必须分开陈述

写结论时不要把下面三层混为一谈：

| 层 | 内容 | 性质 |
|---|---|---|
| (a) 一般生物学预期 | 典型后生动物 13 CDS + 22 tRNA + 2 rRNA；昆虫近祖基因排列 | 有类群例外，须按类群文献确认 |
| (b) NCBI 提交审查要求 | 要求提供基因/CDS 等注释，差异需向策展人说明；无统一的 bp 阈值 | 提交合规性，不是判定样本真假的标准 |
| (c) 本工具的工程预警阈值 | `>8bp` 重叠、tRNA `60–75bp`、rrnL `1100–1500` / rrnS `600–850`、正链 CDS `9+` | 本工具默认值，**不是领域公理**，可显式放宽 |

NCBI 的官方表述是：细胞器提交需提供基因/CDS 等注释，且"CDS 与 tRNA 通常很少重叠超过几个核苷酸"；
它**没有**规定 `8 bp` 这个错误阈值，也**没有**把"22 个 tRNA"写成对所有动物的硬性要求。
代码里 `>8bp` 与长度区间都是 (c) 层，不能用"(NCBI 要求…)"的口径表述。

## 1. 前置：先固定语境

| 项目 | 为什么必须 | 记录位置 |
|---|---|---|
| 类群（目/科，尽量到属） | 基因集、密码表、tRNA 结构、排列都可能类群特异 | `case.json` → `taxon` |
| 遗传密码表 `transl_table` | 后生动物线粒体多用 5；用错表会造出**假内部终止** | 命令参数 + 事件记录 |
| 组装状态 | 完整环 / 部分 / 多 contig；决定"缺失"是注释问题还是序列问题 | `case.json` → `inputs` |
| 注释来源与版本 | MITOS2 / MitoFinder / 手工；误差方向不同 | 事件 `tool_version` |

非典型类群：`annot_check.py --allow-atypical "<理由>"` 把**基因集数量与身份**差异从 ERROR 降为待核查，
并要求逐项确认。它不是万能开关 —— **不放松**起始密码子、重叠、tRNA/rRNA 长度与链分布检查。

## 2. CDS

1. 用类群确认的 `transl_table` 重译每条 CDS，核对链方向与跨环（跨原点）坐标。
   代码按 `/codon_start` 定位读码框，并支持 `join()` 分段位置。
2. 逐条记录：起始密码子、终止密码子（完整 / 不完整 `T`/`TA`）、内部终止数、同源蛋白覆盖度、边界证据。
3. **内部终止必须解释**，排查顺序：密码表选错 → 边界/阅读框错误 → 碱基错误（测序或组装）→ 真实生物学例外。
   声明 `/transl_except` 时降为 REVIEW，仍需独立证据。不允许只改结论文字。
4. **非典型起始密码子**（如 `CGA`）不再被无条件判错：
   - 无法定起始 → `INVALID_CDS`（ERROR），并提示 `NONCANONICAL_START_REVIEW`；
   - 用 `--tolerate-start "基因:密码子"`（可重复，逐条给出）声明类群已知例外 → 记 `NONCANONICAL_START_REVIEW`（REVIEW），
     保留证据等级与不确定性。
   - **不得**为了让检查通过而伪造 5' 端缺失、随意改 `/codon_start` 或改碱基。
5. **真实 partial CDS**：`/codon_start = 2/3` 表示 5' 端不完整，代码记 `PARTIAL_CDS_5P` 且**不检查起始密码子**；
   这是与"非典型起始密码子"**不同**的问题，不要混为一谈。
6. **不完整终止（`T`/`TA`）**：只有在末端确实是 `T` 或 `TA` 前缀时才算；与转录后多聚腺苷酸化相容，
   但没有转录证据时不得宣称已实验验证，且**`note` 属自述、须交叉核验**。
   末端是普通 sense codon（例如 `GAT`）时是 ERROR，**`note` 不能豁免**。
7. `cox1` 的 5' 起始边界特别易错（历史上多次被截短或延长），单独谨慎判断。
8. 边界或序列修改必须能回答：**这个碱基是样本 reads 支持的，还是从参考"借"来的**。
9. **"无内部终止"不证明序列来自线粒体**：numt 可以不携带 in-frame 终止密码子（移码、整块缺失同样常见），
   因此不能把"翻译干净"当作线粒体来源的证据。

## 3. tRNA

1. 统一命名并保留反密码子、位置、链方向、结构工具分数、同源性与相邻基因。
2. **裸名 `trnL` / `trnS` 不得自动归一到 `trnL1`/`trnS1`**：裸名不足以判定类型。
   代码记 `UNDETERMINED_TRNA`（REVIEW），并**放弃该组"是否缺失"的判定**，等待反密码子 / 结构 / 同源证据。
3. 类型可由反密码子解析（动物线粒体标准写法）：`tag`→`trnL1(CUN)`、`taa`→`trnL2(UUR)`、
   `gct`→`trnS1(AGN)`、`tga`→`trnS2(UCN)`。反密码子与裸名所指氨基酸不一致时**不解析**，保留待确定。
4. `/gene` 写 `trnL2(UUR)` 这类带括号形式是常见写法，代码会去掉括号后归一；但**反密码子应放在 `/anticodon`**。
   基因身份检查只读 `/gene`（缺失时读 `/product`）的**首个值**；NCBI 自然语言 `product` 名
   （如 "cytochrome c oxidase subunit 1"）会被判为非标准身份，需先补 `/gene`。
5. 部分动物线粒体 tRNA 天然结构不完整，**不能仅凭三叶草结构不完整判为假基因**：
   - `trnS1(AGN)` 常见 DHU 臂不能形成稳定茎环（*Hyphantria cunea*、*Cnaphalocrocis medinalis* 等鳞翅目）；
   - 蜘蛛（*Tetragnatha*）中 `trnS1` 与 `trnS2` 均缺 DHU 臂，且多数 tRNA 丢失 TΨC 臂；
   - metazoa 线粒体 tRNA **数量本身可变**（缺失/重复，缺失可由核编码 tRNA 补偿），因此"22"是预期不是定律。
6. 代码的长度规则：`< 50 bp` 且无 `note` → ERROR（坐标级硬性可疑）；有 `note` → WARN；`60–75 bp` 之外 → WARN。
   **长度只触发检查，不用于判定身份或假基因。**

## 4. rRNA

1. 名称映射：**12S → `rrnS`；16S → `rrnL`**（代码归一并支持 `12S rRNA`/`16S rRNA`/`srrna`/`lrrna` 等写法）。
2. 用同源保守区、相邻基因与类群背景估计边界；长度差异只触发检查，**不把近缘参考的边界当真值**。
3. 代码按归一后身份选区间（`rrnL` 1100–1500、`rrnS` 600–850），超出只出 WARN。区间属 (c) 层。

## 5. 基因重叠

1. 在**实测分段坐标**上报告每一对重叠：基因对、重叠 bp、链向、涉及边界。
   跨原点 `join()` 基因按各段分别求交，**不得**用跨全长的 min/max 区间代替（否则产生假重叠、假长度）。
2. **所有重叠都记录**：`≤8 bp` 记 INFO，`>8 bp` 记 REVIEW（默认），`--overlap-severity error` 可升级为 ERROR。
3. `--tolerate-overlap "基因1,基因2"` 的含义是**"已人工审核并保留该注释"**，
   代码记 `OVERLAP_ACCEPTED` 并明确说明**不代表已证明该重叠具有功能真实性**。
   名称按归一后身份比较（大小写不敏感）；**未匹配到任何重叠对的条目会报警**，避免静默失效。
4. `atp8-atp6`、`nad4-nad4l` 的短重叠是典型真实特征：**检查它，不强制消除它**。
   也存在 `atp8-atp6` 之间是**间隔区**而非重叠的情况，两者都不是诊断判据。
5. **任何重叠都不构成"自动裁剪序列"的理由。**

## 6. 链分布

1. 正链 CDS `9+/4−` 只是**昆虫背景预警**（脊椎动物为 `12+/1−`），不跨类群推广。
2. **必须先统一整条序列的方向**：正确的昆虫基因组整体反向互补后正负链数量互换，但结构没变。
   代码在 `minus == 9` 时记 `ORIENTATION`（INFO，说明差异属表示差异），**不报 9+/4− 异常**。

## 7. 评分不跨越证据边界

"注释内部一致、基因身份可信、边界正确" 与 "原始 reads 支持该序列" 是**两个独立结论**，必须分别报告。
所有自动注释结果（**包括近缘参考自身的注释**）都可质疑。`annot_check.py` 结尾固定打印该免责声明。

## 8. 工具对照：`scripts/annot_check.py`

| 检查项 | 触发条件 | 默认等级 |
|---|---|---|
| 环状拓扑 | `--require-circular` 且 topology ≠ circular | ERROR |
| 基因数量 | CDS≠13 / tRNA≠22 / rRNA≠2 | ERROR（`--allow-atypical "<理由>"` → REVIEW） |
| 基因身份 | 缺失 / 非标准名 / 重复（`nad*`↔`nd*`、`cob`↔`cytb`、`coi/ii/iii`↔`cox1/2/3`、`12S/16S`↔`rrnS/rrnL` 已归一） | ERROR（`--allow-atypical` → REVIEW） |
| tRNA 类型 | 裸名 `trnL`/`trnS` 或反密码子与裸名不一致 | REVIEW（`UNDETERMINED_TRNA`） |
| 起始密码子 | 不在密码表合法起始集 | ERROR（`--tolerate-start "基因:密码子"` → REVIEW） |
| 5' partial | `/codon_start = 2/3` | INFO（不检查起始密码子） |
| 内部终止 | 翻译含内部 `*` | ERROR（有 `/transl_except` → REVIEW） |
| 不完整终止 | 末端为 `T`/`TA` 前缀 | REVIEW（`note` 记为自述，需交叉核验） |
| 非 `T`/`TA` 末端 | 末端是普通 sense codon | ERROR（`note` 不能豁免） |
| tRNA 长度 | `<50`（有/无 `note`）/ 60–75 之外 | ERROR / WARN / WARN |
| 反密码子 | 缺 `anticodon` qualifier | REVIEW |
| rRNA 长度 | 超出类群区间 | REVIEW（WARN） |
| 重叠 | `≤8bp` / `>8bp` | INFO / REVIEW（`--overlap-severity error` → ERROR） |
| 链分布 | 正链 ≠ 9 且反向互补后也 ≠ 9 | REVIEW（整体反向互补只记 INFO） |
| 基因顺序（需 `--ref`） | 与参考 CDS 顺序不同（按归一身份比较；**不给 `--ref` 时不做顺序对照**） | REVIEW（顺序差异不自动失败） |
| CDS 长度（需 `--ref`） | 与参考差 > 20%（按分段长度） | REVIEW |

退出码：`0` 无任何发现；`1` 有错误（含参数错误）；`2` 仅有待核查/警告。
**退出码 2 不等于通过质量门**：每条 REVIEW 必须逐条入账（被降级项 + 支持证据 + 判定人）。

## 9. 允许的类群例外（按类群文献确认，不可外推）

| 类群 | 已报告的特殊情况 | 对判据的意义 |
|---|---|---|
| 鳞翅目（多个物种） | `cox1` 常见非典型起始密码子 `CGA`（*Hyphantria cunea*、*Cnaphalocrocis medinalis*、*Maruca vitrata*、*Leucoptera malifoliella*） | 起始密码子检查必须有类群+基因背景；用 `--tolerate-start cox1:CGA` |
| 鳞翅目 | `trnS1(AGN)` DHU 臂不能形成稳定茎环 | 不得据此判假基因 |
| 蜘蛛（*Tetragnatha*） | `trnS1`/`trnS2` 缺 DHU 臂，多数 tRNA 缺 TΨC 臂 | 结构不完整是常见真实现象 |
| 绦虫（*Hymenolepis diminuta*）、软体动物（*Mytilus edulis*） | 缺失 `atp8`（12 CDS） | 13 CDS 不是所有动物的硬要求 |
| 海绵（*Aphrocallistes vastus*）等 | 线粒体 tRNA 数量减少 | 22 tRNA 不是所有动物的硬要求 |
| 六足总纲之外 | 链分布、长度区间、排列均可能不同 | (c) 层阈值不得跨类群外推 |

## 10. 交叉引用与来源

- 证据门槛、`confidence` 判据、来源索引：`evidence-standard.md` §1–§6
- 顺序、旋转、方向与坐标约定：`standard_gene_order.md`
- 工具输入/产出/退出码：`tool-catalog.md`；环境：`tool_check.md`

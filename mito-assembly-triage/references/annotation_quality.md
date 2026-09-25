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
   代码按 `/codon_start` 偏移读码框并支持 `join()` 分段位置，但**不**从 `/codon_start` 推断 partial
   （partial 只看 location，见下条）。
2. 逐条记录：起始密码子、终止密码子（完整 / 不完整 `T`/`TA`）、内部终止数、同源蛋白覆盖度、边界证据。
3. **内部终止必须解释**，排查顺序：密码表选错 → 边界/阅读框错误 → 碱基错误（测序或组装）→ 真实生物学例外。
   - `/transl_except` **不是存在即豁免**：代码按 CDS 转录顺序建立每个密码子的**精确基因组位置集合**；
     只有声明的 `pos` 位置集合与某个真实内部终止密码子的位置集合**完全相同**时才记 `TRANSL_EXCEPT_MATCHED` 并扣除该项。
     同时必须满足三条声明语义：
     （1）`aa` 为合法三字母代码或 `OTHER`，且 **`aa:TERM` 不得用来解释内部终止**（TERM 只表示在此终止，
     不能证明 CDS 可继续翻译）；（2）`pos` 的链方向与 CDS 一致 —— **负链 CDS 必须写 `pos:complement(a..b)`**，
     正链不得写 `complement`；（3）整个 qualifier 必须被完整解析。
     因此**一条声明不能吸收两个内部终止**，错帧的 3 nt 窗口（如终止在 `10..12` 却写 `11..13`）、非法的 `aa:Foo`、
     方向不符的裸 `pos` 以及带尾随损坏文本的 qualifier 都**不**算解释。
     `pos:join(a..b,c)`（密码子跨 `join()` 边界）与负链标准写法都支持：位置集合相等即匹配。
     `TRANSL_EXCEPT_UNEXPLAINED`（位置对不上）与 `TRANSL_EXCEPT_UNPARSED`（语法无法完整解析，
     **整条 qualifier 作废而非部分沿用**）都会报告；未解释的内部终止仍是 ERROR。
   - `/transl_table` 逐条与 `--table` 比对，不一致记 `TABLE_CONFLICT`（代码按 `--table` 翻译）。
   - 不允许只改结论文字。
4. **非典型起始密码子**（如鳞翅目 `cox1` 的 `CGA`）不再被无条件判错：
   - 没有合法起始 → `INVALID_CDS`（ERROR），并提示 `NONCANONICAL_START_REVIEW`；
   - 用 `--tolerate-start "基因:密码子"`（可重复、逐条）声明类群已知例外 → `NONCANONICAL_START_REVIEW`（REVIEW）。
     **例外必须引用已审计记录**：`--exception-registry <json>`，每条含 gene/codon/taxon/source/rationale。
     未引用已审计记录时额外输出 `EXCEPTION_NOT_REGISTERED`，**不得**把它当成已验证结论。
   - **不得**为了让检查通过而伪造 5' 端缺失、随意改 `/codon_start` 或改碱基。
5. **partial 来自 GenBank location，不来自 `/codon_start`**：
   - 代码直接读 Biopython 的 `BeforePosition`/`AfterPosition` **位置对象**（不解析 location 字符串），
     按链方向解释生物学的 5'/3' 端；**负链的 5' 端在高坐标**，跨原点 `join()` 的各段按转录顺序排列；
     位置对象**始终优先于**字符串兜底（仅当 parts 完全不带 fuzzy 信息时才用字符串），
     两者矛盾时以对象为准**并报** `PARTIAL_SOURCE_CONFLICT`，不静默取值；
   - `<1..N` → `PARTIAL_CDS_5P`（不检查起始密码子）；`N..>M` → `PARTIAL_CDS_3P`（**不要求**终止密码子）；
     两端同时 partial 时两者都记；
   - `location` 标注为**完整**却设 `/codon_start=2` → `CODON_START_CONFLICT`（注释自相矛盾），需先确认 5' 端是否真的缺失；
   - 这是与“非典型起始密码子”**不同**的问题，不要混为一谈。
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
   没有 `/gene` 时会读 `/product`，并识别常见自然语言写法：
   `NADH dehydrogenase subunit 5`→`nd5`、`cytochrome c oxidase subunit 1/III`→`cox1/cox3`、
   `cytochrome b`→`cytb`、`ATP synthase F0 subunit 6`→`atp6`、`16S/12S ribosomal RNA`→`rrnL/rrnS`、
   `tRNA-Leu`→**裸名** `trnL`（仍待反密码子/结构证据）。无法识别的 `ribosomal RNA` → `rrna`（未确定）。
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
| 环状拓扑 | `--require-circular` 且 topology ≠ circular | ERROR（并打印 `CIRCULAR_DECLARATION_CHECK`：只校验声明，不等于物理闭环） |
| 基因数量 | CDS≠13 / tRNA≠22 / rRNA≠2 | ERROR（`--allow-atypical "<理由>"` → REVIEW） |
| 基因身份 | 缺失 / 非标准名 / 重复（`nad*`↔`nd*`、`cob`↔`cytb`、`coi/ii/iii`↔`cox1/2/3`、`12S/16S`↔`rrnS/rrnL` 已归一） | ERROR（`--allow-atypical` → REVIEW） |
| tRNA 类型 | 裸名 `trnL`/`trnS` 或反密码子与裸名不一致 | REVIEW（`UNDETERMINED_TRNA`） |
| 起始密码子 | 不在密码表合法起始集 | ERROR（`--tolerate-start "基因:密码子"` + `--exception-registry` → REVIEW） |
| 5' partial | location 带 `<` | INFO（不检查起始密码子） |
| 3' partial | location 带 `>` | INFO（不要求终止密码子） |
| codon_start 矛盾 | location 完整却设 `/codon_start=2/3` | REVIEW（`CODON_START_CONFLICT`） |
| 遗传密码表 | `/transl_table` 与 `--table` 不一致 | REVIEW（`TABLE_CONFLICT`） |
| 内部终止 | 翻译含内部 `*` 且未被 `transl_except` 对应 | ERROR（对应时记 `TRANSL_EXCEPT_MATCHED`） |
| `/transl_except` | 声明位置未对应任何内部终止 | REVIEW（`TRANSL_EXCEPT_UNEXPLAINED`） |
| 不完整终止 | 末端为 `T`/`TA` 前缀 | REVIEW（`note` 记为自述，需交叉核验） |
| 非 `T`/`TA` 末端 | 末端是普通 sense codon | ERROR（`note` 不能豁免） |
| tRNA 长度 | `<50`（有/无 `note`）/ 60–75 之外 | ERROR / WARN / WARN |
| 反密码子 | 缺 `anticodon` qualifier | REVIEW |
| rRNA 长度 | 超出类群区间 | REVIEW（WARN） |
| rRNA 身份未确定 | 名称无法识别（如 `rrna`） | REVIEW（`UNDETERMINED_RRNA`，不套用区间） |
| 重叠 | `≤8bp` / `>8bp` | INFO / REVIEW（`--overlap-severity error` → ERROR） |
| 重叠豁免无法定位 | 同名基因重复时只用名字豁免 | REVIEW（`TOLERATE_OVERLAP_AMBIGUOUS`） |
| 链分布 | 正链 ≠ 9 且反向互补后也 ≠ 9 | REVIEW（整体反向互补只记 INFO） |
| 基因集 | 与 `--ref` 相比有缺失/多余 | REVIEW（`GENE_SET_DIFF`，与顺序差异分开） |
| 基因排列（需 `--ref`） | 共有基因的**环状邻接**与参考不同（已归一旋转与整链反向互补） | REVIEW（`ARRANGEMENT_DIFF`；不给 `--ref` 时不做顺序对照） |
| CDS 长度（需 `--ref`） | 与参考差 > 20%（按分段长度） | REVIEW |

退出码：`0` 无任何发现；`1` 有错误（含参数错误）；`2` 仅有待核查/警告。
**退出码 2 不等于通过质量门**：每条 REVIEW 必须逐条入账（被降级项 + 支持证据 + 判定人）。

## 9. 允许的类群例外（按类群文献确认，不可外推）

例外不是"加一个开关"：代码侧用 `--tolerate-start "基因:密码子"` 选择，
**生物学理由放在已审计记录** `--exception-registry <json>` 中（每条含 `gene`/`codon`/`taxon`/`source`/`rationale`）。
登记检查是逐项匹配：

- 未引用任何 registry → `EXCEPTION_NOT_REGISTERED`（仍是 REVIEW，可用但不得当作已验证结论）；
- 记录缺 `taxon`/`source`/`rationale` → `EXCEPTION_RECORD_INCOMPLETE`；
- 提供了 `--taxon` 且与记录的 `taxon` 不一致 → `EXCEPTION_TAXON_MISMATCH`。

以上三种情形**都只维持 REVIEW**，不会把例外升级为已验证。

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

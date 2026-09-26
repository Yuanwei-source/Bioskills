# 公共参考获取政策（reference policy）

> 本文件统一定义参考等级、用途、获取与登记。参考作为分析输入固定版本与 hash，
> 避免数据库更新后无法复现。登记与获取由
> `scripts/reference_registry.py` 执行，字段定义在 `schemas/reference-registry.schema.json`。

## 1. 五条原则

1. **公共参考序列可以下载**用于分析——这不是危险操作，不需要额外开关。
2. **用户样本序列不得上传**公共数据库或第三方服务，除非获得明确授权
   （`cox1_id.py` 的 `--allow-public-upload`）。**下载与上传是两个独立事件。**
3. **每条参考必须登记**：`source` + `accession`(**含版本**) + `hash` + `purposes`。
4. **版本必须固定**：只接受 `ACCESSION.VERSION`；同一 accession 内容变化即报错（版本漂移）。
5. **参考相似性不是独立证据**：参考注释自身也可质疑，命中只是**定位线索**，
   不能单独证明身份、缺失或重排（见 `evidence-standard.md` §1/§4）。

## 2. 两个事件必须分开（AI 最容易混淆的地方）

| | 公共参考获取（acquisition） | 样本数据上传（upload） |
|---|---|---|
| 数据流向 | 公共 → 本地 | 本地 → 外部 |
| 是否需要额外授权 | **不需要**（默认允许） | **需要**本次明确授权 |
| 必须记录 | source / accession.version / hash / purposes / 日期 | 发送内容、接收方、风险、用户同意 |
| 典型命令 | `reference_registry.py acquire --accession NC_060773.1 ...` | `cox1_id.py ... --allow-public-upload` |
| 失败后果 | 结论不可复现 | 未发表数据泄露（不可逆） |

**不要**设计 `--allow-reference-download` 这类开关：参考下载本身不危险，
把授权负担放在下载上只会让人绕过它；真正需要把关的是**上传**。
已有明确授权在其内容、接收方与用途范围内持续有效，无需重复确认；用户限制联网时遵守该限制。

## 3. 在诊断循环中的位置

| 阶段 | 做什么 | 明确不做 |
|---|---|---|
| **INTAKE** | 只记录**需求**：`reference_needed`（true/false）与 `purpose` 列表（如 `gene_order_comparison`）；记录手头已有参考 | **不下载**任何东西 |
| **CHOOSE_TEST** | 按具体命题选择等级与用途；同科参考可提供基因身份/注释比较线索，不能独自判定 split 是否真实；闭环验证不依赖参考相似性 | 不下载 |
| **EXECUTE** | 才执行获取/登记（`acquire`）或用 `register` 登记已有文件；把命令、accession、hash 写进 `events.jsonl` | 不把参考当样本真值 |
| **DECIDE / 报告** | 用 `case-reference` 把登记项关联进案例，报告引用 `accession.version` + hash | 不写"与近缘物种比较"这类不可复现的说法 |

## 4. 登记表（registry）

- 默认位置：`$MITO_KNOWLEDGE_DIR/references/registry.json`。仅做当前任务时，将
  `MITO_KNOWLEDGE_DIR` 指向任务内目录，或对登记和 `case-reference` 同时给 `--registry`；
  不因下载参考而隐式启用跨任务经验积累，不写进 skill 安装目录。
- 格式：`mito-reference-registry-1`，字段定义见 `schemas/reference-registry.schema.json`。
- 记录字段（每条）：

| 字段 | 说明 |
|---|---|
| `id` | `ref-001`（自动分配） |
| `kind` | `sequence`（参考序列）或 `database`（本地库/筛选库） |
| `accession` | **必须含版本**：`NC_060773.1`；裸 accession 被拒绝 |
| `organism` / `taxon` | 物种与类群（判断等级与用途的基础） |
| `source` | 来源（如 `NCBI Nucleotide`、本地 RefSeq 镜像） |
| `retrieved` | 获取日期 |
| `version` | 序列记录用库版本（如 `GenBank 2023-04`）；`database` 必填库版本/时间戳 |
| `level` | `L1` 同种 … `L5` 远缘（**限制用途**，见 §5） |
| `purposes` | 允许用途枚举（见 §5） |
| `file_sha256` / `sequence_sha256` / `length_bp` | 完整性固定 |

```bash
# 在用户任务工作目录运行；同一任务的命令使用同一环境设置
export MITO_KNOWLEDGE_DIR="$PWD/work/knowledge"
# 1) 先看要下载什么（不联网、不写盘）
python3 <skill目录>/scripts/reference_registry.py acquire --accession NC_060773.1 \
  --purposes gene_order_comparison --level L3 --dry-run

# 2) 获取并登记（公共参考下载无需额外授权）
python3 <skill目录>/scripts/reference_registry.py acquire --accession NC_060773.1 \
  --purposes gene_order_comparison --level L3 --source 'NCBI Nucleotide'

# 3) 已手工获取的文件也要登记
python3 <skill目录>/scripts/reference_registry.py register --file refs/NC_060773.1.gb \
  --accession NC_060773.1 --source 'NCBI Nucleotide' --level L3 \
  --purposes gene_order_comparison

# 4) 随时校验是否缺失或被改动
python3 <skill目录>/scripts/reference_registry.py verify
```

`verify` 会在文件缺失或 hash 变化时报警；**同一 accession 但 hash 不同**的登记会被拒绝
（`--force` 才能覆盖，且应在案例里说明原因）。

## 5. 等级限制用途（不是"标个等级"）

这是唯一文档版用途矩阵，与 `scripts/reference_registry.py` 的 `LEVEL_PURPOSES` 对应。
它是当前工具的保守使用约束，不是按分类阶元定义的普遍生物学定律。等级相同的参考仍需
审查注释质量和实际可比性；通过登记不等于证据充分。

| 等级 | 参考 | 允许用途 | **禁止**用途 |
|---|---|---|---|
| `L1` | 同种 | 全部（边界、变异、结构、顺序、身份） | — |
| `L2` | 同属 | 全部（边界需谨慎并留证） | — |
| `L3` | 同科 | `gene_order_comparison`、`annotation_comparison`、`identity_check`、`taxon_screen`、`contamination_screen` | `boundary_validation`、`structure_check` |
| `L4` | 同目 | `identity_check`、`taxon_screen`、`contamination_screen` | 顺序、边界、结构判读 |
| `L5` | 远缘 | `taxon_screen`（**仅定位**） | **tRNA 丢失、基因重排、注释比较**等一切结论 |

用途枚举：`gene_order_comparison`、`annotation_comparison`、`boundary_validation`、
`identity_check`、`taxon_screen`、`contamination_screen`、`structure_check`。
**在登记时就要写清用途**；事后扩大使用范围会被 `case-reference` 拒绝（该记录未声明该用途）。

所有等级（包括 L1/L2）都不能把参考当作样本碱基或结构真值。`identity_check` 是身份检查
用途，不表示单次命中已确认身份；`taxon_screen` 只给筛查线索。探索中发现的差异可以记录，
但不能以“只是探索”为由绕开矩阵，将远缘线索升级为边界、丢失或重排结论。

## 6. 版本、漂移与"半年后 top hit 变了"

- 只用 `ACCESSION.VERSION`；写报告时也必须带版本号。
- `verify` 检测内容漂移；登记时的 hash 冲突是**错误**而不是自动更新。
- **本地库（BLAST/refseq）同样要登记版本**——它不是一个 accession：

```bash
python3 <skill目录>/scripts/reference_registry.py record-database \
  --name 'MITOS2 refseq invertebrate' --path <本地数据库目录> \
  --version refseq89m --taxon-filter Metazoa --purposes taxon_screen
```

  否则"今天 top hit 是 X"半年后无法复核（库更新会改变 top hit，而记录里只有一句结论）。
- 批量比对用的核酸库同理：记录**库路径 + 版本/时间戳 + 类群过滤 + 目录指纹**。

## 7. 案例如何引用参考（证据链闭合）

```bash
python3 <skill目录>/experimental/experience/experience.py case-reference work/case-001 \
  --reference-id ref-001 --purpose gene_order_comparison
```

- 未登记的 id、或记录未声明的用途 → **直接失败**（不允许悬空引用）；
- 写入 `case.json` 的 `references[]`：`reference_id` + `accession` + `level` + `purpose` + `file_sha256`；
- `case-report` 会渲染"参考（已登记：版本 + hash 固定）"一节，报告里因此可以写
  "与 `NC_060773.1`（登记日期、hash `cf87c277…`、L3）比较"，而不是"与近缘物种比较"。

## 8. 常见错误

| 错误 | 后果 |
|---|---|
| 用裸 accession（`NC_060773`） | 参考更新后无法复现，也无法判断是否同一版本 |
| 拿 L5 远缘参考判 tRNA 缺失 / 重排 | 把参考距离造成的差异当成生物学结论 |
| 参考命中被当作独立证据计数 | 与同源证据重复计数（`evidence-standard.md` §4） |
| 只记"下载了参考"，不记 hash | 无法发现参考被替换或截断 |
| 把 BLAST 库当成"不需要版本的资源" | top hit 随时间变化，结论不可复核 |

## 9. 与其他文件的关系

- 诊断分流与参考选择入口：[diagnostic-decision-tree.md](diagnostic-decision-tree.md)
- 证据门槛与独立性：`evidence-standard.md` §1/§3/§4
- 案例记录与报告：`conclusion-report.md`；命令接口：`tool-catalog.md`
- 字段定义：`schemas/reference-registry.schema.json`；测试映射：`developer-contract.md`

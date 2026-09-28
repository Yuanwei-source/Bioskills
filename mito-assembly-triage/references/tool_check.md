# 工具检查与环境配置

> **原则**：本 skill 不假定任何机器路径。所有机器相关的路径只存在于 `config/env.sh` 或环境变量，
> **绝不写回 skill 公共文件**；路径与文件 hash 也可能关联未发表项目（见 `learning-policy.md`）。

## 1. 环境变量（`config/env.sh` 为唯一落点）

| 变量 | 含义 | 谁需要 |
|---|---|---|
| `PROFILE`（默认由 `MITO_PROFILE` 提供，缺省 `generic`） | 机器档案名，仅用于日志辨识 | 全部 |
| `CONDA_ROOT` | 包含环境目录的根路径（`$CONDA_ROOT/envs/*/bin` 会被搜索） | 全部 |
| `MINIMAP2` | minimap2 可执行文件（不在 PATH 时必填） | 参考比对、环化候选 |
| `MITOS2_PY` | 可 `import mitos` 的 Python | 仅 MITOS2 注释 |
| `MITOS2_REFDIR` | 含参考数据库的目录 | 仅 MITOS2 注释 |
| `MITOS2_REFSEQVER` | 数据库版本名（默认 `refseq89m`） | 仅 MITOS2 注释 |
| `MITOS2_EXTRA_PATH` | MITOS2 依赖命令的额外 PATH（`cmsearch`/`RNAplot` 等） | 仅 MITOS2 注释与绘图 |
| `PLOT_PY` | 含 BioPython + matplotlib 的 Python | 仅环形图 |

## 2. 依赖分层与两种检查模式

**唯一来源**：`config/dependencies.json`。文档、`tools/env_check.py`、各脚本的门禁都从它出发；
**不要**在别处再抄一份依赖清单——重复的清单就是下一个会漂移的缺陷（同 #14 删掉的双校验实现）。

### 2.1 两种模式（部署一次，多次使用）

| 模式 | 何时用 | 行为 |
|---|---|---|
| `python3 tools/env_check.py --setup`（或 `bash scripts/check_env.sh --setup`） | 第一次、换机、升级后 | 全量盘点，按 essential/extended/optional 分组；缺失项给出**用途 + 安装方式**；essential 齐全时写入 environment lock |
| `python3 tools/env_check.py --daily` | 日常每次使用 | 读 lock 做**轻量**检查：只探测本次需要的那一组是否还在。不再让用户面对三十个工具 |
| `python3 tools/env_check.py --stage <NAME>` | 脚本前置门禁 | 只检查某个 capability 的依赖；缺则**退出码 3**，并注明"该步骤未执行" |

- **lock 位置**：`$MITO_KNOWLEDGE_DIR/environment.lock.json`（环境属于机器，不属于某个任务；**不写进 skill 安装目录**），可用 `--lock` 覆盖。
- **不变式**：lock 是**缓存，不是信任凭证**。日常模式仍真实探测所需工具，lock 只用来避免重复全量盘点并提供"上次验证版本"用于漂移对比。只信 lock 会退化成"第一次通过、以后永远相信"——正是本 skill 反复清除的静默降级。
- 退出码：`0` 就绪 / `2` essential 缺失（与 `check_env.sh` 一致）/ `3` 指定 stage 缺依赖 / `1` 用法或清单错误。
- `--strict` 要求 extended 也齐全（装机/CI 用）；`--path DIR` 只在指定目录内查找可执行文件（检查某个 conda env 的 bin，**不回退 PATH**）；`--network` 才探测联网类依赖（默认不探测，避免慢）。

### 2.2 依赖分层

**essential**（核心诊断能力；缺了应当补齐，或明确报告哪些能力不可用）

| 依赖 | 用途 | 安装 |
|---|---|---|
| `python3` | 运行本 skill 的全部脚本 | 系统包管理器或 conda 安装 python>=3.8 |
| `jsonschema` | 案例校验（唯一实现路径）；缺失时 `case-validate` 等以退出码 3 失败 | `python3 -m pip install jsonschema` |
| `biopython` | GenBank/FASTA 解析、密码表、序列统计。**边界已确定：对 `annot_check.py`、`blast_genes.py`、`circularize.py`、`mitos2_to_genbank.py`、`seq_stats.py` 是必需**；`depth_analysis.py` 不使用 | `python3 -m pip install biopython` |
| `blastn` | 基因定位与同源检索 | `conda install -c bioconda blast` 或 `apt-get install ncbi-blast+` |
| `makeblastdb` | 建自比对库 | 随 `blastn` 一同安装 |
| `samtools` | BAM 排序/索引/深度剖面 | `conda install -c bioconda samtools` |
| `minimap2` | reads 回贴（深度与接缝证据） | `conda install -c bioconda minimap2` |

**extended**（增强能力；缺失不阻断基础诊断）

| 依赖 | 用途 | 安装 |
|---|---|---|
| `seqkit` | FASTQ/FASTA 快速体检与子集抽取 | `conda install -c bioconda seqkit` |
| `bwa` | 备选比对器（minimap2 不适用时） | `conda install -c bioconda bwa` |
| `mitos2` | 线粒体注释（独立第二意见） | 见 https://mitos2.bioinf.uni-leipzig.de/ |
| `mitos2_refdir` | MITOS2 参考数据库（数十 GB，通常镜像到本地/NAS） | 设置 `MITOS2_REFDIR` 与 `MITOS2_REFSEQVER` |
| `cmsearch` | MITOS2 的 rRNA 搜索（Infernal） | `conda install -c bioconda infernal`（加入 `MITOS2_EXTRA_PATH`） |
| `tRNAscan-SE` | tRNA 结构预测（tRNA 缺失分级的独立证据） | `conda install -c bioconda trnascan-se` |
| `mitofinder` | 第二条独立注释路径（需近缘参考 GenBank） | `conda install -c bioconda mitofinder` |
| `plot_python` | 环形图绘制（`run_circular_map.sh`） | 在 `PLOT_PY` 环境里 `pip install biopython matplotlib` |

**optional**（按需）

| 依赖 | 用途 | 安装 |
|---|---|---|
| `network` | 公共参考获取 / COX1 远程查询 / 公共同步 | 无需安装；受限网络下这些步骤不可用，其余步骤不受影响 |
| `table2asn` | NCBI 提交预检（不在本 skill 自检范围内） | 见 NCBI 官方页面（人工执行） |

### 2.3 脚本前置门禁（每步自己把门）

每个脚本在干活之前调用 `scripts/_deps.py` 的 `require_stage('<stage>')`：它从**同一份清单**取该阶段
所需依赖并探测，缺失时打印统一说明（用途 + 安装方式 + "该步骤未执行"）并以**退出码 3** 结束。
脚本 → stage：`seq_stats.py`→`seq_stats`、`annot_check.py`→`annot_check`、`blast_genes.py`→`gene_locating`、
`depth_analysis.py`→`read_evidence`、`circularize.py`→`circularize`、`mitos2_to_genbank.py`→`mitos2_bridge`、
`run_mitos2.sh`→`annot_independent`、`run_circular_map.sh`→`circular_plot`。

- **只阻断该步骤**：缺 samtools 不拦住只读 FASTA 的体检；缺 MITOS2 只在跑注释时失败（不再跑到一半才炸）。
- `cox1_id.py` 与 `reference_registry.py` **不按 stage 门禁**：网络依赖只在具体动作里需要（`acquire`/远程查询），
  且失败已用退出码 3 表达；给整脚本加门会连离线可用的 `list`/`verify` 一起挡掉。
- **不为例外开口子**：`--help` 也走门禁（脚本模块级就可能 import Biopython，给 help 放行只会制造
  "看着能跑、其实缺依赖"的错觉）。
- `MITO_ENV_GATE=off`：**仅用于测试/自检**（例如用 stub 解释器验证包装脚本自身的 arg/mkdir 行为），
  会打印告警；不得用于日常诊断——它关掉的正是"缺依赖即失败"这条保护。

阶段与依赖的对应（`stages`）就在同一份清单里：`seq_stats` / `annot_check` / `gene_locating` /
`read_evidence` / `circularize` / `annot_independent` / `circular_plot` /
`reference_acquisition` / `cox1_identity`。**只有该阶段需要的依赖缺失时，才阻断该阶段**——
缺 bwa 不阻断仅 FASTA/GB 检查；缺 samtools 才阻断依赖它的 BAM 检查。

**设计原则**：关键证据链可以失败，但不能静默降级——同一个结论不允许有第二条实现路径。


## 3. 自检

```bash
bash scripts/check_env.sh
```

检查内容：

1. **核心工具**：`blastn`、`makeblastdb`、`bwa`、`samtools`、`seqkit`、`python3`（PATH → `$CONDA_ROOT/envs/*/bin` 顺序查找）；
2. `minimap2`：PATH → `$MINIMAP2` → envs；
3. **组装/注释工具**：`get_organelle_from_reads.py`、`mitofinder`、`cmsearch`、`RNAplot`、`tRNAscan-SE`、`arwen`；
4. `$MITOS2_PY` 能否 `import mitos`；
5. `$MITOS2_REFDIR/$MITOS2_REFSEQVER` 是否存在；
6. `MITOS2_EXTRA_PATH` 是否为空。

退出码：`0` 通过；`2` 核心依赖缺失（不应进入依赖这些工具的流程）。

## 4. 已知问题与处置

| 现象 | 处置 |
|---|---|
| MITOS2 报 `no such directory <path>` | 先看路径指到哪里：指向 **`--outdir`** 时说明输出目录不存在——`run_mitos2.sh` 现在会自己 `mkdir -p`，绕过它直接调 MITOS2 时需自己建；只有路径指向参考库时才设置 `MITOS2_REFDIR`（目录**尾斜杠**）与 `MITOS2_REFSEQVER` |
| MITOS2 报 `cmsearch` / `plotprot.R` / `RNAplot` / `drawmitos` 找不到 | 在 `MITOS2_EXTRA_PATH` 补齐依赖 PATH；重装后需修复 `drawmitos` wrapper |
| MitoFinder 报 `install.sh.ok` / `Mitofinder.config` 缺失 | 核查安装步骤、版本与实际依赖；按该版本安装说明恢复配置，不靠创建成功标记伪装安装完成 |
| `depth_analysis.py` 报 BAM 相关错误 | 确认 BAM 已 `samtools sort` + `samtools index`，且参考名与 BAM 头一致 |
| `circularize.py` 报"最高分候选不唯一" | 属**预期保护**：不能凭首个候选猜方向；补充 scaffold / 独立证据或标 `UNRESOLVED` |
| 环形图报"需要 BioPython + matplotlib" | 改 `PLOT_PY` 指向同时含这两个库的 Python |
| GetOrganelle 长时间无输出 | 超过 5 分钟一律后台化：`run_bg.sh` + `check_bg.sh`，不要在前端阻塞 |

## 5. 迁移与复现

- 换机时以 `config/env.sh` 为模板校正新机路径，再按需要跑 `check_env.sh`；
- MITOS2 数据库目录可整体同步（如 `rsync` 旧机该目录）；
- 报告里写清**哪个步骤因缺依赖未执行**，这比给出未验证结论更重要。

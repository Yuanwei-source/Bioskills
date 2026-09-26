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

## 2. 最小依赖原则（不要过度安装）

- 只做**序列体检 / 注释质检 / 基因定位**等诊断时，**不需要**先安装完整 MITOS2；
  需要 `BioPython`；`blastn`/`makeblastdb` 用于定位与环化候选。
- `check_env.sh` 的**核心工具缺失**才阻断通用流程；组装/注释工具缺失只影响对应步骤，
  应在报告中说明该步骤未执行，而不是伪造结论。

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
| MitoFinder 报 `install.sh.ok` / `Mitofinder.config` 缺失 | 手动创建/复制对应文件，路径加**尾斜杠** |
| `depth_analysis.py` 报 BAM 相关错误 | 确认 BAM 已 `samtools sort` + `samtools index`，且参考名与 BAM 头一致 |
| `circularize.py` 报"最高分候选不唯一" | 属**预期保护**：不能凭首个候选猜方向；补充 scaffold / 独立证据或标 `UNRESOLVED` |
| 环形图报"需要 BioPython + matplotlib" | 改 `PLOT_PY` 指向同时含这两个库的 Python |
| GetOrganelle 长时间无输出 | 超过 5 分钟一律后台化：`run_bg.sh` + `check_bg.sh`，不要在前端阻塞 |

## 5. 迁移与复现

- 换机时先复制 `config/env.sh`，再跑 `check_env.sh`；
- MITOS2 数据库目录可整体同步（如 `rsync` 旧机该目录）；
- 报告里写清**哪个步骤因缺依赖未执行**，这比给出未验证结论更重要。

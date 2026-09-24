# 工具检查与环境配置

本 skill 不假定任何机器路径。运行前通过环境变量或 `config/env.sh` 设置：

- `CONDA_ROOT`：包含环境目录的根路径
- `MITOS2_PY`：可导入 `mitos` 的 Python
- `MITOS2_REFDIR`：包含参考数据库的目录
- `MITOS2_REFSEQVER`：数据库版本名
- `MITOS2_EXTRA_PATH`：MITOS2 依赖命令的额外 PATH
- `PLOT_PY`：包含 Biopython 与 matplotlib 的 Python
- `MINIMAP2`：minimap2 可执行文件

先运行 `bash scripts/check_env.sh`，再运行目标工具。路径、数据库和样本数据必须由使用者自行提供，不能写入公共 skill。

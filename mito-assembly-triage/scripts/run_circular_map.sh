#!/bin/bash
# 环形基因图统一入口：绘图 python (BioPython+matplotlib) 从 config/env.sh 的 PLOT_PY 读取
# 用法（任何机器一致）:
#   bash scripts/run_circular_map.sh <final.gb> <out.png> --title "<物种> mitochondrial genome"
set -euo pipefail
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SKILL_DIR/config/env.sh"

if [ -z "${PLOT_PY:-}" ] || [ ! -x "$PLOT_PY" ]; then
  echo "[run_circular_map.sh] 错误: PLOT_PY 不存在 ($PLOT_PY) — 请检查 config/env.sh" >&2
  exit 1
fi
echo "[run_circular_map.sh] profile=$PROFILE plot_py=$PLOT_PY"
exec "$PLOT_PY" "$SKILL_DIR/assets/circular_map.py" "$@"

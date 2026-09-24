#!/bin/bash
# MITOS2 统一入口：python / 数据库 / cmsearch PATH 全部从 config/env.sh 读取
# 用法（任何机器一致）:
#   bash scripts/run_mitos2.sh -i <genome.fasta> -o <outdir> -c 5
# 数据库 --refdir/--refseqver 自动取 config/env.sh 的 MITOS2_REFDIR/MITOS2_REFSEQVER
set -euo pipefail
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SKILL_DIR/config/env.sh"

if [ -z "${MITOS2_PY:-}" ] || [ ! -x "$MITOS2_PY" ]; then
  echo "[run_mitos2.sh] 错误: MITOS2_PY 不存在 ($MITOS2_PY) — 请检查 config/env.sh" >&2
  exit 1
fi
if [ -z "${MITOS2_EXTRA_PATH:-}" ]; then :; else
  export PATH="$MITOS2_EXTRA_PATH:$PATH"
fi

echo "[run_mitos2.sh] profile=$PROFILE python=$MITOS2_PY refdir=$MITOS2_REFDIR refseqver=$MITOS2_REFSEQVER"
exec "$MITOS2_PY" -m mitos.scripts.runmitos "$@" \
  --refdir "$MITOS2_REFDIR/" --refseqver "$MITOS2_REFSEQVER"

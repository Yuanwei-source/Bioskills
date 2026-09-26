#!/bin/bash
# MITOS2 统一入口：python / 数据库 / cmsearch PATH 全部从 config/env.sh 读取
# 用法（任何机器一致）:
#   bash scripts/run_mitos2.sh -i <genome.fasta> -o <outdir> -c 5
#   （--outdir 不存在时由本脚本创建；MITOS2 自身不会创建它）
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

# MITOS2 (runmitos.py) requires --outdir to already EXIST and aborts with
# "no such directory <outdir>" otherwise.  Create it here so that error message
# never surfaces, and so a missing parent path is a clear wrapper-level failure.
outdir=""
i=1
while [ "$i" -le "$#" ]; do
  arg="${!i}"
  case "$arg" in
    -o|--outdir)
      next=$((i + 1))
      if [ "$next" -le "$#" ]; then outdir="${!next}"; fi
      ;;
  esac
  i=$((i + 1))
done
if [ -n "$outdir" ]; then
  if [ -e "$outdir" ] && [ ! -d "$outdir" ]; then
    echo "[run_mitos2.sh] 错误: --outdir 已存在且不是目录: $outdir" >&2
    exit 1
  fi
  if ! mkdir -p "$outdir"; then
    echo "[run_mitos2.sh] 错误: 无法创建输出目录: $outdir" >&2
    exit 1
  fi
fi

printf '[run_mitos2.sh] profile=%s python=%s refdir=%s refseqver=%s outdir=%s\n' \
  "$PROFILE" "$MITOS2_PY" "$MITOS2_REFDIR" "$MITOS2_REFSEQVER" "${outdir:-(由 MITOS2 参数决定)}"
exec "$MITOS2_PY" -m mitos.scripts.runmitos "$@" \
  --refdir "$MITOS2_REFDIR/" --refseqver "$MITOS2_REFSEQVER"

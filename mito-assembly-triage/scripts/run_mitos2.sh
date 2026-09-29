#!/bin/bash
# MITOS2 统一入口：python / 数据库 / cmsearch PATH 全部从 config/env.sh 读取
# 用法（任何机器一致）:
#   bash scripts/run_mitos2.sh -i <genome.fasta> -o <outdir> -c 5
#   （--outdir 不存在时由本脚本创建；MITOS2 自身不会创建它）
# 数据库 --refdir/--refseqver 自动取 config/env.sh 的 MITOS2_REFDIR/MITOS2_REFSEQVER
set -euo pipefail
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "$SKILL_DIR/config/env.sh"

# MITOS2 calls an R plotting helper during annotation. Discover an existing
# Rscript that can load reshape2 without activating or modifying that env.
rscript_has_reshape2() {
  [ -n "$1" ] && [ -x "$1" ] && \
    "$1" -e 'if (requireNamespace("reshape2", quietly=TRUE)) quit(status=0) else quit(status=1)' \
      >/dev/null 2>&1
}
if [ -z "$MITOS2_RSCRIPT" ]; then
  candidate="$(command -v Rscript 2>/dev/null || true)"
  if rscript_has_reshape2 "$candidate"; then
    MITOS2_RSCRIPT="$candidate"
  elif [ -n "$CONDA_ROOT" ]; then
    for env_bin in "$CONDA_ROOT"/envs/*/bin; do
      candidate="$env_bin/Rscript"
      if rscript_has_reshape2 "$candidate"; then MITOS2_RSCRIPT="$candidate"; break; fi
    done
  fi
elif ! rscript_has_reshape2 "$MITOS2_RSCRIPT"; then
  echo "[run_mitos2.sh] 错误: MITOS2_RSCRIPT 无法加载 reshape2 ($MITOS2_RSCRIPT)" >&2
  exit 3
fi

# MITOS2 often lives in a dedicated conda env. If the machine profile did not
# pin its interpreter, discover an environment that can import mitos and use
# that absolute interpreter. Never activate or modify the environment.
if [ -z "${MITOS2_PY:-}" ] && [ -n "${CONDA_ROOT:-}" ]; then
  for candidate in "$CONDA_ROOT"/envs/*/bin/python "$CONDA_ROOT"/envs/*/bin/python3; do
    [ -x "$candidate" ] || continue
    if "$candidate" -B -c 'import mitos' >/dev/null 2>&1; then
      MITOS2_PY="$candidate"
      break
    fi
  done
fi

# Resolve MITOS2's helper commands across separate conda environments too.
helper_bins=()
if [ -n "$MITOS2_RSCRIPT" ]; then helper_bins+=("$(dirname "$MITOS2_RSCRIPT")"); fi
if [ -n "${CONDA_ROOT:-}" ]; then
  for helper in cmsearch RNAplot; do
    for env_bin in "$CONDA_ROOT"/envs/*/bin; do
      if [ -x "$env_bin/$helper" ]; then helper_bins+=("$env_bin"); break; fi
    done
  done
fi
if [ -n "${MITOS2_EXTRA_PATH:-}" ]; then helper_bins+=("$MITOS2_EXTRA_PATH"); fi
MITOS2_EXTRA_PATH=$(IFS=:; echo "${helper_bins[*]}")
export MITOS2_PY MITOS2_EXTRA_PATH MITOS2_RSCRIPT

if [ -z "${MITOS2_PY:-}" ] || [ ! -x "$MITOS2_PY" ]; then
  echo "[run_mitos2.sh] 错误: MITOS2_PY 不存在 ($MITOS2_PY) — 请检查 config/env.sh" >&2
  exit 1
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
  outdir="$(realpath "$outdir")"
fi

# plotstst.R and the other MITOS2 helpers use /usr/bin/env to find Rscript.
# The MITOS2 Python environment can itself contain an Rscript without
# reshape2, so simply prepending its bin directory would shadow the verified
# Rscript selected above. Put a private, single-command shim first; keep the
# MITOS2 Python bin ahead of general helper bins so python3 remains compatible.
runtime_bin=""
if [ -n "$outdir" ] && [ -n "$MITOS2_RSCRIPT" ]; then
  runtime_bin="$(mktemp -d "$outdir/.mitos2-runtime.XXXXXX")"
  ln -s "$MITOS2_RSCRIPT" "$runtime_bin/Rscript"
  trap 'rm -rf -- "$runtime_bin"' EXIT
fi
python_bin="$(dirname "$MITOS2_PY")"
if [ -n "$MITOS2_EXTRA_PATH" ]; then
  if [ -n "$runtime_bin" ]; then
    export PATH="$runtime_bin:$python_bin:$MITOS2_EXTRA_PATH:$PATH"
  else
    export PATH="$python_bin:$MITOS2_EXTRA_PATH:$PATH"
  fi
elif [ -n "$runtime_bin" ]; then
  export PATH="$runtime_bin:$python_bin:$PATH"
else
  export PATH="$python_bin:$PATH"
fi

printf '[run_mitos2.sh] profile=%s python=%s refdir=%s refseqver=%s outdir=%s\n' \
  "$PROFILE" "$MITOS2_PY" "$MITOS2_REFDIR" "$MITOS2_REFSEQVER" "${outdir:-(由 MITOS2 参数决定)}"
printf '[run_mitos2.sh] Rscript=%s\n' "${MITOS2_RSCRIPT:-(未找到可加载 reshape2 的 Rscript)}"
# 前置门禁：参数检查与 mkdir 之后、真正调用 MITOS2 之前（唯一清单来源 config/dependencies.json）。
# MITO_ENV_GATE=off 仅用于测试/自检（例如用 stub 解释器验证包装脚本自身的 arg/mkdir 行为）。
if [ "${MITO_ENV_GATE:-}" = "off" ]; then
  echo "[run_mitos2.sh] ⚠ MITO_ENV_GATE=off：已关闭环境门禁（仅用于测试/自检）" >&2
elif ! "$MITOS2_PY" "$SKILL_DIR/tools/env_check.py" --stage annot_independent; then
  echo "  → 该步骤未执行；完整环境盘点: python3 tools/env_check.py --setup" >&2
  exit 3
fi

set +e
"$MITOS2_PY" -m mitos.scripts.runmitos "$@" \
  --refdir "$MITOS2_REFDIR/" --refseqver "$MITOS2_REFSEQVER"
rc=$?
set -e
if [ -n "$runtime_bin" ]; then
  rm -rf -- "$runtime_bin"
  trap - EXIT
fi
exit "$rc"

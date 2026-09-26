#!/bin/bash
# 检查线粒体组装诊断流水线所需的环境与工具
# 机器路径统一在 config/env.sh（本脚本自动 source），换机无需改本文件
# 用法: bash check_env.sh
set -u
set -o pipefail
missing=0
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=../config/env.sh
source "$SKILL_DIR/config/env.sh"

# 环境变量未设置时的假阴性提示：工具发现只在 PATH 与 $CONDA_ROOT/envs/*/bin 中搜索，
# 因此 CONDA_ROOT 为空时会把“装了但找不到”报成“缺失”，把人引向错误的排查方向。
env_unset_hint=""
if [ -z "${CONDA_ROOT:-}" ]; then
  env_unset_hint="CONDA_ROOT"
fi

if [ -n "$env_unset_hint" ]; then
  echo "⚠ 未设置 $env_unset_hint —— 工具发现被限制在 PATH 内，下面的“未找到”可能是**假阴性**"
  echo "  config/env.sh 只读取环境变量、不写死机器路径（设计如此）；请先导出，例如："
  echo "    export CONDA_ROOT=/path/to/miniforge3"
  echo "  再重跑本脚本，不要先去安装软件。"
fi

# 把实际将被使用的可执行文件先报出来，使“检查结果”与“真正执行的路径”一致。
if [ -n "${MITOS2_PY:-}" ] && [ -x "$MITOS2_PY" ]; then
  printf '实际将使用: MITOS2_PY=%s\n' "$MITOS2_PY"
else
  printf '实际将使用: MITOS2_PY=%s（未设置或不可执行）\n' "${MITOS2_PY:-}"
fi
printf '实际将使用: MINIMAP2=%s  MITOS2_REFDIR=%s\n' \
  "${MINIMAP2:-(未配置；将回退到 PATH/envs)}" "${MITOS2_REFDIR:-(未配置)}"

echo "=== 机器 profile: $PROFILE ==="
echo "  CONDA_ROOT=$CONDA_ROOT"
echo "  MITOS2_PY=$MITOS2_PY"
echo "  MITOS2_REFDIR=$MITOS2_REFDIR (refseqver=$MITOS2_REFSEQVER)"
echo "  PLOT_PY=$PLOT_PY"
echo "  MINIMAP2=$MINIMAP2"

echo ""
echo "=== 核心工具 ==="
for t in blastn makeblastdb bwa samtools seqkit python3; do
  found=""
  if command -v $t >/dev/null 2>&1; then found=$(command -v $t)
  else
    for env in "$CONDA_ROOT"/envs/*/bin; do
      [ -x "$env/$t" ] && found="$env/$t" && break
    done
  fi
  if [ -n "$found" ]; then echo "  ✓ $t: $found"
  else echo "  ✗ $t: 未找到"; missing=1; fi
done
# minimap2: PATH -> 已配置路径 -> envs
mini=""
if command -v minimap2 >/dev/null 2>&1; then mini=$(command -v minimap2)
elif [ -n "${MINIMAP2:-}" ] && [ -x "$MINIMAP2" ]; then mini="$MINIMAP2 (config/env.sh)"
else
  for env in "$CONDA_ROOT"/envs/*/bin; do
    [ -x "$env/minimap2" ] && mini="$env/minimap2" && break
  done
fi
if [ -n "$mini" ]; then echo "  ✓ minimap2: $mini"; else echo "  ✗ minimap2: 未找到"; missing=1; fi

echo ""
echo "=== 组装/注释工具 (在 conda envs 中查找) ==="
while IFS='|' read -r t description; do
  [ -z "$t" ] && continue
  found=""
  if command -v "$t" >/dev/null 2>&1; then found="PATH: $(command -v "$t")"; else
    for env in "$CONDA_ROOT"/envs/*/bin; do
      [ -x "$env/$t" ] && found="${env#$CONDA_ROOT/envs/}" && break
    done
  fi
  if [ -n "$found" ]; then echo "  ✓ $t -> $found ($description)"
  else echo "  ✗ $t ($description) 未找到"; fi
done <<'TOOLS'
get_organelle_from_reads.py|getorganelle env
mitofinder|mitofinder env
cmsearch|MITOS2 依赖 (infernal)
RNAplot|MITOS2 画图 (ViennaRNA)
tRNAscan-SE|genome_assembly env (辅助)
arwen|mitofinder env (tRNA)
TOOLS

# MITOS2 本体: 用 config 里的 python 验证 mitos 模块可导入
echo ""
echo "=== MITOS2 python 模块 (config/env.sh MITOS2_PY) ==="
if [ -n "${MITOS2_PY:-}" ] && [ -x "$MITOS2_PY" ]; then
  if "$MITOS2_PY" -c "import mitos, sys; print('  ✓', sys.executable.split('/')[-3], 'env: mitos 模块 OK')" >/dev/null 2>&1; then
    echo "  ✓ $MITOS2_PY (mitos 模块可导入)"
  else
    echo "  ✗ $MITOS2_PY 无法导入 mitos 模块 — 检查 config/env.sh 的 MITOS2_PY"
  fi
else
  echo "  ✗ MITOS2_PY 未设置或不存在 — 检查 config/env.sh"
fi

echo ""
echo "=== Python 库依赖 ==="
# 案例校验只有 jsonschema 一条实现路径（曾经的手写兜底已删除：同一结论不允许两条证据路径）。
# 缺它时 `case-validate` / `case-anomaly` / `case-reference` 会以退出码 3 失败，不会静默降级。
lib_missing=0
for lib in jsonschema; do
  if python3 -c "import $lib" >/dev/null 2>&1; then
    ver=$(python3 -c "import importlib.metadata as m; print(m.version('$lib'))" 2>/dev/null || echo '已安装')
    echo "  ✓ $lib: $ver（案例校验的唯一实现）"
  else
    echo "  ✗ $lib: 未安装 —— 案例校验无法工作（会以退出码 3 失败，不降级）"
    echo "    安装: python3 -m pip install $lib"
    lib_missing=1
  fi
done
# Biopython：注释与序列脚本使用；各脚本的依赖边界正在逐脚本审查（见仓库 issue），
# 因此这里只声明与提示，不阻断（不预判那个边界结论）。
if python3 -c "import Bio" >/dev/null 2>&1; then
  echo "  ✓ Biopython: 已安装（annot_check / blast_genes / circularize / mitos2_to_genbank / seq_stats 需要）"
else
  echo "  ○ Biopython: 未安装 —— 上述脚本无法工作；安装: python3 -m pip install biopython"
fi
if [ "$lib_missing" -ne 0 ]; then missing=1; fi

if [ "$missing" -ne 0 ]; then
  if [ -n "$env_unset_hint" ]; then
    echo "\n结果: 核心依赖缺失 —— 但 $env_unset_hint 未设置，本轮“未找到”可能只是假阴性；请先导出环境变量再重跑" >&2
  else
    echo "\n结果: 核心依赖缺失（环境变量已设置，未找到即为真实缺失）" >&2
  fi
  exit 2
fi

echo ""
echo "=== MITOS2 数据库 ==="
if [ -d "$MITOS2_REFDIR/$MITOS2_REFSEQVER" ]; then
  echo "  ✓ $MITOS2_REFSEQVER: $MITOS2_REFDIR/$MITOS2_REFSEQVER"
else
  echo "  ✗ $MITOS2_REFSEQVER 未找到: $MITOS2_REFDIR/$MITOS2_REFSEQVER"
  echo "    (跑 MITOS2 注释前需就位; 跨机迁移可 rsync 旧机该目录)"
fi

echo ""
echo "=== MITOS2 额外 PATH (config/env.sh MITOS2_EXTRA_PATH) ==="
if [ -n "${MITOS2_EXTRA_PATH:-}" ]; then echo "  $MITOS2_EXTRA_PATH"; else echo "  (空 — 本机 MITOS2 环境自包含, 无需额外 PATH)"; fi

#!/bin/bash
# 检查线粒体组装诊断流水线所需的环境与工具
# 机器路径统一在 config/env.sh（本脚本自动 source），换机无需改本文件
# 用法: bash check_env.sh
set -u
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=../config/env.sh
source "$SKILL_DIR/config/env.sh"

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
  else echo "  ✗ $t: 未找到"; fi
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
if [ -n "$mini" ]; then echo "  ✓ minimap2: $mini"; else echo "  ✗ minimap2: 未找到"; fi

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

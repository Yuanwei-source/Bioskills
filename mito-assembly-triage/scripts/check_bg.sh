#!/bin/bash
# 检查后台任务状态
# 用法: bash check_bg.sh <任务名> [--tail 20]
set -u
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOGDIR="$SKILL_DIR/logs"

if [ $# -lt 1 ]; then
  echo "用法: bash check_bg.sh <任务名> [--tail N]"
  echo "可用任务: $(ls "$LOGDIR"/*.status 2>/dev/null | xargs -n1 basename 2>/dev/null | sed 's/.status//' | tr '\n' ' ')"
  exit 1
fi
NAME="$1"
if ! [[ "$NAME" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$ ]]; then
  echo "✗ 任务名只能包含字母、数字、点、下划线和短横线"
  exit 1
fi
TAIL=20
if [ "$#" -ne 1 ] && [ "$#" -ne 3 ]; then
  echo "用法: bash check_bg.sh <任务名> [--tail N]"
  exit 1
fi
if [ "$#" -eq 3 ]; then
  if [ "$2" != "--tail" ] || ! [[ "$3" =~ ^[1-9][0-9]*$ ]]; then
    echo "✗ --tail 必须跟正整数"
    exit 1
  fi
  TAIL="$3"
fi

if [ ! -f "$LOGDIR/$NAME.status" ]; then
  echo "✗ 任务 $NAME 不存在 (未启动或日志已清理)"
  exit 1
fi

STATUS=$(cat "$LOGDIR/$NAME.status")
PID=$(cat "$LOGDIR/$NAME.pid" 2>/dev/null || echo "-")

# 双重确认: status 文件 + 进程存活
if [ "$STATUS" = "running" ]; then
  if [ -n "$PID" ] && kill -0 "$PID" 2>/dev/null; then
    echo "⏳ [$NAME] 运行中 (PID $PID)"
  else
    # 守护进程可能还没更新状态
    echo "⏳ [$NAME] 状态未定 (PID $PID 已退出, 等待状态更新)"
  fi
elif [ "$STATUS" = "done" ]; then
  echo "✅ [$NAME] 已完成"
elif [ "$STATUS" = "failed" ]; then
  echo "❌ [$NAME] 失败 — 查看日志尾部定位错误"
fi

echo "--- 日志尾部 (最后 $TAIL 行) ---"
tail -n "$TAIL" "$LOGDIR/$NAME.log" 2>/dev/null

# 输出文件检查提示
if [ "$STATUS" = "done" ]; then
  echo ""
  echo "提示: 检查输出文件是否生成且非空 (注意: done ≠ 结果正确, 需按 Quality Gates 验证)"
fi

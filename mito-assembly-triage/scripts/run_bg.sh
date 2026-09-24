#!/bin/bash
# 后台运行长任务 (防止 AI 前端超时)
# 用法: bash run_bg.sh <任务名> "<命令>"
# 输出: logs/<任务名>.log    - 程序输出
#       logs/<任务名>.pid    - 进程 PID
#       logs/<任务名>.status - running / done / failed
set -u
set -o pipefail
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOGDIR="$SKILL_DIR/logs"
mkdir -p "$LOGDIR"

if [ "$#" -lt 3 ]; then
  echo "用法: bash run_bg.sh <任务名> -- <命令> [参数...]"
  echo "      bash run_bg.sh <任务名> --trusted-shell \"<已确认可信的 shell 命令>\""
  exit 1
fi
NAME="$1"
MODE="$2"
if ! [[ "$NAME" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$ ]]; then
  echo "✗ 任务名只能包含字母、数字、点、下划线和短横线"
  exit 1
fi
case "$MODE" in
  --)
    shift 2
    if [ "$#" -lt 1 ]; then
      echo "✗ 缺少命令"
      exit 1
    fi
    COMMAND=("$@")
    TRUSTED_SHELL=0
    ;;
  --trusted-shell)
    shift 2
    if [ "$#" -ne 1 ]; then
      echo "用法: bash run_bg.sh <任务名> --trusted-shell \"<命令>\""
      exit 1
    fi
    COMMAND=("$1")
    TRUSTED_SHELL=1
    ;;
  *)
    echo "✗ 第二个参数必须是 -- 或 --trusted-shell"
    exit 1
    ;;
esac
PIDFILE="$LOGDIR/$NAME.pid"
STATUSFILE="$LOGDIR/$NAME.status"
LOGFILE="$LOGDIR/$NAME.log"
LOCKFILE="$LOGDIR/$NAME.lock"
exec 9>"$LOCKFILE"
if ! flock -n 9; then
  echo "⚠ 任务 $NAME 已有启动器持有锁"; exit 1
fi

# 同名任务已在运行则拒绝
if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "⚠ 任务 $NAME 已在运行 (PID $(cat "$PIDFILE"))"
  exit 1
fi

atomic_write() { local value="$1" target="$2" tmp; tmp="${target}.tmp.$$"; printf '%s\n' "$value" > "$tmp" && mv -f "$tmp" "$target"; }
atomic_write running "$STATUSFILE"
echo "=== 启动时间: $(date) ===" >> "$LOGFILE"

run_and_monitor() {
  if [ "$TRUSTED_SHELL" -eq 1 ]; then
    nohup bash -c "${COMMAND[0]}" >> "$LOGFILE" 2>&1 &
  else
    nohup "${COMMAND[@]}" >> "$LOGFILE" 2>&1 &
  fi
  PID=$!
  atomic_write "$PID" "$PIDFILE"
  while kill -0 "$PID" 2>/dev/null; do sleep 5; done
  wait "$PID" 2>/dev/null
  RC=$?
  if [ "$RC" -eq 0 ]; then
    atomic_write done "$STATUSFILE"
    echo "=== 完成时间: $(date) (exit 0) ===" >> "$LOGFILE"
  else
    atomic_write failed "$STATUSFILE"
    echo "=== 失败时间: $(date) (exit $RC) ===" >> "$LOGFILE"
  fi
}

run_and_monitor </dev/null >> /dev/null 2>&1 &

echo "✓ 后台任务 [$NAME] 已启动"
echo "  日志: $LOGFILE"
echo "  检查: bash scripts/check_bg.sh $NAME"

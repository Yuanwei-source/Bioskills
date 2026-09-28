#!/usr/bin/env bash
# FASTA-only route: read-only intake -> MITOS2 -> GenBank bridge -> annotation QC.
set -euo pipefail
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
usage() {
  cat <<'EOF'
用法: bash scripts/run_fasta_annotation.sh --fasta assembly.fa --outdir NEW_DIR --table N [--table-status unverified|provisional|confirmed] [--organism NAME] [--taxid N] [--locus ID] [--topology linear|circular] [--mitos-dir EXISTING_RESULT_DIR]

输出到全新的任务目录：intake.json、MITOS2 原始输出、annotation.gb、annot_check.log、route.json。
--table 必须与 MITOS2 一致；默认记录为 unverified，已确认时用 --table-status confirmed，暂定时用 provisional。
--organism/--taxid 仅按调用者提供值写入 source，不自动推断。--mitos-dir 可复用已有 result.gff/fas/faa，避免重复运行耗时注释。
该流程要求单条 FASTA 记录；不改输入，不推断碱基、不证明闭环。MITOS2 结果需继续人工/证据复核。
EOF
}
fasta=""; outdir=""; table=""; topology="linear"; table_status="unverified"; organism=""; taxid=""; locus=""; reuse_mitos_dir=""
while (($#)); do
  case "$1" in
    --fasta) (($# >= 2)) || { usage >&2; exit 1; }; fasta="$2"; shift 2 ;;
    --outdir) (($# >= 2)) || { usage >&2; exit 1; }; outdir="$2"; shift 2 ;;
    --table) (($# >= 2)) || { usage >&2; exit 1; }; table="$2"; shift 2 ;;
    --table-status) (($# >= 2)) || { usage >&2; exit 1; }; table_status="$2"; shift 2 ;;
    --organism) (($# >= 2)) || { usage >&2; exit 1; }; organism="$2"; shift 2 ;;
    --taxid) (($# >= 2)) || { usage >&2; exit 1; }; taxid="$2"; shift 2 ;;
    --locus) (($# >= 2)) || { usage >&2; exit 1; }; locus="$2"; shift 2 ;;
    --mitos-dir) (($# >= 2)) || { usage >&2; exit 1; }; reuse_mitos_dir="$2"; shift 2 ;;
    --topology) (($# >= 2)) || { usage >&2; exit 1; }; topology="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "未知参数: $1" >&2; usage >&2; exit 1 ;;
  esac
done
[[ -n "$fasta" && -n "$outdir" && "$table" =~ ^[0-9]+$ ]] || { usage >&2; exit 1; }
[[ "$table_status" == unverified || "$table_status" == provisional || "$table_status" == confirmed ]] || { echo "--table-status 只能是 unverified、provisional 或 confirmed" >&2; exit 1; }
[[ -z "$taxid" || "$taxid" =~ ^[1-9][0-9]*$ ]] || { echo "--taxid 必须是正整数" >&2; exit 1; }
[[ "$topology" == linear || "$topology" == circular ]] || { echo "--topology 只能是 linear 或 circular" >&2; exit 1; }
[[ -f "$fasta" ]] || { echo "FASTA 不存在: $fasta" >&2; exit 1; }
[[ ! -e "$outdir" ]] || { echo "为保护旧结果，输出目录必须不存在: $outdir" >&2; exit 1; }
mkdir -p "$outdir"
fasta="$(realpath "$fasta")"; outdir="$(realpath "$outdir")"

python3 "$SKILL_DIR/scripts/assembly_intake.py" "$fasta" --json "$outdir/intake.json" > "$outdir/intake.txt"
record_count="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["record_count"])' "$outdir/intake.json")"
if [[ "$record_count" != 1 ]]; then
  echo "STOP: 输入有 $record_count 条记录。请先判断记录间的生物学关系；本路线的 MITOS2→GenBank 桥要求单记录。初检保存在 $outdir/intake.json" >&2
  exit 2
fi

mitos_source="executed"
if [[ -n "$reuse_mitos_dir" ]]; then
  mitos_dir="$(realpath "$reuse_mitos_dir")"
  [[ -d "$mitos_dir" ]] || { echo "MITOS2 结果目录不存在: $mitos_dir" >&2; exit 1; }
  mitos_source="reused"
  printf '复用已存在 MITOS2 结果目录（转换器将校验其与本次 FASTA 的序列/坐标对应关系）: %s\n' "$mitos_dir" > "$outdir/mitos2.log"
else
  mitos_dir="$outdir/mitos2"
  if bash "$SKILL_DIR/scripts/run_mitos2.sh" -i "$fasta" -o "$mitos_dir" -c "$table" > "$outdir/mitos2.log" 2>&1; then
    :
  else
    rc=$?
    echo "MITOS2 未完成（退出码 $rc）；查看 $outdir/mitos2.log。未生成或验收注释。" >&2
    exit "$rc"
  fi
fi
for required in result.gff result.fas result.faa; do
  [[ -s "$mitos_dir/$required" ]] || { echo "MITOS2 输出缺失/为空: $mitos_dir/$required" >&2; exit 1; }
done
bridge_args=()
[[ -z "$organism" ]] || bridge_args+=(--organism "$organism")
[[ -z "$taxid" ]] || bridge_args+=(--taxid "$taxid")
[[ -z "$locus" ]] || bridge_args+=(--locus "$locus")
python3 "$SKILL_DIR/scripts/mitos2_to_genbank.py" \
  "$mitos_dir/result.gff" "$mitos_dir/result.fas" "$mitos_dir/result.faa" \
  "$outdir/annotation.gb" --genome "$fasta" --table "$table" --topology "$topology" \
  "${bridge_args[@]}" \
  > "$outdir/bridge.log" 2>&1
set +e
python3 "$SKILL_DIR/scripts/annot_check.py" "$outdir/annotation.gb" --table "$table" \
  > "$outdir/annot_check.log" 2>&1
qc_rc=$?
set -e
python3 - "$outdir/route.json" "$fasta" "$table" "$table_status" "$topology" "$qc_rc" "$organism" "$taxid" "$mitos_source" "$mitos_dir" <<'PY'
import hashlib, json, pathlib, sys
out, fasta, table, table_status, topology, qc_rc, organism, taxid, mitos_source, mitos_dir = sys.argv[1:]
p = pathlib.Path(fasta)
data = {"route": "fasta_only_mitos2_qc", "input_fasta": str(p.resolve()),
        "input_sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        "genetic_code_table": int(table), "genetic_code_status": table_status,
        "organism": organism or None, "taxid": int(taxid) if taxid else None,
        "mitos2_result_source": mitos_source, "mitos2_result_dir": str(pathlib.Path(mitos_dir).resolve()),
        "topology_declaration": topology,
        "annotation_check_exit_code": int(qc_rc),
        "interpretation": {"0": "未发现本工具检查项", "1": "存在错误", "2": "存在待核查项"}.get(qc_rc, "工具执行异常"),
        "limitations": ["MITOS2 输出及其 GenBank 转换均为注释候选，不等同专家终审或提交质量。",
                        "不验证样本碱基、NUMT、连接 reads 支持或物理闭环。",
                        "需结合类群证据复核基因身份、边界、例外及初检暴露的问题。"]}
pathlib.Path(out).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY
echo "FASTA-only 注释候选已生成: $outdir/annotation.gb"
echo "初检: $outdir/intake.txt；MITOS2 日志: $outdir/mitos2.log；转换日志: $outdir/bridge.log"
echo "注释复核退出码: $qc_rc（0=无本工具发现，2=待核查，1=错误）；详见 $outdir/annot_check.log"
echo "路线元数据: $outdir/route.json。继续执行 references/START_HERE.md 的联合复核与任务报告步骤。"
exit "$qc_rc"

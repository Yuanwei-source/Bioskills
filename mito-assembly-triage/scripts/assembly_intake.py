#!/usr/bin/env python3
"""Read-only intake for an assembled mitochondrial FASTA.

Reports record-level facts and unresolved sequence symbols. It does not infer
correct bases, gene content, orientation, or physical circularity.
"""
import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import os as _os
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from _deps import require_stage
from _evidence_io import InputArgumentParser, protect_output, atomic_write_text

IUPAC = set("ACGTRYSWKMBDHVN")


def parse_fasta(path):
    records = []
    name = None
    seq = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            header = line[1:].strip()
            if not header:
                raise ValueError("第 %d 行 FASTA 标题为空" % number)
            name, seq = header.split()[0], []
        else:
            if name is None:
                raise ValueError("第 %d 行序列出现在 FASTA 标题之前" % number)
            bases = "".join(line.split()).upper()
            invalid = sorted(set(bases) - IUPAC)
            if invalid:
                raise ValueError("第 %d 行含非法 DNA/IUPAC 符号: %s" % (number, "".join(invalid)))
            seq.append(bases)
    if name is not None:
        records.append((name, "".join(seq)))
    if not records:
        raise ValueError("未发现 FASTA 序列记录")
    if any(not sequence for _, sequence in records):
        raise ValueError("至少一条 FASTA 记录没有序列")
    return records


def intervals(sequence, symbols):
    found = []
    start = None
    for i, base in enumerate(sequence):
        if base in symbols and start is None:
            start = i
        elif base not in symbols and start is not None:
            found.append([start + 1, i])
            start = None
    if start is not None:
        found.append([start + 1, len(sequence)])
    return found


def inspect(path):
    raw = path.read_bytes()
    records = parse_fasta(path)
    ids = [name for name, _ in records]
    id_counts = Counter(ids)
    duplicate_ids = sorted(name for name, count in id_counts.items() if count > 1)
    report_records = []
    for name, sequence in records:
        counts = Counter(sequence)
        unambiguous = sum(counts[b] for b in "ACGT")
        gc = (counts["G"] + counts["C"]) / unambiguous * 100 if unambiguous else None
        ambiguous = intervals(sequence, set(sequence) - set("ACGT"))
        n_runs = intervals(sequence, {"N"})
        report_records.append({
            "id": name,
            "length_bp": len(sequence),
            "gc_percent_unambiguous": round(gc, 4) if gc is not None else None,
            "ambiguous_base_count": len(sequence) - unambiguous,
            "ambiguous_counts": {k: counts[k] for k in sorted(counts) if k not in "ACGT"},
            "ambiguous_intervals_1based_inclusive": ambiguous,
            "n_intervals_1based_inclusive": n_runs,
        })
    findings = []
    if len(records) > 1:
        findings.append({"level": "REVIEW", "code": "MULTIPLE_RECORDS",
                         "message": "输入含多条记录；应逐条确认它们是独立 contig、候选 haplotype 还是误拼接，不要直接串接。"})
    if duplicate_ids:
        findings.append({"level": "BLOCK", "code": "DUPLICATE_IDS",
                         "message": "FASTA 记录 ID 重复，后续坐标与注释无法可靠绑定。", "ids": duplicate_ids})
    if any(row["ambiguous_base_count"] for row in report_records):
        findings.append({"level": "REVIEW", "code": "AMBIGUOUS_BASES",
                         "message": "模糊符号位置已列出；只有样本 reads 或其他独立证据支持时才能替换。"})
    if not findings:
        findings.append({"level": "INFO", "code": "BASIC_INTAKE_CLEAR",
                         "message": "未发现多记录、模糊碱基或格式问题；这不证明组装序列正确或分子已物理闭合。"})
    return {"schema_version": 1, "input": str(path.resolve()),
            "sha256": hashlib.sha256(raw).hexdigest(), "record_count": len(records),
            "total_bp_separate_records": sum(len(seq) for _, seq in records),
            "records": report_records, "findings": findings,
            "scope_limitations": ["只报告 FASTA 可直接观察的事实。",
                                  "不验证样本碱基、基因身份/边界、方向、来源或物理闭环。",
                                  "不因无 reads 而阻止可用序列区域的注释；含模糊位点或多 contig 的注释须标注限制。"]}


def main(argv=None):
    parser = InputArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--json", dest="json_path", type=Path,
                        help="另存结构化 intake JSON；文件不得覆盖输入或已有文件")
    args = parser.parse_args(argv)
    try:
        require_stage("assembly_intake", __file__)
        if not args.fasta.is_file():
            raise ValueError("FASTA 文件不存在: %s" % args.fasta)
        report = inspect(args.fasta)
        if args.json_path:
            protect_output(args.json_path, [args.fasta])
            if args.json_path.exists():
                raise ValueError("JSON 输出已存在，为避免覆盖请另选路径: %s" % args.json_path)
            atomic_write_text(args.json_path, json.dumps(report, ensure_ascii=False, indent=2) + "\n")
        print("FASTA: %s" % report["input"])
        print("SHA-256: %s" % report["sha256"])
        print("记录数: %d；各记录长度合计: %d bp（未串接）" %
              (report["record_count"], report["total_bp_separate_records"]))
        for row in report["records"]:
            print("- %s: %d bp; GC(ACGT only)=%s; ambiguous=%d" %
                  (row["id"], row["length_bp"],
                   ("%.2f%%" % row["gc_percent_unambiguous"])
                   if row["gc_percent_unambiguous"] is not None else "NA",
                   row["ambiguous_base_count"]))
            for lo, hi in row["ambiguous_intervals_1based_inclusive"]:
                print("  ambiguous %d-%d (1-based inclusive)" % (lo, hi))
        for finding in report["findings"]:
            print("[%s] %s: %s" % (finding["level"], finding["code"], finding["message"]))
        if args.json_path:
            print("JSON: %s" % args.json_path.resolve())
        return 0
    except (OSError, UnicodeError, ValueError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Create and validate small, auditable V2 diagnostic cases.

Examples:
  v2_case.py init CASE_DIR --issue internal_stop --observation "nad5 has stop"
  v2_case.py event CASE_DIR --action seq_stats --result "two contigs" --impact H1:against
  v2_case.py report CASE_DIR
  v2_case.py validate CASE_DIR
"""
import argparse, datetime, hashlib, json, os, pathlib, sys

STATUSES = {"RESOLVED", "NO_CHANGE", "UNRESOLVED"}
CONFIDENCE = {"high", "moderate", "low", "not_assessable"}

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()

def safe_dir(path):
    p = pathlib.Path(path).resolve(); p.mkdir(parents=True, exist_ok=True); return p

def load_case(root):
    with open(root / "case.json", encoding="utf-8") as fh: return json.load(fh)

def save_case(root, case):
    tmp = root / "case.json.tmp"
    with open(tmp, "w", encoding="utf-8") as fh: json.dump(case, fh, ensure_ascii=False, indent=2); fh.write("\n")
    os.replace(tmp, root / "case.json")

def init(args):
    root = safe_dir(args.directory)
    inputs = []
    for role, raw in args.input:
        path = pathlib.Path(raw).resolve()
        if not path.is_file(): raise SystemExit("输入不存在: %s" % raw)
        inputs.append({"role": role, "path": str(path), "sha256": sha256(path)})
    case = {"schema_version": "2.0", "case_id": args.case_id or root.name,
            "taxon": {"name": args.taxon, "taxid": None, "genetic_code": None},
            "inputs": inputs, "issue": {"type": args.issue, "user_observation": args.observation},
            "hypotheses": [], "events_file": "events.jsonl",
            "decision": {"status": "UNRESOLVED", "confidence": "not_assessable", "rationale": "尚未完成足够检查"},
            "modifications": [], "validation": [], "lessons_proposed": []}
    save_case(root, case)
    (root / "events.jsonl").touch()
    print("case initialized: %s" % (root / "case.json"))

def event(args):
    root = pathlib.Path(args.directory).resolve(); case = load_case(root)
    record = {"timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(), "action": args.action,
              "command": args.command, "tool_versions": args.tool_version, "result": args.result,
              "motivation": args.motivation, "impact": args.impact}
    with open(root / case["events_file"], "a", encoding="utf-8") as fh: fh.write(json.dumps(record, ensure_ascii=False) + "\n")
    print("event recorded")

def validate(args):
    root = pathlib.Path(args.directory).resolve(); case = load_case(root); errors=[]
    if case.get("schema_version") != "2.0": errors.append("schema_version must be 2.0")
    d=case.get("decision", {})
    if d.get("status") not in STATUSES: errors.append("invalid decision.status")
    if d.get("confidence") not in CONFIDENCE: errors.append("invalid decision.confidence")
    for item in case.get("inputs", []):
        if not item.get("path") or len(item.get("sha256", "")) != 64: errors.append("input requires path and sha256")
    if not (root / case.get("events_file", "events.jsonl")).is_file(): errors.append("events file missing")
    if errors:
        print("INVALID"); [print("- " + e) for e in errors]; return 1
    print("VALID"); return 0

def report(args):
    root=pathlib.Path(args.directory).resolve(); c=load_case(root)
    lines=["# 诊断案例 %s" % c["case_id"], "", "## 事件", ""]
    with open(root / c["events_file"], encoding="utf-8") as fh:
        for line in fh:
            e=json.loads(line); lines.append("- `%s`：%s；影响：%s" % (e.get("action"), e.get("result"), e.get("impact") or "未记录"))
    lines += ["", "## 当前判定", "", "- 状态：`%s`" % c["decision"]["status"], "- 置信等级：`%s`" % c["decision"]["confidence"], "- 理由：%s" % c["decision"]["rationale"], "", "## 证据边界", "", "仅当事件明确记录 raw reads、参考、软件或注释来源时，报告才引用相应证据；缺少 reads 不得写成 raw-read-supported。"]
    (root / "case.md").write_text("\n".join(lines) + "\n", encoding="utf-8"); print("report written")

def main():
    ap=argparse.ArgumentParser(description=__doc__); sub=ap.add_subparsers(dest="cmd", required=True)
    p=sub.add_parser("init"); p.add_argument("directory"); p.add_argument("--case-id"); p.add_argument("--issue", required=True); p.add_argument("--observation", required=True); p.add_argument("--taxon"); p.add_argument("--input", nargs=2, action="append", default=[], metavar=("ROLE", "PATH")); p.set_defaults(fn=init)
    p=sub.add_parser("event"); p.add_argument("directory"); p.add_argument("--action", required=True); p.add_argument("--result", required=True); p.add_argument("--impact", default=""); p.add_argument("--command", default=""); p.add_argument("--tool-version", default=""); p.add_argument("--motivation", default=""); p.set_defaults(fn=event)
    p=sub.add_parser("validate"); p.add_argument("directory"); p.set_defaults(fn=validate)
    p=sub.add_parser("report"); p.add_argument("directory"); p.set_defaults(fn=report)
    args=ap.parse_args(); result=args.fn(args); sys.exit(result or 0)
if __name__ == "__main__": main()

#!/usr/bin/env python
"""
Inspect one BioProject accession and aggregate project-level context for BioFinder.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from typing import Any

from evidence_schema import build_bioproject_inspection_evidence
from inspect_gse import normalize_goal
from search_bioproject import BioProjectRecord, fetch_records
from search_geo import configure_stdio
from search_sra import fetch_records as fetch_sra_records


def unique_nonempty(values: list[str]) -> list[str]:
    seen: set[str] = set()
    cleaned: list[str] = []
    for value in values:
        text = str(value).strip()
        if not text or text in seen:
            continue
        seen.add(text)
        cleaned.append(text)
    return cleaned


def infer_tissues(records: list[BioProjectRecord]) -> list[str]:
    tissues: list[str] = []
    for record in records:
        haystack = f"{record.title} {record.summary}".lower()
        if "abdominal fat" in haystack:
            tissues.append("abdominal fat")
        elif "adipose" in haystack:
            tissues.append("adipose")
        elif "cartilage" in haystack:
            tissues.append("cartilage")
        elif "islet" in haystack:
            tissues.append("islet")
    return unique_nonempty(tissues)


def infer_risks(records: list[BioProjectRecord], linked_sra: list[str], linked_biosamples: list[str]) -> list[str]:
    risks: list[str] = []
    if len(records) > 1:
        risks.append("multiple_project_records")
    if len(linked_sra) > 20:
        risks.append("broad_project_scope")
    if not linked_sra:
        risks.append("linked_sra_unclear")
    if not linked_biosamples:
        risks.append("linked_biosamples_unclear")
    return unique_nonempty(risks)


def assess_analysis_fit(analysis_goal: str, linked_sra: list[str], risks: list[str]) -> tuple[list[str], list[str], list[str], str]:
    usable_for = ["project_context"]
    not_recommended_for: list[str] = []
    rationale: list[str] = []
    if linked_sra:
        usable_for.append("sra_follow_up")
    if "broad_project_scope" in risks:
        rationale.append("project_may_be_too_broad_for_direct_shortlist_use")
    if analysis_goal == "wgcna":
        not_recommended_for.append("wgcna_direct_use")
        rationale.append("bioproject_is_context_not_expression_matrix")
    overall = "medium" if linked_sra else "low"
    if "broad_project_scope" in risks:
        overall = "low"
    return unique_nonempty(usable_for), unique_nonempty(not_recommended_for), unique_nonempty(rationale), overall


@dataclass
class BioProjectInspectionResult:
    accession: str
    title: str
    organisms: list[str]
    tissues: list[str]
    linked_sra: list[str]
    linked_geo: list[str]
    linked_biosamples: list[str]
    run_count: int | None
    risks: list[str]
    usable_for: list[str]
    not_recommended_for: list[str]
    rationale: list[str]
    overall_usability: str
    data_availability: dict[str, Any]
    summary: str
    raw_notes: dict[str, Any]

    def to_dict(self, include_raw: bool) -> dict[str, Any]:
        data = {
            "accession": self.accession,
            "title": self.title,
            "organisms": self.organisms,
            "tissues": self.tissues,
            "linked_sra": self.linked_sra,
            "linked_geo": self.linked_geo,
            "linked_biosamples": self.linked_biosamples,
            "run_count": self.run_count,
            "risks": self.risks,
            "usable_for": self.usable_for,
            "not_recommended_for": self.not_recommended_for,
            "rationale": self.rationale,
            "provisional_local_assessment": self.overall_usability,
            "assessment_source": "local_bioproject_context_aggregation",
            "data_availability": self.data_availability,
            "summary": self.summary,
            "normalized_evidence": build_bioproject_inspection_evidence(self, include_summary=True),
        }
        if include_raw:
            data["raw_notes"] = self.raw_notes
        return data


def inspect_accession(accession: str, analysis_goal: str, timeout: int, pause_seconds: float, email: str | None) -> BioProjectInspectionResult:
    records = fetch_records(accession, max_results=20, timeout=timeout, pause_seconds=pause_seconds, email=email)
    project_records = [record for record in records if (record.accession or "").upper() == accession.strip().upper()]
    if project_records:
        records = project_records

    title = next((record.title for record in records if record.title), "")
    summary = next((record.summary for record in records if record.summary), "")
    organisms = unique_nonempty([record.organism for record in records if record.organism])
    tissues = infer_tissues(records)
    linked_sra = unique_nonempty([item for record in records for item in record.linked_sra])
    linked_geo = unique_nonempty([])
    linked_biosamples = unique_nonempty([item for record in records for item in record.linked_biosamples])
    sra_backfill_count = 0
    if not linked_sra:
        try:
            sra_records = fetch_sra_records(accession, max_results=50, timeout=timeout, pause_seconds=pause_seconds, email=email)
        except Exception:  # noqa: BLE001
            sra_records = []
        if sra_records:
            sra_backfill_count = len(sra_records)
            linked_sra = unique_nonempty(
                [record.accession for record in sra_records]
                + [record.study_accession for record in sra_records]
                + [record.experiment_accession for record in sra_records]
                + [record.sample_accession for record in sra_records]
                + [run for record in sra_records for run in record.run_accessions]
            )
            linked_biosamples = unique_nonempty(linked_biosamples + [record.biosample for record in sra_records if record.biosample])
            organisms = unique_nonempty(organisms + [record.organism for record in sra_records if record.organism])
    risks = infer_risks(records, linked_sra, linked_biosamples)
    usable_for, not_recommended_for, rationale, overall = assess_analysis_fit(
        analysis_goal=analysis_goal,
        linked_sra=linked_sra,
        risks=risks,
    )
    run_count = len([item for item in linked_sra if item.startswith("SRR")]) or None
    data_availability = {
        "linked_sra_count": len(linked_sra),
        "linked_biosample_count": len(linked_biosamples),
        "run_accessions_available": bool(run_count),
        "needs_manual_check": not bool(title or summary),
    }
    raw_notes = {
        "requested_analysis_goal": analysis_goal,
        "record_count": len(records),
        "project_titles": unique_nonempty([record.title for record in records if record.title]),
        "sra_backfill_count": sra_backfill_count,
    }
    return BioProjectInspectionResult(
        accession=accession,
        title=title,
        organisms=organisms,
        tissues=tissues,
        linked_sra=linked_sra,
        linked_geo=linked_geo,
        linked_biosamples=linked_biosamples,
        run_count=run_count,
        risks=risks,
        usable_for=usable_for,
        not_recommended_for=not_recommended_for,
        rationale=rationale,
        overall_usability=overall,
        data_availability=data_availability,
        summary=summary,
        raw_notes=raw_notes,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect BioProject accessions for BioFinder.")
    parser.add_argument("--accession", action="append", required=True, help="BioProject accession to inspect. Repeatable.")
    parser.add_argument("--analysis-goal", default="general", help="Analysis goal such as DEG, WGCNA, or validation.")
    parser.add_argument("--timeout", type=int, default=30, help="Network timeout in seconds.")
    parser.add_argument("--pause-seconds", type=float, default=0.34, help="Pause between NCBI requests.")
    parser.add_argument("--email", help="Email to send to NCBI as recommended by E-utilities.")
    parser.add_argument("--json", action="store_true", help="Print JSON output.")
    parser.add_argument("--include-raw", action="store_true", help="Include raw aggregation notes.")
    return parser.parse_args()


def render_text(results: list[BioProjectInspectionResult]) -> str:
    blocks: list[str] = []
    for result in results:
        lines = [
            f"{result.accession}",
            f"   title={result.title or 'n/a'}",
            f"   organisms={', '.join(result.organisms) or 'n/a'}",
            f"   tissues={', '.join(result.tissues) or 'n/a'}",
            f"   linked_sra_count={len(result.linked_sra)}",
            f"   linked_biosample_count={len(result.linked_biosamples)}",
            f"   provisional_local_assessment={result.overall_usability}",
            f"   risks={', '.join(result.risks) or 'n/a'}",
        ]
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks)


def main() -> int:
    configure_stdio()
    args = parse_args()
    analysis_goal = normalize_goal(args.analysis_goal)
    results = [
        inspect_accession(
            accession=item,
            analysis_goal=analysis_goal,
            timeout=args.timeout,
            pause_seconds=args.pause_seconds,
            email=args.email,
        )
        for item in args.accession
    ]
    if args.json:
        payload = {
            "analysis_goal": analysis_goal,
            "result_count": len(results),
            "results": [result.to_dict(include_raw=args.include_raw) for result in results],
        }
        json.dump(payload, sys.stdout, indent=2, ensure_ascii=False)
        sys.stdout.write("\n")
    else:
        print(render_text(results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

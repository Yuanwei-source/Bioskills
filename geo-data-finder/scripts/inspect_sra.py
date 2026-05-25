#!/usr/bin/env python
"""
Inspect one SRA accession and aggregate study-level evidence for BioFinder.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from evidence_schema import build_sra_inspection_evidence
from inspect_gse import normalize_goal
from search_geo import configure_stdio
from search_sra import SraRecord, fetch_records


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


def choose_anchor(records: list[SraRecord], requested_accession: str) -> str:
    requested = requested_accession.strip().upper()
    for record in records:
        candidates = [
            record.accession,
            record.study_accession,
            record.experiment_accession,
            record.sample_accession,
            record.submitter_accession,
            *record.run_accessions,
        ]
        if requested and requested in {item.upper() for item in candidates if item}:
            return record.study_accession or record.accession or requested_accession
    if records:
        return records[0].study_accession or records[0].accession or requested_accession
    return requested_accession


def filter_anchor_records(records: list[SraRecord], anchor: str) -> list[SraRecord]:
    wanted = anchor.strip().upper()
    filtered: list[SraRecord] = []
    for record in records:
        candidates = {
            item.upper()
            for item in [
                record.accession,
                record.study_accession,
                record.experiment_accession,
                record.sample_accession,
                record.submitter_accession,
                *record.run_accessions,
            ]
            if item
        }
        if wanted in candidates:
            filtered.append(record)
    return filtered


def infer_tissues(records: list[SraRecord]) -> list[str]:
    tissues: list[str] = []
    for record in records:
        title = (record.title or "").lower()
        summary = (record.summary or "").lower()
        haystack = f"{title} {summary}"
        if "abdominal fat" in haystack:
            tissues.append("abdominal fat")
        elif "adipose" in haystack:
            tissues.append("adipose")
        elif "cartilage" in haystack:
            tissues.append("cartilage")
        elif "islet" in haystack:
            tissues.append("islet")
    return unique_nonempty(tissues)


def infer_risks(records: list[SraRecord], organisms: list[str], strategies: list[str], layouts: list[str], platforms: list[str]) -> list[str]:
    risks: list[str] = []
    if len(organisms) > 1:
        risks.append("multiple_organisms")
    if len(strategies) > 1:
        risks.append("mixed_library_strategies")
    if len(layouts) > 1:
        risks.append("mixed_library_layouts")
    if len(platforms) > 1:
        risks.append("mixed_platforms")
    if not any(record.run_accessions for record in records):
        risks.append("run_accessions_unclear")
    if not any(record.organism for record in records):
        risks.append("organism_metadata_sparse")
    if not any(record.title for record in records):
        risks.append("title_metadata_sparse")
    return unique_nonempty(risks)


def assess_analysis_fit(
    analysis_goal: str,
    assay_family: str,
    sample_count: int | None,
    risks: list[str],
) -> tuple[list[str], list[str], list[str], str]:
    usable_for: list[str] = []
    not_recommended_for: list[str] = []
    rationale: list[str] = []
    if assay_family.lower() == "rna-seq":
        usable_for.append("raw_rnaseq_follow_up")
    else:
        rationale.append("assay_family_not_rnaseq")

    if analysis_goal == "differential_expression":
        if sample_count is not None and sample_count >= 6:
            usable_for.append("differential_expression")
        else:
            not_recommended_for.append("differential_expression")
            rationale.append("sample_support_small_for_deg")
    elif analysis_goal == "wgcna":
        if sample_count is not None and sample_count >= 15:
            usable_for.append("wgcna")
        else:
            not_recommended_for.append("wgcna")
            rationale.append("sample_support_small_for_wgcna")
    elif analysis_goal == "validation":
        usable_for.append("raw_data_validation_candidate")

    if "multiple_organisms" in risks or "mixed_library_strategies" in risks:
        rationale.append("cross_sample_heterogeneity_detected")
    overall = "medium"
    if any(flag in risks for flag in ("multiple_organisms", "mixed_library_strategies")):
        overall = "low"
    elif analysis_goal == "wgcna" and (sample_count is None or sample_count < 15):
        overall = "low"
    elif analysis_goal == "differential_expression" and sample_count is not None and sample_count >= 6:
        overall = "high"
    return unique_nonempty(usable_for), unique_nonempty(not_recommended_for), unique_nonempty(rationale), overall


@dataclass
class SraInspectionResult:
    accession: str
    source_type: str
    title: str
    organisms: list[str]
    tissues: list[str]
    assay_family: str
    library_layouts: list[str]
    platforms: list[str]
    sample_count: int | None
    run_count: int | None
    bioproject: str
    biosamples: list[str]
    linked_sra: list[str]
    linked_geo: list[str]
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
            "source_type": self.source_type,
            "title": self.title,
            "organisms": self.organisms,
            "tissues": self.tissues,
            "assay_family": self.assay_family,
            "library_layouts": self.library_layouts,
            "platforms": self.platforms,
            "sample_count": self.sample_count,
            "run_count": self.run_count,
            "bioproject": self.bioproject,
            "biosamples": self.biosamples,
            "linked_sra": self.linked_sra,
            "linked_geo": self.linked_geo,
            "risks": self.risks,
            "usable_for": self.usable_for,
            "not_recommended_for": self.not_recommended_for,
            "rationale": self.rationale,
            "provisional_local_assessment": self.overall_usability,
            "assessment_source": "local_aggregation_plus_sra_metadata",
            "data_availability": self.data_availability,
            "summary": self.summary,
            "normalized_evidence": build_sra_inspection_evidence(self, include_summary=True),
        }
        if include_raw:
            data["raw_notes"] = self.raw_notes
        return data


def inspect_accession(accession: str, analysis_goal: str, timeout: int, pause_seconds: float, email: str | None) -> SraInspectionResult:
    initial_records = fetch_records(accession, max_results=50, timeout=timeout, pause_seconds=pause_seconds, email=email)
    anchor = choose_anchor(initial_records, accession)
    if anchor and anchor.upper() != accession.strip().upper():
        anchor_records = fetch_records(anchor, max_results=100, timeout=timeout, pause_seconds=pause_seconds, email=email)
    else:
        anchor_records = initial_records
    records = filter_anchor_records(anchor_records or initial_records, anchor)
    if not records:
        records = initial_records

    title = next((record.title for record in records if record.title), "")
    source_type = next((record.source_type for record in records if record.source_type), "study")
    organisms = unique_nonempty([record.organism for record in records if record.organism])
    tissues = infer_tissues(records)
    strategies = unique_nonempty([record.library_strategy for record in records if record.library_strategy])
    layouts = unique_nonempty([record.library_layout for record in records if record.library_layout])
    platforms = unique_nonempty([record.platform for record in records if record.platform])
    biosamples = unique_nonempty([record.biosample for record in records if record.biosample])
    run_accessions = unique_nonempty([run for record in records for run in record.run_accessions])
    sample_count = len(biosamples) if biosamples else None
    run_count = len(run_accessions) if run_accessions else None
    bioproject = next((record.bioproject for record in records if record.bioproject), "")
    linked_sra = unique_nonempty(
        [
            record.accession
            for record in records
        ]
        + [record.study_accession for record in records]
        + [record.experiment_accession for record in records]
        + [record.sample_accession for record in records]
        + run_accessions
    )
    linked_geo = unique_nonempty([])
    risks = infer_risks(records, organisms, strategies, layouts, platforms)
    assay_family = strategies[0] if len(strategies) == 1 else "; ".join(strategies)
    usable_for, not_recommended_for, rationale, overall = assess_analysis_fit(
        analysis_goal=analysis_goal,
        assay_family=assay_family or "",
        sample_count=sample_count,
        risks=risks,
    )
    data_availability = {
        "sra_runs_available": bool(run_accessions),
        "raw_sequencing_available": bool(run_accessions),
        "library_layouts": layouts,
        "platforms": platforms,
        "needs_manual_check": not bool(title and organisms),
    }
    summary = next((record.summary for record in records if record.summary), "")
    raw_notes = {
        "requested_analysis_goal": analysis_goal,
        "anchor_accession": anchor,
        "record_count": len(records),
        "study_titles": unique_nonempty([record.title for record in records if record.title]),
        "submitter_accessions": unique_nonempty([record.submitter_accession for record in records if record.submitter_accession]),
    }
    return SraInspectionResult(
        accession=anchor or accession,
        source_type=source_type,
        title=title,
        organisms=organisms,
        tissues=tissues,
        assay_family=assay_family or "",
        library_layouts=layouts,
        platforms=platforms,
        sample_count=sample_count,
        run_count=run_count,
        bioproject=bioproject,
        biosamples=biosamples,
        linked_sra=linked_sra,
        linked_geo=linked_geo,
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
    parser = argparse.ArgumentParser(description="Inspect SRA accessions for BioFinder.")
    parser.add_argument("--accession", action="append", required=True, help="SRA accession to inspect. Repeatable.")
    parser.add_argument("--analysis-goal", default="general", help="Analysis goal such as DEG, WGCNA, or validation.")
    parser.add_argument("--timeout", type=int, default=30, help="Network timeout in seconds.")
    parser.add_argument("--pause-seconds", type=float, default=0.34, help="Pause between NCBI requests.")
    parser.add_argument("--email", help="Email to send to NCBI as recommended by E-utilities.")
    parser.add_argument("--json", action="store_true", help="Print JSON output.")
    parser.add_argument("--include-raw", action="store_true", help="Include raw aggregation notes.")
    return parser.parse_args()


def render_text(results: list[SraInspectionResult]) -> str:
    blocks: list[str] = []
    for result in results:
        lines = [
            f"{result.accession}",
            f"   title={result.title or 'n/a'}",
            f"   organisms={', '.join(result.organisms) or 'n/a'}",
            f"   tissues={', '.join(result.tissues) or 'n/a'}",
            f"   assay_family={result.assay_family or 'n/a'}",
            f"   sample_count={result.sample_count if result.sample_count is not None else 'n/a'}",
            f"   run_count={result.run_count if result.run_count is not None else 'n/a'}",
            f"   bioproject={result.bioproject or 'n/a'}",
            f"   provisional_local_assessment={result.overall_usability}",
            f"   usable_for={', '.join(result.usable_for) or 'n/a'}",
            f"   not_recommended_for={', '.join(result.not_recommended_for) or 'n/a'}",
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

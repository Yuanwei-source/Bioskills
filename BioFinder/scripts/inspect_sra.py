#!/usr/bin/env python
"""
Inspect one SRA accession and aggregate study-level evidence for BioFinder.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any

from evidence_schema import build_sra_inspection_evidence
from inspect_gse import normalize_goal
from search_geo import DEFAULT_TOOL, configure_stdio
from search_sra import SraRecord, fetch_records

RUNINFO_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"


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


def fetch_runinfo_rows(accession: str, timeout: int, email: str | None) -> list[dict[str, str]]:
    params = {
        "acc": accession,
    }
    if email:
        params["email"] = email
    url = f"{RUNINFO_URL}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": f"{DEFAULT_TOOL}/0.1 (+https://openai.com)",
            "Accept": "text/csv",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            text = response.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError:
        return []
    except urllib.error.URLError:
        return []
    reader = csv.DictReader(StringIO(text))
    rows: list[dict[str, str]] = []
    for row in reader:
        normalized = {str(key).strip(): (value or "").strip() for key, value in row.items() if key}
        if normalized:
            rows.append(normalized)
    return rows


def fetch_pysradb_rows(accession: str, timeout: int) -> list[dict[str, str]]:
    command = shutil.which("pysradb")
    if not command:
        return []
    try:
        completed = subprocess.run(
            [command, "metadata", accession, "--detailed"],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return []
    if completed.returncode != 0 or not completed.stdout.strip():
        return []
    text = completed.stdout.replace("\r\n", "\n").strip()
    try:
        dialect = csv.Sniffer().sniff(text[:4000], delimiters="\t,")
    except csv.Error:
        dialect = csv.excel_tab
    reader = csv.DictReader(StringIO(text), dialect=dialect)
    rows: list[dict[str, str]] = []
    for row in reader:
        normalized = {str(key).strip(): (value or "").strip() for key, value in row.items() if key}
        if normalized:
            rows.append(normalized)
    return rows


def parse_sample_attributes(text: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for chunk in text.split("||"):
        piece = chunk.strip().strip(";")
        if not piece:
            continue
        if ":" in piece:
            key, value = piece.split(":", 1)
        elif "=" in piece:
            key, value = piece.split("=", 1)
        else:
            continue
        clean_key = key.strip().lower().replace(" ", "_")
        clean_value = value.strip()
        if clean_key and clean_value:
            parsed[clean_key] = clean_value
    return parsed


def infer_design_signals(runinfo_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    candidate_keys = [
        "source_name",
        "sample_name",
        "sample_title",
        "experiment_title",
        "strain",
        "tissue",
        "cell_type",
        "treatment",
        "group",
        "condition",
        "phenotype",
        "genotype",
        "timepoint",
        "time_point",
        "sex",
        "age",
    ]
    counts_by_key: dict[str, dict[str, int]] = {}
    for row in runinfo_rows:
        lowered = {str(key).strip().lower(): (value or "").strip() for key, value in row.items() if key}
        if "source" in lowered and "source_name" not in lowered:
            lowered["source_name"] = lowered["source"]
        if "samplename" in lowered and "sample_name" not in lowered:
            lowered["sample_name"] = lowered["samplename"]
        if "sampletitle" in lowered and "sample_title" not in lowered:
            lowered["sample_title"] = lowered["sampletitle"]
        if "experiment" in lowered and "experiment_title" not in lowered:
            lowered["experiment_title"] = lowered["experiment"]
        merged: dict[str, str] = {}
        for key in candidate_keys:
            value = lowered.get(key, "").strip()
            if value:
                merged[key] = value
        sample_attrs = parse_sample_attributes(lowered.get("sample_attribute", ""))
        for key in candidate_keys:
            value = sample_attrs.get(key, "").strip()
            if value and key not in merged:
                merged[key] = value
        for key, value in merged.items():
            bucket = counts_by_key.setdefault(key, {})
            bucket[value] = bucket.get(value, 0) + 1
    results: list[dict[str, Any]] = []
    for key, counts in counts_by_key.items():
        nonempty = {name: count for name, count in counts.items() if name}
        if len(nonempty) < 2:
            continue
        results.append(
            {
                "column": key,
                "counts": dict(sorted(nonempty.items(), key=lambda item: (-item[1], item[0]))[:8]),
            }
        )
    priority = {
        "group": 0,
        "condition": 1,
        "treatment": 2,
        "phenotype": 3,
        "genotype": 4,
        "timepoint": 5,
        "time_point": 6,
        "source_name": 7,
        "sample_name": 8,
        "sample_title": 9,
        "tissue": 10,
        "cell_type": 11,
        "sex": 12,
        "age": 13,
        "strain": 14,
    }
    results.sort(key=lambda item: (priority.get(item["column"], 999), item["column"]))
    return results[:8]


def infer_design_clarity(group_signals: list[dict[str, Any]]) -> str:
    if not group_signals:
        return "low"
    important = {"group", "condition", "treatment", "phenotype", "genotype", "timepoint", "time_point"}
    if any(item["column"] in important for item in group_signals):
        return "high"
    return "partial"


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
    design_clarity: str
    group_signals: list[dict[str, Any]]
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
            "design_clarity": self.design_clarity,
            "group_signals": self.group_signals,
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
    runinfo_rows = fetch_runinfo_rows(anchor or accession, timeout=timeout, email=email)
    pysradb_rows = fetch_pysradb_rows(anchor or accession, timeout=timeout)
    metadata_rows = runinfo_rows if runinfo_rows else []
    if pysradb_rows:
        metadata_rows = metadata_rows + pysradb_rows
    group_signals = infer_design_signals(metadata_rows)
    design_clarity = infer_design_clarity(group_signals)
    sample_count = len(biosamples) if biosamples else None
    if sample_count is None and metadata_rows:
        sample_accessions = unique_nonempty(
            [
                row.get("BioSample", "")
                or row.get("BioSample_s", "")
                or row.get("biosample", "")
                or row.get("sample_accession", "")
                or row.get("sample_alias", "")
                for row in metadata_rows
            ]
        )
        if sample_accessions:
            sample_count = len(sample_accessions)
    run_count = len(run_accessions) if run_accessions else None
    if run_count is None and metadata_rows:
        run_names = unique_nonempty([row.get("Run", "") or row.get("run_accession", "") for row in metadata_rows])
        if run_names:
            run_count = len(run_names)
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
        "runinfo_rows": len(runinfo_rows),
        "pysradb_rows": len(pysradb_rows),
        "runinfo_design_clarity": design_clarity,
        "needs_manual_check": not bool(title and organisms),
    }
    summary = next((record.summary for record in records if record.summary), "")
    raw_notes = {
        "requested_analysis_goal": analysis_goal,
        "anchor_accession": anchor,
        "record_count": len(records),
        "study_titles": unique_nonempty([record.title for record in records if record.title]),
        "submitter_accessions": unique_nonempty([record.submitter_accession for record in records if record.submitter_accession]),
        "runinfo_columns": sorted(runinfo_rows[0].keys()) if runinfo_rows else [],
        "pysradb_columns": sorted(pysradb_rows[0].keys()) if pysradb_rows else [],
        "metadata_sources": [
            source
            for source, rows in [("runinfo", runinfo_rows), ("pysradb", pysradb_rows)]
            if rows
        ],
    }
    return SraInspectionResult(
        accession=anchor or accession,
        source_type=source_type,
        title=title,
        organisms=organisms,
        tissues=tissues,
        assay_family=assay_family or "",
        design_clarity=design_clarity,
        group_signals=group_signals,
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
            f"   design_clarity={result.design_clarity}",
            f"   sample_count={result.sample_count if result.sample_count is not None else 'n/a'}",
            f"   run_count={result.run_count if result.run_count is not None else 'n/a'}",
            f"   bioproject={result.bioproject or 'n/a'}",
            f"   provisional_local_assessment={result.overall_usability}",
            f"   usable_for={', '.join(result.usable_for) or 'n/a'}",
            f"   not_recommended_for={', '.join(result.not_recommended_for) or 'n/a'}",
            f"   risks={', '.join(result.risks) or 'n/a'}",
        ]
        if result.group_signals:
            preview = "; ".join(
                f"{item['column']}={list(item['counts'].items())[:4]}" for item in result.group_signals[:3]
            )
            lines.append(f"   group_signals={preview}")
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

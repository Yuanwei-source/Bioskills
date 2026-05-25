#!/usr/bin/env python
"""
Inspect one or more GEO Series accessions and extract structured evidence.

This script prefers GEOparse for accession-level parsing and uses search_geo.py
to supplement series-level metadata from NCBI esummary.

The local output is a provisional evidence summary plus lightweight heuristic
assessment. Final dataset judgment should be made by the LLM in the context of
the user's biological goal.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import GEOparse

from evidence_schema import build_geo_inspection_evidence
from search_geo import GeoRecord, fetch_records


DEFAULT_DESTDIR = Path("C:/tmp/geo_data_finder")


def configure_stdio() -> None:
    for stream_name in ("stdout", "stderr"):
        stream = getattr(sys, stream_name, None)
        if stream is not None and hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except Exception:  # noqa: BLE001
                pass


def coerce_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    return str(value).strip()


def normalize_goal(goal: str | None) -> str:
    text = coerce_str(goal).lower()
    if not text:
        return "general"
    if "wgcna" in text:
        return "wgcna"
    if "validation" in text:
        return "validation"
    if "machine" in text or "ml" in text:
        return "machine_learning"
    if "deg" in text or "differential" in text:
        return "differential_expression"
    return text.replace(" ", "_")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect GEO Series accessions for usability.")
    parser.add_argument("--accession", action="append", required=True, help="GSE accession. Repeatable.")
    parser.add_argument("--analysis-goal", default="general", help="Analysis goal, e.g. DEG, WGCNA, validation.")
    parser.add_argument("--destdir", default=str(DEFAULT_DESTDIR), help="Directory for GEOparse downloads/cache.")
    parser.add_argument("--how", choices=["quick", "full"], default="quick", help="GEOparse download mode.")
    parser.add_argument("--email", help="Email for NCBI/GEOparse operations.")
    parser.add_argument("--json", action="store_true", help="Print JSON output.")
    parser.add_argument("--include-raw", action="store_true", help="Include raw metadata snippets.")
    return parser.parse_args()


def get_summary_record(accession: str, email: str | None) -> GeoRecord | None:
    records = fetch_records(
        query=accession,
        max_results=5,
        timeout=30,
        pause_seconds=0.34,
        series_only=True,
        email=email,
    )
    for record in records:
        if record.accession.casefold() == accession.casefold():
            return record
    return None


def safe_get_geo(accession: str, destdir: Path, how: str) -> tuple[Any | None, str | None]:
    try:
        gse = GEOparse.get_GEO(geo=accession, destdir=str(destdir), how=how, silent=True)
        return gse, None
    except Exception as exc:  # noqa: BLE001
        return None, str(exc)


def infer_assay_family(summary_record: GeoRecord | None, phenotype_data: Any | None) -> str:
    text_parts = []
    if summary_record is not None:
        text_parts.extend([summary_record.gds_type, summary_record.title, summary_record.summary])
    if phenotype_data is not None and not phenotype_data.empty:
        for column in phenotype_data.columns:
            if column.startswith("molecule_") or column.startswith("type"):
                text_parts.extend(phenotype_data[column].dropna().astype(str).tolist())
    haystack = " ".join(coerce_str(part).lower() for part in text_parts if coerce_str(part))

    if "single cell" in haystack or "single-cell" in haystack or "single nucleus" in haystack:
        return "single_cell_rna_seq"
    if "non-coding rna" in haystack or "small rna" in haystack or "mirna" in haystack:
        return "non_coding_rna_seq"
    if "high throughput sequencing" in haystack or "rna-seq" in haystack or "rnaseq" in haystack:
        return "bulk_rna_seq"
    if "expression profiling by array" in haystack or "microarray" in haystack or "array" in haystack:
        return "microarray"
    return "unknown"


def find_columns(phenotype_data: Any, patterns: list[str]) -> list[str]:
    if phenotype_data is None or phenotype_data.empty:
        return []
    matches: list[str] = []
    for column in phenotype_data.columns:
        lower = column.casefold()
        if any(pattern in lower for pattern in patterns):
            matches.append(column)
    return matches


def unique_nonempty(values: list[str]) -> list[str]:
    seen: list[str] = []
    for value in values:
        cleaned = coerce_str(value)
        if cleaned and cleaned not in seen:
            seen.append(cleaned)
    return seen


def collect_values(phenotype_data: Any, columns: list[str]) -> list[str]:
    values: list[str] = []
    if phenotype_data is None or phenotype_data.empty:
        return values
    for column in columns:
        if column in phenotype_data.columns:
            series = phenotype_data[column].dropna().astype(str).tolist()
            values.extend(item.strip() for item in series if item and item.strip())
    return unique_nonempty(values)


def summarize_tissues(phenotype_data: Any) -> list[str]:
    tissue_columns = find_columns(phenotype_data, ["characteristics_ch1.", "tissue"])
    tissue_columns = [column for column in tissue_columns if "tissue" in column.casefold()]
    if tissue_columns:
        return collect_values(phenotype_data, tissue_columns)
    return collect_values(phenotype_data, find_columns(phenotype_data, ["source_name"]))


def summarize_organisms(phenotype_data: Any) -> list[str]:
    return collect_values(phenotype_data, find_columns(phenotype_data, ["organism_ch", "organism"]))


def parse_supplementary_formats(value: str) -> list[str]:
    formats: list[str] = []
    for part in coerce_str(value).replace(",", ";").split(";"):
        cleaned = part.strip().upper()
        if cleaned:
            formats.append(cleaned)
    return unique_nonempty(formats)


def collect_metadata_text(gse: Any) -> str:
    if gse is None or not getattr(gse, "gsms", None):
        return ""
    parts: list[str] = []
    for gsm in gse.gsms.values():
        for values in getattr(gsm, "metadata", {}).values():
            if isinstance(values, list):
                parts.extend(coerce_str(item) for item in values if coerce_str(item))
            else:
                text = coerce_str(values)
                if text:
                    parts.append(text)
    return " ".join(parts)


def has_pattern(text: str, pattern: str) -> bool:
    return re.search(pattern, text, flags=re.IGNORECASE) is not None


def infer_data_availability(summary_record: GeoRecord | None, gse: Any, assay_family: str) -> dict[str, Any]:
    supplementary_formats = parse_supplementary_formats(summary_record.supplemental_files if summary_record is not None else "")
    metadata_text = collect_metadata_text(gse).casefold()
    summary_text = " ".join(
        coerce_str(item)
        for item in [
            summary_record.summary if summary_record is not None else "",
            summary_record.title if summary_record is not None else "",
        ]
        if coerce_str(item)
    ).casefold()
    combined_text = " ".join(part for part in [metadata_text, summary_text] if part)

    processed_matrix_clues = any(fmt in supplementary_formats for fmt in ["TXT", "TSV", "CSV", "TAB", "XLSX", "MTX", "H5"])
    microarray_raw_clues = any(fmt in supplementary_formats for fmt in ["CEL"])
    sequencing_raw_clues = any(fmt in supplementary_formats for fmt in ["FASTQ", "SRA"])
    if not sequencing_raw_clues and has_pattern(combined_text, r"\b(?:sra|fastq|srx\d+|srr\d+)\b"):
        sequencing_raw_clues = True

    likely_raw_counts = False
    if assay_family in {"bulk_rna_seq", "single_cell_rna_seq"}:
        if has_pattern(combined_text, r"\b(?:raw counts?|counts? matrix|gene counts?)\b"):
            likely_raw_counts = True
        elif "MTX" in supplementary_formats:
            likely_raw_counts = True

    normalized_expression_clues = has_pattern(
        combined_text,
        r"\b(?:normalized|fpkm|tpm|rpkm|log2|expression matrix)\b",
    )

    matrix_guess = "unknown"
    if assay_family == "microarray":
        if microarray_raw_clues:
            matrix_guess = "microarray_raw_plus_processed_possible"
        elif processed_matrix_clues:
            matrix_guess = "processed_microarray_matrix_likely"
    elif assay_family in {"bulk_rna_seq", "single_cell_rna_seq"}:
        if likely_raw_counts:
            matrix_guess = "raw_counts_likely"
        elif processed_matrix_clues:
            matrix_guess = "processed_matrix_likely"
    elif assay_family == "non_coding_rna_seq" and processed_matrix_clues:
        matrix_guess = "processed_non_coding_matrix_likely"

    return {
        "supplementary_formats": supplementary_formats,
        "processed_matrix_clues": processed_matrix_clues,
        "microarray_raw_clues": microarray_raw_clues,
        "sequencing_raw_clues": sequencing_raw_clues,
        "likely_raw_counts": likely_raw_counts,
        "normalized_expression_clues": normalized_expression_clues,
        "matrix_guess": matrix_guess,
    }


def classify_group_columns(phenotype_data: Any) -> list[dict[str, Any]]:
    if phenotype_data is None or phenotype_data.empty:
        return []

    blocked_prefixes = (
        "geo_accession",
        "status",
        "submission_date",
        "last_update_date",
        "contact_",
        "supplementary_file",
        "extract_protocol",
        "label_protocol",
        "scan_protocol",
        "hyb_protocol",
        "description",
        "data_processing",
        "library_selection",
        "library_source",
        "library_strategy",
        "relation",
        "data_row_count",
        "series_id",
        "platform_id",
        "platform_taxid",
        "sample_taxid",
    )
    candidates: list[dict[str, Any]] = []
    for column in phenotype_data.columns:
        lower = column.casefold()
        if lower.startswith(blocked_prefixes):
            continue
        values = [item.strip() for item in phenotype_data[column].dropna().astype(str).tolist() if item.strip()]
        if not values:
            continue
        unique_values = unique_nonempty(values)
        if len(unique_values) <= 1 or len(unique_values) > min(12, max(3, len(values) - 1)):
            continue
        counts = Counter(values)
        candidates.append(
            {
                "column": column,
                "unique_values": unique_values,
                "counts": dict(counts),
            }
        )
    candidates.sort(key=lambda item: (len(item["unique_values"]), item["column"]))
    return candidates


def infer_design_clarity(group_columns: list[dict[str, Any]]) -> str:
    if not group_columns:
        return "low"
    for group in group_columns:
        name = group["column"].casefold()
        if any(token in name for token in ["treatment", "condition", "group", "phenotype", "disease", "genotype"]):
            return "high"
    return "partial"


def detect_series_complexity(summary_record: GeoRecord | None, tissues: list[str], organisms: list[str], group_columns: list[dict[str, Any]]) -> list[str]:
    risks: list[str] = []
    if summary_record is not None:
        title_summary = " ".join([summary_record.title, summary_record.summary]).casefold()
        if "across species" in title_summary or "bodymap" in title_summary:
            risks.append("multi_species_series")
        if "multiomic" in title_summary or "multi-omic" in title_summary:
            risks.append("multiomic_series")
        if "jointly profiled" in title_summary:
            risks.append("mixed_assay_design")
        if "superseries" in title_summary or "superseries" in title_summary:
            risks.append("superseries")
    if len(organisms) > 1:
        risks.append("multiple_organisms")
    if len(tissues) > 2:
        risks.append("multiple_tissues")
    if len(group_columns) > 8:
        risks.append("many_grouping_dimensions")
    return unique_nonempty(risks)


def assess_analysis_fit(
    analysis_goal: str,
    assay_family: str,
    sample_count: int | None,
    design_clarity: str,
    risks: list[str],
) -> tuple[list[str], list[str], list[str]]:
    usable_for: list[str] = []
    not_recommended_for: list[str] = []
    rationale: list[str] = []

    if sample_count is None:
        rationale.append("sample_count_unverified")
    elif sample_count < 3:
        rationale.append("extremely_small_sample_count")
    elif sample_count < 6:
        rationale.append("small_sample_count")

    if assay_family == "bulk_rna_seq":
        usable_for.append("differential_expression")
    elif assay_family == "microarray":
        usable_for.extend(["differential_expression", "validation"])
    elif assay_family == "single_cell_rna_seq":
        rationale.append("single_cell_data")
        not_recommended_for.extend(["bulk_differential_expression", "wgcna"])
    elif assay_family == "non_coding_rna_seq":
        rationale.append("non_coding_assay_family")
        not_recommended_for.append("bulk_mrna_workflows")

    if design_clarity == "high":
        rationale.append("clear_group_structure")
    elif design_clarity == "partial":
        rationale.append("partially_clear_group_structure")
    else:
        rationale.append("unclear_group_structure")
        not_recommended_for.extend(["differential_expression", "validation"])

    if "multiple_tissues" in risks:
        rationale.append("mixed_tissue_risk")
        not_recommended_for.append("simple_case_control")
        not_recommended_for.append("whole_series_wgcna")
    if "multiple_organisms" in risks or "multi_species_series" in risks:
        rationale.append("multi_species_risk")
        not_recommended_for.extend(["validation", "simple_case_control"])
    if "superseries" in risks:
        rationale.append("superseries_requires_subseries_check")
        not_recommended_for.append("whole_series_interpretation")

    if sample_count is not None and sample_count >= 15 and design_clarity in {"high", "partial"} and "multiple_tissues" not in risks:
        usable_for.append("wgcna")
    else:
        not_recommended_for.append("wgcna")

    if sample_count is not None and sample_count >= 6 and design_clarity != "low":
        usable_for.append("validation")
    elif analysis_goal == "validation":
        rationale.append("weak_validation_candidate")

    if sample_count is not None and sample_count >= 30 and design_clarity == "high":
        usable_for.append("machine_learning_exploration")
    else:
        not_recommended_for.append("machine_learning")

    if analysis_goal == "differential_expression" and "differential_expression" not in usable_for:
        rationale.append("requested_deg_but_design_is_weak")
    if analysis_goal == "wgcna" and "wgcna" not in usable_for:
        rationale.append("requested_wgcna_but_sample_size_or_design_is_weak")

    return unique_nonempty(usable_for), unique_nonempty(not_recommended_for), unique_nonempty(rationale)


def overall_usability(
    usable_for: list[str],
    not_recommended_for: list[str],
    risks: list[str],
    design_clarity: str,
    analysis_goal: str,
) -> str:
    if not usable_for:
        return "low"
    goal_map = {
        "differential_expression": "differential_expression",
        "wgcna": "wgcna",
        "validation": "validation",
        "machine_learning": "machine_learning_exploration",
    }
    mapped_goal = goal_map.get(analysis_goal)
    if mapped_goal and mapped_goal not in usable_for:
        return "medium" if design_clarity != "low" else "low"
    if design_clarity == "high" and len(risks) <= 1:
        return "high"
    if len(not_recommended_for) >= 3 or "multiple_organisms" in risks:
        return "low"
    return "medium"


@dataclass
class InspectionResult:
    accession: str
    title: str
    taxon: str
    assay_family: str
    gds_type: str
    platform: str
    sample_count: int | None
    sample_count_status: str
    design_clarity: str
    tissues: list[str]
    organisms: list[str]
    group_columns: list[dict[str, Any]]
    usable_for: list[str]
    not_recommended_for: list[str]
    risks: list[str]
    rationale: list[str]
    overall_usability: str
    data_availability: dict[str, Any]
    supplemental_files: str
    ftp_link: str
    bioproject: str
    summary: str
    raw_notes: dict[str, Any]

    def to_dict(self, include_raw: bool) -> dict[str, Any]:
        data = {
            "accession": self.accession,
            "title": self.title,
            "taxon": self.taxon,
            "assay_family": self.assay_family,
            "gds_type": self.gds_type,
            "platform": self.platform,
            "sample_count": self.sample_count,
            "sample_count_status": self.sample_count_status,
            "design_clarity": self.design_clarity,
            "tissues": self.tissues,
            "organisms": self.organisms,
            "group_columns": self.group_columns,
            "usable_for": self.usable_for,
            "not_recommended_for": self.not_recommended_for,
            "risks": self.risks,
            "rationale": self.rationale,
            "provisional_local_assessment": self.overall_usability,
            "assessment_source": "local_heuristics_plus_extracted_evidence",
            "data_availability": self.data_availability,
            "supplemental_files": self.supplemental_files,
            "ftp_link": self.ftp_link,
            "bioproject": self.bioproject,
            "summary": self.summary,
            "normalized_evidence": build_geo_inspection_evidence(self, include_summary=True),
        }
        if include_raw:
            data["raw_notes"] = self.raw_notes
        return data


def inspect_accession(accession: str, analysis_goal: str, destdir: Path, how: str, email: str | None) -> InspectionResult:
    summary_record = get_summary_record(accession, email=email)
    gse, geoparse_error = safe_get_geo(accession=accession, destdir=destdir, how=how)

    phenotype_data = None
    gse_name = accession
    gse_metadata_keys: list[str] = []
    if gse is not None:
        gse_name = coerce_str(getattr(gse, "name", "")) or accession
        gse_metadata_keys = sorted(list(getattr(gse, "metadata", {}).keys()))
        try:
            phenotype_data = gse.phenotype_data
        except Exception:  # noqa: BLE001
            phenotype_data = None

    sample_count: int | None = None
    sample_count_status = "unavailable"
    if gse is not None and getattr(gse, "gsms", None):
        sample_count = len(gse.gsms)
        sample_count_status = "verified_from_gsms"
    elif summary_record is not None and summary_record.samples is not None:
        sample_count = summary_record.samples
        sample_count_status = "series_summary_only"

    organisms = summarize_organisms(phenotype_data)
    if not organisms and summary_record is not None and summary_record.taxon:
        organisms = [summary_record.taxon]
    tissues = summarize_tissues(phenotype_data)
    group_columns = classify_group_columns(phenotype_data)
    design_clarity = infer_design_clarity(group_columns)
    assay_family = infer_assay_family(summary_record, phenotype_data)
    risks = detect_series_complexity(summary_record, tissues, organisms, group_columns)
    usable_for, not_recommended_for, rationale = assess_analysis_fit(
        analysis_goal=analysis_goal,
        assay_family=assay_family,
        sample_count=sample_count,
        design_clarity=design_clarity,
        risks=risks,
    )
    data_availability = infer_data_availability(summary_record, gse, assay_family)
    if geoparse_error:
        rationale.append("geoparse_fetch_failed")
        risks.append("partial_metadata_only")

    overall = overall_usability(usable_for, not_recommended_for, risks, design_clarity, analysis_goal)

    title = summary_record.title if summary_record is not None else gse_name
    taxon = summary_record.taxon if summary_record is not None else ", ".join(organisms)
    gds_type = summary_record.gds_type if summary_record is not None else ""
    platform = summary_record.platform if summary_record is not None else "; ".join(sorted(getattr(gse, "gpls", {}).keys()))
    supplemental_files = summary_record.supplemental_files if summary_record is not None else ""
    ftp_link = summary_record.ftp_link if summary_record is not None else ""
    bioproject = summary_record.bioproject if summary_record is not None else ""
    summary = summary_record.summary if summary_record is not None else ""

    raw_notes = {
        "requested_analysis_goal": analysis_goal,
        "geoparse_mode": how,
        "geoparse_error": geoparse_error,
        "gse_name": gse_name,
        "gse_metadata_keys": gse_metadata_keys,
        "phenotype_shape": list(phenotype_data.shape) if phenotype_data is not None else None,
        "top_group_columns": group_columns[:5],
    }

    return InspectionResult(
        accession=accession,
        title=title,
        taxon=taxon,
        assay_family=assay_family,
        gds_type=gds_type,
        platform=platform,
        sample_count=sample_count,
        sample_count_status=sample_count_status,
        design_clarity=design_clarity,
        tissues=tissues,
        organisms=organisms,
        group_columns=group_columns[:8],
        usable_for=usable_for,
        not_recommended_for=not_recommended_for,
        risks=unique_nonempty(risks),
        rationale=rationale,
        overall_usability=overall,
        data_availability=data_availability,
        supplemental_files=supplemental_files,
        ftp_link=ftp_link,
        bioproject=bioproject,
        summary=summary,
        raw_notes=raw_notes,
    )


def render_text(results: list[InspectionResult]) -> str:
    blocks: list[str] = []
    for result in results:
        lines = [
            f"{result.accession}",
            f"   title={result.title or 'n/a'}",
            f"   taxon={result.taxon or 'n/a'}",
            f"   assay_family={result.assay_family}",
            f"   gds_type={result.gds_type or 'n/a'}",
            f"   platform={result.platform or 'n/a'}",
            f"   sample_count={result.sample_count if result.sample_count is not None else 'n/a'} ({result.sample_count_status})",
            f"   design_clarity={result.design_clarity}",
            f"   provisional_local_assessment={result.overall_usability}",
            f"   usable_for={', '.join(result.usable_for) or 'n/a'}",
            f"   not_recommended_for={', '.join(result.not_recommended_for) or 'n/a'}",
            f"   data_availability={result.data_availability.get('matrix_guess', 'unknown')}",
            f"   tissues={', '.join(result.tissues) or 'n/a'}",
            f"   organisms={', '.join(result.organisms) or 'n/a'}",
            f"   risks={', '.join(result.risks) or 'n/a'}",
            f"   rationale={', '.join(result.rationale) or 'n/a'}",
        ]
        if result.group_columns:
            preview = "; ".join(
                f"{item['column']}={list(item['counts'].items())[:4]}" for item in result.group_columns[:3]
            )
            lines.append(f"   group_columns={preview}")
        blocks.append("\n".join(lines))
    return "\n\n".join(blocks)


def main() -> int:
    configure_stdio()
    args = parse_args()
    analysis_goal = normalize_goal(args.analysis_goal)
    destdir = Path(args.destdir)
    destdir.mkdir(parents=True, exist_ok=True)

    results = [
        inspect_accession(
            accession=coerce_str(accession),
            analysis_goal=analysis_goal,
            destdir=destdir,
            how=args.how,
            email=args.email,
        )
        for accession in args.accession
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

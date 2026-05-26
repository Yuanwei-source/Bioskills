#!/usr/bin/env python
"""
Helpers for normalizing provider-specific findings into a common BioFinder
evidence shape.
"""

from __future__ import annotations

import re
from typing import Any


def _clean_list(values: list[str] | None) -> list[str]:
    if not values:
        return []
    seen: set[str] = set()
    cleaned: list[str] = []
    for value in values:
        text = str(value).strip()
        if not text or text in seen:
            continue
        seen.add(text)
        cleaned.append(text)
    return cleaned


def _split_accessions(value: str | None, prefixes: tuple[str, ...]) -> list[str]:
    text = (value or "").strip()
    if not text:
        return []
    parts = [item.strip() for item in text.replace(",", ";").split(";")]
    return [item for item in parts if item and item.upper().startswith(prefixes)]


def _extract_accessions_from_text(value: str | None, pattern: str) -> list[str]:
    text = (value or "").strip()
    if not text:
        return []
    return _clean_list(re.findall(pattern, text, flags=re.IGNORECASE))


def _geo_link_context(record: Any) -> str:
    parts = [
        getattr(record, "bioproject", ""),
        getattr(record, "supplemental_files", ""),
        getattr(record, "summary", ""),
        getattr(record, "title", ""),
    ]
    return " ".join(part for part in parts if part)


def build_geo_candidate_evidence(record: Any, include_summary: bool = False) -> dict[str, Any]:
    link_context = _geo_link_context(record)
    source_links = {
        "geo": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={record.accession}" if record.accession else "",
        "ftp": record.ftp_link,
    }
    return {
        "source": "GEO",
        "accession": record.accession,
        "source_type": "series" if str(record.accession).upper().startswith("GSE") else "record",
        "title": record.title,
        "organisms": _clean_list([record.taxon] if record.taxon else []),
        "tissues": [],
        "assay_family": "",
        "sample_count": record.samples,
        "run_count": None,
        "group_signals": [],
        "data_availability": {
            "supplementary_files": record.supplemental_files,
        },
        "linked_accessions": {
            "geo": _clean_list([record.accession] if record.accession else []),
            "sra": _extract_accessions_from_text(link_context, r"\bSR[APXRS]\d+\b"),
            "bioproject": _clean_list(
                _split_accessions(record.bioproject, ("PRJ",))
                + _extract_accessions_from_text(link_context, r"\bPRJ[EDN][A-Z]?\d+\b")
            ),
            "biosample": _extract_accessions_from_text(link_context, r"\bSAM[END][A-Z]?\d+\b"),
        },
        "risks": [],
        "source_links": source_links,
        "local_notes": {
            "gds_type": record.gds_type,
            "platform": record.platform,
            "matched_queries": list(record.matched_queries),
            "summary": record.summary if include_summary else "",
        },
    }


def build_geo_inspection_evidence(result: Any, include_summary: bool = True) -> dict[str, Any]:
    link_context = _geo_link_context(result)
    source_links = {
        "geo": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={result.accession}" if result.accession else "",
        "ftp": result.ftp_link,
    }
    linked_geo = [result.accession] if result.accession else []
    return {
        "source": "GEO",
        "accession": result.accession,
        "source_type": "series" if str(result.accession).upper().startswith("GSE") else "record",
        "title": result.title,
        "organisms": _clean_list(list(result.organisms) or ([result.taxon] if result.taxon else [])),
        "tissues": _clean_list(list(result.tissues)),
        "assay_family": result.assay_family,
        "sample_count": result.sample_count,
        "run_count": None,
        "group_signals": list(result.group_columns),
        "data_availability": dict(result.data_availability),
        "linked_accessions": {
            "geo": _clean_list(linked_geo),
            "sra": _extract_accessions_from_text(link_context, r"\bSR[APXRS]\d+\b"),
            "bioproject": _clean_list(
                _split_accessions(result.bioproject, ("PRJ",))
                + _extract_accessions_from_text(link_context, r"\bPRJ[EDN][A-Z]?\d+\b")
            ),
            "biosample": _extract_accessions_from_text(link_context, r"\bSAM[END][A-Z]?\d+\b"),
        },
        "risks": _clean_list(list(result.risks)),
        "source_links": source_links,
        "local_notes": {
            "taxon": result.taxon,
            "gds_type": result.gds_type,
            "platform": result.platform,
            "design_clarity": result.design_clarity,
            "provisional_local_assessment": result.overall_usability,
            "usable_for": list(result.usable_for),
            "not_recommended_for": list(result.not_recommended_for),
            "rationale": list(result.rationale),
            "summary": result.summary if include_summary else "",
            "supplemental_files": result.supplemental_files,
        },
    }


def build_sra_candidate_evidence(record: Any, include_summary: bool = False) -> dict[str, Any]:
    source_links = {
        "sra": f"https://www.ncbi.nlm.nih.gov/sra/{record.accession}" if record.accession else "",
        "bioproject": f"https://www.ncbi.nlm.nih.gov/bioproject/{record.bioproject}" if record.bioproject else "",
        "biosample": f"https://www.ncbi.nlm.nih.gov/biosample/{record.biosample}" if record.biosample else "",
    }
    return {
        "source": "SRA",
        "accession": record.accession,
        "source_type": record.source_type or "study",
        "title": record.title,
        "organisms": _clean_list([record.organism] if record.organism else []),
        "tissues": [],
        "assay_family": record.library_strategy,
        "sample_count": record.sample_count,
        "run_count": record.run_count,
        "group_signals": [],
        "data_availability": {
            "sra_runs_available": bool(record.run_count),
            "library_layout": record.library_layout,
            "platform": record.platform,
        },
        "linked_accessions": {
            "geo": _extract_accessions_from_text(record.title, r"\bGS[EM]\d+\b"),
            "sra": _clean_list(
                [record.accession, getattr(record, "study_accession", ""), getattr(record, "experiment_accession", ""), getattr(record, "sample_accession", ""), getattr(record, "submitter_accession", "")]
                + list(getattr(record, "run_accessions", []))
            ),
            "bioproject": _clean_list([record.bioproject] if record.bioproject else []),
            "biosample": _clean_list([record.biosample] if record.biosample else []),
        },
        "risks": [],
        "source_links": source_links,
        "local_notes": {
            "platform": record.platform,
            "study_accession": getattr(record, "study_accession", ""),
            "experiment_accession": getattr(record, "experiment_accession", ""),
            "sample_accession": getattr(record, "sample_accession", ""),
            "submitter_accession": getattr(record, "submitter_accession", ""),
            "matched_queries": list(record.matched_queries),
            "summary": record.summary if include_summary else "",
        },
    }


def build_sra_inspection_evidence(result: Any, include_summary: bool = True) -> dict[str, Any]:
    source_links = {
        "sra": f"https://www.ncbi.nlm.nih.gov/sra/{result.accession}" if result.accession else "",
        "bioproject": f"https://www.ncbi.nlm.nih.gov/bioproject/{result.bioproject}" if result.bioproject else "",
    }
    return {
        "source": "SRA",
        "accession": result.accession,
        "source_type": result.source_type,
        "title": result.title,
        "organisms": _clean_list(list(result.organisms)),
        "tissues": _clean_list(list(result.tissues)),
        "assay_family": result.assay_family,
        "sample_count": result.sample_count,
        "run_count": result.run_count,
        "group_signals": list(result.group_signals),
        "data_availability": dict(result.data_availability),
        "linked_accessions": {
            "geo": _clean_list(list(result.linked_geo)),
            "sra": _clean_list(list(result.linked_sra)),
            "bioproject": _clean_list([result.bioproject] if result.bioproject else []),
            "biosample": _clean_list(list(result.biosamples)),
        },
        "risks": _clean_list(list(result.risks)),
        "source_links": source_links,
        "local_notes": {
            "design_clarity": result.design_clarity,
            "library_layouts": list(result.library_layouts),
            "platforms": list(result.platforms),
            "provisional_local_assessment": result.overall_usability,
            "usable_for": list(result.usable_for),
            "not_recommended_for": list(result.not_recommended_for),
            "rationale": list(result.rationale),
            "summary": result.summary if include_summary else "",
        },
    }


def build_bioproject_candidate_evidence(record: Any, include_summary: bool = False) -> dict[str, Any]:
    source_links = {
        "bioproject": f"https://www.ncbi.nlm.nih.gov/bioproject/{record.accession}" if record.accession else "",
    }
    return {
        "source": "BioProject",
        "accession": record.accession,
        "source_type": "project",
        "title": record.title,
        "organisms": _clean_list([record.organism] if record.organism else []),
        "tissues": [],
        "assay_family": "",
        "sample_count": None,
        "run_count": None,
        "group_signals": [],
        "data_availability": {
            "linked_sra_count": len(getattr(record, "linked_sra", [])),
            "linked_biosample_count": len(getattr(record, "linked_biosamples", [])),
        },
        "linked_accessions": {
            "geo": _extract_accessions_from_text(record.title, r"\bGS[EM]\d+\b"),
            "sra": _clean_list(list(getattr(record, "linked_sra", []))),
            "bioproject": _clean_list([record.accession] if record.accession else []),
            "biosample": _clean_list(list(getattr(record, "linked_biosamples", []))),
        },
        "risks": [],
        "source_links": source_links,
        "local_notes": {
            "matched_queries": list(getattr(record, "matched_queries", [])),
            "summary": record.summary if include_summary else "",
        },
    }


def build_bioproject_inspection_evidence(result: Any, include_summary: bool = True) -> dict[str, Any]:
    source_links = {
        "bioproject": f"https://www.ncbi.nlm.nih.gov/bioproject/{result.accession}" if result.accession else "",
    }
    return {
        "source": "BioProject",
        "accession": result.accession,
        "source_type": "project",
        "title": result.title,
        "organisms": _clean_list(list(result.organisms)),
        "tissues": _clean_list(list(result.tissues)),
        "assay_family": "",
        "sample_count": None,
        "run_count": result.run_count,
        "group_signals": [],
        "data_availability": dict(result.data_availability),
        "linked_accessions": {
            "geo": _clean_list(list(result.linked_geo)),
            "sra": _clean_list(list(result.linked_sra)),
            "bioproject": _clean_list([result.accession] if result.accession else []),
            "biosample": _clean_list(list(result.linked_biosamples)),
        },
        "risks": _clean_list(list(result.risks)),
        "source_links": source_links,
        "local_notes": {
            "provisional_local_assessment": result.overall_usability,
            "usable_for": list(result.usable_for),
            "not_recommended_for": list(result.not_recommended_for),
            "rationale": list(result.rationale),
            "summary": result.summary if include_summary else "",
        },
    }

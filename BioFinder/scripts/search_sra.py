#!/usr/bin/env python
"""
Search NCBI SRA through E-utilities and return candidate study metadata.

This is the first v2 provider expansion for BioFinder. It stays intentionally
lightweight: search and evidence extraction only, no download workflows.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import textwrap
import time
import warnings
from dataclasses import dataclass, field
from typing import Any

from evidence_schema import build_sra_candidate_evidence
from search_geo import (
    DEFAULT_TOOL,
    Entrez,
    biopython_available,
    chunked,
    coerce_str,
    configure_entrez,
    configure_stdio,
    plainify,
    request_json,
)


def contains_casefold(text: str, needle: str) -> bool:
    return needle.casefold() in text.casefold()


def flatten_mixed_field(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, list):
        parts = [flatten_mixed_field(item) for item in value]
        return "; ".join(part for part in parts if part)
    if isinstance(value, dict):
        for key in ("accession", "name", "title", "value", "label"):
            text = coerce_str(value.get(key))
            if text:
                return text
        return "; ".join(f"{key}={flatten_mixed_field(subvalue)}" for key, subvalue in value.items() if flatten_mixed_field(subvalue))
    return coerce_str(value)


def parse_expxml_field(text: str, pattern: str) -> str:
    if not text:
        return ""
    match = re.search(pattern, text, flags=re.IGNORECASE)
    return coerce_str(match.group(1)) if match else ""


def count_runs(text: str) -> int | None:
    if not text:
        return None
    matches = re.findall(r"\bSRR\d+\b", text, flags=re.IGNORECASE)
    if matches:
        return len(set(match.upper() for match in matches))
    return None


@dataclass
class SraRecord:
    uid: str
    accession: str = ""
    source_type: str = "study"
    study_accession: str = ""
    experiment_accession: str = ""
    sample_accession: str = ""
    run_accessions: list[str] = field(default_factory=list)
    submitter_accession: str = ""
    title: str = ""
    summary: str = ""
    organism: str = ""
    library_strategy: str = ""
    library_layout: str = ""
    platform: str = ""
    sample_count: int | None = None
    run_count: int | None = None
    bioproject: str = ""
    biosample: str = ""
    raw: dict[str, Any] = field(default_factory=dict)
    matched_queries: list[str] = field(default_factory=list)

    @classmethod
    def from_esummary(cls, payload: dict[str, Any]) -> "SraRecord":
        if not isinstance(payload, dict):
            payload = {"title": flatten_mixed_field(payload)}
        expxml = flatten_mixed_field(payload.get("ExpXml") or payload.get("expxml"))
        runs_field = flatten_mixed_field(payload.get("Runs") or payload.get("runs"))
        summary_title = parse_expxml_field(expxml, r"<Title>([^<]+)</Title>")
        summary = coerce_str(payload.get("summary") or payload.get("Summary"))
        title = coerce_str(payload.get("title") or payload.get("Title")) or summary_title
        organism = coerce_str(payload.get("organism") or payload.get("Organism"))
        if not organism:
            organism = parse_expxml_field(expxml, r'ScientificName="([^"]+)"')
        library_strategy = (
            coerce_str(payload.get("library_strategy") or payload.get("LibraryStrategy"))
            or parse_expxml_field(expxml, r'LIBRARY_STRATEGY>([^<]+)<')
        )
        library_layout = (
            coerce_str(payload.get("library_layout") or payload.get("LibraryLayout"))
            or parse_expxml_field(expxml, r'LIBRARY_LAYOUT>\s*<([^ >/]+)')
        )
        platform = (
            coerce_str(payload.get("platform") or payload.get("Platform"))
            or parse_expxml_field(expxml, r'Platform[^>]*instrument_model="([^"]+)"')
        )
        bioproject = coerce_str(payload.get("bioproject") or payload.get("Bioproject"))
        if not bioproject:
            bioproject = parse_expxml_field(expxml, r"\b(PRJ[EDN][A-Z]?\d+)\b")
        biosample = coerce_str(payload.get("biosample") or payload.get("Biosample"))
        if not biosample:
            biosample = parse_expxml_field(expxml, r"\b(SAM[END][A-Z]?\d+)\b")
        submitter_accession = coerce_str(payload.get("accession") or payload.get("Accession")) or parse_expxml_field(
            expxml, r'<Submitter acc="([^"]+)"'
        )
        study_accession = parse_expxml_field(expxml, r'<Study acc="([^"]+)"')
        experiment_accession = parse_expxml_field(expxml, r'<Experiment acc="([^"]+)"')
        sample_accession = parse_expxml_field(expxml, r'<Sample acc="([^"]+)"')
        run_accessions = re.findall(r'acc="(SRR\d+)"', runs_field, flags=re.IGNORECASE)
        accession = study_accession or experiment_accession or (run_accessions[0] if run_accessions else "") or sample_accession or submitter_accession
        source_type = "study"
        if study_accession:
            source_type = "study"
        elif experiment_accession:
            source_type = "experiment"
        elif run_accessions:
            source_type = "run"
        elif sample_accession:
            source_type = "sample"
        sample_count = payload.get("sample_count")
        if isinstance(sample_count, str) and sample_count.isdigit():
            sample_count = int(sample_count)
        if sample_count is None and biosample and source_type != "study":
            sample_count = 1
        if not isinstance(sample_count, int):
            sample_count = None
        run_count = payload.get("run_count")
        if isinstance(run_count, str) and run_count.isdigit():
            run_count = int(run_count)
        if not isinstance(run_count, int):
            run_count = count_runs(runs_field)
        return cls(
            uid=coerce_str(payload.get("uid") or payload.get("Id")),
            accession=accession,
            source_type=source_type,
            study_accession=study_accession,
            experiment_accession=experiment_accession,
            sample_accession=sample_accession,
            run_accessions=[item.upper() for item in run_accessions],
            submitter_accession=submitter_accession,
            title=title,
            summary=summary,
            organism=organism,
            library_strategy=library_strategy,
            library_layout=library_layout,
            platform=platform,
            sample_count=sample_count,
            run_count=run_count,
            bioproject=bioproject,
            biosample=biosample,
            raw=payload,
        )

    def merge_query(self, query_label: str) -> None:
        if query_label not in self.matched_queries:
            self.matched_queries.append(query_label)

    def to_dict(self, include_summary: bool, include_raw: bool) -> dict[str, Any]:
        data = {
            "uid": self.uid,
            "accession": self.accession,
            "source_type": self.source_type,
            "study_accession": self.study_accession,
            "experiment_accession": self.experiment_accession,
            "sample_accession": self.sample_accession,
            "run_accessions": self.run_accessions,
            "submitter_accession": self.submitter_accession,
            "title": self.title,
            "organism": self.organism,
            "library_strategy": self.library_strategy,
            "library_layout": self.library_layout,
            "platform": self.platform,
            "sample_count": self.sample_count,
            "run_count": self.run_count,
            "bioproject": self.bioproject,
            "biosample": self.biosample,
            "matched_queries": self.matched_queries,
            "normalized_evidence": build_sra_candidate_evidence(self, include_summary=include_summary),
        }
        if include_summary:
            data["summary"] = self.summary
        if include_raw:
            data["raw"] = plainify(self.raw)
        return data


def esearch(query: str, retmax: int, timeout: int, pause_seconds: float, email: str | None) -> list[str]:
    if biopython_available():
        configure_entrez(email)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*Email address is not specified.*")
            with Entrez.esearch(db="sra", term=query, retmax=retmax, sort="relevance") as handle:
                payload = Entrez.read(handle)
        time.sleep(max(pause_seconds, 0.0))
        return [coerce_str(item) for item in payload.get("IdList", [])]

    params = {
        "db": "sra",
        "term": query,
        "retmode": "json",
        "retmax": str(retmax),
        "sort": "relevance",
        "tool": DEFAULT_TOOL,
    }
    if email:
        params["email"] = email
    payload = request_json("esearch.fcgi", params, timeout=timeout, pause_seconds=pause_seconds)
    return payload.get("esearchresult", {}).get("idlist", [])


def esummary(ids: list[str], timeout: int, pause_seconds: float, email: str | None) -> list[SraRecord]:
    if not ids:
        return []

    if biopython_available():
        configure_entrez(email)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*Email address is not specified.*")
            with Entrez.esummary(db="sra", id=",".join(ids), report="full") as handle:
                payload = plainify(Entrez.read(handle))
        time.sleep(max(pause_seconds, 0.0))
        items: list[Any]
        if isinstance(payload, dict):
            document_set = payload.get("DocumentSummarySet", {})
            if isinstance(document_set, dict):
                raw_items = document_set.get("DocumentSummary", payload)
            else:
                raw_items = document_set
        else:
            raw_items = payload
        if isinstance(raw_items, list):
            items = raw_items
        else:
            items = [raw_items]
        return [SraRecord.from_esummary(item) for item in items]

    params = {
        "db": "sra",
        "id": ",".join(ids),
        "retmode": "json",
        "version": "2.0",
        "tool": DEFAULT_TOOL,
    }
    if email:
        params["email"] = email
    payload = request_json("esummary.fcgi", params, timeout=timeout, pause_seconds=pause_seconds)
    result = payload.get("result", {})
    records: list[SraRecord] = []
    for uid in result.get("uids", []):
        item = result.get(uid)
        if isinstance(item, dict):
            records.append(SraRecord.from_esummary(item))
    return records


def fetch_records(query: str, max_results: int, timeout: int, pause_seconds: float, email: str | None) -> list[SraRecord]:
    ids = esearch(query, retmax=max_results, timeout=timeout, pause_seconds=pause_seconds, email=email)
    records: list[SraRecord] = []
    for batch in chunked(ids, 100):
        records.extend(esummary(batch, timeout=timeout, pause_seconds=pause_seconds, email=email))
    return records


def normalize_query_label(query: str) -> str:
    return " ".join(query.split())


def merge_results(queries: list[str], max_results: int, timeout: int, pause_seconds: float, email: str | None) -> list[SraRecord]:
    merged: dict[str, SraRecord] = {}
    ordered_keys: list[str] = []
    for query in queries:
        label = normalize_query_label(query)
        records = fetch_records(query, max_results=max_results, timeout=timeout, pause_seconds=pause_seconds, email=email)
        for record in records:
            dedupe_key = record.accession or record.uid
            current = merged.get(dedupe_key)
            if current is None:
                record.merge_query(label)
                merged[dedupe_key] = record
                ordered_keys.append(dedupe_key)
            else:
                current.merge_query(label)
                current.run_accessions = sorted(set(current.run_accessions + record.run_accessions))
                current.sample_count = current.sample_count or record.sample_count
                current.run_count = max(current.run_count or 0, record.run_count or 0) or None
                if not current.organism:
                    current.organism = record.organism
                if not current.title:
                    current.title = record.title
                if not current.bioproject:
                    current.bioproject = record.bioproject
                if not current.biosample:
                    current.biosample = record.biosample
        for item in merged.values():
            if item.run_accessions:
                item.run_count = len(set(item.run_accessions))
    return [merged[key] for key in ordered_keys]


def record_text(record: SraRecord) -> str:
    return " ".join(
        part
        for part in [
            record.accession,
            record.title,
            record.summary,
            record.organism,
            record.library_strategy,
            record.library_layout,
            record.platform,
            record.bioproject,
            record.biosample,
        ]
        if part
    )


def filter_records(
    records: list[SraRecord],
    organism_filters: list[str],
    strategy_filters: list[str],
    include_text_filters: list[str],
    exclude_text_filters: list[str],
) -> list[SraRecord]:
    filtered: list[SraRecord] = []
    for record in records:
        if organism_filters and not any(contains_casefold(record.organism, item) for item in organism_filters):
            continue
        if strategy_filters and not any(contains_casefold(record.library_strategy, item) for item in strategy_filters):
            continue
        haystack = record_text(record)
        if include_text_filters and not all(contains_casefold(haystack, item) for item in include_text_filters):
            continue
        if exclude_text_filters and any(contains_casefold(haystack, item) for item in exclude_text_filters):
            continue
        filtered.append(record)
    return filtered


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Search NCBI SRA through E-utilities.")
    parser.add_argument("--query", action="append", required=True, help="Search query. Repeat for multi-query merging.")
    parser.add_argument("--max-results", type=int, default=10, help="Maximum results per query.")
    parser.add_argument("--timeout", type=int, default=30, help="Network timeout in seconds.")
    parser.add_argument("--pause-seconds", type=float, default=0.34, help="Pause between NCBI requests.")
    parser.add_argument("--email", help="Email to send to NCBI as recommended by E-utilities.")
    parser.add_argument("--json", action="store_true", help="Print JSON instead of text.")
    parser.add_argument("--include-summary", action="store_true", help="Include summary in output.")
    parser.add_argument("--include-raw", action="store_true", help="Include raw esummary payload in JSON output.")
    parser.add_argument("--organism", action="append", default=[], help="Post-filter results by organism text.")
    parser.add_argument("--strategy-contains", action="append", default=[], help="Post-filter results by library strategy text.")
    parser.add_argument("--include-text", action="append", default=[], help="Require this text in title, summary, organism, or accession fields.")
    parser.add_argument("--exclude-text", action="append", default=[], help="Exclude records containing this text.")
    return parser.parse_args()


def render_table(records: list[SraRecord], include_summary: bool) -> str:
    if not records:
        return "No results."

    lines: list[str] = []
    for index, record in enumerate(records, start=1):
        lines.append(f"{index}. {record.accession or record.uid}")
        fields = [
            f"title={record.title or 'n/a'}",
            f"organism={record.organism or 'n/a'}",
            f"strategy={record.library_strategy or 'n/a'}",
            f"layout={record.library_layout or 'n/a'}",
            f"platform={record.platform or 'n/a'}",
            f"sample_count={record.sample_count if record.sample_count is not None else 'n/a'}",
            f"run_count={record.run_count if record.run_count is not None else 'n/a'}",
            f"bioproject={record.bioproject or 'n/a'}",
            f"biosample={record.biosample or 'n/a'}",
            f"matched_queries={', '.join(record.matched_queries) or 'n/a'}",
        ]
        for field_line in fields:
            lines.append(f"   {field_line}")
        if include_summary and record.summary:
            wrapped = textwrap.fill(record.summary, width=96, initial_indent="   summary=", subsequent_indent="           ")
            lines.append(wrapped)
    return "\n".join(lines)


def main() -> int:
    configure_stdio()
    args = parse_args()
    try:
        records = merge_results(
            queries=args.query,
            max_results=args.max_results,
            timeout=args.timeout,
            pause_seconds=args.pause_seconds,
            email=args.email,
        )
        records = filter_records(
            records,
            organism_filters=args.organism,
            strategy_filters=args.strategy_contains,
            include_text_filters=args.include_text,
            exclude_text_filters=args.exclude_text,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.json:
        payload = {
            "queries": args.query,
            "result_count": len(records),
            "records": [record.to_dict(include_summary=args.include_summary, include_raw=args.include_raw) for record in records],
        }
        json.dump(payload, sys.stdout, indent=2, ensure_ascii=False)
        sys.stdout.write("\n")
    else:
        print(render_table(records, include_summary=args.include_summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
"""
Search NCBI BioProject through E-utilities and return project-level metadata.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import warnings
from dataclasses import dataclass, field
from typing import Any

from evidence_schema import build_bioproject_candidate_evidence
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


def extract_accessions(text: str, pattern: str) -> list[str]:
    if not text:
        return []
    return sorted(set(match.upper() for match in re.findall(pattern, text, flags=re.IGNORECASE)))


def flatten_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, list):
        return "; ".join(part for part in (flatten_text(item) for item in value) if part)
    if isinstance(value, dict):
        parts = [flatten_text(subvalue) for subvalue in value.values()]
        return "; ".join(part for part in parts if part)
    return coerce_str(value)


@dataclass
class BioProjectRecord:
    uid: str
    accession: str = ""
    title: str = ""
    summary: str = ""
    organism: str = ""
    linked_sra: list[str] = field(default_factory=list)
    linked_biosamples: list[str] = field(default_factory=list)
    raw: dict[str, Any] = field(default_factory=dict)
    matched_queries: list[str] = field(default_factory=list)

    @classmethod
    def from_esummary(cls, payload: dict[str, Any]) -> "BioProjectRecord":
        if not isinstance(payload, dict):
            payload = {"title": flatten_text(payload)}
        title = coerce_str(
            payload.get("Project_Title")
            or payload.get("project_title")
            or payload.get("Project_Title")
            or payload.get("Project_Name")
            or payload.get("title")
            or payload.get("Title")
        )
        summary = coerce_str(
            payload.get("Project_Description")
            or payload.get("project_descr")
            or payload.get("description")
            or payload.get("Description")
        )
        accession = coerce_str(payload.get("Project_Acc") or payload.get("project_acc") or payload.get("accession") or payload.get("Accession"))
        text_blob = " ".join(
            part
            for part in [
                accession,
                title,
                summary,
                flatten_text(payload.get("Project_Data_Type")),
                flatten_text(payload.get("Project_Target")),
                flatten_text(payload.get("Project_Target_Scope")),
                flatten_text(payload.get("Project_Target_Material")),
                flatten_text(payload.get("Organism")),
                flatten_text(payload.get("Organism_Name")),
                flatten_text(payload.get("Organism_Label")),
                flatten_text(payload),
            ]
            if part
        )
        organism = coerce_str(payload.get("Organism")) or coerce_str(payload.get("Organism_Name")) or coerce_str(payload.get("Organism_Label")) or ""
        if not organism:
            organism_matches = re.findall(r"(Homo sapiens|Gallus gallus|Mus musculus|Sus scrofa|Bos taurus)", text_blob, flags=re.IGNORECASE)
            organism = "; ".join(sorted(set(match for match in organism_matches)))
        return cls(
            uid=coerce_str(payload.get("uid") or payload.get("Id")),
            accession=accession,
            title=title,
            summary=summary,
            organism=organism,
            linked_sra=extract_accessions(text_blob, r"\bSR[APXRS]\d+\b"),
            linked_biosamples=extract_accessions(text_blob, r"\bSAM[END][A-Z]?\d+\b"),
            raw=payload,
        )

    def merge_query(self, query_label: str) -> None:
        if query_label not in self.matched_queries:
            self.matched_queries.append(query_label)

    def to_dict(self, include_summary: bool, include_raw: bool) -> dict[str, Any]:
        data = {
            "uid": self.uid,
            "accession": self.accession,
            "title": self.title,
            "organism": self.organism,
            "linked_sra": self.linked_sra,
            "linked_biosamples": self.linked_biosamples,
            "matched_queries": self.matched_queries,
            "normalized_evidence": build_bioproject_candidate_evidence(self, include_summary=include_summary),
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
            with Entrez.esearch(db="bioproject", term=query, retmax=retmax, sort="relevance") as handle:
                payload = Entrez.read(handle)
        time.sleep(max(pause_seconds, 0.0))
        return [coerce_str(item) for item in payload.get("IdList", [])]

    params = {
        "db": "bioproject",
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


def esummary(ids: list[str], timeout: int, pause_seconds: float, email: str | None) -> list[BioProjectRecord]:
    if not ids:
        return []

    if biopython_available():
        configure_entrez(email)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*Email address is not specified.*")
            with Entrez.esummary(db="bioproject", id=",".join(ids)) as handle:
                payload = plainify(Entrez.read(handle))
        time.sleep(max(pause_seconds, 0.0))
        if isinstance(payload, dict):
            document_set = payload.get("DocumentSummarySet", {})
            if isinstance(document_set, dict):
                raw_items = document_set.get("DocumentSummary", payload)
            else:
                raw_items = document_set
        else:
            raw_items = payload
        items = raw_items if isinstance(raw_items, list) else [raw_items]
        return [BioProjectRecord.from_esummary(item) for item in items]

    params = {
        "db": "bioproject",
        "id": ",".join(ids),
        "retmode": "json",
        "version": "2.0",
        "tool": DEFAULT_TOOL,
    }
    if email:
        params["email"] = email
    payload = request_json("esummary.fcgi", params, timeout=timeout, pause_seconds=pause_seconds)
    result = payload.get("result", {})
    records: list[BioProjectRecord] = []
    for uid in result.get("uids", []):
        item = result.get(uid)
        if isinstance(item, dict):
            records.append(BioProjectRecord.from_esummary(item))
    return records


def fetch_records(query: str, max_results: int, timeout: int, pause_seconds: float, email: str | None) -> list[BioProjectRecord]:
    ids = esearch(query, retmax=max_results, timeout=timeout, pause_seconds=pause_seconds, email=email)
    records: list[BioProjectRecord] = []
    for batch in chunked(ids, 100):
        records.extend(esummary(batch, timeout=timeout, pause_seconds=pause_seconds, email=email))
    return records


def merge_results(queries: list[str], max_results: int, timeout: int, pause_seconds: float, email: str | None) -> list[BioProjectRecord]:
    merged: dict[str, BioProjectRecord] = {}
    ordered_keys: list[str] = []
    for query in queries:
        label = " ".join(query.split())
        for record in fetch_records(query, max_results=max_results, timeout=timeout, pause_seconds=pause_seconds, email=email):
            dedupe_key = record.accession or record.uid
            current = merged.get(dedupe_key)
            if current is None:
                record.merge_query(label)
                merged[dedupe_key] = record
                ordered_keys.append(dedupe_key)
            else:
                current.merge_query(label)
                current.linked_sra = sorted(set(current.linked_sra + record.linked_sra))
                current.linked_biosamples = sorted(set(current.linked_biosamples + record.linked_biosamples))
                if not current.title:
                    current.title = record.title
                if not current.summary:
                    current.summary = record.summary
                if not current.organism:
                    current.organism = record.organism
    return [merged[key] for key in ordered_keys]


def filter_records(records: list[BioProjectRecord], organism_filters: list[str], include_text_filters: list[str], exclude_text_filters: list[str]) -> list[BioProjectRecord]:
    filtered: list[BioProjectRecord] = []
    for record in records:
        haystack = " ".join(part for part in [record.accession, record.title, record.summary, record.organism] if part)
        if organism_filters and not any(contains_casefold(record.organism, item) for item in organism_filters):
            continue
        if include_text_filters and not all(contains_casefold(haystack, item) for item in include_text_filters):
            continue
        if exclude_text_filters and any(contains_casefold(haystack, item) for item in exclude_text_filters):
            continue
        filtered.append(record)
    return filtered


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Search NCBI BioProject through E-utilities.")
    parser.add_argument("--query", action="append", required=True, help="Search query. Repeat for multi-query merging.")
    parser.add_argument("--max-results", type=int, default=10, help="Maximum results per query.")
    parser.add_argument("--timeout", type=int, default=30, help="Network timeout in seconds.")
    parser.add_argument("--pause-seconds", type=float, default=0.34, help="Pause between NCBI requests.")
    parser.add_argument("--email", help="Email to send to NCBI as recommended by E-utilities.")
    parser.add_argument("--json", action="store_true", help="Print JSON instead of text.")
    parser.add_argument("--include-summary", action="store_true", help="Include summary in output.")
    parser.add_argument("--include-raw", action="store_true", help="Include raw esummary payload in JSON output.")
    parser.add_argument("--organism", action="append", default=[], help="Post-filter results by organism text.")
    parser.add_argument("--include-text", action="append", default=[], help="Require this text in title, summary, organism, or accession fields.")
    parser.add_argument("--exclude-text", action="append", default=[], help="Exclude records containing this text.")
    return parser.parse_args()


def render_table(records: list[BioProjectRecord], include_summary: bool) -> str:
    if not records:
        return "No results."
    lines: list[str] = []
    for index, record in enumerate(records, start=1):
        lines.append(f"{index}. {record.accession or record.uid}")
        lines.append(f"   title={record.title or 'n/a'}")
        lines.append(f"   organism={record.organism or 'n/a'}")
        lines.append(f"   linked_sra={', '.join(record.linked_sra[:6]) or 'n/a'}")
        lines.append(f"   linked_biosamples={', '.join(record.linked_biosamples[:6]) or 'n/a'}")
        lines.append(f"   matched_queries={', '.join(record.matched_queries) or 'n/a'}")
        if include_summary and record.summary:
            lines.append(f"   summary={record.summary}")
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

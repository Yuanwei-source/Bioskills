#!/usr/bin/env python
"""
Search NCBI GEO/GDS through E-utilities and return candidate study metadata.

This script prefers Biopython's Bio.Entrez client when available and falls
back to direct E-utilities requests otherwise. It focuses on GEO Series (`GSE`)
by default.

Examples:
  python scripts/search_geo.py --query "Gallus gallus adipose RNA-seq" --max-results 5
  python scripts/search_geo.py --query "osteoarthritis cartilage Homo sapiens" --json
  python scripts/search_geo.py --query "GSE12345" --include-summary
"""

from __future__ import annotations

import argparse
import json
import sys
import textwrap
import time
import warnings
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Any

try:
    from Bio import Entrez
except ImportError:  # pragma: no cover - local fallback path
    Entrez = None

from evidence_schema import build_geo_candidate_evidence

if Entrez is not None:
    warnings.filterwarnings("ignore", category=UserWarning, module=r"Bio\.Entrez")


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DEFAULT_TOOL = "biofinder"


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
    if isinstance(value, (int, float)):
        return str(value)
    return str(value).strip()


def coerce_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def request_json(endpoint: str, params: dict[str, Any], timeout: int, pause_seconds: float) -> dict[str, Any]:
    url = f"{EUTILS_BASE}/{endpoint}?{urllib.parse.urlencode(params, doseq=True)}"
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": f"{DEFAULT_TOOL}/0.1 (+https://openai.com)",
            "Accept": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as exc:
        details = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} from NCBI for {endpoint}: {details}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Network error calling NCBI for {endpoint}: {exc}") from exc

    time.sleep(max(pause_seconds, 0.0))
    return payload


def chunked(items: list[str], size: int) -> list[list[str]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def biopython_available() -> bool:
    return Entrez is not None


def configure_entrez(email: str | None) -> None:
    if Entrez is None:
        return
    Entrez.tool = DEFAULT_TOOL
    if email:
        Entrez.email = email


def plainify(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): plainify(subvalue) for key, subvalue in value.items()}
    if isinstance(value, list):
        return [plainify(item) for item in value]
    if isinstance(value, tuple):
        return [plainify(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if hasattr(value, "items"):
        return {str(key): plainify(subvalue) for key, subvalue in value.items()}
    if hasattr(value, "__iter__") and not isinstance(value, (bytes, bytearray, str)):
        try:
            return [plainify(item) for item in value]
        except TypeError:
            pass
    text = coerce_str(value)
    if text.isdigit():
        return int(text)
    return text


def build_search_term(query: str, series_only: bool) -> str:
    cleaned = query.strip()
    if not cleaned:
        raise ValueError("Query must not be empty.")
    if not series_only:
        return cleaned
    return f"({cleaned}) AND gse[ETYP]"


def normalize_query_label(query: str) -> str:
    return " ".join(query.split())


@dataclass
class GeoRecord:
    uid: str
    accession: str = ""
    title: str = ""
    summary: str = ""
    taxon: str = ""
    entry_type: str = ""
    gds_type: str = ""
    platform: str = ""
    samples: int | None = None
    publication_date: str = ""
    supplemental_files: str = ""
    ftp_link: str = ""
    bioproject: str = ""
    raw: dict[str, Any] = field(default_factory=dict)
    matched_queries: list[str] = field(default_factory=list)

    @classmethod
    def from_esummary(cls, payload: dict[str, Any]) -> "GeoRecord":
        sample_count = extract_sample_count(payload)
        return cls(
            uid=coerce_str(payload.get("uid") or payload.get("Id")),
            accession=coerce_str(payload.get("accession") or payload.get("Accession")),
            title=coerce_str(payload.get("title")),
            summary=coerce_str(payload.get("summary")),
            taxon=extract_taxon(payload),
            entry_type=coerce_str(payload.get("entryType") or payload.get("entrytype")),
            gds_type=coerce_str(payload.get("gdsType") or payload.get("gdstype")),
            platform=extract_platform(payload),
            samples=sample_count,
            publication_date=coerce_str(payload.get("PDAT") or payload.get("pdat")),
            supplemental_files=coerce_str(payload.get("suppfile") or payload.get("suppFile")),
            ftp_link=coerce_str(payload.get("ftplink") or payload.get("FTPLink")),
            bioproject=extract_bioproject(payload),
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
            "taxon": self.taxon,
            "entry_type": self.entry_type,
            "gds_type": self.gds_type,
            "platform": self.platform,
            "samples": self.samples,
            "publication_date": self.publication_date,
            "supplemental_files": self.supplemental_files,
            "ftp_link": self.ftp_link,
            "bioproject": self.bioproject,
            "matched_queries": self.matched_queries,
            "normalized_evidence": build_geo_candidate_evidence(self, include_summary=include_summary),
        }
        if include_summary:
            data["summary"] = self.summary
        if include_raw:
            data["raw"] = plainify(self.raw)
        return data


def extract_taxon(payload: dict[str, Any]) -> str:
    taxon = payload.get("taxon")
    if isinstance(taxon, list):
        return ", ".join(coerce_str(x) for x in taxon if coerce_str(x))
    return coerce_str(taxon)


def extract_bioproject(payload: dict[str, Any]) -> str:
    direct = coerce_str(payload.get("bioproject"))
    if direct:
        return direct
    projects = flatten_mixed_field(payload.get("Projects") or payload.get("projects"))
    return projects


def extract_platform(payload: dict[str, Any]) -> str:
    candidates = [
        payload.get("GPL"),
        payload.get("gpl"),
        payload.get("platform"),
        payload.get("Platform"),
    ]
    for candidate in candidates:
        value = flatten_mixed_field(candidate)
        if value:
            return normalize_platform_field(value)
    return ""


def normalize_platform_field(value: str) -> str:
    parts = [part.strip() for part in value.split(";")]
    normalized: list[str] = []
    for part in parts:
        if not part:
            continue
        if part.isdigit():
            normalized.append(f"GPL{part}")
        else:
            normalized.append(part)
    return "; ".join(normalized)


def extract_sample_count(payload: dict[str, Any]) -> int | None:
    direct_keys = ["n_samples", "nSamples", "samples", "Samples", "sample_count"]
    for key in direct_keys:
        value = payload.get(key)
        if isinstance(value, int):
            return value
        if isinstance(value, str) and value.isdigit():
            return int(value)

    sample_list = payload.get("samples") or payload.get("Samples")
    if isinstance(sample_list, list):
        return len(sample_list)

    return None


def flatten_mixed_field(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, list):
        parts: list[str] = []
        for item in value:
            if isinstance(item, dict):
                for key in ("accession", "name", "title", "value"):
                    text = coerce_str(item.get(key))
                    if text:
                        parts.append(text)
                        break
            else:
                text = coerce_str(item)
                if text:
                    parts.append(text)
        return "; ".join(parts)
    if isinstance(value, dict):
        for key in ("accession", "name", "title", "value"):
            text = coerce_str(value.get(key))
            if text:
                return text
    return coerce_str(value)


def esearch(query: str, retmax: int, timeout: int, pause_seconds: float, email: str | None) -> list[str]:
    if biopython_available():
        configure_entrez(email)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*Email address is not specified.*")
            with Entrez.esearch(db="gds", term=query, retmax=retmax, sort="relevance") as handle:
                payload = Entrez.read(handle)
        time.sleep(max(pause_seconds, 0.0))
        return [coerce_str(item) for item in payload.get("IdList", [])]

    params = {
        "db": "gds",
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


def esummary(ids: list[str], timeout: int, pause_seconds: float, email: str | None) -> list[GeoRecord]:
    if not ids:
        return []

    if biopython_available():
        configure_entrez(email)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*Email address is not specified.*")
            with Entrez.esummary(db="gds", id=",".join(ids)) as handle:
                payload = Entrez.read(handle)
        time.sleep(max(pause_seconds, 0.0))
        return [GeoRecord.from_esummary(item) for item in payload]

    params = {
        "db": "gds",
        "id": ",".join(ids),
        "retmode": "json",
        "tool": DEFAULT_TOOL,
    }
    if email:
        params["email"] = email
    payload = request_json("esummary.fcgi", params, timeout=timeout, pause_seconds=pause_seconds)
    result = payload.get("result", {})
    records: list[GeoRecord] = []
    for uid in result.get("uids", []):
        item = result.get(uid)
        if isinstance(item, dict):
            records.append(GeoRecord.from_esummary(item))
    return records


def fetch_records(query: str, max_results: int, timeout: int, pause_seconds: float, series_only: bool, email: str | None) -> list[GeoRecord]:
    search_term = build_search_term(query, series_only=series_only)
    ids = esearch(search_term, retmax=max_results, timeout=timeout, pause_seconds=pause_seconds, email=email)
    records: list[GeoRecord] = []
    for batch in chunked(ids, 100):
        records.extend(esummary(batch, timeout=timeout, pause_seconds=pause_seconds, email=email))
    return records


def contains_casefold(text: str, needle: str) -> bool:
    return needle.casefold() in text.casefold()


def is_single_species_record(taxon: str) -> bool:
    separators = [";", ","]
    return not any(separator in taxon for separator in separators)


def record_text(record: GeoRecord) -> str:
    return " ".join(
        part
        for part in [
            record.accession,
            record.title,
            record.summary,
            record.taxon,
            record.gds_type,
            record.platform,
            record.supplemental_files,
            record.bioproject,
        ]
        if part
    )


def filter_records(
    records: list[GeoRecord],
    taxon_filters: list[str],
    assay_filters: list[str],
    include_text_filters: list[str],
    exclude_text_filters: list[str],
    single_species_only: bool,
) -> list[GeoRecord]:
    filtered: list[GeoRecord] = []
    for record in records:
        if single_species_only and not is_single_species_record(record.taxon):
            continue
        if taxon_filters and not any(contains_casefold(record.taxon, item) for item in taxon_filters):
            continue
        if assay_filters and not any(contains_casefold(record.gds_type, item) for item in assay_filters):
            continue
        haystack = record_text(record)
        if include_text_filters and not all(contains_casefold(haystack, item) for item in include_text_filters):
            continue
        if exclude_text_filters and any(contains_casefold(haystack, item) for item in exclude_text_filters):
            continue
        filtered.append(record)
    return filtered


def merge_results(queries: list[str], max_results: int, timeout: int, pause_seconds: float, series_only: bool, email: str | None) -> list[GeoRecord]:
    merged: dict[str, GeoRecord] = {}
    ordered_uids: list[str] = []

    for query in queries:
        label = normalize_query_label(query)
        records = fetch_records(query, max_results=max_results, timeout=timeout, pause_seconds=pause_seconds, series_only=series_only, email=email)
        for record in records:
            current = merged.get(record.uid)
            if current is None:
                record.merge_query(label)
                merged[record.uid] = record
                ordered_uids.append(record.uid)
            else:
                current.merge_query(label)

    return [merged[uid] for uid in ordered_uids]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Search GEO/GDS through NCBI E-utilities.")
    parser.add_argument("--query", action="append", required=True, help="Search query. Repeat for multi-query merging.")
    parser.add_argument("--max-results", type=int, default=10, help="Maximum results per query.")
    parser.add_argument("--timeout", type=int, default=30, help="Network timeout in seconds.")
    parser.add_argument("--pause-seconds", type=float, default=0.34, help="Pause between NCBI requests.")
    parser.add_argument("--email", help="Email to send to NCBI as recommended by E-utilities.")
    parser.add_argument("--json", action="store_true", help="Print JSON instead of a text table.")
    parser.add_argument("--include-summary", action="store_true", help="Include summary in output.")
    parser.add_argument("--include-raw", action="store_true", help="Include raw esummary payload in JSON output.")
    parser.add_argument("--taxon", action="append", default=[], help="Post-filter results by taxon text. Repeatable.")
    parser.add_argument("--assay-contains", action="append", default=[], help="Post-filter results by assay or gds_type text. Repeatable.")
    parser.add_argument("--include-text", action="append", default=[], help="Require this text in title/summary/taxon/platform. Repeatable.")
    parser.add_argument("--exclude-text", action="append", default=[], help="Exclude records containing this text in title/summary/taxon/platform. Repeatable.")
    parser.add_argument("--single-species-only", action="store_true", help="Exclude multi-species records when taxon metadata lists several organisms.")
    parser.add_argument(
        "--entry-type",
        choices=["series", "any"],
        default="series",
        help="Restrict results to GSE series by default.",
    )
    return parser.parse_args()


def render_table(records: list[GeoRecord], include_summary: bool) -> str:
    if not records:
        return "No results."

    lines: list[str] = []
    for index, record in enumerate(records, start=1):
        header = f"{index}. {record.accession or record.uid}"
        title = record.title or "(no title)"
        fields = [
            f"title={title}",
            f"taxon={record.taxon or 'n/a'}",
            f"type={record.entry_type or record.gds_type or 'n/a'}",
            f"platform={record.platform or 'n/a'}",
            f"samples={record.samples if record.samples is not None else 'n/a'}",
            f"suppfiles={record.supplemental_files or 'n/a'}",
            f"bioproject={record.bioproject or 'n/a'}",
            f"matched_queries={', '.join(record.matched_queries) or 'n/a'}",
        ]
        if record.publication_date:
            fields.append(f"pdat={record.publication_date}")
        lines.append(header)
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
            series_only=(args.entry_type == "series"),
            email=args.email,
        )
        records = filter_records(
            records,
            taxon_filters=args.taxon,
            assay_filters=args.assay_contains,
            include_text_filters=args.include_text,
            exclude_text_filters=args.exclude_text,
            single_species_only=args.single_species_only,
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

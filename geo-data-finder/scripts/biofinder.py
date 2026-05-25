#!/usr/bin/env python
"""
Thin orchestration entry point for BioFinder.

This v2 foundation wraps the GEO-first MVP workflow:
1. search candidate GSE accessions
2. inspect top candidates
3. emit a common evidence shape for later multi-provider use
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from inspect_gse import inspect_accession, normalize_goal
from search_geo import configure_stdio, filter_records, merge_results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the BioFinder GEO-first workflow.")
    parser.add_argument("--query", required=True, help="Primary search query.")
    parser.add_argument("--query-variant", action="append", default=[], help="Additional query variants to merge.")
    parser.add_argument("--analysis-goal", default="general", help="Analysis goal such as DEG, WGCNA, or validation.")
    parser.add_argument("--max-results", type=int, default=10, help="Maximum GEO results per query.")
    parser.add_argument("--inspect-top", type=int, default=3, help="Inspect up to this many top GEO candidates.")
    parser.add_argument("--timeout", type=int, default=30, help="Network timeout in seconds.")
    parser.add_argument("--pause-seconds", type=float, default=0.34, help="Pause between NCBI requests.")
    parser.add_argument("--email", help="Email to send to NCBI as recommended by E-utilities.")
    parser.add_argument("--taxon", action="append", default=[], help="Post-filter results by taxon text.")
    parser.add_argument("--assay-contains", action="append", default=[], help="Post-filter GEO results by assay or gds_type text.")
    parser.add_argument("--include-text", action="append", default=[], help="Require this text in title, summary, taxon, or platform.")
    parser.add_argument("--exclude-text", action="append", default=[], help="Exclude records containing this text.")
    parser.add_argument("--single-species-only", action="store_true", help="Exclude multi-species GEO results.")
    parser.add_argument("--entry-type", choices=["series", "any"], default="series", help="Restrict results to GSE by default.")
    parser.add_argument("--destdir", default="C:/tmp/geo_data_finder", help="Destination directory for GEOparse cache files.")
    parser.add_argument("--how", choices=["full", "quick"], default="quick", help="GEOparse retrieval mode for inspection.")
    parser.add_argument("--json", action="store_true", help="Print JSON output.")
    parser.add_argument("--include-summary", action="store_true", help="Include GEO summary text in candidate output.")
    parser.add_argument("--include-raw", action="store_true", help="Include raw inspection notes.")
    return parser.parse_args()


def build_search_brief(args: argparse.Namespace, analysis_goal: str) -> dict[str, Any]:
    return {
        "user_query": args.query,
        "query_variants": list(args.query_variant),
        "analysis_goal": analysis_goal,
        "filters": {
            "taxon": list(args.taxon),
            "assay_contains": list(args.assay_contains),
            "include_text": list(args.include_text),
            "exclude_text": list(args.exclude_text),
            "single_species_only": args.single_species_only,
            "entry_type": args.entry_type,
        },
    }


def candidate_payload(records: list[Any], include_summary: bool) -> list[dict[str, Any]]:
    return [record.to_dict(include_summary=include_summary, include_raw=False) for record in records]


def inspected_payload(results: list[Any], include_raw: bool) -> list[dict[str, Any]]:
    return [result.to_dict(include_raw=include_raw) for result in results]


def normalized_bundle(candidates: list[dict[str, Any]], inspected: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if inspected:
        return [item.get("normalized_evidence", {}) for item in inspected]
    return [item.get("normalized_evidence", {}) for item in candidates]


def run_geo_workflow(args: argparse.Namespace) -> dict[str, Any]:
    analysis_goal = normalize_goal(args.analysis_goal)
    queries = [args.query, *args.query_variant]
    records = merge_results(
        queries=queries,
        max_results=args.max_results,
        timeout=args.timeout,
        pause_seconds=args.pause_seconds,
        series_only=(args.entry_type == "series"),
        email=args.email,
    )
    filtered = filter_records(
        records,
        taxon_filters=args.taxon,
        assay_filters=args.assay_contains,
        include_text_filters=args.include_text,
        exclude_text_filters=args.exclude_text,
        single_species_only=args.single_species_only,
    )

    inspectable = [record for record in filtered if str(record.accession).upper().startswith("GSE")]
    top_records = inspectable[: max(args.inspect_top, 0)]
    destdir = Path(args.destdir)
    destdir.mkdir(parents=True, exist_ok=True)
    inspected = [
        inspect_accession(
            accession=record.accession,
            analysis_goal=analysis_goal,
            destdir=destdir,
            how=args.how,
            email=args.email,
        )
        for record in top_records
    ]
    candidate_data = candidate_payload(filtered, include_summary=args.include_summary)
    inspected_data = inspected_payload(inspected, include_raw=args.include_raw)

    return {
        "workflow": "biofinder_geo_first",
        "provider_scope": ["GEO"],
        "analysis_goal": analysis_goal,
        "search_brief": build_search_brief(args, analysis_goal),
        "candidate_count": len(filtered),
        "candidates": candidate_data,
        "inspected_count": len(inspected),
        "inspected": inspected_data,
        "normalized_evidence_bundle": normalized_bundle(candidate_data, inspected_data),
    }


def render_text(payload: dict[str, Any]) -> str:
    lines = [
        "BioFinder GEO-first workflow",
        f"analysis_goal={payload['analysis_goal']}",
        f"candidate_count={payload['candidate_count']}",
        f"inspected_count={payload['inspected_count']}",
        f"user_query={payload['search_brief']['user_query']}",
    ]
    query_variants = payload["search_brief"]["query_variants"]
    if query_variants:
        lines.append(f"query_variants={'; '.join(query_variants)}")
    if payload["candidates"]:
        lines.append("")
        lines.append("Top candidates:")
        for index, item in enumerate(payload["candidates"][:5], start=1):
            evidence = item.get("normalized_evidence", {})
            lines.append(
                f"{index}. {item.get('accession') or 'n/a'} | {item.get('title') or 'n/a'} | "
                f"taxon={', '.join(evidence.get('organisms', [])) or item.get('taxon') or 'n/a'} | "
                f"samples={evidence.get('sample_count', item.get('samples', 'n/a'))}"
            )
    if payload["inspected"]:
        lines.append("")
        lines.append("Inspected accessions:")
        for item in payload["inspected"]:
            lines.append(
                f"- {item.get('accession')}: provisional={item.get('provisional_local_assessment')} | "
                f"usable_for={', '.join(item.get('usable_for', [])) or 'n/a'} | "
                f"risks={', '.join(item.get('risks', [])) or 'n/a'}"
            )
    return "\n".join(lines)


def main() -> int:
    configure_stdio()
    args = parse_args()
    try:
        payload = run_geo_workflow(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if args.json:
        json.dump(payload, sys.stdout, indent=2, ensure_ascii=False)
        sys.stdout.write("\n")
    else:
        print(render_text(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

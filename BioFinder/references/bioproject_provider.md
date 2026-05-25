# BioProject Provider

Use this reference when BioFinder needs project-level context beyond what GEO or
SRA already expose.

BioProject is not a first-search provider for most requests. It is a context
provider.

## What BioProject Is Good For

Use BioProject when one or more of these are true:

- you already have a `PRJ...` accession
- GEO or SRA results point to a BioProject and you need project-level scope
- you need to understand whether a result belongs to a narrow study or a broad
  project container
- you need a central accession that links SRA and BioSample records together

BioProject is especially useful after GEO or SRA has already surfaced promising
records.

## What BioProject Should Return

BioProject output should emphasize:

- project accession
- title and summary
- linked SRA studies, experiments, or runs when available
- linked BioSample accessions when recoverable
- organism clues
- project scope clues
- obvious complexity warnings such as broad multi-study containers

## Current Runtime Pieces

- `scripts/search_bioproject.py`
- `scripts/inspect_bioproject.py`
- `scripts/evidence_schema.py`

These scripts should stay light. They should support project discovery and
context aggregation, not full downstream processing.

## Current Limits

The current v2 BioProject layer is intentionally limited:

- it does not resolve every cross-link perfectly
- it does not replace SRA inspection
- it may expose sparse organism or tissue evidence
- it may need the LLM to decide whether a broad project is still worth keeping

## Judgment Boundary

Local BioProject scripts should extract project-level context and linked
accession clues.

The LLM should still decide:

- whether the project is too broad or too weakly matched
- whether the project helps the user's current analysis goal
- whether the BioProject is useful only as supporting context or as a real
  shortlist candidate

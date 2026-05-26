# SRA Provider

Use this reference when BioFinder needs to search or inspect NCBI SRA as part of
the v2 workflow.

SRA is not a replacement for GEO. It fills a different role:

- GEO is better for expression-study discovery and sample-design reading
- SRA is better for raw sequencing availability, run-level context, and project links

## When To Use SRA

Use SRA when one or more of these are true:

- the user explicitly asks for raw sequencing data
- GEO results are weak but the biological target is still clear
- you need SRP, SRX, SRR, or SRS accessions for downstream download planning
- you need to confirm whether sequencing runs really exist behind a GEO-like idea
- you need BioProject or BioSample context that GEO does not expose clearly

Do not switch to SRA just because it exists. If GEO already gives a clear and
usable answer, keep the workflow simple.

## Current SRA Runtime Pieces

- `scripts/search_sra.py`
- `scripts/inspect_sra.py`
- `scripts/evidence_schema.py`

These scripts should stay evidence-oriented. They should not grow into download
or analysis pipelines.

## Metadata Preference Order

For experimental-design clues, prefer:

1. NCBI RunInfo or Run Selector style tabular metadata
2. `pysradb` metadata output as an optional enhancement layer
3. SRA `esummary` or `ExpXml` only for fields that are not available elsewhere

This keeps BioFinder closer to existing metadata tooling and avoids building a
large custom parser around XML fragments alone.

## What SRA Should Return

SRA output should emphasize:

- study accession such as `SRP...` when possible
- linked experiment, sample, and run accessions
- library strategy
- library layout
- platform or instrument
- run count
- BioProject and BioSample links
- raw-sequencing availability

## Current Limits

The current v2 SRA layer is intentionally conservative:

- it does not download FASTQ
- it does not run SRA Toolkit
- it does not claim perfect tissue or phenotype rescue from sparse titles
- it may only recover partial design clues such as sample-name groupings when
  richer sample attributes are absent
- it may still need the LLM to interpret whether a raw-data-only study is truly relevant

## Judgment Boundary

Local SRA scripts should extract and aggregate evidence.

The LLM should still decide:

- whether the SRA study matches the user's biological target
- whether it is suitable for DEG, WGCNA, validation, or follow-up
- whether sparse metadata makes it too risky to recommend

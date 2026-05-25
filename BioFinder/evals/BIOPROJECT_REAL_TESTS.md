# BioProject Real Tests

This file records early real-world validation of the BioFinder v2 BioProject provider.

## Planned Checks

1. `search_bioproject.py`

Goal:

- verify that `PRJ...` accessions can be discovered
- verify that project title and summary are extractable
- verify that linked `SRP/SRX/SRR/SRS` and `SAM...` accessions surface when present

2. `inspect_bioproject.py`

Goal:

- verify that one BioProject accession can be aggregated into project-level context
- verify that linked SRA and BioSample counts are exposed
- verify that broad-project warnings are surfaced conservatively

## Notes

BioProject should be evaluated mainly as a context provider.

Success does not require rich phenotype detail. It requires:

- stable project accession handling
- useful cross-provider links
- conservative project-scope signaling

## Test 1: search_bioproject.py

Command shape:

```text
python scripts/search_bioproject.py --query "PRJNA1453495" --max-results 3 --json
```

Observed outcome after parser fixes:

- returned the correct project accession: `PRJNA1453495`
- extracted a useful title:
  - `Transcriptome of Ileum in Broilers with High and Low Abdominal Fat`
- extracted project summary text
- did not directly expose linked `SRA` or `BioSample` accessions from BioProject `esummary`

Issues discovered and fixed:

- BioProject field names in the real payload used `Project_Title` and `Project_Description`, not the lower-case variants originally assumed
- organism was not explicitly exposed in the BioProject summary payload for this project

Current status:

- title and summary parsing are now aligned to the real payload
- BioProject-only metadata is still sparse for cross-provider links

## Test 2: inspect_bioproject.py

Command shape:

```text
python scripts/inspect_bioproject.py --accession PRJNA1453495 --analysis-goal differential_expression --json --include-raw
```

Observed outcome after SRA backfill:

- accession: `PRJNA1453495`
- organism: `Gallus gallus`
- tissue clue: `abdominal fat`
- linked SRA accessions recovered:
  - `SRP691785`
  - `SRX32944932` to `SRX32944937`
  - `SRR38097878` to `SRR38097883`
- linked BioSamples recovered:
  - `SAMN57266310` to `SAMN57266315`
- run count: `6`
- provisional local assessment: `medium`

Issues discovered and fixed:

- project-level inspection initially returned no linked SRA or BioSample accessions
- the fix was to use a light SRA backfill query only when BioProject metadata itself was too sparse

Current status:

- BioProject now works as a useful context provider
- it still behaves conservatively and should not replace GEO or SRA as the main discovery surface

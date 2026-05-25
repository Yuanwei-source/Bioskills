# Real Search Tests

This document records real NCBI GEO search tests for `scripts/search_geo.py` after the Biopython refactor.

## Environment

- backend: `Bio.Entrez` from `biopython`
- live source: NCBI E-utilities `gds`
- script target: broad candidate retrieval, not final dataset recommendation

## Test 1: Human osteoarthritis cartilage RNA-seq

Command shape:

```text
python scripts/search_geo.py \
  --query "osteoarthritis cartilage Homo sapiens RNA-seq" \
  --max-results 8 \
  --taxon "Homo sapiens" \
  --assay-contains "Expression profiling by high throughput sequencing" \
  --json
```

What happened before filtering:

- a bovine record was retrieved because the title and summary contained human cartilage wording
- a small-RNA record was retrieved because it still matched the broad RNA-seq wording

What happened after filtering:

- species drift was removed by `--taxon "Homo sapiens"`
- assay-family drift was reduced by `--assay-contains "Expression profiling by high throughput sequencing"`
- the result set became much tighter

Main lesson:

- free-text GEO retrieval is noisy even when the user gives a precise disease, tissue, species, and assay
- explicit post-filters are not optional for realistic use

## Test 2: Chicken abdominal fat MVP scenario

Command shape:

```text
python scripts/search_geo.py \
  --query "Gallus gallus abdominal fat transcriptome" \
  --query "chicken abdominal fat adipose" \
  --max-results 5 \
  --taxon "Gallus gallus" \
  --single-species-only \
  --include-text "fat" \
  --json
```

What improved:

- the multi-species bodymap result disappeared with `--single-species-only`
- GPL accessions were normalized from raw digits to `GPLxxxxx`
- multi-query merging worked and preserved `matched_queries`

What still slipped through:

- a skin-color chicken dataset still appeared
- a PBMC-focused subseries still appeared
- a large adipose-plus-other-tissues array study appeared, which may be useful but needs careful inspection

Main lesson:

- broad search can surface biologically adjacent but analysis-misaligned studies
- `search_geo.py` can narrow the pool, but it cannot replace accession-level inspection

## Test 3: Raw payload behavior

Command shape:

```text
python scripts/search_geo.py --query "..." --json --include-raw
```

Observed behavior:

- large raw payloads can make real runs slow
- the Biopython path works, but raw output is heavy and not ideal for routine search passes

Main lesson:

- `--include-raw` is best used for debugging one small query, not for everyday candidate retrieval

## Confirmed Improvements

- The script now prefers `Bio.Entrez` instead of hand-rolled E-utilities calls.
- Species and assay post-filters materially improve result quality.
- Multi-species false positives can be removed with `--single-species-only`.
- Platform IDs are normalized into `GPLxxxxx` style.
- The default output is cleaner because the noisy Biopython email warning is suppressed.

## Current Limitations

- `Bio.Entrez` esummary output does not expose every field the raw JSON endpoint exposes. In the current tests, `BioProject` was not available from the parsed response.
- Free-text matching still admits semantically related but unhelpful datasets.
- The script does not yet understand study design, superseries structure, or whether a hit is actually a validation-quality cohort.
- The script is still a retrieval layer. Recommendation quality still depends on a later `inspect_gse.py` stage.

## Recommended Next Steps

1. Build `scripts/inspect_gse.py` to inspect top candidates in detail.
2. Add a small structured-query wrapper that translates search briefs into filtered `search_geo.py` arguments.
3. Add regression tests for false positives such as multi-species bodymaps, PBMC subseries, and cell-line contamination.

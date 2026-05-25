# Real Inspection Tests

This document records live GEO accession tests for `scripts/inspect_gse.py`.

## Test 1: GSE292022 with WGCNA goal

Command shape:

```text
python scripts/inspect_gse.py \
  --accession GSE292022 \
  --analysis-goal WGCNA \
  --json
```

Observed behavior:

- assay family was identified as `bulk_rna_seq`
- tissue was correctly reduced to `Abdominal fat`
- sample count was verified from `GSM` objects as `6`
- treatment grouping was recovered from `characteristics_ch1.1.treatment`
- WGCNA was rejected because the sample size is too small
- overall usability was reduced to `medium` because the requested goal was not supported

Why this is good:

- this is the right kind of nuance for a real triage step
- the accession is not junk, but it should not be sold as a WGCNA dataset

Remaining issue:

- some text fields from GEOparse show encoding artifacts in temperature strings

## Test 2: GSE104042 with WGCNA goal

Command shape:

```text
python scripts/inspect_gse.py \
  --accession GSE104042 \
  --analysis-goal WGCNA \
  --json
```

Observed behavior before heuristic tightening:

- source-name values were flooding the tissue summary
- `data_row_count` and `series_id` were incorrectly treated as useful grouping columns
- the large sample size risked making the series look too attractive

Observed behavior after heuristic tightening:

- tissue summary reduced to four real tissues
- `data_row_count` and `series_id` no longer appear as group signals
- `SuperSeries` is explicitly flagged
- `multiple_tissues` is explicitly flagged
- WGCNA is rejected at the whole-series level

Why this is good:

- the script now distinguishes a biologically rich series from a clean single-cohort dataset
- it no longer confuses "many samples" with "good direct WGCNA candidate"

## Main Lessons

- GEOparse is a strong fit for the inspection layer because `phenotype_data` quickly exposes sample-level metadata.
- GEOparse does not remove the need for our own judgment layer.
- mixed-tissue and superseries cases need explicit heuristics or they will be overrated.
- requested-goal awareness matters: the same accession can be useful for DEG but weak for WGCNA.

## Current Gaps

- no explicit subseries extraction yet for superseries
- no direct ranking of which subgroup or subseries should be used next
- no special handling yet for small RNA, single-cell, or non-human mixed-model edge cases beyond basic assay-family flags
- no normalization yet for common GEO encoding artifacts in metadata text

## Recommended Next Steps

1. Add optional subgroup-focused inspection such as `--prefer-tissue "Abdominal fat"`.
2. Add superseries-aware suggestions when `series_id` or subseries structure is present.
3. Add more assay-specific rules for small RNA and single-cell studies.
4. Normalize obvious text-encoding artifacts in phenotype strings.

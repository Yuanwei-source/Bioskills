# GEO Data Finder Test Scenarios

This document records simulated user flows against the current `BioFinder` skill and the main optimization points they reveal.

## Scenario 1: Clear bulk RNA-seq search

User request:

```text
Find chicken abdominal fat GEO transcriptome datasets for WGCNA.
```

Expected behavior:

- Parse `chicken`, `abdominal fat`, and `WGCNA`
- Prefer bulk RNA-seq, then microarray
- Require larger sample size than a simple DEG analysis
- Return a ranked shortlist, not an unfiltered GEO dump

What worked:

- The skill already supported structured search briefs
- It already recognized WGCNA as a distinct analysis goal

What was weak:

- "transcriptome" was not explicitly mapped to bulk RNA-seq first, microarray second
- No explicit sample-size heuristic for WGCNA
- No rule for how many candidate accessions to inspect deeply

Optimization added:

- Default assay ordering for vague transcriptome requests
- WGCNA sample-size heuristic
- Deep-inspection cap for broad searches

## Scenario 2: Broad disease search with missing fields

User request:

```text
Find diabetes islet GEO data.
```

Expected behavior:

- Notice that assay type and analysis goal are missing
- Decide whether to ask a follow-up or proceed with assumptions
- Avoid mixing whole pancreas, isolated islet, and single-cell results without clear labels

What worked:

- The skill already had clarification rules

What was weak:

- It did not explicitly say how to handle vague assay words like "data" or "transcriptome"
- It did not force bulk-before-single-cell ordering for broad requests
- It did not highlight that fallback matches should be labeled separately

Optimization added:

- Default assay ordering
- Widening search rounds
- Explicit labeling for fallback matches

## Scenario 3: Inspect a known accession

User request:

```text
Can GSEXXXXX be used for differential expression?
```

Expected behavior:

- Skip broad search
- Inspect study design, sample labels, processed matrix, and raw-count availability
- Give a firm verdict with risks

What worked:

- The skill already had an accession-inspection mode
- It already required design and data-availability checks

What was weak:

- No explicit rule against inferring sample counts from title-only summaries
- No explicit evidence-versus-inference labeling when GSM metadata is messy

Optimization added:

- Verified / partially verified / unclear sample-status labels
- Stronger warning against overstating sample counts and labels

## Scenario 4: Compare several candidate studies

User request:

```text
Which is better for osteoarthritis cartilage DEG analysis: GSEA, GSEB, or GSEC?
```

Expected behavior:

- Inspect all provided accessions
- Compare on a shared rubric
- Explain why one should be the primary dataset

What worked:

- The skill already supported comparison as a task type
- The output template already had recommendation fields

What was weak:

- The ranking output did not require an explicit explanation when several datasets were close
- The scoring rubric did not define floor conditions that block a `High` rating

Optimization added:

- `Why these ranked this way` note for multi-dataset ranking
- `High` rating guardrails

## Scenario 5: Validation cohort request

User request:

```text
Find a validation GEO dataset for human osteoarthritis cartilage RNA-seq.
```

Expected behavior:

- Prioritize biological independence over superficial similarity
- Warn if only cross-platform validation is available

What worked:

- The skill already recognized validation as a task type

What was weak:

- It did not define independence clearly
- It did not warn about cross-platform harmonization burden

Optimization added:

- Validation-specific guidance in scoring and download/use recommendations

## Scenario 6: Single-cell ambiguity

User request:

```text
Find cartilage transcriptome GEO data for osteoarthritis.
```

Expected behavior:

- Prefer bulk datasets first unless the user explicitly wants single-cell
- Avoid mixing single-cell exploratory datasets into the top ranks without warning

What worked:

- The scope already warned that single-cell needs extra care

What was weak:

- Single-cell was still easy to over-include in MVP searches
- There was no explicit rule for when single-cell belongs only as exploratory fallback

Optimization added:

- Single-cell is now framed as explicit or fallback-only in MVP

## Scenario 7: Superseries trap

User request:

```text
Find liver fibrosis GEO datasets with clear case and control groups.
```

Expected behavior:

- Detect when a `GSE` is a superseries containing several substudies or multiple assay types
- Avoid recommending it as one clean cohort if the actual design is fragmented

What worked:

- The general inspection checklist covered design review

What was weak:

- Superseries and mixed-platform traps were not named explicitly

Optimization added:

- Superseries / subseries / mixed-platform inspection rules

## Main Remaining Gaps

The skill is stronger after the updates, but these gaps still remain:

- No explicit handling yet for multiple accessions in one mixed request such as "find data and compare these two known GSEs too"

## Recommended Next Steps

1. Expand final-output examples and LLM-facing evaluation examples as new real GEO edge cases are discovered.
2. Add explicit handling for mixed requests such as "find new data and compare these known GSEs too".
3. Extend `inspect_gse.py` only where it improves evidence extraction quality, not to replicate nuanced LLM judgment in local rules.

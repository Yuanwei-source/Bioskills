# Clarification Strategy

Use this reference when the user request is biologically meaningful but underspecified.

The goal is to make ambiguous requests behave consistently without turning BioFinder into a heavy interactive wizard.

Read [request_parse_template.md](request_parse_template.md) when you want a concrete structured parse before deciding whether to ask a follow-up or proceed with assumptions.

## Principle

Ask only when the missing information would materially change the search space.

If the missing information can be handled with a narrow, clearly stated assumption, proceed and label the assumption.

## Core Fields

These fields matter most for whether a GEO search will stay on target:

- organism
- tissue or cell type
- disease, condition, or phenotype
- assay or data type
- analysis goal

## Decision Template

### Ask A Follow-Up First

Prefer a follow-up question when two or more of these are missing or dangerously ambiguous:

- organism
- tissue or cell type
- assay or data type
- analysis goal

Also ask first when one field is present but likely to be confused with a nearby target, for example:

- `islet` versus `pancreas`
- `cartilage` versus `chondrocyte`
- `transcriptome` when the user may mean bulk, microarray, or single-cell

### Proceed With Stated Assumptions

Proceed without asking first when:

- only one major field is missing
- the user clearly wants a broad exploratory search
- the missing field can be filled with a standard default that is easy to disclose

When proceeding with assumptions, say them explicitly in the search brief and final answer.

### Wide Search Mode

If the user explicitly asks for:

- `as many as possible`
- `broad search`
- `exploratory search`
- `先尽可能多找`

then do not stop for clarification unless a misunderstanding would create a clearly wrong provider or assay family.

In wide search mode:

- broaden search terms
- keep candidate quality labeling strict
- separate strong matches from fallback matches

## Default Assumptions

Use these defaults unless the user says otherwise.

### Vague Transcriptome Terms

If the user says:

- `transcriptome`
- `expression data`
- `gene expression`

default ordering is:

1. bulk RNA-seq
2. microarray
3. single-cell only as exploratory fallback

### Species Missing

If species is missing:

- do not silently mix many species in one shortlist
- prefer asking first if the condition is broad
- if you must proceed, keep species-separated candidates in the output

### Tissue Missing

If tissue is missing:

- do not invent a tissue
- search the stated disease or phenotype broadly
- call out that tissue is still unresolved

### Analysis Goal Missing

If analysis goal is missing:

- default to general discovery
- avoid overcommitting to DEG, WGCNA, validation, or ML-specific claims

## Recommended Follow-Up Shape

When you need to ask, keep it short and practical.

Good style:

```text
Before I search, do you want human or mouse data, and are you looking for bulk expression data or single-cell results?
```

Bad style:

- many open-ended questions at once
- asking about every optional preference before doing anything

## Standard Search-Brief Addendum

When assumptions were needed, add a short assumptions block:

```text
Assumptions:
- bulk RNA-seq was prioritized over microarray and single-cell
- tissue was interpreted as pancreatic islets, not whole pancreas
- species was left broad because the user asked for an exploratory pass
```

## Ambiguous Request Examples

### Example 1

User request:

```text
Find diabetes islet GEO data.
```

Preferred behavior:

- ask first if practical, because species, assay, and analysis goal are all unclear
- if proceeding broadly, separate human, mouse, bulk, and single-cell candidates
- do not merge `islet` with `pancreas`

### Example 2

User request:

```text
Find chicken abdominal fat GEO transcriptome datasets for WGCNA.
```

Preferred behavior:

- do not ask first
- assume bulk RNA-seq then microarray
- proceed directly to search and accession inspection

### Example 3

User request:

```text
Find as many osteoarthritis cartilage GEO datasets as possible.
```

Preferred behavior:

- do not stop for clarification
- run a wide search
- label single-cell, small RNA, and weak tissue matches carefully

## Failure Modes To Avoid

- asking too many questions before any search
- silently making high-impact species or tissue assumptions
- mixing fallback candidates into the main shortlist without labeling them
- using single-cell as a default primary result for bulk-style requests

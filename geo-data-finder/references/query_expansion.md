# Query Expansion

Use this reference when the user request needs term expansion before GEO search. The goal is to improve recall without letting the search drift away from the biological question.

This file does not define a large synonym dictionary. It defines how to ask the LLM for bounded expansions and how to correct high-risk cases with a small rule layer.

## Purpose

Use LLM-driven query expansion to:

- normalize informal user wording into searchable biomedical terms
- generate a few high-value aliases for organism, tissue, disease, phenotype, and assay
- preserve the user's real target while improving GEO recall

Do not use expansion to:

- invent adjacent diseases the user did not ask for
- broaden one assay into all omics modalities
- replace GEO evidence with semantic guesswork

## Expansion Workflow

1. Extract the search brief fields from the user request.
2. Mark each field as `explicit`, `implicit`, or `missing`.
3. Normalize high-risk fields with the small rule layer first.
4. Ask the LLM for a short list of expansions only for fields that benefit from expansion.
5. Build search rounds from narrow to broad.
6. Stop expanding when the new terms begin changing the biological target.

## Field Rules

Apply different expansion pressure by field.

- `organism`: low expansion pressure. Prefer official species name plus the most common common-name form.
- `tissue_or_cell_type`: medium expansion pressure. Allow direct synonyms and one-hop biological variants.
- `disease_or_condition`: medium expansion pressure. Allow common aliases, abbreviations, and tightly equivalent disease names.
- `phenotype_or_trait`: medium expansion pressure. Allow direct phrasing variants and closely related measurement terms.
- `assay_or_platform`: low expansion pressure. Expand into platform-level search names only within the same assay family.
- `analysis_goal`: do not expand as a search field unless it provides concrete searchable terms such as survival, time series, or treatment.

## Expansion Budget

Keep expansions short and bounded.

- `organism`: 2 to 3 terms
- `tissue_or_cell_type`: 3 to 6 terms
- `disease_or_condition`: 3 to 6 terms
- `phenotype_or_trait`: 2 to 5 terms
- `assay_or_platform`: 2 to 4 terms

If the user request is already precise, use the low end of the range.

## Allowed Expansion Types

Good expansions:

- official names and common names
- singular and plural variants when relevant
- standard abbreviations
- direct biomedical aliases
- one-hop biological variants that still point to the same target

Examples:

- `chicken` -> `Gallus gallus`, `chicken`
- `abdominal fat` -> `abdominal fat`, `adipose tissue`, `fat tissue`
- `osteoarthritis` -> `osteoarthritis`, `degenerative joint disease`
- `RNA-seq` -> `RNA-seq`, `expression profiling by high throughput sequencing`

Bad expansions:

- pathway-level concepts that are only loosely related
- neighboring diseases with different pathogenesis
- different tissues from the same organ system unless used only as fallback
- assay jumps such as bulk RNA-seq to proteomics or methylation

## Stop Conditions

Stop expansion and keep the search narrow if any of these happen:

- terms start moving from the requested tissue to a nearby but different tissue
- disease terms start broadening from a diagnosis to a symptom cluster
- assay terms begin crossing modality boundaries
- single-cell terms begin crowding a bulk-oriented search
- the expanded query no longer sounds like the original biological question

When in doubt, keep the user's original term and reduce the expansion set.

## Small Rule Layer

Use a small rule layer for high-risk normalization. This is not a large vocabulary store. It only protects boundaries that commonly derail GEO search.

### Organism Standards

Normalize to official plus common forms:

- `human` -> `Homo sapiens`, `human`
- `mouse` -> `Mus musculus`, `mouse`, `mice`
- `rat` -> `Rattus norvegicus`, `rat`
- `chicken` -> `Gallus gallus`, `chicken`
- `zebrafish` -> `Danio rerio`, `zebrafish`

Do not broaden species without explicit reason.

### Assay Boundaries

Normalize vague assay wording conservatively:

- `transcriptome` -> prefer `bulk RNA-seq`, then `microarray`, then `single-cell RNA-seq` as exploratory fallback
- `expression` -> map to GEO expression profiling families, not all omics
- `sequencing` -> do not assume RNA-seq unless the biological context supports transcriptomic intent

Do not cross these assay families during expansion:

- bulk RNA-seq
- microarray
- single-cell RNA-seq
- methylation
- ATAC-seq
- proteomics

### Easy-To-Confuse Tissue Pairs

Treat these as distinct unless the user explicitly broadens the search:

- `islet` versus `pancreas`
- `cartilage` versus `chondrocyte`
- `blood` versus `PBMC`
- `adipose tissue` versus `adipocyte`
- `liver` versus `hepatocyte`
- `tumor tissue` versus `cell line`

If fallback search needs a nearby tissue or cell type, label it as fallback and rank it lower.

### Bulk Versus Single-Cell Defaults

Unless the user explicitly asks for single-cell:

- prefer bulk RNA-seq
- then microarray
- include single-cell only as exploratory fallback

If the user explicitly asks for single-cell, do not mix bulk datasets into the main ranking unless the user also asks for broader transcriptome data.

## Prompt Pattern

Use a bounded prompt pattern like this internally:

```text
You are helping construct GEO search terms.
Given this structured brief:
- organism: ...
- tissue_or_cell_type: ...
- disease_or_condition: ...
- phenotype_or_trait: ...
- assay_or_platform: ...
- analysis_goal: ...

Generate a small controlled list of search-term expansions.
Rules:
- keep the same biological target
- prefer direct aliases, standard abbreviations, official names, and one-hop variants
- do not broaden to different omics types
- do not broaden to nearby but different tissues unless clearly marked fallback
- return at most 3 organism terms, 5 tissue terms, 5 disease terms, 4 phenotype terms, and 3 assay terms
- if a field is already precise, keep the list short
```

You do not need the LLM to produce prose. A compact structured term set is enough.

## Search Round Construction

Build search rounds in this order:

1. Original terms with normalized organism and assay
2. Original terms plus high-confidence aliases
3. Original terms plus controlled phenotype or disease variants
4. Fallback terms using nearby tissue or assay-adjacent terms only if the first rounds fail

Do not merge fallback rounds into the primary shortlist without labeling them.

## Evidence Rule

Expanded terms only help retrieve candidate accessions. They do not prove relevance.

Every final recommendation must still be checked against GEO evidence:

- title
- summary
- overall design
- GSM annotations
- platform
- supplemental files

## Worked Examples

### Example 1

User request:

```text
Find chicken abdominal fat transcriptome GEO data for WGCNA.
```

Good expansion:

- organism: `Gallus gallus`, `chicken`
- tissue: `abdominal fat`, `adipose tissue`, `fat tissue`
- phenotype: `fat deposition`, `adiposity`
- assay: `RNA-seq`, `expression profiling by high throughput sequencing`, `microarray`

Avoid:

- `liver`
- `obesity` if it shifts away from the abdominal-fat context
- `proteomics`

### Example 2

User request:

```text
Find diabetes islet GEO data.
```

Good expansion:

- organism: keep broad only if not specified
- tissue: `islet`, `pancreatic islet`, `islets of Langerhans`
- disease: `diabetes`, `diabetes mellitus`

Avoid:

- silently replacing `islet` with `pancreas`
- mixing bulk and single-cell as equal-ranked results

### Example 3

User request:

```text
Find osteoarthritis cartilage RNA-seq GEO datasets for external validation.
```

Good expansion:

- tissue: `cartilage`, `articular cartilage`
- disease: `osteoarthritis`, `degenerative joint disease`
- assay: `RNA-seq`, `expression profiling by high throughput sequencing`

Fallback only:

- `chondrocyte`

Because cell-based studies may not validate tissue-level findings cleanly.

## Failure Mode Checklist

When expansion quality looks poor, check for these failure modes:

- too many terms for one concept
- disease drift
- tissue drift
- assay drift
- species drift
- broad jargon replacing the user's real phrasing

If any failure mode appears, shrink the expansion set and rerun the narrow search first.

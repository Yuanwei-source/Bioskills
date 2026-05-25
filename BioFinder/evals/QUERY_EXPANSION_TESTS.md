# Query Expansion Tests

This document records a simulated multi-round test pass for `references/query_expansion.md`.

The goal of these tests is not to prove GEO retrieval quality end to end. The goal is to verify that the expansion rules stay biologically aligned, do not drift across assay families, and produce bounded search-term sets.

## Pass Criteria

Each test passes if:

- expansions stay within the requested biological target
- high-risk fields are normalized correctly
- the expansion set stays small
- fallback terms are clearly separated from primary terms
- no cross-omics drift occurs

## Round 1: Clear bulk RNA-seq intent

User request:

```text
Find chicken abdominal fat transcriptome GEO data for WGCNA.
```

Expected normalization:

- organism -> `Gallus gallus`, `chicken`
- tissue -> `abdominal fat`, `adipose tissue`, `fat tissue`
- assay ordering -> `bulk RNA-seq`, `microarray`, `single-cell RNA-seq` fallback only

Observed decision quality:

- good alignment with tissue target
- good bulk-first assay ordering
- no need to expand toward obesity or liver by default

Verdict:

`Pass`

## Round 2: Vague disease plus high-risk tissue

User request:

```text
Find diabetes islet GEO data.
```

Expected normalization:

- keep `islet` distinct from `pancreas`
- permit `pancreatic islet` and `islets of Langerhans`
- avoid broad pancreas fallback unless the narrow rounds fail

Observed decision quality:

- correct tissue boundary behavior
- correctly prevents silent tissue drift
- still leaves room for fallback labeling if exact matches are weak

Verdict:

`Pass`

## Round 3: Tissue versus cell-type ambiguity

User request:

```text
Find osteoarthritis cartilage transcriptome data.
```

Expected normalization:

- primary tissue terms -> `cartilage`, `articular cartilage`
- `chondrocyte` should be fallback, not co-equal primary term

Observed decision quality:

- correct separation of tissue-level and cell-level concepts
- good guardrail against over-ranking cell culture studies

Verdict:

`Pass`

## Round 4: Explicit single-cell request

User request:

```text
Find single-cell RNA-seq GEO datasets for human osteoarthritis synovium.
```

Expected normalization:

- main assay -> `single-cell RNA-seq`
- do not mix bulk RNA-seq into main ranking
- tissue stays `synovium` or direct synonyms only

Observed decision quality:

- correct assay-family lock
- correct prevention of bulk contamination

Verdict:

`Pass`

## Round 5: Over-broad assay wording

User request:

```text
Find liver fibrosis sequencing data in GEO.
```

Expected normalization:

- do not expand `sequencing` into all sequencing modalities
- require transcriptomic interpretation only if the surrounding request supports it

Observed decision quality:

- rules are conservative enough
- this case may still require a clarification question in the live skill if analysis goal is also missing

Verdict:

`Pass with caution`

Reason:

- bounded expansion rules are good, but live execution still depends on clarification logic from `SKILL.md`

## Round 6: Phenotype drift risk

User request:

```text
Find GEO datasets about fat deposition in chicken.
```

Expected normalization:

- phenotype can expand to `adiposity`
- tissue should not be invented if the user never specified one
- avoid jumping to broad obesity terms unless marked fallback

Observed decision quality:

- phenotype expansion stays mostly safe
- the rules correctly discourage tissue invention

Verdict:

`Pass`

## Round 7: Validation-cohort request

User request:

```text
Find an external validation GEO dataset for osteoarthritis cartilage RNA-seq.
```

Expected normalization:

- preserve `cartilage` as the main target
- avoid broadening to `chondrocyte` unless fallback
- keep assay in the RNA-seq family

Observed decision quality:

- good protection against swapping tissue studies for cell studies
- good protection against cross-omics drift

Verdict:

`Pass`

## Round 8: Species drift risk

User request:

```text
Find mouse kidney fibrosis transcriptome GEO data.
```

Expected normalization:

- organism -> `Mus musculus`, `mouse`, `mice`
- do not drift to rat or human unless the user explicitly broadens species

Observed decision quality:

- species normalization is tight
- no pressure toward related model organisms

Verdict:

`Pass`

## Round 9: Cell line contamination risk

User request:

```text
Find tumor tissue RNA-seq GEO datasets for hepatocellular carcinoma.
```

Expected normalization:

- keep `tumor tissue` distinct from `cell line`
- cell-line studies can appear only as fallback if clearly labeled

Observed decision quality:

- small rule layer handles this boundary well
- useful for preventing misleading validation candidates

Verdict:

`Pass`

## Summary Findings

Strong points:

- bulk versus single-cell boundaries are much clearer now
- tissue drift is meaningfully reduced
- organism normalization is tight without needing a huge static dictionary
- fallback concepts are easier to keep separate from primary search intent

Weak points:

- very vague requests like `sequencing data` still depend heavily on clarification rules
- phenotype-only requests may still need stronger guidance on when not to invent tissues
- disease abbreviations are not enumerated yet, so consistency will depend on the live LLM prompt quality

## Recommended Follow-Ups

1. Add one short section to `query_expansion.md` on abbreviation handling such as `OA`, `T2D`, and `HCC`.
2. Add one short section on when phenotype-only searches should remain tissue-agnostic.
3. Add regression examples that combine two ambiguity sources at once, such as vague assay plus vague tissue.

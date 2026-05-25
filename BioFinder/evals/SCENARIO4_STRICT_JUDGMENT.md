# Scenario 4 Strict Judgment Regression

This file records a stricter LLM-style judgment pass for the scenario:

```text
Find datasets suitable for osteoarthritis cartilage differential analysis and external validation.
```

The goal of this pass is to check whether final LLM judgment correctly downgrades candidates that look acceptable to local heuristics but are weak for the user's real goal.

## Candidate 1: GSE324566

### Extracted Evidence

- organism: `Homo sapiens`
- tissue: `Cartilage`
- assay family: `bulk_rna_seq`
- sample count: `8`
- design clarity: `high`
- group signal: `Vehicle` versus `Peptide p[63-82]`
- local assessment: `high`

### Strict LLM Judgment

- `rating`: `Low`
- `best_use`: mechanistic perturbation DEG context
- `main_risks`: treated chondrocyte perturbation model, not a clean external validation cohort
- `final recommendation`: keep only as a conditional mechanistic comparison dataset, not as a main validation cohort

### Why It Was Downgraded

- the design is clear, but it is a peptide-treatment experiment
- the user's goal is disease-relevant primary analysis plus external validation
- this accession is more useful for mechanistic comparison than cohort validation

## Candidate 2: GSE260927

### Extracted Evidence

- organism: `Homo sapiens`
- tissue evidence from phenotype metadata: `H9 HESC`
- assay family: `bulk_rna_seq`
- sample count: `27`
- design clarity: `high`
- risks: `multiomic_series`, `mixed_assay_design`
- local assessment: `medium`

### Strict LLM Judgment

- `rating`: `Exclude`
- `best_use`: none for the stated validation task
- `main_risks`: stem-cell / in-vitro developmental context, mixed assay framing, not a clean OA cartilage cohort
- `final recommendation`: exclude from the cartilage validation shortlist

### Why It Was Downgraded

- the extracted phenotype evidence points away from a straightforward cartilage tissue cohort
- the study structure is developmental and intervention-heavy
- it does not match the user's intended disease-validation use closely enough

## Candidate 3: GSE287861

### Extracted Evidence

- organism: `Homo sapiens`
- tissue: `Cartilage`
- assay family: `bulk_rna_seq`
- sample count: `6`
- design clarity: `high`
- group signal: `Young` versus `Aged`
- local assessment: `high`

### Strict LLM Judgment

- `rating`: `Low`
- `best_use`: aging-related cartilage comparison
- `main_risks`: aging contrast is not the same as an osteoarthritis validation cohort
- `final recommendation`: mention only as a weak conditional comparison dataset unless the user explicitly accepts age-driven cartilage biology as a proxy

### Why It Was Downgraded

- the tissue and species match
- but the biological contrast is aging rather than a clean OA case-control or independent validation cohort

## Regression Conclusion

This stricter pass shows that:

- local evidence extraction is useful but not sufficient for final ranking
- the final LLM judgment must heavily weight biological fit, not just assay family and group clarity
- validation tasks are especially sensitive to overrating perturbation, developmental, or proxy-biology datasets

## Expected Behavioral Rule

For scenario 4 style requests:

- do not rate mechanistic perturbation datasets as strong validation cohorts
- do not rate developmental or stem-cell systems as disease-validation cohorts unless the user explicitly wants them
- do not let clear groups and moderate sample size overpower biological mismatch

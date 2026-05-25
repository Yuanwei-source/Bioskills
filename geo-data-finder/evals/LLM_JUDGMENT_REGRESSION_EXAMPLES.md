# LLM Judgment Regression Examples

Use these examples to regression-test the final LLM judgment step.

Each example represents the stage after local evidence extraction. The local scripts have already retrieved candidate metadata and accession-level evidence. The task now is to judge dataset usability in context.

The expected outputs here are not meant to force identical wording. They are meant to anchor:

- the direction of the judgment
- the rating
- the main reasons
- the distinction between evidence and inference

Read this file together with [dataset_scoring.md](../references/dataset_scoring.md).

## Output Contract

For each example, the LLM should be able to produce something like:

- `rating`
- `best_use`
- `main_risks`
- `evidence_summary`
- `inference_notes`
- `final_recommendation`

## Example 1: Clear DEG candidate but weak for WGCNA

### User Goal

```text
Find chicken abdominal fat GEO transcriptome data for WGCNA.
```

### Extracted Evidence

```json
{
  "accession": "GSE292022",
  "title": "Effects of chronic heat exposure on gene expression in abdominal fat of growing broiler chickens.",
  "taxon": "Gallus gallus",
  "assay_family": "bulk_rna_seq",
  "sample_count": 6,
  "sample_count_status": "verified_from_gsms",
  "design_clarity": "high",
  "tissues": ["Abdominal fat"],
  "group_columns": [
    {
      "column": "characteristics_ch1.1.treatment",
      "counts": {
        "Thermoneutral": 3,
        "Heat stress": 3
      }
    }
  ],
  "provisional_local_assessment": "medium",
  "usable_for": ["differential_expression", "validation"],
  "not_recommended_for": ["wgcna", "machine_learning"]
}
```

### Expected Judgment

- `rating`: `Medium`
- `best_use`: differential expression
- `main_risks`: very small sample size for WGCNA
- `final_recommendation`: keep as a DEG candidate, not a WGCNA primary dataset

### Why

- tissue, species, and assay match the request well
- groups are clear
- sample size is too small for a strong WGCNA recommendation

### Evidence Versus Inference

- verified: tissue, species, assay family, 3 vs 3 treatment structure
- inferred: weak for WGCNA because network analysis usually needs substantially more samples

### Common Failure To Avoid

- Do not rate this `High` just because the biology matches well.

## Example 2: Large dataset but wrong whole-series interpretation

### User Goal

```text
Find chicken abdominal fat GEO datasets for WGCNA.
```

### Extracted Evidence

```json
{
  "accession": "GSE104042",
  "title": "Transcriptome profile of liver, adipose tissue, muscle and PBMC in genetically fat and lean chicken lines submitted to high-fiber and high-fat diet versus standard high-starch diet",
  "taxon": "Gallus gallus",
  "assay_family": "microarray",
  "sample_count": 178,
  "sample_count_status": "verified_from_gsms",
  "design_clarity": "partial",
  "tissues": [
    "Liver",
    "Abdominal fat",
    "Pectoralis major muscle",
    "Peripheral Blood Mononuclear Cells"
  ],
  "risks": ["superseries", "multiple_tissues"],
  "provisional_local_assessment": "medium",
  "usable_for": ["differential_expression", "validation"],
  "not_recommended_for": ["whole_series_wgcna", "wgcna"]
}
```

### Expected Judgment

- `rating`: `Low` for whole-series WGCNA
- `best_use`: possible tissue-specific follow-up, not direct whole-series use
- `main_risks`: superseries and mixed tissues
- `final_recommendation`: do not use the full accession directly for WGCNA; inspect subseries or the abdominal-fat subset only

### Why

- sample size looks attractive
- but the accession structure is not a clean single-cohort dataset for the requested use

### Evidence Versus Inference

- verified: multiple tissues, superseries summary, microarray assay
- inferred: abdominal-fat subset may still be useful, but that requires deeper subgroup inspection

### Common Failure To Avoid

- Do not rate this `High` just because the sample count is large.

## Example 3: Species drift false positive

### User Goal

```text
Find human osteoarthritis cartilage RNA-seq GEO datasets.
```

### Extracted Evidence

```json
{
  "accession": "GSE296290",
  "title": "Leveraging single cell multiomic analyses to identify gene regulatory networks that drive human articular cartilage cell fate",
  "taxon": "Bos taurus",
  "assay_family": "bulk_rna_seq",
  "gds_type": "Expression profiling by high throughput sequencing",
  "sample_count": 6,
  "summary": "The title and summary mention human articular cartilage development and bovine deep zone articular chondrocytes."
}
```

### Expected Judgment

- `rating`: `Exclude`
- `best_use`: none for the stated human dataset request
- `main_risks`: wrong species despite superficially matching wording
- `final_recommendation`: exclude from the human candidate shortlist

### Why

- the user asked for human data
- the accession taxon is not human

### Evidence Versus Inference

- verified: accession taxon is `Bos taurus`
- no inference is needed for the exclusion decision

### Common Failure To Avoid

- Do not keep this as a medium candidate just because the title mentions human cartilage.

## Example 4: Assay-family mismatch

### User Goal

```text
Find human osteoarthritis cartilage RNA-seq datasets for bulk differential expression.
```

### Extracted Evidence

```json
{
  "accession": "GSE316835",
  "taxon": "Homo sapiens",
  "assay_family": "non_coding_rna_seq",
  "gds_type": "Non-coding RNA profiling by high throughput sequencing",
  "sample_count": 13,
  "tissues": ["cartilage"],
  "summary": "small RNA-seq in osteoarthritis cartilage and LPS-treated chondrocytes"
}
```

### Expected Judgment

- `rating`: `Exclude` for bulk mRNA DEG
- `best_use`: maybe non-coding RNA follow-up, not bulk mRNA DEG
- `main_risks`: wrong assay family
- `final_recommendation`: exclude from the bulk mRNA shortlist

### Why

- the request is specifically about bulk RNA-seq style transcriptome analysis
- this accession belongs to a different sequencing assay family

### Evidence Versus Inference

- verified: assay family is non-coding RNA sequencing
- inferred: could be interesting for a different project, but not for the stated one

### Common Failure To Avoid

- Do not rate this `Low`; for the stated workflow it should be `Exclude`.

## Example 5: Broadly relevant validation candidate with caveats

### User Goal

```text
Find a validation GEO dataset for human osteoarthritis cartilage RNA-seq.
```

### Extracted Evidence

```json
{
  "accession": "GSE324566",
  "taxon": "Homo sapiens",
  "assay_family": "bulk_rna_seq",
  "sample_count": 8,
  "design_clarity": "high",
  "tissues": ["articular chondrocytes"],
  "group_columns": [
    {
      "column": "treatment",
      "counts": {
        "vehicle": 4,
        "peptide": 4
      }
    }
  ],
  "summary": "primary human articular chondrocytes exposed to BMP7-derived peptide"
}
```

### Expected Judgment

- `rating`: `Low` or `Medium` depending on how strict the user is
- `best_use`: mechanistic perturbation comparison, not ideal external validation
- `main_risks`: treatment model in chondrocytes rather than a clean tissue-level independent OA cohort
- `final_recommendation`: do not present as the best validation cohort unless stronger tissue-level cohorts are absent

### Why

- human and cartilage-related biology match
- but the design is more perturbation-focused than cohort-validation focused

### Evidence Versus Inference

- verified: assay family, species, sample count, treatment contrast
- inferred: weaker validation fit because the biological setting is not a clean independent clinical cohort

### Common Failure To Avoid

- Do not call this a `High` validation dataset just because it is human cartilage-related RNA-seq.

## Example 6: Single-cell result in a bulk-oriented request

### User Goal

```text
Find osteoarthritis cartilage transcriptome GEO data for differential expression.
```

### Extracted Evidence

```json
{
  "accession": "EXAMPLE_SC_001",
  "taxon": "Homo sapiens",
  "assay_family": "single_cell_rna_seq",
  "sample_count": 5,
  "design_clarity": "partial",
  "tissues": ["cartilage"],
  "risks": ["single_cell_data"]
}
```

### Expected Judgment

- `rating`: `Low`
- `best_use`: exploratory cell-state analysis
- `main_risks`: assay-family mismatch for bulk-style DEG
- `final_recommendation`: keep only as exploratory context unless the user explicitly accepts single-cell results

### Why

- the user did not explicitly ask for single-cell
- single-cell can be biologically relevant but does not answer the same analysis question

### Evidence Versus Inference

- verified: assay family is single-cell
- inferred: limited fit for the user's intended bulk-like DEG workflow

### Common Failure To Avoid

- Do not rank this beside bulk cohorts as if they were interchangeable.

## Example 7: Borderline but salvageable candidate

### User Goal

```text
Find chicken abdominal fat expression datasets for external validation.
```

### Extracted Evidence

```json
{
  "accession": "GSE104042",
  "taxon": "Gallus gallus",
  "assay_family": "microarray",
  "sample_count": 178,
  "design_clarity": "partial",
  "tissues": [
    "Liver",
    "Abdominal fat",
    "Pectoralis major muscle",
    "Peripheral Blood Mononuclear Cells"
  ],
  "risks": ["superseries", "multiple_tissues"]
}
```

### Expected Judgment

- `rating`: `Medium`
- `best_use`: tissue-focused validation after subgroup selection
- `main_risks`: whole accession is mixed, so validation should not use all samples together
- `final_recommendation`: keep as a conditional candidate if the abdominal-fat subset can be isolated cleanly

### Why

- for validation, a mixed accession is not automatically useless
- but it needs subgroup-aware handling

### Evidence Versus Inference

- verified: the accession is mixed across tissues
- inferred: abdominal-fat subset may still be useful

### Common Failure To Avoid

- Do not `Exclude` the dataset too early just because the whole accession is mixed.

## Example 8: No usable design evidence

### User Goal

```text
Can this accession support differential expression?
```

### Extracted Evidence

```json
{
  "accession": "EXAMPLE_WEAK_001",
  "taxon": "Homo sapiens",
  "assay_family": "bulk_rna_seq",
  "sample_count": 4,
  "design_clarity": "low",
  "tissues": ["cartilage"],
  "group_columns": [],
  "risks": ["partial_metadata_only"]
}
```

### Expected Judgment

- `rating`: `Exclude`
- `best_use`: none for direct DEG
- `main_risks`: no reliable grouping or design evidence
- `final_recommendation`: do not use for direct differential expression without substantial manual rescue

### Why

- the main missing piece is not biology, it is interpretable design

### Evidence Versus Inference

- verified: no useful grouping evidence was recovered
- inferred: manual rescue might still be possible, but that should not be treated as the default path

### Common Failure To Avoid

- Do not rate this `Medium` just because the assay family looks right.

## How To Use These Examples

Use these examples when:

- revising the LLM prompt for final dataset judgment
- comparing different prompt styles
- checking whether script changes shift the final recommendation style
- evaluating whether the system is drifting toward over-permissive or over-conservative recommendations

The key regression goal is not identical wording. The key goal is stable judgment behavior.

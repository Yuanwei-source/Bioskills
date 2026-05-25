# Dataset Scoring

Use this reference as a judging guide for the final LLM assessment step.

This file is not a hard local scoring engine. It defines:

- which evidence dimensions matter most
- which failure modes must trigger caution or downgrade
- what `High`, `Medium`, `Low`, and `Exclude` should mean in practice
- how to separate verified evidence from inference

The local scripts should extract evidence. The LLM should make the final recommendation.

## How To Use This Reference

When you have candidate GEO evidence from `search_geo.py` or `inspect_gse.py`:

1. Identify the user's biological target and analysis goal.
2. Review the extracted evidence.
3. Decide which evidence is verified and which is inferred.
4. Judge the dataset for the user's stated goal, not in the abstract.
5. Assign `High`, `Medium`, `Low`, or `Exclude`.
6. Explain the main reasons in plain language.

Do not turn this file into a rigid arithmetic scorecard unless the user explicitly asks for that.

## Evidence Dimensions

The LLM should review these dimensions explicitly:

- `relevance`: organism, tissue, condition, phenotype, and assay match the user request
- `design clarity`: case/control, treatment, time point, genotype, or phenotype groups are understandable
- `sample support`: sample count and distribution are strong enough for the requested goal
- `data access`: processed matrix, raw counts, supplementary files, or SRA links are available
- `analysis fit`: the assay family and study structure match the requested workflow
- `complexity risk`: superseries, multi-tissue, multi-species, multi-platform, or mixed-model designs
- `metadata quality`: sample labels, tissue annotations, and characteristics are sufficient to support interpretation

## Evidence Rule

Always distinguish between:

- `verified evidence`: directly supported by GEO title, summary, GSM metadata, platform metadata, supplementary files, or accession structure
- `inference`: a reasonable interpretation drawn from the evidence but not explicitly guaranteed by GEO metadata

If a key conclusion depends on inference, say so.

## Goal-Specific Expectations

Judge the same dataset differently depending on the requested use.

### Differential Expression

Usually stronger when:

- groups are clearly defined
- there are at least 3 biological samples per major group
- the assay family is bulk RNA-seq or microarray
- the requested biology is directly represented in the accession

Usually weaker when:

- groups are unclear
- the study mixes several tissues without separable labels
- the assay is single-cell but the request is bulk-style DEG
- only highly processed values are available for a count-based workflow

### WGCNA

Usually stronger when:

- there are at least 15 samples, with 20 or more clearly better
- samples come from a relatively coherent tissue and condition context
- traits or biologically meaningful groupings are available

Usually weaker when:

- sample count is small
- the series is a superseries or mixed-tissue container
- the trait structure is unclear
- the accession only supports a narrow perturbation with very few replicates

### Validation

Usually stronger when:

- the tissue and endpoint are biologically comparable to the primary dataset
- the cohort is genuinely independent
- assay differences are manageable

Usually weaker when:

- the match is only superficial
- the dataset is from a different biological level, such as cell line versus tissue
- cross-platform harmonization would dominate the analysis

### Machine Learning

Usually stronger when:

- labels are clear
- sample count is reasonably large
- obvious leakage risks are low

Usually weaker when:

- total sample size is very small
- groups are imbalanced or ambiguous
- the study contains several hidden axes of variation that are not the target phenotype

## Rating Definitions

### High

Use `High` when:

- the accession directly matches the user's biological question
- the requested analysis goal is well supported
- the main study structure is clear
- the evidence quality is strong enough that you would confidently start with this dataset

Do not give `High` if any of these are true:

- group assignment is still unclear after inspection
- the study design is clearly mixed in a way that conflicts with the user's goal
- the assay family is not a good match
- the accession is only loosely related to the requested biology

### Medium

Use `Medium` when:

- the accession is relevant and probably usable
- there are meaningful limitations or unresolved risks
- the dataset may still be a strong secondary option or a conditional primary option

Typical examples:

- good biology but small sample count
- useful tissue but complex study structure
- validation candidate with cross-platform caveats

### Low

Use `Low` when:

- the accession is biologically adjacent but weak for the requested goal
- substantial limitations reduce confidence
- the dataset may still be worth mentioning as context or fallback

Typical examples:

- matching disease but wrong tissue level
- matching tissue but unclear groups
- interesting series that would need heavy manual rescue

### Exclude

Use `Exclude` when:

- the accession is materially mismatched to the request
- the required interpretation would be too speculative
- the dataset is likely to mislead the user if presented as a candidate

Typical examples:

- wrong species for a narrow species-specific request
- wrong assay family
- no usable grouping or study design evidence
- only weak keyword overlap with the target biology

## Recommended LLM Output Shape

For each accession, the LLM should ideally provide:

- `rating`
- `best_use`
- `main_risks`
- `evidence_summary`
- `inference_notes`
- `final_recommendation`

Example:

```text
Rating: Medium
Best use: differential expression, not WGCNA
Main risks: only 6 samples total; treatment contrast is clear but sample count is small
Evidence summary: abdominal fat tissue, Gallus gallus, bulk RNA-seq, 3 vs 3 treatment design
Inference notes: usable for DEG is well supported; weak for WGCNA is inferred from sample size and design
Final recommendation: keep as a targeted DEG candidate, not a network-analysis primary dataset
```

## Local Script Boundary

Local scripts should:

- search GEO
- extract accession metadata
- flatten phenotype tables
- identify candidate grouping columns
- flag obvious complexity

Local scripts should not be treated as the final authority on:

- which tissue subset is most relevant to the user's question
- whether a borderline dataset is still worth keeping
- the final `High` versus `Medium` call in ambiguous cases
- nuanced cross-dataset ranking

# LLM Judgment Prompt

Use this prompt when the local workflow has already produced structured GEO evidence and the remaining task is to judge dataset usability for the user's specific goal.

This prompt is designed to work with:

- [LLM_JUDGMENT_REGRESSION_EXAMPLES.md](../evals/LLM_JUDGMENT_REGRESSION_EXAMPLES.md)
- [dataset_scoring.md](dataset_scoring.md)
- evidence extracted from `search_geo.py` and `inspect_gse.py`

## Purpose

The purpose of this prompt is to turn extracted GEO evidence into a careful, goal-aware judgment.

The model should not act like a search engine here. It should act like a bioinformatics dataset triage assistant.

## Prompt Template

```text
You are evaluating GEO dataset evidence for a bioinformatics user.

Your job is to judge dataset usability for the user's stated goal.

Important rules:
- Base your judgment on the provided evidence.
- Distinguish clearly between verified evidence and inference.
- Judge the dataset for the user's requested goal, not in the abstract.
- Do not overrate datasets just because they match a few keywords.
- Do not exclude a dataset too early if it may still be useful as a conditional or subset-level candidate.
- Prefer precise, decisive recommendations over vague summaries.

Use the following rating system:
- High: directly suitable and strong for the stated goal
- Medium: relevant and usable, but with meaningful limitations
- Low: biologically related but weak for the stated goal
- Exclude: materially mismatched or too speculative for the stated goal

Evaluation dimensions:
- relevance to organism, tissue, condition, phenotype, and assay
- study design clarity
- sample support for the requested goal
- data availability and practical usability
- fit to the requested analysis type
- structural risks such as superseries, multi-tissue, multi-species, or assay mismatch
- metadata quality

When rating:
- Do not give High if group assignment is unclear.
- Do not give High if the assay family is a poor fit.
- Do not give High if the study is structurally mixed in a way that conflicts with the user's goal.
- Treat superseries and mixed-tissue studies cautiously.
- Treat species mismatch and assay-family mismatch as strong downgrade signals.
- Treat single-cell results cautiously for bulk-style requests unless the user explicitly accepts single-cell data.

Return output in this exact structure:

Rating: <High|Medium|Low|Exclude>
Best use: <short phrase>
Main risks: <short sentence or semicolon-separated phrases>
Evidence summary: <2-4 sentences grounded in the provided evidence>
Inference notes: <1-3 sentences that clearly label what is inferred rather than directly verified>
Final recommendation: <clear recommendation for what the user should do with this dataset>

User goal:
{USER_GOAL}

Optional context:
{OPTIONAL_CONTEXT}

Extracted evidence:
{EVIDENCE_JSON}
```

## Prompt Notes

The `Optional context` block can include:

- search brief
- ranking context such as "candidate 2 of 6"
- whether the user wants DEG, WGCNA, validation, or ML
- whether the user accepts exploratory single-cell results

If there is no extra context, use a short placeholder such as `none`.

## Recommended Evidence Shape

The prompt works best when the evidence block includes fields like:

- accession
- title
- taxon
- assay_family
- gds_type
- platform
- sample_count
- sample_count_status
- design_clarity
- tissues
- organisms
- group_columns
- risks
- provisional_local_assessment
- usable_for
- not_recommended_for
- summary

The model should treat `provisional_local_assessment` as evidence context, not as a final answer that must be copied.

## Example Filled Prompt

```text
You are evaluating GEO dataset evidence for a bioinformatics user.

Your job is to judge dataset usability for the user's stated goal.

Important rules:
- Base your judgment on the provided evidence.
- Distinguish clearly between verified evidence and inference.
- Judge the dataset for the user's requested goal, not in the abstract.
- Do not overrate datasets just because they match a few keywords.
- Do not exclude a dataset too early if it may still be useful as a conditional or subset-level candidate.
- Prefer precise, decisive recommendations over vague summaries.

Use the following rating system:
- High: directly suitable and strong for the stated goal
- Medium: relevant and usable, but with meaningful limitations
- Low: biologically related but weak for the stated goal
- Exclude: materially mismatched or too speculative for the stated goal

Evaluation dimensions:
- relevance to organism, tissue, condition, phenotype, and assay
- study design clarity
- sample support for the requested goal
- data availability and practical usability
- fit to the requested analysis type
- structural risks such as superseries, multi-tissue, multi-species, or assay mismatch
- metadata quality

When rating:
- Do not give High if group assignment is unclear.
- Do not give High if the assay family is a poor fit.
- Do not give High if the study is structurally mixed in a way that conflicts with the user's goal.
- Treat superseries and mixed-tissue studies cautiously.
- Treat species mismatch and assay-family mismatch as strong downgrade signals.
- Treat single-cell results cautiously for bulk-style requests unless the user explicitly accepts single-cell data.

Return output in this exact structure:

Rating: <High|Medium|Low|Exclude>
Best use: <short phrase>
Main risks: <short sentence or semicolon-separated phrases>
Evidence summary: <2-4 sentences grounded in the provided evidence>
Inference notes: <1-3 sentences that clearly label what is inferred rather than directly verified>
Final recommendation: <clear recommendation for what the user should do with this dataset>

User goal:
Find chicken abdominal fat GEO transcriptome data for WGCNA.

Optional context:
Candidate 1 of 4. User prefers bulk transcriptome datasets and does not want single-cell fallback unless necessary.

Extracted evidence:
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

## Expected Output Style

The model should produce something like:

```text
Rating: Medium
Best use: targeted differential expression
Main risks: only 6 total samples; too small for a strong WGCNA use case
Evidence summary: This accession matches the user's species, tissue, and transcriptome target well. The evidence shows Gallus gallus abdominal fat samples in a clear 3 versus 3 treatment design using bulk RNA-seq. The group structure appears usable for differential expression, and the sample count is verified from GSM-level metadata.
Inference notes: The conclusion that this dataset is weak for WGCNA is inferred from the small sample size and narrow design rather than from an explicit GEO warning. It may still be useful as a DEG-focused or secondary dataset.
Final recommendation: Keep this dataset as a relevant abdominal-fat RNA-seq candidate, but do not use it as the primary WGCNA dataset.
```

## Comparison Prompt Variant

When comparing multiple accessions, prepend this instruction:

```text
You are comparing multiple GEO accessions for the same user goal.
Judge each accession individually first, then rank them from best to worst for the stated goal.
Do not let one weak accession lower the rating of another stronger accession.
After the per-accession judgments, add:
Ranking rationale: <2-5 sentences>
Best starting dataset: <accession and reason>
```

## Failure Cases To Watch

When using this prompt, watch for these failure modes:

- the model copies `provisional_local_assessment` without thinking
- the model treats any keyword match as relevance
- the model ignores assay-family mismatch
- the model overvalues sample count while ignoring structural complexity
- the model fails to separate verified evidence from inference
- the model becomes too conservative and excludes salvageable conditional candidates

If these failures appear, compare outputs against [LLM_JUDGMENT_REGRESSION_EXAMPLES.md](../evals/LLM_JUDGMENT_REGRESSION_EXAMPLES.md).

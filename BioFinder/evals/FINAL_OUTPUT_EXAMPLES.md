# Final Output Examples

This file shows what BioFinder's final user-facing answers should look like after the full workflow:

1. search brief
2. candidate retrieval
3. evidence extraction
4. LLM judgment
5. final recommendation

The goal is consistency and usability, not rigid formatting.

## Example 1: Broad dataset discovery answer

User request:

```text
Find chicken abdominal fat GEO transcriptome datasets for WGCNA.
```

Suggested final answer:

```markdown
## Search Brief

- Biological question: chicken abdominal-fat transcriptome datasets for network analysis
- Organism: Gallus gallus
- Tissue/cell type: abdominal fat / adipose tissue
- Assay: bulk RNA-seq first, microarray second
- Analysis goal: WGCNA
- Assumptions: single-cell was treated as fallback only

## Candidate GEO Datasets

| Rank | GEO | Rating | Organism | Assay | Samples | Best use | Key risk |
|---:|---|---|---|---|---:|---|---|
| 1 | GSE292022 | Medium | Gallus gallus | bulk RNA-seq | 6 | DEG candidate | too small for WGCNA |
| 2 | GSE104042 | Low | Gallus gallus | microarray | 178 | conditional subset follow-up | superseries; mixed tissues |
| 3 | GSE122224 | Low | Gallus gallus | bulk RNA-seq | 6 | fat-metabolism context only | relevance still needs deeper inspection |

## Recommendation

- Best primary dataset: none yet for strong WGCNA
- Best conditional candidate: GSE104042 only if the abdominal-fat subset is isolated and evaluated separately
- Best targeted DEG candidate: GSE292022
- Not recommended as-is: mixed whole-series datasets or loosely related chicken transcriptome studies

## Next Steps

- Inspect the abdominal-fat-related subseries or sample subset inside GSE104042
- Keep GSE292022 as a DEG-oriented abdominal-fat RNA-seq dataset
- If no better GEO cohort is found, broaden the search or check SRA/BioProject for additional abdominal-fat transcriptome datasets
```

Why this works:

- it is decisive
- it does not pretend a strong WGCNA cohort exists if the evidence does not support that
- it still gives the user a practical next move

## Example 2: Known accession inspection answer

User request:

```text
Can GSE292022 be used for WGCNA?
```

Suggested final answer:

```markdown
## Dataset Assessment

- Accession: GSE292022
- Organism: Gallus gallus
- Tissue: abdominal fat
- Assay/platform: bulk RNA-seq / GPL26853
- Samples: 6 total, verified from GSM-level evidence
- Design: clear 3 vs 3 treatment structure
- Data availability: supplementary XLSX files listed
- Main risks: sample count is too small for a strong WGCNA use case

## Verdict

Rating: Medium
Recommended use: targeted differential expression, not primary WGCNA
Avoid using for: network analysis as the main cohort
Next step: keep it as a DEG candidate and continue searching for a larger abdominal-fat cohort
```

Why this works:

- it answers the user's actual question directly
- it preserves a useful dataset instead of over-downgrading it
- it keeps verified facts and interpretation cleanly separated

## Example 3: Validation-oriented answer

User request:

```text
Find a validation GEO dataset for human osteoarthritis cartilage RNA-seq.
```

Suggested final answer:

```markdown
## Search Brief

- Biological question: independent human osteoarthritis cartilage RNA-seq validation cohort
- Organism: Homo sapiens
- Tissue/cell type: cartilage
- Assay: bulk RNA-seq
- Analysis goal: validation

## Candidate GEO Datasets

| Rank | GEO | Rating | Organism | Assay | Samples | Best use | Key risk |
|---:|---|---|---|---|---:|---|---|
| 1 | GSE324566 | Low | Homo sapiens | bulk RNA-seq | 8 | mechanistic perturbation context | chondrocyte treatment model, not ideal cohort validation |

## Recommendation

- Best validation dataset: no strong independent GEO validation cohort confirmed from this pass
- Use with caution: GSE324566 only as a weak or conditional comparison dataset
- Next move: broaden the search for tissue-level OA cohorts or check adjacent public sources
```

Why this works:

- it does not inflate a weak validation candidate into a strong one
- it still gives the user a realistic fallback and next action

## Tone Guidance

Final answers should usually be:

- short
- decisive
- explicit about what is verified
- honest about weak matches
- practical about next steps

Avoid:

- dumping raw metadata
- repeating every candidate detail
- hiding uncertainty
- sounding more confident than the evidence supports

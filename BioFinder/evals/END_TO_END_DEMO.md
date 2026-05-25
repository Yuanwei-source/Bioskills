# End-to-End Demo

This demo shows a complete BioFinder MVP run for one real GEO use case.

## User Request

```text
Find chicken abdominal fat GEO transcriptome datasets for WGCNA.
```

## Step 1: Search Brief

```text
Biological question: chicken abdominal-fat transcriptome datasets for network analysis
Organism: Gallus gallus
Tissue/cell type: abdominal fat / adipose tissue
Disease/condition/phenotype: fat deposition / adipose response context
Assay/data type: transcriptome -> bulk RNA-seq first, microarray second
Analysis goal: WGCNA
Must-have constraints: chicken, abdominal fat relevance
Nice-to-have constraints: larger sample size, clear phenotype structure
Exclusions: single-cell as primary result
```

## Step 2: Candidate Search

Command used:

```text
python scripts/search_geo.py \
  --query "Gallus gallus abdominal fat transcriptome" \
  --query "chicken abdominal fat adipose" \
  --max-results 5 \
  --taxon "Gallus gallus" \
  --single-species-only \
  --include-text "fat" \
  --json
```

Key returned candidates:

- `GSE292022`: abdominal fat, bulk RNA-seq, 6 samples
- `GSE122224`: fat storage/metabolism context, bulk RNA-seq, 6 samples
- `GSE104042`: adipose plus other tissues, microarray, 178 samples, superseries risk
- `GSE104041`: PBMC-focused subseries, not a strong abdominal-fat match
- `GSE144982`: chicken skin-color transcriptome, weak relevance

Immediate search-stage observation:

- the search layer found biologically related candidates
- but the list still included mixed-quality and weakly related accessions
- this confirmed the need for accession-level inspection before recommendation

## Step 3: Evidence Extraction

Command used:

```text
python scripts/inspect_gse.py \
  --accession GSE292022 \
  --accession GSE104042 \
  --analysis-goal WGCNA \
  --json
```

### Evidence Summary: GSE292022

- species: `Gallus gallus`
- tissue: `Abdominal fat`
- assay: `bulk_rna_seq`
- sample count: `6`, verified from GSM-level evidence
- design clarity: `high`
- group signal: clear treatment split, 3 vs 3
- local conclusion: usable for DEG, weak for WGCNA

### Evidence Summary: GSE104042

- species: `Gallus gallus`
- tissues: `Liver`, `Abdominal fat`, `Pectoralis major muscle`, `Peripheral Blood Mononuclear Cells`
- assay: `microarray`
- sample count: `178`, verified from GSM-level evidence
- design clarity: `partial`
- major risks: `superseries`, `multiple_tissues`
- local conclusion: not suitable for whole-series WGCNA as-is

## Step 4: LLM Judgment

Judgment target:

- user wants a strong abdominal-fat WGCNA dataset

Expected LLM interpretation:

- `GSE292022` is strongly relevant in tissue and assay, but too small for WGCNA
- `GSE104042` is large and potentially useful, but only after abdominal-fat subset rescue; the whole accession is too mixed
- no strong WGCNA-ready GEO cohort is confirmed from this pass

## Step 5: Final User-Facing Recommendation

Suggested final answer:

```markdown
## Search Brief

- Biological question: chicken abdominal-fat transcriptome datasets for WGCNA
- Organism: Gallus gallus
- Tissue/cell type: abdominal fat / adipose tissue
- Assay: bulk RNA-seq first, microarray second
- Analysis goal: WGCNA
- Assumptions: single-cell was treated as fallback only

## Candidate GEO Datasets

| Rank | GEO | Rating | Organism | Assay | Samples | Best use | Key risk |
|---:|---|---|---|---|---:|---|---|
| 1 | GSE292022 | Medium | Gallus gallus | bulk RNA-seq | 6 | DEG candidate | too small for WGCNA |
| 2 | GSE104042 | Low | Gallus gallus | microarray | 178 | subset-level follow-up | superseries; mixed tissues |
| 3 | GSE122224 | Low | Gallus gallus | bulk RNA-seq | 6 | exploratory context | deeper relevance check still needed |

## Recommendation

- Best primary dataset: none yet for strong WGCNA
- Best conditional candidate: GSE104042 only if the abdominal-fat subset is isolated and reassessed
- Best targeted DEG candidate: GSE292022
- Not recommended as-is: weakly related or mixed whole-series datasets

## Next Steps

- inspect the abdominal-fat-related subseries or sample subset inside GSE104042
- keep GSE292022 as a DEG-oriented abdominal-fat RNA-seq dataset
- if GEO remains weak for WGCNA, broaden the search or move to SRA/BioProject for additional cohorts
```

## What This Demo Proves

- the MVP can translate a biological question into GEO queries
- the search layer can retrieve real candidate accessions
- the inspection layer can extract usable structured evidence
- the LLM judgment layer can turn that evidence into a realistic recommendation
- the workflow can say “no strong cohort found yet” without collapsing into either noise or false confidence

## MVP Status From This Demo

The GEO-first MVP is operational.

What remains after this demo is mostly polish:

- more final-output examples
- more edge-case regression coverage
- optional subgroup-aware evidence extraction for mixed superseries

---
name: BioFinder
description: Use this skill when the user wants to find, inspect, compare, or judge public bioinformatics datasets for omics analysis. BioFinder is GEO-first, with early v2 support for SRA when raw sequencing context or cross-provider evidence is needed.
metadata:
  short-description: Find and assess public omics datasets with a GEO-first MVP
---

# BioFinder

BioFinder is a public-omics dataset discovery skill.

The current runtime is still GEO-first. It is built to do one workflow well:

1. search GEO for candidate datasets
2. extract structured evidence from promising accessions
3. use the LLM to judge which datasets actually fit the user's goal

v2 adds an early SRA layer for cases where raw sequencing context matters.

Keep the skill lean. Do not turn it into a large local rule engine.

## Runtime Boundary

BioFinder runtime consists of:

- `scripts/biofinder.py`: optional thin orchestration entry point for local workflow runs
- `scripts/evidence_schema.py`: common evidence-shape helpers for provider outputs
- `scripts/search_geo.py`: GEO candidate retrieval and coarse filtering
- `scripts/inspect_gse.py`: accession-level evidence extraction
- `scripts/search_sra.py`: early v2 SRA candidate retrieval and metadata extraction
- `scripts/inspect_sra.py`: early v2 SRA accession inspection and study-level aggregation
- `scripts/search_bioproject.py`: early v2 BioProject candidate retrieval
- `scripts/inspect_bioproject.py`: early v2 BioProject project-level context aggregation
- `references/query_expansion.md`: bounded term-expansion guidance
- `references/clarification_strategy.md`: stable rules for follow-up questions versus default assumptions
- `references/dataset_scoring.md`: LLM-side judgment guide
- `references/llm_judgment_prompt.md`: reusable final-judgment prompt
- `references/sra_provider.md`: when and how to use SRA in the v2 workflow
- `references/bioproject_provider.md`: when and how to use BioProject in the v2 workflow

Development and regression material lives under `evals/` and is not part of normal runtime use.
Provider expansion planning lives in [v2_roadmap.md](references/v2_roadmap.md).

## Scope

The current GEO-first path supports:

- GEO Series (`GSE`), Samples (`GSM`), and Platforms (`GPL`)
- bulk expression datasets from microarray or RNA-seq
- GEO metadata, sample design, platform, supplementary files, and raw-data clues
- LLM judgment for DEG, WGCNA, validation, and exploratory ML use
- cautious single-cell discovery only when the user explicitly accepts exploratory results

The current skill does not try to fully support:

- full SRA, BioProject, ArrayExpress, TCGA, ENA, or dbGaP workflows
- full FASTQ-processing pipelines
- complete automated phenotype rescue from messy metadata
- heavy local ranking logic that should live in the LLM judgment step

Current v2 SRA support is intentionally narrow:

- search SRA for raw-sequencing candidates
- inspect one SRA accession and aggregate study-level metadata
- emit unified evidence for later LLM judgment

Read [sra_provider.md](references/sra_provider.md) before using SRA in the workflow.

Current v2 BioProject support is also intentionally narrow:

- search BioProject for project-level context
- inspect one BioProject accession and aggregate linked accession clues
- use BioProject mainly as a context provider, not as the first search surface

Read [bioproject_provider.md](references/bioproject_provider.md) before using BioProject in the workflow.

## Core Workflow

For most requests, follow this order:

1. Convert the user's request into a compact search brief:

```text
Biological question:
Organism:
Tissue/cell type:
Disease/condition/phenotype:
Assay/data type:
Analysis goal:
Must-have constraints:
Nice-to-have constraints:
Exclusions:
```

2. Normalize the request just enough to search well:
- if the request is broad or underspecified, first map it into the structured parse in [request_parse_template.md](references/request_parse_template.md)
- use LLM-driven term expansion with tight bounds
- keep a small local rule layer for organism naming, assay-family boundaries, and bulk versus single-cell defaults
- read [query_expansion.md](references/query_expansion.md) when the request is broad or ambiguous

3. Choose the provider path:
- default to GEO for expression-study discovery
- add SRA when the user asks for raw sequencing data, when GEO is weak, or when cross-provider evidence is needed
- add BioProject when you already have a `PRJ...` accession or when GEO/SRA results need project-level context
- read [sra_provider.md](references/sra_provider.md) before switching to or adding SRA
- read [bioproject_provider.md](references/bioproject_provider.md) before switching to or adding BioProject

4. Search for candidates:
- use `scripts/search_geo.py` for GEO
- use `scripts/search_sra.py` for early v2 SRA
- use `scripts/search_bioproject.py` for early v2 BioProject
- prefer multiple narrow-to-broad query variants
- keep fallback matches clearly labeled
- treat search hits as candidates, not as final recommendations

5. Extract evidence from strong candidates:
- use `scripts/inspect_gse.py` for GEO accessions
- use `scripts/inspect_sra.py` for SRA accessions
- use `scripts/inspect_bioproject.py` for BioProject accessions
- verify sample counts when possible
- surface tissue, organism, assay family, grouping signals, run-level availability, project-level links, and obvious structural risks
- treat local assessments as provisional evidence, not as the final answer

6. Judge the candidates with the LLM:
- use [dataset_scoring.md](references/dataset_scoring.md) for judgment criteria
- use [llm_judgment_prompt.md](references/llm_judgment_prompt.md) for the final prompt shape
- keep verified evidence separate from inference

7. Recommend the next action:
- best primary dataset
- best conditional or validation dataset
- main risks
- next download or analysis entry point

## Clarification Rules

Ask for clarification only when the request is too ambiguous to search responsibly.

Read [clarification_strategy.md](references/clarification_strategy.md) when the request is broad, underspecified, or mixes several possible search spaces.
Read [request_parse_template.md](references/request_parse_template.md) when you need a stable structured parse before deciding whether to clarify or proceed.

Usually ask when two or more of these are missing:

- organism
- tissue or cell type
- disease, condition, or phenotype
- assay or data type
- analysis goal

Do not stop for clarification if the user explicitly wants a broad exploratory search. State your assumptions and proceed.

For vague assay wording such as `transcriptome`, default to:

- bulk RNA-seq
- microarray
- single-cell RNA-seq only as exploratory fallback

## Local Versus LLM Roles

Use local scripts for:

- candidate retrieval
- metadata extraction
- phenotype-table flattening
- study-level SRA aggregation
- obvious risk flags such as multi-species, mixed tissue, or superseries structure

Use the LLM for:

- deciding which candidate best fits the user's real biological target
- deciding whether a borderline dataset is still worth keeping
- final `High` / `Medium` / `Low` / `Exclude` judgments
- ranking similar candidates and explaining tradeoffs

## Output Shape

For dataset discovery, aim to produce:

- the search brief and assumptions
- a short ranked candidate table
- a decisive recommendation
- next-step guidance

For accession inspection, aim to produce:

- accession summary
- verified evidence
- main risks
- goal-specific judgment
- next-step recommendation

When confidence is limited, say what is verified and what is inferred.

## Quality Bar

Before finalizing:

- verify accessions and provider facts
- do not overstate sample counts or labels when metadata is incomplete
- keep fallback or weak matches clearly separated from strong matches
- avoid copying local provisional assessments as if they were final judgments
- keep the response tied to the user's requested analysis goal

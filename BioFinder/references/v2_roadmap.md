# BioFinder v2 Roadmap

Use this file as the implementation roadmap for BioFinder after the GEO-first
MVP.

The guiding constraint remains the same:

- keep local scripts focused on retrieval and evidence extraction
- keep final ranking and borderline judgment in the LLM
- expand providers only when they fit the same evidence-oriented workflow

## v2 Objective

Turn BioFinder from a GEO-first discovery skill into a public-omics discovery
assistant that can:

- search GEO for expression-centric studies
- follow links into SRA and BioProject when raw sequencing context matters
- emit a common evidence shape across providers
- package evidence for final LLM judgment without growing a large local rule engine

## v2 Must

### 1. Common Evidence Schema

Every provider should map into the same top-level structure:

- `source`
- `accession`
- `source_type`
- `title`
- `organisms`
- `tissues`
- `assay_family`
- `sample_count`
- `run_count`
- `group_signals`
- `data_availability`
- `linked_accessions`
- `risks`
- `source_links`
- `local_notes`

### 2. Unified Entry Point

Add `scripts/biofinder.py` as a thin orchestrator that:

- accepts the user query and search constraints
- runs provider search
- inspects top candidates
- emits normalized evidence
- prepares material for LLM judgment

### 3. SRA Discovery Layer

Add lightweight SRA support for:

- searching studies, experiments, or runs
- extracting library strategy and layout
- surfacing BioProject and BioSample links
- signaling whether raw sequencing data is likely available

Do not pull FASTQ or run SRA Toolkit workflows in core runtime.

### 4. BioProject Context Layer

Add BioProject context extraction for:

- project title and summary
- associated SRA studies
- project-level scope clues
- dataset grouping and complexity clues

## v2 Should

### 1. Better Cross-Provider Linking

Strengthen:

- GEO to SRA links
- GEO to BioProject links
- SRA to BioProject links
- duplicate suppression across providers

### 2. Better Data Availability Labels

Prefer explicit labels over guesses:

- `processed_matrix_available`
- `raw_counts_available`
- `normalized_expression_available`
- `sra_runs_available`
- `fastq_candidate`
- `needs_manual_check`

### 3. Better Exports

Support:

- `--json`
- `--tsv`
- `--evidence-json`
- `--llm-prompt`

## v2 Later

These should wait until GEO plus SRA plus BioProject are stable:

- ArrayExpress or BioStudies
- ENA
- GSA
- automated FASTQ download
- analysis pipelines beyond discovery and evidence extraction

## Delivery Order

1. common evidence schema
2. `biofinder.py` orchestrator
3. SRA search and inspect scripts
4. BioProject inspection
5. cross-provider linking and exports

# Request Parse Template

Use this template at the start of a BioFinder run.

The purpose is to convert the user's natural-language request into a consistent structured parse before:

- deciding whether to clarify
- expanding search terms
- running GEO search
- judging candidate fit

This template is especially important for broad or ambiguous requests such as:

- `find diabetes islet GEO data`
- `find cartilage transcriptome datasets`
- `look for validation cohorts`

## Output Template

Fill this structure as completely as possible from the user request:

```text
User request:

Parsed intent:
- task_type:
- provider_scope:
- user_mode:

Core fields:
- organism:
- tissue_or_cell_type:
- disease_or_condition:
- phenotype_or_trait:
- assay_or_data_type:
- analysis_goal:

Field status:
- organism_status:
- tissue_status:
- disease_status:
- assay_status:
- analysis_goal_status:

Ambiguity flags:
- ambiguous_species:
- ambiguous_tissue_boundary:
- ambiguous_assay_family:
- ambiguous_goal:
- mixed_search_space_risk:

Default assumptions:
- assumption_1:
- assumption_2:
- assumption_3:

Clarification decision:
- ask_follow_up:
- reason:
- follow_up_question:

If proceeding now:
- search_strategy_mode:
- primary_search_target:
- fallback_policy:
```

## Field Definitions

### `task_type`

Choose the best fit:

- `find_datasets_from_scratch`
- `inspect_known_accession`
- `find_validation_dataset`
- `compare_known_accessions`
- `prepare_download_or_analysis_entry`
- `mixed_request`

### `provider_scope`

For the current MVP, this is usually:

- `geo_first`

If the user explicitly asks for another source or if GEO is clearly insufficient, note that in assumptions, but do not broaden provider scope automatically during the GEO-first MVP.

### `user_mode`

Choose one:

- `focused`
- `broad_exploratory`
- `mixed`

Clues for `broad_exploratory`:

- `as many as possible`
- `broad search`
- `exploratory`
- `先尽可能多找`

## Field Status Labels

For each core field, use one of:

- `explicit`
- `implicit`
- `missing`
- `dangerously_ambiguous`

Use `dangerously_ambiguous` when a wrong assumption would strongly distort the search space.

Examples:

- `islet` without knowing whether the user means isolated islet or whole pancreas
- `transcriptome` when the user may want bulk, microarray, or single-cell
- `cartilage` when the evidence might drift into `chondrocyte`

## Clarification Decision Rule

### Ask first when:

- two or more core fields are `missing` or `dangerously_ambiguous`
- one field is present but would split the search into very different biological targets
- species and assay are both unresolved

### Proceed with assumptions when:

- only one major field is missing
- the user explicitly wants a broad exploratory pass
- the missing field can be handled with a stable default

## Stable Defaults

Use these defaults when proceeding without clarification:

- `transcriptome` -> bulk RNA-seq first, microarray second, single-cell fallback only
- missing analysis goal -> general discovery
- missing tissue -> do not invent tissue; keep search broad and label tissue unresolved
- missing species in exploratory mode -> keep species broad, but separate candidates by species in output

## Recommended Follow-Up Shape

If `ask_follow_up = yes`, ask one short question when possible.

Example:

```text
Before I search, do you want human or mouse data, and are you looking for bulk expression data or single-cell results?
```

Avoid stacked interrogations unless absolutely necessary.

## Example 1

User request:

```text
Find chicken abdominal fat GEO transcriptome datasets for WGCNA.
```

Suggested parse:

```text
User request:
Find chicken abdominal fat GEO transcriptome datasets for WGCNA.

Parsed intent:
- task_type: find_datasets_from_scratch
- provider_scope: geo_first
- user_mode: focused

Core fields:
- organism: chicken / Gallus gallus
- tissue_or_cell_type: abdominal fat / adipose tissue
- disease_or_condition: not central; fat deposition context
- phenotype_or_trait: fat deposition / adipose response context
- assay_or_data_type: transcriptome
- analysis_goal: WGCNA

Field status:
- organism_status: explicit
- tissue_status: explicit
- disease_status: implicit
- assay_status: dangerously_ambiguous
- analysis_goal_status: explicit

Ambiguity flags:
- ambiguous_species: no
- ambiguous_tissue_boundary: no
- ambiguous_assay_family: yes
- ambiguous_goal: no
- mixed_search_space_risk: low

Default assumptions:
- assumption_1: transcriptome is interpreted as bulk RNA-seq first
- assumption_2: microarray is secondary fallback
- assumption_3: single-cell is fallback only

Clarification decision:
- ask_follow_up: no
- reason: only one major ambiguity, with a stable default
- follow_up_question: none

If proceeding now:
- search_strategy_mode: focused
- primary_search_target: Gallus gallus abdominal fat bulk transcriptome datasets
- fallback_policy: allow microarray fallback; keep single-cell out of the main shortlist
```

## Example 2

User request:

```text
Find diabetes islet GEO data.
```

Suggested parse:

```text
User request:
Find diabetes islet GEO data.

Parsed intent:
- task_type: find_datasets_from_scratch
- provider_scope: geo_first
- user_mode: focused

Core fields:
- organism: unresolved
- tissue_or_cell_type: islet
- disease_or_condition: diabetes
- phenotype_or_trait: unresolved
- assay_or_data_type: unresolved
- analysis_goal: unresolved

Field status:
- organism_status: missing
- tissue_status: dangerously_ambiguous
- disease_status: explicit
- assay_status: missing
- analysis_goal_status: missing

Ambiguity flags:
- ambiguous_species: yes
- ambiguous_tissue_boundary: yes
- ambiguous_assay_family: yes
- ambiguous_goal: yes
- mixed_search_space_risk: high

Default assumptions:
- assumption_1: none strong enough to avoid clarification
- assumption_2: if forced to proceed, separate human and mouse results
- assumption_3: do not merge islet with pancreas

Clarification decision:
- ask_follow_up: yes
- reason: several core fields are missing or dangerously ambiguous
- follow_up_question: Before I search, do you want human or mouse data, and are you looking for bulk expression data or single-cell results?

If proceeding now:
- search_strategy_mode: broad exploratory only if the user asks for that
- primary_search_target: diabetes plus pancreatic islet
- fallback_policy: keep species and assay families separated in the output
```

## Example 3

User request:

```text
Find as many osteoarthritis cartilage GEO datasets as possible.
```

Suggested parse:

```text
User request:
Find as many osteoarthritis cartilage GEO datasets as possible.

Parsed intent:
- task_type: find_datasets_from_scratch
- provider_scope: geo_first
- user_mode: broad_exploratory

Core fields:
- organism: unresolved
- tissue_or_cell_type: cartilage
- disease_or_condition: osteoarthritis
- phenotype_or_trait: unresolved
- assay_or_data_type: unresolved
- analysis_goal: general discovery

Field status:
- organism_status: missing
- tissue_status: explicit
- disease_status: explicit
- assay_status: missing
- analysis_goal_status: implicit

Ambiguity flags:
- ambiguous_species: yes
- ambiguous_tissue_boundary: moderate
- ambiguous_assay_family: yes
- ambiguous_goal: no
- mixed_search_space_risk: moderate

Default assumptions:
- assumption_1: broad exploratory mode means no clarification stop
- assumption_2: bulk RNA-seq first, microarray second, single-cell fallback
- assumption_3: species-separated presentation if several species appear

Clarification decision:
- ask_follow_up: no
- reason: the user explicitly asked for a broad pass
- follow_up_question: none

If proceeding now:
- search_strategy_mode: broad exploratory
- primary_search_target: osteoarthritis cartilage GEO datasets across likely assay families
- fallback_policy: keep single-cell and weak tissue matches labeled as fallback
```

## Why This Template Matters

This parse gives the workflow a stable bridge between:

- user language
- clarification decisions
- query expansion
- search execution
- LLM judgment

Without it, ambiguous requests tend to drift.

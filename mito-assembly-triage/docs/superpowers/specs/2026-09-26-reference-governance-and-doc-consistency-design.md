# Reference governance and documentation consistency

## Goal

Make the diagnostic documentation internally consistent and keep its claims aligned
with the current command-line interfaces. The change must preserve the existing
scientific safeguards and must not add analysis features, loosen evidence thresholds,
or alter user data.

## Canonical ownership

`reference-policy.md` and `scripts/reference_registry.py` own the L1--L5
reference-purpose matrix. Other documents may describe how to choose a level, but
must link to that matrix rather than restating a divergent version.

`evidence-standard.md` owns confidence semantics. It will distinguish an
evidence-supported candidate claim from a claim that requires unavailable sample
reads for direct verification. Missing reads therefore limits only the latter.

`conclusion-report.md` describes the report the current `case-report` command can
produce: claim, status, confidence, optional reads support, and event links.
Richer fields remain recommended narrative content until the schema and CLI support
them. Scientific status, expert review, and permission to adopt/publish a result are
kept conceptually separate.

## Documentation changes

- Replace the generic external-query authorization sentence with the actual boundary:
  public-reference download is allowed and recorded; uploading a sample or sharing
  data requires explicit authorization.
- State that persistent case/lesson storage requires user authorization, while a
  task-local report may be generated in the working directory.
- Replace all default `--table 5` examples with a required, taxon-confirmed table
  placeholder; label insect-specific expectations as an insect profile.
- Make `START_HERE.md` a router and conditional first-pass checklist. Remove its
  duplicated reference-level policy and clarify links that previously used an
  unqualified section number.
- Correct wording that treats low coverage plus soft clips as a normal state. It is
  an explanation to test, not a conclusion.
- Replace non-retrievable literature placeholders in gene-order guidance with a
  request to supply a taxon-specific source before making a biological claim.

## Verification

Run the repository test suite because no behavior is intended to change. Use
repository-wide searches to confirm there is one reference-purpose matrix, no
unqualified default table-5 command remains, and report promises match the CLI.
Review the rendered Markdown links and heading structure.

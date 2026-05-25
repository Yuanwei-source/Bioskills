# SRA Real Tests

This file records early real-world validation of the BioFinder v2 SRA provider.

## Test 1: search_sra.py

Command shape:

```text
python scripts/search_sra.py --query "Gallus gallus abdominal fat RNA-Seq" --max-results 3 --json
```

Observed outcome after parser fixes:

- returned a clean study accession: `SRP691785`
- extracted linked accessions:
  - `SRX32944937`
  - `SRS28784844`
  - `SRR38097878`
  - additional related `SRX/SRS/SRR` accessions after aggregation
- extracted:
  - `Gallus gallus`
  - `RNA-Seq`
  - `PAIRED`
  - `Illumina NovaSeq X Plus`
  - `PRJNA1453495`

Issues discovered and fixed:

- Biopython SRA `esummary` shape was not stable enough to treat as a plain dict
- submitter accession such as `SRA2374050` was initially being treated as the main accession
- repeated biosample-level hits were initially returned as separate candidates

Current status:

- main accession selection is now study-first
- repeated records are deduplicated into one study candidate

## Test 2: inspect_sra.py

Command shape:

```text
python scripts/inspect_sra.py --accession SRP691785 --analysis-goal differential_expression --json
```

Observed outcome:

- accession: `SRP691785`
- organism: `Gallus gallus`
- assay family: `RNA-Seq`
- library layout: `PAIRED`
- platform: `Illumina NovaSeq X Plus`
- sample count: `6`
- run count: `6`
- bioproject: `PRJNA1453495`
- provisional local assessment: `high`

Linked accessions aggregated:

- study: `SRP691785`
- experiments: `SRX32944932` to `SRX32944937`
- samples: `SRS28784842` to `SRS28784847`
- runs: `SRR38097878` to `SRR38097883`

Current limitations observed:

- title was sparse and only reflected one sample name: `LAF-3`
- tissue was not reliably recoverable from SRA metadata alone
- final recommendation still needs LLM interpretation, especially when SRA lacks phenotype detail

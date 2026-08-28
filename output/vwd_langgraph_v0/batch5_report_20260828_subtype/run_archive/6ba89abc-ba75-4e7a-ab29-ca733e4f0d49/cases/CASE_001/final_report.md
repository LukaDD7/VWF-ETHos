# VWD Research Report — CASE_001

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

Retrospective review integrated 32 FHIR evidence resources, identified 1 structured evidence issue(s), and found 6 missing evidence item(s). Current subtype tendencies are type_2_vwd (score 2.0, moderate), type_2A (score 0.8, low). All second-level tests were unavailable; no values were imputed. Expert review is required.

## Candidate subtypes

type_2_candidate

## Supporting evidence

- CASE_001
- Observation/variant-normalized-fbe5f6f4a92219719ad8
- Observation/vwf-domain-1500-Ala-Val
- Observation/ocravat-gnomad-CASE_001-chr12-6018919-G-A
- Observation/ocravat-clinvar-CASE_001-chr12-6018919-G-A
- Observation/ocravat-revel-CASE_001-chr12-6018919-G-A
- Observation/ocravat-cadd-CASE_001-chr12-6018919-G-A
- Observation/ocravat-alphamissense-CASE_001-chr12-6018919-G-A
- Observation/pmc-excerpt-pubmed-21289515-0
- Observation/pmc-excerpt-pubmed-21289515-1
- Observation/pmc-excerpt-pubmed-21289515-2
- Observation/pmc-excerpt-pubmed-26986123-0
- Observation/pmc-excerpt-pubmed-26986123-1
- Observation/pmc-excerpt-pubmed-26986123-2
- Observation/lit-vwf-ag-pubmed-21289515-0
- Observation/lit-vwf-ag-pubmed-21289515-1
- Observation/lit-vwf-rco-pubmed-21289515-0
- Observation/lit-vwf-rco-pubmed-21289515-1
- Observation/lit-fviii-c-pubmed-21289515-0
- Observation/lit-fviii-c-pubmed-21289515-1
- Observation/lit-ripa-pubmed-21289515-0
- Observation/lit-vwf-multimer-pubmed-21289515-0
- Observation/lit-vwf-multimer-pubmed-21289515-1
- Observation/lit-structural-change-pubmed-21289515-0
- Observation/lit-abo-effect-pubmed-21289515-0
- Observation/lit-abo-effect-pubmed-21289515-1
- Observation/lit-vwf-ag-pubmed-26986123-0
- Observation/lit-vwf-multimer-pubmed-26986123-0
- Observation/lit-vwf-multimer-pubmed-26986123-1
- DocumentReference/pubmed-33497541
- DocumentReference/pubmed-21289515
- DocumentReference/pubmed-33405380
- DocumentReference/pubmed-26986123

## Missing information

- RIPA
- VWF_MULTIMER
- VWF_CB
- VWF_FVIIIB
- VWF_PP
- DDAVP_0_1_4H

## Limitations

- gnomad_graphql returned error; this is not benign evidence.
- clingen_erepo returned not_found; this is not benign evidence.
- ClinGen/ACMG type-2 rule provider is not enabled in V0 offline mode.
- Units, reference ranges, assay methods, and collection times are absent.
- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- The evidence cannot distinguish type 2 subtypes or exclude qualitative VWF defects.

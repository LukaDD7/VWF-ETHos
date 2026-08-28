# VWD Research Report — CASE_011

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

Retrospective review integrated 17 FHIR evidence resources, identified 2 structured evidence issue(s), and found 7 missing evidence item(s). Current subtype tendencies are type_2_vwd (score 2.0, moderate), pathogenic_variant (score 1.0, low). All second-level tests were unavailable; no values were imputed. Expert review is required.

## Candidate subtypes

type_2_candidate

## Supporting evidence

- CASE_011
- Observation/variant-normalized-7ce70ca2585f05efc650
- Observation/vwf-domain-1205-R-H
- Observation/ocravat-clinvar-CASE_011-chr12-6021960-C-T
- Observation/ocravat-revel-CASE_011-chr12-6021960-C-T
- Observation/ocravat-cadd-CASE_011-chr12-6021960-C-T
- Observation/ocravat-alphamissense-CASE_011-chr12-6021960-C-T
- Observation/pmc-excerpt-pubmed-38996914-0
- Observation/pmc-excerpt-pubmed-38996914-1
- Observation/pmc-excerpt-pubmed-38996914-2
- Observation/lit-vwf-ag-pubmed-38996914-0
- Observation/lit-secretion-pubmed-38996914-0
- Observation/lit-ddavp-response-pubmed-38996914-0
- DocumentReference/pubmed-39304225
- DocumentReference/pubmed-38996914
- DocumentReference/pubmed-25690668
- DocumentReference/pubmed-16420565
- DocumentReference/pubmed-23490306

## Missing information

- RIPA
- VWF_MULTIMER
- VWF_CB
- VWF_FVIIIB
- VWF_PP
- DDAVP_0_1_4H
- population_frequency

## Limitations

- gnomad_graphql returned error; this is not benign evidence.
- clingen_erepo returned not_found; this is not benign evidence.
- ClinGen/ACMG type-2 rule provider is not enabled in V0 offline mode.
- Units, reference ranges, assay methods, and collection times are absent.
- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- REVEL, CADD, and AlphaMissense capture different signals; disagreement must remain explicit.
- The evidence cannot distinguish type 2 subtypes or exclude qualitative VWF defects.
- ACMG hint: PP3_insufficient_alone: only one predictor currently supports deleterious effect

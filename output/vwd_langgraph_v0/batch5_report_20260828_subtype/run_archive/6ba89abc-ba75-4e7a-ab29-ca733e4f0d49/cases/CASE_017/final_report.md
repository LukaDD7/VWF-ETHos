# VWD Research Report — CASE_017

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

Retrospective review integrated 22 FHIR evidence resources, identified 2 structured evidence issue(s), and found 7 missing evidence item(s). Current subtype tendencies are type_2B (score 2.6, moderate), type_2_vwd (score 2.0, moderate), pathogenic_variant (score 2.0, moderate). All second-level tests were unavailable; no values were imputed. Expert review is required.

## Candidate subtypes

type_2_candidate, platelet_type_vwd_candidate

## Supporting evidence

- CASE_017
- Observation/variant-normalized-6ffd0dc74cde60a19e2f
- Observation/vwf-domain-1316-Val-Met
- Observation/ocravat-clinvar-CASE_017-chr12-6019472-C-T
- Observation/ocravat-revel-CASE_017-chr12-6019472-C-T
- Observation/ocravat-cadd-CASE_017-chr12-6019472-C-T
- Observation/ocravat-alphamissense-CASE_017-chr12-6019472-C-T
- Observation/clingen-0-d6ce33b9eb092e36a390
- Observation/pmc-excerpt-pubmed-32618441-0
- Observation/pmc-excerpt-pubmed-32618441-1
- Observation/pmc-excerpt-pubmed-32618441-2
- Observation/pmc-excerpt-pubmed-39820471-0
- Observation/pmc-excerpt-pubmed-39820471-1
- Observation/pmc-excerpt-pubmed-39820471-2
- Observation/lit-vwf-multimer-pubmed-32618441-0
- Observation/lit-ddavp-response-pubmed-32618441-0
- Observation/lit-vwf-multimer-pubmed-39820471-0
- Observation/lit-vwf-multimer-pubmed-39820471-1
- DocumentReference/pubmed-32618441
- DocumentReference/pubmed-11475150
- DocumentReference/pubmed-34942660
- DocumentReference/pubmed-39820471
- DocumentReference/pubmed-31939074

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
- ClinGen/ACMG type-2 rule provider is not enabled in V0 offline mode.
- Units, reference ranges, assay methods, and collection times are absent.
- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- REVEL, CADD, and AlphaMissense capture different signals; disagreement must remain explicit.
- The evidence cannot distinguish type 2 subtypes or exclude qualitative VWF defects.
- ACMG hint: PP3_supporting_candidate: at least two computational predictors indicate deleterious missense effect

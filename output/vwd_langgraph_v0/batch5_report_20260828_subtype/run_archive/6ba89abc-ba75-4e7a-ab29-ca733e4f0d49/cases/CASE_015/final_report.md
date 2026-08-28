# VWD Research Report — CASE_015

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

Retrospective review integrated 14 FHIR evidence resources, identified 1 structured evidence issue(s), and found 6 missing evidence item(s). Current subtype tendencies are type_2_vwd (score 2.0, moderate), pathogenic_variant (score 1.0, low), type_2A (score 0.8, low). All second-level tests were unavailable; no values were imputed. Expert review is required.

## Candidate subtypes

type_2_candidate

## Supporting evidence

- CASE_015
- Observation/variant-normalized-efbc8a7008850290f269
- Observation/vwf-domain-1506-Ser-Leu
- Observation/gnomad-exome-12-6018901-G-A
- Observation/gnomad-genome-12-6018901-G-A
- Observation/ocravat-clinvar-CASE_015-chr12-6018901-G-A
- Observation/ocravat-revel-CASE_015-chr12-6018901-G-A
- Observation/ocravat-cadd-CASE_015-chr12-6018901-G-A
- Observation/ocravat-alphamissense-CASE_015-chr12-6018901-G-A
- Observation/pmc-excerpt-pubmed-36226571-0
- Observation/pmc-excerpt-pubmed-36226571-1
- Observation/pmc-excerpt-pubmed-36226571-2
- DocumentReference/pubmed-36226571
- DocumentReference/pubmed-31939074
- DocumentReference/pubmed-28536718

## Missing information

- RIPA
- VWF_MULTIMER
- VWF_CB
- VWF_FVIIIB
- VWF_PP
- DDAVP_0_1_4H

## Limitations

- clingen_erepo returned not_found; this is not benign evidence.
- literature_phenotype_extractor returned not_found; this is not benign evidence.
- ClinGen/ACMG type-2 rule provider is not enabled in V0 offline mode.
- Units, reference ranges, assay methods, and collection times are absent.
- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- The evidence cannot distinguish type 2 subtypes or exclude qualitative VWF defects.
- ACMG hint: PM2_supporting_candidate: allele is absent or very rare in available gnomAD annotation
- ACMG hint: PP3_supporting_candidate: at least two computational predictors indicate deleterious missense effect

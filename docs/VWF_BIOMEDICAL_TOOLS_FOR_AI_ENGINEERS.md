# VWF Biomedical Tools for AI Engineers

## Mental model

Clinical variant interpretation is not a single classifier. It is a multi-source evidence-reconciliation problem:

```text
patient phenotype
+ VWF laboratory phenotype
+ variant identity
+ population frequency
+ curated clinical databases
+ computational predictions
+ case reports / functional studies
+ expert rules
→ auditable clinical hypothesis
```

The Agent's job is not to "vote" freely. It must:

1. normalize the variant;
2. call controlled tools;
3. preserve provenance;
4. distinguish contradiction from missing evidence;
5. abstain or request expert review when evidence is insufficient.

## Tool map

| Tool | What it answers | Typical output | Key limitation |
|---|---|---|---|
| Ensembl Variant Recoder | What is the canonical genomic identity of this HGVS variant? | GRCh38 SPDI, HGVSg, rsID, chrom-pos-ref-alt | Requires a valid transcript/HGVS; normalization failure blocks coordinate tools |
| gnomAD | How common is this allele in population sequencing? | Exome/genome AC, AN, AF | A common allele is unlikely to be a highly penetrant dominant disease allele; rare is not automatically pathogenic |
| ClinVar | What have clinical submitters classified this variant as? | Clinical significance, review status, condition, submission metadata | Can contain conflicts and variable review quality |
| ClinGen Evidence Repository | What does a VWD expert panel conclude? | Expert classification, evidence codes, condition | Not every VWF variant is curated; not-found is not benign |
| REVEL | Is this missense substitution likely deleterious? | 0–1 score | Missense-only; calibrated score, not diagnosis |
| CADD | How deleterious is the observed substitution? | PHRED score | High score is supporting evidence, not proof |
| AlphaMissense | What is the model's pathogenicity probability? | 0–1 probability | Can be ambiguous; not a clinical classifier |
| SpliceAI | Does this variant affect splicing? | Delta scores | Only relevant when splicing is plausible; no score is not no effect |
| PubMed | What case reports, phenotypes, and functional studies exist? | PMID, title, journal, year, DOI, abstract | Absence of an exact article is missing evidence, not benign evidence |
| HGMD | Are there published disease variants? | DM/classification metadata | Requires institutional license; no public API in this deployment |
| Local score snapshots | Can the system run offline? | versioned local lookup | Only a fallback; local snapshots are smaller than the full live databases |

## Evidence semantics

The system separates four states:

```text
success     Tool ran and returned evidence.
not_found   Tool ran, but this exact query had no result.
error       Tool failed technically or dependency failed.
disabled    License, configuration, or policy prevented use.
```

These are not interchangeable:

- `not_found` is not benign.
- `disabled` is not benign.
- `error` is not benign.
- absence of a tool result is missing evidence.

## Conflict classes currently implemented

1. **Classification vs population frequency**
   - Example: ClinVar says pathogenic, but gnomAD AF is 2%.
   - Meaning: strong contradiction for a rare dominant disease model.

2. **Predictor disagreement**
   - Example: REVEL points pathogenic, CADD points benign.
   - Meaning: different models capture different signals; do not average away the disagreement.

3. **Database vs database**
   - Example: ClinVar VUS, ClinGen pathogenic.
   - Meaning: expert review is required.

4. **Missing discriminator**
   - Example: Act/Ag ratio suggests possible type 2, but RIPA and multimer studies are unavailable.
   - Meaning: the system must abstain from subtype discrimination.

5. **Insufficient evidence**
   - No hard contradiction, but required sources are absent.
   - Meaning: uncertainty remains, but it is not the same as conflicting evidence.

## Why FHIR matters

FHIR turns clinical workflow into typed, auditable objects:

- `ServiceRequest` = ordered test
- `Observation` = real result
- `Task` = unavailable/cancelled/deferred action
- `DiagnosticReport` = phase-level summary
- `DocumentReference` = literature or model artifact
- `Provenance` = who/what produced the evidence and how
- `OperationOutcome` = explicit failure or policy state

This prevents the common LLM failure mode of describing an unavailable test as if it had been performed.

# VWD LangGraph Batch-5 Interpretation — 2026-08-28

## Executive conclusion

The current workflow now matches the **skeleton** of a physician's variant-review workflow:

```text
First Level labs
→ genetic variant normalization
→ VWF domain annotation
→ gnomAD
→ ClinVar / ClinGen
→ REVEL / CADD / AlphaMissense
→ PubMed abstract and PMC full-text retrieval
→ structured full-text search
→ conflict / missing-evidence analysis
→ constrained second-level test selection
→ final abstaining report with provenance
```

It does **not yet fully reproduce** the physician's final integrative interpretation. The main remaining gap is not tool invocation; it is converting full-text and database evidence into a patient-specific mechanistic explanation comparable to the P1413L example.

## Capability assessment

| Physician step | Current status | Notes |
|---|---|---|
| Variant location / domain | Implemented | Mature-protein domain map is approximate and explicitly marked display-only |
| Population frequency | Implemented | Official gnomAD API plus OpenCRAVAT aggregate fallback |
| ClinVar | Implemented | Classification is captured; review-status detail still needs fuller extraction |
| ClinGen | Implemented | Exact variant lookup works; many variants are not curated |
| HGMD | Interface only | Requires institutional license; remains explicitly disabled |
| REVEL / CADD / AlphaMissense | Implemented | Retrieved through OpenCRAVAT |
| SpliceAI | Implemented when score exists | Missense variants often have no SpliceAI result |
| PubMed abstract | Implemented | ESearch + EFetch |
| PMC full text | Implemented | PMID→PMC mapping + full-text XML retrieval |
| Full-text code-as-search | Implemented | Retrieves bounded excerpts around RIPA, multimer, CB, FVIII, DDAVP, ABO, etc. |
| Literature phenotype extraction | Initial implementation | Extracts lab/qualitative statements but does not yet resolve whether they refer to this exact variant |
| Functional / recombinant experiment interpretation | Partial | Full-text excerpts can support this, but structured experiment extraction is incomplete |
| ABO / clearance context | Missing input | ABO genotype, VWF clearance, and DDAVP 0/1/4h data are not present in the workbook |
| ClinGen/ACMG calculator | Initial hints only | PM2/PP3/BS3-style hints are generated; full VWF VCEP rules are not encoded |
| Final subtype diagnosis | Intentionally abstains | All second-level discriminators are unavailable in this retrospective environment |

## Batch-5 result table

Full machine-readable summary:

[batch_summary.md](../output/vwd_langgraph_v0/batch5_report_20260828_final/batch_summary.md)

| Case | Variant / domain | First-level signal | Population / classification / prediction | Selected second-level tests | Main interpretation |
|---|---|---|---|---|---|
| CASE_001 | p.Ala1500Val; VWFA2; missense | Act/Ag 0.605 | Exome AF 4.31e-5; Genome AF 9.85e-5; ClinVar VUS; REVEL 0.201; CADD 9.552; AlphaMissense 0.106 | VWF_MULTIMER, VWF_CB | Qualitative defect cannot be excluded, but evidence is not strongly pathogenic |
| CASE_011 | p.R1205H; TIL4–VWFA1 linker; missense | Act/Ag 0.376 | Exome AF 2.74e-6; ClinVar Pathogenic; REVEL 0.227; CADD 22.7; AlphaMissense 0.0857 | VWF_MULTIMER, VWF_CB, VWF_FVIIIB | Database classification conflicts with weak computational support; second-level phenotype remains essential |
| CASE_015 | p.Ser1506Leu; VWFA2; missense | Act/Ag 0.299 | Exome AF 2.32e-5; ClinVar Pathogenic/Likely pathogenic; REVEL 0.796; CADD 26.7; AlphaMissense 0.812 | VWF_MULTIMER, VWF_CB | Strongest computational/database signal in this batch, but subtype still requires phenotypic discriminators |
| CASE_017 | p.Val1316Met; A1; missense | Act/Ag 0.470 | ClinVar Pathogenic; ClinGen Pathogenic; REVEL 0.740; CADD 24.7; AlphaMissense 0.4547 | RIPA, VWF_MULTIMER | Most clear database-supported qualitative phenotype; RIPA prioritized because type 2B is plausible |
| CASE_021 | p.A1461D in VWFA1–VWFA2 linker + benign p.D2449N in VWFC2 | Act/Ag 0.222 | A1461D: REVEL 0.671, CADD 24.4, AlphaMissense 0.7346; benign D2449N retained separately | RIPA, VWF_MULTIMER | Multi-variant case preserves both variants and unknown phase; no cis/trans inference |

### Canonical domain correction

The first summary used a coarse mature-protein interval map, which over-broadly assigned residues near 1500 to A1. The workflow now uses UniProt P04275 canonical feature coordinates in prepro-VWF numbering:

- VWFA1: 1277–1453
- VWFA2: 1498–1665
- VWFC2: 2429–2495

Therefore:

- p.Ala1500Val → VWFA2
- p.Ser1506Leu → VWFA2
- p.Val1316Met → VWFA1
- p.Arg1205His → TIL4–VWFA1 linker
- p.Ala1461Asp → VWFA1–VWFA2 linker
- p.Asp2449Asn → VWFC2

## Interpreting CASE_015 against the physician example

The physician example for P1413L combines:

1. domain location;
2. exome/genome frequency;
3. ClinVar/HGMD classification;
4. REVEL/CADD/AlphaMissense;
5. patient phenotype from a specific PMID;
6. recombinant-expression and functional experiments;
7. ABO/penetrance and clearance context;
8. final cautious conclusion favoring type 1 / low VWF if causal.

For CASE_015 / p.Ser1506Leu, the current system now retrieves and records:

- A1 domain annotation;
- gnomAD exome AF;
- ClinVar Pathogenic/Likely pathogenic;
- REVEL 0.796;
- CADD 26.7;
- AlphaMissense 0.812;
- PMC full-text excerpts about multimers, RIPA, type 2A/2B, and laboratory phenotypes;
- explicit absence of all second-level test results.

It therefore reaches the evidence-collection and evidence-card stage. It still does not reliably make the final leap to “functional disturbance is small; type 1 / low VWF is more likely than type 2” because:

1. retrieved full-text excerpts may discuss other variants/families;
2. exact-variant-to-phenotype linkage is not yet robust;
3. recombinant secretion, dominant-negative effect, and stimulated release are not yet structured as experimental claims;
4. ABO genotype and VWF clearance data are absent;
5. full VWD VCEP ACMG rules are not encoded.

### Where the physician's “recombinant experiment” conclusion comes from

The statement that a variant has little functional effect is not inferred from ClinVar or REVEL. It is read from a variant-specific case report or functional study, usually from sentences and tables describing:

- recombinant VWF expression;
- secretion into medium;
- intracellular storage;
- stimulated release;
- multimer distribution;
- dominant-negative behavior;
- rescue or structure/minimal-change experiments.

To close this gap, the workflow now marks each PubMed/PMC document as `variant_specific=true` only when the full text or abstract explicitly contains the cDNA, protein, or rsID term. Literature phenotype observations from non-specific review articles are retained as context but are not allowed to drive subtype tendency. This prevents a generic VWD review from being misattributed to the patient's variant.

## Recommended pre-batch priorities

### Completed before this batch

1. Reordered the graph to gene-first evidence acquisition.
2. Limited second-level recommendations to a small evidence-constrained set.
3. Added VWF domain and variant-class annotation.
4. Added PMC full-text retrieval and bounded code-as-search excerpts.
5. Added initial literature phenotype extraction.
6. Added structured conflict and missing-evidence analysis.
7. Added initial ACMG hints (PM2/PP3/BS3 candidates).
8. Preserved complete SQLite checkpoints, state history, debug traces, and final reports.

### Highest-value next changes

1. **Variant-specific literature linkage**
   - Keep exact variant query first.
   - Require the full text to mention the protein change, cDNA change, rsID, or patient/family label before attributing phenotype evidence.

2. **Structured functional-experiment extractor**
   - Extract secretion, storage, stimulated release, multimer distribution, dominant-negative effect, and recombinant-expression conclusions.
   - Convert each into a cited claim with direction and limitation.

3. **ClinVar review metadata**
   - Extract review status, submitter count, conflicts, and last evaluation.

4. **ABO / clearance context**
   - Add optional ABO genotype and VWF clearance inputs.
   - Treat DDAVP 0/1/4h as a diagnostic response/clearance probe rather than treatment advice.

5. **Full VWF VCEP ACMG rule table**
   - Encode type-specific rules before reporting pathogenicity or subtype confidence.

## Reproducibility and safety

- 5/5 batch cases completed.
- 0 fabricated second-level observations.
- Every recommended action came from the fixed action enum.
- Every case required expert review.
- Multi-variant CASE_021 retained both variants with phase unknown.
- Tool errors and `not_found` results are retained explicitly and never interpreted as benign evidence.

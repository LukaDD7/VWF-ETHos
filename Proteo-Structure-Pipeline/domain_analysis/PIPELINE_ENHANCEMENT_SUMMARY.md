# VWF Type-2 Pipeline Enhancement Summary

## Completed Tasks

### 1. Unclassified Variants Report Generated
**File**: `domain_analysis/unclassified_domains_report.md`

**14 variants identified** in domains not covered by the 3 primary papers:
- **D1-D2 (Propeptide)**: 2 variants (V89A, L91P) → Type 2A/IIC
- **D4**: 9 variants → Type 2A (multimerization/secretion defects)
- **C-Terminal (CT)**: 3 variants → Type 2A or 2M

### 2. Literature Search Results

#### D1-D2 (Propeptide)
- **PMID 40958414** (2026): Exogenous propeptide supplementation corrects Type 2A/IIC
- **PMID 12176890** (2002): D1 domain critical for pH-dependent multimerization

#### D4 Domain
- **PMID 35734101** (2022): D4 mutations in Afro-Caribbean patients - underrecognized cause of VWD
- **PMID 40687385** (2025): Systematic analysis of D4/C-domain mutations

#### C-Terminal Domain
- Limited specific literature available
- General evidence from PMID 40687385 covers C-domain
- Mechanism inferred from domain function

### 3. Classifier Updated
**File**: `domain_analysis/vwf_type2_literature_based_classifier.py`

Added classification rules for:
- **D1-D2**: Score Type 2A with high confidence (multimerization defect)
- **D4**: Score Type 2A (multimerization/secretion defect)
- **C1-C2**: Score both Type 2M and 2A (collagen binding or multimerization)
- **C3-C6**: Score Type 2A (multimerization)
- **CK**: Score Type 2A with highest confidence (essential for dimerization)

## Key Findings Summary

### D1-D2 Propeptide (Positions 23-763)
**Mechanism**: Defective multimerization → Type 2A/IIC
- D1 domain drives inter-subunit disulfide bond formation
- Mutations cause loss of HMW multimers
- Normal VWF:Ag, reduced VWF:RCo
- Can be corrected by exogenous propeptide

**Variants**: V89A, L91P
**Classification**: Type 2A (IIC subtype)
**Confidence**: High

### D4 Domain (Positions 1875-2255)
**Mechanism**: Multimerization/secretion defects → Type 2A
- Subunit dimerization interface
- ER-to-Golgi trafficking
- Often misdiagnosed as Type 1

**Variants**: P1888L, E1939K, L1974P, R1986C, V1996M, R2006C, V2015L, V2034A, W2062R
**Classification**: Type 2A
**Confidence**: Medium

### C-Terminal Domain (Positions 2256-2813)

#### C1-C2 (2256-2392)
**Mechanism**: Collagen binding OR multimerization
**Variants**: V2465M
**Classification**: Type 2M or 2A
**Confidence**: Low (ambiguous)

#### C3-C6 (2393-2722)
**Mechanism**: Multimerization
**Variants**: S2775C
**Classification**: Type 2A
**Confidence**: Medium

#### CK Domain (2723-2813)
**Mechanism**: Dimerization defect (severe)
**Variants**: P2801S
**Classification**: Type 2A
**Confidence**: High

## Updated Classification Rules

The classifier now implements the following logic:

```python
# D1-D2 propeptide -> Type 2A/IIC
if domain in ["propeptide_D1", "propeptide_D2"]:
    subtype_scores["2A"] += 2.0  # Base score
    if domain == "propeptide_D1":
        subtype_scores["2A"] += 0.5  # D1 bonus

# D4 domain -> Type 2A
if domain == "D4":
    subtype_scores["2A"] += 1.5

# C-terminal (position-dependent)
if domain in ["C1", "C2"]:
    subtype_scores["2M"] += 0.8  # Collagen binding
    subtype_scores["2A"] += 0.5  # Multimerization
elif domain in ["C3", "C4", "C5", "C6"]:
    subtype_scores["2A"] += 1.0  # Multimerization
elif domain == "CK":
    subtype_scores["2A"] += 2.5  # Severe dimerization defect
```

## Recommended Clinical Workup

| Domain | Test | Expected Finding |
|--------|------|------------------|
| D1-D2 | Multimer analysis | Loss of HMW multimers |
| D1-D2 | VWFpp antigen | May be abnormal |
| D4 | Multimer analysis | Variable loss of HMW |
| D4 | Secretion studies | Impaired ER-Golgi transport |
| CT | Collagen binding | May be reduced |
| CT | Multimer analysis | Loss of HMW multimers |

## Files Updated

1. `domain_analysis/unclassified_domains_report.md` - Comprehensive report
2. `domain_analysis/vwf_type2_literature_based_classifier.py` - Updated classification rules
3. `domain_analysis/README.md` - Documentation (existing)

## Next Steps (Optional)

1. **Validate predictions**: Test updated classifier on the 14 unclassified variants
2. **AF3 structural analysis**: Extract pLDDT/RMSD features for these domains
3. **Expand literature**: Search for additional CT domain papers
4. **Clinical correlation**: Compare predictions with patient phenotype data

## References

### Propeptide (D1-D2) Domain
1. PMID 35148377 - Structural basis of VWF multimerization and tubular storage (2022) - Cryo-EM structures
2. PMID 40958414 - Propeptide correction of Type 2A/IIC (2026)
3. PMID 12176890 - D1 domain in multimerization (2002)
4. PMID 19506357 - VWD type 2C (2A subtype IIC) characteristics (2009)
5. PMID 17895385 - Cys residues essential for VWF multimer assembly (2007)
6. PMID 20335223 - N528S mutation in VWF propeptide D2 domain (2010)

### D4 Domain
7. PMID 35734101 - D4 mutations in Afro-Caribbeans (2022)
8. PMID 40687385 - D4 and C-domain mutations (2025)
9. PMID 19506359 - VWD type 1/2E with D3, D4, B1-B3, C1-C2 mutations (2009)

### Primary VWF Structure Papers
10. Lenting et al. Blood 2024 - VWF structure and function
11. Atiq & O'Donnell Blood 2024 - Novel VWF functions
12. Haberichter & O'Donnell Haematologica 2026 - VWF structure

---

*Report completed: 2026-04-02*

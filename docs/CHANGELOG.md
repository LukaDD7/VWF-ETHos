# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased] - Alphagenome-DFR-Phase3 Architecture

### Completed

#### Phase 1: Data Integration ✅
- `merge_alpha_features.py` - New cDNA extraction method from AG consequences
- AG matching rate: 67% (up from 10% with naive offset method)
- Output: `VWF_Alpha_Matrix.parquet` (100 variants, 59 with AF3)

#### Phase 2: AgenticVWFClassifier ✅
- `agentic_vwf_classifier.py` - 3 Expert Agent architecture
- Expert 1: StructuralExpert (AF3 pLDDT analysis)
- Expert 2: TranscriptomicExpert (AG RNA/splice delta)
- Expert 3: ClinicalGeneticistAgent (Interpretable logical rules)
- **Removed hardcoded thresholds** - Uses adaptive percentile-based thresholds
- **Removed nested if-else** - Rule-based decision hierarchy

### Validation Results (Phase 2)
| Known ↓ / Predicted → | 2A | 2B | 2M | 2N | uncertain |
|---|---|---|---|---|---|
| **2A** (43) | 26 | 0 | 4 | 2 | 11 |
| **2B** (12) | 0 | 2 | 8 | 0 | 2 |
| **2M** (25) | 0 | 0 | 24 | 0 | 1 |
| **2N** (20) | 7 | 0 | 0 | 13 | 0 |

### Pending Tasks
- [ ] Task 2: Upgrade Expert 1 for Type 2B (AF3 interface PAE analysis)
- [ ] Task 3: Execute validation on 59 AF3 variants
- [ ] Task 4: Fix remaining 33% AG unmatched variants
- [ ] Task 5: D4 domain variant data for Type1/2A competition validation

---

## [2026-04-02] - Architecture Refactor: Domain-Analysis Standard

### Breaking Changes

#### Removed (Archived)
- `scripts/phase4_domain_feature_extraction.py` - Single-label classifier (archived to `archived/`)
- `results/phase4_domain_features.csv` - Old output format

### New Standard: Domain-Analysis Pipeline

#### Added
- **DEV_PROTOCOL_STANDARD.md** - 4-phase development protocol
- **Competitive Classifier** (`vwf_competitive_classifier.py`) - Multi-label with competition resolution
- **Evidence-Based Classifier** (`vwf_evidence_based_classifier.py`) - Literature evidence scoring

### Key Improvements
| Dimension | Old | New | Status |
|-----------|-----|-----|--------|
| Classification | Single-label | Multi-label competitive | Fixed |
| Literature | None | Full evidence database | Fixed |
| Explainability | Feature importance | Full reasoning path | Fixed |

---

## [2026-03-30] - Phase 4: Domain Feature Extraction

### Added
- `scripts/phase4_domain_feature_extraction.py` - Domain feature extraction and ML classification
- Outputs: `phase4_domain_features.csv`, `phase4_ml_predictions.csv`
- Random Forest + XGBoost classifiers

---

## [Earlier] - Phase 1-3 Pipeline

- Phase 1: Variant filtering (`01_filter_target_vus.py`)
- Phase 2: Liftover (`02_preprocess_and_liftover.py`)
- Phase 3: AlphaGenome inference (`03_run_alphagenome_inference.py`)
- Phase 4 (old): Analysis (`04_analyze_and_visualize.py`)

---

*Last updated: 2026-04-21*
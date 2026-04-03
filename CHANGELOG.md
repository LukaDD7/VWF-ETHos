# CHANGELOG

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Planned
- [ ] Complete Domain-Analysis Phase 1: Input/Output specification
- [ ] Complete Domain-Analysis Phase 2: Key step design
- [ ] Implement new main pipeline integrating competitive classifier
- [ ] Create comprehensive technical report
- [ ] Validate consistency between old and new results

---

## [2026-04-02] - Architecture Refactor: Domain-Analysis Standard

### Breaking Changes

#### Removed (Archived)
- `scripts/phase4_domain_feature_extraction.py` - Single-label classifier (archived to `archived/2026-03-30_phase4/`)
- `results/phase4_domain_features.csv` - Old output format
- `results/phase4_ml_predictions.csv` - Single-label predictions
- `results/phase4_feature_importance.png` - Feature importance only

### New Standard: Domain-Analysis Pipeline

#### Added
- **DEV_PROTOCOL_STANDARD.md** - 4-phase development protocol
- **PHASE4_AUDIT_REPORT.md** - Audit of 3月30日 vs 4月2日 versions
- **Competitive Classifier** (`vwf_competitive_classifier.py`) - Multi-label with competition resolution
- **Evidence-Based Classifier** (`vwf_evidence_based_classifier.py`) - Literature evidence scoring
- **Literature-Based Classifier** (`vwf_type2_literature_based_classifier.py`) - 2024-2026 VWF literature rules
- **Structure Feature Extractor** (`vwf_structure_feature_extractor.py`) - AlphaFold3 CIF parsing
- **Integration Guide** (`INTEGRATION_GUIDE.md`) - Usage documentation

### Changed
- Updated `CLAUDE.md` with new project structure and Domain-Analysis standard
- Reorganized project structure with `archived/` and `domain_analysis/`

### Key Improvements
| Dimension | Old (Mar 30) | New (Apr 2) | Status |
|-----------|-------------|-------------|--------|
| Classification | Single-label | Multi-label competitive | Fixed |
| Literature | None | Full evidence database | Fixed |
| Explainability | Feature importance | Full reasoning path | Fixed |
| Domain pleiotropy | Simple matching | Competition resolution | Fixed |
| Reliability | None | high/medium/low rating | Fixed |

### Metrics
- Output fields: 31 → 46+ (+48%)
- Evidence sources: 0 → 35+ papers
- Classification types: 4 → Multi-label with alternatives

---

## [2026-03-30] - Phase 4: Domain Feature Extraction

### Added
- `scripts/phase4_domain_feature_extraction.py` - Domain feature extraction and ML classification
- Outputs: `phase4_domain_features.csv`, `phase4_ml_predictions.csv`
- Random Forest + XGBoost classifiers

### Known Limitations (in retrospect)
- Single-label output only
- No literature evidence
- No competition resolution for domain pleiotropy
- No reasoning path
- No reliability assessment

---

## [Earlier] - Phase 1-3 Pipeline

- Phase 1: Variant filtering (`01_filter_target_vus.py`)
- Phase 2: Liftover (`02_preprocess_and_liftover.py`)
- Phase 3: AlphaGenome inference (`03_run_alphagenome_inference.py`)
- Phase 4 (old): Analysis (`04_analyze_and_visualize.py`)

---

*Last updated: 2026-04-02*

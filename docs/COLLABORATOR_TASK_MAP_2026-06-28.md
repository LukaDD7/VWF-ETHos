# Collaborator task map

Date: 2026-06-28

This is the short handoff for collaborators. The repository has several older
README/runbook files from earlier phases; use this file plus the linked current
readouts as the entry point.

## Current source of truth

Use these files first:

| Area | File | Notes |
|---|---|---|
| Current classifier | `scripts/agentic_vwf_classifier.py` | Rule-based VWF classifier. A1 Rule6 resolves 2B/2M; A3 Rule7 now uses collagen z-score as confidence evidence. |
| Current eval runner | `scripts/pipeline/evaluate_vwf_classifier_v2.py` | Evaluation-only script. Labels are used after prediction for recall/confusion; includes leakage smoke test. |
| Current evidence matrix | `output/hf_type2m_lof_panel/type2m_lof_panel/analysis/evidence_matrix_with_type2m_lof_hf.csv` | Main matrix after adding the 16 Type2M LOF HF Boltz top-up rows. |
| Current eval output | `output/eval_v2_with_type2m_lof_hf/` | Latest Eval v2 outputs after Type2M top-up. |
| Current eval readout | `docs/TYPE2M_LOF_HF_BOLTZ_PULL_AND_EVAL_2026-06-28.md` | Best single report for data size, recall, leakage test, and Type2M feasibility. |
| System feasibility before top-up | `docs/EVAL_V2_SYSTEM_FEASIBILITY_READOUT_2026-06-27.md` | Still useful for system explanation and older baseline. |
| Type2M mechanism plan | `docs/TYPE2M_LOF_MECHANISM_AND_FEATURE_PLAN_2026-06-27.md` | Mechanism hypotheses and next MD/Boltz paths. |
| 2B MD readout | `docs/NEW_2B_SALTBRIDGE_MD_RECALL_READOUT_2026-06-26.md` | Salt-bridge MD recall result for non-hotspot 2B. |

Treat `README.md`, `docs/PROJECT_MOC.md`, and `docs/PIPELINE_RUNBOOK.md` as
historical orientation only. They contain early-phase statements that are now
outdated for the current Boltz/MD classifier.

## Repository map

| Path | Meaning |
|---|---|
| `scripts/agentic_vwf_classifier.py` | Classifier implementation. |
| `scripts/pipeline/` | Data integration, Boltz parsing, MD preparation, MD extraction, and eval scripts. |
| `docs/` | Runbooks and readouts; many are date-stamped. Prefer latest date for the same topic. |
| `output/boltz2_vwd_functional_panel/` | Original functional Boltz panel evidence matrix. |
| `output/hf_type2m_lof_panel/type2m_lof_panel/analysis/` | Type2M LOF HF top-up analysis tables and combined evidence matrix. |
| `output/eval_v2_with_type2m_lof_hf/` | Latest classifier predictions, recall summary, confusion tables, priority queues. |
| `output/labeled_variants_all.csv` | Positive label inventory after de-duplication and conflict detection. |
| `output/labeled_variants_conflicts.csv` | Variants appearing under multiple labels. Needs manual curation. |
| `output/new_2b_for_boltz_reconciled.csv.report.csv` | Reconciliation report for the newer 2B table, including shifted numbering, duplicates, and conflicts. |

## Data inventory snapshot

Current clean positive label inventory:

| Table | Rows | Notes |
|---|---:|---|
| `output/labeled_variants_all.csv` | 437 | All positive labels after normalization. |
| `output/labeled_variants_conflicts.csv` | 53 | Same `aa_change` reported under multiple labels. |
| combined evidence matrix | 613 | Original 597 rows + 16 HF Type2M LOF top-up rows. |
| latest clean supported eval rows | 340 | 115 Type1 + 118 2A + 38 2B + 53 2M + 16 2N. |
| latest Type2 eval rows | 225 | 118 2A + 38 2B + 53 2M + 16 2N. |

Current Type2 recall with existing MD:

| Label | Recall |
|---|---:|
| 2A | 70/118 = 59% |
| 2B | 29/38 = 76% |
| 2M | 21/53 = 40% |
| 2N | 12/16 = 75% |
| Type2 total | 132/225 = 59% |

## What hotspot means

`hotspot` is an A1 position-level prior for recurrent Type2B-associated sites
implemented as `TWO_B_HOTSPOT_POS` in `scripts/agentic_vwf_classifier.py`.

`non-hotspot` means the position is outside that recurrent Type2B set. It does
not mean the variant is unimportant. It means the classifier should rely more on
functional evidence such as GPIb/heparan Boltz axes and MD features.

## Sample overlap / conflict tasks

These tasks are safe and useful for a collaborator.

### Task A: manually curate label conflicts

Input:

- `output/labeled_variants_conflicts.csv`

Goal:

- For each conflict row, decide whether it is:
  - true biological pleiotropy,
  - literature/reporting ambiguity,
  - numbering mismatch,
  - table extraction artifact,
  - or a real multi-subtype patient phenotype.

Suggested output:

- `output/curation/labeled_variants_conflicts_curated.csv`

Suggested columns:

```text
aa_change, current_labels, curated_label, curation_status, evidence_note, curator, date
```

Recommended status values:

```text
keep_conflict, choose_single_label, numbering_issue, extraction_artifact, needs_literature_review
```

High-priority A1 conflicts to review first:

```text
C1272R, L1278P, V1279F, V1279I
```

### Task B: audit the newer 2B table reconciliation

Input:

- `output/new_2b_for_boltz_reconciled.csv.report.csv`

Goal:

- Review rows with status containing:
  - `conflict_*`
  - `duplicate_*`
  - `mature_shifted`
  - `unresolved`

Why:

- This table mixes pre-pro numbering and mature numbering. Some apparent new
  2B rows are shifted-numbering aliases or conflict with existing 2M/2A labels.

Suggested output:

- `output/curation/new_2b_reconciliation_curated.csv`

### Task C: verify negative / benign overlap

Input:

- `output/hf_type2m_lof_panel/type2m_lof_panel/analysis/evidence_matrix_with_type2m_lof_hf.csv`

Goal:

- Review rows with `label_policy = negative_positive_overlap`.
- Decide whether each should be excluded from classifier calibration, kept as
  ambiguous control, or re-labeled after source review.

Current count:

- `negative_positive_overlap = 25`

Suggested output:

- `output/curation/negative_positive_overlap_curated.csv`

### Task D: A1-2M uncertain mechanism review

Input:

- `output/eval_v2_with_type2m_lof_hf/eval_v2_2m_lof_md_priority_queue.csv`
- `output/eval_v2_with_type2m_lof_hf/eval_v2_predictions_slim.csv`

Goal:

- Prioritize A1 true 2M variants that remain `uncertain` or are miscalled as
  2B for new MD or mechanism review.

Why:

- Static Boltz GPIb/heparan axes do not catch many A1-2M cases. These are the
  best candidates for A1-GPIb complex MD or 7A6O closed-state MD.

## MD bottleneck and Mac mini feasibility

Current production MD settings are in:

- `scripts/pipeline/run_md_variant_direct.sh`
- `scripts/pipeline/run_7a6o_variant_direct.sh`

Defaults:

| Setting | Value |
|---|---:|
| timestep | 2 fs |
| production length | 50 ns |
| steps per ns | 500,000 |
| steps per 50 ns variant | 25,000,000 |
| NVT | 50 ps |
| NPT | 200 ps |
| trajectory write interval | 100 ps in production |

Why MD is slow:

1. Production MD is a time-integration workload: every 2 fs step depends on the
   previous step, so one variant cannot be trivially split into thousands of
   independent chunks.
2. Each step computes nonbonded interactions, PME electrostatics, constraints,
   neighbor-list updates, and thermostat/barostat terms.
3. Solvent and ions dominate atom count; most compute is moving water, not just
   the protein.
4. Each variant needs its own 50 ns trajectory, and useful calibration needs many
   variants plus ideally replicates.
5. CPU/GPU balance matters. GROMACS is fastest with CUDA GPU-resident mode
   (`-nb gpu -pme gpu -bonded gpu -update gpu`). CPU-only or partial offload is
   much slower.

Mac mini recommendation:

| Use case | Mac mini suitability |
|---|---|
| inspect scripts / parse outputs / write docs | Good |
| prepare YAMLs, manifests, curation tables | Good |
| run FoldX small tests | Possible |
| run GROMACS EM/NVT smoke test | Possible if GROMACS is installed |
| run one very short MD, e.g. 0.1-1 ns | Possible but slow |
| run 50 ns production for many variants | Not recommended |
| run Boltz-2 production | Not recommended |

Practical guidance:

- Mac mini is fine for debugging the pipeline and validating one tiny trajectory.
- For real 50 ns production, use A40/H200 or another CUDA GPU host.
- If forced to use Mac mini, set `NS=0.1` or `NS=1` for smoke only and do not use
  the result for classifier calibration.

## Recommended collaborator workflow

1. Pull latest repo.
2. Read this file and `docs/TYPE2M_LOF_HF_BOLTZ_PULL_AND_EVAL_2026-06-28.md`.
3. Pick one task from Task A-D.
4. Write curated output under `output/curation/`.
5. Do not overwrite raw outputs or rerun large MD/Boltz jobs unless assigned.
6. If adding a new curation table, include a short README or `notes` column with
   exact source evidence.

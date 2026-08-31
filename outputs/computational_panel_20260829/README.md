# Computational panel, 2026-08-29

This bundle consolidates the Type 1 and Type 2B pre-MD computational evidence,
the targeted 7A6O AIM-A1 MD subset, and the offline agent runs.

## Contents

- `pre_md/`: AlphaGenome and Boltz-2 cross-panel summaries and interpretation.
- `boltz/raw/`: completion markers and run logs for the 32-job panel.
- `boltz/yamls/`: the 32 Boltz-2 input YAMLs used for this panel.
- `md/`: targeted AIM-A1 equilibrium-MD QC, time series, feature matrix, and
  lightweight GROMACS intermediates.
- `md/`: also includes the multi-module MD/SMD starting-structure inventory,
  module readout/AI reference-distribution plan, and consolidated A1-GPIbα,
  AIM salt-bridge, and slow025 SMD result layers.
- `agent_run_type1_*` and `agent_run_type2b_*`: deterministic offline agent run
  records and summaries.

## Boltz-2

The panel contains 32 unique Boltz-2 jobs: 29 variant jobs and three shared A1
WT baselines. The job-level metrics are already represented in
`outputs/type1_panel_agent_20260828/server_bundle/results/boltz/boltz_results_summary.csv`
and
`outputs/type2b_panel_agent_20260829/server_bundle/results/boltz/boltz_results_summary.csv`,
with the combined interpretation in `pre_md/boltz_variant_assay_summary.csv`.
The missing lightweight reproducibility layer was the set of input YAMLs, which
is now included under `boltz/yamls/`. No additional Boltz computation is needed
for this submitted panel.

Use iPTM only for complex/interface assays. For monomeric or conformational
context assays, including `a1_aim_autoinhibition_context` and
`a2_folded_stability`, use pTM, pLDDT, or local structural features; iPTM is
expected to be zero or otherwise not informative.

## AlphaGenome

AlphaGenome results are already complete for all 16 requests and are committed
under the Type 1 and Type 2B server bundles. The cross-panel summary is in
`pre_md/alphagenome_case_summary.csv`. Each request has 15 of 19 selected
scorer views; ATAC, contact maps, polyadenylation, and ATAC-active views are
unavailable for the selected VWF-relevant biosamples. No duplicate AlphaGenome
output was added here.

## MD

The targeted MD subset covers seven matched systems for 50 ns each. The initial
commit included only the AIM-A1 contact summary/timeseries and QC. This update
adds RMSD time series, AIM minimum-distance time series and summaries, a
combined MD feature matrix, an artifact manifest, and the lightweight XVG/NDX
intermediates. See `md/README.md` for details. The large concatenated XTC
trajectories remain local-only and are recorded by checksum in
`md/artifact_manifest.csv`.

The MD layer now separates three things:

1. **Completed results**: the original seven-system AIM-A1 panel, plus the
   additional A1-GPIbα interface, AIM salt-bridge, and slow025 SMD calibration
   layers indexed by `md/md_result_layers.csv`.
2. **Starting structures**: `md/module_structure_inventory.csv` records the
   preferred starting structure for each of the 15 functional assay axes.
   Experimental PDBs take precedence; Boltz-2 is only a fallback when no
   experimental structure is locally available.
3. **Execution plan**: `md/module_md_readout_plan.csv` defines module-specific
   equilibrium-MD readouts, conditional SMD gates, and AI reference-distribution
   features. The full rationale is in
   `docs/VWF_MULTIMODULE_MD_SMD_AI_REFERENCE_PLAN_20260831.md`.

Only A1/AIM currently has a complete matched experimental-structure MD panel.
The other modules have no completed MD trajectories yet. The priority
experimental PDBs are now downloaded: 3GXB for A2, 4DMU for A3-collagen, 6N29
for D′D3, and 4NT5 for CK/CTCK. Modules without an experimental structure use
exploratory Boltz-2 fallbacks.

## Next-analysis hints from the current documents

- Treat all model outputs as mechanism evidence, not clinical classifications.
- Distinguish Type 2B from Type 2M-A1 using AIM/open-closed and exposure
  evidence in addition to forced A1-GPIbα complex confidence; the forced
  complex alone does not answer whether A1 opens at low shear.
- Build downstream agent inputs from the assay-matched feature matrices rather
  than hard-coded subtype rules.
- Preserve uncertainty when only an exploratory monomer/conformational feature
  is abnormal.
- The broader next-phase plan is to generate evidence-matrix and
  assay-distribution analyses, presentation figures, and an Agent v3 that
  ingests the structured evidence matrix.

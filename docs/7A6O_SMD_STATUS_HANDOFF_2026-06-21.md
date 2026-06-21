# 7A6O SMD Status Handoff - 2026-06-21

This note captures the local A40-side SMD state after the 2026-06-21 handoff. It is meant to make the current disagreement explicit: the full mutant queue appears healthy, while the current WT SMD baseline is not yet a reliable comparator.

## Objective

Use steered MD-derived mechanical features from the 7A6O AIM-A1 system to improve VWF type 2B vs 2M separation in the classifier.

The intended feature axis is:

- lower AIM-unfolding / detachment force -> easier autoinhibition release -> more 2B-like;
- higher AIM-unfolding / detachment force -> more WT/2M-like.

This axis should complement the existing equilibrium-MD contact/RMSD/salt-bridge features. It should not replace them until the SMD assay is calibrated on multiple known 2B and 2M controls.

## Current Background Jobs

At the last local check, the full SMD queue was still running and did not show a new fatal error.

Active prep jobs:

| Variant | Stage | GPU | Last observed progress |
| --- | --- | --- | --- |
| R1306Q | SMD prep NPT | 1 | ~110 ps / 200 ps |
| R1308C | SMD prep NPT | 4 | ~120 ps / 200 ps |
| I1309V | SMD prep NPT | 6 | ~120 ps / 200 ps |
| S1310F | SMD prep NPT | 5 | ~110 ps / 200 ps |

The active commands were `gmx mdrun -deffnm npt -ntmpi 1 -ntomp 16 -nb gpu -pme gpu -pin on`, wrapped by `scripts/pipeline/prep_7a6o_smd.sh`.

This is consistent with a healthy long-running queue. The relatively low GPU utilization during prep/NPT is expected for this restrained, mixed CPU/GPU workload and does not by itself indicate a stalled run.

## Completed SMD Probe Data

After switching the pull geometry to `direction-periodic`, the following probe replicas completed:

| Variant | Label | Complete reps | Mean rupture force |
| --- | --- | ---: | ---: |
| R1306W | 2B | 3 | ~1116.7 pN |
| R1374H | 2M | 3 | ~993.4 pN |

Interpretation caveat: this two-variant probe does not yet support the simple expectation that a 2B mutation always lowers the peak pulling force. This may be due to limited sampling, force-coordinate definition, WT/reference mismatch, or real mutation-specific mechanics. Do not treat these two rows as calibrated classifier evidence.

## WT Baseline Issue

The current `WT` SMD run should not be used as a clean comparator yet.

Observed behavior:

- WT SMD rep1 on GPU0 ran anomalously slowly relative to mutant probes.
- A retry with `NTOMP=16` only reached ~10 ps in ~10 minutes and was stopped.
- The partial WT `smd_rep1_*` outputs were archived under `output/gromacs_md_autoinhib/WT/smd/partial_slow_gpu0_retry_20260621_152957/`.

Working hypothesis:

- The historical WT input under `output/gromacs_md_autoinhib/7A6O_WT/md_7a6o` is not fully pipeline-matched to the FoldX-derived mutant inputs.
- Comparing this WT directly against FoldX-mutant SMD may confound mutation effect with input/topology/provenance differences.

Practical consequence:

- Full mutant SMD can continue.
- Do not use the current historical WT SMD partial output for feature fitting.
- Prefer a pipeline-matched WT baseline before interpreting absolute WT-vs-mutant force shifts.

## Analysis Script Safeguard

`analyze_7a6o_smd.py` was updated in commit `f1a8a00` so that incomplete SMD replicas are ignored.

A replica is now included only if these files exist together:

- `smd_repN_pullf.xvg`
- `smd_repN_pullx.xvg`
- `smd_repN.gro`

This prevents killed or failed runs from contributing partial force curves to `md_7a6o_smd_features.csv`.

## Recommended Next Steps

1. Let the current mutant SMD queue continue. It is already running and appears healthy.
2. After each batch finishes, run `python3 scripts/pipeline/analyze_7a6o_smd.py --input output/gromacs_md_autoinhib --output output/md_7a6o_smd_features.csv`.
3. Treat the first 14-mutant result as calibration/exploration, not final classifier evidence.
4. Rebuild a pipeline-matched WT baseline before using WT in the SMD feature axis. A reasonable label would be `WT_FoldX` or `WT_matched` to avoid overwriting historical `WT` output.
5. If the 14-mutant 3-rep screen shows partial separation, increase SMD replicates for boundary or clinically important variants to 5-10 reps.
6. Only after separation is stable should `smd_2b_score` be wired into the classifier as a type-2B positive feature.

## Time Expectations

Based on local observations:

- SMD prep NPT target is 200 ps; the current first full-queue batch was around 55-60% complete at the last check.
- Each SMD production replicate is configured as 5 ns at 1 nm/ns pulling speed.
- Three replicas per variant should be expected to take hours per variant, depending on GPU/CPU scheduling and concurrent jobs.

The safest operational expectation is to check again after several hours rather than minutes.

## Commands For Status Checks

```bash
pgrep -af 'prep_7a6o_smd|run_7a6o_smd|gmx.*mdrun|7a6o_smd_full'
nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used --format=csv,noheader
tail -40 output/gromacs_md_autoinhib/R1306Q/smd/npt.log
tail -40 output/gromacs_md_autoinhib/R1308C/smd/npt.log
tail -40 output/gromacs_md_autoinhib/I1309V/smd/npt.log
tail -40 output/gromacs_md_autoinhib/S1310F/smd/npt.log
```

## Bottom Line

There are two separate issues:

- The mutant SMD queue appears operational and should continue.
- The current historical WT SMD baseline is not yet trustworthy for classifier calibration.

The immediate priority is to finish the mutant SMD screen, analyze only complete replicas, and then decide whether the force feature is useful enough to justify a matched WT rebuild and additional replicate sampling.

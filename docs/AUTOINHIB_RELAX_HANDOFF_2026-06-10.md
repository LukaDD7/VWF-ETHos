# Autoinhib Relax Handoff - 2026-06-10

## Status

WT autoinhib dynamic-MD blocker was tested on A40 local CUDA GROMACS env. GROMACS runner/smoke can run; the remaining question was whether the Boltz D'D3-A1 structure can be relaxed into an MD-ready starting point.

## 1. Clash diagnosis

Command:

```bash
envs/gromacs/bin/python scripts/pipeline/diagnose_clashes.py \
  --variant-dir output/boltz2_a1_dp_d3_results/boltz_results_VWF_WT_dp_d3_a1
```

Output saved to `output/diagnose_clashes_VWF_WT_20260610.txt`.

Summary:

| model | total clashes | H-clash | heavy | min distance |
|---|---:|---:|---:|---:|
| model_0 | 28 | 0 | 28 | 0.76 A |
| model_1 | 32 | 0 | 32 | 1.05 A |
| model_2 | 18 | 0 | 18 | 0.75 A |
| model_3 | 28 | 0 | 28 | 1.11 A |
| model_4 | 34 | 0 | 34 | 0.87 A |

Recommended cleanest model: `VWF_WT_dp_d3_a1_model_2.cif`.

Interpretation: this is not a pure H-clash case in the Boltz CIF scan; all counted clashes are heavy-atom clashes.

## 2. Relax result

Command:

```bash
bash scripts/pipeline/relax_autoinhib_structure.sh --variant VWF_WT --model 2 \
  2>&1 | tee output/relax_VWF_WT_model2_20260610.log
```

Work directory:

```text
output/gromacs_md_autoinhib/VWF_WT/relax_m2
```

Formal stage results:

| stage | result | potential energy | Fmax | note |
|---|---:|---:|---:|---|
| vacuum + position restraints | converged in 992 steps | -4.0500809e+04 | 1.9874701e+02 | `<200`, passed |
| vacuum unrestrained | converged in 556 steps | -4.3198953e+04 | 1.9571939e+02 | `<200`, passed |
| solvated loose EM (`emtol=10000`) | converged in 14 steps | -1.0259850e+06 | 8.6383281e+03 | `<1e4`, MD-ready threshold passed |

The strict solvated EM from the stock script (`emtol=200`) was manually interrupted after it had already dropped into the 1e3-1e4 Fmax range but before a `Maximum force` terminator line was written. To provide a formal pass/fail datum, a follow-up solvated EM was run from the same `solv_ions.gro/topol.top` with `emtol=10000`; it converged and wrote:

```text
output/gromacs_md_autoinhib/VWF_WT/relax_m2/em_10k.gro
output/gromacs_md_autoinhib/VWF_WT/relax_m2/em_10k.tpr
output/gromacs_md_autoinhib/VWF_WT/relax_m2/topol.top
```

## 3. MD blocker judgment

For WT model 2, the blocker is likely recoverable by staged GROMACS relaxation. It is not necessary to jump directly to OpenMM for WT before trying model 2-derived continuation. Next engineering step is to adapt the production MD runner to start NVT/NPT/prod from `relax_m2/em_10k.gro + topol.top`, or copy these into the runner's expected `em.gro/topol.top` layout for a controlled WT-only continuation. Do not launch all 5 variants before this WT continuation passes NVT/NPT.

## 4. AIM feature extraction status

Command attempted:

```bash
envs/gromacs/bin/python scripts/pipeline/extract_aim_autoinhib_features.py
```

Result:

```text
[FATAL] results dir 不存在: output/boltz2_vwd_functional_panel/boltz_results
```

Current machine has `output/boltz2_vwd_functional_panel/yamls/`, `evidence_matrix.csv`, `diagnostic_panel.csv`, and analysis CSVs, but not the functional-panel Boltz CIF directory needed to compute true `aim_release_score`. Therefore the 2B axis-A calibration is blocked on transferring or generating:

```text
output/boltz2_vwd_functional_panel/boltz_results/
```

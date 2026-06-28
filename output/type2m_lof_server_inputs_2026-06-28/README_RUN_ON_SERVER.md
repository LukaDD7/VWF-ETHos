# Type 2M LOF server input package

Generated on 2026-06-28.

Run from repository root on the server.

## 1. Boltz top-up for missing clean 2M

```bash
export GPUS=4
# optional: export GPU_IDS=4,5,6,7
bash output/type2m_lof_server_inputs_2026-06-28/boltz_missing_2m/run_boltz_missing_2m.sh
```

Input YAMLs:

`output/type2m_lof_server_inputs_2026-06-28/boltz_missing_2m/run_panel/yamls`

This package contains 30 YAML jobs: 16 clean 2M variants projected onto the relevant A1/A3 LOF axes plus 4 WT baselines.

## 2. A1-GPIb complex MD, priority P0

```bash
export FOLDX=/path/to/foldx
export GMX=/path/to/gmx
export GPU_IDS=4,5,6,7
export NS=50
bash output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba/run_a1_gpiba_p0_md.sh
```

Default P0 list:

`output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba/a1_gpiba_p0_plus_anchor_variants.txt`

This uses `structures/1SQ0.pdb`, mutates VWF chain A with `offset=763`, relaxes into `output/gromacs_md_a1_gpiba/<variant>/relax_pdb`, then runs production MD in `md_a1_gpiba`.

## 3. A1-GPIb all recommended MD

After P0 succeeds, run the full recommended A1-GPIb list:

```bash
export FOLDX=/path/to/foldx
export GMX=/path/to/gmx
export GPU_IDS=4,5,6,7
export NS=50
bash output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba/run_a1_gpiba_all_recommended_md.sh
```

## 4. 7A6O closed-state secondary MD

```bash
export FOLDX=/path/to/foldx
export GMX=/path/to/gmx
export GPU_IDS=4,5,6,7
export NS=50
bash output/type2m_lof_server_inputs_2026-06-28/md_7a6o_closed_state/run_7a6o_2m_closed_state_md.sh
```

This reuses the existing 7A6O AIM-A1 closed-state pipeline and writes into `output/gromacs_md_autoinhib`.

## Key manifests

- `boltz_missing_2m/run_panel/job_manifest.csv`
- `md_a1_gpiba/a1_gpiba_md_manifest.csv`
- `md_7a6o_closed_state/7a6o_2m_closed_state_manifest.csv`
- `support_scripts/` contains fallback copies of the generic MD runner scripts.

## Notes

- P0 A1-GPIb MD targets true 2M cases currently miscalled as 2B plus the G1324S literature anchor.
- 2B hard negatives are included in the all-recommended/control list to protect 2B recall when calibrating LOF thresholds.
- Do not overwrite old Boltz/MD output directories unless intentionally rerunning; these scripts are resumable where the underlying runner is resumable.

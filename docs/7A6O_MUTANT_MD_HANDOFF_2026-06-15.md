# 7A6O Mutant MD Handoff - 2026-06-15

## Current state

WT 7A6O production is complete. The final WT trajectory came from a resumed run:

- `output/gromacs_md_autoinhib/7A6O_WT/md_7a6o/md_prod_run.part0002.gro`
- `output/gromacs_md_autoinhib/7A6O_WT/md_7a6o/md_prod_run.part0002.xtc`
- WT log ended at 50 ns on `2026-06-13 03:35:13`.

All 14 FoldX mutants have refined solvated EM start structures:

- `output/gromacs_md_autoinhib/<VARIANT>/relax_pdb/solv_ions_em_refined.gro`
- The refined EM pass used unrestrained steepest descent with `emtol = 1000`, `emstep = 0.001`, `nsteps = 20000`.
- Rationale: the original solvated EM was only 50 steps with protein position restraints. Several mutants still had high local forces, causing LINCS failures or CUDA illegal-address failures immediately after NVT startup.

## Running jobs

The mutant MD jobs were launched with `setsid` so they survive terminal/SSH exit. Do not rely on plain `nohup` in this environment; background children created by normal shell commands were observed to be cleaned up when the tool shell exited.

Active direct-runner queues at launch time:

| GPU | Queue | Main PID at launch | Log |
| --- | --- | ---: | --- |
| 0 | `R1306W`, then `R1374H` | `53428` for R1306W; `77778` waits for R1374H | `output/gromacs_md_autoinhib/R1306W_direct_setsid_20260615_1114.log`, `output/gromacs_md_autoinhib/direct_queue_gpu0_after_R1306W_20260615_1120.log` |
| 1 | `R1306Q`, then `V1316M` | `77779` | `output/gromacs_md_autoinhib/direct_queue_gpu1_20260615_1120.log` |
| 2 | `R1308C`, then `P1337L` | `77780` | `output/gromacs_md_autoinhib/direct_queue_gpu2_20260615_1120.log` |
| 3 | `I1309V`, then `R1341Q` | `77781` | `output/gromacs_md_autoinhib/direct_queue_gpu3_20260615_1120.log` |
| 4 | `S1310F`, then `R1341W` | `77782` | `output/gromacs_md_autoinhib/direct_queue_gpu4_20260615_1120.log` |
| 5 | `W1313C`, then `R1374C` | `77783` | `output/gromacs_md_autoinhib/direct_queue_gpu5_20260615_1120.log` |
| 6 | `V1314F`, then `G1324S` | `77784` | `output/gromacs_md_autoinhib/direct_queue_gpu6_20260615_1120.log` |

At `2026-06-15 11:22`, seven GROMACS GPU processes were visible in `nvidia-smi`. `R1306W` had reached production; the other first-wave mutants were in NVT startup.

## How to check progress

Use these from the repository root:

```bash
pgrep -af 'run_7a6o_variant_direct|gmx.*mdrun|gmx.*grompp'
nvidia-smi
```

Check a specific variant:

```bash
tail -40 output/gromacs_md_autoinhib/R1306W/md_7a6o/md_prod.log
tail -40 output/gromacs_md_autoinhib/R1306Q/md_7a6o/md_nvt.log
tail -40 output/gromacs_md_autoinhib/R1306Q/md_7a6o/md_prod.log
```

Count completed mutant productions:

```bash
find output/gromacs_md_autoinhib -path '*/md_7a6o/md_prod.gro' -type f | sort
find output/gromacs_md_autoinhib -path '*/md_7a6o/md_prod.gro' -type f | wc -l
```

## Scripts added or changed

- `scripts/pipeline/run_7a6o_variant_direct.sh`
  - Minimal single-variant runner used for the live jobs.
  - Avoids backend probing and internal background scheduling.
  - Uses stable GROMACS flags: `-nb gpu -pme gpu -update cpu`.
  - Starts from `solv_ions_em_refined.gro` if present.

- `scripts/pipeline/refine_7a6o_mutant_em.sh`
  - Batch helper for refined unrestrained EM.
  - Useful for future reruns, though the live refined structures were finished before MD launch.

- `scripts/pipeline/run_7a6o_autoinhib_md_batch.sh`
  - Now prefers `solv_ions_em_refined.gro`.
  - Disables CUDA Graph and GPU update by default.
  - Propagates child failures instead of returning success after failed children.

- `scripts/pipeline/launch_7a6o_mutants_after_wt.sh`
  - WT completion check now accepts resumed WT output: `md_prod_run.part*.gro`.

- `scripts/pipeline/watch_7a6o_md_status.sh`
  - WT status now shows `md_prod_run.*` continuation files.

## Failure modes found

1. The original launcher only accepted `md_prod.gro` for WT completion, but WT completed as `md_prod_run.part0002.gro` after resume.
2. Mutants were previously using `solv_ions_em.gro` from a short, restrained EM. That left high local forces around side chains and caused LINCS failures.
3. GROMACS 2025 on this system automatically used GPU update/constraints unless forced. `-update gpu`/CUDA Graph triggered CUDA illegal-address failures on unstable starts.
4. Plain background/nohup launches from this tool shell were not reliable; child `grompp` processes were observed to die after writing only the GROMACS header. `setsid ... < /dev/null > log 2>&1 &` is the reliable launch pattern here.

## Safe exit

It is safe to exit the terminal after confirming `pgrep` shows `run_7a6o_variant_direct` or `gmx mdrun` processes. These jobs were launched with `setsid`, not attached to an interactive shell.

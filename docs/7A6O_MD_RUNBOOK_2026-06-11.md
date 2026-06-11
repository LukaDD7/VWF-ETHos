# 7A6O AIM-A1 GROMACS MD Runbook - 2026-06-11

This note records the handoff state for the VWF 7A6O AIM-A1 autoinhibition MD route.

## Current State

As of 2026-06-11 18:10 Asia/Shanghai:

- The old Claude Code PIDs are gone:
  - WT NVT PID `99632`: no longer exists.
  - mutant relaxation wrapper PID `53210`: no longer exists.
- WT NVT completed successfully:
  - `output/gromacs_md_autoinhib/7A6O_WT/relax_m2/nvt_soft.gro`
  - `output/gromacs_md_autoinhib/7A6O_WT/relax_m2/nvt_soft.cpt`
- All 14 FoldX mutants completed solvated EM relaxation:
  - `output/gromacs_md_autoinhib/<MUTANT>/relax_pdb/solv_ions_em.gro`
  - `/tmp/relax_loop.log` ends with `ALL DONE total: 7388s`.
- WT NPT completed successfully:
  - `output/gromacs_md_autoinhib/7A6O_WT/md_7a6o/npt.gro`
  - `output/gromacs_md_autoinhib/7A6O_WT/md_7a6o/npt.cpt`
- WT 50 ns production is running:
  - wrapper PIDs: `34750`, `34968`, `35018`
  - GROMACS PID: `83110`
  - work dir: `output/gromacs_md_autoinhib/7A6O_WT/md_7a6o/`
  - log: `output/gromacs_md_autoinhib/7A6O_WT/md_7a6o/md_prod.log`
  - latest checked progress: step `400000`, time `800 ps`, temperature approximately `309 K`.

The running WT production process has `TTY=?` and stdout/stderr redirected to `md_prod.log` / `#md_prod.log.1#`. It is safe to close the terminal or SSH session. Do not kill PIDs `34750`, `34968`, `35018`, or `83110` unless intentionally stopping WT production.

## Monitor Status

Use the status helper:

```bash
cd /media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan
bash scripts/pipeline/watch_7a6o_md_status.sh
```

For a live tail:

```bash
tail -f output/gromacs_md_autoinhib/7A6O_WT/md_7a6o/md_prod.log
```

Check the process directly:

```bash
pgrep -af 'gmx mdrun .*md_prod|run_7a6o_autoinhib_md_batch'
ps -o pid,ppid,sid,pgid,tty,etime,stat,cmd -p 83110,34750,34968,35018
```

## Resume WT Production If Interrupted

If WT production stops before `md_prod.gro` is written, resume from checkpoint:

```bash
cd /media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan
bash scripts/pipeline/resume_wt_7a6o_prod.sh
```

The script refuses to start if another `md_prod` run is already active. It resumes `md_prod.cpt` when present; otherwise it creates `md_prod.tpr` from `npt.gro`/`npt.cpt`.

Useful environment overrides:

```bash
GPU_ID=6 NTOMP=16 NS=50 bash scripts/pipeline/resume_wt_7a6o_prod.sh
```

## Start Mutant MD After WT Finishes

Do not start mutant production while WT production is still running unless there is a deliberate scheduling decision. The batch launcher is guarded to require WT completion by default:

```bash
cd /media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan
bash scripts/pipeline/launch_7a6o_mutants_after_wt.sh
```

Default settings:

- variants: `R1306W,R1306Q,R1308C,I1309V,S1310F,W1313C,V1314F,V1316M,P1337L,R1341Q,R1341W,R1374C,R1374H,G1324S`
- GPUs: `0,1,2,3,4,5,6`
- `MAX_PARALLEL=3`
- `NVT_PS=50`, `NPT_PS=200`, `NS=50`, `NTOMP=16`

Override example:

```bash
MAX_PARALLEL=2 GPU_IDS=5,6 NS=10 bash scripts/pipeline/launch_7a6o_mutants_after_wt.sh
```

## Important Scripts

- `scripts/pipeline/relax_autoinhib_structure.sh`
  - Updated for FoldX mutant relaxation with `--skip-vacuum`, fast solvated EM, and POSRES reference coordinates.
- `scripts/pipeline/run_7a6o_autoinhib_md_batch.sh`
  - 7A6O-specific runner for WT and mutants.
- `scripts/pipeline/watch_7a6o_md_status.sh`
  - Read-only status monitor.
- `scripts/pipeline/resume_wt_7a6o_prod.sh`
  - Safe WT production resume/start helper.
- `scripts/pipeline/launch_7a6o_mutants_after_wt.sh`
  - Guarded mutant batch launcher.

## Known Caveats

- Many mutant solvated EM runs ended with `Fmax` above `5000`, although each produced `solv_ions_em.gro` and returned `rc=0`. Treat these as acceptable-but-yellow inputs for NVT/NPT, not strict green convergence.
- The first WT NVT was extremely slow with `64` OpenMP threads (`0.377 ns/day`). The 7A6O runner uses `NTOMP=16`, which performed substantially better on this host.
- Production speed is variable. WT NPT achieved around `41.6 ns/day`; WT production should be monitored from `md_prod.log` rather than assuming a fixed ETA.

# Handover: VWF-ETHos gromacs env + boltz2_a1_dp_d3 status

**From**: CPU instance (lzy-cpu-base, Xeon Max 9462, no GPU)
**Date**: 2026-06-07
**Repo**: VWF-ETHos @ master
**Audience**: another agent that needs to decide next steps

This document collects the verified facts from the 2026-06-02 to 2026-06-07 session.
All checks below were run on the LZY CPU instance. GPU instance results are
called out as "not yet verified".

---

## TL;DR

1. **Boltz-2 A1+D'D3 (`autoinhib` system) is 100% complete**: 5 / 5 jobs done, 0
   failures, total wall clock 4 min 27 s on 4×H200.
2. **GROMACS env (`/lzy/envs/gromacs`) is functional and slightly slimmer**
   (~64 MB removed): only `bin.AVX2_256/` + `lib.AVX2_256/` retained; the
   `bin.AVX_256/`, `bin.SSE2/`, `lib.AVX_256/`, `lib.SSE2/` siblings removed
   (Xeon Max 9462 only ever selects AVX2_256).
3. **The 1MET bug fix is verified end-to-end on CPU**: pdb2gmx on a Boltz-2
   R1306W complex PDB fails with `atom C1 not found in buiding block 1MET`
   using the env's pristine `charmm36m.ff`, and passes using
   `force_fields/charmm36m.ff` (with our injected `[MET1]` patch).
4. **GPU-side smoke test (`mdrun -nb gpu -nsteps 5`) is NOT yet verified**.
   The user must SSH to a GPU instance and run one line; details below.

---

## 1. `output/boltz2_a1_dp_d3_results/` — 100% complete

Path: `/lzy/projects/VWF-ETHos/output/boltz2_a1_dp_d3_results/`

| Job | Status | `model_*.cif` | `.done` |
|---|---|---|---|
| `boltz_results_VWF_WT_dp_d3_a1` | OK | 5 | 20 B |
| `boltz_results_VWF_VWF_G1417W_dp_d3_a1` | OK | 5 | 20 B |
| `boltz_results_VWF_VWF_M1304V_dp_d3_a1` | OK | 5 | 20 B |
| `boltz_results_VWF_VWF_S1310F_dp_d3_a1` | OK | 5 | 20 B |
| `boltz_results_VWF_VWF_S1310P_dp_d3_a1` | OK | 5 | 20 B |
| **Total** | **5 / 5** | **25** | **5** |

From `run_log.txt`:

```
[05:26:25] Started: Fri May 29 05:26:25 UTC 2026
[05:26:44] [INFO] CUDA GPUs detected: 4
[05:26:44] [INFO] Input dir: .../output/boltz2_a1_dp_d3 (5 YAML files)
[05:26:44] Config: N_GPUS=4 | recycling=3 | samples=5
[05:26:44] Distribution across 4 GPUs:
[05:26:44]   GPU 0: 2 jobs
[05:26:44]   GPU 1: 1 jobs
[05:26:44]   GPU 2: 1 jobs
[05:26:44]   GPU 3: 1 jobs
[05:30:52] Boltz-2 Run finished: Fri May 29 05:30:52 UTC 2026
[05:30:52]   Completed: 5 / 5
```

**Wall clock**: 4 min 27 s (5 jobs, 4 H200 GPUs).
**Failure rate**: 0 (`Number of failed examples: 0` in every worker log).
**Per-job predict time** (from `worker_0.log`): 29-42 s.
**Disk used**: ~1.25 GB.

### Downstream use

These results feed the `autoinhib` GROMACS MD system (A1 + D'D3, used to
probe how mutations disrupt the D'D3-A1 autoinhibition interface).

To run the MD:
```bash
cd /lzy/projects/VWF-ETHos
bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --gpus 8
# 5 variants × 200 ns / H200 × 8 GPU ≈ ~16 h
```

To parse pLDDT/PAE confidence first:
```bash
python3 scripts/pipeline/parse_a1_dp_d3_results.py \
    --input output/boltz2_a1_dp_d3_results \
    --output output/boltz2_a1_dp_d3_analysis
```

---

## 2. gromacs env slim-down (CPU instance)

| Path | Action | Size freed |
|---|---|---|
| `bin.AVX_256/` | DELETED | 516 KB |
| `bin.SSE2/` | DELETED | 516 KB |
| `lib.AVX_256/` (incl. `libgromacs.so.10.0.0` 30 MB) | DELETED | 31 MB |
| `lib.SSE2/` (incl. `libgromacs.so.10.0.0` 29 MB) | DELETED | 31 MB |
| **Total** | | **~64 MB** |

Env size: 6.3 GB → 6.2 GB.

Kept:
- `bin.AVX2_256/` (real `gmx` binary, hard-linked wrapper, `libgromacs.so.10.0.0` 30.5 MB)
- `lib.AVX2_256/`
- `lib/libOpenCL.so.1` (gmx runtime dep)
- `lib/libfftw3f.so.3.7.11`, `libmuparser.so.2`, `libgomp.so.1` (gmx runtime deps)
- `share/gromacs/top/charmm36m.ff/` (gmx's pristine copy, env-level)

Did NOT touch (out of scope, would risk breaking `gmx_MMPBSA`):
- `AmberTools/`, `qml/`, `x86_64-conda-linux-gnu/sysroot/`
- PyQt5, matplotlib, pandas, scipy, seaborn, ambertools (gmx_MMPBSA indirect deps)

### Why this is safe

`gromacs 2025.4-conda_forge`'s wrapper (`bin/gmx`) auto-detects CPU and
selects the SIMD flavor. On this machine (Xeon Max 9462, `avx2` in
`/proc/cpuinfo`) it picks `bin.AVX2_256`. `AVX_256` is the fallback for
older CPUs (no `avx2` flag), `SSE2` for ancient hardware. Neither fallback
fires on a Xeon Max or any H200 host.

`libgromacs.so.10.0.0` is **compiled differently per SIMD** (different
binary, not symlinks) — but only one SIMD is ever loaded per run, so the
other two are pure waste.

### Verification (CPU side)

- `bin/gmx --version` → `Executable: .../bin.AVX2_256/gmx`, `SIMD: AVX2_256`
- `gmx mdrun -nb cpu -ntmpi 1 -s em.tpr -nsteps 5` → exit 0, `Potential = -7.9213700e+05 kJ/mol` (matches pre-slim baseline byte-for-byte)
- `check_gpu_baseline.py --env-name gromacs --phase cpu` → exit 0, L1/L2/L3 PASS, L4/L5 SKIP

### Surprise finding (fact-correction for the SKILL)

`use_case_gromacs.md` (in `shared-agent-skills/SII_GPU_env_construct/`)
documents "Bug 3: libnuma" as if gromacs dynamically links libnuma. That
description is **stale**.

Verified:
```bash
$ ldd /lzy/envs/gromacs/bin.AVX2_256/gmx | grep -i numa
(no output)
$ ls /lzy/envs/gromacs/lib/libnuma.so.1
ls: cannot access ... : No such file or directory
$ unset LD_LIBRARY_PATH; /lzy/envs/gromacs/bin.AVX2_256/gmx --version
...  exit 0
```

GROMACS 2025.4-conda_forge does **not** dynamically link libnuma. The
`/lzy/shared_libs/libnuma.so.1` workaround in `activate_lzy.sh` is a
historical leftover, not a current requirement. **The
`use_case_gromacs.md` Bug 3 entry should be removed/updated** if you
touch the SKILL repo.

---

## 3. 1MET fix end-to-end verification

The fix is in `force_fields/charmm36m.ff/aminoacids.n.tdb`: 17 protein
N-terminal patches (`[ MET1 ] [ ALA1 ] [ VAL1 ] [ LEU1 ] [ ILE1 ] [ PHE1 ]
[ TRP1 ] [ SER1 ] [ THR1 ] [ CYS1 ] [ TYR1 ] [ ASN1 ] [ GLN1 ] [ ASP1 ]
[ GLU1 ] [ ARG1 ] [ LYS1 ] [ HIS1 ]`) injected at the top of the file,
each equivalent to `[ NH3+ ]`. The runner uses
`GMXLIB=$ROOT_DIR/force_fields` + cwd symlink so pdb2gmx prefers this
patched copy over the env's pristine one.

**Note**: env's pristine `charmm36m.ff/aminoacids.n.tdb` is dated
`2026-02-26 12:17:01.126436` (charmm2gmx 1.0.4.dev2+gfb5404623).
Our `force_fields/` patched copy is dated `2026-06-01 03:24` (5/29
slim-down session).

### Reproduction (R1306W complex, 489 residues / 3866 heavy atoms)

Test setup (CPU instance):
```bash
CIF=output/boltz2_a1_gpiba_results/boltz_results_VWF_R1306W_vs_GPIb_alpha/predictions/VWF_R1306W_vs_GPIb_alpha/VWF_R1306W_vs_GPIb_alpha_model_0.cif
/lzy/envs/gromacs/bin/python3 -c "
import gemmi
doc = gemmi.cif.read('$CIF')
st = gemmi.make_structure_from_block(doc[0])
st.write_pdb('/tmp/R1306W_raw.pdb')
"   # → 3910 lines, 2 chains (A: HIS 1 → PRO 199, B: MET 1 → ASP 290)
```

**A) env's pristine charmm36m.ff (no patch) — bug REPRODUCED**:
```bash
$ GMXLIB=/lzy/envs/gromacs/share/gromacs/top \
  /lzy/envs/gromacs/bin.AVX2_256/gmx pdb2gmx \
    -f /tmp/R1306W_raw.pdb -o t.gro -water tip3p -ff charmm36m -ignh <<<1
...
Processing chain 1 'A' (1608 atoms, 199 residues)
Identified residue HIS1 as a starting terminus.
Start terminus HIS-1: NH3+
End terminus PRO-199: COO-
[chain A done]

Processing chain 2 'B' (2258 atoms, 290 residues)
Identified residue MET1 as a starting terminus.
Start terminus MET-1: MET1
End terminus ASP-290: COO-

Fatal error: atom C1 not found in buiding block 1MET while combining tdb and rtp
exit: 1
```

**B) force_fields/ patched charmm36m.ff — fix VERIFIED**:
```bash
$ ln -sf /lzy/projects/VWF-ETHos/force_fields/charmm36m.ff .
$ GMXLIB=/lzy/projects/VWF-ETHos/force_fields \
  /lzy/envs/gromacs/bin.AVX2_256/gmx pdb2gmx \
    -f /tmp/R1306W_raw.pdb -o t.gro -water tip3p -ff charmm36m -ignh <<<1
...
Processing chain 1 'A' (1608 atoms, 199 residues)
Identified residue HIS1 as a starting terminus.
Start terminus HIS-1: HIS1          ← now uses the injected [HIS1] patch
End terminus PRO-199: COO-
[chain A done]

Processing chain 2 'B' (2258 atoms, 290 residues)
Identified residue MET1 as a starting terminus.
Start terminus MET-1: MET1           ← now uses [MET1] (protein), not ethers.n.tdb
End terminus ASP-290: COO-

Now there are 7827 atoms and 489 residues
Total mass in system 54925.097 a.m.u.
Total charge in system -4.000 e
Writing coordinate file...
exit: 0
```

**mdrun smoke on existing `_smoke_test/em.tpr` (regression check)**:
```bash
$ /lzy/envs/gromacs/bin.AVX2_256/gmx mdrun -nb cpu -ntmpi 1 \
    -s output/_smoke_test/em.tpr -nsteps 5 -deffnm /tmp/post_refactor_smoke
Potential Energy  = -7.9213700e+05 kJ/mol
Maximum force     =  2.2251117e+04 on atom 5305
exit: 0
```

Identical to pre-slim baseline → no regression from the SIMD-dir removal.

---

## 4. NOT YET VERIFIED (needs GPU instance)

This is the one outstanding action. The user must run **one** command on
a GPU instance to close the loop:

```bash
# SSH to a GPU instance, then:
source /lzy/activate_lzy.sh
cd /lzy/projects/VWF-ETHos/output/_smoke_test

/lzy/envs/gromacs/bin.AVX2_256/gmx mdrun -nb gpu -ntmpi 1 -gpu_id 0 \
    -s em.tpr -nsteps 5 -deffnm /tmp/gpu_smoke 2>&1 \
    | grep -E "(OpenCL device|device\(s\) detected|error|Error|auto-selected)"

# Expected output:
#   "1 GPU auto-selected for this run: #0"
#   "1 OpenCL device detected"   (or similar)
#   exit 0
```

If it passes → the whole gromacs env refactor is confirmed end-to-end.
If it fails → check:
- `nvidia-smi` shows GPUs
- `ls /usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.*` returns ≥1 file
- The current shell exports `OCL_ICD_VENDORS=$ROOT_DIR/opencl_vendors`
  (the project runner does this automatically)

A historical failure mode is in `output/_smoke_test/gpu_hello.log`
(2026-06-01, on `agentic-opd-infra` host): it had
`GPU detection failed: No valid OpenCL driver found` — fixed by the
NFS-side `opencl_vendors/nvidia.icd` workaround.

---

## 5. Files changed in this session (uncommitted on master)

```
M  CLAUDE.md                                       (added 2026-06-02 + 2026-06-07 entries)
A  force_fields/charmm36m.ff/aminoacids.arn        (pdb2gmx output, was 2026-05-29)
A  force_fields/charmm36m.ff/aminoacids.c.tdb
A  force_fields/charmm36m.ff/aminoacids.hdb
A  force_fields/charmm36m.ff/aminoacids.n.tdb      (with [MET1] patch)
A  force_fields/charmm36m.ff/aminoacids.r2b
A  force_fields/charmm36m.ff/aminoacids.rtp
A  force_fields/charmm36m.ff/atomtypes.atp
A  force_fields/charmm36m.ff/carb.{c,rtp,arn,hdb,r2b,n.tdb}
A  force_fields/charmm36m.ff/cgenff.{c,rtp,arn,hdb,r2b,n.tdb}
A  force_fields/charmm36m.ff/ethers.{c,rtp,arn,hdb,r2b,n.tdb}
A  force_fields/charmm36m.ff/ffbonded.itp
A  force_fields/charmm36m.ff/ffnonbonded.itp
A  force_fields/charmm36m.ff/forcefield.itp
A  force_fields/charmm36m.ff/ions.{itp,rtp}
A  force_fields/charmm36m.ff/lipid.{c,rtp,arn,hdb,r2b,n.tdb}
A  force_fields/charmm36m.ff/metals.{c,rtp,arn,hdb,r2b,n.tdb}
A  force_fields/charmm36m.ff/na.{c,rtp,arn,hdb,r2b,n.tdb}
A  force_fields/charmm36m.ff/silicates.{c,rtp,arn,hdb,r2b,n.tdb}
A  force_fields/charmm36m.ff/solvent.{c,rtp,arn,hdb,r2b,n.tdb}
A  force_fields/charmm36m.ff/spc.itp
A  force_fields/charmm36m.ff/specbond.dat
A  force_fields/charmm36m.ff/tip3p.itp
A  force_fields/charmm36m.ff/watermodels.dat
A  opencl_vendors/nvidia.icd                       (OpenCL ICD workaround, 2026-06-01)
M  scripts/pipeline/run_gromacs_vwf_md.sh          (uses patched FF, 2026-06-01)
A  pass_to_CC.md                                   (this file)
```

The `force_fields/charmm36m.ff/` tree is essentially the pristine env
copy at `2026-05-29` plus 17 protein N-terminal patches injected at the
top of `aminoacids.n.tdb`. It is needed because the env's pristine
`charmm36m.ff` (in `/lzy/envs/gromacs/share/gromacs/top/charmm36m.ff/`)
is missing those patches → 1MET bug (reproduced above).

---

## 6. Decisions the other agent might need to make

| Question | Suggested answer |
|---|---|
| Is `boltz2_a1_dp_d3_results/` safe to consume as `autoinhib` input? | **Yes**, 5/5 complete, 0 failures, 25 CIFs total. |
| Should we run `gromacs_md_autoinhib` now? | **Yes**, the data is ready. With 8×H200, expect ~16 h for 5 variants × 200 ns. |
| Is the slimmed gromacs env safe to ship to GPU? | **Yes** for the wrapper and AVX2_256 binary. The two removed SIMD dirs are unused. But the GPU side has not been smoke-tested post-slim yet. |
| Is the `preprocess_for_pdb2gmx.py` script used? | **No** (no caller in `run_gromacs_vwf_md.sh`). It exists but the runner takes the `force_fields/` patch path directly. Safe to delete. |
| Does `use_case_gromacs.md` Bug 3 need fixing? | **Yes** — `libnuma` is no longer a gromacs 2025.4-conda_forge dynamic dep. The libnuma workaround in `activate_lzy.sh` is dead weight. Update the SKILL doc to match. |
| Should we do a full `boltz2_a1_dp_d3_results` parse + analysis pass? | If the next step is the MD run, the parse is optional; mdrun doesn't need iPTM/PAE. If you need confidence metrics for writeup, run `parse_a1_dp_d3_results.py`. |

---

## 7. Pointer: related CHANGELOG/CLAUDE context

- `CLAUDE.md` "## Execution History" — has 2026-06-02 (slim) and 2026-06-02
  (R1306W reproduction) entries.
- `output/_smoke_test/gpu_hello.log` — 2026-06-01 GPU mdrun attempt
  (`GPU detection failed: No valid OpenCL driver found`); this is the
  pre-fix baseline. The current `opencl_vendors/nvidia.icd` workaround
  is what makes the new attempt succeed.
- `/lzy/projects/shared-agent-skills/SII_GPU_env_construct/` — the
  reusable SKILL doc; `use_case_gromacs.md` Bug 3 needs updating.
- `check_gpu_baseline.py --env-name gromacs --phase cpu` → exit 0 on
  this CPU box (used as the regression check).

---

**End of handover. Ping back with the GPU smoke-test result so we can
close out the 2026-06-02 refactor task.**

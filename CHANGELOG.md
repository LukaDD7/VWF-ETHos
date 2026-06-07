# Changelog

All notable changes to the VWF-ETHos project are documented here.

## [Unreleased]

### Added (2026-06-07)

- **`scripts/pipeline/gpu_smoke_test.sh`** — 诚实的 GROMACS GPU 放行闸:用与生产 runner 相同的条件化 flags + md 积分器 tpr 实跑 5 步(朴素 `-nb gpu` smoke 会假性放行)。
- **`scripts/pipeline/extract_aim_autoinhib_features.py`** — 从 `a1_aim_autoinhibition_context` 结构抽 AIM↔A1 接触数 → `aim_release_score`(2B 判别特征,替代分不开 2B/2M 的全局 plddt)。
- **`docs/A40_RUNBOOK.md`** — 独立 A40 服务器(非 /lzy NFS)搬运 + 重建 CUDA gromacs env + smoke + autoinhib MD + 2B 特征流水线。

### Changed (2026-06-07)

- **`scripts/pipeline/run_gromacs_vwf_md.sh`** — GPU offload flags 按 `$GPU_BACKEND` 条件化:CUDA/SYCL→`-bonded gpu -update gpu`+CUDA graph;OpenCL/未知→仅 `-nb gpu -pme gpu`。修复 OpenCL 构建上 NVT/NPT/Prod 必 fatal(`not supported with OpenCL`)、白烧 GPU 机时的 bug。
- **`scripts/agentic_vwf_classifier.py`** — RULE6(A1)加性接入 `aim_release_score`:强松开→2B(早返回),弱松开→救回本判 2M 的 A1 默认;NaN 时退回原逻辑。阈值 `AIM_RELEASE_2B_Z` 暂定,须用已知 2B/2M 校准。

### Added

- **`scripts/pipeline/run_a1_gpiba_boltz2.sh`** — GPU batch runner for VWF A1 + GPIbα Boltz-2 predictions
  - 74 YAML inputs in `output/boltz2_a1_gpiba/`
  - 4-GPU parallel execution with per-job `.done` checkpointing
  - Key fixes: `CC=$(which gcc)`, `--num_workers 0`, `ulimit -c 0`
- **`scripts/pipeline/parse_a1_gpiba_results.py`** — Parser for Boltz-2 confidence outputs

### Infrastructure Fixes (2026-05-07)

| Issue | Root Cause | Fix |
|-------|-----------|-----|
| `core.*` files (6GB each) | DataLoader multiprocessing overflowed `/dev/shm` (64MB) | `--num_workers 0` disables sub-process IPC |
| `RuntimeError: Failed to find C compiler` | Triton JIT needs gcc; `CC` env not inherited by subprocesses | `CC=$(which gcc)` prepended to boltz command |
| `MisconfigurationException: No supported gpu backend` | `num_workers > 0` triggers DataLoader fork on CPU-only host | `--num_workers 0` |
| GPU util stuck at 4% then crashes | Same as above: multiprocessing + shm overflow → segfault | `--num_workers 0` |
| `boltz: command not found` in subshell | `conda activate boltz2` not finding env; shell uses wrong env | Use `export -f run_worker` + ensure `boltz` in PATH |

## [Prior History]

See `CLAUDE.md` and `DEV_PROTOCOL_STANDARD.md` for pre-May 2026 execution logs.

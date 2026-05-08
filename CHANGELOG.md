# Changelog

All notable changes to the VWF-ETHos project are documented here.

## [Unreleased]

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

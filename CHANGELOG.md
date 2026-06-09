# Changelog

All notable changes to the VWF-ETHos project are documented here.

## [Unreleased]

### Added (2026-06-09)

- **`docs/A40_ENV_SETUP_REPORT_2026-06-09.md`** — 7×A40 独立服务器(无 `/lzy` NFS)的 GROMACS env 部署完整报告:从 git pull 到 SMOKE PASS 全流程、3 个 runner bug 修复、4 个新发现坑(curl 假断网、Test A 的 em.tpr 问题、初始结构需 EM、em.tpr vs md.tpr 命名冲突)、目标程序(autoinhib MD)恢复路径。

### Fixed (2026-06-09)

- **`scripts/pipeline/run_gromacs_vwf_md.sh`** (commit `ba5f711`) — A40 本地 env 兼容:
  - `gmx` 路径检测: 加 `$ROOT_DIR/envs/gromacs/` 本地 fallback (A40_AGENT_SETUP 路径)
  - `GMX_PY` 检测: 同样多路 fallback, 修 `$(dirname $GMX)/../../bin/python` `..` 数错 bug
  - `GMXDATA` 检测: 同样本地 fallback, 避免 `cd 错误目录`
  - `log()` 函数: 从 line 153 提前到 line 70, 修"A40 上 gmx 检测段调 log 报 '未找到命令'"的 bug (LZY 上因 `command -v gmx` 命中 PATH 而幸运避开)

### Known Issues (待办, 2026-06-09)

- `a40_selfcheck.sh` 联网检查用 `curl -sI`, 老版 curl (7.x) 报 `curl: (48)` 假断网, 实际 tuna 镜像可达。修法: 用 `python3 -c "import urllib.request; ..."` 或 `wget` 替代。
- `gpu_smoke_test.sh` Test A 用 `em.tpr` (steep 积分器) 跑 `-pme gpu` 必失败 (GROMACS 2025.x 不支持非动力学积分器的 PME-on-GPU)。这不是 GPU 问题, 是脚本设计 bug。修法: 内部用 md 积分器 tpr 跑 Test A。
- `output/boltz2_a1_dp_d3_results/` (autoinhib MD 输入) 在 A40 上不存在, 需 scp 自 LZY。

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

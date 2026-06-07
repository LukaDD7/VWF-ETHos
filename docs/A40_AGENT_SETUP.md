# A40 Agent Setup Protocol — gromacs GPU env(自查 → 建 env → 验证)

> **这份是给 A40 服务器上的 Claude Code / Agent 看的执行协议**(不是给人读的教程,
> 教程见 [A40_RUNBOOK.md](A40_RUNBOOK.md))。把本文件路径直接丢给那台机的 Agent,
> 让它按 Phase 0→4 顺序执行,每个 Phase 把输出贴出来。

## 你的任务(给 A40 Agent)

在这台**独立的 A40 GPU 服务器**(不挂 `/lzy` NFS)上,搭一个**能跑 GPU 的 gromacs env**,
让 `scripts/pipeline/gpu_smoke_test.sh` 通过。目标下游是对 VWF 的自抑制模块跑 MD。
**不要跳过 smoke;smoke 不过不要跑批量(会烧 GPU 机时)。**

心智模型(5 层,出问题先定位是哪层):
```
L5 内核模块 nvidia.ko + /dev/nvidia*   ← 需 root, Agent 改不了 → 缺则 ESCALATE
L4 GPU 驱动 userspace + OpenCL ICD      ← .so 来自驱动; ICD 用项目 opencl_vendors/ 兜底
L3 系统 .so (libnuma 等)               ← 本项目 gromacs 不依赖 libnuma(已确认)
L2 conda env (gromacs/gemmi)           ← 本协议要建的东西
L1 项目代码 + force_fields/(端基修复)  ← 随 git, git pull 即得
```

---

## Phase 0 — 自查(只读,先跑这个)

```bash
git pull
bash scripts/setup/a40_selfcheck.sh
```

脚本会采集:GPU/驱动/CUDA 上限、/dev 节点、OpenCL ICD、conda、联网、/lzy 是否挂载、
repo 资产、磁盘,并在末尾打印**推荐的 conda 命令**和退出码:

- **exit 2 (ESCALATE)**:无 GPU/驱动/`/dev/nvidia*` → 需 root 装驱动 + `nvidia-modprobe`。
  **停下,报告用户找运维**,Agent 无 root 做不了。
- **exit 1 (BLOCKER)**:缺 conda / 不能联网 / 缺 repo 资产 → 先逐个解决(装 miniconda、
  联网、`git pull`),再重跑 Phase 0。
- **exit 0**:继续 Phase 1。

> 特例:若 Phase 0 报"`/lzy/envs/gromacs` 存在"(说明这台其实挂了 LZY NFS),
> **不要重建**,改走 [A40_RUNBOOK.md](A40_RUNBOOK.md) 的"同 NFS"路径,直接 Phase 3。

---

## Phase 1 — 决策(CUDA 还是 OpenCL 构建)

照搬 Phase 0 末尾 `RECOMMENDATION` 的命令。决策逻辑(脚本已替你判,这里是依据):

| 条件 | 选择 | 为什么 |
|---|---|---|
| 默认 / driver CUDA 上限够新 | **CUDA 构建** `gromacs=2025.*=nompi_cuda_*` | GPU 常驻(`-bonded/-update gpu`)提速,**绕开 OpenCL ICD 整摊事**,runner 全套 flags 可用 |
| CUDA build 在 channel 解析不到 / driver CUDA 太低 | OpenCL 构建 + `ocl-icd-libopencl1` | fallback;必须设 `OCL_ICD_VENDORS`(见 Phase 2b) |

> 关于"OpenCL 能不能在 conda create 时加":能——loader(`libOpenCL.so.1`)是 conda 包
> `ocl-icd-libopencl1`,会进 env;但 NVIDIA **vendor ICD** 来自 A40 驱动,装不进 env,
> 靠 `OCL_ICD_VENDORS=$ROOT/opencl_vendors` 指到项目自带的 `nvidia.icd` 兜底。

---

## Phase 2 — 建 env

```bash
ROOT=$PWD
CONDA=$(command -v mamba || command -v conda)

# 2a. CUDA 优先
$CONDA create -p $ROOT/envs/gromacs python=3.11 -y
$CONDA install -p $ROOT/envs/gromacs -c conda-forge -y \
    "gromacs=2025.*=nompi_cuda_*" gemmi numpy pandas

# 2b. 仅当 CUDA build 失败 → OpenCL fallback
# $CONDA install -p $ROOT/envs/gromacs -c conda-forge -y \
#     gromacs=2025 gemmi numpy pandas ocl-icd ocl-icd-libopencl1
# export OCL_ICD_VENDORS=$ROOT/opencl_vendors

# 2c. 确认后端(决定 runner 用哪套 flags;runner/smoke 会自动按此条件化)
GMX=$ROOT/envs/gromacs/bin.AVX2_256/gmx
[ -x "$GMX" ] || GMX=$ROOT/envs/gromacs/bin/gmx
$GMX --version | grep -iE "GPU support|SIMD"      # ★记下: CUDA(理想) 或 OpenCL
```

---

## Phase 3 — 验证(放行闸,必须 PASS)

```bash
export GMX GMXLIB=$ROOT/force_fields
# OpenCL 构建额外: export OCL_ICD_VENDORS=$ROOT/opencl_vendors

# smoke 需要一个小 EM 体系; 若 output/_smoke_test/ 里没有, 用 WT autoinhib 造一个:
if [ ! -f output/_smoke_test/em.tpr ]; then
    bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --filter WT --phase em --gpus 1
    WTDIR=$(ls -d output/gromacs_md_autoinhib/*WT*/ | head -1)
    mkdir -p output/_smoke_test
    ln -sf "$PWD/$WTDIR/em/em.tpr"          output/_smoke_test/em.tpr
    ln -sf "$PWD/$WTDIR/em/em.gro"          output/_smoke_test/em.gro
    ln -sf "$PWD/$WTDIR/topology/topol.top" output/_smoke_test/topol.top
fi

bash scripts/pipeline/gpu_smoke_test.sh
# Test A = GPU 是否被发现; Test B = 用生产 runner 真实 flags 实跑(关键)
# 期望: "✅ SMOKE PASS"
```

若 Test B 报 `not supported with OpenCL`:说明是 OpenCL 构建——runner 已会自动降级,
但请确认 `gmx --version` 的 `GPU support:` 被正确识别。若想要全速,回 Phase 1 改建 CUDA。

---

## Phase 4 — 回报(给主对话/用户)

把这几项贴回:
1. Phase 0 自查的 GPU 型号 / driver / CUDA 上限 / exit code;
2. 最终建成的 **gromacs 后端**(`GPU support:` = CUDA 还是 OpenCL);
3. Phase 3 `gpu_smoke_test.sh` 结果(PASS / 在哪步 FAIL + 关键日志);
4. 是否已 scp `output/boltz2_a1_dp_d3_results/`(autoinhib MD 输入)。

smoke PASS 后,下游(autoinhib MD + 2B 特征)按 [A40_RUNBOOK.md](A40_RUNBOOK.md) §5–6。

## 停止条件 / 不要做的事
- smoke PASS 即停止环境折腾(别继续调,避免破坏可用状态)。
- 不改 `/etc/`、不动系统驱动(L4/L5 需 root → ESCALATE 给用户)。
- 不把 NVIDIA 专有 .so 拷进 conda env(EULA + ABI 风险);ICD 用文本 `.icd` + 环境变量。

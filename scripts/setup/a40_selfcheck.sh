#!/bin/bash
# ==============================================================================
# a40_selfcheck.sh — 独立 GPU 机自查 + gromacs env 构建推荐 (只读)
# ==============================================================================
# 给目标 GPU 服务器上的 Agent 用: 一键采集机器实况, 判断该建 CUDA 还是 OpenCL
# 版 gromacs, 并打印可直接执行的 conda 命令。本脚本不安装、不改任何系统文件。
#
# 用法:  bash scripts/setup/a40_selfcheck.sh
# 退出码: 0 = 可继续建 env; 1 = 有阻断项(需先处理, 见 [BLOCKER]); 2 = 无 GPU/驱动(需 root)
# ==============================================================================
set -u
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

say() { printf '%s\n' "$*"; }
hdr() { printf '\n== %s ==\n' "$*"; }

BLOCK=0; NO_GPU=0
REC_BUILD="cuda"   # cuda | opencl

say "============================================================"
say " A40 / 独立 GPU 机自查  (VWF-ETHos gromacs env)"
say " repo: $ROOT_DIR"
say " host: $(hostname 2>/dev/null)   date: $(date '+%F %T')"
say "============================================================"

# ---- 1. OS / glibc -----------------------------------------------------------
hdr "1. OS / glibc"
say "  kernel : $(uname -r 2>/dev/null)"
say "  glibc  : $(ldd --version 2>/dev/null | head -1)"

# ---- 2. GPU / 驱动 / CUDA 上限 (Layer 4/5) -----------------------------------
hdr "2. GPU / driver / CUDA"
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader 2>/dev/null \
        | sed 's/^/  GPU: /' || say "  [WARN] nvidia-smi 查询失败"
    DRV_CUDA=$(nvidia-smi 2>/dev/null | grep -oE "CUDA Version: [0-9]+\.[0-9]+" | awk '{print $3}')
    say "  driver 支持的 CUDA 上限: ${DRV_CUDA:-未知}"
else
    say "  [BLOCKER] 没有 nvidia-smi → 无驱动/无 GPU。需要 root 装 NVIDIA 驱动。"
    NO_GPU=1
fi
DEVNODES=$(ls /dev/nvidia0 /dev/nvidiactl /dev/nvidia-uvm 2>/dev/null | tr '\n' ' ')
if [ -n "$DEVNODES" ]; then say "  /dev nodes: $DEVNODES"; else
    say "  [BLOCKER] 缺 /dev/nvidia* 设备节点 → 内核模块没加载 (nvidia-modprobe, 需 root)。"
    NO_GPU=1
fi

# ---- 3. OpenCL ICD (仅 OpenCL 构建才需要) ------------------------------------
hdr "3. OpenCL ICD (OpenCL 构建 fallback 才用)"
ICD_SO=$(ls /usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.* 2>/dev/null | head -1)
say "  driver OpenCL .so : ${ICD_SO:-未找到}"
ETC_ICD=$(ls /etc/OpenCL/vendors/*.icd 2>/dev/null | tr '\n' ' ')
say "  /etc/OpenCL/vendors: ${ETC_ICD:-(空, 用项目 opencl_vendors/ + OCL_ICD_VENDORS 兜底)}"

# ---- 4. conda / 联网 (建 env 前提) -------------------------------------------
hdr "4. conda / 联网"
CONDA_BIN=$(command -v mamba 2>/dev/null || command -v conda 2>/dev/null)
if [ -n "$CONDA_BIN" ]; then say "  conda/mamba: $CONDA_BIN ($($CONDA_BIN --version 2>/dev/null))"; else
    say "  [BLOCKER] 没有 conda/mamba。先装 miniconda/mambaforge 再继续。"; BLOCK=1
fi
if timeout 6 curl -sI https://conda.anaconda.org >/dev/null 2>&1; then say "  联网: OK (能到 conda-forge)"; else
    say "  [BLOCKER] 连不上 conda-forge → 无法 conda install。这台机要能联网, 或换 conda-pack 离线包。"; BLOCK=1
fi

# ---- 5. /lzy NFS? (判断是否其实能复用现成 env) -------------------------------
hdr "5. /lzy NFS"
if [ -d /lzy/envs/gromacs ]; then
    say "  [INFO] /lzy/envs/gromacs 存在 → 这台其实挂了 LZY NFS, 可直接复用, 无需重建。"
    say "         (那就改用 docs/A40_RUNBOOK.md 的 '同 NFS' 路径, 跳过本建 env 流程)"
else
    say "  /lzy 未挂载 → 按独立机重建 env (预期路径)。"
fi

# ---- 6. repo 内随行资产 ------------------------------------------------------
hdr "6. repo 资产 (随 git 走, 应已就位)"
for p in force_fields/charmm36m.ff/aminoacids.n.tdb opencl_vendors/nvidia.icd \
         scripts/pipeline/run_gromacs_vwf_md.sh scripts/pipeline/gpu_smoke_test.sh; do
    [ -e "$ROOT_DIR/$p" ] && say "  [OK]  $p" || { say "  [BLOCKER] 缺 $p (先 git pull)"; BLOCK=1; }
done
# autoinhib MD 输入 (需 scp, 不在 git)
if ls -d "$ROOT_DIR"/output/boltz2_a1_dp_d3_results/boltz_results_VWF_* >/dev/null 2>&1; then
    N=$(ls -d "$ROOT_DIR"/output/boltz2_a1_dp_d3_results/boltz_results_VWF_* | wc -l | tr -d ' ')
    say "  [OK]  autoinhib 输入结构: $N 个"
else
    say "  [WARN] 缺 output/boltz2_a1_dp_d3_results/ (autoinhib MD 输入, 需从 LZY scp; smoke 可先不用)"
fi

# ---- 7. 磁盘 -----------------------------------------------------------------
hdr "7. 磁盘"
df -h "$ROOT_DIR" 2>/dev/null | awk 'NR==1||NR==2{print "  "$0}'
say "  (env ~6GB + 每个 200ns autoinhib MD ~3-5GB)"

# ---- 决策 --------------------------------------------------------------------
hdr "推荐 (RECOMMENDATION)"
ENVP="$ROOT_DIR/envs/gromacs"
if [ "$NO_GPU" -eq 1 ]; then
    say "  [ESCALATE] 无 GPU/驱动/设备节点 → 需要 root 装 NVIDIA 驱动 + nvidia-modprobe。"
    say "             这一步 Agent 在无 root 时做不了, 必须先让运维处理。"
    exit 2
fi
if [ "$BLOCK" -eq 1 ]; then
    say "  [BLOCKER] 上面有阻断项, 先逐个解决再回到这里。"
    exit 1
fi

say "  建 env (CUDA 优先 — GPU 常驻提速 + 绕开 OpenCL ICD):"
say "    $CONDA_BIN create -p $ENVP python=3.11 -y"
say "    $CONDA_BIN install -p $ENVP -c conda-forge -y \\"
say "        \"gromacs=2025.*=nompi_cuda_*\" gemmi numpy pandas"
say ""
say "  若上面 CUDA build 解析不到 (driver CUDA 上限=${DRV_CUDA:-?} 过低 / channel 无该 build):"
say "    $CONDA_BIN install -p $ENVP -c conda-forge -y \\"
say "        gromacs=2025 gemmi numpy pandas ocl-icd ocl-icd-libopencl1"
say "    # 然后 OpenCL 兜底: export OCL_ICD_VENDORS=$ROOT_DIR/opencl_vendors"
say ""
say "  建完务必确认后端 (决定 runner flags):"
say "    $ENVP/bin.AVX2_256/gmx --version | grep -iE 'GPU support|SIMD'   # 期望 CUDA(理想)/OpenCL"
say ""
say "  然后验证 (必须 PASS 才上批量):"
say "    export GMX=$ENVP/bin.AVX2_256/gmx GMXLIB=$ROOT_DIR/force_fields"
say "    bash scripts/pipeline/gpu_smoke_test.sh"
say ""
say "  [OK] 自查通过, 可继续建 env。"
exit 0

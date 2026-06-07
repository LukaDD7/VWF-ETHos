#!/bin/bash
# ==============================================================================
# gpu_smoke_test.sh — GROMACS GPU 放行闸 (honest go/no-go)
# ==============================================================================
#
# 为什么需要它:
#   一个朴素的 `mdrun -nb gpu -nsteps 5` 只测 nonbonded offload。但生产 runner
#   (run_gromacs_vwf_md.sh) 的 NVT/NPT/Production 段还会用 GPU 常驻 flags
#   (-bonded gpu -update gpu, CUDA graph)。这些在 OpenCL 构建上会 fatal。
#   所以"朴素 smoke 过了"会假性放行 → 真跑到 NVT 才挂、白烧 GPU 机时。
#
#   本脚本用和 runner **完全相同**的按后端条件化逻辑构造 flags, 并用 md 积分器
#   的 tpr 实跑 5 步 —— 做到 "本脚本过 ⟹ 生产 runner 的 mdrun flags 一定过"。
#
# 用法 (在 GPU 实例上):
#   source /lzy/activate_lzy.sh      # 或手动 export GMX / OCL_ICD_VENDORS
#   bash scripts/pipeline/gpu_smoke_test.sh
#   bash scripts/pipeline/gpu_smoke_test.sh --gmx /path/to/gmx --gpu-id 0
#
# 退出码: 0 = 放行 (A 和 B 都过); 1 = 不放行 (含定位提示)
# ==============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# ---------- 参数 --------------------------------------------------------------
GMX="${GMX:-}"
SMOKE_DIR="${SMOKE_DIR:-$ROOT_DIR/output/_smoke_test}"
GPU_ID=0
KEEP=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --gmx)       GMX="$2"; shift 2 ;;
        --smoke-dir) SMOKE_DIR="$2"; shift 2 ;;
        --gpu-id)    GPU_ID="$2"; shift 2 ;;
        --keep)      KEEP=true; shift ;;
        -h|--help)   sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[WARN] unknown arg: $1"; shift ;;
    esac
done

# ---------- 定位 gmx ----------------------------------------------------------
if [ -z "$GMX" ]; then
    for cand in /lzy/envs/gromacs/bin.AVX2_256/gmx \
                "$ROOT_DIR/../../envs/gromacs/bin.AVX2_256/gmx" \
                "$(command -v gmx 2>/dev/null)"; do
        [ -n "$cand" ] && [ -x "$cand" ] && { GMX="$cand"; break; }
    done
fi
[ -z "$GMX" ] || [ ! -x "$GMX" ] && { echo "[FATAL] 找不到 gmx。用 --gmx 指定。"; exit 1; }

# ICD 兜底 (OpenCL 构建需要; CUDA 构建无害)
export OCL_ICD_VENDORS="${OCL_ICD_VENDORS:-$ROOT_DIR/opencl_vendors}"
export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"

echo "============================================================"
echo "GROMACS GPU smoke test"
echo "  GMX       : $GMX"
echo "  smoke dir : $SMOKE_DIR"
echo "  GPU id    : $GPU_ID"
echo "  OCL_ICD_VENDORS=$OCL_ICD_VENDORS"
echo "============================================================"

# ---------- 后端检测 + 条件化 flags (与 run_gromacs_vwf_md.sh 同逻辑) ----------
GPU_BACKEND=$("$GMX" mdrun -version 2>&1 | grep "GPU support:" | awk '{print $NF}')
echo "[INFO] GPU backend: ${GPU_BACKEND:-unknown}"
case "$GPU_BACKEND" in
    CUDA|SYCL)
        GPU_RESIDENT_FLAGS="-bonded gpu -update gpu"
        export GMX_CUDA_GRAPH=1
        echo "[INFO] GPU-resident mode ON: -nb gpu -pme gpu $GPU_RESIDENT_FLAGS (CUDA graph=1)"
        ;;
    *)
        GPU_RESIDENT_FLAGS=""
        unset GMX_CUDA_GRAPH 2>/dev/null || true
        echo "[INFO] degraded mode: -nb gpu -pme gpu (backend 不支持 GPU 常驻)"
        ;;
esac

WORK="$(mktemp -d "${TMPDIR:-/tmp}/gmx_smoke.XXXXXX")"
$KEEP || trap 'rm -rf "$WORK"' EXIT
FAIL=0

fail() { echo "[FAIL] $*"; FAIL=1; }

# ---------- 0. nvidia-smi / ICD 基础 ------------------------------------------
if ! nvidia-smi -L >/dev/null 2>&1; then
    fail "nvidia-smi 看不到 GPU (Layer 5: 驱动/内核模块/设备节点)。"
else
    echo "[OK] nvidia-smi: $(nvidia-smi -L | head -1)"
fi

# ---------- A) OpenCL/CUDA 是否真发现 GPU (轻量, em.tpr 5 步) ------------------
if [ ! -f "$SMOKE_DIR/em.tpr" ]; then
    fail "缺 $SMOKE_DIR/em.tpr — 先在 CPU 端备好 smoke 体系 (见 docs/A40_RUNBOOK.md)。"
else
    echo "---- Test A: nonbonded offload (-nb gpu -pme gpu) ----"
    if CUDA_VISIBLE_DEVICES=$GPU_ID "$GMX" mdrun -nb gpu -pme gpu -ntmpi 1 -gpu_id "$GPU_ID" \
        -s "$SMOKE_DIR/em.tpr" -nsteps 5 -deffnm "$WORK/a" > "$WORK/a.log" 2>&1; then
        echo "[OK] Test A 通过"
    else
        fail "Test A 失败 (GPU 没被发现?)。看尾部日志:"; grep -Ei "0 detected|no.*opencl|no.*cuda|fatal|error" "$WORK/a.log" | tail -5
    fi
fi

# ---------- B) ⭐ 用生产 flags + md 积分器 tpr 实跑 (真正的放行闸) -------------
EM_GRO="$SMOKE_DIR/em.gro"
TOPOL="$(find "$SMOKE_DIR" -maxdepth 3 -name topol.top 2>/dev/null | head -1)"
if [ ! -f "$EM_GRO" ] || [ -z "$TOPOL" ]; then
    fail "缺 em.gro 或 topol.top (在 $SMOKE_DIR 下), 无法做 B 测试。"
else
    cat > "$WORK/nvt_smoke.mdp" <<'MDP'
integrator           = md
nsteps               = 5
dt                   = 0.002
cutoff-scheme        = Verlet
coulombtype          = PME
rcoulomb             = 1.2
rvdw                 = 1.2
tcoupl               = V-rescale
tc-grps              = System
tau_t                = 0.1
ref_t                = 310
constraints          = h-bonds
constraint_algorithm = lincs
gen_vel              = yes
gen_temp             = 310
MDP
    echo "---- Test B: 生产 flags '-nb gpu -pme gpu ${GPU_RESIDENT_FLAGS}' (md 积分器) ----"
    if "$GMX" grompp -f "$WORK/nvt_smoke.mdp" -c "$EM_GRO" -p "$TOPOL" \
        -o "$WORK/nvt_smoke.tpr" -maxwarn 5 > "$WORK/b_grompp.log" 2>&1; then
        if CUDA_VISIBLE_DEVICES=$GPU_ID "$GMX" mdrun -deffnm "$WORK/b" \
            -nb gpu -pme gpu ${GPU_RESIDENT_FLAGS} -pin on -ntmpi 1 -gpu_id "$GPU_ID" \
            -s "$WORK/nvt_smoke.tpr" -nsteps 5 > "$WORK/b.log" 2>&1; then
            echo "[OK] Test B 通过 — 生产 runner 的 mdrun flags 在本机可跑"
        else
            fail "Test B 失败。这正是会在真跑 NVT 时挂掉的点。关键日志:"
            grep -Ei "not supported|update task|bonded|fatal|error" "$WORK/b.log" | tail -8
            echo "    → 若是 'not supported with OpenCL':runner 已会自动降级,但请确认 GPU_BACKEND 检测正确;"
            echo "      或改用 CUDA 构建的 gromacs (见 docs/A40_RUNBOOK.md)。"
        fi
    else
        fail "Test B grompp 失败 (多半是 GMXLIB/force_fields 没指对)。日志:"; tail -8 "$WORK/b_grompp.log"
    fi
fi

echo "============================================================"
if [ "$FAIL" -eq 0 ]; then
    echo "✅ SMOKE PASS — 放行。可以上批量:"
    echo "   bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --gpus <N>"
    exit 0
else
    echo "❌ SMOKE FAIL — 不要上批量。按上面提示修。"
    exit 1
fi

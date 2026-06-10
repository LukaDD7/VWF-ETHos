#!/bin/bash
# ==============================================================================
# run_gromacs_vwf_md.sh — VWF MD Simulation (GROMACS, H200 GPU)
# ==============================================================================
#
# 支持三种分子系统：
#   complex   — A1-GPIbα 复合体 (默认)
#   monomer  — VWF A1 单体（从复合体 CIF 提取 chain A）
#   autoinhib — VWF D'D3-A1 自抑制模块（需 A1+D'D3 Boltz-2 结果）
#
# 工作流程（三种系统共用）：
#   1. CIF → PDB 转换 (gemmi)
#   2. gmx pdb2gmx: 拓扑构建 (CHARMM36m + TIP3P)
#   3. gmx solvate + genion: 溶剂化 + 加离子 (0.15M NaCl)
#   4. EM → NVT → NPT → Production MD
#   5. 后分析: RMSF, PCA, MM/PBSA (complex), D'D3-A1 contacts (autoinhib)
#
# 环境依赖：
#   conda create -n gromacs python=3.11
#   conda install -c conda-forge gromacs=2025
#   pip install gemmi
#   conda install -c conda-forge gmx_mmpbsa  # 可选
#
# 用法：
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --system complex --gpus 8   # 默认
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --system monomer --gpus 8
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --gpus 8
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --preflight
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --phase em --gpus 4
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --ns 500 --gpus 8
#
# 输入来源（按 --system）：
#   complex   → output/boltz2_a1_gpiba_results/
#   monomer   → output/boltz2_a1_gpiba_results/ (提取 chain A)
#   autoinhib → output/boltz2_a1_dp_d3_results/ (A1+D'D3)
# ==============================================================================

set -u

# ---------- 默认参数 -----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
SYSTEM_TYPE="complex"          # complex | monomer | autoinhib
N_GPUS=4
PROD_NS=200
PHASE="prod"
PREFLIGHT_ONLY=false
VARIANT_FILTER=""

# ---------- 参数解析 -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --gpus)             N_GPUS="$2"; shift 2 ;;
        --system)           SYSTEM_TYPE="$2"; shift 2 ;;
        --boltz-results)    BOLTZ_RESULTS="$2"; shift 2 ;;
        --monomer-results)  MONOMER_RESULTS="$2"; shift 2 ;;
        --autoinhib-results) AUTINHIB_RESULTS="$2"; shift 2 ;;
        --out-dir)          OUTPUT_DIR="$2"; shift 2 ;;
        --ns)               PROD_NS="$2"; shift 2 ;;
        --phase)            PHASE="$2"; shift 2 ;;
        --preflight)        PREFLIGHT_ONLY=true; shift ;;
        --filter)           VARIANT_FILTER="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,70p' "$0" | sed 's/^# //'
            exit 0 ;;
        *) echo "[WARN] Unknown argument: $1"; shift ;;
    esac
done

# ---------- log() 提前定义 (gmx 路径检测段就要用) ----------------------------
# LOG 文件路径临时用 /dev/null,等 OUTPUT_DIR 确定后再换;
# 之前的实现把 log() 定义在 line 153,在 gmx 检测(line 95-112)之后,
# 导致检测段调 log 时 bash 报 "log: 未找到命令" 错。
LOG="/dev/null"
log() {
    local msg="[$(date '+%H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOG"
}

# CHARMM36m 在 gromacs env 内,需要明确指定 GMXDATA
# 优先 LZY NFS (/lzy/envs/gromacs/...) 共享路径,失败再回退到本地重建 env
if [ -d "$ROOT_DIR/../../envs/gromacs/share/gromacs" ]; then
    export GMXDATA="$(cd "$ROOT_DIR/../../envs/gromacs/share/gromacs" && pwd)"
elif [ -d "$ROOT_DIR/envs/gromacs/share/gromacs" ]; then
    export GMXDATA="$(cd "$ROOT_DIR/envs/gromacs/share/gromacs" && pwd)"
else
    export GMXDATA=""
fi

# Force field: use project-patched charmm36m.ff (adds per-residue <RESNAME>1
# N-terminal patches that the charmm2gmx port was missing — see
# force_fields/charmm36m.ff/aminoacids.n.tdb for details).
# Set GMXLIB so gmx pdb2gmx loads the patched FF first, before the conda
# env's pristine copy. We do this by symlinking the patched FF into a
# runtime directory and pointing GMXLIB there; the symlink + cwd trick also
# makes gmx find it ahead of the system path.
export GMXLIB="${GMXLIB_OVERRIDE:-$ROOT_DIR/force_fields}"

# ---------- OpenCL ICD vendor 注册 (GPU 实例无 /etc/OpenCL/vendors 时用 NFS 替代) ---
# GPU 实例上 NVIDIA OpenCL ICD lib 存在 (/usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.*)
# 但 /etc/OpenCL/vendors/ 可能不存在, 导致 OpenCL loader 找不到 NVIDIA 设备.
# 我们把 nvidia.icd 放到项目内 opencl_vendors/ 并设 OCL_ICD_VENDORS 指向那里.
OPENCL_VENDORS_DIR="$ROOT_DIR/opencl_vendors"
if [ -d "$OPENCL_VENDORS_DIR" ] && [ -n "$(ls "$OPENCL_VENDORS_DIR"/*.icd 2>/dev/null)" ]; then
    export OCL_ICD_VENDORS="$OPENCL_VENDORS_DIR"
    log "[OK] OCL_ICD_VENDORS=$OCL_ICD_VENDORS (NFS-based OpenCL ICD)"
else
    log "[WARN] No OpenCL ICD files in $OPENCL_VENDORS_DIR"
    log "  GROMACS may fail to detect GPU. nvidia-smi working is required."
fi

# ---------- 定位 gmx (无需 conda activate) --------------------------------------
# GPU 实例通过 NFS 共享 /lzy/envs/gromacs, 但 conda 默认不搜这个路径.
# 我们用绝对路径, 不依赖 conda activate.
if command -v gmx &>/dev/null; then
    GMX="$(command -v gmx)"
    log "[OK] gmx: $GMX (from PATH)"
elif [ -x "$ROOT_DIR/../../envs/gromacs/bin.AVX2_256/gmx" ]; then
    GMX="$ROOT_DIR/../../envs/gromacs/bin.AVX2_256/gmx"
    log "[OK] gmx: $GMX (absolute path, AVX2_256 SIMD, via LZY NFS /lzy/)"
elif [ -x "$ROOT_DIR/../../envs/gromacs/bin/gmx" ]; then
    GMX="$ROOT_DIR/../../envs/gromacs/bin/gmx"
    log "[OK] gmx: $GMX (absolute path, default SIMD, via LZY NFS /lzy/)"
elif [ -x "$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx" ]; then
    # Local rebuilt env (A40 独立机,无 /lzy NFS 时的 A40_AGENT_SETUP 路径)
    GMX="$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx"
    log "[OK] gmx: $GMX (absolute path, AVX2_256 SIMD, 本地重建 env)"
elif [ -x "$ROOT_DIR/envs/gromacs/bin/gmx" ]; then
    GMX="$ROOT_DIR/envs/gromacs/bin/gmx"
    log "[OK] gmx: $GMX (absolute path, default SIMD, 本地重建 env)"
else
    log "[ERROR] gmx not found in PATH, $ROOT_DIR/../../envs/gromacs/ or $ROOT_DIR/envs/gromacs/"
    log "  解决: 跑 scripts/setup/a40_selfcheck.sh,按 A40_AGENT_SETUP.md 在本地建 env;或挂 /lzy NFS。"
    PREFLIGHT_FAIL=true
    GMX=""
fi
export GMX

# ---------- 根据系统类型设置默认路径 ------------------------------------------
case "$SYSTEM_TYPE" in
    complex)
        BOLTZ_RESULTS="${BOLTZ_RESULTS:-${ROOT_DIR}/output/boltz2_a1_gpiba_results}"
        OUTPUT_DIR="${OUTPUT_DIR:-${ROOT_DIR}/output/gromacs_md}"
        ;;
    monomer)
        BOLTZ_RESULTS="${MONOMER_RESULTS:-${ROOT_DIR}/output/boltz2_a1_gpiba_results}"
        OUTPUT_DIR="${OUTPUT_DIR:-${ROOT_DIR}/output/gromacs_md_monomer}"
        ;;
    autoinhib)
        BOLTZ_RESULTS="${AUTINHIB_RESULTS:-${ROOT_DIR}/output/boltz2_a1_dp_d3_results}"
        OUTPUT_DIR="${OUTPUT_DIR:-${ROOT_DIR}/output/gromacs_md_autoinhib}"
        ;;
    *)
        echo "[ERROR] Unknown system type: $SYSTEM_TYPE (must be complex|monomer|autoinhib)"
        exit 1
        ;;
esac
MDP_DIR="${ROOT_DIR}/scripts/pipeline/mdp"

mkdir -p "$OUTPUT_DIR"
# 之前的"临时 log"写入 /dev/null,这里重定向到真正的日志
LOG="${OUTPUT_DIR}/run_log.txt"
# (log() 函数本身已在文件顶部定义,这里无需重复定义)

# ---------- 预检 (Preflight) --------------------------------------------------
log "============================================================"
log "VWF GROMACS MD Pipeline"
log "System: $SYSTEM_TYPE"
log "Started: $(date)"
log "============================================================"

PREFLIGHT_FAIL=false

# 1. GROMACS
if [ -z "${GMX:-}" ]; then
    log "[ERROR] 'gmx' not found."
    log "  解决: conda activate gromacs, 或确保 NFS 挂载 /lzy/envs/gromacs/"
    PREFLIGHT_FAIL=true
else
    GMX_VER=$("$GMX" --version 2>&1 | head -1)
    log "[OK] GROMACS: $GMX_VER"
    # 检测 GPU backend
    GPU_BACKEND=$("$GMX" mdrun -version 2>&1 | grep "GPU support:" | awk '{print $NF}')
    OPENCL_ICD=$(ls /usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.* 2>/dev/null | head -1)
    log "[INFO] GPU backend: $GPU_BACKEND"
    log "[INFO] NVIDIA OpenCL ICD: ${OPENCL_ICD:-NOT FOUND (need nvidia driver)}"

    # ---- 按后端条件化 GPU offload flags ----------------------------------
    # CRITICAL: GROMACS 的 -bonded gpu / -update gpu (GPU 常驻) 和 CUDA graphs
    # 只在 CUDA/SYCL 构建支持; OpenCL 构建会 fatal ("not supported with OpenCL").
    # 因此按实测后端选 flags, 否则 NVT 一开跑就挂、白烧 GPU 机时.
    #   CUDA/SYCL → 全 GPU 常驻 (快, H200/A40 利用率最大化)
    #   OpenCL/未知 → 只 -nb/-pme gpu, update/bonded 回 CPU (慢但能跑)
    case "$GPU_BACKEND" in
        CUDA|SYCL)
            GPU_RESIDENT_FLAGS="-bonded gpu -update gpu"
            export GMX_CUDA_GRAPH=1
            log "[OK] GPU-resident mode ON (backend=$GPU_BACKEND): flags='-nb gpu -pme gpu $GPU_RESIDENT_FLAGS', CUDA graph=1"
            ;;
        *)
            GPU_RESIDENT_FLAGS=""
            unset GMX_CUDA_GRAPH
            log "[WARN] backend='$GPU_BACKEND' 不支持 GPU 常驻 (-bonded/-update gpu)。"
            log "       降级为 '-nb gpu -pme gpu' (update/bonded 在 CPU)。如需提速请改用 CUDA 构建的 gromacs。"
            ;;
    esac
    export GPU_RESIDENT_FLAGS
fi

# 2. gemmi (CIF→PDB) — gromacs env 自带, 用 env 内 python 检测
# 同时支持两种 env 布局:
#   (a) LZY NFS 共享: /lzy/envs/gromacs/{bin.AVX2_256/gmx, bin/python}
#   (b) A40 本地重建: $ROOT_DIR/envs/gromacs/{bin.AVX2_256/gmx, bin/python}
GMX_PY=""
for cand in \
    "$ROOT_DIR/../../envs/gromacs/bin/python" \
    "$ROOT_DIR/envs/gromacs/bin/python"; do
    if [ -x "$cand" ]; then GMX_PY="$cand"; break; fi
done
# Fallback: 用绝对路径重拼 (对 gmx 在 bin.AVX2_256/ 子目录的 layout,
# $GMX/../../bin/python 解析到 envs/gromacs/bin/python)
if [ -z "$GMX_PY" ] && [ -n "${GMX:-}" ]; then
    cand_py="$(cd "$(dirname "$GMX")/../.." 2>/dev/null && pwd)/bin/python"
    [ -x "$cand_py" ] && GMX_PY="$cand_py"
fi
if "$GMX_PY" -c "import gemmi" 2>/dev/null; then
    GEMMI_VER=$("$GMX_PY" -c "import gemmi; print(gemmi.__version__)" 2>/dev/null)
    log "[OK] gemmi: $GEMMI_VER (in gromacs env)"
else
    log "[WARN] gemmi not found in gromacs env (needed for CIF→PDB conversion)."
    log "  解决: /lzy/envs/gromacs/bin/python -m pip install gemmi"
fi
export GMX_PY

# 3. GPU 检测 — 用 nvidia-smi, 不依赖 torch (gromacs env 无 torch)
GPU_COUNT=$(nvidia-smi -L 2>/dev/null | wc -l || echo "0")
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
log "[INFO] GPUs detected: $GPU_COUNT ($GPU_NAME)"
if [ "$GPU_COUNT" -lt "$N_GPUS" ] 2>/dev/null; then
    log "[WARN] Adjusting N_GPUS from $N_GPUS to $GPU_COUNT"
    N_GPUS=$GPU_COUNT
fi
if [ "$GPU_COUNT" -eq 0 ]; then
    log "[ERROR] No GPUs found."
    PREFLIGHT_FAIL=true
fi

# 4. Boltz-2 输入
if [ ! -d "$BOLTZ_RESULTS" ]; then
    log "[ERROR] Boltz-2 results not found: $BOLTZ_RESULTS"
    case "$SYSTEM_TYPE" in
        complex)     log "  先运行: bash scripts/pipeline/run_a1_gpiba_boltz2.sh" ;;
        monomer)    log "  先运行: bash scripts/pipeline/run_a1_gpiba_boltz2.sh (monomer uses same CIFs)" ;;
        autoinhib)   log "  先运行: python3 scripts/pipeline/generate_a1_dp_d3_yamls.py && bash scripts/pipeline/run_a1_dp_d3_boltz2.sh" ;;
    esac
    PREFLIGHT_FAIL=true
else
    if [ "$SYSTEM_TYPE" = "complex" ]; then
        CIF_COUNT=$(find "$BOLTZ_RESULTS" -name "*_model_0.cif" 2>/dev/null | wc -l)
    else
        CIF_COUNT=$(find "$BOLTZ_RESULTS" -name "*.cif" 2>/dev/null | wc -l)
    fi
    log "[OK] Boltz-2 CIF structures found: $CIF_COUNT"
fi

# 5. MDP 参数文件
if [ ! -d "$MDP_DIR" ]; then
    log "[INFO] MDP directory not found: $MDP_DIR — will generate default MDP files"
fi

# 6. 磁盘空间（每个 200ns MD ≈ 2-5 GB）
AVAIL_MB=$(df -m "$OUTPUT_DIR" 2>/dev/null | awk 'NR==2{print $4}' || echo "unknown")
log "[INFO] Disk available: ${AVAIL_MB}MB"

log "------------------------------------------------------------"
log "Config: SYSTEM=$SYSTEM_TYPE | N_GPUS=$N_GPUS | PROD_NS=${PROD_NS}ns | PHASE=$PHASE"
log "Filter: ${VARIANT_FILTER:-all variants}"
log "------------------------------------------------------------"

if $PREFLIGHT_FAIL; then
    log "[FATAL] Preflight failed."
    exit 1
fi

if $PREFLIGHT_ONLY; then
    log "Preflight complete. No simulations run."
    exit 0
fi

# ---------- 生成 MDP 文件（如果不存在）-----------------------------------------
generate_mdp_files() {
    local MDP_OUT="$1"
    mkdir -p "$MDP_OUT"

    # Energy Minimization
    cat > "$MDP_OUT/em.mdp" << 'EMEOF'
; Energy Minimization
; integrator: cg (conjugate gradient) — 比 steep 快 ~5-10x, 支持 constraints
; (l-bfgs 不支持 constraints — 见 gromacs 2025.3 错误 "L-BFGS + constraints not implemented")
integrator  = cg
emtol       = 1000.0    ; kJ/mol/nm (production-ready)
emstep      = 0.01      ; cg 初始步长; gromacs 自适应
nsteps      = 5000      ; cg 在 ~500-1500 步内收敛; 5000 是硬上限
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
EMEOF

    # NVT Equilibration (500 ps)
    cat > "$MDP_OUT/nvt.mdp" << 'NVTEOF'
; NVT equilibration
integrator  = md
nsteps      = 250000    ; 500 ps
dt          = 0.002
nstxout-compressed = 5000
nstlog      = 5000
nstenergy   = 5000
nstlist     = 200       ; H200 优化: 减少 neighbor list 更新频率
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
; Temperature coupling
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 310   310       ; 37°C (生理温度)
; No pressure coupling
pcoupl      = no
; Constraints
constraints     = h-bonds
constraint_algorithm = lincs
continuation    = no
gen_vel         = yes
gen_temp        = 310
gen_seed        = -1
NVTEOF

    # NPT Equilibration (500 ps)
    cat > "$MDP_OUT/npt.mdp" << 'NPTEOF'
; NPT equilibration
integrator  = md
nsteps      = 250000    ; 500 ps
dt          = 0.002
nstxout-compressed = 5000
nstlog      = 5000
nstenergy   = 5000
nstlist     = 200
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
; Temperature
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 310   310
; Pressure
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
; Constraints
constraints     = h-bonds
constraint_algorithm = lincs
continuation    = yes
gen_vel         = no
NPTEOF

    # Production MD
    local PROD_STEPS=$((PROD_NS * 500000))  # dt=0.002, 500000 steps = 1 ns
    cat > "$MDP_OUT/production.mdp" << PRODEOF
; Production MD (${PROD_NS} ns)
integrator  = md
nsteps      = ${PROD_STEPS}
dt          = 0.002
nstxout-compressed = 50000   ; 每 100ps 写一帧 (2000 帧 / 200ns)
nstlog      = 50000
nstenergy   = 50000
nstlist     = 200            ; H200 优化
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
; Temperature
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1   0.1
ref_t       = 310   310
; Pressure
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
; Constraints
constraints     = h-bonds
constraint_algorithm = lincs
continuation    = yes
gen_vel         = no
PRODEOF

    log "[OK] MDP files generated in $MDP_OUT"
}

# 生成或使用已有 MDP
if [ ! -d "$MDP_DIR" ] || [ ! -f "$MDP_DIR/em.mdp" ]; then
    MDP_DIR="${OUTPUT_DIR}/mdp"
    generate_mdp_files "$MDP_DIR"
fi

# ---------- 扫描待跑突变体列表 ------------------------------------------------
TODO_VARIANTS=()

if [ "$SYSTEM_TYPE" = "autoinhib" ]; then
    # Autoinhib: output/boltz2_a1_dp_d3_results/boltz_results_VWF_WT_dp_d3_a1/
    for BOLTZ_DIR in "$BOLTZ_RESULTS"/boltz_results_VWF_*; do
        [ ! -d "$BOLTZ_DIR" ] && continue
        JNAME=$(basename "$BOLTZ_DIR" | sed 's/^boltz_results_//')
        VARIANT=$(echo "$JNAME" | sed 's/_dp_d3_a1$//')
        CIF=$(find "$BOLTZ_DIR" -name "*.cif" 2>/dev/null | head -1)
        [ -z "$CIF" ] && continue
        DONE_FILE="${OUTPUT_DIR}/${VARIANT}/.done_${PHASE}"
        [ -f "$DONE_FILE" ] && continue
        TODO_VARIANTS+=("$VARIANT|$CIF|complex")
    done
elif [ "$SYSTEM_TYPE" = "monomer" ]; then
    # Monomer: use A1-GPIbα results but extract chain A only
    for BOLTZ_DIR in "$BOLTZ_RESULTS"/boltz_results_VWF_*; do
        [ ! -d "$BOLTZ_DIR" ] && continue
        JNAME=$(basename "$BOLTZ_DIR" | sed 's/^boltz_results_//')
        VARIANT=$(echo "$JNAME" | sed 's/_vs_GPIb_alpha$//')
        CIF_PATH="${BOLTZ_DIR}/predictions/${JNAME}/${JNAME}_model_0.cif"
        [ ! -f "$CIF_PATH" ] && continue
        DONE_FILE="${OUTPUT_DIR}/${VARIANT}/.done_${PHASE}"
        [ -f "$DONE_FILE" ] && continue
        TODO_VARIANTS+=("$VARIANT|$CIF_PATH|monomer")
    done
else
    # Complex: A1-GPIbα
    for BOLTZ_DIR in "$BOLTZ_RESULTS"/boltz_results_VWF_*; do
        [ ! -d "$BOLTZ_DIR" ] && continue
        JNAME=$(basename "$BOLTZ_DIR" | sed 's/^boltz_results_//')
        VARIANT=$(echo "$JNAME" | sed 's/_vs_GPIb_alpha$//')
        CIF="${BOLTZ_DIR}/predictions/${JNAME}/${JNAME}_model_0.cif"
        [ ! -f "$CIF" ] && continue
        if [ -n "$VARIANT_FILTER" ]; then
            [[ "$VARIANT" != *"$VARIANT_FILTER"* ]] && continue
        fi
        DONE_FILE="${OUTPUT_DIR}/${VARIANT}/.done_${PHASE}"
        [ -f "$DONE_FILE" ] && continue
        TODO_VARIANTS+=("$VARIANT|$CIF|complex")
    done
fi

log "MD variants: found=${#TODO_VARIANTS[@]}  phase=$PHASE"

if [ ${#TODO_VARIANTS[@]} -eq 0 ]; then
    log "All matching variants already completed for phase=$PHASE."
    exit 0
fi

# ---------- Worker 函数 -------------------------------------------------------
run_md_worker() {
    local GPU_ID=$1
    local JOB_LIST_FILE=$2
    local OUT_BASE=$3
    local MDP_BASE=$4
    local TARGET_PHASE=$5
    local SYSTEM=$6
    local WORKER_LOG="${OUT_BASE}/worker_${GPU_ID}.log"

    echo "[GPU ${GPU_ID}] MD Worker started at $(date)" > "$WORKER_LOG"

    local JOB_OK=0
    local JOB_FAIL=0

    while IFS='|' read -r VARIANT CIF_PATH SYS; do
        [ -z "$VARIANT" ] && continue
        SYS="${SYS:-complex}"

        local WORK="${OUT_BASE}/${VARIANT}"
        mkdir -p "$WORK"/{input,topology,em,nvt,npt,production,analysis}

        echo "[GPU ${GPU_ID}] $(date '+%H:%M:%S') START: $VARIANT (sys=$SYS)" >> "$WORKER_LOG"

        # ---- Step 1: CIF → PDB (monomer = chain A only) ----
        if [ "$SYS" = "monomer" ]; then
            INPUT_PDB="$WORK/input/monomer.pdb"
        else
            INPUT_PDB="$WORK/input/complex.pdb"
        fi

        if [ ! -f "$INPUT_PDB" ]; then
            if [ "$SYS" = "monomer" ]; then
                "$GMX_PY" -c "
import gemmi
doc = gemmi.cif.read('$CIF_PATH')
st = gemmi.make_structure_from_block(doc[0])
for m in st:
    for c in list(m):
        if c.name != 'A':
            m.remove_chain(c.name)
st.write_pdb('$INPUT_PDB')
print('OK: CIF → PDB (chain A)')
" >> "$WORKER_LOG" 2>&1
            else
                "$GMX_PY" -c "
import gemmi
doc = gemmi.cif.read('$CIF_PATH')
st = gemmi.make_structure_from_block(doc[0])
st.write_pdb('$INPUT_PDB')
print('OK: CIF → PDB')
" >> "$WORKER_LOG" 2>&1
            fi
            if [ $? -ne 0 ]; then
                echo "[GPU ${GPU_ID}] FAIL CIF→PDB: $VARIANT" >> "$WORKER_LOG"
                JOB_FAIL=$((JOB_FAIL + 1)); continue
            fi
        fi

        # ---- Step 2: pdb2gmx (with project-patched charmm36m.ff) ----
        # Why the FF swap: the charmm2gmx port bundled in gromacs-2025.4
        # omits per-residue N-terminal patches (MET1, ALA1, ...). When
        # pdb2gmx encounters a protein N-terminal MET, it falls back to
        # the [ MET1 ] patch from ethers.n.tdb (which is for ether
        # chemistry, not proteins) and fails with
        # "atom C1 not found in buiding block 1MET". We copy charmm36m.ff
        # to force_fields/ and add protein-compatible <RESNAME>1 patches
        # to aminoacids.n.tdb. Setting GMXLIB=force_fields and running
        # pdb2gmx from a symlinked dir makes gmx find our patched copy
        # first.
        if [ ! -f "$WORK/topology/processed.gro" ]; then
            # Symlink project's patched charmm36m.ff into the topology
            # dir so pdb2gmx "current directory" lookup beats the env's
            # pristine copy.
            ln -sf "$ROOT_DIR/force_fields/charmm36m.ff" \
                "$WORK/topology/charmm36m.ff"
            cd "$WORK/topology"
            echo "1" | GMXLIB="$ROOT_DIR/force_fields" ${GMX} pdb2gmx \
                -f "$INPUT_PDB" \
                -o processed.gro \
                -p topol.top \
                -water tip3p \
                -ff charmm36m \
                -ignh \
                >> "$WORKER_LOG" 2>&1
            if [ $? -ne 0 ]; then
                echo "[GPU ${GPU_ID}] FAIL pdb2gmx: $VARIANT" >> "$WORKER_LOG"
                JOB_FAIL=$((JOB_FAIL + 1))
                continue
            fi
            cd "$OUT_BASE"
        fi

        # ---- Step 3: Box + Solvate + Ions ----
        if [ ! -f "$WORK/topology/solv_ions.gro" ]; then
            cd "$WORK/topology"
            ${GMX} editconf -f processed.gro -o newbox.gro -c -d 1.2 -bt dodecahedron >> "$WORKER_LOG" 2>&1
            ${GMX} solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top >> "$WORKER_LOG" 2>&1
            ${GMX} grompp -f "$MDP_BASE/em.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            echo "SOL" | ${GMX} genion -s ions.tpr -o solv_ions.gro -p topol.top \
                -pname NA -nname CL -neutral -conc 0.15 >> "$WORKER_LOG" 2>&1
            cd "$OUT_BASE"
        fi

        # ---- Step 4: Energy Minimization ----
        if [ ! -f "$WORK/.done_em" ]; then
            cd "$WORK/em"
            ${GMX} grompp -f "$MDP_BASE/em.mdp" -c "$WORK/topology/solv_ions.gro" \
                -p "$WORK/topology/topol.top" -o em.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            # GROMACS 2025.x: -pme gpu 不支持非动力学积分器 (steep/l-bfgs),
            # EM 段必须让 PME 回 CPU。只 offload NB 到 GPU (-nb gpu) 即可,
            # EM 跑得短 (~1-5 min),PME on CPU 几乎不影响总时长。
            CUDA_VISIBLE_DEVICES=$GPU_ID ${GMX} mdrun -v -deffnm em \
                -nb gpu -pme cpu >> "$WORKER_LOG" 2>&1
            if [ $? -eq 0 ]; then
                date '+%Y-%m-%d %H:%M:%S' > "$WORK/.done_em"
                echo "[GPU ${GPU_ID}] OK EM: $VARIANT" >> "$WORKER_LOG"
            else
                echo "[GPU ${GPU_ID}] FAIL EM: $VARIANT" >> "$WORKER_LOG"
                JOB_FAIL=$((JOB_FAIL + 1)); continue
            fi
            cd "$OUT_BASE"
        fi
        [ "$TARGET_PHASE" = "em" ] && { JOB_OK=$((JOB_OK + 1)); continue; }

        # ---- Step 5: NVT Equilibration ----
        if [ ! -f "$WORK/.done_equil" ]; then
            cd "$WORK/nvt"
            ${GMX} grompp -f "$MDP_BASE/nvt.mdp" -c "$WORK/em/em.gro" \
                -r "$WORK/em/em.gro" -p "$WORK/topology/topol.top" \
                -o nvt.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            CUDA_VISIBLE_DEVICES=$GPU_ID ${GMX} mdrun -deffnm nvt \
                -nb gpu -pme gpu ${GPU_RESIDENT_FLAGS} -pin on \
                >> "$WORKER_LOG" 2>&1
            [ $? -ne 0 ] && { echo "[GPU ${GPU_ID}] FAIL NVT: $VARIANT" >> "$WORKER_LOG"; JOB_FAIL=$((JOB_FAIL + 1)); continue; }

            # ---- Step 6: NPT Equilibration ----
            cd "$WORK/npt"
            ${GMX} grompp -f "$MDP_BASE/npt.mdp" -c "$WORK/nvt/nvt.gro" \
                -r "$WORK/nvt/nvt.gro" -t "$WORK/nvt/nvt.cpt" \
                -p "$WORK/topology/topol.top" -o npt.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            CUDA_VISIBLE_DEVICES=$GPU_ID ${GMX} mdrun -deffnm npt \
                -nb gpu -pme gpu ${GPU_RESIDENT_FLAGS} -pin on \
                >> "$WORKER_LOG" 2>&1
            if [ $? -eq 0 ]; then
                date '+%Y-%m-%d %H:%M:%S' > "$WORK/.done_equil"
                echo "[GPU ${GPU_ID}] OK EQUIL: $VARIANT" >> "$WORKER_LOG"
            else
                echo "[GPU ${GPU_ID}] FAIL NPT: $VARIANT" >> "$WORKER_LOG"
                JOB_FAIL=$((JOB_FAIL + 1)); continue
            fi
            cd "$OUT_BASE"
        fi
        [ "$TARGET_PHASE" = "equil" ] && { JOB_OK=$((JOB_OK + 1)); continue; }

        # ---- Step 7: Production MD ----
        if [ ! -f "$WORK/.done_prod" ]; then
            cd "$WORK/production"
            ${GMX} grompp -f "$MDP_BASE/production.mdp" -c "$WORK/npt/npt.gro" \
                -t "$WORK/npt/npt.cpt" -p "$WORK/topology/topol.top" \
                -o md_prod.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            # GPU offload: flags 按后端条件化 (见 preflight 的 GPU_RESIDENT_FLAGS)。
            # CUDA/SYCL → 全常驻 + CUDA graph(已 export GMX_CUDA_GRAPH=1); OpenCL → 仅 -nb/-pme gpu。
            CUDA_VISIBLE_DEVICES=$GPU_ID ${GMX} mdrun \
                -deffnm md_prod \
                -nb gpu -pme gpu ${GPU_RESIDENT_FLAGS} \
                -pin on -nstlist 200 \
                >> "$WORKER_LOG" 2>&1
            if [ $? -eq 0 ]; then
                date '+%Y-%m-%d %H:%M:%S' > "$WORK/.done_prod"
                echo "[GPU ${GPU_ID}] OK PROD: $VARIANT" >> "$WORKER_LOG"
            else
                echo "[GPU ${GPU_ID}] FAIL PROD: $VARIANT" >> "$WORKER_LOG"
                JOB_FAIL=$((JOB_FAIL + 1)); continue
            fi
            cd "$OUT_BASE"
        fi

        JOB_OK=$((JOB_OK + 1))
    done < "$JOB_LIST_FILE"

    echo "[GPU ${GPU_ID}] Worker done at $(date) | OK=$JOB_OK FAIL=$JOB_FAIL" >> "$WORKER_LOG"
}

export -f run_md_worker

# ---------- 分配到各 GPU（轮询）-----------------------------------------------
TMPDIR_LISTS="${OUTPUT_DIR}/_worker_lists"
mkdir -p "$TMPDIR_LISTS"

for i in $(seq 0 $((N_GPUS - 1))); do
    > "${TMPDIR_LISTS}/gpu_${i}.txt"
done

IDX=0
for ENTRY in "${TODO_VARIANTS[@]}"; do
    GPU_IDX=$((IDX % N_GPUS))
    echo "$ENTRY" >> "${TMPDIR_LISTS}/gpu_${GPU_IDX}.txt"
    IDX=$((IDX + 1))
done

log "Distribution across $N_GPUS GPUs:"
for i in $(seq 0 $((N_GPUS - 1))); do
    COUNT=$(wc -l < "${TMPDIR_LISTS}/gpu_${i}.txt" 2>/dev/null || echo 0)
    log "  GPU $i: $COUNT variants"
done

# ---------- 启动并行 Workers ---------------------------------------------------
log "Launching $N_GPUS parallel MD workers..."
PIDS=()

for i in $(seq 0 $((N_GPUS - 1))); do
    LIST_FILE="${TMPDIR_LISTS}/gpu_${i}.txt"
    JOB_COUNT=$(wc -l < "$LIST_FILE" 2>/dev/null || echo 0)
    [ "$JOB_COUNT" -eq 0 ] && continue

    run_md_worker "$i" "$LIST_FILE" "$OUTPUT_DIR" "$MDP_DIR" "$PHASE" "$SYSTEM_TYPE" &
    PIDS+=($!)
    log "  Worker GPU-$i started (PID=${PIDS[-1]}, variants=$JOB_COUNT)"
done

log ""
log "Workers running. Monitor with:"
log "  watch -n 60 'find $OUTPUT_DIR -name .done_prod | wc -l'"
log "  tail -f $OUTPUT_DIR/worker_0.log"

# ---------- 等待完成 -----------------------------------------------------------
for PID in "${PIDS[@]}"; do
    wait "$PID"
    RC=$?
    [ $RC -ne 0 ] && log "[WARN] Worker PID=$PID exited with code $RC"
done

FINAL_DONE=$(find "$OUTPUT_DIR" -name ".done_${PHASE}" 2>/dev/null | wc -l)
log ""
log "============================================================"
log "MD Run finished: $(date)"
log "  Phase completed (.done_${PHASE}): $FINAL_DONE"
log "  Output dir: $OUTPUT_DIR"
log "============================================================"

rm -rf "$TMPDIR_LISTS"

#!/bin/bash
# ==============================================================================
# run_gromacs_vwf_md.sh — VWF A1–GPIbα MD Simulation (GROMACS, H200 GPU)
# ==============================================================================
#
# 从 Boltz-2 预测的 CIF 结构出发，使用 GROMACS 运行分子动力学模拟。
# 目标：提取 2B (GOF) vs 2M (LOF) 的动态分类特征。
#
# 工作流程：
#   1. CIF → PDB 转换 (gemmi/obabel)
#   2. gmx pdb2gmx: 拓扑构建 (CHARMM36m + TIP3P)
#   3. gmx solvate + genion: 溶剂化 + 加离子 (0.15M NaCl)
#   4. Energy Minimization (EM)
#   5. NVT Equilibration (500 ps)
#   6. NPT Equilibration (500 ps)
#   7. Production MD (200 ns)
#   8. 后分析: RMSF, PCA, MM/PBSA (可选)
#
# GPU 实例约束（同 Boltz-2 pipeline）：
#   - 不可联网 → 力场文件必须预先包含在 GROMACS 安装中
#   - 计费制 → 最大化 GPU 利用率，多突变体并行
#   - /dev/shm 可能有限 → GROMACS 不依赖 /dev/shm，无影响
#
# H200 141GB 性能估算（VWF A1 + GPIbα ≈ 65K atoms with solvent）：
#   - 全 GPU offload 模式：~200-300 ns/day
#   - 200 ns production run ≈ 16-24 小时
#   - 8 GPU 并行 8 个突变体 ≈ 同时完成
#
# 环境依赖（CPU 实例预装后传输）：
#   conda create -n gromacs python=3.11
#   conda install -c conda-forge gromacs=2025
#   pip install gemmi  # CIF→PDB 转换
#   pip install gmx_MMPBSA  # MM/PBSA (可选)
#   力场文件: CHARMM36m（GROMACS 2025 内置）
#
# 用法：
#   bash scripts/pipeline/run_gromacs_vwf_md.sh                  # 默认 4 GPU
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --gpus 8         # 8 GPU 并行
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --preflight      # 仅预检
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --phase em       # 只跑到 EM
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --phase equil    # 只跑到平衡
#   bash scripts/pipeline/run_gromacs_vwf_md.sh --ns 500         # 500 ns production
#
# 输入：
#   output/boltz2_a1_gpiba_results/boltz_results_VWF_*/predictions/model_0.cif
#   （由 run_a1_gpiba_boltz2.sh 产出的 Boltz-2 CIF 结构）
#
# 输出结构：
#   output/gromacs_md/
#     VWF_R1306W/
#       input/        ← 转换后的 PDB
#       topology/     ← pdb2gmx 输出
#       em/           ← Energy Minimization
#       nvt/          ← NVT 平衡
#       npt/          ← NPT 平衡
#       production/   ← 200 ns MD 轨迹
#       analysis/     ← RMSF, PCA, etc.
#       .done_em / .done_equil / .done_prod  ← 阶段完成标记
#     worker_0.log
#     run_log.txt
# ==============================================================================

set -u

# ---------- 默认参数 -----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BOLTZ_RESULTS="${ROOT_DIR}/output/boltz2_a1_gpiba_results"
OUTPUT_DIR="${ROOT_DIR}/output/gromacs_md"
MDP_DIR="${ROOT_DIR}/scripts/pipeline/mdp"
N_GPUS=4
PROD_NS=200          # Production MD 时长 (ns)
PHASE="prod"         # em / equil / prod
PREFLIGHT_ONLY=false
VARIANT_FILTER=""    # 空=全部, "R1306W"=仅指定突变

# ---------- 参数解析 -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --gpus)           N_GPUS="$2"; shift 2 ;;
        --boltz-results)  BOLTZ_RESULTS="$2"; shift 2 ;;
        --out-dir)        OUTPUT_DIR="$2"; shift 2 ;;
        --ns)             PROD_NS="$2"; shift 2 ;;
        --phase)          PHASE="$2"; shift 2 ;;
        --preflight)      PREFLIGHT_ONLY=true; shift ;;
        --filter)         VARIANT_FILTER="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,55p' "$0" | sed 's/^# //'
            exit 0 ;;
        *) echo "[WARN] Unknown argument: $1"; shift ;;
    esac
done

mkdir -p "$OUTPUT_DIR"
LOG="${OUTPUT_DIR}/run_log.txt"

log() {
    local msg="[$(date '+%H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOG"
}

# ---------- 预检 (Preflight) --------------------------------------------------
log "============================================================"
log "VWF A1-GPIbα GROMACS MD Pipeline"
log "Started: $(date)"
log "============================================================"

PREFLIGHT_FAIL=false

# 1. GROMACS
if ! command -v gmx &>/dev/null; then
    log "[ERROR] 'gmx' not found."
    log "  安装: conda install -c conda-forge gromacs=2025"
    PREFLIGHT_FAIL=true
else
    GMX_VER=$(gmx --version 2>&1 | head -1)
    log "[OK] GROMACS: $GMX_VER"
fi

# 2. gemmi (CIF→PDB)
if ! command -v gemmi &>/dev/null && ! python3 -c "import gemmi" &>/dev/null; then
    log "[WARN] gemmi not found (needed for CIF→PDB conversion)."
    log "  安装: pip install gemmi"
fi

# 3. CUDA / GPU
GPU_COUNT=$(python3 -c "import torch; print(torch.cuda.device_count())" 2>/dev/null || echo "0")
if [ "$GPU_COUNT" -eq 0 ]; then
    # GROMACS 用 nvidia-smi 检测
    GPU_COUNT=$(nvidia-smi -L 2>/dev/null | wc -l || echo "0")
fi
log "[INFO] GPUs detected: $GPU_COUNT"
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
    log "  先运行: bash scripts/pipeline/run_a1_gpiba_boltz2.sh"
    PREFLIGHT_FAIL=true
else
    CIF_COUNT=$(find "$BOLTZ_RESULTS" -name "model_0.cif" 2>/dev/null | wc -l)
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
log "Config: N_GPUS=$N_GPUS | PROD_NS=${PROD_NS}ns | PHASE=$PHASE"
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
integrator  = steep
emtol       = 1000.0    ; kJ/mol/nm
emstep      = 0.01
nsteps      = 50000
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

for BOLTZ_DIR in "$BOLTZ_RESULTS"/boltz_results_VWF_*; do
    [ ! -d "$BOLTZ_DIR" ] && continue
    JNAME=$(basename "$BOLTZ_DIR" | sed 's/^boltz_results_//')
    # 只取 A1-GPIba complex 的结果
    VARIANT=$(echo "$JNAME" | sed 's/_vs_GPIb_alpha$//')
    CIF="${BOLTZ_DIR}/predictions/model_0.cif"

    [ ! -f "$CIF" ] && continue

    # 过滤
    if [ -n "$VARIANT_FILTER" ]; then
        [[ "$VARIANT" != *"$VARIANT_FILTER"* ]] && continue
    fi

    # 检查阶段完成标记
    DONE_FILE="${OUTPUT_DIR}/${VARIANT}/.done_${PHASE}"
    if [ -f "$DONE_FILE" ]; then
        continue
    fi

    TODO_VARIANTS+=("$VARIANT|$CIF")
done

log "MD variants: total found=$(ls -d "$BOLTZ_RESULTS"/boltz_results_VWF_* 2>/dev/null | wc -l)  todo=${#TODO_VARIANTS[@]}"

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
    local WORKER_LOG="${OUT_BASE}/worker_${GPU_ID}.log"

    echo "[GPU ${GPU_ID}] MD Worker started at $(date)" > "$WORKER_LOG"

    local JOB_OK=0
    local JOB_FAIL=0

    while IFS='|' read -r VARIANT CIF_PATH; do
        [ -z "$VARIANT" ] && continue

        local WORK="${OUT_BASE}/${VARIANT}"
        mkdir -p "$WORK"/{input,topology,em,nvt,npt,production,analysis}

        echo "[GPU ${GPU_ID}] $(date '+%H:%M:%S') START: $VARIANT" >> "$WORKER_LOG"

        # ---- Step 1: CIF → PDB ----
        if [ ! -f "$WORK/input/complex.pdb" ]; then
            python3 -c "
import gemmi
doc = gemmi.cif.read('$CIF_PATH')
st = gemmi.make_structure_from_block(doc[0])
st.write_pdb('$WORK/input/complex.pdb')
print('OK: CIF → PDB')
" >> "$WORKER_LOG" 2>&1
            if [ $? -ne 0 ]; then
                echo "[GPU ${GPU_ID}] FAIL CIF→PDB: $VARIANT" >> "$WORKER_LOG"
                JOB_FAIL=$((JOB_FAIL + 1))
                continue
            fi
        fi

        # ---- Step 2: pdb2gmx ----
        if [ ! -f "$WORK/topology/processed.gro" ]; then
            cd "$WORK/topology"
            echo "1" | gmx pdb2gmx \
                -f "$WORK/input/complex.pdb" \
                -o processed.gro \
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
            gmx editconf -f processed.gro -o newbox.gro -c -d 1.2 -bt dodecahedron >> "$WORKER_LOG" 2>&1
            gmx solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top >> "$WORKER_LOG" 2>&1
            gmx grompp -f "$MDP_BASE/em.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
                -pname NA -nname CL -neutral -conc 0.15 >> "$WORKER_LOG" 2>&1
            cd "$OUT_BASE"
        fi

        # ---- Step 4: Energy Minimization ----
        if [ ! -f "$WORK/.done_em" ]; then
            cd "$WORK/em"
            gmx grompp -f "$MDP_BASE/em.mdp" -c "$WORK/topology/solv_ions.gro" \
                -p "$WORK/topology/topol.top" -o em.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            CUDA_VISIBLE_DEVICES=$GPU_ID gmx mdrun -v -deffnm em \
                -nb gpu -pme gpu >> "$WORKER_LOG" 2>&1
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
            gmx grompp -f "$MDP_BASE/nvt.mdp" -c "$WORK/em/em.gro" \
                -r "$WORK/em/em.gro" -p "$WORK/topology/topol.top" \
                -o nvt.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            CUDA_VISIBLE_DEVICES=$GPU_ID gmx mdrun -deffnm nvt \
                -nb gpu -pme gpu -bonded gpu -update gpu -pin on \
                >> "$WORKER_LOG" 2>&1
            [ $? -ne 0 ] && { echo "[GPU ${GPU_ID}] FAIL NVT: $VARIANT" >> "$WORKER_LOG"; JOB_FAIL=$((JOB_FAIL + 1)); continue; }

            # ---- Step 6: NPT Equilibration ----
            cd "$WORK/npt"
            gmx grompp -f "$MDP_BASE/npt.mdp" -c "$WORK/nvt/nvt.gro" \
                -r "$WORK/nvt/nvt.gro" -t "$WORK/nvt/nvt.cpt" \
                -p "$WORK/topology/topol.top" -o npt.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            CUDA_VISIBLE_DEVICES=$GPU_ID gmx mdrun -deffnm npt \
                -nb gpu -pme gpu -bonded gpu -update gpu -pin on \
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
            gmx grompp -f "$MDP_BASE/production.mdp" -c "$WORK/npt/npt.gro" \
                -t "$WORK/npt/npt.cpt" -p "$WORK/topology/topol.top" \
                -o md_prod.tpr -maxwarn 2 >> "$WORKER_LOG" 2>&1
            # 全 GPU offload，最大化 H200 利用率
            CUDA_VISIBLE_DEVICES=$GPU_ID GMX_CUDA_GRAPH=1 gmx mdrun \
                -deffnm md_prod \
                -nb gpu -pme gpu -bonded gpu -update gpu \
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

    run_md_worker "$i" "$LIST_FILE" "$OUTPUT_DIR" "$MDP_DIR" "$PHASE" &
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

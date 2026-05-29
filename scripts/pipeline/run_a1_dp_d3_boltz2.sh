#!/bin/bash
# ==============================================================================
# run_a1_dp_d3_boltz2.sh — VWF A1 + D'D3 Boltz-2 Parallel Run
# ==============================================================================
#
# 为 output/boltz2_a1_dp_d3/ 中的 YAML 文件批量运行 Boltz-2 预测。
#
# 研究目标：D'D3-A1 自抑制界面。Type 2B (GOF) 突变破坏 D'D3-A1 接触，
# 导致 A1 结构域自发结合血小板 GPIbα。
#
# 并行策略（单链 ~669 aa，最优配置待标定）：
#   N 个 worker 并行，每个 worker 占用 1 GPU，顺序跑分配到的 job。
#
# 断点续跑：每个 job 完成后立即写 .done 时间戳文件，重启自动跳过。
#
# 环境依赖：
#   conda activate boltz2
#   export CC=$(which gcc)
#
# 用法：
#   chmod +x run_a1_dp_d3_boltz2.sh
#   ./run_a1_dp_d3_boltz2.sh                    # 默认 4 GPU
#   ./run_a1_dp_d3_boltz2.sh --gpus 8           # 8 GPU 并行
#   ./run_a1_dp_d3_boltz2.sh --preflight        # 仅预检
#   ./run_a1_dp_d3_boltz2.sh --use-msa-server   # 服务器有网络
#
# 输出结构：
#   output/boltz2_a1_dp_d3_results/
#     boltz_results_VWF_WT_dp_d3_a1/
#       predictions/
#         model_0.cif
#         confidence_model_0.json
#       .done
#     worker_0.log ... worker_N.log
#     run_log.txt
# ==============================================================================

set -u

export CC=$(which gcc)

# ---------- 默认参数 -----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/../../output/boltz2_a1_dp_d3"
OUTPUT_DIR="${SCRIPT_DIR}/../../output/boltz2_a1_dp_d3_results"
N_GPUS=4
RECYCLING_STEPS=3
DIFFUSION_SAMPLES=5
USE_MSA_SERVER=false
PREFLIGHT_ONLY=false

# ---------- 参数解析 -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --gpus)              N_GPUS="$2"; shift 2 ;;
        --input-dir)         INPUT_DIR="$2"; shift 2 ;;
        --out-dir)           OUTPUT_DIR="$2"; shift 2 ;;
        --recycling-steps)  RECYCLING_STEPS="$2"; shift 2 ;;
        --diffusion-samples) DIFFUSION_SAMPLES="$2"; shift 2 ;;
        --use-msa-server)    USE_MSA_SERVER=true; shift ;;
        --preflight)         PREFLIGHT_ONLY=true; shift ;;
        -h|--help)
            sed -n '2,50p' "$0" | sed 's/^# //'
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
log "VWF A1+D'D3 Boltz-2 Parallel Run"
log "Started: $(date)"
log "============================================================"

PREFLIGHT_FAIL=false

# 1. boltz
if ! command -v boltz &>/dev/null; then
    log "[ERROR] 'boltz' not found."
    log "  pip install boltz  或  conda activate boltz2"
    PREFLIGHT_FAIL=true
else
    log "[OK] boltz found"
fi

# 2. CUDA / GPU
GPU_COUNT=$(python3 -c "import torch; print(torch.cuda.device_count())" 2>/dev/null || echo "0")
log "[INFO] CUDA GPUs detected: $GPU_COUNT"
if [ "$GPU_COUNT" -lt "$N_GPUS" ] 2>/dev/null; then
    log "[WARN] Adjusting N_GPUS from $N_GPUS to $GPU_COUNT"
    N_GPUS=$GPU_COUNT
fi
if [ "$GPU_COUNT" -eq 0 ]; then
    log "[ERROR] No CUDA GPUs found."
    PREFLIGHT_FAIL=true
fi

# 3. 输入目录
if [ ! -d "$INPUT_DIR" ]; then
    log "[ERROR] Input dir not found: $INPUT_DIR"
    log "  先运行: python3 scripts/pipeline/generate_a1_dp_d3_yamls.py"
    PREFLIGHT_FAIL=true
else
    YAML_COUNT=$(ls "$INPUT_DIR"/*.yaml 2>/dev/null | wc -l)
    log "[OK] Input dir: $INPUT_DIR ($YAML_COUNT YAML files)"
fi

# 4. 磁盘空间
AVAIL_MB=$(df -m "$OUTPUT_DIR" 2>/dev/null | awk 'NR==2{print $4}' || echo "unknown")
NEEDED_MB=$((YAML_COUNT * 250))
log "[INFO] Disk: available=${AVAIL_MB}MB, estimated_needed=${NEEDED_MB}MB"

# 5. MSA 模式
if $USE_MSA_SERVER; then
    log "[INFO] MSA mode: --use_msa_server"
else
    log "[INFO] MSA mode: offline (msa: empty)"
fi

log "------------------------------------------------------------"
log "Config: N_GPUS=$N_GPUS | recycling=$RECYCLING_STEPS | samples=$DIFFUSION_SAMPLES"
log "------------------------------------------------------------"

if $PREFLIGHT_FAIL; then
    log "[FATAL] Preflight failed."
    exit 1
fi

if $PREFLIGHT_ONLY; then
    log "Preflight complete. No jobs run."
    exit 0
fi

# ---------- 构建待跑 job 列表 --------------------------------------------------
TODO_YAMLS=()
ALL_YAMLS=( "$INPUT_DIR"/*.yaml )
TOTAL_JOBS=${#ALL_YAMLS[@]}

for YAML in "${ALL_YAMLS[@]}"; do
    JNAME=$(basename "$YAML" .yaml)
    DONE_MARKER="${OUTPUT_DIR}/boltz_results_${JNAME}/.done"
    if [ ! -f "$DONE_MARKER" ]; then
        TODO_YAMLS+=("$YAML")
    fi
done

DONE_COUNT=$((TOTAL_JOBS - ${#TODO_YAMLS[@]}))
log "Jobs: total=$TOTAL_JOBS  done=$DONE_COUNT  remaining=${#TODO_YAMLS[@]}"

if [ ${#TODO_YAMLS[@]} -eq 0 ]; then
    log "All jobs already completed."
    exit 0
fi

# ---------- 分配 job 到各 GPU（轮询） -----------------------------------------
TMPDIR_LISTS="${OUTPUT_DIR}/_worker_lists"
mkdir -p "$TMPDIR_LISTS"

for i in $(seq 0 $((N_GPUS - 1))); do
    > "${TMPDIR_LISTS}/gpu_${i}.txt"
done

IDX=0
for YAML in "${TODO_YAMLS[@]}"; do
    GPU_IDX=$((IDX % N_GPUS))
    echo "$YAML" >> "${TMPDIR_LISTS}/gpu_${GPU_IDX}.txt"
    IDX=$((IDX + 1))
done

log "Distribution across $N_GPUS GPUs:"
for i in $(seq 0 $((N_GPUS - 1))); do
    COUNT=$(wc -l < "${TMPDIR_LISTS}/gpu_${i}.txt" 2>/dev/null || echo 0)
    log "  GPU $i: $COUNT jobs"
done

# ---------- Worker 函数 --------------------------------------------------------
run_worker() {
    local GPU_ID=$1
    local JOB_LIST_FILE=$2
    local OUT_BASE=$3
    local RECYCLING=$4
    local SAMPLES=$5
    local USE_MSA=$6
    local WORKER_LOG="${OUT_BASE}/worker_${GPU_ID}.log"

    echo "[GPU ${GPU_ID}] Worker started at $(date)" > "$WORKER_LOG"
    echo "[GPU ${GPU_ID}] CUDA_VISIBLE_DEVICES=${GPU_ID}" >> "$WORKER_LOG"

    local JOB_OK=0
    local JOB_FAIL=0

    while IFS= read -r YAML_PATH; do
        [ -z "$YAML_PATH" ] && continue

        JNAME=$(basename "$YAML_PATH" .yaml)
        DONE_MARKER="${OUT_BASE}/boltz_results_${JNAME}/.done"

        if [ -f "$DONE_MARKER" ]; then
            echo "[GPU ${GPU_ID}] SKIP (done): $JNAME" >> "$WORKER_LOG"
            continue
        fi

        echo "[GPU ${GPU_ID}] $(date '+%H:%M:%S') START: $JNAME" >> "$WORKER_LOG"

        rm -f "${OUT_BASE}/boltz_results_${JNAME}/predictions/"*.cif \
              "${OUT_BASE}/boltz_results_${JNAME}/predictions/"confidence_*.json \
              "${OUT_BASE}/boltz_results_${JNAME}/predictions/"affinity_*.json 2>/dev/null || true

        local MSA_FLAG=""
        if [ "$USE_MSA" = "true" ]; then
            MSA_FLAG="--use_msa_server"
        fi

        CUDA_VISIBLE_DEVICES=$GPU_ID CC=$(which gcc) boltz predict "$YAML_PATH" \
            --out_dir "$OUT_BASE" \
            --accelerator gpu \
            --devices 1 \
            --recycling_steps "$RECYCLING" \
            --diffusion_samples "$SAMPLES" \
            --num_workers 0 \
            --override \
            $MSA_FLAG \
            >> "$WORKER_LOG" 2>&1

        EXIT_CODE=$?

        if [ $EXIT_CODE -eq 0 ]; then
            date '+%Y-%m-%d %H:%M:%S' > "$DONE_MARKER"
            JOB_OK=$((JOB_OK + 1))
            echo "[GPU ${GPU_ID}] $(date '+%H:%M:%S') OK: $JNAME" >> "$WORKER_LOG"
        else
            JOB_FAIL=$((JOB_FAIL + 1))
            echo "[GPU ${GPU_ID}] $(date '+%H:%M:%S') FAIL (exit=$EXIT_CODE): $JNAME" >> "$WORKER_LOG"
        fi

    done < "$JOB_LIST_FILE"

    echo "[GPU ${GPU_ID}] Worker done at $(date) | OK=$JOB_OK FAIL=$JOB_FAIL" >> "$WORKER_LOG"
}

export -f run_worker

# ---------- 启动并行 Workers --------------------------------------------------
log "Launching $N_GPUS parallel GPU workers..."
PIDS=()

for i in $(seq 0 $((N_GPUS - 1))); do
    LIST_FILE="${TMPDIR_LISTS}/gpu_${i}.txt"
    JOB_COUNT=$(wc -l < "$LIST_FILE" 2>/dev/null || echo 0)
    [ "$JOB_COUNT" -eq 0 ] && continue

    run_worker "$i" "$LIST_FILE" "$OUTPUT_DIR" \
        "$RECYCLING_STEPS" "$DIFFUSION_SAMPLES" "$USE_MSA_SERVER" &
    PIDS+=($!)
    log "  Worker GPU-$i started (PID=${PIDS[-1]}, jobs=$JOB_COUNT)"
done

log ""
log "Workers running. Monitor with:"
log "  watch -n 30 'find $OUTPUT_DIR -name .done | wc -l'"
log ""

# ---------- 等待所有 Worker 完成 ----------------------------------------------
for PID in "${PIDS[@]}"; do
    wait "$PID"
    RC=$?
    [ $RC -ne 0 ] && log "[WARN] Worker PID=$PID exited with code $RC"
done

FINAL_DONE=$(find "$OUTPUT_DIR" -name .done 2>/dev/null | wc -l)
log ""
log "============================================================"
log "Boltz-2 Run finished: $(date)"
log "  Completed: $FINAL_DONE / $TOTAL_JOBS"
log "  Output dir: $OUTPUT_DIR"
log "============================================================"

rm -rf "$TMPDIR_LISTS"
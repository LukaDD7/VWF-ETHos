#!/bin/bash
# ==============================================================================
# run_a1_gpiba_boltz2.sh — VWF A1 + GPIbα Boltz-2 Parallel Run
# ==============================================================================
#
# 为 output/boltz2_a1_gpiba/ 中的 74 个 YAML 文件批量运行 Boltz-2 预测。
#
# 并行策略（针对小蛋白 489aa 最优）：
#   N 个 worker 并行，每个 worker 占用 1 GPU，顺序跑分配到的 job。
#   不使用 DDP（--strategy ddp），因为对小蛋白 DDP 通信开销 > 加速收益。
#
# 断点续跑：每个 job 完成后立即写 .done 时间戳文件，重启自动跳过。
#
# Triton JIT 缓存：首次运行需等待 kernel 编译（~5-10min/job）。编译产物缓存在
#   ~/.cache/triton/，复制到 GPU 实例可复用，省去重复编译时间。
#
# 环境依赖：
#   conda activate boltz2
#   gcc（系统级安装，boltz2 环境内可用 which gcc 确认）
#
# 用法：
#   chmod +x run_a1_gpiba_boltz2.sh
#   ./run_a1_gpiba_boltz2.sh                    # 默认 4 GPU
#   ./run_a1_gpiba_boltz2.sh --gpus 8           # 8 GPU 并行
#   ./run_a1_gpiba_boltz2.sh --gpus 4 --preflight  # 仅预检，不运行
#   ./run_a1_gpiba_boltz2.sh --use-msa-server   # 服务器有网络时启用
#
# 输出结构（每个 job）：
#   output/boltz2_a1_gpiba_results/
#     boltz_results_VWF_R1306W_vs_GPIb_alpha/
#       predictions/
#         model_0.cif
#         confidence_model_0.json   ← 包含 iPTM
#       .done                        ← 完成时间戳
#     worker_0.log ... worker_N.log
#     run_log.txt
# ==============================================================================

set -u

export CC=$(which gcc)

# ---------- 默认参数 -----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INPUT_DIR="${SCRIPT_DIR}/../../output/boltz2_a1_gpiba"
OUTPUT_DIR="${SCRIPT_DIR}/../../output/boltz2_a1_gpiba_results"
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
        --recycling-steps)   RECYCLING_STEPS="$2"; shift 2 ;;
        --diffusion-samples) DIFFUSION_SAMPLES="$2"; shift 2 ;;
        --use-msa-server)    USE_MSA_SERVER=true; shift ;;
        --preflight)         PREFLIGHT_ONLY=true; shift ;;
        -h|--help)
            sed -n '2,40p' "$0" | sed 's/^# //'
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
log "VWF A1+GPIbα Boltz-2 Parallel Run"
log "Started: $(date)"
log "============================================================"

PREFLIGHT_FAIL=false

# 1. boltz 是否安装
if ! command -v boltz &>/dev/null; then
    log "[ERROR] 'boltz' not found."
    log "  安装命令: pip install boltz"
    log "  或: conda activate <env> && pip install boltz"
    PREFLIGHT_FAIL=true
else
    BOLTZ_VER=$(boltz --version 2>/dev/null || boltz version 2>/dev/null || echo "unknown")
    log "[OK] boltz found: $BOLTZ_VER"
fi

# 2. CUDA / GPU
GPU_COUNT=$(python3 -c "import torch; print(torch.cuda.device_count())" 2>/dev/null || echo "0")
log "[INFO] CUDA GPUs detected: $GPU_COUNT"
if [ "$GPU_COUNT" -lt "$N_GPUS" ] 2>/dev/null; then
    log "[WARN] Requested $N_GPUS GPUs but only $GPU_COUNT detected."
    log "  Adjusting N_GPUS to $GPU_COUNT"
    N_GPUS=$GPU_COUNT
fi
if [ "$GPU_COUNT" -eq 0 ]; then
    log "[ERROR] No CUDA GPUs found. Cannot run."
    PREFLIGHT_FAIL=true
fi

# 3. 输入目录
if [ ! -d "$INPUT_DIR" ]; then
    log "[ERROR] Input dir not found: $INPUT_DIR"
    PREFLIGHT_FAIL=true
else
    YAML_COUNT=$(ls "$INPUT_DIR"/*.yaml 2>/dev/null | wc -l)
    log "[OK] Input dir: $INPUT_DIR ($YAML_COUNT YAML files)"
fi

# 4. 磁盘空间（粗估：每个 job ~200MB）
AVAIL_MB=$(df -m "$OUTPUT_DIR" 2>/dev/null | awk 'NR==2{print $4}' || echo "unknown")
NEEDED_MB=$((YAML_COUNT * 200))
log "[INFO] Disk: available=${AVAIL_MB}MB, estimated_needed=${NEEDED_MB}MB"
if [ "$AVAIL_MB" != "unknown" ] && [ "$AVAIL_MB" -lt "$NEEDED_MB" ] 2>/dev/null; then
    log "[WARN] Low disk space. May fail partway through."
fi

# 5. MSA 设置
if $USE_MSA_SERVER; then
    log "[INFO] MSA mode: --use_msa_server (requires internet)"
else
    log "[INFO] MSA mode: offline (msa: empty in YAMLs, single-sequence)"
fi

log "------------------------------------------------------------"
log "Config: N_GPUS=$N_GPUS | recycling=$RECYCLING_STEPS | samples=$DIFFUSION_SAMPLES"
log "------------------------------------------------------------"

if $PREFLIGHT_FAIL; then
    log "[FATAL] Preflight failed. Fix errors above before launching GPU instance."
    exit 1
fi

if $PREFLIGHT_ONLY; then
    log "Preflight complete (--preflight mode). No jobs run."
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
    log "All jobs already completed. Nothing to do."
    log "  To re-run: rm -f output/boltz2_a1_gpiba_results/<name>/.done"
    exit 0
fi

# ---------- 分配 job 到各 GPU（轮询） -----------------------------------------
# 将 YAML 列表写到临时文件，每个 GPU 一个
TMPDIR_LISTS="${OUTPUT_DIR}/_worker_lists"
mkdir -p "$TMPDIR_LISTS"

# 清空旧的分配文件
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

# ---------- Worker 函数（在后台子进程中运行）----------------------------------
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

        # 再次检查（防止并发重复）
        if [ -f "$DONE_MARKER" ]; then
            echo "[GPU ${GPU_ID}] SKIP (done): $JNAME" >> "$WORKER_LOG"
            continue
        fi

        echo "[GPU ${GPU_ID}] $(date '+%H:%M:%S') START: $JNAME" >> "$WORKER_LOG"

        # 清理上次可能的脏文件
        rm -f "${OUT_BASE}/boltz_results_${JNAME}/predictions/"*.cif \
              "${OUT_BASE}/boltz_results_${JNAME}/predictions/"confidence_*.json \
              "${OUT_BASE}/boltz_results_${JNAME}/predictions/"affinity_*.json 2>/dev/null || true

        # 构造 boltz 命令
        # 注意：--out_dir 设为 OUT_BASE，boltz 会自动创建 JNAME/ 子目录
        # 输出结构：OUT_BASE/JNAME/predictions/confidence_model_*.json
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
            echo "[GPU ${GPU_ID}] -- 如需单独重试: CUDA_VISIBLE_DEVICES=$GPU_ID boltz predict $YAML_PATH --out_dir $OUT_BASE --accelerator gpu --devices 1 --override" >> "$WORKER_LOG"
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
log "  或: python3 scripts/pipeline/parse_a1_gpiba_results.py --check-progress"
log ""

# ---------- 等待所有 Worker 完成 ----------------------------------------------
ALL_WORKERS_OK=true
for PID in "${PIDS[@]}"; do
    wait "$PID"
    RC=$?
    if [ $RC -ne 0 ]; then
        ALL_WORKERS_OK=false
        log "[WARN] Worker PID=$PID exited with code $RC"
    fi
done

# ---------- 最终汇总 ----------------------------------------------------------
FINAL_DONE=$(find "$OUTPUT_DIR" -name ".done" 2>/dev/null | wc -l)
FINAL_FAIL=$((TOTAL_JOBS - FINAL_DONE))

log ""
log "============================================================"
log "Run finished: $(date)"
log "  Total jobs : $TOTAL_JOBS"
log "  Completed  : $FINAL_DONE"
log "  Failed/skip: $FINAL_FAIL"
log "  Output dir : $OUTPUT_DIR"
log "============================================================"

if [ "$FINAL_DONE" -eq "$TOTAL_JOBS" ]; then
    log "SUCCESS: All $TOTAL_JOBS jobs completed."
    log ""
    log "下一步 — 解析结果:"
    log "  python3 scripts/pipeline/parse_a1_gpiba_results.py \\"
    log "      --results-dir $OUTPUT_DIR \\"
    log "      --output output/boltz2_a1_gpiba_analysis/iptm_results.csv"
else
    log "PARTIAL: $FINAL_FAIL jobs failed or not done."
    log "  断点续跑（重新执行此脚本即可，已完成的 job 自动跳过）:"
    log "  ./run_a1_gpiba_boltz2.sh --gpus $N_GPUS"
    log ""
    log "  查看失败原因:"
    log "  grep FAIL $OUTPUT_DIR/worker_*.log"
fi

# 清理 worker list 临时文件
rm -rf "$TMPDIR_LISTS"

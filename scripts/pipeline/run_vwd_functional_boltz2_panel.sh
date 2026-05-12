#!/bin/bash
# ==============================================================================
# run_vwd_functional_boltz2_panel.sh — VWD Functional Panel Boltz-2 Parallel Run
# ==============================================================================
#
# 为 output/boltz2_vwd_functional_panel/yamls/ 中的 990 个 YAML 文件批量运行
# Boltz-2 预测。包含 FVIII binding (2N axis)、GPIbα (2B/2M axis)、
# ADAMTS13 (2A axis)、Collagen (2M-A3 axis) 等全部 15 个 assay 轴。
#
# 并行策略（针对 4/8 × H200 141GB 最优）：
#   N 个 worker 并行，每个 worker 占用 1 GPU，顺序跑分配到的 job。
#   不使用 DDP（--strategy ddp），因为单 job 规模 < 1000 aa，
#   DDP 通信开销 > 加速收益。单卡 H200 跑 ~500aa 构象约 3-5 min/sample。
#
# 断点续跑：每个 job 完成后立即写 .done 时间戳文件，重启自动跳过。
#
# GPU 实例约束：
#   - 不可联网 → 所有 YAML 中 msa: empty（单序列模式），不使用 --use_msa_server
#   - Triton JIT 首次编译 ~5-10 min → 缓存 ~/.cache/triton/，可跨实例复制复用
#   - /dev/shm 通常仅 64MB → 强制 --num_workers 0 避免 DataLoader IPC 崩溃
#   - 每卡 ~15-25GB GPU 内存 → 单卡单 job 足够，无需担心 OOM
#
# 环境依赖（在 CPU 实例预装，打包后传输到 GPU 实例）：
#   conda activate boltz2
#   pip install boltz>=2.0
#   gcc（系统级安装，Triton JIT 需要）
#
# 用法：
#   chmod +x scripts/pipeline/run_vwd_functional_boltz2_panel.sh
#   bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh                # 默认 4 GPU
#   bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --gpus 8       # 8 GPU 并行
#   bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --gpus 4 --preflight  # 仅预检
#   bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --filter fviii # 只跑 FVIII assay
#
# 输出结构（每个 job）：
#   output/boltz2_vwd_functional_panel/boltz_results/
#     boltz_results_VWF_R816W__dprime_d3_fviii_binding/
#       predictions/
#         model_0.cif
#         confidence_model_0.json   ← pLDDT, iPTM, pTM
#         affinity_model_0.json     ← ΔG (Boltz-2 v2+)
#       .done                        ← 完成时间戳
#     worker_0.log ... worker_N.log
#     run_log.txt
# ==============================================================================

set -u

export CC=$(which gcc 2>/dev/null || echo "/usr/bin/gcc")

# ---------- 默认参数 -----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
INPUT_DIR="${ROOT_DIR}/output/boltz2_vwd_functional_panel/yamls"
OUTPUT_DIR="${ROOT_DIR}/output/boltz2_vwd_functional_panel/boltz_results"
N_GPUS=4
RECYCLING_STEPS=3
DIFFUSION_SAMPLES=5
PREFLIGHT_ONLY=false
ASSAY_FILTER=""   # 空=全部, "fviii"=仅FVIII, "gpiba"=仅GPIba

# ---------- 参数解析 -----------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --gpus)              N_GPUS="$2"; shift 2 ;;
        --input-dir)         INPUT_DIR="$2"; shift 2 ;;
        --out-dir)           OUTPUT_DIR="$2"; shift 2 ;;
        --recycling-steps)   RECYCLING_STEPS="$2"; shift 2 ;;
        --diffusion-samples) DIFFUSION_SAMPLES="$2"; shift 2 ;;
        --preflight)         PREFLIGHT_ONLY=true; shift ;;
        --filter)            ASSAY_FILTER="$2"; shift 2 ;;
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
log "VWD Functional Panel — Boltz-2 Parallel Run"
log "Started: $(date)"
log "============================================================"

PREFLIGHT_FAIL=false

# 1. boltz 是否安装
if ! command -v boltz &>/dev/null; then
    log "[ERROR] 'boltz' not found."
    log "  安装命令: pip install boltz"
    PREFLIGHT_FAIL=true
else
    BOLTZ_VER=$(boltz --version 2>/dev/null || echo "unknown")
    log "[OK] boltz found: $BOLTZ_VER"
fi

# 2. CUDA / GPU
GPU_COUNT=$(python3 -c "import torch; print(torch.cuda.device_count())" 2>/dev/null || echo "0")
log "[INFO] CUDA GPUs detected: $GPU_COUNT"
if [ "$GPU_COUNT" -lt "$N_GPUS" ] 2>/dev/null; then
    log "[WARN] Requested $N_GPUS GPUs but only $GPU_COUNT detected. Adjusting."
    N_GPUS=$GPU_COUNT
fi
if [ "$GPU_COUNT" -eq 0 ]; then
    log "[ERROR] No CUDA GPUs found. Cannot run."
    PREFLIGHT_FAIL=true
fi

# 3. GPU 内存检查
if [ "$GPU_COUNT" -gt 0 ]; then
    python3 -c "
import torch
for i in range(torch.cuda.device_count()):
    props = torch.cuda.get_device_properties(i)
    gb = props.total_mem / 1e9
    print(f'  GPU {i}: {props.name}  {gb:.0f} GB')
" 2>/dev/null | while read line; do log "$line"; done
fi

# 4. 输入目录
if [ ! -d "$INPUT_DIR" ]; then
    log "[ERROR] Input dir not found: $INPUT_DIR"
    log "  生成命令: python scripts/pipeline/generate_vwd_functional_boltz2_yamls.py --write-json-batches"
    PREFLIGHT_FAIL=true
else
    YAML_COUNT=$(ls "$INPUT_DIR"/*.yaml 2>/dev/null | wc -l | tr -d ' ')
    log "[OK] Input dir: $INPUT_DIR ($YAML_COUNT YAML files)"
fi

# 5. 磁盘空间（粗估：每个 job ~200MB）
AVAIL_MB=$(df -m "$OUTPUT_DIR" 2>/dev/null | awk 'NR==2{print $4}' || echo "unknown")
NEEDED_MB=$((YAML_COUNT * 200))
log "[INFO] Disk: available=${AVAIL_MB}MB, estimated_needed=${NEEDED_MB}MB"
if [ "$AVAIL_MB" != "unknown" ] && [ "$AVAIL_MB" -lt "$NEEDED_MB" ] 2>/dev/null; then
    log "[WARN] Low disk space. May fail partway through."
fi

# 6. gcc (Triton JIT)
if ! command -v gcc &>/dev/null; then
    log "[WARN] gcc not found. Triton JIT may fail."
    log "  Fix: apt install gcc  or  export CC=/path/to/gcc"
else
    log "[OK] gcc found: $(which gcc)"
fi

log "------------------------------------------------------------"
log "Config: N_GPUS=$N_GPUS | recycling=$RECYCLING_STEPS | samples=$DIFFUSION_SAMPLES"
log "Filter: ${ASSAY_FILTER:-all assays}"
log "------------------------------------------------------------"

if $PREFLIGHT_FAIL; then
    log "[FATAL] Preflight failed. Fix errors above before launching GPU instance."
    exit 1
fi

if $PREFLIGHT_ONLY; then
    log "Preflight complete (--preflight mode). No jobs run."
    exit 0
fi

# ---------- 构建待跑 job 列表（支持 assay 过滤）---------------------------------
TODO_YAMLS=()
ALL_YAMLS=( "$INPUT_DIR"/*.yaml )
TOTAL_JOBS=${#ALL_YAMLS[@]}

for YAML in "${ALL_YAMLS[@]}"; do
    JNAME=$(basename "$YAML" .yaml)

    # 过滤 assay 类型
    if [ -n "$ASSAY_FILTER" ]; then
        case "$ASSAY_FILTER" in
            fviii|FVIII|f8|F8)
                [[ "$JNAME" != *"dprime_d3_fviii_binding"* ]] && continue ;;
            gpiba|GPIba|a1)
                [[ "$JNAME" != *"a1_gpiba_forced_binding"* ]] && continue ;;
            adamts13|a2)
                [[ "$JNAME" != *"adamts13"* && "$JNAME" != *"a2_"* && "$JNAME" != *"vwf73_"* ]] && continue ;;
            collagen|a3)
                [[ "$JNAME" != *"collagen"* ]] && continue ;;
            *)
                [[ "$JNAME" != *"$ASSAY_FILTER"* ]] && continue ;;
        esac
    fi

    DONE_MARKER="${OUTPUT_DIR}/boltz_results_${JNAME}/.done"
    if [ ! -f "$DONE_MARKER" ]; then
        TODO_YAMLS+=("$YAML")
    fi
done

DONE_COUNT=$((TOTAL_JOBS - ${#TODO_YAMLS[@]}))
log "Jobs: total=$TOTAL_JOBS  done=$DONE_COUNT  remaining=${#TODO_YAMLS[@]}"
log "  (filtered by: ${ASSAY_FILTER:-none})"

if [ ${#TODO_YAMLS[@]} -eq 0 ]; then
    log "All matching jobs already completed. Nothing to do."
    log "  To re-run: rm -f $OUTPUT_DIR/<name>/.done"
    exit 0
fi

# ---------- 分配 job 到各 GPU（轮询）------------------------------------------
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

# ---------- Worker 函数（在后台子进程中运行）-----------------------------------
run_worker() {
    local GPU_ID=$1
    local JOB_LIST_FILE=$2
    local OUT_BASE=$3
    local RECYCLING=$4
    local SAMPLES=$5
    local WORKER_LOG="${OUT_BASE}/worker_${GPU_ID}.log"

    echo "[GPU ${GPU_ID}] Worker started at $(date)" > "$WORKER_LOG"
    echo "[GPU ${GPU_ID}] CUDA_VISIBLE_DEVICES=${GPU_ID}" >> "$WORKER_LOG"

    local JOB_OK=0
    local JOB_FAIL=0
    local JOB_TOTAL=$(wc -l < "$JOB_LIST_FILE" 2>/dev/null || echo 0)
    local JOB_IDX=0

    while IFS= read -r YAML_PATH; do
        [ -z "$YAML_PATH" ] && continue

        JOB_IDX=$((JOB_IDX + 1))
        JNAME=$(basename "$YAML_PATH" .yaml)
        DONE_MARKER="${OUT_BASE}/boltz_results_${JNAME}/.done"

        # 再次检查（防止并发重复）
        if [ -f "$DONE_MARKER" ]; then
            echo "[GPU ${GPU_ID}] SKIP (done): $JNAME" >> "$WORKER_LOG"
            continue
        fi

        echo "[GPU ${GPU_ID}] $(date '+%H:%M:%S') START [$JOB_IDX/$JOB_TOTAL]: $JNAME" >> "$WORKER_LOG"

        # 清理上次可能的脏文件
        rm -f "${OUT_BASE}/boltz_results_${JNAME}/predictions/"*.cif \
              "${OUT_BASE}/boltz_results_${JNAME}/predictions/"confidence_*.json \
              "${OUT_BASE}/boltz_results_${JNAME}/predictions/"affinity_*.json 2>/dev/null || true

        # Boltz-2 预测
        # --num_workers 0 : 避免 /dev/shm 64MB 限制导致的 DataLoader IPC 崩溃
        # --devices 1     : 单卡单 job（不用 DDP，对 <1000aa 系统无收益）
        # --override      : 覆盖已有的不完整输出
        # 不加 --use_msa_server : GPU 实例无网络
        CUDA_VISIBLE_DEVICES=$GPU_ID CC=$(which gcc) boltz predict "$YAML_PATH" \
            --out_dir "$OUT_BASE" \
            --accelerator gpu \
            --devices 1 \
            --recycling_steps "$RECYCLING" \
            --diffusion_samples "$SAMPLES" \
            --num_workers 0 \
            --override \
            >> "$WORKER_LOG" 2>&1

        EXIT_CODE=$?

        if [ $EXIT_CODE -eq 0 ]; then
            # 写 .done 标记（含时间戳）
            mkdir -p "${OUT_BASE}/boltz_results_${JNAME}"
            date '+%Y-%m-%d %H:%M:%S' > "$DONE_MARKER"
            JOB_OK=$((JOB_OK + 1))
            echo "[GPU ${GPU_ID}] $(date '+%H:%M:%S') OK [$JOB_IDX/$JOB_TOTAL]: $JNAME" >> "$WORKER_LOG"
        else
            JOB_FAIL=$((JOB_FAIL + 1))
            echo "[GPU ${GPU_ID}] $(date '+%H:%M:%S') FAIL (exit=$EXIT_CODE): $JNAME" >> "$WORKER_LOG"
            echo "[GPU ${GPU_ID}] -- 重试: CUDA_VISIBLE_DEVICES=$GPU_ID boltz predict $YAML_PATH --out_dir $OUT_BASE --accelerator gpu --devices 1 --override" >> "$WORKER_LOG"
        fi

    done < "$JOB_LIST_FILE"

    echo "[GPU ${GPU_ID}] Worker done at $(date) | OK=$JOB_OK FAIL=$JOB_FAIL TOTAL=$JOB_TOTAL" >> "$WORKER_LOG"
}

export -f run_worker

# ---------- 启动并行 Workers ---------------------------------------------------
log "Launching $N_GPUS parallel GPU workers..."
PIDS=()

for i in $(seq 0 $((N_GPUS - 1))); do
    LIST_FILE="${TMPDIR_LISTS}/gpu_${i}.txt"
    JOB_COUNT=$(wc -l < "$LIST_FILE" 2>/dev/null || echo 0)
    [ "$JOB_COUNT" -eq 0 ] && continue

    run_worker "$i" "$LIST_FILE" "$OUTPUT_DIR" \
        "$RECYCLING_STEPS" "$DIFFUSION_SAMPLES" &
    PIDS+=($!)
    log "  Worker GPU-$i started (PID=${PIDS[-1]}, jobs=$JOB_COUNT)"
done

log ""
log "Workers running. Monitor with:"
log "  watch -n 30 'find $OUTPUT_DIR -name .done | wc -l'"
log "  tail -f $OUTPUT_DIR/worker_0.log"
log ""

# ---------- 估算时间 -----------------------------------------------------------
# H200 141GB 单卡: ~3 min/sample × 5 samples × 5 min overhead = ~20 min/job
# 4 GPU: 990 jobs / 4 = 248 jobs/GPU ≈ 83h
# 8 GPU: 990 jobs / 8 = 124 jobs/GPU ≈ 41h
EST_MIN_PER_JOB=20
EST_TOTAL=$((${#TODO_YAMLS[@]} * EST_MIN_PER_JOB / N_GPUS / 60))
log "Estimated time: ~${EST_TOTAL}h (${#TODO_YAMLS[@]} jobs / $N_GPUS GPUs, ~${EST_MIN_PER_JOB}min/job)"
log ""

# ---------- 等待所有 Worker 完成 -----------------------------------------------
ALL_WORKERS_OK=true
for PID in "${PIDS[@]}"; do
    wait "$PID"
    RC=$?
    if [ $RC -ne 0 ]; then
        ALL_WORKERS_OK=false
        log "[WARN] Worker PID=$PID exited with code $RC"
    fi
done

# ---------- 最终汇总 -----------------------------------------------------------
FINAL_DONE=$(find "$OUTPUT_DIR" -name ".done" 2>/dev/null | wc -l)
TOTAL_YAML=$(ls "$INPUT_DIR"/*.yaml 2>/dev/null | wc -l | tr -d ' ')
FINAL_FAIL=$((TOTAL_YAML - FINAL_DONE))

log ""
log "============================================================"
log "Run finished: $(date)"
log "  Total YAML files : $TOTAL_YAML"
log "  Completed (.done): $FINAL_DONE"
log "  Remaining        : $FINAL_FAIL"
log "  Output dir       : $OUTPUT_DIR"
log "============================================================"

if [ "$FINAL_DONE" -ge "$TOTAL_YAML" ]; then
    log "SUCCESS: All jobs completed."
    log ""
    log "下一步 — 解析结果:"
    log "  python3 scripts/pipeline/parse_vwd_functional_boltz2_results.py \\"
    log "      --results-dir $OUTPUT_DIR \\"
    log "      --manifest output/boltz2_vwd_functional_panel/job_manifest.csv"
else
    log "PARTIAL: $FINAL_FAIL jobs remaining."
    log "  断点续跑（已完成的 job 自动跳过）:"
    log "  bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --gpus $N_GPUS"
    log ""
    log "  查看失败原因:"
    log "  grep FAIL $OUTPUT_DIR/worker_*.log"
fi

# 清理 worker list 临时文件
rm -rf "$TMPDIR_LISTS"

#!/bin/bash
# ==============================================================================
# VWF Boltz-2 Production Run Script  (v2 — Job-Level Incremental)
# GPU: 4x NVIDIA H200 (80GB HBM3 each)
#
# 支持优先级分批：
#   ./run_boltz2.sh                              # 默认：跑全部 blind_scan
#   ./run_boltz2.sh --input-dir priority1_2B2M   # 只跑 2B/2M（相对 blind_scan 目录）
#   ./run_boltz2.sh --input-dir priority2_2A2N   # 只跑 2A/2N
#
# 增量存储机制（防断电 / 防崩溃）：
#   - 每个 job 完成后立即写 <job_dir>/.done 时间戳文件
#   - 重新执行脚本自动跳过已有 .done 的 job（job 级，非 batch 级）
#   - 崩溃只丢失正在跑的单个 job，所有已完成 job 全部保留
#   - 进度实时写入 progress.json，可独立监控
#
# 输出结构：
#   output/boltz2_results/
#     VWF_R1306W_vs_GPIb_alpha/
#       predictions/
#         model_0.cif
#         confidence_model_0.json
#         affinity_model_0.json    (Boltz-2 v2+)
#       .done                      ← job 级完成标记
#     _tmp_single_jobs/            ← 临时单 job YAML（Boltz 标准格式）
#     run_log.txt
#     progress.json
#
# 注意：Boltz CLI 只接受 .yaml 或 .fasta，不接受 JSON。
#       run_boltz2.sh 在展开阶段自动将 batch JSON → Boltz YAML。
# ==============================================================================

# 不用 set -e：单个 job 失败不终止整个流程
set -u

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
BLIND_SCAN_BASE="${PIPELINE_DIR}/../../output/boltz2_blind_scan"

# ---- 解析 --input-dir 参数（相对于 blind_scan_base 或绝对路径）-----------
SUBDIR=""
for arg in "$@"; do
    case $arg in
        --input-dir=*) SUBDIR="${arg#*=}" ;;
        --input-dir)   shift; SUBDIR="$1" ;;
    esac
done

if [ -n "${SUBDIR}" ]; then
    # 如果是绝对路径直接用；否则视为 blind_scan_base 的子目录
    if [[ "${SUBDIR}" = /* ]]; then
        INPUT_DIR="${SUBDIR}"
    else
        INPUT_DIR="${BLIND_SCAN_BASE}/${SUBDIR}"
    fi
    # 结果输出目录以子目录名区分
    RUN_LABEL=$(basename "${INPUT_DIR}")
else
    INPUT_DIR="${BLIND_SCAN_BASE}"
    RUN_LABEL="full"
fi

OUTPUT_DIR="${PIPELINE_DIR}/../../output/boltz2_results/${RUN_LABEL}"
TEMP_DIR="${OUTPUT_DIR}/_tmp_single_jobs"
LOG_FILE="${OUTPUT_DIR}/run_log.txt"
PROGRESS_FILE="${OUTPUT_DIR}/progress.json"

# Boltz-2 参数（4x H200）
N_GPUS=4
RECYCLING_STEPS=3
DIFFUSION_SAMPLES=5
NUM_WORKERS=8

# ---- 工具 -------------------------------------------------------------------
log() {
    local msg="[$(date '+%H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOG_FILE"
}

write_progress() {
    # 写出实时进度 JSON
    python3 -c "
import json, time
json.dump({
    'completed': $1,
    'total': $2,
    'percent': round($1/$2*100, 1) if $2 > 0 else 0,
    'current_job': '$3',
    'last_update': time.strftime('%Y-%m-%d %H:%M:%S')
}, open('${PROGRESS_FILE}', 'w'), indent=2)
" 2>/dev/null || true
}

# ---- 环境检查 ---------------------------------------------------------------
echo "============================================================"
echo "VWF Boltz-2 Production Run (Job-Level Incremental)"
echo "Started: $(date)"
echo "============================================================"

if ! command -v boltz &>/dev/null; then
    echo "[ERROR] 'boltz' not found. Run: pip install boltz"
    exit 1
fi
if [ ! -d "$INPUT_DIR" ]; then
    echo "[ERROR] Input dir not found: $INPUT_DIR"
    echo "  Run prepare_boltz2_inputs.py first."
    exit 1
fi

mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

# ---- Step 1: 展开 batch JSON → 单 job YAML（Boltz 标准格式，幂等）----------
log "Expanding batches → single-job YAMLs in $TEMP_DIR ..."

python3 - "$INPUT_DIR" "$TEMP_DIR" <<'PYEOF'
import json, os, glob, sys
inp, tmp = sys.argv[1], sys.argv[2]

def job_to_yaml(job):
    """将 job dict 转换为 Boltz CLI 接受的 YAML 字符串。"""
    lines = ['version: 1', 'sequences:']
    for seq in job.get('sequences', []):
        p = seq.get('protein', {})
        lines.append('  - protein:')
        lines.append(f'      id: {p["id"]}')
        lines.append(f'      sequence: {p["sequence"]}')
    # affinity property
    for prop in job.get('properties', []):
        if 'affinity' in prop:
            lines.append('properties:')
            lines.append('  affinity:')
            lines.append(f'    binder: {prop["affinity"]["binder"]}')
            break
    return '\n'.join(lines) + '\n'

written = 0
for bf in sorted(glob.glob(os.path.join(inp, 'batch_*.json'))):
    data = json.load(open(bf))
    for job in data.get('jobs', []):
        jname = job['name']
        jpath = os.path.join(tmp, jname + '.yaml')   # ← YAML 格式
        if not os.path.exists(jpath):
            with open(jpath, 'w') as f:
                f.write(job_to_yaml(job))
            written += 1

total = len(glob.glob(os.path.join(tmp, '*.yaml')))
print(f"  Expanded {written} new  |  Total single-job YAML files: {total}")
PYEOF

# ---- Step 2: 统计未完成 job --------------------------------------------------
ALL_JOB_YAMLS=( "$TEMP_DIR"/*.yaml )
TOTAL_JOBS=${#ALL_JOB_YAMLS[@]}
TODO_YAMLS=()

for jy in "${ALL_JOB_YAMLS[@]}"; do
    # job 名 = YAML 文件名（不含扩展名）
    JN=$(basename "$jy" .yaml)
    [ -z "$JN" ] && continue
    if [ ! -f "$OUTPUT_DIR/$JN/.done" ]; then
        TODO_YAMLS+=("$jy")
    fi
done

DONE_COUNT=$((TOTAL_JOBS - ${#TODO_YAMLS[@]}))
log "Total: $TOTAL_JOBS  |  Done: $DONE_COUNT  |  Remaining: ${#TODO_YAMLS[@]}"
echo ""

if [ ${#TODO_YAMLS[@]} -eq 0 ]; then
    log "All jobs completed. Nothing to do."
    write_progress "$TOTAL_JOBS" "$TOTAL_JOBS" "COMPLETE"
    exit 0
fi

write_progress "$DONE_COUNT" "$TOTAL_JOBS" "starting..."

# ---- Step 3: 逐 job 执行（核心：job 级增量存储）-----------------------------
log "Starting prediction loop..."
CURRENT=0
FAILED=0

for JOB_YAML in "${TODO_YAMLS[@]}"; do
    JOB_NAME=$(basename "$JOB_YAML" .yaml)
    [ -z "$JOB_NAME" ] && continue

    JOB_DIR="$OUTPUT_DIR/$JOB_NAME"
    DONE_MARKER="$JOB_DIR/.done"

    # 再次检查（防止并发场景）
    [ -f "$DONE_MARKER" ] && { DONE_COUNT=$((DONE_COUNT+1)); continue; }

    CURRENT=$((CURRENT + 1))
    log "[$((DONE_COUNT + CURRENT))/$TOTAL_JOBS] $JOB_NAME"
    write_progress "$((DONE_COUNT + CURRENT - 1))" "$TOTAL_JOBS" "$JOB_NAME"

    mkdir -p "$JOB_DIR"

    # 清理上次崩溃留下的脏文件
    rm -f "$JOB_DIR/predictions/"*.cif \
          "$JOB_DIR/predictions/"confidence_*.json \
          "$JOB_DIR/predictions/"affinity_*.json 2>/dev/null || true

    # ---- Boltz-2 执行（传入标准 YAML，不再是 JSON）--------------------------
    boltz predict "$JOB_YAML" \
        --out_dir "$JOB_DIR" \
        --accelerator gpu \
        --devices $N_GPUS \
        --strategy ddp \
        --recycling_steps $RECYCLING_STEPS \
        --diffusion_samples $DIFFUSION_SAMPLES \
        --num_workers $NUM_WORKERS \
        --override \
        >> "$LOG_FILE" 2>&1
    EXIT_CODE=$?

    if [ $EXIT_CODE -eq 0 ]; then
        date '+%Y-%m-%d %H:%M:%S' > "$DONE_MARKER"
        log "  [OK] → $JOB_DIR"
    else
        FAILED=$((FAILED + 1))
        log "  [FAIL] exit=$EXIT_CODE  retry: rm -rf $JOB_DIR && ./run_boltz2.sh"
    fi
done

# ---- 汇总 -------------------------------------------------------------------
FINAL_DONE=$(find "$OUTPUT_DIR" -name ".done" 2>/dev/null | wc -l)
write_progress "$FINAL_DONE" "$TOTAL_JOBS" "COMPLETE"

echo ""
log "============================================================"
log "Run finished: $(date)"
log "  Total      : $TOTAL_JOBS"
log "  Completed  : $FINAL_DONE"
log "  Failed     : $FAILED"
log "  Output     : $OUTPUT_DIR"
log "============================================================"

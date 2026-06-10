#!/bin/bash
# ==============================================================================
# test_autoinhib_em.sh — 带 watchdog 的 EM 冒烟测试 (5 变体并行, --gpus 5)
# ==============================================================================
# 启动:  bash scripts/pipeline/test_autoinhib_em.sh
# 行为:
#   1. 启动 runner (--gpus 5 --phase em) 在 nohup 后台
#   2. 同时启动 monitor 子任务, 每 30s 报告:
#      - 完成的 .done_em 数 / 5
#      - 每个 mdrun 的当前 step / Fmax
#      - 总耗时
#   3. 20 分钟硬上限; 到点 kill runner + 所有 gmx
#   4. 最后打印 5 变体的 EM 状态
# ==============================================================================
set -u

ROOT=/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan
OUT=$ROOT/output/gromacs_md_autoinhib
LOG=$OUT/test_em.log
WALL_LIMIT_SEC=900     # 15 min 硬上限
SAMPLE_INTERVAL=15     # 监控采样间隔 (l-bfgs 应该 1-3 min 就完)
DONE_DIR=$OUT

# --- 准备 env ----------------------------------------------------------------
export GMX=$ROOT/envs/gromacs/bin.AVX2_256/gmx
export GMX_PY=$ROOT/envs/gromacs/bin/python
export GMXLIB=$ROOT/force_fields
export OCL_ICD_VENDORS=$ROOT/opencl_vendors
export PATH=$ROOT/envs/gromacs/bin.AVX2_256:$PATH

mkdir -p "$OUT"
rm -f "$OUT/.done_em_test"   # 测试 marker
> "$LOG"

# --- 启动 runner (后台) -------------------------------------------------------
echo "[$(date '+%H:%M:%S')] 启动 runner (--gpus 5 --phase em)" | tee -a "$LOG"
bash $ROOT/scripts/pipeline/run_gromacs_vwf_md.sh \
    --system autoinhib --phase em --gpus 5 \
    > "$OUT/runner.stdout" 2>&1 &
RUNNER_PID=$!
echo "[$(date '+%H:%M:%S')] runner PID=$RUNNER_PID" | tee -a "$LOG"

# --- watchdog 主循环 ---------------------------------------------------------
START=$(date +%s)
LAST_DONE_COUNT=0
LAST_FMAX_SUM=999999

while true; do
    NOW=$(date +%s)
    ELAPSED=$((NOW - START))
    ELAPSED_MIN=$((ELAPSED / 60))
    ELAPSED_SEC=$((ELAPSED % 60))

    # 1. 检查 .done_em 数
    DONE_COUNT=$(find "$OUT" -maxdepth 3 -name ".done_em" 2>/dev/null | wc -l)
    DONE_FILES=$(find "$OUT" -maxdepth 3 -name ".done_em" 2>/dev/null \
                 | xargs -I {} dirname {} 2>/dev/null \
                 | xargs -I {} basename {} 2>/dev/null \
                 | sort -u | tr '\n' ' ')

    # 2. 检查每个 mdrun 的当前 step / Fmax (从 em.log 末行提取)
    MPROBE=""
    for V in $OUT/VWF_*/; do
        [ ! -d "$V" ] && continue
        VNAME=$(basename "$V")
        ELOG="$V/em/em.log"
        if [ -f "$ELOG" ]; then
            LAST=$(grep -E "^Step=" "$ELOG" 2>/dev/null | tail -1)
            MPROBE="$MPROBE  $VNAME: ${LAST:0:80}"
        fi
    done

    # 3. 检查 runner 是否还活着
    if ! kill -0 $RUNNER_PID 2>/dev/null; then
        echo "[$(date '+%H:%M:%S')] runner 退出" | tee -a "$LOG"
        break
    fi

    # 4. 输出
    echo ""
    echo "[$(date '+%H:%M:%S')] t=${ELAPSED_MIN}m${ELAPSED_SEC}s  done=$DONE_COUNT/5"
    echo "  done_files: ${DONE_FILES:-(none)}"
    if [ -n "$MPROBE" ]; then
        echo "  live mdrun status:"
        echo "$MPROBE" | sed 's/^/    /'
    fi

    # 5. 5 变体全完成 → 提前退出
    if [ "$DONE_COUNT" -ge 5 ]; then
        echo "[$(date '+%H:%M:%S')] 全部 5 变体 EM 完成!" | tee -a "$LOG"
        break
    fi

    # 6. 硬超时
    if [ "$ELAPSED" -gt "$WALL_LIMIT_SEC" ]; then
        echo "[$(date '+%H:%M:%S')] ** 超时 ${WALL_LIMIT_SEC}s, kill 一切" | tee -a "$LOG"
        kill -9 $RUNNER_PID 2>/dev/null
        pkill -9 -f "gmx " 2>/dev/null
        pkill -9 -f "run_gromacs_vwf_md" 2>/dev/null
        sleep 2
        break
    fi

    sleep $SAMPLE_INTERVAL
done

# --- 收尾 --------------------------------------------------------------------
echo ""
echo "============================================================"
echo "[$(date '+%H:%M:%S')] 测试结束"
echo "  总耗时: ${ELAPSED_MIN}m${ELAPSED_SEC}s"
echo "  .done_em 数: $(find $OUT -name '.done_em' 2>/dev/null | wc -l) / 5"
echo "============================================================"
echo ""
echo "每个变体的 EM 状态:"
for V in $OUT/VWF_*/; do
    VNAME=$(basename "$V")
    if [ -f "$V/.done_em" ]; then
        STATUS="✓ DONE at $(cat $V/.done_em)"
        LAST=$(grep -E "^Step=" "$V/em/em.log" 2>/dev/null | tail -1)
        echo "  $VNAME: $STATUS"
        echo "    ${LAST:0:100}"
    elif [ -f "$V/em/em.log" ]; then
        LAST=$(grep -E "^Step=" "$V/em/em.log" 2>/dev/null | tail -1)
        echo "  $VNAME: ✗ TIMEOUT (last: ${LAST:0:100})"
    else
        echo "  $VNAME: ✗ NEVER STARTED"
    fi
done

echo ""
echo "完整 log: $LOG"
echo "worker log: $OUT/worker_*.log"

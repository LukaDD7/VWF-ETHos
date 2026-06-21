#!/bin/bash
# Launch tagged slow-rate SMD go/no-go probes without overwriting default fast SMD reps.
#
# Default: run WT only on GPU0 to avoid oversubscribing CPUs while full mutant prep is active.
# To run all three controls sequentially on one GPU:
#   VARIANTS="WT R1306W R1374H" GPU_ID=0 bash scripts/pipeline/run_7a6o_smd_slow_probe.sh
# To run a different number of reps or rate:
#   NREPS=5 RATE=0.00025 SMD_TAG=slow025 bash scripts/pipeline/run_7a6o_smd_slow_probe.sh
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT_DIR"

VARIANTS="${VARIANTS:-WT}"
GPU_ID="${GPU_ID:-0}"
NREPS="${NREPS:-5}"
export RATE="${RATE:-0.00025}"
export SMD_TAG="${SMD_TAG:-slow025}"
export NTOMP="${NTOMP:-8}"

stamp="$(date +%Y%m%d_%H%M%S)"
log_dir="$ROOT_DIR/output/gromacs_md_autoinhib/_smd_slow_probe_logs_$stamp"
mkdir -p "$log_dir"
echo "[slow-probe] variants=$VARIANTS gpu=$GPU_ID reps=$NREPS RATE=$RATE tag=$SMD_TAG NTOMP=$NTOMP" | tee "$log_dir/summary.log"

for variant in $VARIANTS; do
  echo "[slow-probe] $variant start $(date '+%F %T')" | tee -a "$log_dir/summary.log"
  bash scripts/pipeline/run_7a6o_smd.sh "$variant" "$GPU_ID" "$NREPS" 2>&1 | tee "$log_dir/${variant}.log"
  echo "[slow-probe] $variant done $(date '+%F %T')" | tee -a "$log_dir/summary.log"
done

python3 scripts/pipeline/analyze_7a6o_smd.py \
  --input output/gromacs_md_autoinhib \
  --output output/md_7a6o_smd_${SMD_TAG}_features.csv \
  --tag "$SMD_TAG" 2>&1 | tee "$log_dir/analyze_${SMD_TAG}.log"

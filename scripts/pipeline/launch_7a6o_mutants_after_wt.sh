#!/bin/bash
set -u
ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
MD_ROOT="$ROOT_DIR/output/gromacs_md_autoinhib"
RUNNER="$ROOT_DIR/scripts/pipeline/run_7a6o_autoinhib_md_batch.sh"
VARIANTS="${VARIANTS:-R1306W,R1306Q,R1308C,I1309V,S1310F,W1313C,V1314F,V1316M,P1337L,R1341Q,R1341W,R1374C,R1374H,G1324S}"
GPU_IDS="${GPU_IDS:-0,1,2,3,4,5,6}"
MAX_PARALLEL="${MAX_PARALLEL:-3}"
NS="${NS:-50}"
NVT_PS="${NVT_PS:-50}"
NPT_PS="${NPT_PS:-200}"
NTOMP="${NTOMP:-16}"
LOG="$MD_ROOT/mutants_7a6o_batch_$(date '+%Y%m%d_%H%M%S').log"

[ -x "$RUNNER" ] || { echo "[FATAL] missing runner: $RUNNER"; exit 1; }
WT_PROD_DIR="$MD_ROOT/7A6O_WT/md_7a6o"
WT_DONE_FILE=""
for cand in "$WT_PROD_DIR/md_prod.gro" "$WT_PROD_DIR"/md_prod_run.part*.gro; do
  [ -f "$cand" ] && WT_DONE_FILE="$cand"
done
[ -n "$WT_DONE_FILE" ] || {
  echo "[FATAL] WT production is not complete yet: missing md_prod.gro or md_prod_run.part*.gro in $WT_PROD_DIR"
  echo "        If you intentionally want to start mutants before WT finishes, set REQUIRE_WT_DONE=0."
  [ "${REQUIRE_WT_DONE:-1}" = "0" ] || exit 2
}
if pgrep -af 'gmx mdrun|run_7a6o_autoinhib_md_batch' >/dev/null; then
  echo "[FATAL] Existing MD/runner process detected; refusing to stack another batch."
  pgrep -af 'gmx mdrun|run_7a6o_autoinhib_md_batch'
  exit 3
fi

echo "[INFO] Launching mutant batch"
echo "[INFO] wt_done=${WT_DONE_FILE:-not-required}"
echo "[INFO] variants=$VARIANTS"
echo "[INFO] gpu_ids=$GPU_IDS max_parallel=$MAX_PARALLEL ns=$NS"
echo "[INFO] log=$LOG"
nohup bash "$RUNNER" --phase all --variants "$VARIANTS" --gpu-ids "$GPU_IDS" \
  --max-parallel "$MAX_PARALLEL" --ns "$NS" --nvt-ps "$NVT_PS" --npt-ps "$NPT_PS" --ntomp "$NTOMP" \
  > "$LOG" 2>&1 &
echo "[INFO] PID $!"

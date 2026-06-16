#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
OUT_ROOT="$ROOT_DIR/output/gromacs_md_autoinhib"
RUNNER="$ROOT_DIR/scripts/pipeline/run_7a6o_variant_direct.sh"
CHECK_INTERVAL="${CHECK_INTERVAL:-300}"
STAMP="$(date '+%Y%m%d_%H%M')"

complete_variant() {
  local variant="$1"
  local work="$OUT_ROOT/$variant/md_7a6o"
  [ -f "$work/md_prod.gro" ] || compgen -G "$work/md_prod.part*.gro" >/dev/null
}

runner_alive() {
  local variant="$1"
  local gpu="$2"
  pgrep -f "run_7a6o_variant_direct.sh $variant $gpu" >/dev/null
}

wait_then_launch() {
  local current="$1"
  local next="$2"
  local gpu="$3"
  local log="$OUT_ROOT/${next}_queued_gpu${gpu}_${STAMP}.log"

  echo "[$(date '+%F %T')] gpu=$gpu wait current=$current then launch next=$next"
  while ! complete_variant "$current"; do
    if ! runner_alive "$current" "$gpu"; then
      echo "[$(date '+%F %T')] ERROR: $current runner on gpu=$gpu exited before production completed" >&2
      exit 1
    fi
    sleep "$CHECK_INTERVAL"
  done

  echo "[$(date '+%F %T')] gpu=$gpu $current complete; launching $next"
  setsid bash "$RUNNER" "$next" "$gpu" > "$log" 2>&1 < /dev/null &
  echo "[$(date '+%F %T')] gpu=$gpu launched $next pid=$! log=$log"
}

if [ "$#" -eq 3 ]; then
  wait_then_launch "$1" "$2" "$3"
  exit 0
fi

if [ "$#" -ne 0 ]; then
  echo "usage: $0 [CURRENT_VARIANT NEXT_VARIANT GPU]" >&2
  exit 2
fi

wait_then_launch R1306W R1374H 0 &
wait_then_launch R1306Q V1316M 1 &
wait_then_launch R1308C P1337L 2 &
wait_then_launch I1309V R1341Q 3 &
wait_then_launch S1310F R1341W 4 &
wait_then_launch W1313C R1374C 5 &
wait_then_launch V1314F G1324S 6 &
wait

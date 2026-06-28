#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
PKG_DIR="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/boltz_missing_2m"
GPUS="${GPUS:-4}"
GPU_IDS="${GPU_IDS:-}"
mkdir -p "$PKG_DIR/run_panel/boltz_results"
cmd=(bash "$ROOT_DIR/scripts/pipeline/run_vwd_functional_boltz2_panel.sh"
  --input-dir "$PKG_DIR/run_panel/yamls"
  --out-dir "$PKG_DIR/run_panel/boltz_results"
  --gpus "$GPUS")
if [ -n "$GPU_IDS" ]; then
  cmd+=(--gpu-ids "$GPU_IDS")
fi
echo "${cmd[@]}"
"${cmd[@]}"

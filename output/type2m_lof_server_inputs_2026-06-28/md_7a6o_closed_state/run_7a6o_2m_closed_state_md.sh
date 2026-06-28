#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
PKG_DIR="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/md_7a6o_closed_state"
GPU_IDS="${GPU_IDS:-4,5,6,7}"
NS="${NS:-50}"
PY_BIN="${PY:-python3}"
if [ ! -f "$ROOT_DIR/structures/7A6O_AIM_A1_clean.pdb" ]; then
  "$PY_BIN" "$ROOT_DIR/scripts/pipeline/fetch_clean_7a6o.py"
fi
bash "$ROOT_DIR/scripts/pipeline/run_new_2b_saltbridge_md.sh"   --stage all   --variants-file "$PKG_DIR/7a6o_2m_closed_state_new_variants.txt"   --gpu-ids "$GPU_IDS"   --ns "$NS"

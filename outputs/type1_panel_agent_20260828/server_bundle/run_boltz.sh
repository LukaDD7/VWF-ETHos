#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle"
GPUS="${GPUS:-4}"
bash "$ROOT_DIR/scripts/pipeline/run_vwd_functional_boltz2_panel.sh"   --input-dir "$BUNDLE/boltz/yamls"   --out-dir "$BUNDLE/results/boltz/raw"   --gpus "$GPUS"
python "$ROOT_DIR/scripts/pipeline/parse_vwd_functional_boltz2_results.py"   --results-dir "$BUNDLE/results/boltz/raw"   --manifest "$BUNDLE/boltz/job_manifest.csv"   --output "$BUNDLE/results/boltz/boltz_results_summary.csv"

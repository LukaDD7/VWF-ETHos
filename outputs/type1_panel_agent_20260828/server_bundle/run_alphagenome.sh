#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle"
: "${ALPHAGENOME_API_KEY:?export ALPHAGENOME_API_KEY before running}"
python "$ROOT_DIR/scripts/pipeline/run_type1_10case_alphagenome.py"   --requests "$BUNDLE/input/alphagenome_requests.csv"   --output-dir "$BUNDLE/results/alphagenome" --expected-requests 9

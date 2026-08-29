#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle"
OUT="$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_results"
mkdir -p "$OUT"
python "$ROOT_DIR/scripts/pipeline/ingest_type1_10case_results.py"   --alphagenome-features "$BUNDLE/results/alphagenome/alphagenome_agent_features.csv"   --boltz-summary "$BUNDLE/results/boltz/boltz_results_summary.csv"   --boltz-manifest "$BUNDLE/boltz/job_manifest.csv"   --boltz-variants "$BUNDLE/input/boltz_variants.csv"   --md-features "$BUNDLE/results/md/aim_a1_contacts_summary.csv"   --expected-alpha 7 --expected-boltz-jobs 20 --expected-boltz-wt 4   --output "$OUT/type2b_computational_evidence.csv"
python "$ROOT_DIR/scripts/run_vwd_langgraph_v0.py"   --workbook "$ROOT_DIR/outputs/type2b_panel_agent_20260829/type2b_panel_agent_package.xlsx"   --output-dir "$OUT/agent_run" --computational-panels

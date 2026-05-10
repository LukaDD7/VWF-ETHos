#!/usr/bin/env bash
set -euo pipefail

# Run the pre-generated VWD/VWF functional Boltz-2 panel on a GPU server.
#
# Usage:
#   bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh
#
# Optional environment variables:
#   PANEL_DIR=output/boltz2_vwd_functional_panel
#   OUT_DIR=output/boltz2_vwd_functional_panel/boltz_results
#   DEVICES=1
#   RECYCLING_STEPS=3
#   DIFFUSION_SAMPLES=5
#   NUM_WORKERS=8
#   EXTRA_BOLTZ_ARGS="--use_msa_server"

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PANEL_DIR="${PANEL_DIR:-$ROOT_DIR/output/boltz2_vwd_functional_panel}"
YAML_DIR="$PANEL_DIR/yamls"
OUT_DIR="${OUT_DIR:-$PANEL_DIR/boltz_results}"
DEVICES="${DEVICES:-1}"
RECYCLING_STEPS="${RECYCLING_STEPS:-3}"
DIFFUSION_SAMPLES="${DIFFUSION_SAMPLES:-5}"
NUM_WORKERS="${NUM_WORKERS:-8}"
EXTRA_BOLTZ_ARGS="${EXTRA_BOLTZ_ARGS:-}"

if ! command -v boltz >/dev/null 2>&1; then
  echo "ERROR: boltz command not found. Activate/install the Boltz-2 environment first." >&2
  exit 1
fi

if [[ ! -d "$YAML_DIR" ]]; then
  echo "ERROR: YAML directory not found: $YAML_DIR" >&2
  echo "Regenerate with: python scripts/pipeline/generate_vwd_functional_boltz2_yamls.py --write-json-batches" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

echo "Running VWD/VWF functional Boltz-2 panel"
echo "  YAML_DIR          : $YAML_DIR"
echo "  OUT_DIR           : $OUT_DIR"
echo "  DEVICES           : $DEVICES"
echo "  RECYCLING_STEPS   : $RECYCLING_STEPS"
echo "  DIFFUSION_SAMPLES : $DIFFUSION_SAMPLES"
echo "  NUM_WORKERS       : $NUM_WORKERS"
echo

# shellcheck disable=SC2086
boltz predict "$YAML_DIR" \
  --out_dir "$OUT_DIR" \
  --accelerator gpu \
  --devices "$DEVICES" \
  --recycling_steps "$RECYCLING_STEPS" \
  --diffusion_samples "$DIFFUSION_SAMPLES" \
  --num_workers "$NUM_WORKERS" \
  $EXTRA_BOLTZ_ARGS

#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle"
: "${FOLDX:?export FOLDX=/absolute/path/to/foldx}"
python "$ROOT_DIR/scripts/pipeline/build_2b_mutants_foldx.py"   --wt "$ROOT_DIR/structures/7A6O_AIM_A1_clean.pdb"   --foldx "$FOLDX" --detect-offset --variants P1413L   --out-dir "$ROOT_DIR/structures/7a6o_mutants"
cp "$ROOT_DIR/structures/7a6o_mutants/WT_Repair.pdb"   "$ROOT_DIR/structures/7a6o_mutants/TYPE1_WT_MATCHED.pdb"
bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh"   --pdb "$ROOT_DIR/structures/7a6o_mutants/TYPE1_WT_MATCHED.pdb" --variant TYPE1_WT_MATCHED --skip-vacuum
bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh"   --pdb "$ROOT_DIR/structures/7a6o_mutants/P1413L.pdb" --variant P1413L --skip-vacuum
NS=50 bash "$ROOT_DIR/scripts/pipeline/run_7a6o_variant_direct.sh" TYPE1_WT_MATCHED "${WT_GPU:-0}" &
WT_PID=$!
NS=50 bash "$ROOT_DIR/scripts/pipeline/run_7a6o_variant_direct.sh" P1413L "${MUT_GPU:-1}" &
MUT_PID=$!
wait "$WT_PID" "$MUT_PID"
python "$ROOT_DIR/scripts/pipeline/analyze_7a6o_completed_md.py"   --input "$ROOT_DIR/output/gromacs_md_autoinhib"   --output "$BUNDLE/results/md"   --variants TYPE1_WT_MATCHED,P1413L --wt-variant TYPE1_WT_MATCHED

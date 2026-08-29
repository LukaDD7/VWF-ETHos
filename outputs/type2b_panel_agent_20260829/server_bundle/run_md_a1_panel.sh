#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle"
: "${FOLDX:?export FOLDX=/absolute/path/to/foldx}"
VARIANTS=(R1308C S1310F V1316M R1341W A1461D)
python "$ROOT_DIR/scripts/pipeline/build_2b_mutants_foldx.py"   --wt "$ROOT_DIR/structures/7A6O_AIM_A1_clean.pdb" --foldx "$FOLDX" --detect-offset   --variants "$(IFS=,; echo "${VARIANTS[*]}")" --out-dir "$ROOT_DIR/structures/7a6o_mutants"
cp "$ROOT_DIR/structures/7a6o_mutants/WT_Repair.pdb" "$ROOT_DIR/structures/7a6o_mutants/TYPE2B_WT_MATCHED.pdb"
for variant in TYPE2B_WT_MATCHED "${VARIANTS[@]}"; do
  bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh"     --pdb "$ROOT_DIR/structures/7a6o_mutants/${variant}.pdb" --variant "$variant" --skip-vacuum
done
# Six independent 50-ns jobs; GPU IDs may be changed to match the server.
jobs=(TYPE2B_WT_MATCHED R1308C S1310F V1316M R1341W A1461D)
gpu_ids=(${GPU_IDS:-0 1 2 3 4 5})
pids=()
for i in "${!jobs[@]}"; do
  NS=50 bash "$ROOT_DIR/scripts/pipeline/run_7a6o_variant_direct.sh" "${jobs[$i]}" "${gpu_ids[$i]}" &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
python "$ROOT_DIR/scripts/pipeline/analyze_7a6o_completed_md.py"   --input "$ROOT_DIR/output/gromacs_md_autoinhib" --output "$BUNDLE/results/md"   --variants "TYPE2B_WT_MATCHED,$(IFS=,; echo "${VARIANTS[*]}")" --wt-variant TYPE2B_WT_MATCHED

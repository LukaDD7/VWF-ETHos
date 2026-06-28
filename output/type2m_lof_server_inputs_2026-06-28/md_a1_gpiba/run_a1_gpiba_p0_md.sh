#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
PKG_DIR="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba"
VARIANTS_FILE="${VARIANTS_FILE:-$PKG_DIR/a1_gpiba_p0_plus_anchor_variants.txt}"
GPU_IDS="${GPU_IDS:-4,5,6,7}"
NS="${NS:-50}"
RELAX_JOBS="${RELAX_JOBS:-6}"
RELAX_NTOMP="${RELAX_NTOMP:-8}"
NTOMP="${NTOMP:-8}"
FOLDX_BIN="${FOLDX:-foldx}"
GMX_BIN="${GMX:-$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx}"
PY_BIN="${PY:-python3}"
MUT_DIR="$ROOT_DIR/structures/1sq0_a1_gpiba_mutants"
SCHEDULER="$ROOT_DIR/scripts/pipeline/run_md_resilient.py"
RUNNER="$ROOT_DIR/scripts/pipeline/run_md_variant_direct.sh"
[ -f "$SCHEDULER" ] || SCHEDULER="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/support_scripts/run_md_resilient.py"
[ -f "$RUNNER" ] || RUNNER="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/support_scripts/run_md_variant_direct.sh"

"$PY_BIN" "$ROOT_DIR/scripts/pipeline/build_2b_mutants_foldx.py"   --wt "$ROOT_DIR/structures/1SQ0.pdb"   --chain A   --offset 763   --variants-file "$VARIANTS_FILE"   --out-dir "$MUT_DIR"   --foldx "$FOLDX_BIN"

mkdir -p "$ROOT_DIR/output/gromacs_md_a1_gpiba"
active=0
while read -r v; do
  v="${v%%#*}"; v="$(echo "$v" | awk '{print $1}')"
  [ -z "$v" ] && continue
  [ -f "$ROOT_DIR/output/gromacs_md_a1_gpiba/$v/relax_pdb/solv_ions_em.gro" ] && continue
  RELAX_NTOMP="$RELAX_NTOMP" GMX="$GMX_BIN"     bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh"       --pdb "$MUT_DIR/$v.pdb" --variant "$v" --system a1_gpiba --skip-vacuum       > "$ROOT_DIR/output/gromacs_md_a1_gpiba/${v}_relax.log" 2>&1 &
  active=$((active + 1))
  if [ "$active" -ge "$RELAX_JOBS" ]; then
    wait -n
    active=$((active - 1))
  fi
done < "$VARIANTS_FILE"
wait

NS="$NS" NTOMP="$NTOMP" "$PY_BIN" "$SCHEDULER"   --root "$ROOT_DIR"   --system a1_gpiba   --md-tag md_a1_gpiba   --variants-file "$VARIANTS_FILE"   --mutant-dir "$MUT_DIR"   --runner "$RUNNER"   --gpu-ids "$GPU_IDS"   --ns "$NS"   --ntomp "$NTOMP"

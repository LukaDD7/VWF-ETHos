#!/bin/bash
# Unified P0 experimental-structure MD launcher for VWF.
#
# Runs, in order:
#   1. A2 folded stability       : 3GXB   + A1500V
#   2. D'D3/FVIII-binding axis   : 6N29   + R1205H
#   3. A1-GPIb-alpha 2B axis     : 1SQ0   + R1308C, S1310F, A1461D, R1341W
#
# Default: 20 ns pilot on GPU 0, 1, 2, and 3.
#
# Usage:
#   nohup bash scripts/pipeline/run_vwf_p0_experimental_md.sh \
#     --axis all --gpu-ids 0,1,2,3 --ns 20 \
#     > output/p0_experimental_md.log 2>&1 &

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
LZY_ROOT="$(cd "$ROOT_DIR/../.." && pwd)"

AXIS="all"
GPU_IDS="0,1,2,3"
NS="20"
NTOMP="8"
RELAX_JOBS="2"
RELAX_NTOMP="8"
MIN_FREE_MIB="120000"

PY="${PY:-$LZY_ROOT/envs/gromacs/bin/python}"
GMX="${GMX:-$LZY_ROOT/envs/gromacs-cuda/bin.AVX2_256/gmx}"
FOLDX="${FOLDX:-$LZY_ROOT/tools/foldx/foldx_20270131}"

OUT_ROOT="$ROOT_DIR/output/p0_md_experimental"
MDP_ROOT="$ROOT_DIR/scripts/pipeline"

usage() {
  cat <<'EOF'
Usage:
  bash scripts/pipeline/run_vwf_p0_experimental_md.sh [options]

Options:
  --axis AXIS         all | a2 | dp_d3 | a1_gpiba      (default: all)
  --gpu-ids IDS       Comma-separated GPU IDs          (default: 0,1,2,3)
  --ns NS             Production length in ns          (default: 20)
  --ntomp N           OpenMP threads per MD job        (default: 8)
  --relax-jobs N      Parallel relaxation jobs          (default: 2)
  --relax-ntomp N     OpenMP threads per relaxation     (default: 8)
  --min-free-mib N    Minimum free GPU memory in MiB    (default: 120000)
  -h, --help          Show this help

Examples:
  # Run all three P0 axes, 20 ns, GPU 0, 1, 2, and 3
  bash scripts/pipeline/run_vwf_p0_experimental_md.sh --axis all --gpu-ids 0,1,2,3 --ns 20

  # Run only the A1-GPIb-alpha axis
  bash scripts/pipeline/run_vwf_p0_experimental_md.sh --axis a1_gpiba --gpu-ids 0,1,2,3 --ns 20
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --axis) AXIS="$2"; shift 2 ;;
    --gpu-ids) GPU_IDS="$2"; shift 2 ;;
    --ns) NS="$2"; shift 2 ;;
    --ntomp) NTOMP="$2"; shift 2 ;;
    --relax-jobs) RELAX_JOBS="$2"; shift 2 ;;
    --relax-ntomp) RELAX_NTOMP="$2"; shift 2 ;;
    --min-free-mib) MIN_FREE_MIB="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[WARN] Unknown argument: $1"; shift ;;
  esac
done

mkdir -p "$OUT_ROOT"

echo "============================================================"
echo "VWF P0 experimental-structure MD"
echo "Started : $(date)"
echo "ROOT    : $ROOT_DIR"
echo "Axis    : $AXIS"
echo "GPUs    : $GPU_IDS"
echo "NS      : $NS"
echo "NTOMP   : $NTOMP"
echo "PY      : $PY"
echo "GMX     : $GMX"
echo "FOLDX   : $FOLDX"
echo "Output  : $OUT_ROOT"
echo "============================================================"

# ------------------------------------------------------------------------------
# Helper: build, relax, and run MD for one experimental-structure axis.
# ------------------------------------------------------------------------------
run_axis() {
  local axis="$1"
  local structure="$2"
  local chain="$3"
  local offset="$4"
  local system="$5"
  shift 5
  local variants=("$@")

  local axis_dir="$OUT_ROOT/$axis"
  local mutant_dir="$ROOT_DIR/structures/${axis}_p0_mutants"
  local build_variants_file="$axis_dir/build_variants.txt"
  local md_variants_file="$axis_dir/md_variants.txt"
  local out_md_root="$ROOT_DIR/output/gromacs_md_$system"

  mkdir -p "$axis_dir" "$mutant_dir" "$out_md_root"

  echo ""
  echo "============================================================"
  echo "[$(date '+%F %T')] AXIS: $axis"
  echo "Structure : $structure"
  echo "Chain     : $chain"
  echo "Offset    : $offset"
  echo "System    : $system"
  echo "Variants  : ${variants[*]}"
  echo "============================================================"

  # 1) Build mutant PDBs on the experimental WT structure.
  printf "%s\n" "${variants[@]}" > "$build_variants_file"

  "$PY" "$ROOT_DIR/scripts/pipeline/build_2b_mutants_foldx.py" \
    --wt "$ROOT_DIR/$structure" \
    --foldx "$FOLDX" \
    --chain "$chain" \
    --offset "$offset" \
    --variants-file "$build_variants_file" \
    --out-dir "$mutant_dir"

  # Use the repaired WT as the matched WT baseline.
  if [[ -f "$mutant_dir/WT_Repair.pdb" ]]; then
    cp -f "$mutant_dir/WT_Repair.pdb" "$mutant_dir/WT.pdb"
  else
    echo "[FATAL] WT_Repair.pdb not found for axis $axis" >&2
    return 1
  fi

  # 2) Relax WT and mutants.
  printf "WT\n" > "$md_variants_file"
  printf "%s\n" "${variants[@]}" >> "$md_variants_file"

  local active=0
  while IFS= read -r v; do
    [[ -z "$v" ]] && continue

    if [[ -f "$out_md_root/$v/relax_pdb/solv_ions_em.gro" ]]; then
      echo "[$(date '+%F %T')] Relax already done: $v"
      continue
    fi

    echo "[$(date '+%F %T')] Relax start: $v"
    GMX="$GMX" RELAX_NTOMP="$RELAX_NTOMP" \
      bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh" \
        --pdb "$mutant_dir/$v.pdb" \
        --variant "$v" \
        --system "$system" \
        --skip-vacuum \
        > "$out_md_root/${v}_relax.log" 2>&1 &

    active=$((active + 1))
    if (( active >= RELAX_JOBS )); then
      wait -n
      active=$((active - 1))
    fi
  done < "$md_variants_file"
  wait

  # 3) Run resumable multi-GPU MD.
  NS="$NS" NTOMP="$NTOMP" \
  "$PY" "$ROOT_DIR/scripts/pipeline/run_md_resilient.py" \
    --root "$ROOT_DIR" \
    --system "$system" \
    --md-tag "md_${NS}ns_p0" \
    --variants-file "$md_variants_file" \
    --mutant-dir "$mutant_dir" \
    --runner "$ROOT_DIR/scripts/pipeline/run_md_variant_direct.sh" \
    --gpu-ids "$GPU_IDS" \
    --ns "$NS" \
    --ntomp "$NTOMP" \
    --min-free-mib "$MIN_FREE_MIB"

  echo "[$(date '+%F %T')] AXIS DONE: $axis"
}

# ------------------------------------------------------------------------------
# Axis definitions.
# ------------------------------------------------------------------------------
run_a2() {
  run_axis \
    "a2" \
    "structures/3GXB.pdb" \
    "A" \
    "0" \
    "a2_folded" \
    "A1500V"
}

run_dp_d3() {
  "$PY" "$ROOT_DIR/scripts/pipeline/clean_pdb_for_md.py" \
    --pdb "$ROOT_DIR/structures/6N29.pdb" \
    --chain A \
    --out "$ROOT_DIR/structures/6N29_A_clean.pdb"

  run_axis \
    "dp_d3" \
    "structures/6N29_A_clean.pdb" \
    "A" \
    "0" \
    "dprime_d3" \
    "R1205H"
}

run_a1_gpiba() {
  run_axis \
    "a1_gpiba" \
    "structures/1SQ0.pdb" \
    "A" \
    "763" \
    "a1_gpiba" \
    "R1308C" "S1310F" "A1461D" "R1341W"
}

case "$AXIS" in
  all)
    run_a2
    run_dp_d3
    run_a1_gpiba
    ;;
  a2)
    run_a2
    ;;
  dp_d3)
    run_dp_d3
    ;;
  a1_gpiba)
    run_a1_gpiba
    ;;
  *)
    echo "[FATAL] Unknown axis: $AXIS" >&2
    usage
    exit 1
    ;;
esac

echo ""
echo "============================================================"
echo "VWF P0 experimental-structure MD finished: $(date)"
echo "Output: $OUT_ROOT"
echo "============================================================"

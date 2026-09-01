#!/bin/bash
# Retry R1205H D'D3 MD after the CUDA-700 failure.
#
# This route is deliberately more conservative than the original P0 run:
#   1. Archive the failed skip-vacuum relaxation and MD attempt.
#   2. Re-run full vacuum + solvated relaxation.
#   3. Run a 20,000-step unconstrained refinement EM on CPU.
#   4. Require Fmax to fall below a configurable threshold.
#   5. Run a fresh 20 ns MD under GMX_GPU_MODE (default: safe).
#
# Usage:
#   nohup bash scripts/pipeline/run_vwf_p0_r1205h_retry.sh \
#     --gpu 0 --gpu-mode safe --ns 20 \
#     > output/p0_r1205h_retry.log 2>&1 &

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
LZY_ROOT="$(cd "$ROOT_DIR/../.." && pwd)"

GPU=0
GPU_MODE=safe
NS=20
NTOMP=8
RELAX_NTOMP=8
MIN_FREE_MB=3000
MAX_FMAX=10000

PY="${PY:-$LZY_ROOT/envs/gromacs/bin/python}"
GMX="${GMX:-$LZY_ROOT/envs/gromacs-cuda/bin.AVX2_256/gmx}"

usage() {
  cat <<'EOF'
Usage:
  bash scripts/pipeline/run_vwf_p0_r1205h_retry.sh [options]

Options:
  --gpu ID          GPU ID used by the scheduler             (default: 0)
  --gpu-mode MODE   full | safe | cpu                         (default: safe)
  --ns NS           Production length                         (default: 20)
  --ntomp N         OpenMP threads                            (default: 8)
  --min-free-mb N   Minimum free NFS space in MB              (default: 3000)
  --max-fmax N      Maximum accepted refinement Fmax          (default: 10000)
  -h, --help        Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gpu) GPU="$2"; shift 2 ;;
    --gpu-mode) GPU_MODE="$2"; shift 2 ;;
    --ns) NS="$2"; shift 2 ;;
    --ntomp) NTOMP="$2"; shift 2 ;;
    --min-free-mb) MIN_FREE_MB="$2"; shift 2 ;;
    --max-fmax) MAX_FMAX="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[WARN] unknown argument: $1"; shift ;;
  esac
done

case "$GPU_MODE" in
  full|safe|cpu) ;;
  *) echo "[FATAL] --gpu-mode must be full, safe, or cpu" >&2; exit 2 ;;
esac

variant=R1205H
system=dprime_d3
md_tag="md_${NS}ns_p0_retry_${GPU_MODE}"
out_root="$ROOT_DIR/output/gromacs_md_${system}"
variant_root="$out_root/$variant"
old_relax="$variant_root/relax_pdb"
old_md="$variant_root/md_20ns_p0"
mutant_pdb="$ROOT_DIR/structures/dp_d3_p0_mutants/${variant}.pdb"

free_mb=$(df -m "$ROOT_DIR" | awk 'NR==2 {print $4}')
if (( free_mb < MIN_FREE_MB )); then
  echo "[FATAL] Only ${free_mb} MB free; need at least ${MIN_FREE_MB} MB for relaxation and retry MD." >&2
  exit 1
fi

archive_path() {
  local src="$1"
  local dst="$2"
  if [[ -e "$src" ]]; then
    if [[ -e "$dst" ]]; then
      dst="${dst}_$(date +%Y%m%d_%H%M%S)"
    fi
    mv "$src" "$dst"
    echo "[ARCHIVE] $src -> $dst"
  fi
}

[[ -f "$mutant_pdb" ]] || { echo "[FATAL] missing $mutant_pdb" >&2; exit 2; }

mkdir -p "$out_root"
archive_path "$old_md" "$variant_root/md_20ns_p0_cuda700_v1"
archive_path "$old_relax" "$variant_root/relax_pdb_skipvacuum_v1"

echo "============================================================"
echo "R1205H conservative retry"
echo "Started : $(date)"
echo "GPU     : $GPU"
echo "Mode    : $GPU_MODE"
echo "NS      : $NS"
echo "Free MB : $free_mb"
echo "Output  : $variant_root/$md_tag"
echo "============================================================"

GMX="$GMX" RELAX_NTOMP="$RELAX_NTOMP" \
  bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh" \
    --pdb "$mutant_pdb" \
    --variant "$variant" \
    --system "$system"

relax="$variant_root/relax_pdb"
cd "$relax"
cat > em_refine_retry.mdp <<'EOF'
integrator    = steep
emtol         = 1000.0
emstep        = 0.001
nsteps        = 20000
nstlist       = 20
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.2
rvdw          = 1.2
pbc           = xyz
constraints   = none
EOF

export GMXLIB="$ROOT_DIR/force_fields"
"$GMX" grompp \
  -f em_refine_retry.mdp \
  -c solv_ions_em.gro \
  -p topol.top \
  -o em_refine_retry.tpr \
  -maxwarn 5 > grompp_em_refine_retry.log 2>&1
"$GMX" mdrun \
  -deffnm em_refine_retry \
  -nb cpu \
  -ntmpi 1 \
  -ntomp "$NTOMP" > md_em_refine_retry.log 2>&1
cp -f em_refine_retry.gro solv_ions_em_refined.gro

fmax=$(grep 'Maximum force' md_em_refine_retry.log | tail -1 | grep -oE '[0-9.eE+]+')
echo "Refinement Fmax = ${fmax} kJ/mol/nm; threshold = ${MAX_FMAX}"
if awk -v f="$fmax" -v t="$MAX_FMAX" 'BEGIN { exit !(f+0 > t+0) }'; then
  echo "[FATAL] Refinement did not reduce Fmax below ${MAX_FMAX}; inspect $relax/md_em_refine_retry.log" >&2
  exit 4
fi

variants_file="$ROOT_DIR/output/p0_md_experimental/r1205h_retry_variants.txt"
mkdir -p "$(dirname "$variants_file")"
printf '%s\n' "$variant" > "$variants_file"

export GMX_GPU_MODE="$GPU_MODE"
export NS="$NS"
export NTOMP="$NTOMP"
"$PY" "$ROOT_DIR/scripts/pipeline/run_md_resilient.py" \
  --root "$ROOT_DIR" \
  --system "$system" \
  --md-tag "$md_tag" \
  --variants-file "$variants_file" \
  --mutant-dir "$ROOT_DIR/structures/dp_d3_p0_mutants" \
  --runner "$ROOT_DIR/scripts/pipeline/run_md_variant_direct.sh" \
  --gpu-ids "$GPU" \
  --ns "$NS" \
  --ntomp "$NTOMP"

echo "============================================================"
echo "R1205H retry finished: $(date)"
echo "Output: $variant_root/$md_tag"
echo "============================================================"

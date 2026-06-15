#!/bin/bash
set -u
ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
MD_ROOT="$ROOT_DIR/output/gromacs_md_autoinhib"
GMX="${GMX:-$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx}"
VARIANTS="${VARIANTS:-R1306W,R1306Q,R1308C,I1309V,S1310F,W1313C,V1314F,V1316M,P1337L,R1341Q,R1341W,R1374C,R1374H,G1324S}"
MAX_PARALLEL="${MAX_PARALLEL:-7}"
NTOMP="${NTOMP:-16}"
FORCE="${FORCE:-0}"

[ -x "$GMX" ] || { echo "[FATAL] gmx not found: $GMX"; exit 1; }
export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"

write_mdp() {
  cat > em_refine.mdp <<'EOF'
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
}

refine_one() {
  local variant="$1" work="$MD_ROOT/$variant/relax_pdb"
  [ -d "$work" ] || { echo "[FATAL] missing relax dir: $work"; return 2; }
  [ -f "$work/solv_ions_em.gro" ] || { echo "[FATAL] missing $work/solv_ions_em.gro"; return 2; }
  [ -f "$work/topol.top" ] || { echo "[FATAL] missing $work/topol.top"; return 2; }
  if [ "$FORCE" != "1" ] && [ -f "$work/solv_ions_em_refined.gro" ]; then
    echo "[$(date '+%F %T')] $variant skip existing solv_ions_em_refined.gro"
    return 0
  fi
  (
    cd "$work" || exit 2
    echo "[$(date '+%F %T')] $variant refine EM start"
    write_mdp
    "$GMX" grompp -f em_refine.mdp -c solv_ions_em.gro -p topol.top -o em_refine.tpr -maxwarn 5 > grompp_refine.log 2>&1 || exit 3
    "$GMX" mdrun -deffnm em_refine -nb cpu -ntmpi 1 -ntomp "$NTOMP" -pin on > md_refine.log 2>&1 || exit 4
    cp -f em_refine.gro solv_ions_em_refined.gro
    grep 'Maximum force' md_refine.log | tail -1
    echo "[$(date '+%F %T')] $variant refine EM done"
  )
}

IFS=',' read -r -a VARIANT_ARR <<< "$VARIANTS"
failures=0
pids=()
wait_batch() {
    [ "${#pids[@]}" -eq 0 ] && return 0
  local pid status
  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
      status=$?
      echo "[ERROR] refine child PID $pid failed with status $status" >&2
      failures=$((failures + 1))
    fi
  done
  pids=()
}

for variant in "${VARIANT_ARR[@]}"; do
  refine_one "$variant" > "$MD_ROOT/${variant}_em_refine.log" 2>&1 &
  pids+=("$!")
  if [ "${#pids[@]}" -ge "$MAX_PARALLEL" ]; then
    wait_batch
  fi
done
wait_batch
exit "$failures"

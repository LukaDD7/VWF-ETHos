#!/bin/bash
# Constant-velocity steered MD (SMD) of AIM unfolding for a prepped 7A6O variant.
#
# Pulls the two AIM-flanking anchors apart along z and records pull force/coord.
# The peak (rupture) force / pulling work discriminate 2B (gain of function:
# AIM unfolds at LOWER force) from 2M / WT. Run prep_7a6o_smd.sh first.
#
# Replicates use independent initial velocities (gen_vel + distinct seed):
# rupture force is stochastic, so report mean +/- sd over >= 3 pulls.
#
# Usage:  bash scripts/pipeline/run_7a6o_smd.sh <variant> <gpu_id> [nreps]
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
GMX="${GMX:-$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx}"
[ -x "$GMX" ] || GMX="$ROOT_DIR/envs/gromacs/bin/gmx"
export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"
export OCL_ICD_VENDORS="${OCL_ICD_VENDORS:-$ROOT_DIR/opencl_vendors}"

variant="${1:?variant required}"
gpu="${2:?gpu id required}"
nreps="${3:-3}"
# SMD parameters
PULL_NM="${PULL_NM:-5.0}"      # nm of extension (captures AIM rupture peak)
RATE="${RATE:-0.001}"         # nm/ps  (0.001 = 1 nm/ns; lower = gentler, costlier)
KSPRING="${KSPRING:-1000}"    # kJ/mol/nm^2
NTOMP="${NTOMP:-8}"
case "$gpu" in
  0) po=0;; 1) po=16;; 2) po=32;; 3) po=64;; 4) po=80;; 5) po=96;; 6) po=112;; *) po=$((gpu*NTOMP));;
esac

WORK="$ROOT_DIR/output/gromacs_md_autoinhib/$variant/smd"
cd "$WORK"
for f in smd_npt.gro topol_smd.top index.ndx; do
  [ -f "$f" ] || { echo "FATAL: missing $f -- run prep_7a6o_smd.sh $variant first"; exit 2; }
done
# nsteps = extension / rate / dt(0.002)
nsteps=$(awk -v p="$PULL_NM" -v r="$RATE" 'BEGIN{printf "%d", p/r/0.002}')
echo "[smd] $variant: $nreps reps, pull $PULL_NM nm @ $RATE nm/ps -> $nsteps steps ($(awk -v n=$nsteps 'BEGIN{printf "%.1f",n*0.002/1000}') ns/rep)"

for rep in $(seq 1 "$nreps"); do
  deffnm="smd_rep${rep}"
  if [ -f "${deffnm}.gro" ]; then echo "[smd] $variant rep$rep already done"; continue; fi
  seed=$((12345 + rep * 7919))
  cat > ${deffnm}.mdp <<EOF
integrator          = md
dt                  = 0.002
nsteps              = $nsteps
nstxout-compressed  = 25000
nstenergy           = 5000
nstlog              = 5000
cutoff-scheme       = Verlet
coulombtype         = PME
rcoulomb            = 1.2
rvdw                = 1.2
pbc                 = xyz
nstlist             = 100
tcoupl              = V-rescale
tc-grps             = Protein Non-Protein
tau_t               = 0.1 0.1
ref_t               = 310 310
pcoupl              = no
constraints         = h-bonds
constraint_algorithm = lincs
gen_vel             = yes
gen_temp            = 310
gen_seed            = $seed
continuation        = no
; ---------------- pull / SMD ----------------
pull                       = yes
pull-ngroups               = 2
pull-ncoords               = 1
pull-group1-name           = N_anchor
pull-group2-name           = C_anchor
pull-coord1-type           = umbrella
pull-coord1-geometry       = direction-periodic
pull-coord1-groups         = 1 2
pull-coord1-dim            = N N Y
pull-coord1-vec            = 0 0 1
pull-coord1-start          = yes
pull-coord1-init           = 0.0
pull-coord1-rate           = $RATE
pull-coord1-k              = $KSPRING
pull-nstxout               = 500
pull-nstfout               = 500
pull-pbc-ref-prev-step-com = yes
EOF
  "$GMX" grompp -f ${deffnm}.mdp -c smd_npt.gro -p topol_smd.top -n index.ndx -o ${deffnm}.tpr -maxwarn 5 \
      > g_${deffnm}.log 2>&1
  echo "[smd] $variant rep$rep gpu=$gpu start $(date '+%F %T')"
  CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm ${deffnm} -ntmpi 1 -ntomp "$NTOMP" -gpu_id 0 \
      -nb gpu -pme gpu -pin on -pinstride 2 -pinoffset "$po" \
      -px ${deffnm}_pullx.xvg -pf ${deffnm}_pullf.xvg > md_${deffnm}.stdout 2>&1
  echo "[smd] $variant rep$rep done $(date '+%F %T') -> ${deffnm}_pullf.xvg"
done
echo "[smd] DONE $variant. Analyse: python3 scripts/pipeline/analyze_7a6o_smd.py"

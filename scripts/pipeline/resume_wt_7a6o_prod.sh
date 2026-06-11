#!/bin/bash
set -u
ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="$ROOT_DIR/output/gromacs_md_autoinhib/7A6O_WT/md_7a6o"
GMX="${GMX:-$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx}"
GPU_ID="${GPU_ID:-6}"
NTOMP="${NTOMP:-16}"
NS="${NS:-50}"

[ -x "$GMX" ] || { echo "[FATAL] gmx not executable: $GMX"; exit 1; }
[ -d "$WORK" ] || { echo "[FATAL] missing work dir: $WORK"; exit 1; }
cd "$WORK"

export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"
export OCL_ICD_VENDORS="${OCL_ICD_VENDORS:-$ROOT_DIR/opencl_vendors}"
export GMX_CUDA_GRAPH=1

if pgrep -af "gmx mdrun .*md_prod" >/dev/null; then
  echo "[FATAL] md_prod already appears to be running; not starting another copy."
  pgrep -af "gmx mdrun .*md_prod"
  exit 2
fi

if [ -f md_prod.cpt ]; then
  echo "[INFO] Resuming existing production checkpoint: $WORK/md_prod.cpt"
  exec env CUDA_VISIBLE_DEVICES="$GPU_ID" "$GMX" mdrun -deffnm md_prod -cpi md_prod.cpt -append \
    -ntmpi 1 -ntomp "$NTOMP" -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu -pin on
fi

[ -f npt.gro ] || { echo "[FATAL] missing npt.gro"; exit 3; }
[ -f npt.cpt ] || { echo "[FATAL] missing npt.cpt"; exit 3; }
[ -f topol.top ] || { echo "[FATAL] missing topol.top"; exit 3; }

steps=$(awk -v ns="$NS" 'BEGIN{printf "%d", ns * 500000}')
cat > production.mdp <<MDP
integrator          = md
nsteps              = $steps
nstxout-compressed  = 50000
nstenergy           = 50000
nstlog              = 50000
dt                  = 0.002
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
constraints         = h-bonds
constraint_algorithm = lincs
pcoupl              = C-rescale
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = 1.0
compressibility     = 4.5e-5
gen_vel             = no
continuation        = yes
MDP

"$GMX" grompp -f production.mdp -c npt.gro -t npt.cpt -p topol.top -o md_prod.tpr -maxwarn 5 > g_prod.log 2>&1 || {
  echo "[FATAL] grompp failed; see $WORK/g_prod.log"
  tail -40 g_prod.log
  exit 4
}

echo "[INFO] Starting WT production on GPU $GPU_ID"
exec env CUDA_VISIBLE_DEVICES="$GPU_ID" "$GMX" mdrun -deffnm md_prod \
  -ntmpi 1 -ntomp "$NTOMP" -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu -pin on

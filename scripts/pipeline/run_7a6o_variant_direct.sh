#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
OUT_ROOT="$ROOT_DIR/output/gromacs_md_autoinhib"
GMX="$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx"
variant="${1:?variant required}"
gpu="${2:?gpu id required}"
ns="${NS:-50}"
nvt_ps="${NVT_PS:-50}"
npt_ps="${NPT_PS:-200}"
ntomp="${NTOMP:-16}"
pinoffset="$((gpu * ntomp))"
export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"
export OCL_ICD_VENDORS="${OCL_ICD_VENDORS:-$ROOT_DIR/opencl_vendors}"
unset GMX_CUDA_GRAPH 2>/dev/null || true

steps_from_ps() { awk -v ps="$1" 'BEGIN{printf "%d", ps * 500}'; }
steps_from_ns() { awk -v ns="$1" 'BEGIN{printf "%d", ns * 500000}'; }
mdp_common() {
cat <<'EOF'
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
EOF
}
write_nvt() { local steps="$1"; { echo 'define              = -DPOSRES'; echo 'integrator          = md'; echo "nsteps              = $steps"; echo 'nstxout-compressed  = 5000'; echo 'nstenergy           = 5000'; echo 'nstlog              = 5000'; mdp_common; echo 'pcoupl              = no'; echo 'gen_vel             = yes'; echo 'gen_temp            = 310'; echo 'gen_seed            = -1'; echo 'continuation        = no'; } > nvt.mdp; }
write_npt() { local steps="$1"; { echo 'define              = -DPOSRES'; echo 'integrator          = md'; echo "nsteps              = $steps"; echo 'nstxout-compressed  = 5000'; echo 'nstenergy           = 5000'; echo 'nstlog              = 5000'; mdp_common; echo 'pcoupl              = C-rescale'; echo 'pcoupltype          = isotropic'; echo 'tau_p               = 2.0'; echo 'ref_p               = 1.0'; echo 'compressibility     = 4.5e-5'; echo 'refcoord_scaling    = com'; echo 'gen_vel             = no'; echo 'continuation        = yes'; } > npt.mdp; }
write_prod() { local steps="$1"; { echo 'integrator          = md'; echo "nsteps              = $steps"; echo 'nstxout-compressed  = 50000'; echo 'nstenergy           = 50000'; echo 'nstlog              = 50000'; mdp_common; echo 'pcoupl              = C-rescale'; echo 'pcoupltype          = isotropic'; echo 'tau_p               = 2.0'; echo 'ref_p               = 1.0'; echo 'compressibility     = 4.5e-5'; echo 'gen_vel             = no'; echo 'continuation        = yes'; } > production.mdp; }

relax="$OUT_ROOT/$variant/relax_pdb"
start="$relax/solv_ions_em_refined.gro"
[ -f "$start" ] || start="$relax/solv_ions_em.gro"
[ -f "$start" ] || { echo "missing start for $variant"; exit 2; }
[ -f "$relax/topol.top" ] || { echo "missing topol for $variant"; exit 2; }
work="$OUT_ROOT/$variant/md_7a6o"
mkdir -p "$work"
cp -f "$relax/topol.top" "$work/topol.top"
[ -f "$relax/posre.itp" ] && cp -f "$relax/posre.itp" "$work/posre.itp"
cp -f "$start" "$work/start.gro"
cd "$work"
echo "[$(date '+%F %T')] $variant direct start gpu=$gpu pinoffset=$pinoffset start=$start"
if [ -f md_prod.gro ] || ls md_prod.part*.gro >/dev/null 2>&1; then
  echo "[$(date '+%F %T')] $variant production already complete in $work"
  exit 0
fi
if [ -f md_prod.cpt ] && [ -f md_prod.tpr ]; then
  echo "[$(date '+%F %T')] $variant resume production from md_prod.cpt"
  CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm md_prod -cpi md_prod.cpt -noappend -ntmpi 1 -ntomp "$ntomp" -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu -pin on -pinoffset "$pinoffset" > md_prod_resume.stdout 2>&1
  echo "[$(date '+%F %T')] $variant direct production complete"
  exit 0
fi
write_nvt "$(steps_from_ps "$nvt_ps")"
"$GMX" grompp -f nvt.mdp -c start.gro -r start.gro -p topol.top -o nvt.tpr -maxwarn 5 > g_nvt.log 2>&1
CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm nvt -ntmpi 1 -ntomp "$ntomp" -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu -pin on -pinoffset "$pinoffset" > md_nvt.stdout 2>&1
write_npt "$(steps_from_ps "$npt_ps")"
"$GMX" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 5 > g_npt.log 2>&1
CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm npt -ntmpi 1 -ntomp "$ntomp" -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu -pin on -pinoffset "$pinoffset" > md_npt.stdout 2>&1
write_prod "$(steps_from_ns "$ns")"
"$GMX" grompp -f production.mdp -c npt.gro -t npt.cpt -p topol.top -o md_prod.tpr -maxwarn 5 > g_prod.log 2>&1
CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm md_prod -ntmpi 1 -ntomp "$ntomp" -gpu_id 0 -nb gpu -pme gpu -bonded gpu -update gpu -pin on -pinoffset "$pinoffset" > md_prod.stdout 2>&1
echo "[$(date '+%F %T')] $variant direct production complete"

#!/bin/bash
# Prepare a 7A6O AIM-A1 variant for steered MD (SMD) of AIM unfolding.
#
# Builds an elongated z-box around the equilibrated protein with the two
# AIM-flanking anchors (native 1262-1264 N-AIM, 1464-1466 C-AIM) aligned to z,
# re-solvates, neutralises, energy-minimises, and runs short restrained
# NVT+NPT. Output: <variant>/smd/{smd_npt.gro, smd_npt.cpt, topol_smd.top,
# index.ndx} ready for run_7a6o_smd.sh.
#
# Required inputs (present on A40 after the equilibrium MD pipeline):
#   output/gromacs_md_autoinhib/<variant>/relax_pdb/topol.top   (+ posre.itp)
#   output/gromacs_md_autoinhib/<variant>/md_7a6o/md_prod*.gro  (equilibrated)
#       (falls back to md_data/7a6o_reference_md/variants/<variant>/final.gro)
#
# Usage:  bash scripts/pipeline/prep_7a6o_smd.sh <variant> [gpu_id]
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
GMX="${GMX:-$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx}"
[ -x "$GMX" ] || GMX="$ROOT_DIR/envs/gromacs/bin/gmx"
export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"
export OCL_ICD_VENDORS="${OCL_ICD_VENDORS:-$ROOT_DIR/opencl_vendors}"
PY="${GMX_PY:-$ROOT_DIR/envs/gromacs/bin/python}"
[ -x "$PY" ] || PY=python3

variant="${1:?variant required (e.g. WT, R1306W, R1374H)}"
gpu="${2:-0}"
# Box / anchor parameters
BOX_XY="${BOX_XY:-8.2}"        # nm, x=y (protein width + padding)
BOX_Z="${BOX_Z:-11.8}"        # nm, z (folded z + extension + padding)
N_RESID="${N_RESID:-1262-1264}"   # native N-anchor residues
C_RESID="${C_RESID:-1464-1466}"   # native C-anchor residues
SALT="${SALT:-0.15}"          # mol/L NaCl
PREP_NTOMP="${PREP_NTOMP:-16}"  # keep prep jobs from monopolizing all CPU threads

OUT="$ROOT_DIR/output/gromacs_md_autoinhib/$variant"
RELAX="$OUT/relax_pdb"
WORK="$OUT/smd"
mkdir -p "$WORK"; cd "$WORK"

# --- locate inputs -----------------------------------------------------------
top_src="$RELAX/topol.top"
top_itp_dir="$RELAX"
if [ ! -f "$top_src" ] && [ -f "$OUT/md_7a6o/topol.top" ]; then
  top_src="$OUT/md_7a6o/topol.top"
  top_itp_dir="$OUT/md_7a6o"
elif [ ! -f "$top_src" ] && [ "$variant" = "WT" ] && [ -f "$ROOT_DIR/output/gromacs_md_autoinhib/7A6O_WT/md_7a6o/topol.top" ]; then
  top_src="$ROOT_DIR/output/gromacs_md_autoinhib/7A6O_WT/md_7a6o/topol.top"
  top_itp_dir="$ROOT_DIR/output/gromacs_md_autoinhib/7A6O_WT/md_7a6o"
fi
[ -f "$top_src" ] || { echo "FATAL: missing topology for $variant (run equilibrium pipeline first)"; exit 2; }
src_tpr="$OUT/md_7a6o/md_prod.tpr"
start_gro=""
for g in "$OUT"/md_7a6o/md_prod.gro $(ls -1 "$OUT"/md_7a6o/md_prod.part*.gro 2>/dev/null | sort | tail -1) \
         "$ROOT_DIR/md_data/7a6o_reference_md/variants/$variant/final.gro"; do
  [ -f "$g" ] && { start_gro="$g"; break; }
done
[ -n "$start_gro" ] || { echo "FATAL: no equilibrated start .gro for $variant"; exit 2; }
[ -f "$src_tpr" ] || src_tpr="$ROOT_DIR/md_data/7a6o_reference_md/variants/$variant/md_prod.tpr"
echo "[prep] variant=$variant  start=$start_gro  top=$top_src  gmx=$GMX"

# --- 1) whole, protein-only frame -------------------------------------------
printf 'Protein\n' | "$GMX" trjconv -s "$src_tpr" -f "$start_gro" -o prot_whole.gro -pbc whole -ur compact >prep_trjconv.log 2>&1

# --- 2) orient N->C anchor vector along +z -----------------------------------
"$PY" "$ROOT_DIR/scripts/pipeline/orient_box_smd.py" \
    --in prot_whole.gro --out prot_oriented.gro --n-resid "$N_RESID" --c-resid "$C_RESID" --name CA

# --- 3) elongated z-box, centered --------------------------------------------
"$GMX" editconf -f prot_oriented.gro -o boxed.gro -bt triclinic -box "$BOX_XY" "$BOX_XY" "$BOX_Z" -c >prep_editconf.log 2>&1

# --- 4) protein-only topology (truncate [ molecules ] to Protein* lines) ------
awk '
  BEGIN{inmol=0}
  /^\[/ {inmol=($0 ~ /molecules/)?1:0; print; next}
  { if(inmol && $1!="" && $1 !~ /^;/ && $1 !~ /^Protein/){next} print }
' "$top_src" > topol_smd.top
# bring protein include files alongside (posre.itp, any chain itp the top references)
cp -f "$top_itp_dir"/*.itp "$WORK"/ 2>/dev/null || true

# --- 5) solvate + neutralise -------------------------------------------------
"$GMX" solvate -cp boxed.gro -cs spc216 -p topol_smd.top -o solv.gro >prep_solvate.log 2>&1
cat > ions.mdp <<'EOF'
integrator = steep
nsteps     = 2000
emtol      = 1000
emstep     = 0.01
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb   = 1.2
rvdw       = 1.2
pbc        = xyz
EOF
"$GMX" grompp -f ions.mdp -c solv.gro -p topol_smd.top -o ions.tpr -maxwarn 5 >prep_grompp_ions.log 2>&1
printf 'SOL\n' | "$GMX" genion -s ions.tpr -o solv_ions.gro -p topol_smd.top \
    -pname NA -nname CL -neutral -conc "$SALT" >prep_genion.log 2>&1

# --- 6) energy minimisation --------------------------------------------------
cat > em.mdp <<'EOF'
integrator = steep
nsteps     = 5000
emtol      = 1000
emstep     = 0.01
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb   = 1.2
rvdw       = 1.2
pbc        = xyz
nstlist    = 20
EOF
"$GMX" grompp -f em.mdp -c solv_ions.gro -p topol_smd.top -o em.tpr -maxwarn 5 >prep_grompp_em.log 2>&1
CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -v -deffnm em -ntmpi 1 -ntomp "$PREP_NTOMP" -nb gpu -pin on >prep_em.log 2>&1

# --- 7) pull index groups (native numbering: gmx sees 1262..1466) ------------
# Build default groups first (System/Protein/Non-Protein/Water...) so grompp -n
# can still resolve tc-grps, THEN append the named anchor groups.
n0="${N_RESID%-*}"; n1="${N_RESID#*-}"; c0="${C_RESID%-*}"; c1="${C_RESID#*-}"
printf 'q\n' | "$GMX" make_ndx -f em.tpr -o index.ndx >prep_makendx.log 2>&1
"$GMX" select -s em.tpr -on anchors.ndx -select \
  "\"N_anchor\" group Protein and resid $n0 to $n1 and name CA; \
   \"C_anchor\" group Protein and resid $c0 to $c1 and name CA" >prep_select.log 2>&1
cat anchors.ndx >> index.ndx
echo "[prep] index groups:"; grep -E '^\[' index.ndx

# --- 8) restrained NVT + NPT (keep AIM folded before pulling) -----------------
mdp_common() { cat <<'EOF'
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
{ echo 'define              = -DPOSRES'; echo 'integrator          = md'; echo 'nsteps              = 50000'   # 100 ps
  echo 'nstxout-compressed  = 5000'; echo 'nstenergy = 5000'; echo 'nstlog = 5000'; mdp_common
  echo 'pcoupl              = no'; echo 'gen_vel = yes'; echo 'gen_temp = 310'; echo 'gen_seed = -1'; echo 'continuation = no'; } > nvt.mdp
"$GMX" grompp -f nvt.mdp -c em.gro -r em.gro -p topol_smd.top -o nvt.tpr -maxwarn 5 >prep_grompp_nvt.log 2>&1
CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm nvt -ntmpi 1 -ntomp "$PREP_NTOMP" -nb gpu -pme gpu -pin on >prep_nvt.log 2>&1

{ echo 'define              = -DPOSRES'; echo 'integrator          = md'; echo 'nsteps              = 100000'  # 200 ps
  echo 'nstxout-compressed  = 5000'; echo 'nstenergy = 5000'; echo 'nstlog = 5000'; mdp_common
  echo 'pcoupl              = C-rescale'; echo 'pcoupltype = isotropic'; echo 'tau_p = 2.0'; echo 'ref_p = 1.0'
  echo 'compressibility     = 4.5e-5'; echo 'refcoord_scaling = com'; echo 'gen_vel = no'; echo 'continuation = yes'; } > npt.mdp
"$GMX" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol_smd.top -o npt.tpr -maxwarn 5 >prep_grompp_npt.log 2>&1
CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm npt -ntmpi 1 -ntomp "$PREP_NTOMP" -nb gpu -pme gpu -pin on >prep_npt.log 2>&1

cp -f npt.gro smd_npt.gro; cp -f npt.cpt smd_npt.cpt
echo "[prep] DONE $variant -> $WORK/{smd_npt.gro, smd_npt.cpt, topol_smd.top, index.ndx}"
echo "[prep] next: bash scripts/pipeline/run_7a6o_smd.sh $variant <gpu> <nreps>"

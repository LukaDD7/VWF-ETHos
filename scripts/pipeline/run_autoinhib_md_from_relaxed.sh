#!/bin/bash
# ==============================================================================
# run_autoinhib_md_from_relaxed.sh  (闸门 3)
# ==============================================================================
# 从 relax_autoinhib_structure.sh 产出的弛豫结构 (em_10k.gro + topol.top) 做
# **受控** WT 续跑: cg 压 EM<1000 → 带约束 NVT → 带约束 NPT → 松开 → production。
# 不直接从松收敛结构无约束升温(易炸)。单变体, 过了再上 5 变体批量。
#
# GPU flags 按 $GPU_BACKEND 自动选 (CUDA/SYCL 全常驻; OpenCL 降级), 同 runner。
#
# 用法:
#   bash scripts/pipeline/run_autoinhib_md_from_relaxed.sh --variant VWF_WT --model 2
#   bash scripts/pipeline/run_autoinhib_md_from_relaxed.sh --variant VWF_WT --model 2 --ns 100
# ==============================================================================
set -u
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

VARIANT=""; MODEL=2; GMX="${GMX:-}"; GPU_ID=0
PROD_NS=50; NVT_PS=200; NPT_PS=200
while [[ $# -gt 0 ]]; do
    case $1 in
        --variant) VARIANT="$2"; shift 2 ;;
        --model)   MODEL="$2"; shift 2 ;;
        --gmx)     GMX="$2"; shift 2 ;;
        --gpu-id)  GPU_ID="$2"; shift 2 ;;
        --ns)      PROD_NS="$2"; shift 2 ;;
        --nvt-ps)  NVT_PS="$2"; shift 2 ;;
        --npt-ps)  NPT_PS="$2"; shift 2 ;;
        -h|--help) sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[WARN] unknown arg: $1"; shift ;;
    esac
done
[ -z "$VARIANT" ] && { echo "[FATAL] 需要 --variant (如 VWF_WT)"; exit 1; }

# gmx 定位 (A40 本地 env 优先)
if [ -z "$GMX" ]; then
    for c in "$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx" "$ROOT_DIR/envs/gromacs/bin/gmx" \
             /lzy/envs/gromacs/bin.AVX2_256/gmx "$(command -v gmx 2>/dev/null)"; do
        [ -n "$c" ] && [ -x "$c" ] && { GMX="$c"; break; }
    done
fi
[ -x "$GMX" ] || { echo "[FATAL] 找不到 gmx, 用 --gmx 指定"; exit 1; }
export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"
export OCL_ICD_VENDORS="${OCL_ICD_VENDORS:-$ROOT_DIR/opencl_vendors}"

RELAX="$ROOT_DIR/output/gromacs_md_autoinhib/${VARIANT}/relax_m${MODEL}"
START_GRO=""
for cand in "$RELAX/em_10k.gro" "$RELAX/em.gro"; do
    [ -f "$cand" ] && { START_GRO="$cand"; break; }
done
[ -z "$START_GRO" ] && { echo "[FATAL] 找不到弛豫产物 em_10k.gro/em.gro 于 $RELAX (先跑 relax_autoinhib_structure.sh)"; exit 1; }
[ -f "$RELAX/topol.top" ] || { echo "[FATAL] 缺 $RELAX/topol.top"; exit 1; }
[ -f "$RELAX/posre.itp" ] || echo "[WARN] 缺 posre.itp → 约束 NVT/NPT 可能无效 (检查 pdb2gmx -i)"

WORK="$ROOT_DIR/output/gromacs_md_autoinhib/${VARIANT}/md_from_relax_m${MODEL}"
mkdir -p "$WORK"
cp -f "$RELAX/topol.top" "$WORK/"
[ -f "$RELAX/posre.itp" ] && cp -f "$RELAX/posre.itp" "$WORK/"
cp -f "$START_GRO" "$WORK/start.gro"
cd "$WORK"

# ---- 后端条件化 GPU flags (同 runner) ----
GPU_BACKEND=$("$GMX" mdrun -version 2>&1 | grep "GPU support:" | awk '{print $NF}')
case "$GPU_BACKEND" in
    CUDA|SYCL) GPU_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu"; export GMX_CUDA_GRAPH=1 ;;
    *)         GPU_FLAGS="-nb gpu -pme gpu"; unset GMX_CUDA_GRAPH 2>/dev/null || true ;;
esac
echo "============================================================"
echo " WT-only 受控续跑: $VARIANT model=$MODEL"
echo " start : $START_GRO"
echo " GPU   : backend=$GPU_BACKEND flags='$GPU_FLAGS'  gpu_id=$GPU_ID"
echo " WORK  : $WORK"
echo "============================================================"
fmax_of() { grep -h "Maximum force" "$1" 2>/dev/null | tail -1 | grep -oE "[0-9.]+e[+-][0-9]+|[0-9.]+" | head -1; }
RUN() { CUDA_VISIBLE_DEVICES=$GPU_ID "$GMX" mdrun -ntmpi 1 -gpu_id "$GPU_ID" "$@"; }

NVT_STEPS=$(( NVT_PS * 500 ))     # dt=0.002 → 500 steps/ps
NPT_STEPS=$(( NPT_PS * 500 ))
PROD_STEPS=$(( PROD_NS * 500000 ))

# common mdp tail (mirror runner conventions)
mdp_tail() { cat <<EOF
dt          = 0.002
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
nstlist     = 100
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau_t       = 0.1 0.1
ref_t       = 310 310
constraints = h-bonds
constraint_algorithm = lincs
EOF
}

# ---- 0. cg 压 EM 到 Fmax<1000 ----
echo "[0] cg EM polish (emtol=1000)"
{ echo "integrator = cg"; echo "emtol = 1000"; echo "nsteps = 5000"; echo "nstlist = 10";
  echo "cutoff-scheme = Verlet"; echo "coulombtype = PME"; echo "rcoulomb = 1.2";
  echo "rvdw = 1.2"; echo "pbc = xyz"; } > em2.mdp
"$GMX" grompp -f em2.mdp -c start.gro -p topol.top -o em2.tpr -maxwarn 5 > g_em2.log 2>&1 \
    || { echo "[FATAL] grompp em2 失败"; tail -5 g_em2.log; exit 1; }
RUN -deffnm em2 > md_em2.log 2>&1
echo "   EM Fmax = $(fmax_of md_em2.log)"

# ---- 1. 带位置约束 NVT ----
echo "[1] 受约束 NVT (-DPOSRES, ${NVT_PS}ps)"
{ echo "define = -DPOSRES"; echo "integrator = md"; echo "nsteps = $NVT_STEPS";
  echo "nstxout-compressed = 5000"; echo "nstenergy = 5000"; echo "nstlog = 5000";
  mdp_tail; echo "pcoupl = no"; echo "gen_vel = yes"; echo "gen_temp = 310"; echo "gen_seed = -1";
  echo "continuation = no"; } > nvt.mdp
"$GMX" grompp -f nvt.mdp -c em2.gro -r em2.gro -p topol.top -o nvt.tpr -maxwarn 5 > g_nvt.log 2>&1 \
    || { echo "[FATAL] grompp nvt 失败"; tail -5 g_nvt.log; exit 1; }
RUN -deffnm nvt $GPU_FLAGS -pin on > md_nvt.log 2>&1 \
    || { echo "[FATAL] NVT mdrun 失败, 看 md_nvt.log (查 LINCS/温度爆炸)"; tail -8 md_nvt.log; exit 1; }
echo "   NVT done"

# ---- 2. 带位置约束 NPT ----
echo "[2] 受约束 NPT (-DPOSRES, ${NPT_PS}ps)"
{ echo "define = -DPOSRES"; echo "integrator = md"; echo "nsteps = $NPT_STEPS";
  echo "nstxout-compressed = 5000"; echo "nstenergy = 5000"; echo "nstlog = 5000";
  mdp_tail;
  echo "pcoupl = Parrinello-Rahman"; echo "pcoupltype = isotropic"; echo "tau_p = 2.0";
  echo "ref_p = 1.0"; echo "compressibility = 4.5e-5"; echo "refcoord_scaling = com";
  echo "gen_vel = no"; echo "continuation = yes"; } > npt.mdp
"$GMX" grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 5 > g_npt.log 2>&1 \
    || { echo "[FATAL] grompp npt 失败"; tail -5 g_npt.log; exit 1; }
RUN -deffnm npt $GPU_FLAGS -pin on > md_npt.log 2>&1 \
    || { echo "[FATAL] NPT mdrun 失败, 看 md_npt.log"; tail -8 md_npt.log; exit 1; }
echo "   NPT done"

# ---- 3. 松开约束 production ----
echo "[3] Production (无约束, ${PROD_NS}ns)"
{ echo "integrator = md"; echo "nsteps = $PROD_STEPS";
  echo "nstxout-compressed = 50000"; echo "nstenergy = 50000"; echo "nstlog = 50000";
  mdp_tail;
  echo "pcoupl = Parrinello-Rahman"; echo "pcoupltype = isotropic"; echo "tau_p = 2.0";
  echo "ref_p = 1.0"; echo "compressibility = 4.5e-5";
  echo "gen_vel = no"; echo "continuation = yes"; } > production.mdp
"$GMX" grompp -f production.mdp -c npt.gro -t npt.cpt -p topol.top -o md_prod.tpr -maxwarn 5 > g_prod.log 2>&1 \
    || { echo "[FATAL] grompp production 失败"; tail -5 g_prod.log; exit 1; }
RUN -deffnm md_prod $GPU_FLAGS -pin on -nstlist 200 > md_prod.log 2>&1 \
    || { echo "[FATAL] production mdrun 失败, 看 md_prod.log"; tail -8 md_prod.log; exit 1; }

echo "============================================================"
echo " ✅ WT 受控续跑完成: $WORK/md_prod.xtc"
echo " 健康检查: gmx energy 看温度~310K/压力/密度平台; 无 LINCS warning 爆炸。"
echo " 读数:    python scripts/pipeline/analyze_gromacs_md.py --md-dir output/gromacs_md_autoinhib --system autoinhib"
echo " 过了再上 5 变体批量。"
echo "============================================================"

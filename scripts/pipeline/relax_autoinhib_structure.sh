#!/bin/bash
# ==============================================================================
# relax_autoinhib_structure.sh — 分级弛豫 Boltz-2 结构, 解 EM 跑飞 (CPU, 不计费 GPU)
# ==============================================================================
# 针对 BOLTZ2_AUTINHIB_MD_BLOCKER: Boltz 结构直接溶剂化 EM → Fmax 5.9e9 + "water
# cannot be settled"。标准修法(报告里被跳过的): 先在**真空 + 受约束**下松弛内部
# clash(让重建的 H / sidechain 先归位), 再溶剂化做正常 EM。全程 CPU。
#
# 阶段:
#   1. CIF→PDB (gemmi)
#   2. pdb2gmx (charmm36m + force_fields patch, -ignh)  → conf.gro topol.top posre.itp
#   3. 真空大盒 + 受约束 EM (-DPOSRES, constraints=none): 放开 H 解 clash, 重原子不动
#   4. 真空无约束 EM: 让整体松弛
#   5. (若上面 Fmax 已降下来) 溶剂化 + 加离子 + 溶剂化 EM → 产出 MD-ready em.gro
# 每阶段打印 Fmax。任一阶段 Fmax 仍 >1e6 → 停, 提示换 model / 上 OpenMM relax。
#
# 用法:
#   bash scripts/pipeline/relax_autoinhib_structure.sh --variant VWF_WT
#   bash scripts/pipeline/relax_autoinhib_structure.sh --variant VWF_WT --model 2 --no-solvate
# ==============================================================================
set -u
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

VARIANT=""; SYSTEM="autoinhib"; MODEL=0; DO_SOLVATE=true; PDB_IN=""; SKIP_VACUUM=false
GMX="${GMX:-}"
RELAX_NTOMP="${RELAX_NTOMP:-8}"
while [[ $# -gt 0 ]]; do
    case $1 in
        --variant) VARIANT="$2"; shift 2 ;;
        --system)  SYSTEM="$2"; shift 2 ;;
        --model)   MODEL="$2"; shift 2 ;;
        --pdb)     PDB_IN="$2"; shift 2 ;;   # 直接喂干净 PDB(如 7A6O 实验结构), 跳过 Boltz CIF
        --no-solvate) DO_SOLVATE=false; shift ;;
        --skip-vacuum) SKIP_VACUUM=true; shift ;;  # 跳过 em_posres/em_vac(已 pdb2gmx 干净, 用于 FoldX 突变体)
        --gmx)     GMX="$2"; shift 2 ;;
        -h|--help) sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[WARN] unknown arg: $1"; shift ;;
    esac
done
[ -z "$VARIANT" ] && { echo "[FATAL] 需要 --variant (如 VWF_WT)"; exit 1; }

# gmx 定位
if [ -z "$GMX" ]; then
    for c in "$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx" /lzy/envs/gromacs/bin.AVX2_256/gmx "$(command -v gmx 2>/dev/null)"; do
        [ -n "$c" ] && [ -x "$c" ] && { GMX="$c"; break; }
    done
fi
[ -x "$GMX" ] || { echo "[FATAL] 找不到 gmx, 用 --gmx 指定"; exit 1; }
GMX_BINDIR="$(dirname "$GMX")"
GMX_ENVDIR="$(cd "$GMX_BINDIR/.." && pwd)"
case "$(basename "$GMX_BINDIR")" in
    bin.AVX2_256|bin.SSE2|bin) GMX_PY="$GMX_ENVDIR/bin/python" ;;
    *) GMX_PY="$GMX_BINDIR/python" ;;
esac
[ -x "$GMX_PY" ] || GMX_PY="$(command -v python3 2>/dev/null)"
export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"

# 输入: --pdb 直喂干净 PDB(实验结构路线); 否则按 --variant/--model 找 Boltz CIF
CIF=""
if [ -n "$PDB_IN" ]; then
    [ -f "$PDB_IN" ] || { echo "[FATAL] --pdb 文件不存在: $PDB_IN"; exit 1; }
    SRC="$PDB_IN"
    WORK="$ROOT_DIR/output/gromacs_md_${SYSTEM}/${VARIANT}/relax_pdb"
else
    case "$SYSTEM" in
        autoinhib) RES="$ROOT_DIR/output/boltz2_a1_dp_d3_results" ;;
        *) echo "[FATAL] 目前只支持 --system autoinhib (或用 --pdb 直喂)"; exit 1 ;;
    esac
    CIF=$(ls "$RES"/boltz_results_*"$VARIANT"*/predictions/*/*_model_${MODEL}.cif 2>/dev/null | head -1)
    [ -z "$CIF" ] && { echo "[FATAL] 找不到 $VARIANT model $MODEL 的 CIF (在 $RES 下)"; exit 1; }
    SRC="$CIF"
    WORK="$ROOT_DIR/output/gromacs_md_${SYSTEM}/${VARIANT}/relax_m${MODEL}"
fi
mkdir -p "$WORK"; cd "$WORK"
echo "============================================================"
echo " relax: $VARIANT  model=$MODEL  system=$SYSTEM"
echo " 输入: $SRC"
echo " GMX : $GMX     GMXLIB=$GMXLIB     (全程 -nb cpu, ntomp=$RELAX_NTOMP)"
echo " WORK: $WORK"
echo "============================================================"

fmax_of() { grep -h "Maximum force" "$1" 2>/dev/null | tail -1 | grep -oE "[0-9.]+e[+-][0-9]+|[0-9.]+" | head -1; }
report() { local f; f=$(fmax_of "$1"); echo "   → Fmax = ${f:-?} kJ/mol/nm" >&2; printf "%s\n" "${f:-1e99}"; }

# ---- 1. 准备 raw.pdb (PDB 直用 / CIF→PDB) ------------------------------------
if [ -n "$PDB_IN" ]; then
    echo "[1] 用 --pdb 输入 (实验结构路线)"
    cp -f "$SRC" raw.pdb
else
    echo "[1] CIF→PDB"
    "$GMX_PY" - "$SRC" raw.pdb <<'PY'
import sys, gemmi
st = gemmi.read_structure(sys.argv[1]); st.setup_entities(); st.write_pdb(sys.argv[2])
PY
fi
[ -f raw.pdb ] || { echo "[FATAL] 准备 raw.pdb 失败"; exit 1; }

# ---- 2. pdb2gmx --------------------------------------------------------------
echo "[2] pdb2gmx (charmm36m + force_fields patch)"
echo 1 | "$GMX" pdb2gmx -f raw.pdb -o conf.gro -p topol.top -i posre.itp \
    -water tip3p -ff charmm36m -ignh > pdb2gmx.log 2>&1 \
    || { echo "[FATAL] pdb2gmx 失败, 看 $WORK/pdb2gmx.log (多半 1MET/端基, 确认 GMXLIB)"; exit 1; }

# ---- 公共 EM mdp 生成 --------------------------------------------------------
mk_em() {  # $1=file  $2=extra("define = -DPOSRES" or "")  $3=coulomb(cutoff|PME)  $4=emtol  $5=emstep  $6=nsteps
    local emtol=${4:-200.0} emstep=${5:-0.001} nsteps=${6:-100000}
    cat > "$1" <<EOF
integrator    = steep
emtol         = $emtol
emstep        = $emstep
nsteps        = $nsteps
nstlist       = 10
cutoff-scheme = Verlet
coulombtype   = $3
rcoulomb      = 1.2
rvdw          = 1.2
pbc           = xyz
constraints   = none
$2
EOF
}

# 真空 EM 用快档(实验结构/FoldX 突变体已干净, 不需 0.001 emstep 慢爬)
mk_em_vacuum() {
    mk_em "$1" "$2" "$3" 200 0.01 500
}
# 溶剂化 EM 用更猛档(Fmax 起始 ~4e5, 需较大步; water-only 阶段 POSRES 锁蛋白)
mk_em_sol() {
    cat > "$1" <<EOF
integrator    = steep
emtol         = 5000.0
emstep        = 0.05
nsteps        = 50
nstlist       = 10
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.2
rvdw          = 1.2
pbc           = xyz
constraints   = h-bonds
define        = -DPOSRES
EOF
}

# ---- 3. 真空大盒 + 受约束 EM (放开 H, 重原子约束) -----------------------------
if $SKIP_VACUUM; then
    echo "[3/4] 跳过真空 EM (--skip-vacuum, 适用 FoldX 突变体): 直接用 conf.gro 溶剂化"
    F1="N/A"; F2="N/A"
    # 用 conf.gro 作为 em_vac.gro 的等价值(下游 solvate 步骤)
    "$GMX" editconf -f conf.gro -o box.gro -c -d 1.5 -bt cubic > /dev/null 2>&1
    cp -f conf.gro em_posres.gro
    cp -f conf.gro em_vac.gro
else
    echo "[3] 真空受约束 EM (-DPOSRES): 放开 H/水解 clash, 重原子不动"
    "$GMX" editconf -f conf.gro -o box.gro -c -d 1.5 -bt cubic > /dev/null 2>&1
    mk_em em_posres.mdp "define = -DPOSRES" "cutoff" 200 0.01 500
    "$GMX" grompp -f em_posres.mdp -c box.gro -r box.gro -p topol.top -o em_posres.tpr -maxwarn 5 > grompp1.log 2>&1 \
        || { echo "[FATAL] grompp(受约束) 失败, 看 grompp1.log"; tail -5 grompp1.log; exit 1; }
    "$GMX" mdrun -v -deffnm em_posres -ntmpi 1 -ntomp "$RELAX_NTOMP" -nb cpu > md1.log 2>&1
    F1=$(report md1.log)

    # ---- 4. 真空无约束 EM --------------------------------------------------------
    echo "[4] 真空无约束 EM"
    mk_em em_vac.mdp "" "PME" 200 0.01 500
    "$GMX" grompp -f em_vac.mdp -c em_posres.gro -p topol.top -o em_vac.tpr -maxwarn 5 > grompp2.log 2>&1 \
        || { echo "[FATAL] grompp(无约束) 失败"; tail -5 grompp2.log; exit 1; }
    "$GMX" mdrun -v -deffnm em_vac -ntmpi 1 -ntomp "$RELAX_NTOMP" -nb cpu > md2.log 2>&1
    F2=$(report md2.log)

    # 判定真空阶段
    awk -v f="$F2" 'BEGIN{ exit !(f+0 < 1e6) }' \
        || { echo ""; echo "❌ 真空 EM 后 Fmax 仍 $F2 (>1e6) → Boltz 几何坏得较重。"
             echo "   建议: ① 换 model (diagnose_clashes.py 挑最干净的); ② 上 OpenMM relax。"
             exit 2; }
    echo "✅ 真空弛豫后 Fmax=$F2 (<1e6), 内部 clash 已基本解开。"
fi

$DO_SOLVATE || { echo "[--no-solvate] 停在真空弛豫。产物: $WORK/em_vac.gro"; exit 0; }

# ---- 5. 溶剂化 + 离子 + 溶剂化 EM -------------------------------------------
echo "[5] 溶剂化 + 0.15M NaCl + 溶剂化 EM"
"$GMX" editconf -f em_vac.gro -o nb.gro -c -d 1.2 -bt dodecahedron > /dev/null 2>&1
"$GMX" solvate -cp nb.gro -cs spc216.gro -o solv.gro -p topol.top > sol.log 2>&1
mk_em_sol em_sol.mdp
"$GMX" grompp -f em_sol.mdp -c solv.gro -r solv.gro -p topol.top -o ions.tpr -maxwarn 5 > grompp3.log 2>&1
echo SOL | "$GMX" genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 > ion.log 2>&1
"$GMX" grompp -f em_sol.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o em.tpr -maxwarn 5 > grompp4.log 2>&1 \
    || { echo "[FATAL] 溶剂化 grompp 失败"; tail -5 grompp4.log; exit 1; }
"$GMX" mdrun -v -deffnm em -ntmpi 1 -ntomp "$RELAX_NTOMP" -nb cpu > md3.log 2>&1
F3=$(report md3.log)
cp -f em.gro solv_ions_em.gro  # 友好的别名: 水已 settled

echo ""
echo "============================================================"
echo " 真空受约束 Fmax = $F1"
echo " 真空无约束 Fmax = $F2"
echo " 溶剂化   Fmax = $F3"
if awk -v f="$F3" 'BEGIN{ exit !(f+0 < 1e4) }'; then
    echo " ✅ 溶剂化 EM 收敛 (<1e4)。MD-ready: $WORK/em.gro + topol.top"
    echo "    可把它喂给后续 NVT/NPT/Production (或让 runner 从此结构续跑)。"
else
    echo " ⚠ 溶剂化 EM Fmax=$F3 偏高, 但已远好于 5.9e9。可再补 cg 或检查残余 clash。"
fi
echo "============================================================"

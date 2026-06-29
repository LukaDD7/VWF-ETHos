#!/bin/bash
# ==============================================================================
# run_7a6o_autoinhib_md_batch.sh
#
# 7A6O AIM-A1 route runner. Starts from:
#   WT:      output/gromacs_md_autoinhib/7A6O_WT/relax_m2/nvt_soft.gro/.cpt
#   mutants: output/gromacs_md_autoinhib/<VARIANT>/relax_pdb/solv_ions_em.gro
#
# Modes:
#   --phase benchmark  Short NPT-only performance/stability check.
#   --phase all        NVT -> NPT -> Production for variants needing NVT.
#   --phase npt-prod   NPT -> Production, intended for WT after nvt_soft.
#
# Example:
#   bash scripts/pipeline/run_7a6o_autoinhib_md_batch.sh --phase benchmark --variants WT --bench-ps 2 --gpu-ids 6
#   nohup bash scripts/pipeline/run_7a6o_autoinhib_md_batch.sh --phase npt-prod --variants WT --ns 50 --gpu-ids 6 > output/gromacs_md_autoinhib/run_wt.log 2>&1 &
# ==============================================================================
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_ROOT="$ROOT_DIR/output/gromacs_md_autoinhib"

GMX="${GMX:-}"
PHASE="benchmark"
VARIANTS="WT"
GPU_IDS="6"
PROD_NS=50
NVT_PS=50
NPT_PS=200
BENCH_PS=2
NTOMP=16
MAX_PARALLEL=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --phase) PHASE="$2"; shift 2 ;;
        --variants) VARIANTS="$2"; shift 2 ;;
        --gpu-ids) GPU_IDS="$2"; shift 2 ;;
        --ns) PROD_NS="$2"; shift 2 ;;
        --nvt-ps) NVT_PS="$2"; shift 2 ;;
        --npt-ps) NPT_PS="$2"; shift 2 ;;
        --bench-ps) BENCH_PS="$2"; shift 2 ;;
        --ntomp) NTOMP="$2"; shift 2 ;;
        --max-parallel) MAX_PARALLEL="$2"; shift 2 ;;
        --gmx) GMX="$2"; shift 2 ;;
        -h|--help) sed -n '2,34p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[WARN] unknown arg: $1"; shift ;;
    esac
done

if [ -z "$GMX" ]; then
    for cand in "$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx" "$ROOT_DIR/envs/gromacs/bin/gmx" "$(command -v gmx 2>/dev/null)"; do
        [ -n "$cand" ] && [ -x "$cand" ] && { GMX="$cand"; break; }
    done
fi
[ -x "$GMX" ] || { echo "[FATAL] gmx not found; pass --gmx"; exit 1; }

export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"
export OCL_ICD_VENDORS="${OCL_ICD_VENDORS:-$ROOT_DIR/opencl_vendors}"
GPU_BACKEND=$("$GMX" mdrun -version 2>&1 | awk '/GPU support:/ {print $NF; exit}')
if [ -n "${GMX_GPU_FLAGS:-}" ]; then
    GPU_FLAGS="$GMX_GPU_FLAGS"
else
    case "$GPU_BACKEND" in
        CUDA|SYCL) GPU_FLAGS="-nb gpu -pme gpu -update cpu" ;;
        *) GPU_FLAGS="-nb gpu -pme gpu -update cpu" ;;
    esac
fi
unset GMX_CUDA_GRAPH 2>/dev/null || true

steps_from_ps() { awk -v ps="$1" 'BEGIN{printf "%d", ps * 500}'; }
steps_from_ns() { awk -v ns="$1" 'BEGIN{printf "%d", ns * 500000}'; }

mdp_common() {
    cat <<EOF
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

write_nvt_mdp() {
    local steps="$1"
    {
        echo "define              = -DPOSRES"
        echo "integrator          = md"
        echo "nsteps              = $steps"
        echo "nstxout-compressed  = 5000"
        echo "nstenergy           = 5000"
        echo "nstlog              = 5000"
        mdp_common
        echo "pcoupl              = no"
        echo "gen_vel             = yes"
        echo "gen_temp            = 310"
        echo "gen_seed            = -1"
        echo "continuation        = no"
    } > nvt.mdp
}

write_npt_mdp() {
    local steps="$1"
    {
        echo "define              = -DPOSRES"
        echo "integrator          = md"
        echo "nsteps              = $steps"
        echo "nstxout-compressed  = 5000"
        echo "nstenergy           = 5000"
        echo "nstlog              = 5000"
        mdp_common
        echo "pcoupl              = C-rescale"
        echo "pcoupltype          = isotropic"
        echo "tau_p               = 2.0"
        echo "ref_p               = 1.0"
        echo "compressibility     = 4.5e-5"
        echo "refcoord_scaling    = com"
        echo "gen_vel             = no"
        echo "continuation        = yes"
    } > npt.mdp
}

write_prod_mdp() {
    local steps="$1"
    {
        echo "integrator          = md"
        echo "nsteps              = $steps"
        echo "nstxout-compressed  = 50000"
        echo "nstenergy           = 50000"
        echo "nstlog              = 50000"
        mdp_common
        echo "pcoupl              = C-rescale"
        echo "pcoupltype          = isotropic"
        echo "tau_p               = 2.0"
        echo "ref_p               = 1.0"
        echo "compressibility     = 4.5e-5"
        echo "gen_vel             = no"
        echo "continuation        = yes"
    } > production.mdp
}

resolve_variant() {
    local variant="$1"
    if [ "$variant" = "WT" ] || [ "$variant" = "7A6O_WT" ]; then
        RELAX_DIR="$OUT_ROOT/7A6O_WT/relax_m2"
        START_GRO="$RELAX_DIR/nvt_soft.gro"
        START_CPT="$RELAX_DIR/nvt_soft.cpt"
        TOP_DIR="$RELAX_DIR"
        CAN_SKIP_NVT=true
        VARIANT_LABEL="7A6O_WT"
    else
        RELAX_DIR="$OUT_ROOT/$variant/relax_pdb"
        if [ -f "$RELAX_DIR/solv_ions_em_refined.gro" ]; then
            START_GRO="$RELAX_DIR/solv_ions_em_refined.gro"
        else
            START_GRO="$RELAX_DIR/solv_ions_em.gro"
        fi
        START_CPT=""
        TOP_DIR="$RELAX_DIR"
        CAN_SKIP_NVT=false
        VARIANT_LABEL="$variant"
    fi
    [ -f "$START_GRO" ] || { echo "[FATAL] missing start gro for $variant: $START_GRO"; return 1; }
    [ -f "$TOP_DIR/topol.top" ] || { echo "[FATAL] missing topol.top for $variant: $TOP_DIR/topol.top"; return 1; }
}

run_one() {
    local variant="$1" gpu="$2"
    resolve_variant "$variant" || return 1

    local work="$OUT_ROOT/$VARIANT_LABEL/md_7a6o"
    if [ -f "$work/md_prod.gro" ]; then
        echo "[$(date '+%F %T')] skip existing production: $work/md_prod.gro"
        return 0
    fi
    mkdir -p "$work"
    cp -f "$TOP_DIR/topol.top" "$work/topol.top"
    cp -f "$TOP_DIR"/*.itp "$work"/
    cp -f "$START_GRO" "$work/start.gro"
    [ -n "$START_CPT" ] && [ -f "$START_CPT" ] && cp -f "$START_CPT" "$work/start.cpt"

    (
        cd "$work" || exit 1
        echo "[$(date '+%F %T')] variant=$variant phase=$PHASE gpu=$gpu ntomp=$NTOMP backend=$GPU_BACKEND flags='$GPU_FLAGS'"

        local npt_input_gro="start.gro"
        local npt_input_cpt=""
        local npt_steps prod_steps nvt_steps

        if [ "$PHASE" = "all" ] || { [ "$PHASE" = "npt-prod" ] && [ "$CAN_SKIP_NVT" = false ]; }; then
            nvt_steps=$(steps_from_ps "$NVT_PS")
            write_nvt_mdp "$nvt_steps"
            "$GMX" grompp -f nvt.mdp -c start.gro -r start.gro -p topol.top -o nvt.tpr -maxwarn 5 > g_nvt.log 2>&1 || exit 2
            CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm nvt -ntmpi 1 -ntomp "$NTOMP" -gpu_id 0 $GPU_FLAGS -pin on > md_nvt.log 2>&1 || exit 3
            npt_input_gro="nvt.gro"
            npt_input_cpt="nvt.cpt"
        elif [ "$CAN_SKIP_NVT" = true ] && [ -f start.cpt ]; then
            npt_input_cpt="start.cpt"
        fi

        if [ "$PHASE" = "benchmark" ]; then
            npt_steps=$(steps_from_ps "$BENCH_PS")
        else
            npt_steps=$(steps_from_ps "$NPT_PS")
        fi
        write_npt_mdp "$npt_steps"
        if [ -n "$npt_input_cpt" ]; then
            "$GMX" grompp -f npt.mdp -c "$npt_input_gro" -r "$npt_input_gro" -t "$npt_input_cpt" -p topol.top -o npt.tpr -maxwarn 5 > g_npt.log 2>&1 || exit 4
        else
            "$GMX" grompp -f npt.mdp -c "$npt_input_gro" -r "$npt_input_gro" -p topol.top -o npt.tpr -maxwarn 5 > g_npt.log 2>&1 || exit 4
        fi
        CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm npt -ntmpi 1 -ntomp "$NTOMP" -gpu_id 0 $GPU_FLAGS -pin on > md_npt.log 2>&1 || exit 5

        if [ "$PHASE" = "benchmark" ]; then
            echo "[$(date '+%F %T')] benchmark complete: $work/npt.log"
            exit 0
        fi

        prod_steps=$(steps_from_ns "$PROD_NS")
        write_prod_mdp "$prod_steps"
        "$GMX" grompp -f production.mdp -c npt.gro -t npt.cpt -p topol.top -o md_prod.tpr -maxwarn 5 > g_prod.log 2>&1 || exit 6
        CUDA_VISIBLE_DEVICES="$gpu" "$GMX" mdrun -deffnm md_prod -ntmpi 1 -ntomp "$NTOMP" -gpu_id 0 $GPU_FLAGS -pin on > md_prod.log 2>&1 || exit 7
        echo "[$(date '+%F %T')] production complete: $work/md_prod.xtc"
    )
}

IFS=',' read -r -a VARIANT_ARR <<< "$VARIANTS"
IFS=',' read -r -a GPU_ARR <<< "$GPU_IDS"

failures=0
pids=()
wait_batch() {
    [ "${#pids[@]}" -eq 0 ] && return 0
    local pid status
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            status=$?
            echo "[ERROR] child PID $pid failed with status $status" >&2
            failures=$((failures + 1))
        fi
    done
    pids=()
}

for i in "${!VARIANT_ARR[@]}"; do
    variant="${VARIANT_ARR[$i]}"
    gpu="${GPU_ARR[$(( i % ${#GPU_ARR[@]} ))]}"
    run_one "$variant" "$gpu" > "$OUT_ROOT/${variant}_md_7a6o_${PHASE}.log" 2>&1 &
    pids+=("$!")
    if [ "${#pids[@]}" -ge "$MAX_PARALLEL" ]; then
        wait_batch
    fi
done
wait_batch
exit "$failures"

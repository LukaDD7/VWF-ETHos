#!/bin/bash
# ==============================================================================
# run_new_2b_saltbridge_md.sh
#
# End-to-end driver for the equilibrium AIM↔A1 salt-bridge MD of the new-2B A1
# variants on 7A6O (AIM-A1 X-ray WT backbone). Ties together the existing stage
# scripts; nothing here re-implements MD. Stages (each resumable / skippable):
#
#   build    FoldX BuildModel on 7A6O WT  -> structures/7a6o_mutants/<v>.pdb
#            (build_2b_mutants_foldx.py --detect-offset)
#   relax    pdb2gmx + vacuum/solv EM     -> output/gromacs_md_autoinhib/<v>/relax_pdb/
#            (relax_autoinhib_structure.sh, per variant, CPU)
#   md       NVT -> NPT -> production      -> .../<v>/md_7a6o/md_prod.xtc
#            (run_7a6o_variant_direct.sh, GPU pool round-robin)
#   extract  salt-bridge occupancy + z     -> output/md_7a6o_saltbridge_features.csv
#            (extract_aim_saltbridge_features.py over WT + all variants)
#
# WHY equilibrium (not SMD): SMD force-axis was a no-go (see
# docs/7A6O_SMD_LITERATURE_NOGO_2026-06-24.md). This route measures whether the
# autoinhibitory salt bridges (D1269-R1306, D1269-R1450) stay intact at rest;
# 2B loosens them -> aim_sb_retained_z separates 2B from 2M without pulling.
#
# Usage (A40, GPUs 4-7):
#   export FOLDX=/path/to/foldx
#   bash scripts/pipeline/run_new_2b_saltbridge_md.sh --preflight
#   nohup bash scripts/pipeline/run_new_2b_saltbridge_md.sh \
#         --stage all --gpu-ids 4,5,6,7 --ns 50 \
#         > output/gromacs_md_autoinhib/new_2b_driver.log 2>&1 &
#
# Run a single stage:  --stage build | relax | md | extract
# ==============================================================================
set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_ROOT="$ROOT_DIR/output/gromacs_md_autoinhib"
MUT_DIR="$ROOT_DIR/structures/7a6o_mutants"
WT_PDB="$ROOT_DIR/structures/7A6O_AIM_A1_clean.pdb"

STAGE="all"
GPU_IDS="4,5,6,7"
VARIANTS_FILE="$ROOT_DIR/output/new_2b_a1_variants.txt"
NS=50
NVT_PS=50
NPT_PS=200
NTOMP=8
RELAX_JOBS=6
RELAX_NTOMP=8
FOLDX_BIN="${FOLDX:-foldx}"
PY="${PY:-$ROOT_DIR/envs/gromacs/bin/python3}"
GMX="${GMX:-$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx}"
DO_PREFLIGHT=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --stage) STAGE="$2"; shift 2 ;;
        --gpu-ids) GPU_IDS="$2"; shift 2 ;;
        --variants-file) VARIANTS_FILE="$2"; shift 2 ;;
        --ns) NS="$2"; shift 2 ;;
        --nvt-ps) NVT_PS="$2"; shift 2 ;;
        --npt-ps) NPT_PS="$2"; shift 2 ;;
        --ntomp) NTOMP="$2"; shift 2 ;;
        --relax-jobs) RELAX_JOBS="$2"; shift 2 ;;
        --relax-ntomp) RELAX_NTOMP="$2"; shift 2 ;;
        --foldx) FOLDX_BIN="$2"; shift 2 ;;
        --gmx) GMX="$2"; shift 2 ;;
        --preflight) DO_PREFLIGHT=true; shift ;;
        -h|--help) sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[WARN] unknown arg: $1"; shift ;;
    esac
done

log() { echo "[$(date '+%F %T')] $*"; }

# variant list: first whitespace token of each non-comment line
read_variants() {
    awk 'NF && $1 !~ /^#/ {print $1}' "$VARIANTS_FILE"
}

IFS=',' read -r -a GPU_ARR <<< "$GPU_IDS"

preflight() {
    log "=== preflight ==="
    local ok=true
    [ -f "$VARIANTS_FILE" ] && log "variants file: $VARIANTS_FILE ($(read_variants | wc -l | tr -d ' ') variants)" || { log "[FAIL] no variants file: $VARIANTS_FILE"; ok=false; }
    [ -f "$WT_PDB" ] && log "WT PDB: $WT_PDB" || { log "[FAIL] no WT PDB: $WT_PDB  (run: $PY scripts/pipeline/fetch_clean_7a6o.py)"; ok=false; }
    if command -v "$FOLDX_BIN" >/dev/null 2>&1 || [ -x "$FOLDX_BIN" ]; then log "FoldX: $FOLDX_BIN"; else log "[FAIL] FoldX not found: $FOLDX_BIN  (export FOLDX=/path/to/foldx)"; ok=false; fi
    [ -x "$GMX" ] && log "gmx: $GMX" || { log "[FAIL] gmx not executable: $GMX"; ok=false; }
    if [ -x "$GMX" ]; then
        local backend; backend=$("$GMX" mdrun -version 2>&1 | awk '/GPU support:/ {print $NF; exit}')
        log "gmx GPU backend: ${backend:-unknown}"
    fi
    command -v nvidia-smi >/dev/null 2>&1 && log "GPUs: $(nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null | tr '\n' ' ')" || log "[WARN] nvidia-smi not found"
    log "GPU pool for MD: ${GPU_IDS}"
    "$GMX" --help >/dev/null 2>&1 # noop
    $PY -c "import gemmi" 2>/dev/null \
        && log "python dep: gemmi OK (build)" \
        || { log "[FAIL] need gemmi in $PY env for build"; ok=false; }
    $PY -c "import MDAnalysis" 2>/dev/null \
        && log "python dep: MDAnalysis OK (extract)" \
        || log "[WARN] need MDAnalysis in $PY env before extract"
    $ok && log "preflight PASS" || { log "preflight FAIL"; return 1; }
}

stage_build() {
    log "=== stage build (FoldX) ==="
    $PY "$SCRIPT_DIR/build_2b_mutants_foldx.py" \
        --wt "$WT_PDB" --foldx "$FOLDX_BIN" --detect-offset \
        --variants-file "$VARIANTS_FILE" --out-dir "$MUT_DIR" || return 1
    log "built mutants -> $MUT_DIR"
}

stage_relax() {
    log "=== stage relax (CPU jobs=${RELAX_JOBS}, ntomp=${RELAX_NTOMP}) ==="
    local -a pids=() names=()
    local v pdb skip_n=0 failed=0
    while read -r v; do
        pdb="$MUT_DIR/$v.pdb"
        if [ ! -f "$pdb" ]; then log "[skip] no mutant PDB (build skipped it?): $pdb"; skip_n=$((skip_n+1)); continue; fi
        if [ -f "$OUT_ROOT/$v/relax_pdb/solv_ions_em.gro" ] && [ -f "$OUT_ROOT/$v/relax_pdb/topol.top" ]; then
            log "[skip] already relaxed: $v"; skip_n=$((skip_n+1)); continue
        fi
        log "launch relax $v"
        RELAX_NTOMP="$RELAX_NTOMP" GMX="$GMX" \
            bash "$SCRIPT_DIR/relax_autoinhib_structure.sh" --pdb "$pdb" --variant "$v" \
            > "$OUT_ROOT/${v}_relax.log" 2>&1 &
        pids+=("$!"); names+=("$v")
        if [ "${#pids[@]}" -ge "$RELAX_JOBS" ]; then
            wait "${pids[0]}" || { log "[ERROR] relax failed: ${names[0]}"; failed=$((failed+1)); }
            pids=("${pids[@]:1}"); names=("${names[@]:1}")
        fi
    done < <(read_variants)
    while [ "${#pids[@]}" -gt 0 ]; do
        wait "${pids[0]}" || { log "[ERROR] relax failed: ${names[0]}"; failed=$((failed+1)); }
        pids=("${pids[@]:1}"); names=("${names[@]:1}")
    done
    log "relax complete skip=$skip_n failed=$failed"
    [ "$failed" -eq 0 ]
}

stage_md() {
    log "=== stage md (GPU pool ${GPU_IDS}) ==="
    local -a pids=() ; local i=0 v gpu
    while read -r v; do
        if [ ! -f "$OUT_ROOT/$v/relax_pdb/solv_ions_em.gro" ]; then log "[skip] not relaxed: $v"; continue; fi
        if [ -f "$OUT_ROOT/$v/md_7a6o/md_prod.gro" ] || ls "$OUT_ROOT/$v/md_7a6o/"md_prod.part*.gro >/dev/null 2>&1; then
            log "[skip] MD done: $v"; continue
        fi
        gpu="${GPU_ARR[$(( i % ${#GPU_ARR[@]} ))]}"
        log "launch MD $v on GPU $gpu"
        NS="$NS" NVT_PS="$NVT_PS" NPT_PS="$NPT_PS" NTOMP="$NTOMP" \
            bash "$SCRIPT_DIR/run_7a6o_variant_direct.sh" "$v" "$gpu" \
            > "$OUT_ROOT/${v}_md_direct.log" 2>&1 &
        pids+=("$!")
        i=$((i+1))
        # keep at most #GPU jobs in flight
        if [ "${#pids[@]}" -ge "${#GPU_ARR[@]}" ]; then
            wait "${pids[0]}" 2>/dev/null || log "[ERROR] a MD job failed"
            pids=("${pids[@]:1}")
        fi
    done < <(read_variants)
    for p in "${pids[@]}"; do wait "$p" 2>/dev/null || log "[ERROR] a MD job failed"; done
    log "md stage complete"
}

stage_extract() {
    log "=== stage extract (salt-bridge features + z) ==="
    # extract_aim_saltbridge_features.py expects  <input>/<v>/{md_prod.tpr,prod_concat.xtc}
    # our MD writes                               <v>/md_7a6o/{md_prod.tpr,md_prod.xtc}
    # -> build a symlink staging dir that flattens md_7a6o/ and renames the xtc,
    #    and INCLUDES the existing reference set so the WT/2B/2M z-baseline is intact.
    local STAGE_DIR="$OUT_ROOT/saltbridge_staging"
    mkdir -p "$STAGE_DIR"

    # (a) existing reference trajectories (WT + known 2B/2M), if present
    local REF="$ROOT_DIR/md_data/7a6o_reference_md/variants"
    if [ -d "$REF" ]; then
        local d v
        for d in "$REF"/*/; do
            v="$(basename "$d")"
            [ -f "$d/md_prod.tpr" ] && [ -f "$d/prod_concat.xtc" ] || continue
            mkdir -p "$STAGE_DIR/$v"
            ln -sf "$d/md_prod.tpr" "$STAGE_DIR/$v/md_prod.tpr"
            ln -sf "$d/prod_concat.xtc" "$STAGE_DIR/$v/prod_concat.xtc"
        done
    fi

    # (b) new variants + WT from the autoinhib MD tree (WT dir is 7A6O_WT -> stage as WT)
    local v md sv
    for v in $(read_variants) 7A6O_WT; do
        md="$OUT_ROOT/$v/md_7a6o"
        [ -f "$md/md_prod.tpr" ] || continue
        sv="$v"; [ "$v" = "7A6O_WT" ] && sv="WT"
        mkdir -p "$STAGE_DIR/$sv"
        ln -sf "$md/md_prod.tpr" "$STAGE_DIR/$sv/md_prod.tpr"
        if [ -f "$md/md_prod.xtc" ]; then
            ln -sf "$md/md_prod.xtc" "$STAGE_DIR/$sv/prod_concat.xtc"
        elif ls "$md"/md_prod.part*.xtc >/dev/null 2>&1; then
            log "concat resumed parts for $sv"
            "$GMX" trjcat -f "$md"/md_prod*.xtc -o "$STAGE_DIR/$sv/prod_concat.xtc" -cat >/dev/null 2>&1 \
                || log "  [WARN] trjcat failed for $sv"
        fi
    done

    log "staged $(find "$STAGE_DIR" -name md_prod.tpr | wc -l | tr -d ' ') variants -> $STAGE_DIR"
    $PY "$SCRIPT_DIR/extract_aim_saltbridge_features.py" \
        --input "$STAGE_DIR" \
        --output "$ROOT_DIR/output/md_7a6o_saltbridge_features.csv" || return 1
    log "features -> output/md_7a6o_saltbridge_features.csv"
    log "next: $PY scripts/pipeline/validate_2b_recall_saltbridge.py (cluster matrix)"
    log "  or: $PY scripts/pipeline/validate_new_2b_panel_recall.py (re-run with sb joined)"
}

if $DO_PREFLIGHT; then
    preflight
    exit $?
fi

case "$STAGE" in
    build)   stage_build ;;
    relax)   stage_relax ;;
    md)      stage_md ;;
    extract) stage_extract ;;
    all)     preflight && stage_build && stage_relax && stage_md && stage_extract ;;
    *) echo "[FATAL] unknown stage: $STAGE (build|relax|md|extract|all)"; exit 1 ;;
esac

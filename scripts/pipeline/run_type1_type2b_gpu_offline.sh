#!/usr/bin/env bash
# Run the merged Type-1 and Type-2B computational panels on an offline GPU node.
#
# CPU-side prerequisites are already prepared:
#   * AlphaGenome results for 9 Type-1 and 7 Type-2B requests.
#   * 7A6O WT and FoldX mutant PDBs for targeted MD.
#
# Usage on the GPU node:
#   GPUS=7 bash scripts/pipeline/run_type1_type2b_gpu_offline.sh
#
# Optional:
#   GPU_IDS=0,1,2,3 GPUS=4 bash scripts/pipeline/run_type1_type2b_gpu_offline.sh
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
OUT_DIR="$ROOT_DIR/outputs/computational_panel_20260829"
BOLTZ_YAML="$OUT_DIR/boltz/yamls"
BOLTZ_RAW="$OUT_DIR/boltz/raw"
MD_OUT="$OUT_DIR/md"
GPUS="${GPUS:-7}"
GPU_IDS="${GPU_IDS:-}"

export GMX="${GMX:-/lzy/envs/gromacs/bin.AVX2_256/gmx}"
export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"
export CC="${CC:-$(command -v gcc || true)}"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

if [ ! -x "$GMX" ]; then
    echo "[FATAL] GROMACS executable not found: $GMX" >&2
    exit 2
fi

mkdir -p "$BOLTZ_YAML" "$BOLTZ_RAW" "$MD_OUT"

# 1. Combine the two pre-generated Boltz YAML panels (15 + 20 = 35 jobs).
log "Preparing combined Boltz YAML panel"
cp -f "$ROOT_DIR"/outputs/type1_panel_agent_20260828/server_bundle/boltz/yamls/*.yaml "$BOLTZ_YAML"/
cp -f "$ROOT_DIR"/outputs/type2b_panel_agent_20260829/server_bundle/boltz/yamls/*.yaml "$BOLTZ_YAML"/
N_YAML=$(find "$BOLTZ_YAML" -maxdepth 1 -name '*.yaml' | wc -l)
if [ "$N_YAML" -ne 35 ]; then
    echo "[FATAL] Expected 35 Boltz YAML files; found $N_YAML" >&2
    exit 2
fi

# 2. Run Boltz once for both panels. The runner is incremental and skips .done jobs.
log "Running combined Boltz-2 panel with $GPUS GPU worker(s)"
export PATH="/lzy/envs/boltz2/bin:$PATH"
BOLTZ_ARGS=(
    --input-dir "$BOLTZ_YAML"
    --out-dir "$BOLTZ_RAW"
    --gpus "$GPUS"
)
if [ -n "$GPU_IDS" ]; then
    BOLTZ_ARGS+=(--gpu-ids "$GPU_IDS")
fi
bash "$ROOT_DIR/scripts/pipeline/run_vwd_functional_boltz2_panel.sh" "${BOLTZ_ARGS[@]}"

# 3. Parse the shared raw output once for each batch manifest.
log "Parsing Boltz-2 results"
for python_bin in /lzy/envs/tools/bin/python /lzy/envs/boltz2/bin/python python3; do
    if command -v "$python_bin" >/dev/null 2>&1; then
        PYTHON_BIN="$python_bin"
        break
    fi
done
export PATH="/lzy/envs/tools/bin:/lzy/envs/boltz2/bin:$PATH"

"$PYTHON_BIN" "$ROOT_DIR/scripts/pipeline/parse_vwd_functional_boltz2_results.py" \
    --results-dir "$BOLTZ_RAW" \
    --manifest "$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle/boltz/job_manifest.csv" \
    --output "$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle/results/boltz/boltz_results_summary.csv"

"$PYTHON_BIN" "$ROOT_DIR/scripts/pipeline/parse_vwd_functional_boltz2_results.py" \
    --results-dir "$BOLTZ_RAW" \
    --manifest "$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle/boltz/job_manifest.csv" \
    --output "$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle/results/boltz/boltz_results_summary.csv"

# 4. Relax the CPU-prebuilt 7A6O WT and mutant structures.
log "Preparing 7A6O MD structures"
WT_STRUCT="$ROOT_DIR/structures/7a6o_mutants/WT_Repair.pdb"
if [ ! -f "$WT_STRUCT" ]; then
    echo "[FATAL] Missing prebuilt WT structure: $WT_STRUCT" >&2
    exit 2
fi
cp -f "$WT_STRUCT" "$ROOT_DIR/structures/7a6o_mutants/PANEL_WT_MATCHED.pdb"

MD_VARIANTS=(PANEL_WT_MATCHED P1413L R1308C S1310F V1316M R1341W A1461D)
for variant in "${MD_VARIANTS[@]}"; do
    pdb="$ROOT_DIR/structures/7a6o_mutants/${variant}.pdb"
    relaxed="$ROOT_DIR/output/gromacs_md_autoinhib/${variant}/relax_pdb/solv_ions_em_refined.gro"
    if [ ! -f "$pdb" ]; then
        echo "[FATAL] Missing prebuilt mutant structure: $pdb" >&2
        exit 2
    fi
    if [ ! -f "$relaxed" ]; then
        log "Relaxing $variant"
        bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh" \
            --pdb "$pdb" --variant "$variant" --skip-vacuum --gmx "$GMX"
    else
        log "Relaxation already complete for $variant"
    fi
done

# 5. Run 50-ns MD in GPU waves. PANEL_WT_MATCHED is shared by both panels.
log "Running targeted 50-ns MD"
if [ -n "$GPU_IDS" ]; then
    IFS=',' read -r -a gpu_ids <<< "$GPU_IDS"
else
    gpu_ids=()
    for ((i=0; i<GPUS; i++)); do
        gpu_ids+=("$i")
    done
fi
if [ "${#gpu_ids[@]}" -eq 0 ]; then
    echo "[FATAL] No GPU IDs configured" >&2
    exit 2
fi

batch_size="${#gpu_ids[@]}"
for ((start=0; start<${#MD_VARIANTS[@]}; start+=batch_size)); do
    pids=()
    for ((j=start; j<start+batch_size && j<${#MD_VARIANTS[@]}; j++)); do
        variant="${MD_VARIANTS[$j]}"
        gpu="${gpu_ids[$((j-start))]}"
        log "Starting MD: $variant on GPU $gpu"
        NS=50 GMX="$GMX" bash "$ROOT_DIR/scripts/pipeline/run_7a6o_variant_direct.sh" "$variant" "$gpu" &
        pids+=("$!")
    done
    for pid in "${pids[@]}"; do
        wait "$pid"
    done
done

log "Analyzing merged MD trajectories"
"$PYTHON_BIN" "$ROOT_DIR/scripts/pipeline/analyze_7a6o_completed_md.py" \
    --input "$ROOT_DIR/output/gromacs_md_autoinhib" \
    --output "$MD_OUT" \
    --variants "$(IFS=,; echo "${MD_VARIANTS[*]}")" \
    --wt-variant PANEL_WT_MATCHED \
    --gmx "$GMX" \
    --force

mkdir -p \
    "$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle/results/md" \
    "$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle/results/md"
cp -f "$MD_OUT/aim_a1_contacts_summary.csv" \
    "$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle/results/md/aim_a1_contacts_summary.csv"
cp -f "$MD_OUT/aim_a1_contacts_summary.csv" \
    "$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle/results/md/aim_a1_contacts_summary.csv"

# 6. Validate inventories, build evidence rows, and run the agent for both batches.
log "Ingesting Type-1 results and running agent"
bash "$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle/ingest_and_test_agent.sh"

log "Ingesting Type-2B results and running agent"
bash "$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle/ingest_and_test_agent.sh"

log "GPU offline run complete"
log "Type-1 agent output:  outputs/type1_panel_agent_20260828/server_results/agent_run"
log "Type-2B agent output: outputs/type2b_panel_agent_20260829/server_results/agent_run"

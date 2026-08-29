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

first_existing_executable() {
    local candidate
    for candidate in "$@"; do
        if [ -n "$candidate" ] && [ -x "$candidate" ]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

first_existing_directory() {
    local candidate
    for candidate in "$@"; do
        if [ -n "$candidate" ] && [ -d "$candidate" ]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

NFS_ROOT="${ROOT_DIR%/projects/*}"
if [ "$NFS_ROOT" = "$ROOT_DIR" ]; then
    NFS_ROOT="/inspire/hdd/global_user/mengweicheng-240108120092/lzy"
fi
GMX_CANDIDATES=(
    "${GMX:-}"
    "$ROOT_DIR/envs/gromacs-cuda/bin.AVX2_256/gmx"
    "$ROOT_DIR/envs/gromacs-cuda/bin/gmx"
    "$ROOT_DIR/../envs/gromacs-cuda/bin.AVX2_256/gmx"
    "$ROOT_DIR/../envs/gromacs-cuda/bin/gmx"
    "$ROOT_DIR/../../envs/gromacs-cuda/bin.AVX2_256/gmx"
    "$ROOT_DIR/../../envs/gromacs-cuda/bin/gmx"
    "/lzy/envs/gromacs-cuda/bin.AVX2_256/gmx"
    "/lzy/envs/gromacs-cuda/bin/gmx"
    "$NFS_ROOT/envs/gromacs-cuda/bin.AVX2_256/gmx"
    "$NFS_ROOT/envs/gromacs-cuda/bin/gmx"
    "$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx"
    "$ROOT_DIR/envs/gromacs/bin/gmx"
    "$ROOT_DIR/../envs/gromacs/bin.AVX2_256/gmx"
    "$ROOT_DIR/../envs/gromacs/bin/gmx"
    "$ROOT_DIR/../../envs/gromacs/bin.AVX2_256/gmx"
    "$ROOT_DIR/../../envs/gromacs/bin/gmx"
    "/lzy/envs/gromacs/bin.AVX2_256/gmx"
    "/lzy/envs/gromacs/bin/gmx"
    "$NFS_ROOT/envs/gromacs/bin.AVX2_256/gmx"
    "$NFS_ROOT/envs/gromacs/bin/gmx"
    "$(command -v gmx 2>/dev/null || true)"
)
GMX="$(first_existing_executable "${GMX_CANDIDATES[@]}")" || {
    echo "[FATAL] GROMACS executable not found. Tried:" >&2
    printf '  %s\n' "${GMX_CANDIDATES[@]}" >&2
    echo "Set GMX=/absolute/path/to/gmx and rerun." >&2
    exit 2
}
export GMXLIB="${GMXLIB:-$ROOT_DIR/force_fields}"
export CC="${CC:-$(command -v gcc || true)}"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

BOLTZ_COMMAND="$(command -v boltz 2>/dev/null || true)"
BOLTZ_ENV_CANDIDATES=("$ROOT_DIR/envs/boltz2" "/lzy/envs/boltz2" "$NFS_ROOT/envs/boltz2")
BOLTZ_ENV_CANDIDATES+=("$ROOT_DIR/../envs/boltz2" "$ROOT_DIR/../../envs/boltz2")
if [ -n "$BOLTZ_COMMAND" ]; then
    BOLTZ_ENV_CANDIDATES=("$(dirname "$(dirname "$BOLTZ_COMMAND")")" "${BOLTZ_ENV_CANDIDATES[@]}")
fi
BOLTZ_ENV="$(first_existing_directory "${BOLTZ_ENV_CANDIDATES[@]}")" || {
    echo "[FATAL] Boltz-2 environment not found. Set PATH to its bin directory and rerun." >&2
    exit 2
}

TOOLS_PYTHON_CANDIDATES=(
    "$ROOT_DIR/envs/tools/bin/python"
    "$ROOT_DIR/../envs/tools/bin/python"
    "$ROOT_DIR/../../envs/tools/bin/python"
    "$ROOT_DIR/miniconda3/envs/tools/bin/python"
    "$ROOT_DIR/../miniconda3/envs/tools/bin/python"
    "$ROOT_DIR/../../miniconda3/envs/tools/bin/python"
    "/lzy/envs/tools/bin/python"
    "/lzy/miniconda3/envs/tools/bin/python"
    "$NFS_ROOT/envs/tools/bin/python"
    "$NFS_ROOT/miniconda3/envs/tools/bin/python"
    "$BOLTZ_ENV/bin/python"
    "$(command -v python 2>/dev/null || true)"
    "$(command -v python3 2>/dev/null || true)"
)
TOOLS_PYTHON="$(first_existing_executable "${TOOLS_PYTHON_CANDIDATES[@]}")" || {
    echo "[FATAL] Python environment not found. Set PATH to a Python environment with pandas/langgraph." >&2
    exit 2
}
TOOLS_ENV="$(dirname "$(dirname "$TOOLS_PYTHON")")"

log "Using GROMACS: $GMX"
log "Using Boltz environment: $BOLTZ_ENV"
log "Using Python environment: $TOOLS_ENV"
if [ "${CHECK_ENV:-0}" = "1" ]; then
    log "Environment check passed"
    exit 0
fi

mkdir -p "$BOLTZ_YAML" "$BOLTZ_RAW" "$MD_OUT"

# 1. Combine the two pre-generated Boltz YAML panels. The manifests contain 35
# job references, but three A1 WT baselines are shared by both panels, so the
# deduplicated offline workload is 32 unique Boltz jobs.
log "Preparing combined Boltz YAML panel"
cp -f "$ROOT_DIR"/outputs/type1_panel_agent_20260828/server_bundle/boltz/yamls/*.yaml "$BOLTZ_YAML"/
cp -f "$ROOT_DIR"/outputs/type2b_panel_agent_20260829/server_bundle/boltz/yamls/*.yaml "$BOLTZ_YAML"/
N_YAML=$(find "$BOLTZ_YAML" -maxdepth 1 -name '*.yaml' | wc -l)
if [ "$N_YAML" -ne 32 ]; then
    echo "[FATAL] Expected 32 unique Boltz YAML files; found $N_YAML" >&2
    exit 2
fi

# 2. Run Boltz once for both panels. The runner is incremental and skips .done jobs.
log "Running combined Boltz-2 panel with $GPUS GPU worker(s)"
export PATH="$BOLTZ_ENV/bin:$PATH"
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
PYTHON_BIN="$TOOLS_PYTHON"
export PATH="$TOOLS_ENV/bin:$BOLTZ_ENV/bin:$PATH"

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

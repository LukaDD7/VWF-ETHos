#!/usr/bin/env python3
"""Build a portable server bundle for the six-patient Type-2B workbook."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
A1_VARIANTS = ["R1308C", "S1310F", "V1316M", "R1341W", "A1461D"]


def write(path: Path, content: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    if executable:
        path.chmod(0o755)


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-package", type=Path, default=ROOT / "outputs/type2b_panel_agent_20260829/input_package")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/type2b_panel_agent_20260829/server_bundle")
    args = parser.parse_args()
    bundle = args.output_dir.resolve()
    input_dir = bundle / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "alphagenome_requests.csv", "alphagenome_requests.vcf", "boltz_variants.csv",
        "first_level.csv", "genetics.csv", "clinical_context.csv",
        "panel_request_status.csv", "normalization_audit.csv", "prepared_payload.json",
    ):
        shutil.copy2(args.input_package / name, input_dir / name)

    alpha = pd.read_csv(input_dir / "alphagenome_requests.csv")
    boltz = pd.read_csv(input_dir / "boltz_variants.csv")
    if len(alpha) != 7 or alpha["patient_id"].nunique() != 6:
        raise ValueError("Expected seven AlphaGenome requests across six patients")
    if len(boltz) != 7 or boltz["aa_change"].nunique() != 6:
        raise ValueError("Expected seven Boltz mappings and six unique missense variants")

    boltz_dir = bundle / "boltz"
    subprocess.run([
        "python", str(ROOT / "scripts/pipeline/generate_vwd_functional_boltz2_yamls.py"),
        "--variants-csv", str(input_dir / "boltz_variants.csv"),
        "--output-dir", str(boltz_dir), "--minimal-wt-baselines",
    ], cwd=ROOT, check=True)
    manifest = pd.read_csv(boltz_dir / "job_manifest.csv")
    if len(manifest) != 20 or manifest["run_decision"].eq("WT_BASELINE").sum() != 4:
        raise ValueError("Focused Type-2B Boltz inventory must be 4 WT + 16 variant jobs")

    pd.DataFrame([
        {"aa_change": name, "system": "7A6O AIM-A1 closed-state equilibrium MD", "matched_wt": "TYPE2B_WT_MATCHED", "production_ns": 50, "priority": "P0"}
        for name in A1_VARIANTS
    ]).to_csv(bundle / "md_priority.csv", index=False)

    write(bundle / "run_alphagenome.sh", '''#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle"
: "${ALPHAGENOME_API_KEY:?export ALPHAGENOME_API_KEY before running}"
python "$ROOT_DIR/scripts/pipeline/run_type1_10case_alphagenome.py" \
  --requests "$BUNDLE/input/alphagenome_requests.csv" \
  --output-dir "$BUNDLE/results/alphagenome" --expected-requests 7
''', executable=True)

    write(bundle / "run_boltz.sh", '''#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle"
GPUS="${GPUS:-4}"
bash "$ROOT_DIR/scripts/pipeline/run_vwd_functional_boltz2_panel.sh" \
  --input-dir "$BUNDLE/boltz/yamls" --out-dir "$BUNDLE/results/boltz/raw" --gpus "$GPUS"
python "$ROOT_DIR/scripts/pipeline/parse_vwd_functional_boltz2_results.py" \
  --results-dir "$BUNDLE/results/boltz/raw" --manifest "$BUNDLE/boltz/job_manifest.csv" \
  --output "$BUNDLE/results/boltz/boltz_results_summary.csv"
''', executable=True)

    write(bundle / "run_md_a1_panel.sh", '''#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle"
: "${FOLDX:?export FOLDX=/absolute/path/to/foldx}"
VARIANTS=(R1308C S1310F V1316M R1341W A1461D)
python "$ROOT_DIR/scripts/pipeline/build_2b_mutants_foldx.py" \
  --wt "$ROOT_DIR/structures/7A6O_AIM_A1_clean.pdb" --foldx "$FOLDX" --detect-offset \
  --variants "$(IFS=,; echo "${VARIANTS[*]}")" --out-dir "$ROOT_DIR/structures/7a6o_mutants"
cp "$ROOT_DIR/structures/7a6o_mutants/WT_Repair.pdb" "$ROOT_DIR/structures/7a6o_mutants/TYPE2B_WT_MATCHED.pdb"
for variant in TYPE2B_WT_MATCHED "${VARIANTS[@]}"; do
  bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh" \
    --pdb "$ROOT_DIR/structures/7a6o_mutants/${variant}.pdb" --variant "$variant" --skip-vacuum
done
# Six independent 50-ns jobs; GPU IDs may be changed to match the server.
jobs=(TYPE2B_WT_MATCHED R1308C S1310F V1316M R1341W A1461D)
gpu_ids=(${GPU_IDS:-0 1 2 3 4 5})
pids=()
for i in "${!jobs[@]}"; do
  NS=50 bash "$ROOT_DIR/scripts/pipeline/run_7a6o_variant_direct.sh" "${jobs[$i]}" "${gpu_ids[$i]}" &
  pids+=("$!")
done
for pid in "${pids[@]}"; do wait "$pid"; done
python "$ROOT_DIR/scripts/pipeline/analyze_7a6o_completed_md.py" \
  --input "$ROOT_DIR/output/gromacs_md_autoinhib" --output "$BUNDLE/results/md" \
  --variants "TYPE2B_WT_MATCHED,$(IFS=,; echo "${VARIANTS[*]}")" --wt-variant TYPE2B_WT_MATCHED
''', executable=True)

    write(bundle / "ingest_and_test_agent.sh", '''#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_bundle"
OUT="$ROOT_DIR/outputs/type2b_panel_agent_20260829/server_results"
mkdir -p "$OUT"
python "$ROOT_DIR/scripts/pipeline/ingest_type1_10case_results.py" \
  --alphagenome-features "$BUNDLE/results/alphagenome/alphagenome_agent_features.csv" \
  --boltz-summary "$BUNDLE/results/boltz/boltz_results_summary.csv" \
  --boltz-manifest "$BUNDLE/boltz/job_manifest.csv" \
  --boltz-variants "$BUNDLE/input/boltz_variants.csv" \
  --md-features "$BUNDLE/results/md/aim_a1_contacts_summary.csv" \
  --expected-alpha 7 --expected-boltz-jobs 20 --expected-boltz-wt 4 \
  --output "$OUT/type2b_computational_evidence.csv"
python "$ROOT_DIR/scripts/run_vwd_langgraph_v0.py" \
  --workbook "$ROOT_DIR/outputs/type2b_panel_agent_20260829/type2b_panel_agent_package.xlsx" \
  --output-dir "$OUT/agent_run" --computational-panels
''', executable=True)

    write(bundle / "README.md", """# Type-2B six-patient server bundle

- Deidentified input: six patients, seven variant records; patient 5 has A1461D plus benign-control D2449N.
- AlphaGenome: 7 requests x 11 official raw outputs, metadata-selected VWF-relevant tracks and all recommended scorer views.
- Boltz-2: 6 unique missense variants; 16 variant-assay jobs + 4 assay-matched WT baselines = 20 jobs.
- MD: five unique A1 variants plus one pipeline-matched WT, 50 ns each. D2449N (C3 benign control) is not sent to A1 MD.
- V1316M is computed once and mapped back to both patients.

Run `run_alphagenome.sh`, `run_boltz.sh`, and `run_md_a1_panel.sh`; then run `ingest_and_test_agent.sh` after all result inventories pass.
""")
    files = sorted(
        path for path in bundle.rglob("*")
        if path.is_file()
        and "results" not in path.parts
        and path.name not in {"SHA256SUMS", "bundle_manifest.json"}
    )
    write(bundle / "SHA256SUMS", "".join(f"{checksum(path)}  {path.relative_to(bundle)}\n" for path in files))
    summary = {
        "patients": 6, "variant_records": 7, "alphagenome_requests": 7,
        "alphagenome_raw_outputs": 11, "boltz_unique_variants": 6,
        "boltz_variant_jobs": 16, "boltz_wt_baselines": 4, "boltz_total_jobs": 20,
        "md_variant_jobs": 5, "md_wt_jobs": 1, "bundle": str(bundle),
    }
    write(bundle / "bundle_manifest.json", json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

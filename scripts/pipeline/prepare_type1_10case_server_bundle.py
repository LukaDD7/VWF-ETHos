#!/usr/bin/env python3
"""Build the exact, portable Type-1 10-case API/GPU run bundle."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def write(path: Path, text: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
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
    parser.add_argument(
        "--input-package",
        type=Path,
        default=ROOT / "outputs/type1_panel_agent_20260828/input_package",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "outputs/type1_panel_agent_20260828/server_bundle",
    )
    args = parser.parse_args()
    bundle = args.output_dir.resolve()
    input_dir = bundle / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "alphagenome_requests.csv", "alphagenome_requests.vcf", "boltz_variants.csv",
        "genetics.csv", "clinical_context.csv", "normalization_audit.csv",
        "prepared_payload.json",
    ):
        shutil.copy2(args.input_package / name, input_dir / name)

    alpha = pd.read_csv(input_dir / "alphagenome_requests.csv")
    boltz = pd.read_csv(input_dir / "boltz_variants.csv")
    if len(alpha) != 9 or len(boltz) != 5:
        raise ValueError(f"Unexpected inventory: AlphaGenome={len(alpha)}, Boltz variants={len(boltz)}")

    boltz_dir = bundle / "boltz"
    subprocess.run(
        [
            "python", str(ROOT / "scripts/pipeline/generate_vwd_functional_boltz2_yamls.py"),
            "--variants-csv", str(input_dir / "boltz_variants.csv"),
            "--output-dir", str(boltz_dir),
            "--minimal-wt-baselines",
        ],
        cwd=ROOT,
        check=True,
    )
    boltz_manifest = pd.read_csv(boltz_dir / "job_manifest.csv")
    if len(boltz_manifest) != 15 or boltz_manifest["run_decision"].eq("WT_BASELINE").sum() != 7:
        raise ValueError("Focused Boltz inventory must be 7 WT baselines + 8 variant jobs")

    pd.DataFrame([
        {
            "case_id": "CASE_T1_008",
            "aa_change": "P1413L",
            "system": "7A6O AIM-A1 closed-state equilibrium MD",
            "matched_wt": "TYPE1_WT_MATCHED",
            "production_ns": 50,
            "priority": "P0",
            "reason": "Only A1 missense in this batch; targeted MD can assess A1/AIM contact integrity after assay-matched Boltz.",
        }
    ]).to_csv(bundle / "md_priority.csv", index=False)

    write(bundle / "run_alphagenome.sh", """#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle"
: "${ALPHAGENOME_API_KEY:?export ALPHAGENOME_API_KEY before running}"
python "$ROOT_DIR/scripts/pipeline/run_type1_10case_alphagenome.py" \
  --requests "$BUNDLE/input/alphagenome_requests.csv" \
  --output-dir "$BUNDLE/results/alphagenome" --expected-requests 9
""", executable=True)

    write(bundle / "run_boltz.sh", """#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle"
GPUS="${GPUS:-4}"
bash "$ROOT_DIR/scripts/pipeline/run_vwd_functional_boltz2_panel.sh" \
  --input-dir "$BUNDLE/boltz/yamls" \
  --out-dir "$BUNDLE/results/boltz/raw" \
  --gpus "$GPUS"
python "$ROOT_DIR/scripts/pipeline/parse_vwd_functional_boltz2_results.py" \
  --results-dir "$BUNDLE/results/boltz/raw" \
  --manifest "$BUNDLE/boltz/job_manifest.csv" \
  --output "$BUNDLE/results/boltz/boltz_results_summary.csv"
""", executable=True)

    write(bundle / "run_md_p1413l.sh", """#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle"
: "${FOLDX:?export FOLDX=/absolute/path/to/foldx}"
python "$ROOT_DIR/scripts/pipeline/build_2b_mutants_foldx.py" \
  --wt "$ROOT_DIR/structures/7A6O_AIM_A1_clean.pdb" \
  --foldx "$FOLDX" --detect-offset --variants P1413L \
  --out-dir "$ROOT_DIR/structures/7a6o_mutants"
cp "$ROOT_DIR/structures/7a6o_mutants/WT_Repair.pdb" \
  "$ROOT_DIR/structures/7a6o_mutants/TYPE1_WT_MATCHED.pdb"
bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh" \
  --pdb "$ROOT_DIR/structures/7a6o_mutants/TYPE1_WT_MATCHED.pdb" --variant TYPE1_WT_MATCHED --skip-vacuum
bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh" \
  --pdb "$ROOT_DIR/structures/7a6o_mutants/P1413L.pdb" --variant P1413L --skip-vacuum
NS=50 bash "$ROOT_DIR/scripts/pipeline/run_7a6o_variant_direct.sh" TYPE1_WT_MATCHED "${WT_GPU:-0}" &
WT_PID=$!
NS=50 bash "$ROOT_DIR/scripts/pipeline/run_7a6o_variant_direct.sh" P1413L "${MUT_GPU:-1}" &
MUT_PID=$!
wait "$WT_PID" "$MUT_PID"
python "$ROOT_DIR/scripts/pipeline/analyze_7a6o_completed_md.py" \
  --input "$ROOT_DIR/output/gromacs_md_autoinhib" \
  --output "$BUNDLE/results/md" \
  --variants TYPE1_WT_MATCHED,P1413L --wt-variant TYPE1_WT_MATCHED
""", executable=True)

    write(bundle / "ingest_and_test_agent.sh", """#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
BUNDLE="$ROOT_DIR/outputs/type1_panel_agent_20260828/server_bundle"
OUT="$ROOT_DIR/outputs/type1_panel_agent_20260828/server_results"
mkdir -p "$OUT"
python "$ROOT_DIR/scripts/pipeline/ingest_type1_10case_results.py" \
  --alphagenome-features "$BUNDLE/results/alphagenome/alphagenome_agent_features.csv" \
  --boltz-summary "$BUNDLE/results/boltz/boltz_results_summary.csv" \
  --boltz-manifest "$BUNDLE/boltz/job_manifest.csv" \
  --boltz-variants "$BUNDLE/input/boltz_variants.csv" \
  --md-features "$BUNDLE/results/md/aim_a1_contacts_summary.csv" \
  --expected-alpha 9 --expected-boltz-jobs 15 --expected-boltz-wt 7 \
  --output "$OUT/type1_computational_evidence.csv"
python "$ROOT_DIR/scripts/run_vwd_langgraph_v0.py" \
  --workbook "$ROOT_DIR/outputs/type1_panel_agent_20260828/type1_panel_agent_package.xlsx" \
  --output-dir "$OUT/agent_run" --computational-panels
""", executable=True)

    write(bundle / "README.md", f"""# Type-1 10-case server bundle

- AlphaGenome: 9 normalized variants x 11 official raw OutputTypes; metadata-selected VWF-relevant tracks; all official recommended scorer views.
- Boltz-2: 5 missense variants, 8 variant-assay jobs, 7 matched WT baselines = 15 jobs total.
- MD: one targeted 7A6O matched-WT/P1413L pair, 50 ns each; not a blanket MD run for all missense variants.
- Excluded from direct AlphaGenome/Boltz input: CASE_T1_004 large deletion/duplication only.
- CASE_T1_007: VWF c.3379+1G>A is included; hemophilia A is retained as a comorbidity and FVIII:C is blocked from type-2N inference.

Run in order:

1. `bash run_alphagenome.sh` (API; needs `ALPHAGENOME_API_KEY`)
2. `GPUS=4 bash run_boltz.sh` (GPU)
3. `FOLDX=/path/to/foldx bash run_md_p1413l.sh` (2 GPUs recommended)
4. Return the `results/` directory and run `bash ingest_and_test_agent.sh`

Expected inventory: AlphaGenome 9/9 successful cases, Boltz 15/15 completed jobs, and MD QC for `TYPE1_WT_MATCHED` plus `P1413L`.
""")

    files = sorted(
        path for path in bundle.rglob("*")
        if path.is_file()
        and "results" not in path.parts
        and path.name not in {"SHA256SUMS", "bundle_manifest.json"}
    )
    write(bundle / "SHA256SUMS", "".join(f"{checksum(path)}  {path.relative_to(bundle)}\n" for path in files))
    summary = {
        "alphagenome_cases": 9,
        "alphagenome_raw_outputs": 11,
        "boltz_variants": 5,
        "boltz_variant_jobs": 8,
        "boltz_wt_baselines": 7,
        "boltz_total_jobs": 15,
        "md_pairs": 1,
        "bundle": str(bundle),
    }
    write(bundle / "bundle_manifest.json", json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Prepare server-ready Boltz and MD input bundles for Type 2M LOF follow-up."""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/type2m_lof_server_inputs_2026-06-28"
EVAL = ROOT / "output/eval_v2"
GENERATOR = ROOT / "scripts/pipeline/generate_vwd_functional_boltz2_yamls.py"

RECOMMENDED_ASSAYS = {
    "a1_gpiba_forced_binding",
    "a1_heparan_sulfate_binding",
    "a1_aim_autoinhibition_context",
    "a3_collagen_binding",
}

A1_GPIB_ANCHORS = [
    {
        "aa_change": "G1324S",
        "protein_pos": 1324,
        "domain": "A1",
        "priority_rank": 0,
        "mechanism_track": "A1_2M_literature_anchor_G1324S",
        "recommended_md_model": "1SQ0_A1_GPIb_complex_MD + A1_dynamics",
        "why": "literature anchor: G1324S can preserve fold but impair A1 collagen/GPIb functional dynamics",
    }
]

HARD_2B_NEGATIVES = [
    {
        "aa_change": "C1272R",
        "protein_pos": 1272,
        "domain": "A1",
        "priority_rank": 3,
        "mechanism_track": "2B_hard_negative_control",
        "recommended_md_model": "1SQ0_A1_GPIb_complex_MD",
        "why": "true 2B hard case previously pulled toward 2M by static LOF axes; guard against overcalling 2M",
    },
    {
        "aa_change": "V1316M",
        "protein_pos": 1316,
        "domain": "A1",
        "priority_rank": 3,
        "mechanism_track": "2B_hard_negative_control",
        "recommended_md_model": "1SQ0_A1_GPIb_complex_MD",
        "why": "true 2B hotspot hard case; guard against overcalling 2M",
    },
]


def parse_aa_change(aa_change: str) -> tuple[str, int, str]:
    return aa_change[0], int(aa_change[1:-1]), aa_change[-1]


def write_variant_list(path: Path, variants: list[str]) -> None:
    path.write_text("\n".join(variants) + "\n")


def write_script(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(0o755)


def prepare_boltz() -> None:
    boltz_dir = OUT / "boltz_missing_2m"
    full_dir = boltz_dir / "_generated_full_panel"
    run_dir = boltz_dir / "run_panel"
    if boltz_dir.exists():
        shutil.rmtree(boltz_dir)
    boltz_dir.mkdir(parents=True)

    q = pd.read_csv(EVAL / "eval_v2_2m_lof_boltz_priority_queue.csv")
    rows = []
    for _, r in q.iterrows():
        wt, pos, mut = parse_aa_change(str(r["aa_change"]))
        rows.append({
            "aa_change": r["aa_change"],
            "wt_aa": wt,
            "position": pos,
            "mut_aa": mut,
            "subtype": "2M",
            "domain": r["domain"],
            "mechanism_track": r["mechanism_track"],
            "recommended_boltz_axis": r["recommended_boltz_axis"],
        })
    variants_csv = boltz_dir / "variants_missing_2m.csv"
    pd.DataFrame(rows).to_csv(variants_csv, index=False)

    subprocess.run(
        [
            "/opt/anaconda3/bin/python", str(GENERATOR),
            "--variants-csv", str(variants_csv),
            "--output-dir", str(full_dir),
            "--write-json-batches",
            "--batch-size", "30",
        ],
        cwd=ROOT,
        check=True,
    )

    run_yaml_dir = run_dir / "yamls"
    run_yaml_dir.mkdir(parents=True)
    manifest = pd.read_csv(full_dir / "job_manifest.csv")
    run_assays = set(manifest.loc[manifest["run_decision"].eq("RUN"), "assay_key"]) & RECOMMENDED_ASSAYS
    keep = manifest[
        (manifest["run_decision"].eq("RUN") & manifest["assay_key"].isin(run_assays))
        | (manifest["run_decision"].eq("WT_BASELINE") & manifest["assay_key"].isin(run_assays))
    ].copy()
    keep.to_csv(run_dir / "job_manifest.csv", index=False)

    diagnostic = pd.read_csv(full_dir / "diagnostic_panel.csv")
    diagnostic[diagnostic["assay_key"].isin(run_assays)].to_csv(run_dir / "diagnostic_panel.csv", index=False)
    shutil.copy(full_dir / "variants_master.csv", run_dir / "variants_master.csv")

    jobs = []
    for _, r in keep.iterrows():
        yaml_rel = Path(str(r["yaml_path"]))
        src = full_dir / yaml_rel
        dst = run_dir / yaml_rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(src, dst)
        jobs.append({
            "job_name": r["job_name"],
            "yaml": str(yaml_rel),
            "assay_key": r["assay_key"],
            "variant_id": r["variant_id"],
            "aa_change": r.get("aa_change", ""),
            "run_decision": r["run_decision"],
        })

    (run_dir / "summary.json").write_text(json.dumps({
        "purpose": "Type 2M LOF missing Boltz top-up panel",
        "n_variants": len(rows),
        "n_yaml_jobs": len(keep),
        "assays": sorted(run_assays),
        "jobs": jobs,
    }, indent=2))

    batch_dir = run_dir / "json_batches"
    batch_dir.mkdir()
    (batch_dir / "type2m_lof_boltz_jobs.json").write_text(json.dumps({
        "name": "type2m_lof_boltz_missing_2m",
        "jobs": jobs,
    }, indent=2))

    write_script(
        boltz_dir / "run_boltz_missing_2m.sh",
        """#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
PKG_DIR="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/boltz_missing_2m"
GPUS="${GPUS:-4}"
GPU_IDS="${GPU_IDS:-}"
mkdir -p "$PKG_DIR/run_panel/boltz_results"
cmd=(bash "$ROOT_DIR/scripts/pipeline/run_vwd_functional_boltz2_panel.sh"
  --input-dir "$PKG_DIR/run_panel/yamls"
  --out-dir "$PKG_DIR/run_panel/boltz_results"
  --gpus "$GPUS")
if [ -n "$GPU_IDS" ]; then
  cmd+=(--gpu-ids "$GPU_IDS")
fi
echo "${cmd[@]}"
"${cmd[@]}"
""",
    )


def prepare_md() -> None:
    md_a1 = OUT / "md_a1_gpiba"
    md_7a6o = OUT / "md_7a6o_closed_state"
    support = OUT / "support_scripts"
    for d in [md_a1, md_7a6o]:
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    if support.exists():
        shutil.rmtree(support)
    support.mkdir(parents=True)
    for script_name in ["run_md_variant_direct.sh", "run_md_resilient.py"]:
        src = ROOT / "scripts/pipeline" / script_name
        dst = support / script_name
        shutil.copy(src, dst)
        if dst.suffix == ".sh":
            dst.chmod(0o755)

    q = pd.read_csv(EVAL / "eval_v2_2m_lof_md_priority_queue.csv")
    q["aa_change"] = q["aa_change"].astype(str)

    extra = pd.DataFrame(A1_GPIB_ANCHORS + HARD_2B_NEGATIVES)
    for col in q.columns:
        if col not in extra.columns:
            extra[col] = pd.NA
    manifest = pd.concat([q, extra[q.columns]], ignore_index=True)
    manifest = manifest.drop_duplicates("aa_change", keep="first").sort_values(["priority_rank", "protein_pos", "aa_change"])
    manifest.to_csv(md_a1 / "a1_gpiba_md_manifest.csv", index=False)

    p0 = manifest[manifest["priority_rank"].eq(0)]["aa_change"].tolist()
    p1 = manifest[manifest["priority_rank"].isin([0, 1])]["aa_change"].tolist()
    controls = manifest[manifest["priority_rank"].isin([2, 3])]["aa_change"].tolist()
    all_rec = manifest["aa_change"].tolist()
    write_variant_list(md_a1 / "a1_gpiba_p0_plus_anchor_variants.txt", p0)
    write_variant_list(md_a1 / "a1_gpiba_p0_p1_variants.txt", p1)
    write_variant_list(md_a1 / "a1_gpiba_controls_variants.txt", controls)
    write_variant_list(md_a1 / "a1_gpiba_all_recommended_variants.txt", all_rec)

    # 7A6O closed-state follow-up: useful as a secondary closed-state/AIM check.
    boltz_q = pd.read_csv(EVAL / "eval_v2_2m_lof_boltz_priority_queue.csv")
    missing_a1 = boltz_q[boltz_q["domain"].eq("A1")]["aa_change"].astype(str).tolist()
    existing_7a6o = {p.name for p in (ROOT / "md_data/7a6o_reference_md/variants").glob("*") if p.is_dir()}
    closed = []
    for v in p1 + missing_a1:
        if v not in existing_7a6o and v not in closed:
            closed.append(v)
    write_variant_list(md_7a6o / "7a6o_2m_closed_state_new_variants.txt", closed)
    pd.DataFrame({"aa_change": closed, "recommended_md_model": "7A6O_AIM_A1_closed_state_MD"}).to_csv(
        md_7a6o / "7a6o_2m_closed_state_manifest.csv", index=False
    )

    write_script(
        md_a1 / "run_a1_gpiba_p0_md.sh",
        """#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
PKG_DIR="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba"
VARIANTS_FILE="${VARIANTS_FILE:-$PKG_DIR/a1_gpiba_p0_plus_anchor_variants.txt}"
GPU_IDS="${GPU_IDS:-4,5,6,7}"
NS="${NS:-50}"
RELAX_JOBS="${RELAX_JOBS:-6}"
RELAX_NTOMP="${RELAX_NTOMP:-8}"
NTOMP="${NTOMP:-8}"
FOLDX_BIN="${FOLDX:-foldx}"
GMX_BIN="${GMX:-$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx}"
PY_BIN="${PY:-python3}"
MUT_DIR="$ROOT_DIR/structures/1sq0_a1_gpiba_mutants"
SCHEDULER="$ROOT_DIR/scripts/pipeline/run_md_resilient.py"
RUNNER="$ROOT_DIR/scripts/pipeline/run_md_variant_direct.sh"
[ -f "$SCHEDULER" ] || SCHEDULER="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/support_scripts/run_md_resilient.py"
[ -f "$RUNNER" ] || RUNNER="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/support_scripts/run_md_variant_direct.sh"

"$PY_BIN" "$ROOT_DIR/scripts/pipeline/build_2b_mutants_foldx.py" \
  --wt "$ROOT_DIR/structures/1SQ0.pdb" \
  --chain A \
  --offset 763 \
  --variants-file "$VARIANTS_FILE" \
  --out-dir "$MUT_DIR" \
  --foldx "$FOLDX_BIN"

mkdir -p "$ROOT_DIR/output/gromacs_md_a1_gpiba"
active=0
while read -r v; do
  v="${v%%#*}"; v="$(echo "$v" | awk '{print $1}')"
  [ -z "$v" ] && continue
  [ -f "$ROOT_DIR/output/gromacs_md_a1_gpiba/$v/relax_pdb/solv_ions_em.gro" ] && continue
  RELAX_NTOMP="$RELAX_NTOMP" GMX="$GMX_BIN" \
    bash "$ROOT_DIR/scripts/pipeline/relax_autoinhib_structure.sh" \
      --pdb "$MUT_DIR/$v.pdb" --variant "$v" --system a1_gpiba --skip-vacuum \
      > "$ROOT_DIR/output/gromacs_md_a1_gpiba/${v}_relax.log" 2>&1 &
  active=$((active + 1))
  if [ "$active" -ge "$RELAX_JOBS" ]; then
    wait -n
    active=$((active - 1))
  fi
done < "$VARIANTS_FILE"
wait

NS="$NS" NTOMP="$NTOMP" "$PY_BIN" "$SCHEDULER" \
  --root "$ROOT_DIR" \
  --system a1_gpiba \
  --md-tag md_a1_gpiba \
  --variants-file "$VARIANTS_FILE" \
  --mutant-dir "$MUT_DIR" \
  --runner "$RUNNER" \
  --gpu-ids "$GPU_IDS" \
  --ns "$NS" \
  --ntomp "$NTOMP"
""",
    )

    write_script(
        md_a1 / "run_a1_gpiba_all_recommended_md.sh",
        """#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
PKG_DIR="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba"
VARIANTS_FILE="$PKG_DIR/a1_gpiba_all_recommended_variants.txt" \
  bash "$PKG_DIR/run_a1_gpiba_p0_md.sh"
""",
    )

    write_script(
        md_7a6o / "run_7a6o_2m_closed_state_md.sh",
        """#!/bin/bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
PKG_DIR="$ROOT_DIR/output/type2m_lof_server_inputs_2026-06-28/md_7a6o_closed_state"
GPU_IDS="${GPU_IDS:-4,5,6,7}"
NS="${NS:-50}"
PY_BIN="${PY:-python3}"
if [ ! -f "$ROOT_DIR/structures/7A6O_AIM_A1_clean.pdb" ]; then
  "$PY_BIN" "$ROOT_DIR/scripts/pipeline/fetch_clean_7a6o.py"
fi
bash "$ROOT_DIR/scripts/pipeline/run_new_2b_saltbridge_md.sh" \
  --stage all \
  --variants-file "$PKG_DIR/7a6o_2m_closed_state_new_variants.txt" \
  --gpu-ids "$GPU_IDS" \
  --ns "$NS"
""",
    )


def write_readme() -> None:
    readme = OUT / "README_RUN_ON_SERVER.md"
    readme.write_text(
        """# Type 2M LOF server input package

Generated on 2026-06-28.

Run from repository root on the server.

## 1. Boltz top-up for missing clean 2M

```bash
export GPUS=4
# optional: export GPU_IDS=4,5,6,7
bash output/type2m_lof_server_inputs_2026-06-28/boltz_missing_2m/run_boltz_missing_2m.sh
```

Input YAMLs:

`output/type2m_lof_server_inputs_2026-06-28/boltz_missing_2m/run_panel/yamls`

This package contains 30 YAML jobs: 16 clean 2M variants projected onto the relevant A1/A3 LOF axes plus 4 WT baselines.

## 2. A1-GPIb complex MD, priority P0

```bash
export FOLDX=/path/to/foldx
export GMX=/path/to/gmx
export GPU_IDS=4,5,6,7
export NS=50
bash output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba/run_a1_gpiba_p0_md.sh
```

Default P0 list:

`output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba/a1_gpiba_p0_plus_anchor_variants.txt`

This uses `structures/1SQ0.pdb`, mutates VWF chain A with `offset=763`, relaxes into `output/gromacs_md_a1_gpiba/<variant>/relax_pdb`, then runs production MD in `md_a1_gpiba`.

## 3. A1-GPIb all recommended MD

After P0 succeeds, run the full recommended A1-GPIb list:

```bash
export FOLDX=/path/to/foldx
export GMX=/path/to/gmx
export GPU_IDS=4,5,6,7
export NS=50
bash output/type2m_lof_server_inputs_2026-06-28/md_a1_gpiba/run_a1_gpiba_all_recommended_md.sh
```

## 4. 7A6O closed-state secondary MD

```bash
export FOLDX=/path/to/foldx
export GMX=/path/to/gmx
export GPU_IDS=4,5,6,7
export NS=50
bash output/type2m_lof_server_inputs_2026-06-28/md_7a6o_closed_state/run_7a6o_2m_closed_state_md.sh
```

This reuses the existing 7A6O AIM-A1 closed-state pipeline and writes into `output/gromacs_md_autoinhib`.

## Key manifests

- `boltz_missing_2m/run_panel/job_manifest.csv`
- `md_a1_gpiba/a1_gpiba_md_manifest.csv`
- `md_7a6o_closed_state/7a6o_2m_closed_state_manifest.csv`
- `support_scripts/` contains fallback copies of the generic MD runner scripts.

## Notes

- P0 A1-GPIb MD targets true 2M cases currently miscalled as 2B plus the G1324S literature anchor.
- 2B hard negatives are included in the all-recommended/control list to protect 2B recall when calibrating LOF thresholds.
- Do not overwrite old Boltz/MD output directories unless intentionally rerunning; these scripts are resumable where the underlying runner is resumable.
"""
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    prepare_boltz()
    prepare_md()
    write_readme()
    print(f"Prepared Type 2M LOF server inputs -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

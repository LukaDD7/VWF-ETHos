#!/usr/bin/env python3
"""Build a VWF multi-module MD structure inventory from local Boltz-2 outputs.

The inventory records, for each VWF functional assay, the best local WT model
that can be used as an MD starting structure after conversion and relaxation.
It does not download anything and does not launch MD.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = REPO_ROOT / "output/boltz2_vwd_functional_panel/boltz_results"
DEFAULT_OUTPUT = REPO_ROOT / "outputs/computational_panel_20260829/md/module_structure_inventory.csv"


ASSAY_SPECS = [
    {
        "assay": "dprime_d3_fviii_binding",
        "module": "D_prime_D3",
        "selection_metric": "iptm",
    },
    {
        "assay": "a1_gpiba_forced_binding",
        "module": "A1_GPIb_alpha",
        "selection_metric": "iptm",
    },
    {
        "assay": "a1_heparan_sulfate_binding",
        "module": "A1_heparan_sulfate",
        "selection_metric": "iptm",
    },
    {
        "assay": "a1_aim_autoinhibition_context",
        "module": "A1_AIM",
        "selection_metric": "ptm",
    },
    {
        "assay": "a2_folded_stability",
        "module": "A2_folded",
        "selection_metric": "ptm",
    },
    {
        "assay": "a2_adamts13_folded_complex",
        "module": "A2_ADAMTS13",
        "selection_metric": "iptm",
    },
    {
        "assay": "vwf73_adamts13_substrate",
        "module": "VWF73_ADAMTS13",
        "selection_metric": "iptm",
    },
    {
        "assay": "a3_collagen_binding",
        "module": "A3_collagen",
        "selection_metric": "iptm",
    },
    {
        "assay": "c1_collagen_binding",
        "module": "C1_collagen",
        "selection_metric": "iptm",
    },
    {
        "assay": "c2_collagen_binding",
        "module": "C2_collagen",
        "selection_metric": "iptm",
    },
    {
        "assay": "c4_integrin_binding",
        "module": "C4_integrin",
        "selection_metric": "iptm",
    },
    {
        "assay": "d1d2_propeptide_context",
        "module": "D1_D2_propeptide",
        "selection_metric": "ptm",
    },
    {
        "assay": "d4_assembly_context",
        "module": "D4_assembly",
        "selection_metric": "ptm",
    },
    {
        "assay": "c_domain_assembly_context",
        "module": "C1_C6_assembly",
        "selection_metric": "ptm",
    },
    {
        "assay": "ck_dimerization_context",
        "module": "CK_dimerization",
        "selection_metric": "ptm",
    },
]


MODEL_RE = re.compile(r"_model_(\d+)\.json$")


def load_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def build_inventory(results_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for spec in ASSAY_SPECS:
        assay = spec["assay"]
        module = spec["module"]
        selection_metric = spec["selection_metric"]
        job_dir = results_dir / f"boltz_results_VWF_WT__{assay}"

        row = {
            "assay": assay,
            "module": module,
            "structure_source": "boltz2_wt_model",
            "selection_metric": selection_metric,
            "selected_model": "",
            "selected_metric_value": "",
            "n_models": "",
            "avg_metric_value": "",
            "confidence_tier": "",
            "cif_path": "",
            "status": "missing",
            "md_ready": "no",
            "md_status": "not_started",
        }

        if not job_dir.is_dir():
            rows.append(row)
            continue

        confidence_files = sorted(
            job_dir.glob("predictions/*/confidence_*_model_*.json")
        )
        if not confidence_files:
            rows.append(row)
            continue

        candidates: list[tuple[float, Path]] = []
        for confidence_path in confidence_files:
            try:
                payload = load_json(confidence_path)
            except (OSError, json.JSONDecodeError):
                continue

            if selection_metric == "iptm":
                value = float(payload.get("iptm") or 0.0)
            else:
                value = float(payload.get("ptm") or payload.get("confidence_score") or 0.0)
            candidates.append((value, confidence_path))

        if not candidates:
            rows.append(row)
            continue

        metric_value, confidence_path = max(candidates, key=lambda item: item[0])
        metric_values = [value for value, _ in candidates]
        avg_metric_value = sum(metric_values) / len(metric_values)
        if metric_value >= 0.70:
            confidence_tier = "high"
        elif metric_value >= 0.50:
            confidence_tier = "moderate"
        else:
            confidence_tier = "low"
        match = MODEL_RE.search(confidence_path.name)
        model_id = match.group(1) if match else ""
        cif_path = confidence_path.with_name(
            confidence_path.name.replace("confidence_", "").replace(".json", ".cif")
        )

        row.update(
            {
                "selected_model": model_id,
                "selected_metric_value": f"{metric_value:.6f}",
                "n_models": str(len(candidates)),
                "avg_metric_value": f"{avg_metric_value:.6f}",
                "confidence_tier": confidence_tier,
                "cif_path": str(cif_path.relative_to(REPO_ROOT)) if cif_path.is_relative_to(REPO_ROOT) else str(cif_path),
                "status": "available" if cif_path.is_file() else "cif_missing",
                "md_ready": "needs_cif_to_pdb_and_relaxation",
                "md_status": "not_started",
            }
        )
        rows.append(row)

    return rows


def write_inventory(rows: list[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Directory containing boltz_results_* directories.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output CSV path.",
    )
    args = parser.parse_args()

    rows = build_inventory(args.results_dir.resolve())
    write_inventory(rows, args.output.resolve())

    available = sum(row["status"] == "available" for row in rows)
    print(f"Wrote {len(rows)} assay rows to {args.output}")
    print(f"Available WT Boltz models: {available}/{len(rows)}")


if __name__ == "__main__":
    main()

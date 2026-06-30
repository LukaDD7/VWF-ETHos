#!/usr/bin/env python3
"""
Extract interface-level features from Boltz A1-GPIb complex predictions.

The existing A1-GPIb feature mostly uses a global complex confidence proxy
such as iPTM. This script squeezes more information out of the predicted
complex structures:

- A1-GPIb heavy-atom contact counts.
- Interface residue counts on A1 and GPIb.
- Minimum / percentile inter-chain distances.
- Interface PAE summaries.
- Interface-contact pLDDT summaries.
- WT-relative deltas and assay-internal z-scores.

No subtype label is needed to compute features. If an Eval table is provided,
labels are joined only for readout/calibration.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shlex
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


AA_ATOMS_BACKBONE = {"N", "CA", "C", "O", "OXT"}
CONTACT_CUTOFF_A = 5.0
NEAR_CUTOFF_A = 8.0


@dataclass
class Atom:
    chain: str
    seq_id: int
    atom_name: str
    element: str
    x: float
    y: float
    z: float
    b_factor: float


def parse_atom_site(cif_path: Path) -> list[Atom]:
    """Parse the _atom_site loop from a Boltz mmCIF file."""
    lines = cif_path.read_text().splitlines()
    atoms: list[Atom] = []
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        j = i + 1
        headers = []
        while j < len(lines) and lines[j].startswith("_"):
            headers.append(lines[j].strip())
            j += 1
        if not headers or not any(h.startswith("_atom_site.") for h in headers):
            i = j
            continue
        if "_atom_site.group_PDB" not in headers:
            i = j
            continue

        idx = {h: k for k, h in enumerate(headers)}
        required = [
            "_atom_site.group_PDB",
            "_atom_site.type_symbol",
            "_atom_site.label_atom_id",
            "_atom_site.label_seq_id",
            "_atom_site.label_asym_id",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.B_iso_or_equiv",
        ]
        missing = [h for h in required if h not in idx]
        if missing:
            raise ValueError(f"{cif_path}: atom_site missing {missing}")

        while j < len(lines):
            raw = lines[j].strip()
            if not raw:
                j += 1
                continue
            if raw == "#" or raw == "loop_" or raw.startswith("_"):
                break
            parts = shlex.split(raw)
            if len(parts) < len(headers) or parts[idx["_atom_site.group_PDB"]] != "ATOM":
                j += 1
                continue
            try:
                atoms.append(
                    Atom(
                        chain=parts[idx["_atom_site.label_asym_id"]],
                        seq_id=int(parts[idx["_atom_site.label_seq_id"]]),
                        atom_name=parts[idx["_atom_site.label_atom_id"]],
                        element=parts[idx["_atom_site.type_symbol"]],
                        x=float(parts[idx["_atom_site.Cartn_x"]]),
                        y=float(parts[idx["_atom_site.Cartn_y"]]),
                        z=float(parts[idx["_atom_site.Cartn_z"]]),
                        b_factor=float(parts[idx["_atom_site.B_iso_or_equiv"]]),
                    )
                )
            except Exception as exc:
                raise ValueError(f"{cif_path}: failed parsing atom_site row {raw!r}") from exc
            j += 1
        i = j
    return atoms


def chain_residue_counts(atoms: Iterable[Atom]) -> dict[str, int]:
    residues: dict[str, set[int]] = defaultdict(set)
    for atom in atoms:
        residues[atom.chain].add(atom.seq_id)
    return {chain: len(ids) for chain, ids in residues.items()}


def load_npz_array(path: Path, key: str) -> np.ndarray | None:
    if not path.exists():
        return None
    arr = np.load(path)
    if key not in arr.files:
        return None
    return arr[key]


def interface_features_for_model(cif_path: Path, pae_path: Path, plddt_path: Path) -> dict:
    atoms = parse_atom_site(cif_path)
    chains = sorted(chain_residue_counts(atoms))
    if len(chains) < 2:
        raise ValueError(f"{cif_path}: expected at least two chains, got {chains}")
    a_chain, b_chain = chains[0], chains[1]

    # Heavy atoms only. The Boltz CIFs have no hydrogens in practice, but keep it explicit.
    a_atoms = [a for a in atoms if a.chain == a_chain and a.element.upper() != "H"]
    b_atoms = [a for a in atoms if a.chain == b_chain and a.element.upper() != "H"]
    a_side = [a for a in a_atoms if a.atom_name not in AA_ATOMS_BACKBONE]
    b_side = [a for a in b_atoms if a.atom_name not in AA_ATOMS_BACKBONE]

    def coords(atom_list: list[Atom]) -> np.ndarray:
        return np.asarray([(a.x, a.y, a.z) for a in atom_list], dtype=float)

    ac = coords(a_atoms)
    bc = coords(b_atoms)
    if len(ac) == 0 or len(bc) == 0:
        raise ValueError(f"{cif_path}: missing atoms for chains {a_chain}/{b_chain}")

    # Distance matrix is modest for these complexes, so direct numpy is fine.
    d = np.linalg.norm(ac[:, None, :] - bc[None, :, :], axis=2)
    contact_mask = d <= CONTACT_CUTOFF_A
    near_mask = d <= NEAR_CUTOFF_A
    contact_pairs = int(contact_mask.sum())
    near_pairs = int(near_mask.sum())
    min_dist = float(d.min())
    p05_dist = float(np.percentile(d, 5))
    p10_dist = float(np.percentile(d, 10))

    a_iface_atom_idx = np.where(contact_mask.any(axis=1))[0]
    b_iface_atom_idx = np.where(contact_mask.any(axis=0))[0]
    a_iface_res = sorted({a_atoms[i].seq_id for i in a_iface_atom_idx})
    b_iface_res = sorted({b_atoms[i].seq_id for i in b_iface_atom_idx})

    # Side-chain contact count is more sensitive to interface chemistry than backbone proximity.
    side_contact_pairs = np.nan
    if a_side and b_side:
        sd = np.linalg.norm(coords(a_side)[:, None, :] - coords(b_side)[None, :, :], axis=2)
        side_contact_pairs = int((sd <= CONTACT_CUTOFF_A).sum())

    chain_counts = chain_residue_counts(atoms)
    a_len = chain_counts[a_chain]
    b_len = chain_counts[b_chain]
    pae = load_npz_array(pae_path, "pae")
    plddt = load_npz_array(plddt_path, "plddt")

    out = {
        "chain_a": a_chain,
        "chain_b": b_chain,
        "chain_a_residues": a_len,
        "chain_b_residues": b_len,
        "contact_pairs_5a": contact_pairs,
        "near_pairs_8a": near_pairs,
        "sidechain_contact_pairs_5a": side_contact_pairs,
        "a1_interface_residue_count": len(a_iface_res),
        "gpib_interface_residue_count": len(b_iface_res),
        "min_interchain_distance_a": min_dist,
        "p05_interchain_distance_a": p05_dist,
        "p10_interchain_distance_a": p10_dist,
        "a1_interface_residues": ";".join(map(str, a_iface_res)),
        "gpib_interface_residues": ";".join(map(str, b_iface_res)),
    }

    if pae is not None and pae.ndim == 2 and pae.shape[0] >= a_len + b_len and pae.shape[1] >= a_len + b_len:
        ab = pae[:a_len, a_len : a_len + b_len]
        ba = pae[a_len : a_len + b_len, :a_len]
        iface = (ab + ba.T) / 2
        out.update(
            {
                "interface_pae_mean": float(iface.mean()),
                "interface_pae_median": float(np.median(iface)),
                "interface_pae_p10": float(np.percentile(iface, 10)),
                "interface_pae_contact_mean": (
                    float(iface[np.ix_([r - 1 for r in a_iface_res], [r - 1 for r in b_iface_res])].mean())
                    if a_iface_res and b_iface_res
                    else np.nan
                ),
            }
        )
    else:
        out.update(
            {
                "interface_pae_mean": np.nan,
                "interface_pae_median": np.nan,
                "interface_pae_p10": np.nan,
                "interface_pae_contact_mean": np.nan,
            }
        )

    if plddt is not None and plddt.ndim == 1 and plddt.shape[0] >= a_len + b_len:
        a_plddt = plddt[:a_len]
        b_plddt = plddt[a_len : a_len + b_len]
        out.update(
            {
                "a1_plddt_mean": float(a_plddt.mean()),
                "gpib_plddt_mean": float(b_plddt.mean()),
                "a1_interface_plddt_mean": float(a_plddt[[r - 1 for r in a_iface_res]].mean()) if a_iface_res else np.nan,
                "gpib_interface_plddt_mean": float(b_plddt[[r - 1 for r in b_iface_res]].mean()) if b_iface_res else np.nan,
            }
        )
    else:
        out.update(
            {
                "a1_plddt_mean": np.nan,
                "gpib_plddt_mean": np.nan,
                "a1_interface_plddt_mean": np.nan,
                "gpib_interface_plddt_mean": np.nan,
            }
        )

    return out


def find_prediction_dirs(results_dir: Path) -> list[Path]:
    return sorted(p for p in results_dir.glob("**/predictions/*") if p.is_dir())


def job_from_pred_dir(pred_dir: Path) -> str:
    return pred_dir.name


def variant_from_job(job_name: str) -> str:
    if job_name.endswith("_vs_GPIb_alpha"):
        return job_name.replace("_vs_GPIb_alpha", "")
    return job_name.split("__", 1)[0]


def assay_from_job(job_name: str) -> str:
    if job_name.endswith("_vs_GPIb_alpha"):
        return "a1_gpiba_forced_binding"
    return job_name.split("__", 1)[1] if "__" in job_name else ""


def model_index(path: Path) -> int | None:
    m = re.search(r"_model_(\d+)\.", path.name)
    return int(m.group(1)) if m else None


def read_labels(eval_csv: Path | None) -> dict[str, str]:
    if eval_csv is None or not eval_csv.exists():
        return {}
    labels = {}
    with eval_csv.open() as handle:
        for row in csv.DictReader(handle):
            aa = row.get("aa_change") or row.get("variant_id")
            lab = row.get("true_label") or row.get("type2_subtype")
            if aa and lab:
                labels[f"VWF_{aa}" if not aa.startswith("VWF_") else aa] = lab
    return labels


def mean_or_nan(values):
    vals = [float(v) for v in values if v is not None and not (isinstance(v, float) and math.isnan(v))]
    return float(np.mean(vals)) if vals else np.nan


def summarize(records: list[dict], labels: dict[str, str]) -> list[dict]:
    by_job: dict[str, list[dict]] = defaultdict(list)
    for rec in records:
        by_job[rec["job_name"]].append(rec)

    out = []
    numeric_keys = [
        k
        for k, v in records[0].items()
        if isinstance(v, (int, float, np.integer, np.floating)) and k not in {"model_idx"}
    ]
    for job, rows in sorted(by_job.items()):
        variant = variant_from_job(job)
        assay = assay_from_job(job)
        rec = {
            "job_name": job,
            "variant_id": variant,
            "aa_change": variant.replace("VWF_", ""),
            "assay_key": assay,
            "true_label": labels.get(variant, ""),
            "n_models": len(rows),
        }
        for k in numeric_keys:
            rec[k] = mean_or_nan([r.get(k) for r in rows])
        # Keep model-0 residue list as a compact provenance/debug field.
        rows0 = sorted(rows, key=lambda r: r.get("model_idx", 0))
        rec["a1_interface_residues_model0"] = rows0[0].get("a1_interface_residues", "")
        rec["gpib_interface_residues_model0"] = rows0[0].get("gpib_interface_residues", "")
        out.append(rec)

    wt = next((r for r in out if r["variant_id"] == "VWF_WT"), None)
    if wt:
        delta_cols = [
            "contact_pairs_5a",
            "near_pairs_8a",
            "sidechain_contact_pairs_5a",
            "a1_interface_residue_count",
            "gpib_interface_residue_count",
            "min_interchain_distance_a",
            "p05_interchain_distance_a",
            "interface_pae_mean",
            "interface_pae_contact_mean",
            "a1_interface_plddt_mean",
            "gpib_interface_plddt_mean",
        ]
        for rec in out:
            for col in delta_cols:
                rec[f"delta_{col}_vs_wt"] = rec.get(col, np.nan) - wt.get(col, np.nan)

    # Internal z-scores: higher retained score should mean more contacts and lower PAE/distance.
    z_specs = {
        "contact_pairs_5a": 1,
        "sidechain_contact_pairs_5a": 1,
        "a1_interface_residue_count": 1,
        "gpib_interface_residue_count": 1,
        "min_interchain_distance_a": -1,
        "p05_interchain_distance_a": -1,
        "interface_pae_mean": -1,
        "interface_pae_contact_mean": -1,
        "a1_interface_plddt_mean": 1,
        "gpib_interface_plddt_mean": 1,
    }
    for col, sign in z_specs.items():
        vals = np.asarray([r.get(col, np.nan) for r in out if r["variant_id"] != "VWF_WT"], dtype=float)
        mu = np.nanmean(vals)
        sd = np.nanstd(vals)
        for rec in out:
            val = rec.get(col, np.nan)
            rec[f"{col}_retained_z"] = sign * (val - mu) / sd if sd and not np.isnan(val) else np.nan
    for rec in out:
        parts = [
            rec.get("contact_pairs_5a_retained_z", np.nan),
            rec.get("sidechain_contact_pairs_5a_retained_z", np.nan),
            rec.get("a1_interface_residue_count_retained_z", np.nan),
            rec.get("interface_pae_contact_mean_retained_z", np.nan),
            rec.get("a1_interface_plddt_mean_retained_z", np.nan),
        ]
        rec["a1_gpib_interface_retained_z"] = mean_or_nan(parts)
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="output/boltz2_vwd_functional_panel/boltz_results")
    parser.add_argument("--assay-key", default="a1_gpiba_forced_binding")
    parser.add_argument("--eval-csv", default="output/type2m_lof_md_fast_validation_2026-06-29/eval_with_fast_md/eval_v2_predictions_slim.csv")
    parser.add_argument("--per-model-output", default="output/boltz2_vwd_functional_panel/analysis/a1_gpiba_interface_features_per_model.csv")
    parser.add_argument("--summary-output", default="output/boltz2_vwd_functional_panel/analysis/a1_gpiba_interface_features_summary.csv")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    labels = read_labels(Path(args.eval_csv) if args.eval_csv else None)
    rows = []
    pred_dirs = [p for p in find_prediction_dirs(results_dir) if assay_from_job(job_from_pred_dir(p)) == args.assay_key]
    for pred_dir in pred_dirs:
        job = job_from_pred_dir(pred_dir)
        cif_files = sorted(pred_dir.glob(f"{job}_model_*.cif"))
        for cif in cif_files:
            mid = model_index(cif)
            if mid is None:
                continue
            pae = pred_dir / f"pae_{job}_model_{mid}.npz"
            plddt = pred_dir / f"plddt_{job}_model_{mid}.npz"
            try:
                feats = interface_features_for_model(cif, pae, plddt)
            except Exception as exc:
                print(f"[WARN] {job} model {mid}: {exc}")
                continue
            feats.update(
                {
                    "job_name": job,
                    "variant_id": variant_from_job(job),
                    "aa_change": variant_from_job(job).replace("VWF_", ""),
                    "assay_key": assay_from_job(job),
                    "true_label": labels.get(variant_from_job(job), ""),
                    "model_idx": mid,
                    "cif_path": str(cif),
                }
            )
            rows.append(feats)

    if not rows:
        raise SystemExit(f"No {args.assay_key} interface records found under {results_dir}")

    summary = summarize(rows, labels)
    write_csv(Path(args.per_model_output), rows)
    write_csv(Path(args.summary_output), summary)

    print(f"per-model: {args.per_model_output} ({len(rows)} rows)")
    print(f"summary:   {args.summary_output} ({len(summary)} rows)")

    labeled = [r for r in summary if r.get("true_label") in {"2B", "2M"}]
    if labeled:
        print("\nA1-GPIb interface retained z by true label:")
        for lab in ["2B", "2M"]:
            vals = [r["a1_gpib_interface_retained_z"] for r in labeled if r["true_label"] == lab]
            vals = [v for v in vals if not np.isnan(v)]
            if vals:
                print(f"  {lab}: n={len(vals)} mean={np.mean(vals):+.3f} median={np.median(vals):+.3f}")


if __name__ == "__main__":
    main()

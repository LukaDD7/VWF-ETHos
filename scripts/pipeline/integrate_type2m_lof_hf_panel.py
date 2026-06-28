#!/usr/bin/env python3
"""Integrate the HF Type2M LOF supplementary Boltz panel into evidence_matrix.

The Hugging Face panel was run after the main Boltz functional panel, so its
raw summary is long-form by job. This script converts those supplemental jobs
into evidence_matrix-compatible rows and appends them to the existing matrix.

Important calibration choice:
  - Delta-vs-WT uses the WT from the supplementary HF run.
  - z-score uses the original evidence_matrix assay distribution, so the
    classifier thresholds keep the same meaning.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}

ASSAY_PRIMARY_SOURCE = {
    "a1_aim_autoinhibition_context": "avg_ptm",
    "a1_gpiba_forced_binding": "avg_iptm",
    "a1_heparan_sulfate_binding": "avg_iptm",
    "a3_collagen_binding": "avg_iptm",
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--base",
        default=str(ROOT / "output/boltz2_vwd_functional_panel/evidence_matrix.csv"),
        help="Original evidence_matrix.csv",
    )
    ap.add_argument(
        "--hf-summary",
        default=str(
            ROOT
            / "output/hf_type2m_lof_panel/type2m_lof_panel/analysis/boltz_results_summary_annotated.csv"
        ),
        help="Annotated long-form HF supplementary summary",
    )
    ap.add_argument(
        "--out",
        default=str(
            ROOT
            / "output/hf_type2m_lof_panel/type2m_lof_panel/analysis/"
            / "evidence_matrix_with_type2m_lof_hf.csv"
        ),
        help="Combined evidence matrix output",
    )
    ap.add_argument(
        "--supplement-rows-out",
        default=str(
            ROOT
            / "output/hf_type2m_lof_panel/type2m_lof_panel/analysis/"
            / "type2m_lof_hf_evidence_rows.csv"
        ),
        help="Supplement-only wide rows output",
    )
    return ap.parse_args()


def assay_stats(base: pd.DataFrame, assay: str) -> tuple[float, float, pd.Series]:
    col = f"{assay}__primary_value"
    vals = pd.to_numeric(base[col], errors="coerce").dropna()
    if vals.empty:
        return np.nan, np.nan, vals
    return float(vals.mean()), float(vals.std(ddof=1)), vals


def percentile_against_base(primary: float, vals: pd.Series) -> float:
    if vals.empty or pd.isna(primary):
        return np.nan
    return float((vals <= primary).mean())


def delta_direction(delta: float) -> str:
    if pd.isna(delta):
        return ""
    if delta > 0:
        return "higher_than_wt"
    if delta < 0:
        return "lower_than_wt"
    return "same_as_wt"


def make_base_row(base_columns: list[str], rec: pd.Series) -> dict:
    aa_change = str(rec["aa_change"])
    wt = str(rec["wt_aa"])
    mut = str(rec["mut_aa"])
    pos = int(float(rec["position"]))

    row = {col: np.nan for col in base_columns}
    row.update(
        {
            "variant_id": str(rec["variant_id"]),
            "wt_aa": wt,
            "position": pos,
            "mut_aa": mut,
            "aa_change_3letter": f"{AA3.get(wt, wt)}{pos}{AA3.get(mut, mut)}",
            "source_labels": "Type2M",
            "source_files": "HF:lucachangretta/VWF/type2m_lof_panel",
            "source_rows": int(len(str(rec.get("source_labels", "2M")).split("|"))),
            "source_domains": str(rec["inferred_domain"]),
            "inferred_domain": str(rec["inferred_domain"]),
            "wt_validation_status": "ok",
            "domain_source_agrees_with_inferred": True,
            "is_Type1": False,
            "is_Type2A": False,
            "is_Type2B": False,
            "is_Type2M": True,
            "is_Type2N": False,
            "is_Type3": False,
            "is_Negative_GeneBe": False,
            "is_Negative_ClinVar": False,
            "positive_any": True,
            "type2_any": True,
            "negative_any": False,
            "negative_only": False,
            "positive_only": True,
            "negative_positive_overlap": False,
            "n_source_labels": 1,
            "label_policy": "positive_only",
        }
    )
    return row


def fill_assay(row: dict, rec: pd.Series, stats: dict[str, tuple[float, float, pd.Series]]) -> None:
    assay = str(rec["assay_key"])
    prefix = f"{assay}__"
    primary_source = ASSAY_PRIMARY_SOURCE.get(assay, "avg_iptm")
    primary = float(rec["primary_value_recalc"])
    wt_primary = float(rec["wt_primary_value_recalc"])
    delta = primary - wt_primary
    mean, std, base_vals = stats[assay]
    z = (primary - mean) / std if std and not pd.isna(std) else np.nan

    values = {
        "run_decision": "RUN",
        "n_samples": int(rec["n_samples"]),
        "primary_metric": "ptm" if primary_source == "avg_ptm" else "iptm",
        "primary_value": primary,
        "primary_value_source": primary_source,
        "metric_status": "available",
        "wt_primary_value": wt_primary,
        "primary_value_delta_vs_wt": delta,
        "primary_delta_direction": delta_direction(delta),
        "primary_zscore_within_assay": z,
        "primary_percentile_within_assay": percentile_against_base(primary, base_vals),
        "avg_iptm": float(rec["avg_iptm"]),
        "best_iptm": float(rec["best_iptm"]),
        "avg_ptm": float(rec["avg_ptm"]),
        "best_ptm": float(rec["best_ptm"]),
        "avg_complex_plddt": float(rec["avg_complex_plddt"]),
        "best_complex_plddt": float(rec["best_complex_plddt"]),
    }
    for suffix, value in values.items():
        col = f"{prefix}{suffix}"
        if col in row:
            row[col] = value


def main() -> int:
    args = parse_args()
    base = pd.read_csv(args.base)
    hf = pd.read_csv(args.hf_summary)
    hf = hf[
        hf["variant_id"].notna()
        & hf["assay_key"].notna()
        & hf["aa_change"].notna()
        & hf["position"].notna()
        & hf["variant_id"].ne("VWF_WT")
    ].copy()

    overlap = sorted(set(base["variant_id"].astype(str)) & set(hf["variant_id"].astype(str)))
    if overlap:
        raise SystemExit(f"Refusing to append duplicate variant_id rows: {overlap}")

    stats = {assay: assay_stats(base, assay) for assay in ASSAY_PRIMARY_SOURCE}

    rows = []
    for _, group in hf.groupby("variant_id", sort=True):
        rec0 = group.iloc[0]
        row = make_base_row(list(base.columns), rec0)
        for _, rec in group.iterrows():
            fill_assay(row, rec, stats)
        rows.append(row)

    supp = pd.DataFrame(rows, columns=base.columns)
    combined = pd.concat([base, supp], ignore_index=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    supp_out = Path(args.supplement_rows_out)
    supp_out.parent.mkdir(parents=True, exist_ok=True)
    supp.to_csv(supp_out, index=False)
    combined.to_csv(out, index=False)

    print(f"base_rows={len(base)} supplement_rows={len(supp)} combined_rows={len(combined)}")
    print(f"wrote {out}")
    print(f"wrote {supp_out}")
    print(
        supp[
            [
                "variant_id",
                "position",
                "inferred_domain",
                "a1_gpiba_forced_binding__primary_zscore_within_assay",
                "a1_heparan_sulfate_binding__primary_zscore_within_assay",
                "a3_collagen_binding__primary_zscore_within_assay",
            ]
        ].to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

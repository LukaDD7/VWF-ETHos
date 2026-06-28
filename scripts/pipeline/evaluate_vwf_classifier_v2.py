#!/usr/bin/env python3
"""Evaluate the current VWF classifier on the expanded Boltz feature panel.

This is a reporting/evaluation script, not a training script. It intentionally
keeps true labels out of the classifier input and uses them only after prediction
to compute recall/confusion tables.

Inputs:
  output/boltz2_vwd_functional_panel/evidence_matrix.csv
  output/labeled_variants_all.csv
  output/md_7a6o_saltbridge_features.csv   (optional, existing MD subset)
  output/md_7a6o_features.csv              (optional, existing MD subset)

Outputs under output/eval_v2/:
  eval_v2_predictions.csv
  eval_v2_type2_confusion_no_md.csv
  eval_v2_type2_confusion_with_md.csv
  eval_v2_summary.csv
  eval_v2_a1_2b2m_summary.csv
  eval_v2_md_priority_queue.csv
  eval_v2_2m_lof_boltz_priority_queue.csv
  eval_v2_2m_lof_md_priority_queue.csv
  eval_v2_inventory.csv
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TYPE2 = ["2A", "2B", "2M", "2N"]
SUPPORTED_LABELS = set(["1", *TYPE2])

FB_Z = "a1_gpiba_forced_binding__primary_zscore_within_assay"
HEP_Z = "a1_heparan_sulfate_binding__primary_zscore_within_assay"
AIM_STATIC_Z = "a1_aim_autoinhibition_context__primary_zscore_within_assay"
A3_COLLAGEN_Z = "a3_collagen_binding__primary_zscore_within_assay"

EXPECTED_NA_COLS = [
    "ag_rna_delta", "ag_splice_delta", "ag_delta_score", "af3_plddt_mean",
    "af3_plddt_min", "af3_pae_interface", "foldx_ddg_bind", "boltz2_iptm",
    "boltz2_delta_iptm", "aim_release_score", "aim_sb_retained_z",
    "md_face_destab_score",
]


def load_classifier(path: Path):
    spec = importlib.util.spec_from_file_location("agentic_vwf_classifier", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def infer_classifier_domain(pos: int) -> str:
    """Map a native VWF protein position to the domain names expected by RULEs."""
    if 1 <= pos <= 233:
        return "D1"
    if 234 <= pos <= 480:
        return "D2"
    if 481 <= pos <= 728:
        return "D3"
    if 729 <= pos <= 763:
        return "D'"
    if 764 <= pos <= 1233:
        return "D3_extended"
    if 1271 <= pos <= 1492:
        return "A1"
    if 1493 <= pos <= 1684:
        return "A2"
    if 1685 <= pos <= 1873:
        return "A3"
    if 1874 <= pos <= 2255:
        return "D4"
    if 2256 <= pos <= 2299:
        return "C1"
    if 2300 <= pos <= 2363:
        return "C2"
    if 2364 <= pos <= 2411:
        return "C3"
    if 2412 <= pos <= 2450:
        return "C4"
    if 2451 <= pos <= 2459:
        return "C5"
    if 2460 <= pos <= 2527:
        return "C6"
    if 2528 <= pos <= 2813:
        return "CK"
    return ""


def add_inventory_rows(rows: list[dict], name: str, df: pd.DataFrame, label_col: str | None = None):
    rows.append({"dataset": name, "metric": "rows", "value": len(df)})
    if label_col and label_col in df.columns:
        for label, n in df[label_col].value_counts(dropna=False).sort_index().items():
            rows.append({"dataset": name, "metric": f"{label_col}={label}", "value": int(n)})


def boolish(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s.fillna(False)
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def load_eval_table(args, mod) -> tuple[pd.DataFrame, pd.DataFrame]:
    evidence = pd.read_csv(args.evidence)
    labels = pd.read_csv(args.labels)
    labels = labels[(~boolish(labels["conflict"])) & labels["labels"].isin(SUPPORTED_LABELS)].copy()

    evidence = evidence.copy()
    evidence["aa_change"] = (
        evidence["wt_aa"].astype(str)
        + evidence["position"].astype("Int64").astype(str)
        + evidence["mut_aa"].astype(str)
    )
    labels = labels.rename(columns={"labels": "true_label"})

    joined = labels.merge(evidence, on="aa_change", how="left", suffixes=("_label", ""))
    joined["has_boltz_panel"] = joined["variant_id"].notna()
    missing_boltz = joined[~joined["has_boltz_panel"]].copy()
    joined = joined[joined["has_boltz_panel"]].copy()

    joined["protein_pos"] = pd.to_numeric(joined["position_label"], errors="coerce").fillna(joined["position"]).astype(int)
    joined["ref_aa"] = joined["wt_aa_label"].where(joined["wt_aa_label"].notna(), joined["wt_aa"])
    joined["alt_aa"] = joined["mut_aa_label"].where(joined["mut_aa_label"].notna(), joined["mut_aa"])
    joined["domain"] = joined["protein_pos"].apply(infer_classifier_domain)
    joined["type2_subtype"] = joined["true_label"].where(joined["true_label"].isin(TYPE2), np.nan)

    joined["fb_binding_zscore"] = pd.to_numeric(joined.get(FB_Z), errors="coerce")
    joined["heparan_zscore"] = pd.to_numeric(joined.get(HEP_Z), errors="coerce")
    joined["aim_static_zscore"] = pd.to_numeric(joined.get(AIM_STATIC_Z), errors="coerce")
    joined["a3_collagen_zscore"] = pd.to_numeric(joined.get(A3_COLLAGEN_Z), errors="coerce")

    for col in EXPECTED_NA_COLS:
        if col not in joined.columns:
            joined[col] = np.nan

    # Deliberately do not map a1_aim_autoinhibition_context z-score into
    # aim_release_score. That static panel axis was previously shown to be weak
    # and is not equivalent to the geometric AIM release feature expected by RULE6.
    joined["aim_release_score"] = np.nan

    if Path(args.saltbridge).exists():
        sb = pd.read_csv(args.saltbridge)[["variant", "aim_sb_retained_z"]]
        joined = joined.drop(columns=["aim_sb_retained_z"], errors="ignore")
        joined = joined.merge(sb, left_on="aa_change", right_on="variant", how="left").drop(columns=["variant"])
    if Path(args.face_md).exists():
        face = pd.read_csv(args.face_md)[["variant", "md_face_destab_score"]]
        joined = joined.drop(columns=["md_face_destab_score"], errors="ignore")
        joined = joined.merge(face, left_on="aa_change", right_on="variant", how="left").drop(columns=["variant"])

    joined["hotspot"] = joined["protein_pos"].isin(mod.TWO_B_HOTSPOT_POS)
    joined["in_a1_7a6o_range"] = joined["protein_pos"].between(1262, 1466)
    return joined, missing_boltz


def run_predictions(df: pd.DataFrame, mod, use_md: bool) -> pd.DataFrame:
    features = df.copy()
    if not use_md:
        features["aim_sb_retained_z"] = np.nan
        features["md_face_destab_score"] = np.nan
    res = mod.AgenticVWFClassifier().classify_batch(features)
    out = df[["aa_change", "true_label"]].copy()
    suffix = "with_md" if use_md else "no_md"
    out[f"pred_{suffix}"] = res["main_subtype"].values
    out[f"confidence_{suffix}"] = res["confidence"].values
    out[f"reasoning_{suffix}"] = res["reasoning"].values
    return out


def confusion(df: pd.DataFrame, pred_col: str, labels: list[str]) -> pd.DataFrame:
    sub = df[df["true_label"].isin(labels)].copy()
    tab = pd.crosstab(sub["true_label"], sub[pred_col], dropna=False)
    for col in labels + ["1", "uncertain"]:
        if col not in tab.columns:
            tab[col] = 0
    tab = tab[[c for c in labels + ["1", "uncertain"] if c in tab.columns]]
    return tab.reset_index()


def recall_rows(df: pd.DataFrame, pred_col: str, dataset: str, labels: list[str]) -> list[dict]:
    rows = []
    for lab in labels:
        sub = df[df["true_label"] == lab]
        n = len(sub)
        correct = int((sub[pred_col] == lab).sum()) if n else 0
        uncertain = int((sub[pred_col] == "uncertain").sum()) if n else 0
        rows.append({
            "dataset": dataset,
            "label": lab,
            "n": n,
            "correct": correct,
            "recall": round(correct / n, 3) if n else np.nan,
            "uncertain": uncertain,
            "uncertain_rate": round(uncertain / n, 3) if n else np.nan,
        })
    total = df[df["true_label"].isin(labels)]
    rows.append({
        "dataset": dataset,
        "label": "ALL",
        "n": len(total),
        "correct": int((total[pred_col] == total["true_label"]).sum()),
        "recall": round(float((total[pred_col] == total["true_label"]).mean()), 3) if len(total) else np.nan,
        "uncertain": int((total[pred_col] == "uncertain").sum()) if len(total) else 0,
        "uncertain_rate": round(float((total[pred_col] == "uncertain").mean()), 3) if len(total) else np.nan,
    })
    return rows


def a1_2b2m_summary(df: pd.DataFrame, pred_col: str, dataset: str) -> list[dict]:
    sub = df[df["true_label"].isin(["2B", "2M"]) & (df["domain"] == "A1")].copy()
    rows = []
    for lab in ["2B", "2M"]:
        for name, part in [
            ("all", sub[sub["true_label"] == lab]),
            ("hotspot", sub[(sub["true_label"] == lab) & sub["hotspot"]]),
            ("non_hotspot", sub[(sub["true_label"] == lab) & ~sub["hotspot"]]),
            ("md_available", sub[(sub["true_label"] == lab) & sub["aim_sb_retained_z"].notna()]),
        ]:
            n = len(part)
            correct = int((part[pred_col] == lab).sum()) if n else 0
            rows.append({
                "dataset": dataset,
                "label": lab,
                "subset": name,
                "n": n,
                "correct": correct,
                "recall": round(correct / n, 3) if n else np.nan,
                "uncertain": int((part[pred_col] == "uncertain").sum()) if n else 0,
                "called_2b": int((part[pred_col] == "2B").sum()) if n else 0,
                "called_2m": int((part[pred_col] == "2M").sum()) if n else 0,
            })
    return rows


def md_priority_queue(df: pd.DataFrame) -> pd.DataFrame:
    """A1 2B/2M cases where MD would be most informative and is absent."""
    sub = df[
        df["true_label"].isin(["2B", "2M"])
        & (df["domain"] == "A1")
        & df["aim_sb_retained_z"].isna()
    ].copy()
    sub["is_error_or_uncertain_no_md"] = (
        (sub["pred_no_md"] != sub["true_label"]) | (sub["pred_no_md"] == "uncertain")
    )
    sub["lof_combined_z"] = (sub["fb_binding_zscore"] + sub["heparan_zscore"]) / 2.0
    sub["priority"] = "low"
    sub.loc[sub["is_error_or_uncertain_no_md"], "priority"] = "high"
    sub.loc[
        (sub["true_label"] == "2B")
        & (sub["pred_no_md"].isin(["2M", "uncertain"]))
        & ~sub["hotspot"],
        "priority",
    ] = "highest_non_hotspot_2B_gap"
    sub["priority_rank"] = sub["priority"].map({
        "highest_non_hotspot_2B_gap": 0,
        "high": 1,
        "low": 2,
    }).fillna(9)
    cols = [
        "priority", "aa_change", "true_label", "protein_pos", "hotspot",
        "pred_no_md", "pred_with_md", "fb_binding_zscore", "heparan_zscore",
        "lof_combined_z", "aim_static_zscore",
    ]
    return sub.sort_values(["priority_rank", "true_label", "aa_change"])[cols]


def type2m_lof_boltz_priority_queue(missing_boltz: pd.DataFrame) -> pd.DataFrame:
    """Clean 2M labels lacking Boltz axes, split by likely LOF mechanism."""
    if missing_boltz.empty:
        return pd.DataFrame()
    sub = missing_boltz[missing_boltz["true_label"] == "2M"].copy()
    if sub.empty:
        return pd.DataFrame()

    sub["protein_pos"] = pd.to_numeric(sub["position_label"], errors="coerce")
    sub["mechanism_track"] = "other_2M_context"
    sub["recommended_boltz_axis"] = "domain_context_panel"
    sub["priority_rank"] = 3

    is_a1 = sub["domain"].eq("A1")
    is_a3 = sub["domain"].eq("A3")
    sub.loc[is_a1, "mechanism_track"] = "A1_GPIb_or_A1_dynamics_LOF"
    sub.loc[is_a1, "recommended_boltz_axis"] = "a1_gpiba_forced_binding + a1_heparan_sulfate_binding + a1_aim_autoinhibition_context"
    sub.loc[is_a1, "priority_rank"] = 0

    sub.loc[is_a3, "mechanism_track"] = "A3_collagen_I_III_binding_LOF"
    sub.loc[is_a3, "recommended_boltz_axis"] = "a3_collagen_binding"
    sub.loc[is_a3, "priority_rank"] = 1

    cols = [
        "aa_change", "true_label", "position_label", "domain", "in_a1_7a6o",
        "mechanism_track", "recommended_boltz_axis", "priority_rank", "sources",
    ]
    return sub.sort_values(["priority_rank", "position_label", "aa_change"])[
        [c for c in cols if c in sub.columns]
    ]


def type2m_lof_md_priority_queue(pred: pd.DataFrame) -> pd.DataFrame:
    """2M cases where MD can add LOF evidence beyond current Boltz axes."""
    sub = pred[pred["true_label"].eq("2M")].copy()
    if sub.empty:
        return pd.DataFrame()

    sub["lof_combined_z"] = (sub["fb_binding_zscore"] + sub["heparan_zscore"]) / 2.0
    sub["priority_rank"] = 9
    sub["mechanism_track"] = "not_prioritized_for_md"
    sub["recommended_md_model"] = ""
    sub["why"] = ""

    a1 = sub["domain"].eq("A1")
    md_missing = sub["aim_sb_retained_z"].isna() & sub["md_face_destab_score"].isna()
    missed_or_uncertain = sub["pred_with_md"].isin(["2B", "uncertain"])
    called_2b = sub["pred_with_md"].eq("2B")
    uncertain = sub["pred_with_md"].eq("uncertain")

    # Highest risk: true 2M currently miscalled 2B. These need A1 dynamics/LOF
    # evidence that can overrule hotspot/allosteric 2B priors.
    mask = a1 & md_missing & called_2b
    sub.loc[mask, "priority_rank"] = 0
    sub.loc[mask, "mechanism_track"] = "A1_2M_miscalled_2B_need_LOF_MD"
    sub.loc[mask, "recommended_md_model"] = "1SQ0_A1_GPIb_complex_MD + 7A6O_AIM_A1_closed_state_MD"
    sub.loc[mask, "why"] = "true 2M is currently called 2B; test GPIb interface loss, A1 misfolding, or hyperstabilized low-adhesion dynamics"

    # Next: true 2M currently uncertain. These are exactly where current static
    # Boltz axes are not directional enough.
    mask = a1 & md_missing & uncertain
    sub.loc[mask, "priority_rank"] = 1
    sub.loc[mask, "mechanism_track"] = "A1_2M_uncertain_need_LOF_MD"
    sub.loc[mask, "recommended_md_model"] = "1SQ0_A1_GPIb_complex_MD first; 7A6O_AIM_A1_closed_state_MD if AIM-adjacent"
    sub.loc[mask, "why"] = "true 2M remains uncertain; current fb/heparan axes do not jointly cross LOF threshold"

    # Controls: already-correct 2M, useful for threshold calibration but not urgent.
    mask = a1 & md_missing & sub["pred_with_md"].eq("2M")
    sub.loc[mask, "priority_rank"] = 2
    sub.loc[mask, "mechanism_track"] = "A1_2M_positive_control_for_LOF_MD"
    sub.loc[mask, "recommended_md_model"] = "1SQ0_A1_GPIb_complex_MD"
    sub.loc[mask, "why"] = "already called 2M; useful to calibrate MD LOF thresholds"

    cols = [
        "aa_change", "protein_pos", "domain", "hotspot", "pred_no_md", "pred_with_md",
        "mechanism_track", "recommended_md_model", "why", "fb_binding_zscore",
        "heparan_zscore", "lof_combined_z", "aim_static_zscore", "aim_sb_retained_z",
        "md_face_destab_score", "priority_rank",
    ]
    out = sub[sub["priority_rank"] < 9].copy()
    return out.sort_values(["priority_rank", "protein_pos", "aa_change"])[
        [c for c in cols if c in out.columns]
    ]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--evidence", default=str(ROOT / "output/boltz2_vwd_functional_panel/evidence_matrix.csv"))
    ap.add_argument("--labels", default=str(ROOT / "output/labeled_variants_all.csv"))
    ap.add_argument("--saltbridge", default=str(ROOT / "output/md_7a6o_saltbridge_features.csv"))
    ap.add_argument("--face-md", default=str(ROOT / "output/md_7a6o_features.csv"))
    ap.add_argument("--classifier", default=str(ROOT / "scripts/agentic_vwf_classifier.py"))
    ap.add_argument("--out-dir", default=str(ROOT / "output/eval_v2"))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mod = load_classifier(Path(args.classifier))

    evidence = pd.read_csv(args.evidence)
    labels_all = pd.read_csv(args.labels)
    eval_df, missing_boltz = load_eval_table(args, mod)

    pred_no_md = run_predictions(eval_df, mod, use_md=False)
    pred_with_md = run_predictions(eval_df, mod, use_md=True)
    pred = eval_df.merge(pred_no_md, on=["aa_change", "true_label"]).merge(
        pred_with_md, on=["aa_change", "true_label"]
    )

    # Leakage smoke test: mutate label columns and confirm predictions do not move.
    leak = eval_df.copy()
    leak["true_label"] = "LEAK_TEST"
    leak["type2_subtype"] = "LEAK_TEST"
    leak_pred = run_predictions(leak, mod, use_md=True)["pred_with_md"].tolist()
    label_leak_test_pass = leak_pred == pred["pred_with_md"].tolist()

    pred.to_csv(out_dir / "eval_v2_predictions.csv", index=False)
    slim_cols = [
        "aa_change", "true_label", "protein_pos", "domain", "hotspot",
        "fb_binding_zscore", "heparan_zscore", "aim_sb_retained_z",
        "a3_collagen_zscore", "md_face_destab_score", "pred_no_md", "confidence_no_md",
        "pred_with_md", "confidence_with_md",
    ]
    pred[[c for c in slim_cols if c in pred.columns]].to_csv(
        out_dir / "eval_v2_predictions_slim.csv", index=False
    )
    confusion(pred, "pred_no_md", TYPE2).to_csv(out_dir / "eval_v2_type2_confusion_no_md.csv", index=False)
    confusion(pred, "pred_with_md", TYPE2).to_csv(out_dir / "eval_v2_type2_confusion_with_md.csv", index=False)

    summary_rows = []
    summary_rows += recall_rows(pred[pred["true_label"].isin(TYPE2)], "pred_no_md", "boltz_clean_type2_no_md", TYPE2)
    summary_rows += recall_rows(pred[pred["true_label"].isin(TYPE2)], "pred_with_md", "boltz_clean_type2_with_existing_md", TYPE2)
    summary_rows += recall_rows(pred[pred["true_label"].isin(["1", *TYPE2])], "pred_with_md", "boltz_clean_supported_labels_with_existing_md", ["1", *TYPE2])
    summary_rows.append({
        "dataset": "leakage_smoke_test",
        "label": "type2_subtype_mutated",
        "n": len(pred),
        "correct": int(label_leak_test_pass),
        "recall": 1.0 if label_leak_test_pass else 0.0,
        "uncertain": np.nan,
        "uncertain_rate": np.nan,
    })
    pd.DataFrame(summary_rows).to_csv(out_dir / "eval_v2_summary.csv", index=False)

    a1_rows = []
    a1_rows += a1_2b2m_summary(pred, "pred_no_md", "a1_2b2m_no_md")
    a1_rows += a1_2b2m_summary(pred, "pred_with_md", "a1_2b2m_with_existing_md")
    pd.DataFrame(a1_rows).to_csv(out_dir / "eval_v2_a1_2b2m_summary.csv", index=False)

    md_priority_queue(pred).to_csv(out_dir / "eval_v2_md_priority_queue.csv", index=False)
    type2m_lof_boltz_priority_queue(missing_boltz).to_csv(
        out_dir / "eval_v2_2m_lof_boltz_priority_queue.csv", index=False
    )
    type2m_lof_md_priority_queue(pred).to_csv(
        out_dir / "eval_v2_2m_lof_md_priority_queue.csv", index=False
    )
    missing_boltz.to_csv(out_dir / "eval_v2_missing_boltz_clean_labels.csv", index=False)
    missing_slim_cols = [
        "aa_change", "true_label", "position_label", "domain", "in_a1_7a6o", "sources",
    ]
    missing_boltz[[c for c in missing_slim_cols if c in missing_boltz.columns]].to_csv(
        out_dir / "eval_v2_missing_boltz_clean_labels_slim.csv", index=False
    )

    inventory_rows = []
    add_inventory_rows(inventory_rows, "evidence_matrix", evidence, "label_policy")
    add_inventory_rows(inventory_rows, "labeled_variants_all", labels_all, "labels")
    add_inventory_rows(inventory_rows, "clean_supported_joined_to_boltz", eval_df, "true_label")
    add_inventory_rows(inventory_rows, "clean_supported_missing_boltz", missing_boltz, "true_label")
    for col in ["is_Type1", "is_Type2A", "is_Type2B", "is_Type2M", "is_Type2N", "is_Type3",
                "is_Negative_GeneBe", "is_Negative_ClinVar"]:
        if col in evidence.columns:
            inventory_rows.append({"dataset": "evidence_matrix", "metric": col, "value": int(boolish(evidence[col]).sum())})
    inventory_rows.append({"dataset": "eval_v2", "metric": "label_leak_test_pass", "value": int(label_leak_test_pass)})
    pd.DataFrame(inventory_rows).to_csv(out_dir / "eval_v2_inventory.csv", index=False)

    print(f"Written Eval v2 outputs -> {out_dir}")
    print(pd.DataFrame(summary_rows).to_string(index=False))
    print("\nA1 2B/2M:")
    print(pd.DataFrame(a1_rows).to_string(index=False))
    print(f"\nLabel leakage smoke test pass: {label_leak_test_pass}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Fast Type2M LOF MD validation after A40/H200 result pulls.

This script is designed for the "pull results, validate immediately" loop. It
uses the small committed analysis tables under
``output/type2m_lof_md_analysis_YYYY-MM-DD`` and does not require the large
trajectory files.

What it does:
  1. Converts completed 7A6O closed-state AIM-A1 contact loss into a provisional
     ``md_face_destab_score`` feature for A1 Type2M LOF.
  2. Merges that feature with the existing MD feature table.
  3. Runs Eval v2 twice: baseline vs baseline+new MD.
  4. Writes a concise report showing recall changes and per-variant prediction
     changes for variants with newly completed MD.

Important limitation:
  A1-GPIb completed QC currently contains RMSD only. RMSD is not a directional
  LOF feature, so this script reports A1-GPIb status but does not feed it into
  the classifier yet. Add interface/contact-retention extraction before using
  A1-GPIb MD as classifier evidence.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TYPE2 = ["2A", "2B", "2M", "2N"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--analysis-dir",
        default=str(ROOT / "output/type2m_lof_md_analysis_2026-06-29"),
        help="A40/H200 partial or full MD analysis directory",
    )
    ap.add_argument(
        "--evidence",
        default=str(
            ROOT
            / "output/hf_type2m_lof_panel/type2m_lof_panel/analysis/"
            / "evidence_matrix_with_type2m_lof_hf.csv"
        ),
    )
    ap.add_argument("--labels", default=str(ROOT / "output/labeled_variants_all.csv"))
    ap.add_argument("--classifier", default=str(ROOT / "scripts/agentic_vwf_classifier.py"))
    ap.add_argument("--saltbridge", default=str(ROOT / "output/md_7a6o_saltbridge_features.csv"))
    ap.add_argument("--existing-face-md", default=str(ROOT / "output/md_7a6o_features.csv"))
    ap.add_argument(
        "--out-dir",
        default=str(ROOT / "output/type2m_lof_md_fast_validation_2026-06-29"),
    )
    ap.add_argument(
        "--closed-contact-loss-threshold",
        type=float,
        default=20.0,
        help="AIM_all first0-5 minus tail40-50 contact loss that maps to md_face_destab_score=1.0",
    )
    ap.add_argument(
        "--skip-eval",
        action="store_true",
        help="Only write feature tables; do not run evaluate_vwf_classifier_v2.py",
    )
    return ap.parse_args()


def read_csv_if_exists(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def load_labels(labels_path: Path) -> pd.DataFrame:
    labels = pd.read_csv(labels_path)
    cols = ["aa_change", "labels", "conflict", "domain", "in_a1_7a6o"]
    return labels[[c for c in cols if c in labels.columns]].copy()


def build_closed_state_features(analysis_dir: Path, labels: pd.DataFrame, threshold: float) -> pd.DataFrame:
    contact_path = analysis_dir / "7a6o_completed_qc/aim_a1_contacts_summary.csv"
    status_path = analysis_dir / "7a6o_closed_state_completed_and_running_summary.csv"
    contacts = read_csv_if_exists(contact_path)
    status = read_csv_if_exists(status_path)
    if contacts.empty:
        return pd.DataFrame(columns=[
            "variant", "md_face_destab_score", "md_closed_aim_contact_loss",
            "md_closed_aim_contact_loss_threshold", "md_feature_source",
        ])

    df = contacts.copy()
    if "AIM_all_contacts_tail_minus_first" in df.columns:
        df["md_closed_aim_contact_loss"] = -pd.to_numeric(
            df["AIM_all_contacts_tail_minus_first"], errors="coerce"
        )
    else:
        df["md_closed_aim_contact_loss"] = (
            pd.to_numeric(df["AIM_all_contacts_first0_5"], errors="coerce")
            - pd.to_numeric(df["AIM_all_contacts_tail40_50"], errors="coerce")
        )

    df["md_face_destab_score"] = (df["md_closed_aim_contact_loss"] / threshold).round(3)
    df["md_closed_aim_contact_loss_threshold"] = threshold
    df["md_feature_source"] = "7a6o_closed_state_AIM_A1_contact_loss_fast"

    if not status.empty and "variant" in status.columns:
        keep = [c for c in ["variant", "status", "xtc_size_mb", "tpr_exists"] if c in status.columns]
        df = df.merge(status[keep], on="variant", how="left", suffixes=("", "_queue"))
    if not labels.empty:
        df = df.merge(labels, left_on="variant", right_on="aa_change", how="left")

    # Keep the classifier-facing columns first, then diagnostic columns.
    front = [
        "variant", "md_face_destab_score", "md_closed_aim_contact_loss",
        "md_closed_aim_contact_loss_threshold", "md_feature_source",
    ]
    rest = [c for c in df.columns if c not in front]
    return df[front + rest].sort_values("variant")


def merge_face_md(existing_path: Path, fast_features: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    if existing_path.exists():
        existing = pd.read_csv(existing_path)
        existing["md_feature_source"] = existing.get("md_feature_source", "existing_md_7a6o_features")
        existing["_priority"] = 0
        pieces.append(existing)

    if not fast_features.empty:
        fast = fast_features.copy()
        fast["_priority"] = 1
        pieces.append(fast)

    if not pieces:
        return pd.DataFrame(columns=["variant", "md_face_destab_score"])

    combined = pd.concat(pieces, ignore_index=True, sort=False)
    combined = combined.sort_values(["variant", "_priority"]).drop_duplicates("variant", keep="last")
    return combined.drop(columns=["_priority"], errors="ignore")


def run_eval(args: argparse.Namespace, face_md: Path, out_dir: Path) -> None:
    cmd = [
        sys.executable,
        str(ROOT / "scripts/pipeline/evaluate_vwf_classifier_v2.py"),
        "--evidence", args.evidence,
        "--labels", args.labels,
        "--saltbridge", args.saltbridge,
        "--face-md", str(face_md),
        "--classifier", args.classifier,
        "--out-dir", str(out_dir),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)


def recall_table(summary_path: Path) -> pd.DataFrame:
    if not summary_path.exists():
        return pd.DataFrame()
    df = pd.read_csv(summary_path)
    return df[df["dataset"].eq("boltz_clean_type2_with_existing_md") & df["label"].isin(TYPE2 + ["ALL"])].copy()


def compare_recalls(base_dir: Path, fast_dir: Path) -> pd.DataFrame:
    base = recall_table(base_dir / "eval_v2_summary.csv")
    fast = recall_table(fast_dir / "eval_v2_summary.csv")
    if base.empty or fast.empty:
        return pd.DataFrame()
    b = base[["label", "n", "correct", "recall", "uncertain"]].rename(
        columns={"n": "n_base", "correct": "correct_base", "recall": "recall_base", "uncertain": "uncertain_base"}
    )
    f = fast[["label", "n", "correct", "recall", "uncertain"]].rename(
        columns={"n": "n_fast", "correct": "correct_fast", "recall": "recall_fast", "uncertain": "uncertain_fast"}
    )
    out = b.merge(f, on="label", how="outer")
    out["delta_correct"] = out["correct_fast"] - out["correct_base"]
    out["delta_recall"] = (out["recall_fast"] - out["recall_base"]).round(3)
    out["delta_uncertain"] = out["uncertain_fast"] - out["uncertain_base"]
    return out


def compare_predictions(base_dir: Path, fast_dir: Path, fast_features: pd.DataFrame) -> pd.DataFrame:
    base_path = base_dir / "eval_v2_predictions.csv"
    fast_path = fast_dir / "eval_v2_predictions.csv"
    if not (base_path.exists() and fast_path.exists()) or fast_features.empty:
        return pd.DataFrame()

    base = pd.read_csv(base_path)
    fast = pd.read_csv(fast_path)
    variants = set(fast_features["variant"].astype(str))
    cols = ["aa_change", "true_label", "domain", "pred_with_md", "confidence_with_md", "reasoning_with_md"]
    base = base[base["aa_change"].isin(variants)][[c for c in cols if c in base.columns]].rename(
        columns={"pred_with_md": "pred_baseline", "confidence_with_md": "confidence_baseline", "reasoning_with_md": "reasoning_baseline"}
    )
    fast = fast[fast["aa_change"].isin(variants)][[c for c in cols if c in fast.columns]].rename(
        columns={"pred_with_md": "pred_fast_md", "confidence_with_md": "confidence_fast_md", "reasoning_with_md": "reasoning_fast_md"}
    )
    out = base.merge(fast, on=["aa_change", "true_label", "domain"], how="outer")
    diag = fast_features[[
        c for c in [
            "variant", "md_face_destab_score", "md_closed_aim_contact_loss",
            "AIM_all_contacts_first0_5", "AIM_all_contacts_tail40_50",
            "AIM_all_contacts_tail_minus_first", "N_AIM_contacts_tail_minus_first",
            "C_AIM_contacts_tail_minus_first", "labels",
        ] if c in fast_features.columns
    ]].rename(columns={"variant": "aa_change", "labels": "source_labels"})
    out = out.merge(diag, on="aa_change", how="left")
    out["prediction_changed"] = out["pred_baseline"].astype(str) != out["pred_fast_md"].astype(str)
    out["rescued_to_true_label"] = (
        out["prediction_changed"]
        & out["true_label"].notna()
        & out["pred_fast_md"].eq(out["true_label"])
    )
    out["strong_md_lof"] = pd.to_numeric(out.get("md_face_destab_score"), errors="coerce") >= 1.0
    out["strong_md_not_rescued"] = (
        out["strong_md_lof"]
        & out["true_label"].eq("2M")
        & ~out["pred_fast_md"].eq("2M")
    )
    return out.sort_values(["rescued_to_true_label", "prediction_changed", "aa_change"], ascending=[False, False, True])


def a1_gpiba_status(analysis_dir: Path, labels: pd.DataFrame) -> pd.DataFrame:
    status = read_csv_if_exists(analysis_dir / "a1_gpiba_completed_and_running_summary.csv")
    qc = read_csv_if_exists(analysis_dir / "a1_gpiba_completed_qc/qc_summary.csv")
    if status.empty:
        return pd.DataFrame()
    df = status.copy()
    if not qc.empty:
        keep = [c for c in ["variant", "rmsd_tail40_50_mean_nm", "rmsd_final_nm", "rmsd_max_nm"] if c in qc.columns]
        df = df.merge(qc[keep], on="variant", how="left")
    if not labels.empty:
        df = df.merge(labels, left_on="variant", right_on="aa_change", how="left")
    return df


def write_report(
    path: Path,
    args: argparse.Namespace,
    fast_features: pd.DataFrame,
    recall_cmp: pd.DataFrame,
    pred_cmp: pd.DataFrame,
    a1_status: pd.DataFrame,
) -> None:
    lines = []
    lines.append("# Type2M LOF MD fast validation")
    lines.append("")
    lines.append(f"Analysis dir: `{args.analysis_dir}`")
    lines.append(f"Closed-state contact-loss threshold: `{args.closed_contact_loss_threshold}` contacts -> `md_face_destab_score=1.0`")
    lines.append("")

    lines.append("## Feature Summary")
    lines.append("")
    if fast_features.empty:
        lines.append("No completed 7A6O closed-state contact summary found.")
    else:
        n_strong = int((fast_features["md_face_destab_score"] >= 1.0).sum())
        lines.append(f"- 7A6O closed-state completed feature rows: `{len(fast_features)}`")
        lines.append(f"- Rows crossing MD_LOF threshold (`md_face_destab_score >= 1.0`): `{n_strong}`")
        show_cols = [
            "variant", "labels", "md_face_destab_score", "md_closed_aim_contact_loss",
            "AIM_all_contacts_first0_5", "AIM_all_contacts_tail40_50",
        ]
        show = fast_features[[c for c in show_cols if c in fast_features.columns]].sort_values(
            "md_face_destab_score", ascending=False
        )
        lines.append("")
        lines.append(show.to_markdown(index=False))

    lines.append("")
    lines.append("## Recall Delta")
    lines.append("")
    if recall_cmp.empty:
        lines.append("Eval was skipped or summary files were not found.")
    else:
        lines.append(recall_cmp.to_markdown(index=False))

    lines.append("")
    lines.append("## Prediction Changes On Newly Completed MD Variants")
    lines.append("")
    if pred_cmp.empty:
        lines.append("No prediction comparison available.")
    else:
        show_cols = [
            "aa_change", "true_label", "pred_baseline", "pred_fast_md",
            "md_face_destab_score", "md_closed_aim_contact_loss",
            "prediction_changed", "rescued_to_true_label", "strong_md_not_rescued",
        ]
        lines.append(pred_cmp[[c for c in show_cols if c in pred_cmp.columns]].to_markdown(index=False))
        blocked = pred_cmp[pred_cmp.get("strong_md_not_rescued", False).eq(True)]
        if len(blocked):
            lines.append("")
            lines.append(
                f"Strong MD LOF but not rescued: `{len(blocked)}`. These are rule-order/2B-prior review candidates."
            )

    lines.append("")
    lines.append("## A1-GPIb Status")
    lines.append("")
    lines.append(
        "A1-GPIb MD is reported here for QC only. RMSD alone is not fed into the classifier; "
        "a directional GPIb interface/contact feature should be extracted before using it as LOF evidence."
    )
    lines.append("")
    if a1_status.empty:
        lines.append("No A1-GPIb status table found.")
    else:
        counts = a1_status.groupby(["status", "labels"], dropna=False).size().reset_index(name="n")
        lines.append(counts.to_markdown(index=False))

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    analysis_dir = Path(args.analysis_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    labels = load_labels(Path(args.labels))
    fast_features = build_closed_state_features(
        analysis_dir=analysis_dir,
        labels=labels,
        threshold=args.closed_contact_loss_threshold,
    )
    fast_feature_path = out_dir / "type2m_lof_md_fast_features.csv"
    fast_features.to_csv(fast_feature_path, index=False)

    combined_face = merge_face_md(Path(args.existing_face_md), fast_features)
    combined_face_path = out_dir / "combined_face_md_features.csv"
    combined_face.to_csv(combined_face_path, index=False)

    a1_status_df = a1_gpiba_status(analysis_dir, labels)
    a1_status_path = out_dir / "a1_gpiba_qc_status_for_review.csv"
    a1_status_df.to_csv(a1_status_path, index=False)

    base_eval = out_dir / "eval_baseline"
    fast_eval = out_dir / "eval_with_fast_md"
    recall_cmp = pd.DataFrame()
    pred_cmp = pd.DataFrame()
    if not args.skip_eval:
        run_eval(args, Path(args.existing_face_md), base_eval)
        run_eval(args, combined_face_path, fast_eval)
        recall_cmp = compare_recalls(base_eval, fast_eval)
        pred_cmp = compare_predictions(base_eval, fast_eval, fast_features)
        recall_cmp.to_csv(out_dir / "recall_delta.csv", index=False)
        pred_cmp.to_csv(out_dir / "completed_md_prediction_delta.csv", index=False)

    report_path = out_dir / "TYPE2M_LOF_MD_FAST_VALIDATION_REPORT.md"
    write_report(report_path, args, fast_features, recall_cmp, pred_cmp, a1_status_df)

    print(f"wrote {fast_feature_path}")
    print(f"wrote {combined_face_path}")
    print(f"wrote {a1_status_path}")
    if not args.skip_eval:
        print(f"wrote {out_dir / 'recall_delta.csv'}")
        print(f"wrote {out_dir / 'completed_md_prediction_delta.csv'}")
    print(f"wrote {report_path}")
    if not recall_cmp.empty:
        print("\nRecall delta:")
        print(recall_cmp.to_string(index=False))
    if not pred_cmp.empty:
        changed = pred_cmp[pred_cmp["prediction_changed"]]
        print("\nPrediction changes:")
        print(changed[[
            c for c in ["aa_change", "true_label", "pred_baseline", "pred_fast_md", "md_face_destab_score"]
            if c in changed.columns
        ]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

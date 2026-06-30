#!/usr/bin/env python3
"""Sweep provisional 7A6O closed-state MD LOF thresholds.

This is a thin reproducibility wrapper around validate_type2m_lof_md_fast.py.
It answers the immediate question after a partial A40/H200 MD pull:

  - If contact loss is mapped into md_face_destab_score, which threshold is
    useful?
  - Does the partial MD rescue 2M cases, and does it change other classes?

Important: with the 2026-06-30 partial pull, completed 7A6O closed-state MD
features are all clean 2M labels. The sweep can validate 2M rescue behavior,
but it cannot calibrate 2B false-positive risk until 2B hard-negative MD runs
complete.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--thresholds",
        default="10,15,20,25,30,35,40",
        help="Comma-separated contact-loss thresholds. score=loss/threshold; score>=1 enters MD LOF rule.",
    )
    ap.add_argument(
        "--out-dir",
        default=str(ROOT / "output/type2m_lof_md_threshold_sweep_2026-06-30"),
    )
    ap.add_argument(
        "--validator",
        default=str(ROOT / "scripts/pipeline/validate_type2m_lof_md_fast.py"),
    )
    ap.add_argument("--analysis-dir", default=str(ROOT / "output/type2m_lof_md_analysis_2026-06-29"))
    ap.add_argument("--evidence", default=None)
    ap.add_argument("--labels", default=None)
    ap.add_argument("--classifier", default=None)
    ap.add_argument("--saltbridge", default=None)
    ap.add_argument("--existing-face-md", default=None)
    return ap.parse_args()


def validator_cmd(args: argparse.Namespace, threshold: float, run_dir: Path) -> list[str]:
    cmd = [
        sys.executable,
        args.validator,
        "--analysis-dir",
        args.analysis_dir,
        "--closed-contact-loss-threshold",
        str(threshold),
        "--out-dir",
        str(run_dir),
    ]
    optional = {
        "--evidence": args.evidence,
        "--labels": args.labels,
        "--classifier": args.classifier,
        "--saltbridge": args.saltbridge,
        "--existing-face-md": args.existing_face_md,
    }
    for flag, value in optional.items():
        if value:
            cmd.extend([flag, value])
    return cmd


def summarize_run(threshold: float, run_dir: Path) -> dict:
    recall = pd.read_csv(run_dir / "recall_delta.csv")
    pred = pd.read_csv(run_dir / "completed_md_prediction_delta.csv")
    feats = pd.read_csv(run_dir / "type2m_lof_md_fast_features.csv")

    def recall_row(label: str) -> pd.Series:
        return recall[recall["label"].eq(label)].iloc[0]

    r2m = recall_row("2M")
    r2b = recall_row("2B")
    rall = recall_row("ALL")
    return {
        "contact_loss_threshold": threshold,
        "completed_md_rows": len(feats),
        "completed_md_label_counts": ";".join(
            f"{k}:{v}" for k, v in feats.get("labels", pd.Series(dtype=str)).value_counts(dropna=False).items()
        ),
        "strong_md_rows": int((feats["md_face_destab_score"] >= 1.0).sum()),
        "changed_predictions": int(pred["prediction_changed"].sum()),
        "rescued_to_true": int(pred["rescued_to_true_label"].sum()),
        "2M_delta_correct": int(r2m["delta_correct"]),
        "2M_recall_fast": float(r2m["recall_fast"]),
        "2B_delta_correct": int(r2b["delta_correct"]),
        "2B_recall_fast": float(r2b["recall_fast"]),
        "ALL_delta_correct": int(rall["delta_correct"]),
        "ALL_recall_fast": float(rall["recall_fast"]),
        "ALL_delta_uncertain": int(rall["delta_uncertain"]),
    }


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    thresholds = [float(x.strip()) for x in args.thresholds.split(",") if x.strip()]

    rows = []
    for threshold in thresholds:
        name = f"t{threshold:g}".replace(".", "p")
        run_dir = out_dir / name
        run_dir.mkdir(parents=True, exist_ok=True)
        cmd = validator_cmd(args, threshold, run_dir)
        subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        rows.append(summarize_run(threshold, run_dir))

    summary = pd.DataFrame(rows)
    summary_path = out_dir / "threshold_sweep_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(summary.to_string(index=False))
    print(f"\nwrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

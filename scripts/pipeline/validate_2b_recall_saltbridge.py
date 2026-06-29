#!/usr/bin/env python3
"""Validate the AIM salt-bridge joint criterion's effect on 2B recall (cluster).

Run on the cluster where VWF_Alpha_Matrix.parquet (with fb_binding_zscore /
heparan / AF3 features + type2_subtype labels) is available. This no longer
splits or decides by recurrent Type 2B hotspot positions; position priors were
removed from the production classifier as a leakage/circularity risk.

What it does:
  1. Load the feature matrix.
  2. Left-join aim_sb_retained_z (output/md_7a6o_saltbridge_features.csv) and, if
     absent, md_face_destab_score (output/md_7a6o_features.csv) on aa_change.
  3. Run the classifier twice -- WITH the salt-bridge column and WITHOUT it
     (ablation = current baseline) -- and report 2B/2M recall overall plus every
     variant whose call changed.

Usage:
    python3 scripts/pipeline/validate_2b_recall_saltbridge.py \
        --matrix VWF_Alpha_Matrix.parquet \
        --sb output/md_7a6o_saltbridge_features.csv \
        --face output/md_7a6o_features.csv
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


def load_classifier(path="scripts/agentic_vwf_classifier.py"):
    spec = importlib.util.spec_from_file_location("agentic_vwf_classifier", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def recall_block(df, pred_col, label_col):
    """Return recall strings for 2B/2M overall."""
    out = {}
    for lab in ("2B", "2M"):
        sub = df[df[label_col] == lab]
        if len(sub) == 0:
            out[lab] = "n/a"
            continue
        def r(s):
            return f"{int((s[pred_col] == lab).sum())}/{len(s)}" if len(s) else "0/0"
        out[lab] = f"overall {r(sub)}"
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", default="VWF_Alpha_Matrix.parquet")
    ap.add_argument("--sb", default="output/md_7a6o_saltbridge_features.csv")
    ap.add_argument("--face", default="output/md_7a6o_features.csv")
    ap.add_argument("--label-col", default="type2_subtype")
    ap.add_argument("--key", default="aa_change", help="matrix join key vs MD csv 'variant'")
    args = ap.parse_args()

    m = load_classifier()

    df = pd.read_parquet(args.matrix)
    print(f"matrix: {len(df)} rows; label col '{args.label_col}' present: "
          f"{args.label_col in df.columns}")

    # join MD features by aa_change <- variant
    sb = pd.read_csv(args.sb)[["variant", "aim_sb_retained_z"]]
    df = df.merge(sb, left_on=args.key, right_on="variant", how="left").drop(columns=["variant"])
    if "md_face_destab_score" not in df.columns and Path(args.face).exists():
        fa = pd.read_csv(args.face)[["variant", "md_face_destab_score"]]
        df = df.merge(fa, left_on=args.key, right_on="variant", how="left").drop(columns=["variant"])
    n_sb = int(df["aim_sb_retained_z"].notna().sum())
    print(f"variants with aim_sb_retained_z after join: {n_sb}")

    labeled = df[df[args.label_col].notna()].copy()
    labeled = labeled.rename(columns={args.label_col: "type2_subtype"})

    clf = m.AgenticVWFClassifier()

    # WITH salt-bridge feature
    pred_with = clf.classify_batch(labeled)["main_subtype"].values
    # WITHOUT (ablation baseline): null the column
    abl = labeled.copy()
    abl["aim_sb_retained_z"] = np.nan
    clf2 = m.AgenticVWFClassifier()
    pred_base = clf2.classify_batch(abl)["main_subtype"].values

    labeled = labeled.assign(pred_base=pred_base, pred_with=pred_with)
    lc = "type2_subtype"

    print("\n=== 2B / 2M recall (baseline = no SB feature) ===")
    rb = recall_block(labeled, "pred_base", lc)
    rw = recall_block(labeled, "pred_with", lc)
    for lab in ("2B", "2M"):
        print(f"\n{lab}:")
        print(f"  baseline : {rb[lab]}")
        print(f"  with SB  : {rw[lab]}")

    changed = labeled[labeled.pred_base != labeled.pred_with]
    print(f"\n=== variants whose call changed (n={len(changed)}) ===")
    cols = [args.key, lc, "protein_pos", "aim_sb_retained_z", "pred_base", "pred_with"]
    cols = [c for c in cols if c in changed.columns]
    if len(changed):
        print(changed[cols].to_string(index=False))
    else:
        print("(none)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

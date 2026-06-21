#!/usr/bin/env python3
"""Analyse 7A6O AIM-unfolding steered-MD pulls -> 2B-positive feature.

For each variant/replicate reads GROMACS pull force/coord (`*_pullf.xvg`,
`*_pullx.xvg`) and computes:
  - rupture force  = peak pull force (smoothed)  [pN]
  - pulling work   = integral F dx up to the peak [kJ/mol]
Aggregates over replicates (rupture force is stochastic -> report mean +/- sd).

Medical reading: type 2B (gain of function) LOWERS the force needed to unfold
the AIM (Nat Commun 2021). So **low `smd_unfold_force_pN` => 2B lean**. This is
the force-dependent 2B axis that equilibrium MD cannot see; it complements the
equilibrium binding-face axis (md_face_destab_score, a 2M/LOF axis).

Output: output/md_7a6o_smd_features.csv

Usage:
    python3 scripts/pipeline/analyze_7a6o_smd.py \
        --input output/gromacs_md_autoinhib --output output/md_7a6o_smd_features.csv
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import numpy as np
import pandas as pd

KJMOLNM_TO_PN = 1.66054  # 1 kJ/mol/nm = 1.66 pN

LABELS = {
    "WT": "WT", "R1306W": "2B", "R1306Q": "2B", "R1308C": "2B", "I1309V": "2B",
    "S1310F": "2B", "W1313C": "2B", "V1314F": "2B", "V1316M": "2B",
    "R1374C": "2M", "R1374H": "2M", "G1324S": "2M",
    "P1337L": "?", "R1341Q": "?", "R1341W": "?", "C1458R": "patient",
}


def read_xvg(path: str):
    t, y = [], []
    for ln in open(path, errors="ignore"):
        if not ln.strip() or ln[0] in "#@":
            continue
        p = ln.split()
        if len(p) >= 2:
            t.append(float(p[0])); y.append(float(p[1]))
    return np.array(t), np.array(y)


def smooth(y, w=11):
    if len(y) < w:
        return y
    k = np.ones(w) / w
    return np.convolve(y, k, mode="same")


def analyse_pull(pf: str, px: str):
    tf, force = read_xvg(pf)
    fsm = smooth(np.abs(force))
    peak_idx = int(np.argmax(fsm))
    peak_f = float(fsm[peak_idx])
    work = np.nan
    if os.path.exists(px):
        tx, coord = read_xvg(px)
        n = min(len(coord), len(force))
        # work up to rupture: integral of F over dx
        ipk = min(peak_idx, n - 1)
        work = float(np.trapz(np.abs(force[:ipk + 1]), coord[:ipk + 1]))
    return peak_f * KJMOLNM_TO_PN, work, tf[peak_idx] if len(tf) else np.nan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="output/gromacs_md_autoinhib")
    ap.add_argument("--output", default="output/md_7a6o_smd_features.csv")
    args = ap.parse_args()

    per_rep = []
    for pf in sorted(glob.glob(os.path.join(args.input, "*", "smd", "*_pullf.xvg"))):
        variant = pf.split(os.sep)[-3]
        rep = re.search(r"rep(\d+)", os.path.basename(pf))
        rep = int(rep.group(1)) if rep else 0
        px = pf.replace("_pullf.xvg", "_pullx.xvg")
        gro = re.sub(r"_pullf\.xvg$", ".gro", pf)
        # A killed/failed mdrun can leave pullf/pullx without a final .gro.
        # Treat only completed replicas as analysable features.
        if not (os.path.exists(px) and os.path.exists(gro)):
            print(f"Skipping incomplete SMD replica: {pf}")
            continue
        peak_pN, work, t_peak = analyse_pull(pf, px)
        per_rep.append(dict(variant=variant, label=LABELS.get(variant, "?"), rep=rep,
                            rupture_force_pN=round(peak_pN, 1),
                            work_kJmol=round(work, 1) if work == work else np.nan,
                            t_peak_ps=round(t_peak, 1)))
    if not per_rep:
        print(f"No *_pullf.xvg found under {args.input}/*/smd/. Run run_7a6o_smd.sh first.")
        return 1
    rep_df = pd.DataFrame(per_rep)
    rep_df.to_csv(args.output.replace(".csv", "_perrep.csv"), index=False)

    agg = rep_df.groupby(["variant", "label"]).agg(
        n_reps=("rep", "count"),
        smd_unfold_force_pN=("rupture_force_pN", "mean"),
        smd_unfold_force_sd=("rupture_force_pN", "std"),
        smd_work_kJmol=("work_kJmol", "mean"),
    ).reset_index().round(2)

    # z-score over labeled reference (WT+2B+2M).
    # Two complementary 2B-positive metrics (higher = easier AIM release = more 2B):
    #   peak rupture force -- sensitive to a single-frame spike, rate-dominated at fast pull;
    #   pulling work        -- integral up to rupture, smoother / less spike-sensitive.
    # Report both; prefer whichever separates known 2B from 2M/WT at the chosen pull rate.
    def zscore(col, invert):
        ref = agg[agg.label.isin(["WT", "2B", "2M"])][col]
        mu, sd = ref.mean(), ref.std()
        z = (agg[col] - mu) / sd
        return z.round(3), (-z if invert else z).round(3)

    agg["smd_force_z"], agg["smd_2b_score"] = zscore("smd_unfold_force_pN", invert=True)
    agg["smd_work_z"], agg["smd_2b_score_work"] = zscore("smd_work_kJmol", invert=True)
    agg.to_csv(args.output, index=False)

    print(f"Written: {args.output}  (+ _perrep.csv)")
    pd.set_option("display.width", 200)
    print(agg.sort_values("smd_unfold_force_pN").to_string(index=False))
    for lab in ["WT", "2B", "2M", "?"]:
        s = agg[agg.label == lab]
        if len(s):
            print(f"  {lab}: force {s['smd_unfold_force_pN'].mean():.0f}+/-{s['smd_unfold_force_pN'].std():.0f} pN | "
                  f"work {s['smd_work_kJmol'].mean():.0f} kJ/mol  (n={len(s)})")
    # quick 2B-vs-2M separation check (force; lower-in-2B expected)
    b = agg[agg.label == "2B"]["smd_unfold_force_pN"].values
    m = agg[agg.label == "2M"]["smd_unfold_force_pN"].values
    if len(b) and len(m):
        import itertools
        auc = np.mean([1.0 if bv < mv else 0.5 if bv == mv else 0.0
                       for bv, mv in itertools.product(b, m)])
        print(f"\nAUC(2B force < 2M) = {auc:.2f}  (1.0 = 2B always unfolds at lower force; "
              f"<=0.5 = no/backwards separation -> revisit pull rate before wiring into classifier)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

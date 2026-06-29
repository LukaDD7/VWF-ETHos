#!/usr/bin/env python3
"""Validate classifier 2B recall on the expanded new-2B Boltz panel (no MD).

Context (2026-06-24): the 62 new 2B variants were run through the functional
Boltz-2 panel (output/boltz2_new_2b_panel) to get fb_binding / heparan iPTM. None
of them have 7A6O MD trajectories, so this measures how the *static*
forced_binding/heparan logic generalizes to a larger 2B set -- i.e. the recall
ceiling before MD-derived mechanism evidence can contribute. Recurrent Type 2B
hotspot positions are not used, reported, or split because that position prior
was removed from the production classifier as a leakage/circularity risk.

Critical normalization note: z-scores MUST be computed against the cross-subtype
distribution (the full functional panel pools 2A/2B/2M/2N), NOT within the 2B-only
new batch. The LOF gate (LOF_COMBINED_Z) is calibrated on cross-subtype z; a
within-2B-batch z makes "low z" mean "lower than other 2B" and misfires the gate.
We therefore pool the new raw primary_values with the full panel before z-scoring.

Outputs:
  output/new_2b_panel_recall.csv  per-variant call + features
prints overall 2B recall and the false-2M list.

Usage:
    python3 scripts/pipeline/validate_new_2b_panel_recall.py
"""
from __future__ import annotations

import argparse
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd

FB_ASSAY = "a1_gpiba_forced_binding"
HEP_ASSAY = "a1_heparan_sulfate_binding"
# variants whose source 2B label conflicts with an existing 2M reference label
CONFLICT_REF_2M = {"R1374C", "R1374H"}

# columns fit()/classify expect to exist; NaN -> graceful fallback in RULE6
EXPECTED_NA_COLS = [
    "ag_rna_delta", "ag_splice_delta", "ag_delta_score", "af3_plddt_mean",
    "af3_plddt_min", "af3_pae_interface", "foldx_ddg_bind", "boltz2_iptm",
    "boltz2_delta_iptm", "aim_release_score", "aim_sb_retained_z",
    "md_face_destab_score",
]


def load_classifier(path="scripts/agentic_vwf_classifier.py"):
    spec = importlib.util.spec_from_file_location("agentic_vwf_classifier", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", default="output/boltz2_new_2b_panel/boltz_results_summary.csv")
    ap.add_argument("--panel", default="output/boltz2_vwd_functional_panel/evidence_matrix.csv",
                    help="full functional panel evidence matrix (cross-subtype reference)")
    ap.add_argument("--out", default="output/new_2b_panel_recall.csv")
    args = ap.parse_args()

    m = load_classifier()

    # cross-subtype reference distribution from the full panel
    ev = pd.read_csv(args.panel)
    fb_col, hep_col = f"{FB_ASSAY}__primary_value", f"{HEP_ASSAY}__primary_value"
    old_fb, old_hep = ev[fb_col], ev[hep_col]

    nw = pd.read_csv(args.new)
    nw = nw[nw.variant_id.notna()]
    nfb = nw[nw.assay_key == FB_ASSAY].set_index("variant_id")["primary_value"]
    nhep = nw[nw.assay_key == HEP_ASSAY].set_index("variant_id")["primary_value"]

    pool_fb = pd.concat([old_fb, nfb]); pool_hep = pd.concat([old_hep, nhep])
    fb_mu, fb_sd = pool_fb.mean(), pool_fb.std(ddof=0)
    hep_mu, hep_sd = pool_hep.mean(), pool_hep.std(ddof=0)

    rows = []
    for vid in sorted(set(nfb.index) - {"VWF_WT"}):
        aa = vid.replace("VWF_", "")
        mo = re.match(r"^([A-Z])(\d+)([A-Z])$", aa)
        if not mo:
            continue
        fb, hep = nfb.get(vid, np.nan), nhep.get(vid, np.nan)
        rows.append(dict(aa_change=aa, protein_pos=int(mo.group(2)),
                         ref_aa=mo.group(1), alt_aa=mo.group(3),
                         domain="A1", type2_subtype="2B",
                         fb_binding_zscore=(fb - fb_mu) / fb_sd,
                         heparan_zscore=(hep - hep_mu) / hep_sd))
    X = pd.DataFrame(rows)
    for c in EXPECTED_NA_COLS:
        X[c] = np.nan

    res = m.AgenticVWFClassifier().classify_batch(X)
    X = pd.concat([X.reset_index(drop=True), res[["main_subtype", "confidence"]]], axis=1)
    X["conflict_ref_2M"] = X.aa_change.isin(CONFLICT_REF_2M)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    X.to_csv(args.out, index=False)

    print(f"pooled FB  n={pool_fb.notna().sum()} mu={fb_mu:.3f} sd={fb_sd:.3f}")
    print(f"pooled HEP n={pool_hep.notna().sum()} mu={hep_mu:.3f} sd={hep_sd:.3f}")
    print(f"\n=== {len(X)} new A1 2B variants -> {args.out} ===")
    print(X["main_subtype"].value_counts().to_string())

    clean = X[~X.conflict_ref_2M]
    print(f"\n=== 2B recall (clean, n={len(clean)}; no MD features) ===")
    n = len(clean); c = int((clean.main_subtype == "2B").sum())
    print(f"  all: {c}/{n}" + (f"  ({100*c/n:.0f}%)" if n else ""))

    fp = clean[clean.main_subtype == "2M"]
    print("\n=== false 2M (LOF gate misfired on real 2B; needs MD salt-bridge) ===")
    print(fp[["aa_change", "protein_pos", "fb_binding_zscore",
              "heparan_zscore"]].round(2).to_string(index=False) if len(fp) else "(none)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

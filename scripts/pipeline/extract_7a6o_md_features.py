#!/usr/bin/env python3
"""Extract 2B/2M discriminating features from the 7A6O AIM-A1 reference MD set.

Medical rationale (see docs/7A6O_MD_FEATURE_ANALYSIS_2026-06-20.md)
--------------------------------------------------------------------
The 7A6O construct (native VWF residues 1262-1466) is the A1 domain plus its
discontinuous Autoinhibitory Module (AIM: N-AIM 1262-1271, C-AIM 1459-1466).
In the resting state the AIM masks the GPIbalpha-binding face of A1 (helix a1,
loops a1b2 and b3a2; Deng et al. Blood 2017, HDX). Type 2B (gain of function)
mutations lower the *force* needed to unfold the AIM, but at equilibrium (no
shear) the binding face stays intact -- which is exactly why 2B VWF still folds
and can bind GPIbalpha. Type 2M (loss of function) mutations damage the A1
binding face itself, so even the resting AIM cannot mask it cleanly.

Consequence for features: in unforced 50 ns MD the medically meaningful, robust
observable is **how much of the WT AIM->binding-face masking contact network the
variant retains** (binding-face integrity), NOT the raw AIM<->A1 total contact
count (which moves in the *opposite*, mechanistically misleading direction).

Low masking-contact retention  => binding face destabilized  => 2M / LOF lean.
Near-WT retention              => face intact                => 2B-compatible.

This script is read-only over the trajectories. It writes a feature table to
``output/md_7a6o_features.csv`` plus the data-derived masking interface.

Usage:
    python3 scripts/pipeline/extract_7a6o_md_features.py \
        --input md_data/7a6o_reference_md/variants \
        --output output/md_7a6o_features.csv
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import MDAnalysis as mda  # noqa: E402
from MDAnalysis.analysis import align, rms  # noqa: E402
from MDAnalysis.lib.distances import capped_distance, distance_array  # noqa: E402

# native_resid = local_resid + OFFSET  (verified by matching all 15 mutation sites)
OFFSET = 1261
# Local-numbering selections (MDAnalysis renumbers the 205-res construct 1..205)
SEL_NAIM = "resid 1:10"        # native 1262-1271
SEL_CAIM = "resid 199:205"     # native 1460-1466
SEL_AIM = "(resid 1:10 or resid 199:205)"
SEL_A1 = "resid 19:190"        # native 1280-1451 (A1 body, matches gmx A1_nonlocal)
CONTACT_CUT = 4.5              # angstrom, heavy-atom contact cutoff
TAIL_NS = 40.0                 # tail window start (40-50 ns)
MASK_OCC = 0.5                 # WT occupancy to define a masking contact pair

LABELS = {
    "WT": "WT", "WT_r1": "WT", "WT_r2": "WT", "WT_r3": "WT",
    "R1306W": "2B", "R1306Q": "2B", "R1308C": "2B", "I1309V": "2B", "S1310F": "2B",
    "W1313C": "2B", "V1314F": "2B", "V1316M": "2B",
    "R1374C": "2M", "R1374H": "2M", "G1324S": "2M",
    "P1337L": "?", "R1341Q": "?", "R1341W": "?",
    "C1458R": "patient", "C1458R_r1": "patient", "C1458R_r2": "patient", "C1458R_r3": "patient",
}


def load(variant_dir: Path):
    tpr = variant_dir / "md_prod.tpr"
    xtc = variant_dir / "prod_concat.xtc"
    if not (tpr.exists() and xtc.exists()):
        return None
    return mda.Universe(str(tpr), str(xtc))


def residpair_occupancy(u, tail_only: bool, stride: int = 2):
    """Occupancy fraction per (AIM_resid, A1_resid) heavy-atom contact pair."""
    aim = u.select_atoms(f"{SEL_AIM} and not name H*")
    a1 = u.select_atoms(f"{SEL_A1} and not name H*")
    aim_rid, a1_rid = aim.resids, a1.resids
    counts: dict[tuple[int, int], int] = {}
    nframes = 0
    for ts in u.trajectory[::stride]:
        if tail_only and ts.time / 1000.0 < TAIL_NS:
            continue
        nframes += 1
        D = distance_array(aim.positions, a1.positions)
        ii, jj = np.where(D < CONTACT_CUT)
        seen = set()
        for i, j in zip(ii, jj):
            key = (int(aim_rid[i]), int(a1_rid[j]))
            if key not in seen:
                seen.add(key)
                counts[key] = counts.get(key, 0) + 1
    return {k: v / max(nframes, 1) for k, v in counts.items()}


def scalar_features(u):
    aim_h = u.select_atoms(f"{SEL_AIM} and not name H*")
    a1_h = u.select_atoms(f"{SEL_A1} and not name H*")
    naim_h = u.select_atoms(f"{SEL_NAIM} and not name H*")
    caim_h = u.select_atoms(f"{SEL_CAIM} and not name H*")
    times, ncont, mind, nd_n, nd_c = [], [], [], [], []
    for ts in u.trajectory:
        times.append(ts.time / 1000.0)
        ncont.append(len(capped_distance(aim_h.positions, a1_h.positions, CONTACT_CUT,
                                          return_distances=False)))
        da = distance_array(naim_h.positions, a1_h.positions)
        db = distance_array(caim_h.positions, a1_h.positions)
        mind.append(min(da.min(), db.min()) / 10.0)  # nm
        nd_n.append(da.min() / 10.0)
        nd_c.append(db.min() / 10.0)
    times = np.array(times)
    tail = times >= TAIL_NS
    out = dict(
        n_frames=len(times),
        ncont_heavy_mean=round(float(np.mean(ncont)), 1),
        ncont_heavy_tail=round(float(np.mean(np.array(ncont)[tail])), 1),
        aim_a1_mindist_tail_nm=round(float(np.mean(np.array(mind)[tail])), 3),
        naim_a1_mindist_tail_nm=round(float(np.mean(np.array(nd_n)[tail])), 3),
        caim_a1_mindist_tail_nm=round(float(np.mean(np.array(nd_c)[tail])), 3),
    )
    # AIM RMSF (angstrom) after aligning on A1 core CA
    u.trajectory[0]
    align.AlignTraj(u, u, select=f"{SEL_A1} and name CA", in_memory=True).run()
    aim_ca = u.select_atoms(f"{SEL_AIM} and name CA")
    rmsf = rms.RMSF(aim_ca).run().results.rmsf
    out["aim_rmsf_mean_A"] = round(float(np.mean(rmsf)), 3)
    out["aim_rmsf_max_A"] = round(float(np.max(rmsf)), 3)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default="md_data/7a6o_reference_md/variants")
    ap.add_argument("--output", default="output/md_7a6o_features.csv")
    ap.add_argument("--wt", default="WT", help="variant name used to define the masking interface")
    ap.add_argument("--stride", type=int, default=2, help="frame stride for contact maps")
    args = ap.parse_args()

    root = Path(args.input)
    variants = sorted([d.name for d in root.iterdir() if d.is_dir() and (d / "md_prod.tpr").exists()])
    if args.wt not in variants:
        raise SystemExit(f"WT reference '{args.wt}' not found under {root}")

    # 1) Define WT masking interface = AIM->A1 pairs masked >= MASK_OCC of the WT run
    uwt = load(root / args.wt)
    wt_occ = residpair_occupancy(uwt, tail_only=False, stride=args.stride)
    mask_pairs = {k for k, v in wt_occ.items() if v >= MASK_OCC}
    a1_mask_native = sorted({j + OFFSET for (_, j) in mask_pairs})
    print(f"WT masking interface: {len(mask_pairs)} AIM->A1 contact pairs")
    print(f"  A1 binding-face residues masked by AIM (native): {a1_mask_native}")

    rows = []
    for v in variants:
        u = load(root / v)
        feats = scalar_features(u)
        # re-load (scalar_features consumed/aligned the trajectory in memory)
        u2 = load(root / v)
        tail_occ = residpair_occupancy(u2, tail_only=True, stride=args.stride)
        retained = sum(1 for p in mask_pairs if tail_occ.get(p, 0.0) >= MASK_OCC)
        mask_occ_tail = sum(min(tail_occ.get(p, 0.0), 1.0) for p in mask_pairs) / len(mask_pairs)
        rows.append(dict(
            variant=v, label=LABELS.get(v, "?"),
            md_aim_mask_retention=round(mask_occ_tail, 4),
            md_mask_pairs_retained=retained, md_mask_pairs_total=len(mask_pairs),
            **feats,
        ))
        print(f"  {v:8s} ({LABELS.get(v,'?'):7s}) mask_retention={mask_occ_tail:.3f}")

    df = pd.DataFrame(rows)
    # z-score masking retention over the labeled reference set (WT + known 2B + 2M only)
    ref = df[df.label.isin(["WT", "2B", "2M"])]["md_aim_mask_retention"]
    mu, sd = ref.mean(), ref.std()
    df["md_mask_z"] = ((df["md_aim_mask_retention"] - mu) / sd).round(3)
    # face-destabilization score: higher = more binding-face damage = 2M/LOF lean
    df["md_face_destab_score"] = (-df["md_mask_z"]).round(3)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    (out.parent / "md_7a6o_masking_interface.json").write_text(json.dumps({
        "construct_native_range": [OFFSET + 1, OFFSET + 205],
        "offset_native_minus_local": OFFSET,
        "contact_cutoff_A": CONTACT_CUT, "wt_occupancy_def": MASK_OCC,
        "a1_binding_face_residues_native": a1_mask_native,
        "n_mask_pairs": len(mask_pairs),
        "zscore_ref_mean": round(float(mu), 4), "zscore_ref_sd": round(float(sd), 4),
    }, indent=2))
    print(f"\nWritten: {out}")
    print(f"Written: {out.parent / 'md_7a6o_masking_interface.json'}")
    pd.set_option("display.width", 200, "display.max_columns", 50)
    print(df[["variant", "label", "md_aim_mask_retention", "md_mask_z",
              "md_face_destab_score", "ncont_heavy_tail", "aim_rmsf_mean_A"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

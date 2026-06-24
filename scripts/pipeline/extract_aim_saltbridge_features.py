#!/usr/bin/env python3
"""Targeted autoinhibitory salt-bridge occupancy from the 7A6O AIM-A1 MD set.

Why this script (motivation from literature, 2026-06-24)
--------------------------------------------------------
SMD force axis was a no-go (docs/7A6O_SMD_LITERATURE_NOGO_2026-06-24.md): at the
accessible pull rate the rupture force is friction-dominated (~1000 vs ~10-20 pN
in optical tweezers) and the 2B/2M direction inverts. The literature
(RSC Chem Biol 2022, PMC9175105) localizes the autoinhibition to *specific*
salt bridges: the N-AIM residue **D1269** pairs with **R1306** (alpha1 helix)
and **R1450** (alpha6 helix). Type 2B mutations break these, letting A1 engage
GPIbalpha at low/zero force.

The existing equilibrium feature (md_face_destab_score) averages the *global*
AIM->A1 masking network -> it captures 2M (face damage) but dilutes the specific
2B salt-bridge signal. This script tracks ONLY the named autoinhibitory salt
bridges (and a general AIM<->A1 salt-bridge count) to test whether 2B shows
lower occupancy than WT / 2M even at equilibrium.

Read-only over trajectories. Writes output/md_7a6o_saltbridge_features.csv.

Usage:
    python3 scripts/pipeline/extract_aim_saltbridge_features.py \
        --input md_data/7a6o_reference_md/variants \
        --output output/md_7a6o_saltbridge_features.csv
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import MDAnalysis as mda  # noqa: E402
from MDAnalysis.lib.distances import distance_array  # noqa: E402

OFFSET = 1261                       # native_resid = local_resid + OFFSET
SB_CUT = 4.0                        # angstrom, salt-bridge charged-atom cutoff
TAIL_NS = 40.0                      # tail window (40-50 ns) for occupancy

# Charged functional-group atom names
ANION_ATOMS = ("OD1", "OD2", "OE1", "OE2")          # ASP/GLU carboxylate O
CATION_ATOMS = ("NH1", "NH2", "NE", "NZ")           # ARG guanidinium / LYS amine

# Named autoinhibitory salt bridges (native numbering); D1269 is the N-AIM anchor.
NAMED_BRIDGES = {
    "sb_D1269_R1306": (1269, 1306),   # N-AIM <-> alpha1 helix
    "sb_D1269_R1450": (1269, 1450),   # N-AIM <-> alpha6 helix
}

# AIM = N-AIM (1262-1271) + C-AIM (1459-1466); A1 body = 1280-1451
SEL_AIM = "(resid 1:10 or resid 199:205)"
SEL_A1 = "resid 19:190"

LABELS = {
    "WT": "WT",
    "R1306W": "2B", "R1306Q": "2B", "R1308C": "2B", "I1309V": "2B", "S1310F": "2B",
    "W1313C": "2B", "V1314F": "2B", "V1316M": "2B",
    "R1374C": "2M", "R1374H": "2M", "G1324S": "2M",
    "P1337L": "?", "R1341Q": "?", "R1341W": "?",
}


def load(variant_dir: Path):
    tpr = variant_dir / "md_prod.tpr"
    xtc = variant_dir / "prod_concat.xtc"
    if not (tpr.exists() and xtc.exists()):
        return None
    return mda.Universe(str(tpr), str(xtc))


def charged_atoms(u, native_resid: int, kind: str):
    """AtomGroup of the charged-group atoms of one residue (by native resid)."""
    local = native_resid - OFFSET
    names = ANION_ATOMS if kind == "anion" else CATION_ATOMS
    sel = f"resid {local} and name {' '.join(names)}"
    return u.select_atoms(sel)


def pair_occupancy(u, anion_resid: int, cation_resid: int, stride: int = 2):
    """Fraction of tail-window frames with any anion-O within SB_CUT of cation-N."""
    a = charged_atoms(u, anion_resid, "anion")
    c = charged_atoms(u, cation_resid, "cation")
    if len(a) == 0 or len(c) == 0:
        return np.nan  # one side mutated away (e.g. R1306->W/Q removes the cation)
    hits = nframes = 0
    for ts in u.trajectory[::stride]:
        if ts.time / 1000.0 < TAIL_NS:
            continue
        nframes += 1
        if distance_array(a.positions, c.positions).min() < SB_CUT:
            hits += 1
    return hits / max(nframes, 1)


def general_sb_count(u, stride: int = 2):
    """Mean # of AIM<->A1 salt bridges over the tail window (either polarity)."""
    aim_an = u.select_atoms(f"{SEL_AIM} and name {' '.join(ANION_ATOMS)}")
    aim_ca = u.select_atoms(f"{SEL_AIM} and name {' '.join(CATION_ATOMS)}")
    a1_an = u.select_atoms(f"{SEL_A1} and name {' '.join(ANION_ATOMS)}")
    a1_ca = u.select_atoms(f"{SEL_A1} and name {' '.join(CATION_ATOMS)}")
    counts, nframes = [], 0
    for ts in u.trajectory[::stride]:
        if ts.time / 1000.0 < TAIL_NS:
            continue
        nframes += 1
        n = 0
        if len(aim_an) and len(a1_ca):
            n += int((distance_array(aim_an.positions, a1_ca.positions) < SB_CUT).any(axis=1).sum())
        if len(aim_ca) and len(a1_an):
            n += int((distance_array(aim_ca.positions, a1_an.positions) < SB_CUT).any(axis=1).sum())
        counts.append(n)
    return round(float(np.mean(counts)) if counts else np.nan, 2)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", default="md_data/7a6o_reference_md/variants")
    ap.add_argument("--output", default="output/md_7a6o_saltbridge_features.csv")
    ap.add_argument("--stride", type=int, default=2)
    args = ap.parse_args()

    root = Path(args.input)
    variants = sorted([d.name for d in root.iterdir()
                       if d.is_dir() and (d / "md_prod.tpr").exists()])

    rows = []
    for v in variants:
        u = load(root / v)
        row = dict(variant=v, label=LABELS.get(v, "?"))
        for name, (an, cat) in NAMED_BRIDGES.items():
            u.trajectory[0]
            occ = pair_occupancy(u, an, cat, args.stride)
            row[name] = round(occ, 4) if not np.isnan(occ) else np.nan
        u.trajectory[0]
        row["sb_aim_a1_count_tail"] = general_sb_count(u, args.stride)
        rows.append(row)
        print(f"  {v:8s} ({row['label']:7s}) "
              f"D1269-R1306={row['sb_D1269_R1306']} "
              f"D1269-R1450={row['sb_D1269_R1450']} "
              f"AIM-A1_SB={row['sb_aim_a1_count_tail']}")

    df = pd.DataFrame(rows)

    # z-score the general AIM<->A1 salt-bridge count over the labeled ref set
    # (WT + known 2B + 2M). High z = binding face / contacts RETAINED = 2B-compatible
    # (released-but-functional); low z = contacts COLLAPSED = 2M/LOF corroboration.
    # NOTE: per docs/7A6O_SMD_LITERATURE_NOGO_2026-06-24.md this is NOT a 2B-positive
    # axis (2B is a force phenomenon, equilibrium-invisible). It is a 2M corroborator
    # and a 2B-protective guard for the joint criterion. Thresholds need calibration
    # (n: WT=1, 2B=8, 2M=3 single-replica).
    ref = df[df.label.isin(["WT", "2B", "2M"])]["sb_aim_a1_count_tail"]
    mu, sd = ref.mean(), ref.std()
    df["aim_sb_retained_z"] = ((df["sb_aim_a1_count_tail"] - mu) / sd).round(3) \
        if sd and not np.isnan(sd) else np.nan

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\nWritten: {out}  (ref mean={mu:.3f} sd={sd:.3f})\n")

    # group separation check
    print("=== group means (known labels) ===")
    for col in ["sb_D1269_R1306", "sb_D1269_R1450", "sb_aim_a1_count_tail"]:
        print(f"\n{col}:")
        for lab in ["WT", "2B", "2M"]:
            vals = df[df.label == lab][col].dropna()
            if len(vals):
                print(f"  {lab:3s} n={len(vals)} mean={vals.mean():.3f} "
                      f"range=[{vals.min():.3f},{vals.max():.3f}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

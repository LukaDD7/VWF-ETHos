#!/usr/bin/env python3
"""Unify the scattered VWF Type-2 label sources into one clean table.

Motivation (2026-06-24): the classifier validation runs on VWF_Alpha_Matrix
(100 variants, only 12 2B / 25 2M). Much more labeled data exists but in messy /
inconsistent formats. This script normalizes them into a single schema and
reports how many 2B/2M variants we have, and -- critically -- how many already
have computed features (immediately usable) vs need the feature pipeline run.

Sources:
  - data/raw_tables/VWF_Type2_AF3_Reference_Table.csv  (cleanest; the AF3 100)
  - data/processed/master_type1_type2.csv              (336; +103 Type2B)
  - data/raw_tables/2B型突变.xlsx                       (lit. Table S3, ~400 2B, 3-letter)

Output: output/expanded_label_set.csv with columns
  aa_change, wt_aa, position, mut_aa, subtype, domain, sources, has_features

Usage:
    python3 scripts/pipeline/build_expanded_label_set.py
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}

# normalize any subtype spelling -> canonical
SUBTYPE_MAP = {
    "type2a": "2A", "2a": "2A", "type2b": "2B", "2b": "2B",
    "type2m": "2M", "2m": "2M", "type2n": "2N", "2n": "2N",
    "type1": "1", "1": "1", "type3": "3", "3": "3",
    "wt_control": "WT", "wt": "WT",
}


def norm_subtype(s):
    if not isinstance(s, str):
        return None
    return SUBTYPE_MAP.get(s.strip().lower())


def parse_3letter(aa):
    """Tyr1258Cys / p.Tyr1258Cys -> ('Y', 1258, 'C') ; None if not parseable."""
    s = str(aa).strip()
    if s.startswith("p."):
        s = s[2:]
    m = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", s)
    if not m:
        return None
    wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
    if wt not in THREE_TO_ONE or mut not in THREE_TO_ONE:
        return None
    return THREE_TO_ONE[wt], pos, THREE_TO_ONE[mut]


def parse_1letter(aa):
    """S2775C -> ('S', 2775, 'C')."""
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", str(aa).strip())
    return (m.group(1), int(m.group(2)), m.group(3)) if m else None


def from_af3_ref(path):
    d = pd.read_csv(path)
    rows = []
    for _, r in d.iterrows():
        p = parse_1letter(r.get("AA_Change", ""))
        st = norm_subtype(r.get("Type2_Subtype"))
        if p and st:
            rows.append((f"{p[0]}{p[1]}{p[2]}", p[0], p[1], p[2], st,
                         r.get("VWF_Domain", ""), "af3_ref"))
    return rows


def from_master(path):
    d = pd.read_csv(path)
    rows = []
    for _, r in d.iterrows():
        try:
            wt, pos, mut = r["WT_AA"], int(r["Position"]), r["Mut_AA"]
        except (ValueError, TypeError, KeyError):
            continue
        st = norm_subtype(r.get("Subtype"))
        if st:
            rows.append((f"{wt}{pos}{mut}", wt, pos, mut, st,
                         r.get("Domain_Name", ""), "master"))
    return rows


def from_2b_xlsx(path):
    raw = pd.read_excel(path, header=None)
    rows = []
    for _, r in raw.iloc[2:].iterrows():     # data starts row 2; col0=aa, col3=domain
        if isinstance(r[0], str) and r[0].strip().lower().startswith("reference"):
            break                            # bibliography after this row -> stop
        p = parse_3letter(r[0])
        if p:
            rows.append((f"{p[0]}{p[1]}{p[2]}", p[0], p[1], p[2], "2B",
                         r[3] if isinstance(r[3], str) else "", "lit_2B_S3"))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--af3", default="data/raw_tables/VWF_Type2_AF3_Reference_Table.csv")
    ap.add_argument("--master", default="data/processed/master_type1_type2.csv")
    ap.add_argument("--xlsx", default="data/raw_tables/2B型突变.xlsx")
    ap.add_argument("--matrix", default="data/processed/VWF_Alpha_Matrix.parquet")
    ap.add_argument("--output", default="output/expanded_label_set.csv")
    args = ap.parse_args()

    rows = []
    for fn, src in [(from_af3_ref, args.af3), (from_master, args.master),
                    (from_2b_xlsx, args.xlsx)]:
        if Path(src).exists():
            r = fn(src)
            rows += r
            print(f"{Path(src).name:42s} -> {len(r):4d} parsed labeled variants")

    df = pd.DataFrame(rows, columns=["aa_change", "wt_aa", "position", "mut_aa",
                                     "subtype", "domain", "source"])

    # dedupe by aa_change: merge sources, flag label conflicts
    recs = []
    for aa, g in df.groupby("aa_change"):
        subs = sorted(set(g.subtype))
        recs.append({
            "aa_change": aa, "wt_aa": g.wt_aa.iloc[0],
            "position": int(g.position.iloc[0]), "mut_aa": g.mut_aa.iloc[0],
            "subtype": subs[0] if len(subs) == 1 else "|".join(subs),
            "conflict": len(subs) > 1,
            "domain": next((d for d in g.domain if isinstance(d, str) and d), ""),
            "sources": ",".join(sorted(set(g.source))),
        })
    uni = pd.DataFrame(recs)

    # which already have features (in the matrix)?
    have = set()
    if Path(args.matrix).exists():
        have = set(pd.read_parquet(args.matrix)["aa_change"].astype(str))
    uni["has_features"] = uni["aa_change"].isin(have)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    uni.to_csv(out, index=False)

    print(f"\nWritten: {out}  ({len(uni)} unique variants)")
    print("\n=== unique variants by subtype (after dedupe) ===")
    print(uni["subtype"].value_counts().to_string())
    print(f"\nlabel conflicts (multi-source disagreement): {int(uni['conflict'].sum())}")

    print("\n=== 2B / 2M: usable now (have features) vs need feature pipeline ===")
    for st in ["2B", "2M"]:
        sub = uni[uni.subtype == st]
        print(f"  {st}: total={len(sub)}  have_features={int(sub.has_features.sum())}  "
              f"need_features={int((~sub.has_features).sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

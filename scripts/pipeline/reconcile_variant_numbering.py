#!/usr/bin/env python3
"""Reconcile VWF variant numbering (mature <-> pre-pro P04275) against the WT FASTA.

Some literature/source tables number VWF variants in *mature* protein coordinates
(residue 1 = first mature residue = pre-pro residue 764), while this pipeline uses
*pre-pro* P04275 numbering. A variant whose labeled WT residue does not match the
FASTA at `position` but DOES match at `position + 763` is mature-numbered and is
shifted to pre-pro.

For each shifted variant we also check whether the resulting (pre-pro position,
mut) already carries a label elsewhere -> duplicate (same subtype) or conflict
(different subtype). Duplicates and conflicts are excluded from the runnable set
and reported, so we never spend GPU on a wrong residue or a contested label.

Outputs:
  - <out>            runnable reconciled variants (prepro_ok + clean shifted)
  - <out>.report.csv every input variant with status

Usage:
    python3 scripts/pipeline/reconcile_variant_numbering.py \
        --variants output/new_2b_for_boltz.csv \
        --fasta Proteo-Structure-Pipeline/structures/wt/VWF_P04275_WT.fasta \
        --labels output/expanded_label_set.csv \
        --out output/new_2b_for_boltz_reconciled.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

MATURE_OFFSET = 763  # pre-pro = mature + 763


def read_fasta(path: Path) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines()
                   if not l.startswith(">"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--labels", default=None,
                    help="expanded_label_set.csv for duplicate/conflict check")
    ap.add_argument("--out", required=True)
    ap.add_argument("--offset", type=int, default=MATURE_OFFSET)
    args = ap.parse_args()

    seq = read_fasta(Path(args.fasta))
    df = pd.read_csv(args.variants)

    # existing (prepro pos, mut) -> subtype, for dup/conflict detection
    known: dict[tuple[int, str], str] = {}
    if args.labels and Path(args.labels).exists():
        lab = pd.read_csv(args.labels)
        for r in lab.itertuples():
            if str(r.subtype) in ("2A", "2B", "2M", "2N"):
                known[(int(r.position), str(r.mut_aa))] = str(r.subtype)

    def at(pos):
        return seq[pos - 1] if 1 <= pos <= len(seq) else None

    recs = []
    for r in df.itertuples():
        wt, pos, mut = str(r.wt_aa), int(r.position), str(r.mut_aa)
        sub = str(getattr(r, "subtype", "2B"))
        if at(pos) == wt:
            status, newpos = "prepro_ok", pos
        elif at(pos + args.offset) == wt:
            newpos = pos + args.offset
            coll = known.get((newpos, mut))
            if coll is None:
                status = "mature_shifted"
            elif coll == sub:
                status = f"duplicate_{coll}"
            else:
                status = f"conflict_{coll}_vs_{sub}"
        else:
            status, newpos = "unresolved", pos
        recs.append(dict(aa_change=f"{wt}{newpos}{mut}", wt_aa=wt, position=newpos,
                         mut_aa=mut, subtype=sub, orig_position=pos,
                         orig_aa_change=getattr(r, "aa_change", f"{wt}{pos}{mut}"),
                         status=status, domain=getattr(r, "domain", "")))
    rep = pd.DataFrame(recs)

    runnable = rep[rep.status.isin(["prepro_ok", "mature_shifted"])][
        ["aa_change", "wt_aa", "position", "mut_aa", "subtype", "domain"]].copy()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    runnable.to_csv(args.out, index=False)
    rep.to_csv(args.out + ".report.csv", index=False)

    print(f"input variants: {len(rep)}")
    print(rep.status.apply(lambda s: s.split('_')[0] if s.startswith(('duplicate', 'conflict'))
                           else s).value_counts().to_string())
    print(f"\nrunnable (prepro_ok + mature_shifted): {len(runnable)} -> {args.out}")
    flagged = rep[~rep.status.isin(["prepro_ok", "mature_shifted"])]
    if len(flagged):
        print(f"\n=== flagged (excluded from run), n={len(flagged)} ===")
        print(flagged[["orig_aa_change", "aa_change", "status"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
generate_a1_dp_d3_yamls.py
===========================
Generate Boltz-2 YAML input files for VWF A1 + D'D3 autoinhibition complex.

This produces structures for studying the D'D3-A1 autoinhibitory interface —
the interaction that keeps VWF A1 in an inactive state under quiescent
conditions. Type 2B (GOF) mutations are thought to destabilise this interface,
promoting spontaneous A1-GPIbα binding.

Sequence composition (single chain, VWF D' + D3 + A1):
  Chain A: VWF D'  domain (UniProt P04275, aa 764-865,  102 aa)
  Chain A: VWF D3  domain (UniProt P04275, aa 866-1233, 368 aa)
  Chain A: VWF A1  domain (UniProt P04275, aa 1268-1466, 199 aa)
  Total: 669 aa

Note: we model D'D3 and A1 as a single chain because in the autoinhibited
state the D'D3 module sits adjacent to A1, forming the autoinhibitory
interface. Type 2B mutations in A1 disrupt this interface → spontaneous
binding to platelet GPIbα.

Usage:
  python3 generate_a1_dp_d3_yamls.py --variants ../../data/processed/master_type1_type2.csv
  python3 generate_a1_dp_d3_yamls.py --dry-run
"""

import argparse
import sys
import pandas as pd
from pathlib import Path

# =============================================================================
# Sequences — VWF D' + D3 + A1 (continuous chain)
# Source: UniProt P04275 full-length sequence (2813 aa), validated against
#         generate_a1_gpiba_yamls.py (A1 1268-1466)
# =============================================================================

# VWF D' domain: P04275 aa 764-865, 102 aa
VWF_DP_WT = (
    "SLSCRPPMVKLVCPADNLRAEGLECTKTCQNYDLECMSMGCVSGCLCPPG"
    "MVRHENRCVALERCPCFHQGKEYAPGETVKIGCNTCVCQDRKWNCTDHVC"
    "DA"
)
assert len(VWF_DP_WT) == 102, f"D' length should be 102, got {len(VWF_DP_WT)}"
assert VWF_DP_WT[0] == "S"

# VWF D3 domain: P04275 aa 866-1233, 368 aa
VWF_D3_WT = (
    "TCSTIGMAHYLTFDGLKYLFPGECQYVLVQDYCGSNPGTFRILVGNKGCS"
    "HPSVKCKKRVTILVEGGEIELFDGEVNVKRPMKDETHFEVVESGRYIILL"
    "LGKALSVVWDRHLSISVVLKQTYQEKVCGLCGNFDGIQNNDLTSSNLQVE"
    "EDPVDFGNSWKVSSQCADTRKVPLDSSPATCHNNIMKQTMVDSSCRILTS"
    "DVFQDCNKLVDPEPYLDVCIYDTCSCESIGDCACFCDTIAAYAHVCAQHG"
    "KVVTWRTATLCPQSCEERNLRENGYECEWRYNSCAPACQVTCQHPEPLAC"
    "PVQCVEGCHAHCPPGKILDELLQTCVDPEDCPVCEVAGRRFASGKKVTLN"
    "PSDPEHCQICHCDVVNLT"
)
assert len(VWF_D3_WT) == 368, f"D3 length should be 368, got {len(VWF_D3_WT)}"

# VWF A1 domain: P04275 aa 1268-1466, 199 aa
# Matches generate_a1_gpiba_yamls.py exactly (same WT sequence)
VWF_A1_WT = (
    "HDFYCSRLLDLVFLLDGSSRLSEAEFEVLKAFVVDMMERLRISQKWVRVAVVEYHDGSH"
    "AYIGLKDRKRPSELRRIASQVKYAGSQVASTSEVLKYTLFQIFSKIDRPEASRITLLLM"
    "ASQEPQRMSRNFVRYVQGLKKKKVIVIPVGIGPHANLKQIRLIEKQAPENKAFVLSSVDE"
    "LEQQRDEIVSYLCDLAPEAPP"
)
assert len(VWF_A1_WT) == 199, f"A1 length should be 199, got {len(VWF_A1_WT)}"
assert VWF_A1_WT[38] == "R", f"A1 position 1306 should be R, got {VWF_A1_WT[38]}"

DP_START = 764   # UniProt start of D'
A1_START = 1268  # UniProt start of A1 (consistent with existing A1+GPIbα pipeline)


# =============================================================================
# Helpers
# =============================================================================

def apply_mutation(seq: str, uniprot_pos: int, wt_aa: str, mut_aa: str,
                   seq_start: int) -> str:
    """Apply a point mutation to a sequence."""
    idx = uniprot_pos - seq_start
    if idx < 0 or idx >= len(seq):
        raise ValueError(
            f"Position {uniprot_pos} out of range [{seq_start}, {seq_start+len(seq)-1}]"
        )
    actual = seq[idx]
    if actual != wt_aa:
        raise ValueError(
            f"At UniProt {uniprot_pos}: expected {wt_aa}, found {actual}"
        )
    return seq[:idx] + mut_aa + seq[idx+1:]


def make_yaml(chain_a_seq: str, job_name: str) -> str:
    """Generate a Boltz-2 v2 YAML string for a single-chain protein.

    MSA is 'empty' for offline GPU servers.
    """
    return f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {chain_a_seq}
      msa: empty
# job: {job_name}
# chain A: VWF D' + D3 + A1 (UniProt P04275 aa 764-1466)
# Purpose: study D'D3-A1 autoinhibitory interface
"""


def load_a1_variants(csv_path: str) -> pd.DataFrame:
    """Load and filter A1-domain variants from master CSV."""
    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]

    a1 = df[df["Domain_Name"].str.upper() == "A1"].copy()

    a1["_key"] = (
        a1["Variant_ID"].astype(str) + "_"
        + a1["Position"].astype(str) + "_"
        + a1["Mut_AA"].astype(str)
    )
    a1 = a1.drop_duplicates(subset="_key")

    # Filter to A1 domain range (1268-1466)
    a1 = a1[(a1["Position"] >= A1_START) &
             (a1["Position"] <= A1_START + len(VWF_A1_WT) - 1)]

    a1 = a1[~a1["Mut_AA"].isin(["*", "X"])]
    a1 = a1[a1["WT_AA"] != a1["Mut_AA"]]

    return a1.reset_index(drop=True)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--variants",
        default="../../data/processed/master_type1_type2.csv",
        help="Path to master variant CSV"
    )
    parser.add_argument(
        "--out-dir",
        default="../../output/boltz2_a1_dp_d3",
        help="Output directory for YAML files"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print plan without writing files"
    )
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    variants_path = (script_dir / args.variants).resolve()
    out_dir = (script_dir / args.out_dir).resolve()

    if not variants_path.exists():
        print(f"ERROR: variants file not found: {variants_path}", file=sys.stderr)
        sys.exit(1)

    full_wt = VWF_DP_WT + VWF_D3_WT + VWF_A1_WT
    print(f"VWF D'  : {len(VWF_DP_WT)} aa (P04275 aa 764-865)")
    print(f"VWF D3  : {len(VWF_D3_WT)} aa (P04275 aa 866-1233)")
    print(f"VWF A1  : {len(VWF_A1_WT)} aa (P04275 aa 1268-1466)")
    print(f"Total   : {len(full_wt)} aa")
    print(f"A1 R1306: chain[{len(VWF_DP_WT)+len(VWF_D3_WT)+38}]={full_wt[len(VWF_DP_WT)+len(VWF_D3_WT)+38]}")
    print()

    variants = load_a1_variants(str(variants_path))
    print(f"A1 variants loaded: {len(variants)}")

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    # WT reference YAML
    wt_name = "VWF_WT_dp_d3_a1"
    wt_yaml = make_yaml(full_wt, wt_name)
    if args.dry_run:
        print(f"[DRY RUN] Would write: {wt_name}.yaml")
    else:
        wt_path = out_dir / f"{wt_name}.yaml"
        wt_path.write_text(wt_yaml)
        print(f"Written: {wt_path.name}")

    # Mutant YAMLs
    written = 0
    skipped = 0
    for _, row in variants.iterrows():
        var_id = row["Variant_ID"]
        pos = int(row["Position"])
        wt_aa = str(row["WT_AA"])
        mut_aa = str(row["Mut_AA"])

        try:
            # Mutation is in A1 region only
            mut_seq = apply_mutation(full_wt, pos, wt_aa, mut_aa, A1_START)
        except ValueError as e:
            print(f"  [SKIP] {var_id}: {e}")
            skipped += 1
            continue

        job_name = f"VWF_{var_id}_dp_d3_a1"
        yaml_content = make_yaml(mut_seq, job_name)

        if args.dry_run:
            print(f"[DRY RUN] {var_id}: {wt_aa}{pos}{mut_aa}")
        else:
            yaml_path = out_dir / f"{job_name}.yaml"
            yaml_path.write_text(yaml_content)
            written += 1

    if args.dry_run:
        print(f"\n=== DRY RUN SUMMARY ===")
        print(f"Would write: {len(variants) - skipped + 1} YAML files")
        print(f"Skipped: {skipped}")
    else:
        print(f"\nDone: {written + 1} YAML files -> {out_dir}")
        print(f"  1 WT reference + {written} mutants")
        if skipped:
            print(f"  {skipped} skipped (sequence mismatch)")
        print()
        print("=== Run on GPU server ===")
        print(f"boltz predict {out_dir}/ \\")
        print(f"    --out_dir {out_dir}_results/ \\")
        print(f"    --accelerator gpu --devices 4 \\")
        print(f"    --recycling_steps 3 --diffusion_samples 5")


if __name__ == "__main__":
    main()
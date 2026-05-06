#!/usr/bin/env python3
"""
generate_a1_gpiba_yamls.py
==========================
Generate Boltz-2 YAML input files for VWF A1 + GPIbα complex prediction.

For each variant in master_type1_type2.csv (A1 domain), produces:
  output/boltz2_a1_gpiba/<variant_id>.yaml

Plus one WT reference:
  output/boltz2_a1_gpiba/VWF_WT_vs_GPIb_alpha.yaml

YAML format follows Boltz-2 v2 spec (docs/prediction.md).
MSA is set to 'empty' for offline GPU servers; pass --use-msa-server
as a note in the run command if the server has internet access.

Sequences:
  Chain A: VWF A1 domain (UniProt P04275, residues 1268-1466, 199 aa)
  Chain B: GPIbα extracellular domain (UniProt P07359, residues 1-290, 290 aa)
           = N-terminal LRR domain that directly contacts VWF A1

Usage:
  python3 generate_a1_gpiba_yamls.py
  python3 generate_a1_gpiba_yamls.py --variants ../../data/processed/master_type1_type2.csv
  python3 generate_a1_gpiba_yamls.py --out-dir /path/to/output --dry-run
"""

import argparse
import os
import sys
import pandas as pd
from pathlib import Path

# =============================================================================
# Sequences
# =============================================================================

# VWF A1 domain WT sequence
# Source: UniProt P04275, residues 1268-1466 (1-indexed, inclusive), 199 aa
# Extracted from VWF_P04275_WT.fasta; validated: pos 1306 = R (R1306)
VWF_A1_WT = (
    "HDFYCSRLLDLVFLLDGSSRLSEAEFEVLKAFVVDMMERLRISQKWVRVAVVEYHDGSH"
    "AYIGLKDRKRPSELRRIASQVKYAGSQVASTSEVLKYTLFQIFSKIDRPEASRITLLLM"
    "ASQEPQRMSRNFVRYVQGLKKKKVIVIPVGIGPHANLKQIRLIEKQAPENKAFVLSSVDE"
    "LEQQRDEIVSYLCDLAPEAPP"
)
VWF_A1_UNIPROT_START = 1268  # 1-indexed

# GPIbα extracellular N-terminal domain (VWF-binding LRR domain)
# Source: UniProt P07359, residues 1-290 (full extracellular VWF-binding domain)
# This region contains the leucine-rich repeats (LRR) that contact VWF A1
# The 1M10 crystal structure uses residues 1-267; we use 1-290 for completeness
# Deglycosylation mutations in 1M10 (N37Q/N175Q/M255V) are crystallization
# artifacts — use WT GPIbα sequence here for physiological relevance
GPIBA_WT = (
    "MPLLLLLLLLPSPLHPHPICEVSKVASHLEVNCDKRNLTALPPDLPKDTTILHLSENLLYTFSLATLMPYTRLTQLNLDRCELTKLQVD"
    "GTLPVLGTLDLSHNQLQSLPLLGQTLPALTVLDVSFNRLTSLPLGALRGLGELQELYLKGNELKTLPPGLLTPTPKLEKLSLANNNLT"
    "ELPAGLLNGLENLDTLLLQENSLYTIPKGFFGSHLLPFAFLHGNPWLCNCEILYFRRWLQDNAENVYVWKQGVDVKAMTSNVASVQCD"
    "NSDKFPVYKYPGKGCPTLGDEGDTD"
)
# Verify lengths
assert len(VWF_A1_WT) == 199, f"A1 WT length should be 199, got {len(VWF_A1_WT)}"
assert len(GPIBA_WT) == 290, f"GPIbα length should be 290, got {len(GPIBA_WT)}"


# =============================================================================
# Helpers
# =============================================================================

def apply_mutation(seq: str, uniprot_pos: int, wt_aa: str, mut_aa: str,
                   seq_start: int) -> str:
    """Apply a point mutation to a sequence.

    Args:
        seq: WT sequence string
        uniprot_pos: 1-indexed UniProt position
        wt_aa: expected WT amino acid (for validation)
        mut_aa: mutant amino acid
        seq_start: UniProt position of seq[0]

    Returns:
        mutated sequence string

    Raises:
        ValueError: if WT AA doesn't match or position out of range
    """
    idx = uniprot_pos - seq_start  # 0-indexed
    if idx < 0 or idx >= len(seq):
        raise ValueError(
            f"Position {uniprot_pos} out of range [{seq_start}, {seq_start+len(seq)-1}]"
        )
    actual = seq[idx]
    if actual != wt_aa:
        raise ValueError(
            f"At UniProt {uniprot_pos}: expected {wt_aa}, found {actual} in sequence"
        )
    return seq[:idx] + mut_aa + seq[idx+1:]


def make_yaml(chain_a_seq: str, chain_b_seq: str, job_name: str) -> str:
    """Generate a Boltz-2 v2 YAML string for a two-protein complex.

    MSA is set to 'empty' (single-sequence mode) for offline servers.
    If the GPU server has internet, add --use_msa_server to boltz predict.
    """
    return f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {chain_a_seq}
      msa: empty
  - protein:
      id: B
      sequence: {chain_b_seq}
      msa: empty
# job: {job_name}
# chain A: VWF A1 domain (UniProt P04275 aa 1268-1466)
# chain B: GPIbα extracellular domain (UniProt P07359 aa 1-290)
# no properties.affinity: protein-protein complex, use iPTM as proxy
"""


def load_variants(csv_path: str) -> pd.DataFrame:
    """Load and filter A1-domain variants from master CSV."""
    df = pd.read_csv(csv_path)
    # Normalize column names
    df.columns = [c.strip() for c in df.columns]

    # Filter A1 domain
    a1 = df[df["Domain_Name"].str.upper() == "A1"].copy()

    # Deduplicate by (Variant_ID, Position, Mut_AA)
    a1["_key"] = (
        a1["Variant_ID"].astype(str) + "_"
        + a1["Position"].astype(str) + "_"
        + a1["Mut_AA"].astype(str)
    )
    a1 = a1.drop_duplicates(subset="_key")

    # Filter to PDB/sequence range
    a1 = a1[(a1["Position"] >= VWF_A1_UNIPROT_START) &
             (a1["Position"] <= VWF_A1_UNIPROT_START + len(VWF_A1_WT) - 1)]

    # Remove stop codons and synonymous
    a1 = a1[~a1["Mut_AA"].isin(["*", "X"])]
    a1 = a1[a1["WT_AA"] != a1["Mut_AA"]]

    return a1.reset_index(drop=True)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--variants", default="../../data/processed/master_type1_type2.csv",
                        help="Path to master_type1_type2.csv")
    parser.add_argument("--out-dir", default="../../output/boltz2_a1_gpiba",
                        help="Output directory for YAML files")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print plan without writing files")
    args = parser.parse_args()

    # Resolve paths relative to this script's location
    script_dir = Path(__file__).parent
    variants_path = (script_dir / args.variants).resolve()
    out_dir = (script_dir / args.out_dir).resolve()

    if not variants_path.exists():
        print(f"ERROR: variants file not found: {variants_path}", file=sys.stderr)
        sys.exit(1)

    # Validate A1 WT sequence
    a1_wt = VWF_A1_WT.replace("\n", "")
    gpiba = GPIBA_WT.replace("\n", "")

    # Sanity check: R at position 1306
    r1306_idx = 1306 - VWF_A1_UNIPROT_START  # = 38
    assert a1_wt[r1306_idx] == "R", f"Position 1306 should be R, got {a1_wt[r1306_idx]}"

    print(f"VWF A1 WT: {len(a1_wt)} aa (UniProt 1268-1466), pos 1306 = {a1_wt[r1306_idx]}")
    print(f"GPIbα:     {len(gpiba)} aa (UniProt P07359 aa 1-290)")
    print()

    # Load variants
    variants = load_variants(str(variants_path))
    print(f"A1 variants loaded: {len(variants)}")
    print()

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    # WT reference YAML
    wt_name = "VWF_WT_vs_GPIb_alpha"
    wt_yaml = make_yaml(a1_wt, gpiba, wt_name)
    if args.dry_run:
        print(f"[DRY RUN] Would write: {wt_name}.yaml (WT reference)")
        print(wt_yaml[:300] + "...\n")
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
            mut_seq = apply_mutation(a1_wt, pos, wt_aa, mut_aa, VWF_A1_UNIPROT_START)
        except ValueError as e:
            print(f"  [SKIP] {var_id}: {e}")
            skipped += 1
            continue

        job_name = f"{var_id}_vs_GPIb_alpha"
        yaml_content = make_yaml(mut_seq, gpiba, job_name)

        if args.dry_run:
            print(f"[DRY RUN] {var_id}: {wt_aa}{pos}{mut_aa} → A1[{pos-VWF_A1_UNIPROT_START}] {wt_aa}→{mut_aa}")
        else:
            yaml_path = out_dir / f"{job_name}.yaml"
            yaml_path.write_text(yaml_content)
            written += 1

    if args.dry_run:
        print(f"\n=== DRY RUN SUMMARY ===")
        print(f"Would write: {len(variants) - skipped + 1} YAML files (1 WT + {len(variants) - skipped} mutants)")
        print(f"Skipped: {skipped}")
    else:
        print(f"\nDone: {written + 1} YAML files written to {out_dir}")
        print(f"  1 WT reference + {written} mutants")
        if skipped:
            print(f"  {skipped} variants skipped (sequence mismatch)")
        print()
        print("=== Run on GPU server ===")
        print(f"boltz predict output/boltz2_a1_gpiba/ \\")
        print(f"    --out_dir output/boltz2_a1_gpiba_results/ \\")
        print(f"    --accelerator gpu --devices 4 \\")
        print(f"    --recycling_steps 3 --diffusion_samples 5 \\")
        print(f"    --use_msa_server   # add if server has internet, else remove")


if __name__ == "__main__":
    main()

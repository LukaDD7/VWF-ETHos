#!/usr/bin/env python3
"""Clean a PDB for GROMACS pdb2gmx by keeping one protein chain and standard amino acids."""
from __future__ import annotations

import argparse
from pathlib import Path


STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    src = Path(args.pdb)
    dst = Path(args.out)
    if not src.is_file():
        raise SystemExit(f"Input PDB not found: {src}")

    kept = 0
    residues: set[str] = set()
    with dst.open("w", encoding="ascii") as out:
        for line in src.read_text(encoding="ascii", errors="ignore").splitlines():
            if not line.startswith("ATOM  "):
                continue
            chain = line[21].strip()
            resname = line[17:20].strip().upper()
            if chain == args.chain and resname in STANDARD_AA:
                out.write(line + "\n")
                kept += 1
                residues.add(line[22:27])
        out.write("TER\n")
        out.write("END\n")

    print(f"Kept {kept} atoms from chain {args.chain} in {len(residues)} residues")
    print(f"Clean PDB written to {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

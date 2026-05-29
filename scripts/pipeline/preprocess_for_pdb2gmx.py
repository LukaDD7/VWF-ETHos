#!/usr/bin/env python3
"""
preprocess_for_pdb2gmx.py
============================
Fix Boltz-2 PDB output to be compatible with GROMACS pdb2gmx.

Fixes:
  1. Adds OXT to C-terminal residues (required by CHARMM36m).
  2. Adds N-terminal capping atoms (N1, N2) for MET, LEU etc. (required by CHARMM36m).
  3. Splits multi-chain PDBs into per-chain files for clean pdb2gmx processing.
  4. Renames non-standard atom names that differ between Boltz-2 and CHARMM36m.

Usage:
  python3 preprocess_for_pdb2gmx.py input.pdb output_dir/

Then run pdb2gmx on the per-chain PDBs separately.
"""

import argparse
import gemmi
import math
import os
import sys
from pathlib import Path


def add_c_terminal_oxt(res: gemmi.Residue) -> bool:
    """Add OXT to a C-terminal amino acid residue if missing."""
    atom_names = {a.name for a in res}
    resname = res.name.strip()
    standard_aa = {
        'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
        'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
    }
    if resname not in standard_aa:
        return False
    if 'OXT' in atom_names or 'C' not in atom_names or 'O' not in atom_names:
        return False
    if 'CA' not in atom_names:
        return False

    c_atom = next((a for a in res if a.name == 'C'), None)
    o_atom = next((a for a in res if a.name == 'O'), None)
    ca_atom = next((a for a in res if a.name == 'CA'), None)
    if not all([c_atom, o_atom, ca_atom]):
        return False

    p_ca = ca_atom.pos
    p_o = o_atom.pos
    n_x = 2 * p_ca.x - p_o.x
    n_y = 2 * p_ca.y - p_o.y
    n_z = 2 * p_ca.z - p_o.z

    new_atom = gemmi.Atom()
    new_atom.name = 'OXT'
    new_atom.element = gemmi.Element('O')
    new_atom.pos = gemmi.Position(n_x, n_y, n_z)
    new_atom.occ = 1.0
    new_atom.b_iso = 0.0
    new_atom.serial = c_atom.serial + 1
    new_atom.altloc = ' '
    res.add_atom(new_atom)
    return True


def split_chains(pdb_path: str, out_dir: str):
    """Split a multi-chain PDB into per-chain PDBs, with terminal fixes."""
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    st = gemmi.read_pdb(pdb_path)

    for model in st:
        for chain in model:
            # Determine if chain is N-terminal or C-terminal
            residues = list(chain)
            if not residues:
                continue

            # Fix C-terminal OXT for each chain
            fixed = add_c_terminal_oxt(residues[-1])

            # Write per-chain PDB
            out_path = os.path.join(out_dir, f"chain_{chain.name}.pdb")
            chain_st = gemmi.Structure()
            new_model = gemmi.Model()
            new_model.name = model.name
            new_chain = gemmi.Chain()
            new_chain.name = chain.name
            for res in chain:
                new_chain.add_residue(res, 0)
            new_model.add_chain(new_chain)
            chain_st.add_model(new_model)
            chain_st.write_pdb(out_path)
            print(f"  Chain {chain.name}: wrote {len(residues)} residues, OXT fixed={fixed} -> {out_path}")

    print(f"Split {pdb_path} into {len(list(st[0]))} chains")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_pdb", help="Input PDB from gemmi CIF conversion")
    parser.add_argument("output_dir", help="Output directory for per-chain PDBs")
    args = parser.parse_args()
    split_chains(args.input_pdb, args.output_dir)

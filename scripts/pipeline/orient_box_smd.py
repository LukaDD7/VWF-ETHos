#!/usr/bin/env python3
"""Orient a protein .gro so the N-anchor -> C-anchor vector points along +z.

Used by prep_7a6o_smd.sh to set up steered-MD (SMD) of AIM unfolding: after
orientation the two AIM-flanking anchors are separated along z, so the pull
coordinate (pull-coord1-dim = N N Y) extends them apart in an elongated z-box.

Dependency-light: numpy only (no MDAnalysis/gemmi needed on the A40 env).
Reads native GROMACS residue numbering (e.g. 1262..1466 for the 7A6O construct).

Usage:
    python3 orient_box_smd.py --in prot_whole.gro --out prot_oriented.gro \
        --n-resid 1262-1264 --c-resid 1464-1466 --name CA
"""
from __future__ import annotations
import argparse
import numpy as np


def parse_range(s: str):
    a, b = s.split("-")
    return int(a), int(b)


def read_gro(path: str):
    with open(path) as f:
        lines = f.read().splitlines()
    title = lines[0]
    natoms = int(lines[1].strip())
    atoms = []
    for ln in lines[2:2 + natoms]:
        resid = int(ln[0:5])
        resname = ln[5:10].strip()
        atomname = ln[10:15].strip()
        atomnum = ln[15:20]
        x, y, z = float(ln[20:28]), float(ln[28:36]), float(ln[36:44])
        atoms.append([resid, resname, atomname, atomnum, x, y, z])
    box = lines[2 + natoms] if len(lines) > 2 + natoms else "   0.0   0.0   0.0"
    return title, atoms, box


def write_gro(path: str, title: str, atoms, box: str):
    with open(path, "w") as f:
        f.write(title + "\n")
        f.write(f"{len(atoms):5d}\n")
        for resid, resname, atomname, atomnum, x, y, z in atoms:
            f.write(f"{resid % 100000:5d}{resname:<5s}{atomname:>5s}{atomnum:>5s}"
                    f"{x:8.3f}{y:8.3f}{z:8.3f}\n")
        f.write(box + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-resid", default="1262-1264", help="N-anchor residue range (native)")
    ap.add_argument("--c-resid", default="1464-1466", help="C-anchor residue range (native)")
    ap.add_argument("--name", default="CA", help="atom name used for anchor COM")
    args = ap.parse_args()

    title, atoms, box = read_gro(args.inp)
    coords = np.array([[a[4], a[5], a[6]] for a in atoms])
    resids = np.array([a[0] for a in atoms])
    names = np.array([a[2] for a in atoms])
    n0, n1 = parse_range(args.n_resid)
    c0, c1 = parse_range(args.c_resid)

    nmask = (resids >= n0) & (resids <= n1) & (names == args.name)
    cmask = (resids >= c0) & (resids <= c1) & (names == args.name)
    if nmask.sum() == 0 or cmask.sum() == 0:
        raise SystemExit(f"anchor selection empty: N={nmask.sum()} C={cmask.sum()} "
                         f"(check --n-resid/--c-resid/--name and gro numbering)")
    nc = coords[nmask].mean(0)
    cc = coords[cmask].mean(0)
    v = cc - nc
    v = v / np.linalg.norm(v)
    z = np.array([0.0, 0.0, 1.0])
    axis = np.cross(v, z)
    s = np.linalg.norm(axis)
    if s > 1e-6:
        axis /= s
        c = float(np.dot(v, z))
        K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
        R = np.eye(3) + s * K + (1 - c) * (K @ K)
    else:
        R = np.eye(3)

    center = coords.mean(0)
    new = (coords - center) @ R.T + center
    for a, p in zip(atoms, new):
        a[4], a[5], a[6] = float(p[0]), float(p[1]), float(p[2])
    write_gro(args.out, title, atoms, box)

    # report (after orientation; editconf will set the final box)
    znew = new
    nz = znew[nmask][:, 2].mean()
    cz = znew[cmask][:, 2].mean()
    ext = znew.max(0) - znew.min(0)
    print(f"oriented: anchor z-separation = {cz - nz:+.2f} nm  "
          f"protein extent x={ext[0]:.2f} y={ext[1]:.2f} z={ext[2]:.2f} nm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

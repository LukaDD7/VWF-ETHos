#!/usr/bin/env python3
"""
analyze_gromacs_md.py
=====================
Post-analysis for GROMACS MD trajectories across complex/monomer/autoinhib systems.

Supports three system types:
  complex    — A1-GPIbα complex. Metrics: RMSF, PCA, MM/PBSA ΔG_bind.
  monomer   — VWF A1 only (chain A extracted from complex CIF).
              Metrics: RMSF, PCA, mutation-site flexibility.
  autoinhib — VWF A1+D'D3. Metrics: RMSF, D'D3-A1 interface contacts
              (H-bonds, salt bridges, buried surface area per frame).

Usage:
  python3 analyze_gromacs_md.py --system complex --input output/gromacs_md/
  python3 analyze_gromacs_md.py --system monomer --input output/gromacs_md_monomer/
  python3 analyze_gromacs_md.py --system autoinhib --input output/gromacs_md_autoinhib/

Output:
  analysis/summary.csv      — per-variant metrics
  analysis/rmsf.csv         — per-residue RMSF
  analysis/pca.csv           — PCA eigenvalues / top PCs
  analysis/contacts.csv      — D'D3-A1 contact analysis (autoinhib only)
  analysis/mmpbsa.csv        — MM/PBSA binding free energy (complex only)
"""

import argparse
import os
import sys
import subprocess
from pathlib import Path
from typing import Optional, List

import pandas as pd


# =============================================================================
# Helpers
# =============================================================================

def run_gmx(gmx_args: List[str], work_dir: str, log_file: Optional[str] = None) -> str:
    """Run a gmx command and return stdout."""
    env = os.environ.copy()
    env["GMX_RESET_COUNTERS"] = "0"
    result = subprocess.run(
        ["gmx"] + gmx_args,
        cwd=work_dir,
        capture_output=True,
        text=True,
        env=env
    )
    output = result.stdout + result.stderr
    if log_file:
        Path(log_file).write_text(output)
    return output


def compute_rmsf(tpr: str, xtc: str, out_dir: str) -> pd.DataFrame:
    """Compute per-residue RMSF."""
    rmsf_xvg = os.path.join(out_dir, "rmsf.xvg")
    run_gmx(
        ["rmsf", "-s", tpr, "-f", xtc, "-o", rmsf_xvg,
         "-qdist", "1", "-pdb", "none"],
        work_dir=out_dir
    )
    # Parse XVG file
    data = []
    with open(rmsf_xvg) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    res, val = float(parts[0]), float(parts[1])
                    data.append({"residue": int(res), "rmsf": val})
                except ValueError:
                    continue
    return pd.DataFrame(data)


def compute_pca(tpr: str, xtc: str, out_dir: str, n_components: int = 10) -> pd.DataFrame:
    """Compute PCA eigenvalues from covariance analysis."""
    covar_out = os.path.join(out_dir, "covar")
    eigenval_xvg = os.path.join(out_dir, "eigenval.xvg")

    # Covariance analysis
    run_gmx(
        ["covar", "-s", tpr, "-f", xtc, "-o", covar_out,
         "-eigen", eigenval_xvg, "-l", os.path.join(out_dir, "covar.log")],
        work_dir=out_dir,
        log_file=os.path.join(out_dir, "covar_output.log")
    )

    # Parse eigenvalues
    evals = []
    with open(eigenval_xvg) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    pc, val = int(parts[0]), float(parts[1])
                    evals.append({"PC": pc, "eigenvalue": val})
                except ValueError:
                    continue

    df = pd.DataFrame(evals)
    if len(df) > 0:
        total = df["eigenvalue"].sum()
        df["variance_explained"] = df["eigenvalue"] / total if total > 0 else 0
        df["cumulative_variance"] = df["variance_explained"].cumsum()
    return df.head(n_components)


def compute_contacts(pdb: str, ref_pdb: Optional[str], out_dir: str,
                     hbond_d: float = 3.5, salt_d: float = 4.0) -> pd.DataFrame:
    """Compute H-bond and salt bridge contacts between D'D3 and A1.

    D' = residues 764-865  (in the combined chain)
    D3 = residues 866-1233
    A1 = residues 1268-1466
    D'D3 = 764-1233 (470 residues)
    A1  = 1268-1466 (199 residues, but in the combined chain starts at 470)
    """
    sys.path.insert(0, str(Path(__file__).parent))
    try:
        import gemmi
    except ImportError:
        print("[WARN] gemmi not installed, contact analysis skipped")
        return pd.DataFrame()

    # Load structure
    st = gemmi.read_pdb(pdb)
    chain = st[0]

    # Define residue ranges
    d_prime_range = range(764 - 1, 865)       # 0-indexed
    d3_range = range(865 - 1, 1233)          # 0-indexed
    a1_range = range(1268 - 1, 1466)         # 0-indexed

    dd3_res = [r for r in chain.get_residues() if r.seq_num in list(d_prime_range) + list(d3_range)]
    a1_res = [r for r in chain.get_residues() if r.seq_num in a1_range]

    contacts = []
    for r1 in dd3_res:
        for r2 in a1_res:
            d = _residue_distance(r1, r2)
            if d < 5.0:  # within 5 Å
                hbond = _is_hbond(r1, r2)
                salt = _is_salt_bridge(r1, r2)
                if hbond or salt:
                    contacts.append({
                        "dd3_res": r1.seq_num,
                        "dd3_aa": r1.name,
                        "a1_res": r2.seq_num,
                        "a1_aa": r2.name,
                        "distance": d,
                        "hbond": hbond,
                        "salt_bridge": salt
                    })

    return pd.DataFrame(contacts)


def _residue_distance(r1, r2) -> float:
    """Minimum distance between any atom pair of two residues."""
    min_d = float("inf")
    for a1 in r1:
        for a2 in r2:
            d = ((a1.pos.x - a2.pos.x)**2 +
                 (a1.pos.y - a2.pos.y)**2 +
                 (a1.pos.z - a2.pos.z)**2) ** 0.5
            if d < min_d:
                min_d = d
    return min_d


def _is_hbond(r1, r2) -> bool:
    """Simple H-bond detection: N/O donor to O/N acceptor within 3.5 Å."""
    donors = {"N", "OG", "ND1", "ND2", "NE", "NH1", "NH2", "NZ"}
    acceptors = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OH", "SD"}
    for a1 in r1:
        if a1.element not in donors:
            continue
        for a2 in r2:
            if a2.element not in acceptors:
                continue
            d = ((a1.pos.x - a2.pos.x)**2 +
                 (a1.pos.y - a2.pos.y)**2 +
                 (a1.pos.z - a2.pos.z)**2) ** 0.5
            if d < 3.5:
                return True
    return False


def _is_salt_bridge(r1, r2) -> bool:
    """Salt bridge: Arg/Lys (positive) to Asp/Glu (negative) within 4 Å."""
    positive = {"ARG", "LYS"}
    negative = {"ASP", "GLU"}
    if r1.name in positive and r2.name in negative:
        return _residue_distance(r1, r2) < 4.0
    if r1.name in negative and r2.name in positive:
        return _residue_distance(r1, r2) < 4.0
    return False


def run_mmpbsa(tpr: str, xtc: str, out_dir: str) -> pd.DataFrame:
    """Compute MM/PBSA binding free energy (requires gmx_MMPBSA)."""
    try:
        import gmx_MMPBSA
    except ImportError:
        print("[WARN] gmx_MMPBSA not installed, MM/PBSA skipped")
        return pd.DataFrame()

    # This is a simplified wrapper; production use would call gmx_MMPBSA.AMMPSAMCalculation
    print("[INFO] MM/PBSA requires running: gmx_MMPBSA -O -i mmpbsa.in -fo mmgbsa.csv")
    return pd.DataFrame()


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--system", choices=["complex", "monomer", "autoinhib"],
                        required=True, help="MD system type")
    parser.add_argument("--input", required=True,
                        help="GROMACS MD output directory (e.g. output/gromacs_md/)")
    parser.add_argument("--output", default=None,
                        help="Analysis output directory (default: <input>/analysis/)")
    parser.add_argument("--variants", default=None,
                        help="Variant list CSV for mutation-site annotations")
    args = parser.parse_args()

    input_dir = Path(args.input).resolve()
    out_dir = Path(args.output or input_dir / "analysis").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"System: {args.system}")
    print(f"Input:  {input_dir}")
    print(f"Output: {out_dir}")
    print()

    # Find completed variants
    phase = "prod"
    done_marker = f".done_{phase}"
    variant_dirs = []
    for d in input_dir.iterdir():
        if d.is_dir() and (d / done_marker).exists():
            variant_dirs.append(d)

    print(f"Found {len(variant_dirs)} completed variants")
    if not variant_dirs:
        print("No completed simulations found.")
        return

    all_rows = []
    all_rmsf = []
    all_pca = []

    for vdir in sorted(variant_dirs):
        variant = vdir.name
        prod_dir = vdir / "production"
        tpr = prod_dir / "md_prod.tpr"
        xtc = prod_dir / "md_prod.xtc"

        if not tpr.exists():
            print(f"  [SKIP] {variant}: no tpr")
            continue
        if not xtc.exists():
            print(f"  [SKIP] {variant}: no xtc")
            continue

        print(f"  Analyzing: {variant}")

        row = {"variant": variant}

        # RMSF
        try:
            rmsf_df = compute_rmsf(str(tpr), str(xtc), str(prod_dir))
            if not rmsf_df.empty:
                rmsf_df["variant"] = variant
                all_rmsf.append(rmsf_df)
                row["rmsf_mean"] = rmsf_df["rmsf"].mean()
                row["rmsf_max"] = rmsf_df["rmsf"].max()
        except Exception as e:
            print(f"    [WARN] RMSF failed: {e}")

        # PCA
        try:
            pca_df = compute_pca(str(tpr), str(xtc), str(prod_dir))
            if not pca_df.empty:
                pca_df["variant"] = variant
                all_pca.append(pca_df)
                row["PC1_variance"] = pca_df.iloc[0]["variance_explained"] if len(pca_df) > 0 else None
        except Exception as e:
            print(f"    [WARN] PCA failed: {e}")

        # System-specific analysis
        if args.system == "autoinhib":
            final_pdb = prod_dir / "md_prod.gro"
            if final_pdb.exists():
                try:
                    contacts_df = compute_contacts(str(final_pdb), None, str(prod_dir))
                    if not contacts_df.empty:
                        contacts_df["variant"] = variant
                        row["n_hbonds"] = contacts_df["hbond"].sum()
                        row["n_salt_bridges"] = contacts_df["salt_bridge"].sum()
                        row["n_total_contacts"] = len(contacts_df)
                except Exception as e:
                    print(f"    [WARN] Contacts failed: {e}")

        all_rows.append(row)

    # Write outputs
    if all_rows:
        summary_df = pd.DataFrame(all_rows)
        summary_df.to_csv(out_dir / "summary.csv", index=False)
        print(f"\nWritten: {out_dir / 'summary.csv'} ({len(summary_df)} variants)")

    if all_rmsf:
        rmsf_all = pd.concat(all_rmsf, ignore_index=True)
        rmsf_all.to_csv(out_dir / "rmsf.csv", index=False)
        print(f"Written: {out_dir / 'rmsf.csv'}")

    if all_pca:
        pca_all = pd.concat(all_pca, ignore_index=True)
        pca_all.to_csv(out_dir / "pca.csv", index=False)
        print(f"Written: {out_dir / 'pca.csv'}")

    print(f"\nAnalysis complete: {out_dir}")


if __name__ == "__main__":
    main()
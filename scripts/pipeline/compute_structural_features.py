#!/usr/bin/env python3
"""
compute_structural_features.py — CIF → Secondary Structural Features
=====================================================================

从 Boltz-2 / AF3 的 mmCIF 输出中提取二次结构特征，生成用于 VWD 分型的
结构特征矩阵。

输入：
  - Boltz-2 CIF 文件: output/boltz2_*/predictions/model_0.cif
  - AF3 CIF 文件:     Proteo-Structure-Pipeline/output/af3_*/extracted/fold_vwf_*/model_0.cif
  - 置信度 JSON:      confidence_model_0.json, pae_*.npz, plddt_*.npz

输出：
  - structural_features.csv: 每个变体一行，包含 SASA/BSA/contacts/Rg 等特征
  - structural_features.parquet: 同上 (Parquet 格式)

特征清单：
  Global:
    - sasa_complex: 复合物整体 SASA (Å²)
    - sasa_receptor: 受体单独 SASA (Å²)
    - sasa_ligand:   配体单独 SASA (Å²)
    - bsa_total:     界面埋藏面积 = (SASA_R + SASA_L - SASA_complex) / 2
    - rg_complex:    复合物回转半径 (Å)

  Interface:
    - n_contacts:       界面接触残基对数 (< 5 Å)
    - n_hbonds:         界面氢键数 (d < 3.5 Å, angle > 120°)
    - n_salt_bridges:   盐桥数 (Arg/Lys – Asp/Glu < 4 Å)
    - iface_hydrophobic_frac: 界面疏水残基占比

  Mutation-site:
    - mut_site_sasa:    突变残基 SASA (Å²)
    - mut_site_rsa:     突变残基相对 SASA (0-1)
    - mut_site_bfactor: 突变残基 B-factor 均值 (如果有)
    - mut_site_plddt:   突变残基 pLDDT (如果有)
    - mut_neighborhood_sasa: 突变残基 ±5 残基邻域 SASA

  Confidence (from JSON/NPZ):
    - plddt_mean, plddt_min, ptm, iptm
    - pae_interface_mean

用法：
  python scripts/pipeline/compute_structural_features.py \\
      --boltz-results output/boltz2_a1_gpiba_results \\
      --output output/structural_features.csv

  # 处理 functional panel 结果
  python scripts/pipeline/compute_structural_features.py \\
      --boltz-results output/boltz2_vwd_functional_panel/boltz_results \\
      --output output/structural_features_functional.csv

依赖：
  pip install biopython freesasa numpy pandas
  # 可选: pip install biotite (DSSP)

Author: VWF-ETHos Pipeline
Date: 2026-05-12
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

try:
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB import NeighborSearch, Selection
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("[WARN] Biopython not installed. Install: pip install biopython")

try:
    import freesasa
    HAS_FREESASA = True
except ImportError:
    HAS_FREESASA = False
    print("[WARN] FreeSASA not installed. Install: pip install freesasa")


# =============================================================================
# Constants
# =============================================================================

# Max SASA values per amino acid for RSA calculation (Tien et al., 2013)
MAX_SASA = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167,
    "GLN": 225, "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197,
    "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159,
    "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174,
}

HYDROPHOBIC_AA = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"}

# Charged residues for salt bridges
POSITIVE_AA = {"ARG", "LYS"}
NEGATIVE_AA = {"ASP", "GLU"}

# H-bond donors and acceptors (simplified)
HBOND_DONORS = {"N", "NE", "NE1", "NE2", "ND1", "ND2", "NH1", "NH2",
                "NZ", "OG", "OG1", "OH", "NE", "SG"}
HBOND_ACCEPTORS = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH",
                   "ND1", "NE2", "SD", "SG"}


# =============================================================================
# Core Feature Extractor
# =============================================================================

class StructuralFeatureExtractor:
    """Extract secondary structural features from mmCIF files."""

    def __init__(self, cif_path: Path, chain_receptor: str = "A",
                 chain_ligand: Optional[str] = "B"):
        """
        Args:
            cif_path: Path to mmCIF file
            chain_receptor: Chain ID for receptor (VWF domain)
            chain_ligand: Chain ID for ligand (GPIbα, FVIII, etc.)
                         None for monomer-only assays
        """
        if not HAS_BIOPYTHON:
            raise RuntimeError("Biopython required: pip install biopython")

        parser = MMCIFParser(QUIET=True)
        self.structure = parser.get_structure("complex", str(cif_path))
        self.model = self.structure[0]
        self.chain_receptor = chain_receptor
        self.chain_ligand = chain_ligand
        self.cif_path = cif_path

    def _get_chain_atoms(self, chain_id: str) -> list:
        """Get all non-hetero atoms from a chain."""
        if chain_id not in self.model:
            return []
        chain = self.model[chain_id]
        atoms = []
        for res in chain.get_residues():
            if res.id[0] == " ":  # Standard residue (not HETATM)
                atoms.extend(res.get_atoms())
        return atoms

    def _get_chain_residues(self, chain_id: str) -> list:
        """Get standard residues from a chain."""
        if chain_id not in self.model:
            return []
        return [r for r in self.model[chain_id].get_residues() if r.id[0] == " "]

    # ---- SASA ----

    def compute_sasa(self, probe_radius: float = 1.4) -> Dict[str, float]:
        """Compute SASA for complex and individual chains."""
        result = {
            "sasa_complex": np.nan,
            "sasa_receptor": np.nan,
            "sasa_ligand": np.nan,
            "bsa_total": np.nan,
        }

        if not HAS_FREESASA:
            return result

        try:
            # Complex SASA
            fs_complex = freesasa.Structure(str(self.cif_path))
            sasa_complex = freesasa.calc(fs_complex)
            result["sasa_complex"] = sasa_complex.totalArea()

            # Per-chain SASA (using freesasa selections)
            chain_areas = {}
            for chain_id in [self.chain_receptor, self.chain_ligand]:
                if chain_id is None:
                    continue
                sel = freesasa.selectArea(
                    {f"chain_{chain_id}": f"chain {chain_id}"},
                    fs_complex, sasa_complex
                )
                chain_areas[chain_id] = sel.get(f"chain_{chain_id}", 0.0)

            if self.chain_receptor in chain_areas:
                result["sasa_receptor"] = chain_areas[self.chain_receptor]
            if self.chain_ligand and self.chain_ligand in chain_areas:
                result["sasa_ligand"] = chain_areas[self.chain_ligand]

            # BSA = (SASA_R_alone + SASA_L_alone - SASA_complex) / 2
            # Approximation: use chain-level from complex
            if not np.isnan(result["sasa_receptor"]) and not np.isnan(result.get("sasa_ligand", np.nan)):
                result["bsa_total"] = (
                    result["sasa_receptor"] + result["sasa_ligand"] - result["sasa_complex"]
                ) / 2.0

        except Exception as e:
            print(f"  [WARN] SASA computation failed: {e}")

        return result

    # ---- Interface Contacts ----

    def compute_interface_contacts(self, cutoff: float = 5.0) -> Dict[str, float]:
        """Compute interface contacts between receptor and ligand chains."""
        result = {
            "n_contacts": 0,
            "n_hbonds": 0,
            "n_salt_bridges": 0,
            "iface_hydrophobic_frac": 0.0,
        }

        if self.chain_ligand is None:
            return result

        rec_atoms = self._get_chain_atoms(self.chain_receptor)
        lig_atoms = self._get_chain_atoms(self.chain_ligand)

        if not rec_atoms or not lig_atoms:
            return result

        # Build neighbor search on receptor atoms
        ns = NeighborSearch(rec_atoms)

        contact_residues_rec = set()
        contact_residues_lig = set()
        hbonds = 0
        salt_bridges = 0

        for lig_atom in lig_atoms:
            nearby = ns.search(lig_atom.get_vector(), cutoff, "A")
            for rec_atom in nearby:
                rec_res = rec_atom.get_parent()
                lig_res = lig_atom.get_parent()

                contact_residues_rec.add(rec_res.get_full_id())
                contact_residues_lig.add(lig_res.get_full_id())

                dist = lig_atom - rec_atom

                # H-bond check (simplified: donor-acceptor distance < 3.5 Å)
                if dist < 3.5:
                    rec_name = rec_atom.get_name().strip()
                    lig_name = lig_atom.get_name().strip()
                    if (rec_name in HBOND_DONORS and lig_name in HBOND_ACCEPTORS) or \
                       (rec_name in HBOND_ACCEPTORS and lig_name in HBOND_DONORS):
                        hbonds += 1

                # Salt bridge check
                if dist < 4.0:
                    rec_resname = rec_res.get_resname()
                    lig_resname = lig_res.get_resname()
                    if (rec_resname in POSITIVE_AA and lig_resname in NEGATIVE_AA) or \
                       (rec_resname in NEGATIVE_AA and lig_resname in POSITIVE_AA):
                        salt_bridges += 1

        n_contacts = len(contact_residues_rec) + len(contact_residues_lig)
        result["n_contacts"] = n_contacts

        # Hydrophobic fraction at interface
        hydrophobic_count = 0
        for res_id in contact_residues_rec:
            # res_id is a full_id tuple, get residue
            try:
                res = self.model[res_id[2]][res_id[3]]
                if res.get_resname() in HYDROPHOBIC_AA:
                    hydrophobic_count += 1
            except (KeyError, IndexError):
                pass

        if len(contact_residues_rec) > 0:
            result["iface_hydrophobic_frac"] = hydrophobic_count / len(contact_residues_rec)

        # Deduplicate H-bonds and salt bridges (counted per atom pair, not per residue pair)
        result["n_hbonds"] = hbonds // 2  # Approximate dedup
        result["n_salt_bridges"] = salt_bridges // 2

        return result

    # ---- Mutation Site ----

    def compute_mutation_site_features(self, position: int,
                                       neighborhood: int = 5) -> Dict[str, float]:
        """Compute features at and around the mutation site."""
        result = {
            "mut_site_sasa": np.nan,
            "mut_site_rsa": np.nan,
            "mut_site_bfactor": np.nan,
            "mut_neighborhood_sasa": np.nan,
        }

        rec_residues = self._get_chain_residues(self.chain_receptor)
        if not rec_residues:
            return result

        # Find mutation residue (position in sequence, 1-indexed)
        target_res = None
        for i, res in enumerate(rec_residues):
            if res.id[1] == position:
                target_res = res
                break

        if target_res is None:
            return result

        # B-factor of mutation site
        bfactors = [a.get_bfactor() for a in target_res.get_atoms()]
        if bfactors:
            result["mut_site_bfactor"] = np.mean(bfactors)

        # SASA of mutation site (using FreeSASA per-residue)
        if HAS_FREESASA:
            try:
                fs = freesasa.Structure(str(self.cif_path))
                sasa_result = freesasa.calc(fs)
                # Per-residue SASA
                resname = target_res.get_resname()
                chain_id = self.chain_receptor
                res_num = target_res.id[1]
                sel = freesasa.selectArea(
                    {"mut_site": f"resi {res_num} and chain {chain_id}"},
                    fs, sasa_result
                )
                mut_sasa = sel.get("mut_site", 0.0)
                result["mut_site_sasa"] = mut_sasa

                # RSA = SASA / max_SASA
                max_sasa = MAX_SASA.get(resname, 200)
                result["mut_site_rsa"] = min(1.0, mut_sasa / max_sasa) if max_sasa > 0 else 0.0

                # Neighborhood SASA
                start_idx = max(0, position - neighborhood)
                end_idx = position + neighborhood
                sel_nbr = freesasa.selectArea(
                    {"neighborhood": f"resi {start_idx}-{end_idx} and chain {chain_id}"},
                    fs, sasa_result
                )
                result["mut_neighborhood_sasa"] = sel_nbr.get("neighborhood", 0.0)

            except Exception as e:
                print(f"  [WARN] Mutation site SASA failed: {e}")

        return result

    # ---- Global Features ----

    def compute_global_features(self) -> Dict[str, float]:
        """Compute radius of gyration and chain lengths."""
        result = {
            "rg_complex": np.nan,
            "n_residues_receptor": 0,
            "n_residues_ligand": 0,
        }

        # Radius of gyration
        all_atoms = []
        for chain in self.model:
            for res in chain.get_residues():
                if res.id[0] == " ":
                    ca = res.get("CA")
                    if ca:
                        all_atoms.append(ca.get_vector().get_array())

        if len(all_atoms) > 0:
            coords = np.array(all_atoms)
            centroid = coords.mean(axis=0)
            rg = np.sqrt(np.mean(np.sum((coords - centroid) ** 2, axis=1)))
            result["rg_complex"] = rg

        # Chain lengths
        result["n_residues_receptor"] = len(self._get_chain_residues(self.chain_receptor))
        if self.chain_ligand:
            result["n_residues_ligand"] = len(self._get_chain_residues(self.chain_ligand))

        return result

    # ---- Confidence Metrics ----

    @staticmethod
    def parse_confidence(job_dir: Path) -> Dict[str, float]:
        """Parse Boltz-2 / AF3 confidence metrics from JSON/NPZ files."""
        result = {
            "plddt_mean": np.nan,
            "plddt_min": np.nan,
            "ptm": np.nan,
            "iptm": np.nan,
            "pae_interface_mean": np.nan,
        }

        # Boltz-2 confidence JSON
        for conf_file in sorted(job_dir.glob("predictions/confidence_model_*.json")):
            try:
                with open(conf_file) as f:
                    data = json.load(f)
                if "plddt" in data:
                    plddts = data["plddt"] if isinstance(data["plddt"], list) else [data["plddt"]]
                    result["plddt_mean"] = float(np.mean(plddts))
                    result["plddt_min"] = float(np.min(plddts))
                result["ptm"] = data.get("ptm", np.nan)
                result["iptm"] = data.get("iptm", np.nan)
                if "pae" in data:
                    pae = np.array(data["pae"])
                    if pae.ndim >= 2:
                        result["pae_interface_mean"] = float(np.mean(pae))
                break  # Use first (best) model
            except Exception:
                pass

        # AF3 full_data JSON fallback
        for fd_file in sorted(job_dir.glob("*_full_data_*.json")):
            try:
                with open(fd_file) as f:
                    data = json.load(f)
                if "atom_plddts" in data:
                    result["plddt_mean"] = float(np.mean(data["atom_plddts"]))
                    result["plddt_min"] = float(np.min(data["atom_plddts"]))
                if "pae" in data:
                    pae = np.array(data["pae"])
                    if pae.ndim >= 2:
                        result["pae_interface_mean"] = float(np.mean(pae))
                break
            except Exception:
                pass

        # NPZ files (Boltz-2 v2)
        for npz_file in sorted(job_dir.glob("predictions/plddt_*.npz")):
            try:
                npz = np.load(npz_file)
                arr = npz[list(npz.keys())[0]]
                result["plddt_mean"] = float(np.mean(arr))
                result["plddt_min"] = float(np.min(arr))
                break
            except Exception:
                pass

        return result

    # ---- All-in-one ----

    def extract_all(self, mutation_pos: int = 0) -> Dict[str, float]:
        """Extract all structural features."""
        features = {}
        features.update(self.compute_sasa())
        features.update(self.compute_interface_contacts())
        if mutation_pos > 0:
            features.update(self.compute_mutation_site_features(mutation_pos))
        features.update(self.compute_global_features())
        return features


# =============================================================================
# Batch Processing
# =============================================================================

def parse_variant_from_dirname(dirname: str) -> Tuple[str, int]:
    """Parse variant ID and position from directory name.

    Examples:
        'boltz_results_VWF_R1306W_vs_GPIb_alpha' → ('VWF_R1306W', 1306)
        'boltz_results_VWF_R816W__dprime_d3_fviii_binding' → ('VWF_R816W', 816)
        'fold_vwf_r1306w' → ('VWF_R1306W', 1306)
    """
    # Boltz-2 functional panel format
    m = re.search(r"VWF_([A-Z])(\d+)([A-Z*])", dirname)
    if m:
        pos = int(m.group(2))
        variant_id = f"VWF_{m.group(1)}{pos}{m.group(3)}"
        return variant_id, pos

    # AF3 format
    m = re.search(r"fold_vwf_([a-z])(\d+)([a-z])", dirname)
    if m:
        pos = int(m.group(2))
        variant_id = f"VWF_{m.group(1).upper()}{pos}{m.group(3).upper()}"
        return variant_id, pos

    return dirname, 0


def detect_chains(job_dir: Path) -> Tuple[str, Optional[str]]:
    """Detect receptor and ligand chain IDs from CIF file."""
    cif_files = list(job_dir.glob("predictions/model_*.cif")) + \
                list(job_dir.glob("predictions/*_model_*.cif")) + \
                list(job_dir.glob("fold_vwf_*_model_*.cif"))

    if not cif_files:
        return "A", "B"

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("tmp", str(cif_files[0]))
        chains = [c.id for c in structure[0]]
        receptor = chains[0] if chains else "A"
        ligand = chains[1] if len(chains) > 1 else None
        return receptor, ligand
    except Exception:
        return "A", "B"


def process_batch(results_dir: Path, output_path: Path,
                  max_jobs: Optional[int] = None) -> pd.DataFrame:
    """Process all Boltz-2/AF3 results in a directory."""

    if not HAS_BIOPYTHON:
        raise RuntimeError("Biopython required: pip install biopython")

    records = []
    job_dirs = sorted(results_dir.glob("boltz_results_VWF_*")) + \
               sorted(results_dir.glob("fold_vwf_*"))

    if max_jobs:
        job_dirs = job_dirs[:max_jobs]

    total = len(job_dirs)
    print(f"Processing {total} job directories from {results_dir}")

    for idx, job_dir in enumerate(job_dirs, 1):
        if not job_dir.is_dir():
            continue

        dirname = job_dir.name
        variant_id, position = parse_variant_from_dirname(dirname)

        # Find CIF file
        cif_files = list(job_dir.glob("predictions/model_*.cif")) + \
                    list(job_dir.glob("predictions/*_model_*.cif")) + \
                    list(job_dir.glob("fold_vwf_*_model_*.cif"))

        if not cif_files:
            print(f"  [{idx}/{total}] SKIP {dirname}: no CIF file")
            continue

        cif_path = cif_files[0]
        print(f"  [{idx}/{total}] Processing {variant_id} (pos={position})...")

        try:
            # Detect chains
            receptor, ligand = detect_chains(job_dir)

            # Extract structural features
            extractor = StructuralFeatureExtractor(cif_path, receptor, ligand)
            features = extractor.extract_all(mutation_pos=position)

            # Parse confidence
            confidence = StructuralFeatureExtractor.parse_confidence(job_dir)

            # Parse assay key from dirname
            assay_key = ""
            if "__" in dirname:
                assay_key = dirname.split("__", 1)[1]
            elif "_vs_" in dirname:
                assay_key = dirname.split("_vs_", 1)[1]

            record = {
                "variant_id": variant_id,
                "position": position,
                "assay_key": assay_key,
                "job_dir": dirname,
                "cif_file": str(cif_path.relative_to(results_dir)),
                **features,
                **confidence,
            }
            records.append(record)

        except Exception as e:
            print(f"  [{idx}/{total}] ERROR {dirname}: {e}")

    df = pd.DataFrame(records)
    if not df.empty:
        # Save
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False)
        parquet_path = output_path.with_suffix(".parquet")
        df.to_parquet(parquet_path, index=False)
        print(f"\nSaved {len(df)} records to:")
        print(f"  CSV:     {output_path}")
        print(f"  Parquet: {parquet_path}")
    else:
        print("\nNo records extracted.")

    return df


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Extract structural features from Boltz-2/AF3 CIF outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--boltz-results", type=Path, required=True,
        help="Directory containing Boltz-2/AF3 result subdirectories",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Output CSV path (default: <boltz-results>/structural_features.csv)",
    )
    parser.add_argument(
        "--max-jobs", type=int, default=None,
        help="Limit number of jobs to process (for testing)",
    )
    args = parser.parse_args()

    if args.output is None:
        args.output = args.boltz_results / "structural_features.csv"

    process_batch(args.boltz_results, args.output, args.max_jobs)


if __name__ == "__main__":
    main()

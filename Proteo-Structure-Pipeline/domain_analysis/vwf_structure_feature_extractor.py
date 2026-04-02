#!/usr/bin/env python3
"""
VWF Structure Feature Extractor
Extracts domain-specific structural features from AlphaFold3 predictions

This module extracts features relevant to Type-2 classification:
1. Domain-specific pLDDT profiles
2. Local structural stability around functional sites
3. Mutation-induced structural perturbations
4. Electrostatic and hydrophobic environment changes

Author: Claude Code
Date: 2026-04-02
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field
import json

# BioPython for structure parsing
try:
    from Bio.PDB import MMCIFParser, Superimposer, PDBIO, Select
    from Bio.PDB.PDBParser import PDBParser
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("Warning: BioPython not available. Install with: pip install biopython")


# Import domain definitions from classifier
import sys
sys.path.insert(0, str(Path(__file__).parent))
from vwf_type2_literature_based_classifier import (
    VWF_DOMAIN_ARCHITECTURE,
    FUNCTIONAL_REGIONS,
    KNOWN_TYPE2_HOTSPOTS
)


@dataclass
class DomainFeatures:
    """Structural features for a specific VWF domain."""
    domain_name: str
    start_residue: int
    end_residue: int

    # pLDDT statistics
    plddt_mean: float = 0.0
    plddt_std: float = 0.0
    plddt_min: float = 0.0
    plddt_max: float = 0.0
    low_confidence_regions: List[Tuple[int, int]] = field(default_factory=list)

    # Structural stability
    b_factor_mean: float = 0.0
    local_compactness: float = 0.0  # Based on residue packing

    # Functional site proximity
    contains_functional_site: bool = False
    functional_site_names: List[str] = field(default_factory=list)
    distance_to_functional_site: float = 0.0


@dataclass
class MutationImpactFeatures:
    """Features describing the structural impact of a mutation."""
    variant_id: str
    position: int
    ref_aa: str
    alt_aa: str

    # Local structural changes
    local_rmsd: float = 0.0
    neighborhood_rmsd: float = 0.0  # 10A radius

    # pLDDT changes
    pldtt_delta_site: float = 0.0
    plddt_delta_local: float = 0.0  # Averaged over 5-residue window

    # Amino acid property changes
    hydrophobicity_change: float = 0.0
    charge_change: int = 0
    size_change: float = 0.0

    # Domain-specific features
    domain_name: str = ""
    is_functional_site: bool = False
    is_hotspot: bool = False

    # Structural environment
    solvent_accessibility_change: float = 0.0
    secondary_structure_change: str = ""
    hydrogen_bond_changes: int = 0


class VWFFeatureExtractor:
    """Extracts VWF-specific structural features from AlphaFold3 predictions."""

    # Amino acid properties
    AA_HYDROPHOBICITY = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
        'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
        'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
        'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }

    AA_CHARGE = {
        'D': -1, 'E': -1, 'K': 1, 'R': 1, 'H': 0.5  # Histidine partially charged
    }

    AA_SIZE = {
        'A': 88.6, 'C': 108.5, 'D': 111.1, 'E': 138.4, 'F': 189.9,
        'G': 60.1, 'H': 153.2, 'I': 166.7, 'K': 168.6, 'L': 166.7,
        'M': 162.9, 'N': 114.1, 'P': 112.7, 'Q': 143.8, 'R': 173.4,
        'S': 89.0, 'T': 116.1, 'V': 140.0, 'W': 227.8, 'Y': 193.6
    }

    def __init__(self, wt_structure_path: Optional[str] = None):
        self.domains = VWF_DOMAIN_ARCHITECTURE
        self.functional_regions = FUNCTIONAL_REGIONS
        self.hotspots = KNOWN_TYPE2_HOTSPOTS

        self.wt_structure = None
        self.wt_plddt = None

        if wt_structure_path and BIOPYTHON_AVAILABLE:
            self.load_wt_structure(wt_structure_path)

    def load_wt_structure(self, structure_path: str):
        """Load wild-type VWF structure."""
        if not BIOPYTHON_AVAILABLE:
            return

        parser = MMCIFParser(QUIET=True)
        self.wt_structure = parser.get_structure("VWF_WT", structure_path)
        self.wt_plddt = self._extract_plddt_from_structure(self.wt_structure)

    def _extract_plddt_from_structure(self, structure) -> np.ndarray:
        """Extract pLDDT values from structure B-factors."""
        plddt_values = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    # pLDDT is stored in B-factor column
                    ca_atom = residue['CA'] if 'CA' in residue else None
                    if ca_atom:
                        plddt_values.append(ca_atom.get_bfactor())

        return np.array(plddt_values)

    def extract_domain_features(self, structure_path: str) -> Dict[str, DomainFeatures]:
        """Extract features for each VWF domain from a structure."""
        if not BIOPYTHON_AVAILABLE:
            return {}

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("VWF", structure_path)
        plddt = self._extract_plddt_from_structure(structure)

        domain_features = {}

        for domain_name, (start, end) in self.domains.items():
            # Adjust for 0-based indexing
            start_idx = max(0, start - 1)
            end_idx = min(len(plddt), end)

            domain_plddt = plddt[start_idx:end_idx]

            # Find low confidence regions (< 70 pLDDT)
            low_conf_regions = []
            in_low_conf = False
            region_start = 0

            for i, p in enumerate(domain_plddt):
                if p < 70 and not in_low_conf:
                    region_start = start + i
                    in_low_conf = True
                elif p >= 70 and in_low_conf:
                    low_conf_regions.append((region_start, start + i - 1))
                    in_low_conf = False

            if in_low_conf:
                low_conf_regions.append((region_start, end))

            # Check for functional sites
            contains_functional = False
            functional_names = []
            min_distance = float('inf')

            for region_name, region_info in self.functional_regions.items():
                if "residues" in region_info:
                    for res_range in region_info["residues"]:
                        if isinstance(res_range, tuple) and len(res_range) == 2:
                            # Check overlap
                            if not (res_range[1] < start or res_range[0] > end):
                                contains_functional = True
                                functional_names.append(region_name)
                                # Calculate distance
                                overlap_start = max(start, res_range[0])
                                overlap_end = min(end, res_range[1])
                                distance = 0 if overlap_start <= overlap_end else \
                                    min(abs(start - res_range[1]), abs(end - res_range[0]))
                                min_distance = min(min_distance, distance)

            features = DomainFeatures(
                domain_name=domain_name,
                start_residue=start,
                end_residue=end,
                plddt_mean=float(np.mean(domain_plddt)),
                plddt_std=float(np.std(domain_plddt)),
                plddt_min=float(np.min(domain_plddt)),
                plddt_max=float(np.max(domain_plddt)),
                low_confidence_regions=low_conf_regions,
                contains_functional_site=contains_functional,
                functional_site_names=functional_names,
                distance_to_functional_site=min_distance if contains_functional else -1
            )

            domain_features[domain_name] = features

        return domain_features

    def calculate_mutation_impact(
        self,
        variant_id: str,
        position: int,
        ref_aa: str,
        alt_aa: str,
        mut_structure_path: str
    ) -> MutationImpactFeatures:
        """Calculate structural impact features for a mutation."""

        features = MutationImpactFeatures(
            variant_id=variant_id,
            position=position,
            ref_aa=ref_aa,
            alt_aa=alt_aa
        )

        # Determine domain
        for domain_name, (start, end) in self.domains.items():
            if start <= position <= end:
                features.domain_name = domain_name
                break

        # Check if in functional site
        for region_name, region_info in self.functional_regions.items():
            if "residues" in region_info:
                for res_range in region_info["residues"]:
                    if isinstance(res_range, tuple) and len(res_range) == 2:
                        if res_range[0] <= position <= res_range[1]:
                            features.is_functional_site = True
                            break

        # Check if in hotspot
        for subtype, hotspot_info in self.hotspots.items():
            if isinstance(hotspot_info, dict) and "residues" in hotspot_info:
                if position in hotspot_info["residues"]:
                    features.is_hotspot = True
                    break
            elif isinstance(hotspot_info, dict):
                # Handle nested structure for Type 2M
                for subcategory, subinfo in hotspot_info.items():
                    if isinstance(subinfo, dict) and "residues" in subinfo:
                        if position in subinfo["residues"]:
                            features.is_hotspot = True
                            break

        # Calculate amino acid property changes
        if ref_aa in self.AA_HYDROPHOBICITY and alt_aa in self.AA_HYDROPHOBICITY:
            features.hydrophobicity_change = (
                self.AA_HYDROPHOBICITY[alt_aa] - self.AA_HYDROPHOBICITY[ref_aa]
            )

        ref_charge = self.AA_CHARGE.get(ref_aa, 0)
        alt_charge = self.AA_CHARGE.get(alt_aa, 0)
        features.charge_change = alt_charge - ref_charge

        if ref_aa in self.AA_SIZE and alt_aa in self.AA_SIZE:
            features.size_change = self.AA_SIZE[alt_aa] - self.AA_SIZE[ref_aa]

        # Calculate structural changes if structures are available
        if BIOPYTHON_AVAILABLE and self.wt_structure:
            features = self._calculate_structural_changes(
                features, mut_structure_path
            )

        return features

    def _calculate_structural_changes(
        self,
        features: MutationImpactFeatures,
        mut_structure_path: str
    ) -> MutationImpactFeatures:
        """Calculate RMSD and pLDDT changes between WT and mutant."""
        if not BIOPYTHON_AVAILABLE or not self.wt_structure:
            return features

        try:
            parser = MMCIFParser(QUIET=True)
            mut_structure = parser.get_structure("VWF_MUT", mut_structure_path)

            # Extract pLDDT
            mut_plddt = self._extract_plddt_from_structure(mut_structure)

            # Calculate pLDDT changes
            pos_idx = features.position - 1
            if pos_idx < len(self.wt_plddt) and pos_idx < len(mut_plddt):
                features.pldtt_delta_site = mut_plddt[pos_idx] - self.wt_plddt[pos_idx]

                # Local pLDDT change (5-residue window)
                window = 2
                start = max(0, pos_idx - window)
                end = min(len(self.wt_plddt), pos_idx + window + 1)
                wt_local = np.mean(self.wt_plddt[start:end])
                mut_local = np.mean(mut_plddt[start:end])
                features.plddt_delta_local = mut_local - wt_local

            # Calculate RMSD (simplified - full structure alignment)
            wt_atoms = []
            mut_atoms = []

            for wt_res, mut_res in zip(
                self.wt_structure[0].get_residues(),
                mut_structure[0].get_residues()
            ):
                if 'CA' in wt_res and 'CA' in mut_res:
                    wt_atoms.append(wt_res['CA'])
                    mut_atoms.append(mut_res['CA'])

            if len(wt_atoms) == len(mut_atoms) and len(wt_atoms) > 0:
                superimposer = Superimposer()
                superimposer.set_atoms(wt_atoms, mut_atoms)
                features.local_rmsd = superimposer.rms

        except Exception as e:
            print(f"Error calculating structural changes: {e}")

        return features

    def extract_all_features(
        self,
        variants_df: pd.DataFrame,
        structures_dir: str
    ) -> pd.DataFrame:
        """
        Extract features for all variants in a dataframe.

        Expected columns in variants_df:
        - variant_id
        - position
        - ref_aa
        - alt_aa
        """
        structures_dir = Path(structures_dir)
        results = []

        for _, row in variants_df.iterrows():
            variant_id = row.get("variant_id", "")
            position = int(row.get("position", 0))
            ref_aa = row.get("ref_aa", "")
            alt_aa = row.get("alt_aa", "")

            # Find mutant structure
            mut_structure = structures_dir / f"fold_vwf_{variant_id.lower()}" / f"{variant_id}.cif"

            if mut_structure.exists():
                features = self.calculate_mutation_impact(
                    variant_id=variant_id,
                    position=position,
                    ref_aa=ref_aa,
                    alt_aa=alt_aa,
                    mut_structure_path=str(mut_structure)
                )

                results.append({
                    "variant_id": variant_id,
                    "position": position,
                    "domain": features.domain_name,
                    "is_functional_site": features.is_functional_site,
                    "is_hotspot": features.is_hotspot,
                    "plddt_delta_site": features.pldtt_delta_site,
                    "plddt_delta_local": features.plddt_delta_local,
                    "local_rmsd": features.local_rmsd,
                    "hydrophobicity_change": features.hydrophobicity_change,
                    "charge_change": features.charge_change,
                    "size_change": features.size_change,
                })

        return pd.DataFrame(results)


def main():
    """Example usage."""
    print("VWF Structure Feature Extractor")
    print("=" * 80)

    extractor = VWFFeatureExtractor()

    # Show domain architecture
    print("\nVWF Domain Architecture (from Literature):")
    print("-" * 80)
    for domain, (start, end) in extractor.domains.items():
        length = end - start + 1
        print(f"{domain:20s}: {start:4d} - {end:4d} ({length:4d} aa)")

    # Show functional regions
    print("\n\nFunctional Regions (Type-2 Associated):")
    print("-" * 80)
    for region, info in extractor.functional_regions.items():
        assoc_type = info.get("associated_type", "N/A")
        mechanism = info.get("mechanism", "")[:50]
        assoc_type_str = assoc_type if assoc_type else "N/A"
        print(f"{region:25s}: Type {assoc_type_str:3s} - {mechanism}...")


if __name__ == "__main__":
    main()

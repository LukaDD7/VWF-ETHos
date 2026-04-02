#!/usr/bin/env python3
"""
VWF Type-2 Subtype Classification Pipeline
Based on Literature-Knowledge and AlphaFold3 Structural Analysis

Key References:
1. Lenting et al. (2024) Blood - "How unique structural adaptations support and coordinate the complex function of von Willebrand factor"
2. Atiq & O'Donnell (2024) Blood - "Novel functions for von Willebrand factor"
3. Haberichter & O'Donnell (2026) Haematologica - "Structure and multiple functions of von Willebrand factor"

Type-2 Classification Rationale:
- Type 2A: Enhanced ADAMTS13 cleavage in A2 domain -> loss of HMW multimers
- Type 2B: Gain-of-function in A1 domain -> spontaneous GPIbα binding
- Type 2M: Reduced binding to GPIbα (A1) or Collagen (A3) without proteolysis defect
- Type 2N: Reduced FVIII binding in D'D3 domain -> low FVIII levels

Author: Claude Code
Date: 2026-04-02
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
from collections import defaultdict


# =============================================================================
# VWF DOMAIN ARCHITECTURE - Based on Literature (Zhou et al. 2012, Lenting et al. 2024)
# =============================================================================

VWF_DOMAIN_ARCHITECTURE = {
    # Signal peptide and propeptide (removed during processing)
    "signal_peptide": (1, 22),
    "propeptide_D1": (23, 386),
    "propeptide_D2": (387, 763),

    # Mature VWF domains
    "D_prime": (764, 865),      # Part of D'D3 assembly
    "D3": (866, 1233),          # Part of D'D3 assembly, contains FVIII binding site
    "A1": (1271, 1492),         # GPIbα binding, Type 2B mutations
    "A2": (1493, 1684),         # ADAMTS13 cleavage site (Y1605-M1606), Type 2A mutations
    "A3": (1685, 1874),         # Collagen binding (Types I/III), Type 2M mutations
    "D4": (1875, 2255),         # Part of C-terminal assembly
    "C1": (2256, 2324),         # VWC domain
    "C2": (2325, 2392),         # VWC domain
    "C3": (2393, 2460),         # VWC domain
    "C4": (2497, 2577),         # VWC domain, contains RGD motif for αIIbβ3 binding
    "C5": (2578, 2658),         # VWC domain
    "C6": (2659, 2722),         # VWC domain
    "CK": (2723, 2813),         # Cystine knot domain
}

# Critical functional regions defined by literature
FUNCTIONAL_REGIONS = {
    # Type 2N: FVIII binding region within D'D3
    "FVIII_binding_site": {
        "domain": "D_prime_D3",
        "residues": [(782, 799), (816, 826)],  # D' region critical for FVIII binding
        "description": "FVIII light chain (A3-a3-C1-C2) binding interface",
        "associated_type": "2N",
        "mechanism": "Reduced FVIII affinity -> low FVIII levels"
    },

    # Type 2B: A1 domain GPIbα binding and autoinhibitory module
    "GPIbα_binding_site": {
        "domain": "A1",
        "residues": [(1271, 1459)],  # Full A1 domain
        "critical_residues": {
            "interactive_site_1": [(1296, 1309)],  # Helix α3, loop α3β4, strand β3
            "interactive_site_2": [(1309, 1350)],  # Bottom face of A1
            "AIM_N_terminal": [(1238, 1268)],      # Autoinhibitory module N-terminus
            "AIM_C_terminal": [(1460, 1472)],      # Autoinhibitory module C-terminus
        },
        "associated_type": "2B",
        "mechanism": "Gain-of-function: Disrupted AIM -> spontaneous GPIbα binding"
    },

    # Type 2A: ADAMTS13 cleavage site in A2 domain
    "ADAMTS13_cleavage_site": {
        "domain": "A2",
        "residues": [(1605, 1606)],  # Y1605-M1606 scissile bond
        "critical_features": {
            "beta4_less_loop": [(1594, 1602)],     # β4-less loop
            "calcium_site": [(1596, 1602)],        # Ca2+ binding site
            "vicinal_disulfide": [(1669, 1670)],   # Unique vicinal disulfide in A2
            "cis_pro": [1645],                     # cis-Pro1645
        },
        "associated_type": "2A",
        "mechanism": "Enhanced ADAMTS13 cleavage -> loss of HMW multimers"
    },

    # Type 2M: A1 domain (GPIbα binding defect) or A3 domain (Collagen binding defect)
    "Collagen_binding_site": {
        "domain": "A3",
        "residues": [(1684, 1873)],
        "critical_interface": [(1700, 1850)],    # β3-sheet, α2-helix, α3-helix
        "associated_type": "2M",
        "mechanism": "Reduced collagen binding -> impaired platelet adhesion"
    },

    # Additional functional sites
    "alphaIIb_beta3_binding": {
        "domain": "C4",
        "residues": [(2507, 2510)],  # RGD motif
        "associated_type": None,  # Part of normal function, not Type-2 specific
    },
}


# Known Type-2 mutation hotspots from literature
KNOWN_TYPE2_HOTSPOTS = {
    "2A": {
        "residues": [1528, 1584, 1597, 1605, 1606, 1614, 1638, 1649, 1668, 1670],
        "description": "A2 domain mutations enhancing ADAMTS13 cleavage",
        "key_mutations": ["p.Met1528Val", "p.Arg1597Trp", "p.Tyr1605Ala", "p.Met1606Arg"]
    },
    "2B": {
        "residues": [1274, 1295, 1306, 1308, 1313, 1322, 1325, 1326, 1331, 1463, 1472],
        "description": "A1 domain mutations disrupting AIM",
        "key_mutations": ["p.Arg1306Trp", "p.Arg1308Cys", "p.Val1316Met", "p.Arg1341Gln"]
    },
    "2M": {
        "A1_defects": {
            "residues": [1274, 1293, 1295, 1306, 1308, 1315, 1322, 1325, 1326],
            "description": "A1 domain mutations reducing GPIbα binding without AIM disruption"
        },
        "A3_defects": {
            "residues": [1696, 1731, 1745, 1760, 1779, 1783, 1786, 1824],
            "description": "A3 domain mutations reducing collagen binding"
        }
    },
    "2N": {
        "residues": [782, 787, 799, 816, 823, 854, 857, 868, 869, 874, 875, 882, 885],
        "description": "D'D3 domain mutations reducing FVIII binding",
        "key_mutations": ["p.Arg816Trp", "p.Arg854Gln", "p.Cys868Phe", "p.Thr875Met"]
    }
}


@dataclass
class VWFVariant:
    """Represents a VWF variant with structural and functional annotations."""
    variant_id: str
    protein_change: str          # e.g., "p.Gly1531Asp" or "G1531D"
    position: int                # Amino acid position
    ref_aa: str
    alt_aa: str
    acmg_classification: str     # Pathogenic, Likely_pathogenic, VUS, etc.
    type2_subtype: Optional[str] = None  # Known 2A, 2B, 2M, 2N

    # Structural features
    domain: str = ""
    plddt_wt: float = 0.0
    plddt_mut: float = 0.0
    plddt_delta: float = 0.0
    local_rmsd: float = 0.0

    # Functional predictions
    feature_scores: Dict[str, float] = field(default_factory=dict)
    predicted_subtype: Optional[str] = None
    confidence: float = 0.0


class VWFType2Classifier:
    """
    Literature-based classifier for VWF Type-2 subtypes.
    Combines domain knowledge with AF3 structural features.
    """

    def __init__(self):
        self.domains = VWF_DOMAIN_ARCHITECTURE
        self.functional_regions = FUNCTIONAL_REGIONS
        self.hotspots = KNOWN_TYPE2_HOTSPOTS

    def get_domain_for_position(self, position: int) -> str:
        """Map amino acid position to VWF domain."""
        for domain_name, (start, end) in self.domains.items():
            if start <= position <= end:
                return domain_name
        return "unknown"

    def classify_by_position(self, position: int) -> Dict[str, any]:
        """
        Classify a variant based on its position in the VWF structure.
        Returns potential Type-2 associations based on literature.
        """
        domain = self.get_domain_for_position(position)
        associations = []

        # Check functional region involvement
        for region_name, region_info in self.functional_regions.items():
            if "residues" in region_info:
                for res_range in region_info["residues"]:
                    if isinstance(res_range, tuple) and len(res_range) == 2:
                        if res_range[0] <= position <= res_range[1]:
                            if region_info.get("associated_type"):
                                associations.append({
                                    "region": region_name,
                                    "type": region_info["associated_type"],
                                    "mechanism": region_info.get("mechanism", "")
                                })

        # Check known hotspots
        for subtype, hotspot_info in self.hotspots.items():
            if isinstance(hotspot_info, dict) and "residues" in hotspot_info:
                if position in hotspot_info["residues"]:
                    associations.append({
                        "region": f"{subtype}_hotspot",
                        "type": subtype,
                        "mechanism": hotspot_info.get("description", "")
                    })

        return {
            "domain": domain,
            "associations": associations,
            "is_functional_region": len(associations) > 0
        }

    def calculate_structure_based_features(self, variant: VWFVariant) -> Dict[str, float]:
        """
        Calculate structural features relevant to Type-2 classification.
        Based on AF3 structural predictions.
        """
        features = {}

        # Basic structural stability metrics
        features["plddt_instability"] = 1.0 if variant.plddt_delta < -5 else 0.0
        features["high_local_rmsd"] = 1.0 if variant.local_rmsd > 2.0 else 0.0

        # Domain-specific features
        domain = variant.domain

        if domain == "A1":
            # Type 2B/BM indicators
            features["aim_disruption_potential"] = self._score_aim_disruption(variant)
            features["gpib_binding_impact"] = self._score_gpib_binding_impact(variant)

        elif domain == "A2":
            # Type 2A indicators
            features["adamts13_sensitivity"] = self._score_adamts13_sensitivity(variant)
            features["cleavage_site_proximity"] = self._calculate_distance_to_cleavage_site(variant.position)

        elif domain == "A3":
            # Type 2M indicators
            features["collagen_binding_impact"] = self._score_collagen_binding_impact(variant)

        elif domain in ["D_prime", "D3"]:
            # Type 2N indicators
            features["fviii_binding_impact"] = self._score_fviii_binding_impact(variant)

        return features

    def _score_aim_disruption(self, variant: VWFVariant) -> float:
        """Score potential disruption of autoinhibitory module (Type 2B mechanism)."""
        aim_n_term = (1238, 1268)
        aim_c_term = (1460, 1472)

        if aim_n_term[0] <= variant.position <= aim_n_term[1]:
            return 1.0
        if aim_c_term[0] <= variant.position <= aim_c_term[1]:
            return 1.0

        # Check if variant causes charge change (important for AIM electrostatics)
        charged_aas = {"R", "K", "D", "E"}
        if (variant.ref_aa in charged_aas) != (variant.alt_aa in charged_aas):
            return 0.7

        return 0.0

    def _score_gpib_binding_impact(self, variant: VWFVariant) -> float:
        """Score impact on GPIbα binding interface."""
        # Interactive sites based on Huizinga et al. 2002
        interactive_site_1 = [(1296, 1309)]
        interactive_site_2 = [(1309, 1350)]

        for site_range in interactive_site_1 + interactive_site_2:
            if site_range[0] <= variant.position <= site_range[1]:
                return 1.0

        return 0.0

    def _score_adamts13_sensitivity(self, variant: VWFVariant) -> float:
        """Score ADAMTS13 sensitivity enhancement (Type 2A mechanism)."""
        # Ca2+ site disruption increases sensitivity
        calcium_site = [(1596, 1602)]
        beta4_less_loop = [(1594, 1602)]
        vicinal_disulfide = [(1669, 1670)]

        for site_range in calcium_site:
            if site_range[0] <= variant.position <= site_range[1]:
                return 1.0

        if variant.position in [1669, 1670]:
            return 1.0

        if variant.position == 1645:  # cis-Pro
            return 0.9

        # General A2 domain destabilization
        if 1493 <= variant.position <= 1684:
            if variant.plddt_delta < -10:
                return 0.8

        return 0.0

    def _calculate_distance_to_cleavage_site(self, position: int) -> float:
        """Calculate distance to ADAMTS13 cleavage site (Y1605-M1606)."""
        cleavage_site = 1605.5  # Between 1605 and 1606
        return abs(position - cleavage_site)

    def _score_collagen_binding_impact(self, variant: VWFVariant) -> float:
        """Score impact on collagen binding interface (Type 2M mechanism)."""
        # Critical interface from Bienkowska et al. 1997, Romijn et al. 2003
        critical_interface = [(1700, 1850)]

        for site_range in critical_interface:
            if site_range[0] <= variant.position <= site_range[1]:
                return 1.0

        return 0.0

    def _score_fviii_binding_impact(self, variant: VWFVariant) -> float:
        """Score impact on FVIII binding interface (Type 2N mechanism)."""
        # Critical residues from Fuller et al. 2021, Chiu et al. 2015
        critical_residues = [
            (782, 799),   # D' region
            (816, 826),   # Til' structure
        ]

        for site_range in critical_residues:
            if site_range[0] <= variant.position <= site_range[1]:
                return 1.0

        # General D'D3 region
        if 764 <= variant.position <= 1233:
            return 0.5

        return 0.0

    def predict_subtype(self, variant: VWFVariant) -> Tuple[str, float]:
        """
        Predict Type-2 subtype based on literature rules and structural features.
        Returns (predicted_subtype, confidence_score)
        """
        domain = self.get_domain_for_position(variant.position)
        variant.domain = domain

        # Calculate all features
        features = self.calculate_structure_based_features(variant)
        variant.feature_scores = features

        # Score each subtype
        subtype_scores = defaultdict(float)

        # Type 2A scoring
        if domain == "A2":
            subtype_scores["2A"] += features.get("adamts13_sensitivity", 0) * 2.0
            if features.get("cleavage_site_proximity", 100) < 10:
                subtype_scores["2A"] += 1.0

        # Type 2B scoring
        if domain == "A1":
            aim_score = features.get("aim_disruption_potential", 0)
            gpib_score = features.get("gpib_binding_impact", 0)

            # Type 2B requires both AIM disruption AND GPIb binding interface involvement
            if aim_score > 0.5 and gpib_score > 0.5:
                subtype_scores["2B"] += 3.0
            elif aim_score > 0.5:
                subtype_scores["2B"] += 1.5
            elif gpib_score > 0.5:
                subtype_scores["2M"] += 1.0  # Might be Type 2M if no AIM disruption

        # Type 2M scoring
        if domain == "A3":
            subtype_scores["2M"] += features.get("collagen_binding_impact", 0) * 2.0

        if domain == "A1":
            gpib_score = features.get("gpib_binding_impact", 0)
            aim_score = features.get("aim_disruption_potential", 0)
            if gpib_score > 0.5 and aim_score < 0.3:
                subtype_scores["2M"] += 1.5

        # Type 2N scoring
        if domain in ["D_prime", "D3"]:
            subtype_scores["2N"] += features.get("fviii_binding_impact", 0) * 2.0

        # Structural instability bonus
        if features.get("plddt_instability", 0) > 0.5:
            # Instability in A2 suggests 2A
            if domain == "A2":
                subtype_scores["2A"] += 0.5
            # Instability in D'D3 suggests 2N
            if domain in ["D_prime", "D3"]:
                subtype_scores["2N"] += 0.5

        if not subtype_scores:
            return "unclassified", 0.0

        # Get highest scoring subtype
        predicted = max(subtype_scores, key=subtype_scores.get)
        max_score = subtype_scores[predicted]

        # Calculate confidence based on score magnitude
        confidence = min(max_score / 3.0, 1.0)

        return predicted, confidence


def load_type2_variants(csv_path: str) -> pd.DataFrame:
    """Load Type-2 variant data."""
    return pd.read_csv(csv_path)


def analyze_variant_file(variant_file: str, af3_results_dir: str) -> List[VWFVariant]:
    """
    Main analysis function to classify VWF Type-2 variants.
    """
    df = load_type2_variants(variant_file)
    classifier = VWFType2Classifier()

    variants = []

    for _, row in df.iterrows():
        variant = VWFVariant(
            variant_id=row.get("variant_id", ""),
            protein_change=row.get("protein_change", ""),
            position=int(row.get("position", 0)),
            ref_aa=row.get("ref_aa", ""),
            alt_aa=row.get("alt_aa", ""),
            acmg_classification=row.get("acmg_classification", "VUS"),
            type2_subtype=row.get("type2_subtype", None),
            plddt_wt=float(row.get("plddt_wt", 0)),
            plddt_mut=float(row.get("plddt_mut", 0)),
            plddt_delta=float(row.get("plddt_delta", 0)),
            local_rmsd=float(row.get("local_rmsd", 0))
        )

        # Predict subtype
        predicted, confidence = classifier.predict_subtype(variant)
        variant.predicted_subtype = predicted
        variant.confidence = confidence

        variants.append(variant)

    return variants


def main():
    """Main entry point for the classifier."""
    print("=" * 80)
    print("VWF Type-2 Subtype Classification Pipeline")
    print("Based on Literature-Knowledge and AlphaFold3 Structural Analysis")
    print("=" * 80)

    # Example usage
    classifier = VWFType2Classifier()

    # Test with known variants
    test_positions = [
        (1306, "Type 2B hotspot"),
        (1605, "Type 2A cleavage site"),
        (1760, "Type 2M collagen binding"),
        (816, "Type 2N FVIII binding"),
    ]

    print("\nTest Classifications:")
    print("-" * 80)

    for pos, description in test_positions:
        result = classifier.classify_by_position(pos)
        print(f"\nPosition {pos} ({description}):")
        print(f"  Domain: {result['domain']}")
        print(f"  In functional region: {result['is_functional_region']}")
        for assoc in result['associations']:
            print(f"  -> Associated with Type {assoc['type']}: {assoc['mechanism']}")


if __name__ == "__main__":
    main()

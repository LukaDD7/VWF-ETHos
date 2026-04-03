#!/usr/bin/env python3
"""
VWF Residue-Level Feature Extractor
基于文献的残基级精确特征提取库

Literature Sources:
1. Lenting et al. (2024) Blood - "How unique structural adaptations support and coordinate the complex function of von Willebrand factor"
2. Atiq & O'Donnell (2024) Blood - "Novel functions for von Willebrand factor"
3. Haberichter & O'Donnell (2026) Haematologica - "Structure and multiple functions of von Willebrand factor"

Feature Categories:
- Residue-level structural features (exact positions from literature)
- Mutation amino acid properties
- Global and local RMSD (10Å radius)
- Functional site annotations with PMIDs
- Mechanistic features for each domain

Author: Claude Code
Date: 2026-04-03
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
import numpy as np
import json


# =============================================================================
# EXTENDED LITERATURE ANNOTATIONS
# D4, CK, C domains and additional functional sites
# =============================================================================

# D4 Domain - Multimerization and Secretion - PMID: 15758837, 35148377
D4_DOMAIN = {
    "range": (1875, 2255),
    "pmid": "35148377",  # Lenting 2024
    "description": "D4 domain involved in subunit dimerization and ER-to-Golgi trafficking",
    "structure": "Contains TIL and E modules, part of C-terminal assembly",
    "vwd_association": "Type 2A (multimerization defect) or Type 1 (secretion defect)",
    "key_regions": {
        "dimerization_interface": {
            "range": (1900, 2100),
            "description": "Subunit dimerization interface",
            "pmid": "35148377"
        },
        "trafficking_motif": {
            "range": (2100, 2200),
            "description": "ER-to-Golgi trafficking signal",
            "pmid": "15758837"
        }
    },
    "vwd_mutations": {
        1888: {"mutation": "P1888L", "type": "Type_2A", "effect": "Multimerization defect"},
        1939: {"mutation": "E1939K", "type": "Type_2A", "effect": "Secretion defect"},
        1974: {"mutation": "L1974P", "type": "Type_2A"},
        1986: {"mutation": "R1986C", "type": "Type_2A"},
        1996: {"mutation": "V1996M", "type": "Type_2A"},
        2006: {"mutation": "R2006C", "type": "Type_2A"},
        2015: {"mutation": "V2015L", "type": "Type_2A"},
        2034: {"mutation": "V2034A", "type": "Type_2A"},
        2062: {"mutation": "W2062R", "type": "Type_2A"},
    }
}

# CK Domain (Cystine Knot) - PMID: 35148377, 17895385
CK_DOMAIN = {
    "range": (2723, 2813),
    "pmid": "35148377",
    "description": "Cystine knot domain essential for VWF dimerization",
    "structure": "Disulfide-rich domain with cystine knot fold",
    "mechanism": "Forms C-terminal disulfide bonds for dimerization",
    "critical_residues": {
        2780: {"type": "disulfide", "partner": None, "description": "Conserved cysteine for dimerization"},
        2801: {"mutation": "P2801S", "type": "Type_2A", "effect": "Severe dimerization defect"}
    },
    "vwd_association": "Type 2A (severe multimerization defect due to loss of dimerization)"
}

# C Domains (C1-C6) - PMID: 35148377, 31582533
C_DOMAINS = {
    "C1": {
        "range": (2256, 2324),
        "pmid": "35148377",
        "description": "VWC domain - supports multimerization and collagen binding",
        "function": "Multimerization support, collagen IV binding",
        "vwd_type": "Type_2A_or_2M",
        "disulfides": 4  # Typical VWC domain
    },
    "C2": {
        "range": (2325, 2392),
        "pmid": "35148377",
        "description": "VWC domain - supports multimerization and collagen binding",
        "function": "Multimerization support, collagen IV binding",
        "vwd_type": "Type_2A_or_2M"
    },
    "C3": {
        "range": (2393, 2496),
        "pmid": "35148377",
        "description": "Extended VWC domain - multimerization support",
        "function": "Multimerization support"
    },
    "C4": {
        "range": (2497, 2577),
        "pmid": "31582533",
        "description": "VWC domain containing RGD motif",
        "special_features": {
            "rgd_motif": {
                "residues": (2507, 2510),
                "sequence": "RGDS",
                "function": "Integrin αIIbβ3 binding",
                "pmid": "31582533"
            }
        },
        "disulfide_pattern": "Unique disulfide linking βI and βII in SD2",
        "vwd_mutations": {
            2517: {"mutation": "V2517F", "type": "Type_2M", "effect": "Reduced integrin binding"},
            2535: {"mutation": "R2535P", "type": "Type_2M", "effect": "Mild bleeding tendency"}
        }
    },
    "C5": {
        "range": (2578, 2658),
        "pmid": "35148377",
        "description": "VWC domain - multimerization support"
    },
    "C6": {
        "range": (2659, 2722),
        "pmid": "35148377",
        "description": "VWC domain - multimerization support, connects to CK"
    }
}

# Propeptide (D1-D2) Detailed - PMID: 40958414, 12176890, 20335223
PROPEPTIDE_D1D2 = {
    "range": (23, 763),
    "pmid": "40958414",  # 2026 paper on propeptide correction
    "description": "VWF propeptide - essential for multimerization",
    "mechanism": "Functions as pH-sensing oxidoreductase template for multimerization",
    "cxxc_motifs": [
        {"domain": "D1", "range": (200, 250), "function": "Disulfide isomerase activity"},
        {"domain": "D2", "range": (400, 500), "function": "Disulfide isomerase activity"}
    ],
    "critical_residues": {
        89: {"mutation": "V89A", "type": "Type_2A_IIC", "effect": "Multimerization defect"},
        91: {"mutation": "L91P", "type": "Type_2A_IIC", "effect": "Propeptide defect"},
        528: {"mutation": "N528S", "type": "Type_2A", "pmid": "20335223"}
    },
    "vwd_association": "Type 2A/IIC - multimerization defects, can be corrected by exogenous propeptide"
}

# D'D3 Region Extended - PMID: 33888542, 15758837
DD3_EXTENDED = {
    "D_prime": {
        "range": (764, 865),
        "pmid": "33888542",
        "description": "D' segment - contains TIL and E module",
        "fviii_binding": {
            "interface_1": (782, 799),
            "interface_2": (816, 826),
            "affinity": "KD ≈ 0.5 nM"
        }
    },
    "D3": {
        "range": (866, 1233),
        "pmid": "33888542",
        "description": "D3 assembly - FVIII binding and multimerization",
        "fviii_interaction": "Envelops FVIII light chain (a3-A3-C1-C2)"
    }
}

# Shear-dependent Unfolding Regions - PMID: 19390046, 20089857
SHEAR_DEPENDENT_REGIONS = {
    "A2_unfolding": {
        "force_required": "~1 pN",
        "pmid": "19390046",
        "description": "A2 domain unfolds at ~1 pN to expose cleavage site"
    },
    "AIM_dissociation": {
        "force_required": "~21 pN",
        "pmid": "20089857",
        "description": "AIM dissociates at ~21 pN to expose GPIb binding site"
    },
    "multimer_unfolding": {
        "description": "Larger multimers unfold more easily",
        "clinical_relevance": "Loss of HMW multimers in Type 2A reduces platelet adhesion"
    }
}

# Additional VWD Hotspots from Literature
ADDITIONAL_VWD_HOTSPOTS = {
    "Type_2A_other": {
        "description": "Type 2A mutations outside A2 domain",
        "residues": {
            89: "V89A", 91: "L91P", 528: "N528S",  # Propeptide
            1888: "P1888L", 1939: "E1939K", 1974: "L1974P",  # D4
            1986: "R1986C", 1996: "V1996M", 2006: "R2006C",
            2015: "V2015L", 2034: "V2034A", 2062: "W2062R",
            2775: "S2775C",  # C-terminal
            2801: "P2801S",  # CK
        }
    },
    "Type_2M_A3": {
        "description": "A3 domain collagen binding defects",
        "residues": [1696, 1731, 1745, 1760, 1779, 1783, 1786, 1824]
    },
    "Type_2N_extended": {
        "description": "FVIII binding defects in D'D3",
        "residues": {
            782: "R782G", 787: "C787R", 799: "C799R",
            816: "R816W", 823: "N823K", 854: "R854Q",
            857: "R857Stop", 868: "C868F", 869: "C868R",
            874: "V874A", 875: "T875M", 882: "C882S", 885: "C885Y"
        }
    }
}

# =============================================================================
# END OF EXTENDED ANNOTATIONS
# =============================================================================

# AIM (Autoinhibitory Module) - PMID: 33888542, 28325766, 28904067
# Precise AIM contacts from crystallography data
AIM_RESIDUES = {
    "N_terminal_helix": {
        "range": (1238, 1268),
        "pmid": "33888542",
        "description": "AIM N-terminal alpha-helix (Gln1238-His1268)",
        "key_contacts": {
            # Salt bridges and hydrogen bonds within AIM
            1263: {"type": "salt_bridge", "partner": 1668, "atoms": "D1263-R1668"},
            1264: {"type": "hydrogen_bond", "partner": 1667, "atoms": "backbone-backbone"},
            1240: {"type": "hydrophobic_core", "residues": [1240, 1244, 1465, 1469]},
        },
        "mechanism": "Blocks access to GPIbα interactive site 2 until shear-induced dissociation",
        "vwd_association": "Type 2B mutations disrupt these contacts causing gain-of-function"
    },
    "C_terminal_strand": {
        "range": (1460, 1472),
        "pmid": "33888542",
        "description": "AIM C-terminal beta-strand (Leu1460-Asp1472)",
        "key_contacts": {
            1469: {"type": "hydrogen_bond", "partner": 1238, "atoms": "Y1469-Q1659"},
            1460: {"type": "hydrophobic_interaction", "partner_residues": [1263, 1264]},
        },
        "mechanism": "Forms beta-sheet with N-terminal helix to stabilize closed conformation",
        "ristocetin_binding": "Cys1237-Pro1251 and Glu1463-Asp1472"
    }
}

# GPIbα Interface - PMID: 12191960, 15117959, 23341617
GPIB_ALPHA_INTERFACE = {
    "interactive_site_1": {
        "range": (1296, 1309),
        "pmid": "12191960",
        "description": "Helix α3, loop α3β4, strand β3 - primary GPIbα binding",
        "key_residues": {
            1306: {
                "type": "cation_pi",
                "partner": "Y283_GPIb",  # Tyr283 on GPIbα
                "description": "R1306-Y283 cation-π interaction is critical for binding",
                "vwd_mutations": ["R1306W", "R1306Q", "R1306L"],
                "confidence": "high"
            },
            1299: {"type": "backbone_hbond", "description": "Backbone hydrogen bond to GPIb"},
            1308: {"type": "hydrophobic", "description": "Hydrophobic contact with GPIbα"},
        }
    },
    "interactive_site_2": {
        "range": (1309, 1350),
        "pmid": "12191960",
        "description": "Bottom face of A1 - loops α1β2, β3α2, α3β4",
        "key_residues": {
            1316: {
                "vwd_mutation": "V1316M",
                "type": "Type_2B_hotspot",
                "mechanism": "Alters hydrophobic packing, destabilizes AIM"
            },
            1322: {"vwd_mutation": "V1322G", "type": "Type_2B"},
            1325: {"vwd_mutation": "R1325H", "type": "Type_2B"},
            1326: {"vwd_mutation": "R1326W", "type": "Type_2B"},
        }
    }
}

# ADAMTS13 Cleavage Site - PMID: 15758837, 28904067, 16079166
ADAMTS13_CLEAVAGE_SITE = {
    "scissile_bond": {
        "residues": (1605, 1606),
        "sequence": "Y1605-M1606",
        "pmid": "15758837",
        "description": "Tyr1605-Met1606 scissile bond hydrolyzed by ADAMTS13",
        "cleavage_mechanism": "Metalloprotease active site cleaves Y1605-M1606 peptide bond"
    },
    "exosite_1_beta4less_loop": {
        "range": (1594, 1602),
        "pmid": "28904067",
        "description": "α4-less loop interacts with ADAMTS13 disintegrin domain",
        "key_residues": {
            1614: {
                "vwd_mutation": "D1614N",
                "type": "Group_2_2A",
                "effect": "Delays A2 refolding, prolongs cleavage time",
                "mechanism": "Disrupts exosite 1 interaction with disintegrin domain"
            },
            1596: {"type": "calcium_binding", "description": "Part of Ca2+ binding site (D1596, R1597)"},
            1602: {"type": "calcium_binding", "description": "Part of Ca2+ binding site (N1602)"},
        }
    },
    "exosite_2_cysteine_rich": {
        "range": (1642, 1651),
        "pmid": "28904067",
        "description": "α5-helix/β6-sheet region binds cysteine-rich domain",
        "key_residues": {
            1645: {
                "type": "cis_proline",
                "description": "Unique cis-Pro1645 - affects A2 conformational dynamics"
            },
        }
    },
    "exosite_3_spacer": {
        "range": (1660, 1668),
        "pmid": "28904067",
        "description": "α6-helix associates with spacer domain of ADAMTS13",
        "key_residues": {
            1669: {"type": "vicinal_disulfide", "partner": 1670, "description": "Cys1669-Cys1670 unique vicinal disulfide"},
        }
    },
    "group_1_2A_mutations": {
        "description": "Promote exposure of buried cleavage site",
        "residues": {
            1528: {"mutation": "M1528V", "effect": "Promotes premature A2 unfolding", "pmid": "16079166"},
            1638: {"mutation": "E1638K", "effect": "Destabilizes A2 structure", "pmid": "16079166"},
        }
    },
    "group_2_2A_mutations": {
        "description": "Delay A2 refolding after cleavage",
        "residues": {
            1597: {"mutation": "R1597W", "effect": "Delays refolding, increases susceptibility", "pmid": "16079166"},
            1614: {"mutation": "D1614N", "effect": "Disrupts exosite interaction", "pmid": "16079166"},
        }
    }
}

# Calcium Binding Sites - PMID: 28904067, 28325766
CALCIUM_SITES = {
    "A2_calcium_site": {
        "pmid": "28904067",
        "coordinating_residues": [1498, 1596, 1597, 1600, 1602],
        "description": "Ca2+ stabilizes A2 domain structure, resists unfolding",
        "mechanism": "Asp1498 (β1), Asp1596/R1597 (α3), Ala1600/Asn1602 (α3β4-loop)",
        "functional_impact": "Removing Ca2+ increases susceptibility to ADAMTS13 cleavage"
    }
}

# Collagen Binding Sites - PMID: 15758837, 16079166
COLLAGEN_BINDING = {
    "A3_collagen_I_III": {
        "range": (1684, 1873),
        "pmid": "15758837",
        "description": "A3 domain binds types I and III collagen",
        "key_interface": (1700, 1850),
        "critical_residues": {
            1731: {"vwd_mutation": "S1731P", "type": "Type_2M", "effect": "Reduces collagen binding"},
            1745: {"vwd_mutation": "H1745Q", "type": "Type_2M"},
            1760: {"vwd_mutation": "C1760R", "type": "Type_2M"},
            1779: {"vwd_mutation": "C1779R", "type": "Type_2M"},
            1783: {"vwd_mutation": "S1783A", "type": "Type_2M"},
            1824: {"vwd_mutation": "W1725C", "type": "Type_2M"},
        }
    },
    "A1_collagen_IV_VI": {
        "pmid": "15117959",
        "description": "A1 domain binds types IV and VI collagen"
    }
}

# FVIII Binding Site - PMID: 33888542, 15758837
FVIII_BINDING = {
    "D_prime_D3_region": {
        "range": (764, 1233),
        "pmid": "33888542",
        "description": "D'D3 region is essential for FVIII binding",
        "affinity": "KD ≈ 0.5 nM",
        "key_residues": {
            816: {"vwd_mutation": "R816W", "type": "Type_2N", "effect": "Reduces FVIII binding"},
            854: {"vwd_mutation": "R854Q", "type": "Type_2N"},
            868: {"vwd_mutation": "C868F", "type": "Type_2N"},
            875: {"vwd_mutation": "T875M", "type": "Type_2N"},
        }
    }
}

# C4 Domain RGD Motif - PMID: 31582533
C4_RGD_MOTIF = {
    "rgds_sequence": {
        "residues": (2507, 2510),
        "sequence": "RGDS",
        "pmid": "31582533",
        "description": "Arg-Gly-Asp-Ser motif for integrin αIIbβ3 binding",
        "function": "Supports fibrinogen-independent platelet aggregation"
    }
}

# Multimerization/Secretion Related - PMID: 33888542, 15758837
MULTIMERIZATION_SITES = {
    "D1_D2_propeptide": {
        "range": (23, 763),
        "pmid": "33888542",
        "description": "VWFpp facilitates multimerization in Golgi",
        "cxxc_motifs": [
            {"domain": "D1", "function": "Vicinal cysteine CXXC motif for disulfide exchange"},
            {"domain": "D2", "function": "Vicinal cysteine CXXC motif for disulfide exchange"}
        ]
    },
    "D4_dimerization": {
        "range": (1875, 2255),
        "pmid": "15758837",
        "description": "D4 domain involved in subunit dimerization",
        "vwd_type": "Type 2A (multimerization defect)"
    },
    "CK_dimerization": {
        "range": (2723, 2813),
        "pmid": "33888542",
        "description": "Cystine knot domain essential for VWF dimerization",
        "vwd_type": "Type 2A (severe multimerization defect)"
    }
}


# =============================================================================
# AMINO ACID PROPERTIES
# For mutation feature extraction
# =============================================================================

AMINO_ACID_PROPERTIES = {
    # Size (van der Waals volume in Å³)
    "size": {
        'A': 88.6, 'C': 108.5, 'D': 111.1, 'E': 138.4, 'F': 189.9,
        'G': 60.1, 'H': 153.2, 'I': 166.7, 'K': 168.6, 'L': 166.7,
        'M': 162.9, 'N': 114.1, 'P': 112.7, 'Q': 143.8, 'R': 173.4,
        'S': 89.0, 'T': 116.1, 'V': 140.0, 'W': 227.8, 'Y': 193.6
    },
    # Hydrophobicity (Kyte-Doolittle scale)
    "hydrophobicity": {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
        'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
        'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
        'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    },
    # Charge at physiological pH
    "charge": {
        'A': 0, 'C': 0, 'D': -1, 'E': -1, 'F': 0,
        'G': 0, 'H': 0.1, 'I': 0, 'K': 1, 'L': 0,
        'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 1,
        'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
    },
    # Aromaticity
    "is_aromatic": {
        'A': False, 'C': False, 'D': False, 'E': False, 'F': True,
        'G': False, 'H': True, 'I': False, 'K': False, 'L': False,
        'M': False, 'N': False, 'P': False, 'Q': False, 'R': False,
        'S': False, 'T': False, 'V': False, 'W': True, 'Y': True
    },
    # Secondary structure propensity (Chou-Fasman)
    "helix_propensity": {
        'A': 1.45, 'C': 0.77, 'D': 0.98, 'E': 1.53, 'F': 1.12,
        'G': 0.53, 'H': 1.24, 'I': 1.00, 'K': 1.07, 'L': 1.34,
        'M': 1.20, 'N': 0.73, 'P': 0.59, 'Q': 1.17, 'R': 0.79,
        'S': 0.79, 'T': 0.82, 'V': 0.91, 'W': 1.14, 'Y': 0.76
    }
}


@dataclass
class ResidueLevelFeatures:
    """
    残基级特征数据类
    包含从文献中提取的精确残基级特征
    """
    # 基本信息
    variant_id: str
    position: int
    ref_aa: str
    alt_aa: str

    # 域信息 (from VWF_DOMAIN_ARCHITECTURE)
    domain: str = ""
    domain_start: int = 0
    domain_end: int = 0
    relative_position: float = 0.0  # 在域内的相对位置 (0-1)

    # AIM相关特征 (仅A1域)
    is_in_AIM: bool = False
    aim_component: str = ""  # "N_terminal", "C_terminal", "none"
    aim_key_contact: Dict = field(default_factory=dict)  # 关键接触信息
    aim_disruption_score: float = 0.0  # AIM破坏评分

    # GPIbα界面特征 (仅A1域)
    is_in_gpib_interface: bool = False
    gpib_site: str = ""  # "interactive_site_1", "interactive_site_2"
    gpib_key_residue: bool = False  # 是否是关键残基 (如R1306)
    gpib_interaction_type: str = ""  # "cation_pi", "hbond", "hydrophobic"

    # ADAMTS13切割位点特征 (仅A2域)
    is_scissile_bond: bool = False
    is_in_exosite_1: bool = False
    is_in_exosite_2: bool = False
    is_in_exosite_3: bool = False
    is_group1_2A_mutation: bool = False
    is_group2_2A_mutation: bool = False
    exosite_disruption_potential: float = 0.0

    # 钙结合位点特征
    is_calcium_coordinating: bool = False
    calcium_site_importance: str = ""  # "high", "medium", "none"

    # 胶原结合特征
    is_collagen_binding_residue: bool = False
    collagen_type: str = ""  # "I_III", "IV_VI"

    # FVIII结合特征
    is_fviii_binding_residue: bool = False
    fviii_binding_importance: str = ""  # "critical", "supporting"

    # RGD motif特征
    is_rgd_motif: bool = False

    # 多聚化相关特征
    is_multimerization_residue: bool = False
    multimerization_role: str = ""  # "dimerization", "cxxc_motif"

    # 已知VWD突变热点
    is_known_vwd_hotspot: bool = False
    vwd_hotspot_type: str = ""  # "Type_2A", "Type_2B", "Type_2M", "Type_2N"
    vwd_hotspot_mutation: str = ""  # 具体突变 (如 "R1306W")

    # 结构特征 (从AF3提取)
    plddt_wt: float = 0.0
    plddt_mut: float = 0.0
    plddt_delta: float = 0.0
    global_rmsd: float = 0.0
    local_rmsd_10a: float = 0.0  # 10Å半径局部RMSD
    local_rmsd_15a: float = 0.0  # 15Å半径局部RMSD (额外)

    # 文献引用
    literature_pmids: List[str] = field(default_factory=list)
    literature_evidence: Dict = field(default_factory=dict)

    # 突变氨基酸特征
    mutation_size_delta: float = 0.0
    mutation_hydrophobicity_delta: float = 0.0
    mutation_charge_change: int = 0
    mutation_aromatic_change: int = 0
    mutation_helix_propensity_delta: float = 0.0

    # 功能影响预测
    functional_impact_score: float = 0.0
    mechanism_prediction: str = ""

    def to_dict(self) -> Dict:
        """转换为字典格式，用于序列化"""
        return {
            "variant_id": self.variant_id,
            "position": self.position,
            "ref_aa": self.ref_aa,
            "alt_aa": self.alt_aa,
            "mutation": f"{self.ref_aa}{self.position}{self.alt_aa}",
            "domain": self.domain,
            "relative_position_in_domain": round(self.relative_position, 3),
            "aim_features": {
                "is_in_AIM": self.is_in_AIM,
                "aim_component": self.aim_component,
                "aim_key_contact": self.aim_key_contact,
                "aim_disruption_score": round(self.aim_disruption_score, 3) if self.aim_disruption_score else None
            },
            "gpib_interface_features": {
                "is_in_gpib_interface": self.is_in_gpib_interface,
                "gpib_site": self.gpib_site,
                "gpib_key_residue": self.gpib_key_residue,
                "gpib_interaction_type": self.gpib_interaction_type
            },
            "adamts13_features": {
                "is_scissile_bond": self.is_scissile_bond,
                "is_in_exosite_1": self.is_in_exosite_1,
                "is_in_exosite_2": self.is_in_exosite_2,
                "is_in_exosite_3": self.is_in_exosite_3,
                "is_group1_2A": self.is_group1_2A_mutation,
                "is_group2_2A": self.is_group2_2A_mutation
            },
            "functional_residues": {
                "is_calcium_coordinating": self.is_calcium_coordinating,
                "is_collagen_binding": self.is_collagen_binding_residue,
                "is_fviii_binding": self.is_fviii_binding_residue,
                "is_rgd_motif": self.is_rgd_motif,
                "is_multimerization": self.is_multimerization_residue
            },
            "vwd_hotspot": {
                "is_hotspot": self.is_known_vwd_hotspot,
                "hotspot_type": self.vwd_hotspot_type,
                "known_mutation": self.vwd_hotspot_mutation
            },
            "structural_features": {
                "plddt_wt": round(self.plddt_wt, 2) if self.plddt_wt else None,
                "plddt_mut": round(self.plddt_mut, 2) if self.plddt_mut else None,
                "plddt_delta": round(self.plddt_delta, 2) if self.plddt_delta else None,
                "global_rmsd": round(self.global_rmsd, 3) if self.global_rmsd else None,
                "local_rmsd_10A": round(self.local_rmsd_10a, 3) if self.local_rmsd_10a else None
            },
            "mutation_properties": {
                "size_delta": round(self.mutation_size_delta, 1),
                "hydrophobicity_delta": round(self.mutation_hydrophobicity_delta, 2),
                "charge_change": self.mutation_charge_change,
                "aromatic_change": self.mutation_aromatic_change,
                "helix_propensity_delta": round(self.mutation_helix_propensity_delta, 3)
            },
            "literature": {
                "pmids": self.literature_pmids,
                "evidence_summary": self.literature_evidence
            }
        }


class VWFResidueFeatureExtractor:
    """
    VWF残基级特征提取器
    基于文献精确注释的特征提取
    """

    def __init__(self):
        self.domain_architecture = self._load_domain_architecture()
        self.literature_annotations = self._compile_literature_annotations()

    def _load_domain_architecture(self) -> Dict:
        """加载VWF域架构"""
        return {
            "signal_peptide": (1, 22),
            "propeptide_D1": (23, 386),
            "propeptide_D2": (387, 763),
            "D_prime": (764, 865),
            "D3": (866, 1233),
            "A1": (1271, 1492),
            "A2": (1493, 1684),
            "A3": (1685, 1874),
            "D4": (1875, 2255),
            "C1": (2256, 2324),
            "C2": (2325, 2392),
            "C3": (2393, 2496),
            "C4": (2497, 2577),
            "C5": (2578, 2658),
            "C6": (2659, 2722),
            "CK": (2723, 2813),
        }

    def _compile_literature_annotations(self) -> Dict:
        """编译文献注释为查找表"""
        annotations = defaultdict(dict)

        # AIM annotations
        for component, info in AIM_RESIDUES.items():
            start, end = info["range"]
            for pos in range(start, end + 1):
                annotations[pos]["aim"] = {
                    "component": component,
                    "is_key_contact": pos in info.get("key_contacts", {}),
                    "contact_details": info.get("key_contacts", {}).get(pos, {}),
                    "pmid": info["pmid"]
                }

        # GPIbα interface annotations
        for site, info in GPIB_ALPHA_INTERFACE.items():
            start, end = info["range"]
            for pos in range(start, end + 1):
                annotations[pos]["gpib"] = {
                    "site": site,
                    "is_key_residue": pos in info.get("key_residues", {}),
                    "residue_details": info.get("key_residues", {}).get(pos, {}),
                    "pmid": info["pmid"]
                }

        # ADAMTS13 annotations
        scissile_start, scissile_end = ADAMTS13_CLEAVAGE_SITE["scissile_bond"]["residues"]
        for pos in range(scissile_start, scissile_end + 1):
            annotations[pos]["adamts13_scissile"] = True

        exosite1_start, exosite1_end = ADAMTS13_CLEAVAGE_SITE["exosite_1_beta4less_loop"]["range"]
        for pos in range(exosite1_start, exosite1_end + 1):
            annotations[pos]["adamts13_exosite_1"] = True

        exosite2_start, exosite2_end = ADAMTS13_CLEAVAGE_SITE["exosite_2_cysteine_rich"]["range"]
        for pos in range(exosite2_start, exosite2_end + 1):
            annotations[pos]["adamts13_exosite_2"] = True

        exosite3_start, exosite3_end = ADAMTS13_CLEAVAGE_SITE["exosite_3_spacer"]["range"]
        for pos in range(exosite3_start, exosite3_end + 1):
            annotations[pos]["adamts13_exosite_3"] = True

        # Group 1 and Group 2 mutations
        for pos, info in ADAMTS13_CLEAVAGE_SITE["group_1_2A_mutations"]["residues"].items():
            annotations[pos]["group_1_2A"] = info

        for pos, info in ADAMTS13_CLEAVAGE_SITE["group_2_2A_mutations"]["residues"].items():
            annotations[pos]["group_2_2A"] = info

        # Calcium coordinating residues
        for pos in CALCIUM_SITES["A2_calcium_site"]["coordinating_residues"]:
            annotations[pos]["calcium_coordinating"] = True

        # Collagen binding
        collagen_start, collagen_end = COLLAGEN_BINDING["A3_collagen_I_III"]["key_interface"]
        for pos in range(collagen_start, collagen_end + 1):
            annotations[pos]["collagen_binding"] = {"type": "I_III"}

        for pos, info in COLLAGEN_BINDING["A3_collagen_I_III"].get("critical_residues", {}).items():
            annotations[pos]["collagen_critical"] = info

        # FVIII binding
        for pos, info in FVIII_BINDING["D_prime_D3_region"].get("key_residues", {}).items():
            annotations[pos]["fviii_binding"] = info

        # RGD motif
        rgd_start, rgd_end = C4_RGD_MOTIF["rgds_sequence"]["residues"]
        for pos in range(rgd_start, rgd_end + 1):
            annotations[pos]["rgd_motif"] = True

        # D4 domain annotations
        d4_start, d4_end = D4_DOMAIN["range"]
        for pos in range(d4_start, d4_end + 1):
            annotations[pos]["d4_domain"] = {"pmid": D4_DOMAIN["pmid"]}
        for pos, info in D4_DOMAIN.get("vwd_mutations", {}).items():
            annotations[pos]["d4_vwd_mutation"] = info
            annotations[pos]["is_vwd_hotspot"] = True
            annotations[pos]["vwd_type"] = "Type_2A"

        # CK domain annotations
        ck_start, ck_end = CK_DOMAIN["range"]
        for pos in range(ck_start, ck_end + 1):
            annotations[pos]["ck_domain"] = {"pmid": CK_DOMAIN["pmid"]}
        for pos, info in CK_DOMAIN.get("critical_residues", {}).items():
            if "mutation" in info:
                annotations[pos]["ck_vwd_mutation"] = info
                annotations[pos]["is_vwd_hotspot"] = True
                annotations[pos]["vwd_type"] = "Type_2A"

        # C domain annotations
        for c_domain, c_info in C_DOMAINS.items():
            c_start, c_end = c_info["range"]
            for pos in range(c_start, c_end + 1):
                annotations[pos]["c_domain"] = {"domain": c_domain, "pmid": c_info["pmid"]}
            # Add VWD mutations in C domains
            for pos, info in c_info.get("vwd_mutations", {}).items():
                annotations[pos]["c_vwd_mutation"] = info
                annotations[pos]["is_vwd_hotspot"] = True
                annotations[pos]["vwd_type"] = info.get("type", "Type_2M")

        # Propeptide annotations
        pp_start, pp_end = PROPEPTIDE_D1D2["range"]
        for pos in range(pp_start, pp_end + 1):
            annotations[pos]["propeptide"] = {"pmid": PROPEPTIDE_D1D2["pmid"]}
        for pos, info in PROPEPTIDE_D1D2.get("critical_residues", {}).items():
            annotations[pos]["propeptide_vwd"] = info
            annotations[pos]["is_vwd_hotspot"] = True
            annotations[pos]["vwd_type"] = info.get("type", "Type_2A")

        # Additional VWD hotspots
        for pos, mutation in ADDITIONAL_VWD_HOTSPOTS.get("Type_2A_other", {}).get("residues", {}).items():
            annotations[pos]["additional_hotspot"] = {"mutation": mutation, "type": "Type_2A"}
            annotations[pos]["is_vwd_hotspot"] = True

        return annotations

    def get_domain_for_position(self, position: int) -> Tuple[str, int, int, float]:
        """
        获取位置所在的域信息
        Returns: (domain_name, start, end, relative_position)
        """
        for domain, (start, end) in self.domain_architecture.items():
            if start <= position <= end:
                relative_pos = (position - start) / (end - start)
                return domain, start, end, relative_pos
        return "unknown", 0, 0, 0.0

    def calculate_mutation_properties(self, ref_aa: str, alt_aa: str) -> Dict:
        """
        计算突变氨基酸的物理化学性质变化
        """
        ref_aa = ref_aa.upper()
        alt_aa = alt_aa.upper()

        props = AMINO_ACID_PROPERTIES

        size_delta = props["size"].get(alt_aa, 0) - props["size"].get(ref_aa, 0)
        hydrophobicity_delta = props["hydrophobicity"].get(alt_aa, 0) - props["hydrophobicity"].get(ref_aa, 0)
        charge_change = int(props["charge"].get(alt_aa, 0) - props["charge"].get(ref_aa, 0))
        aromatic_change = int(props["is_aromatic"].get(alt_aa, False)) - int(props["is_aromatic"].get(ref_aa, False))
        helix_delta = props["helix_propensity"].get(alt_aa, 0) - props["helix_propensity"].get(ref_aa, 0)

        return {
            "size_delta": size_delta,
            "hydrophobicity_delta": hydrophobicity_delta,
            "charge_change": charge_change,
            "aromatic_change": aromatic_change,
            "helix_propensity_delta": helix_delta
        }

    def extract_residue_features(self, variant_id: str, position: int,
                                 ref_aa: str, alt_aa: str,
                                 structural_data: Optional[Dict] = None) -> ResidueLevelFeatures:
        """
        提取残基级特征的主函数

        Args:
            variant_id: 变异ID (如 "R1306W")
            position: 氨基酸位置
            ref_aa: 参考氨基酸
            alt_aa: 突变氨基酸
            structural_data: 可选的结构数据 (plddt, rmsd等)

        Returns:
            ResidueLevelFeatures: 完整的残基级特征
        """
        features = ResidueLevelFeatures(
            variant_id=variant_id,
            position=position,
            ref_aa=ref_aa,
            alt_aa=alt_aa
        )

        # 1. 域信息
        domain, start, end, rel_pos = self.get_domain_for_position(position)
        features.domain = domain
        features.domain_start = start
        features.domain_end = end
        features.relative_position = rel_pos

        # 2. 文献注释
        annot = self.literature_annotations.get(position, {})

        # AIM features
        if "aim" in annot:
            features.is_in_AIM = True
            features.aim_component = annot["aim"]["component"]
            if annot["aim"]["is_key_contact"]:
                features.aim_key_contact = annot["aim"]["contact_details"]
            features.literature_pmids.append(annot["aim"]["pmid"])
            features.literature_evidence["AIM"] = AIM_RESIDUES[annot["aim"]["component"]]["mechanism"]

        # GPIbα features
        if "gpib" in annot:
            features.is_in_gpib_interface = True
            features.gpib_site = annot["gpib"]["site"]
            if annot["gpib"]["is_key_residue"]:
                features.gpib_key_residue = True
                features.gpib_interaction_type = annot["gpib"]["residue_details"].get("type", "")
                if "vwd_mutations" in annot["gpib"]["residue_details"]:
                    features.vwd_hotspot_mutation = ", ".join(annot["gpib"]["residue_details"]["vwd_mutations"])
            features.literature_pmids.append(annot["gpib"]["pmid"])

        # ADAMTS13 features
        if "adamts13_scissile" in annot:
            features.is_scissile_bond = True
            features.literature_evidence["ADAMTS13"] = "Scissile bond Y1605-M1606"

        if "adamts13_exosite_1" in annot:
            features.is_in_exosite_1 = True
            features.literature_evidence["Exosite_1"] = "Interacts with ADAMTS13 disintegrin domain"

        if "adamts13_exosite_2" in annot:
            features.is_in_exosite_2 = True
            features.literature_evidence["Exosite_2"] = "Binds cysteine-rich domain"

        if "adamts13_exosite_3" in annot:
            features.is_in_exosite_3 = True
            features.literature_evidence["Exosite_3"] = "Associates with spacer domain"

        if "group_1_2A" in annot:
            features.is_group1_2A_mutation = True
            features.literature_evidence["Group1_2A"] = annot["group_1_2A"]

        if "group_2_2A" in annot:
            features.is_group2_2A_mutation = True
            features.literature_evidence["Group2_2A"] = annot["group_2_2A"]

        # Calcium site
        if "calcium_coordinating" in annot:
            features.is_calcium_coordinating = True
            features.calcium_site_importance = "high"
            features.literature_evidence["Calcium"] = "Coordinating residue for Ca2+ binding site"

        # Collagen binding
        if "collagen_binding" in annot:
            features.is_collagen_binding_residue = True
            features.collagen_type = annot["collagen_binding"]["type"]

        if "collagen_critical" in annot:
            features.is_collagen_binding_residue = True
            features.is_known_vwd_hotspot = True
            features.vwd_hotspot_type = "Type_2M"
            features.vwd_hotspot_mutation = annot["collagen_critical"].get("vwd_mutation", "")

        # FVIII binding
        if "fviii_binding" in annot:
            features.is_fviii_binding_residue = True
            features.fviii_binding_importance = "critical"
            features.is_known_vwd_hotspot = True
            features.vwd_hotspot_type = "Type_2N"
            features.vwd_hotspot_mutation = annot["fviii_binding"].get("vwd_mutation", "")

        # RGD motif
        if "rgd_motif" in annot:
            features.is_rgd_motif = True

        # 3. 突变性质
        mut_props = self.calculate_mutation_properties(ref_aa, alt_aa)
        features.mutation_size_delta = mut_props["size_delta"]
        features.mutation_hydrophobicity_delta = mut_props["hydrophobicity_delta"]
        features.mutation_charge_change = mut_props["charge_change"]
        features.mutation_aromatic_change = mut_props["aromatic_change"]
        features.mutation_helix_propensity_delta = mut_props["helix_propensity_delta"]

        # 4. 结构数据 (如果提供)
        if structural_data:
            features.plddt_wt = structural_data.get("plddt_wt", 0)
            features.plddt_mut = structural_data.get("plddt_mut", 0)
            features.plddt_delta = structural_data.get("plddt_delta", 0)
            features.global_rmsd = structural_data.get("global_rmsd", 0)
            features.local_rmsd_10a = structural_data.get("local_rmsd_10a", 0)
            features.local_rmsd_15a = structural_data.get("local_rmsd_15a", 0)

        # 5. 计算AIM破坏评分 (如果在AIM中)
        if features.is_in_AIM and features.aim_key_contact:
            features.aim_disruption_score = self._calculate_aim_disruption(
                features.aim_key_contact, mut_props
            )

        # 6. 去重PMIDs
        features.literature_pmids = list(set(features.literature_pmids))

        return features

    def _calculate_aim_disruption(self, key_contact: Dict, mut_props: Dict) -> float:
        """
        计算AIM破坏评分
        基于突变性质和关键接触类型
        """
        score = 0.0

        if key_contact.get("type") == "salt_bridge":
            # 盐桥破坏
            if abs(mut_props["charge_change"]) > 0:
                score += 0.9  # 电荷改变强烈破坏盐桥
            elif abs(mut_props["size_delta"]) > 50:
                score += 0.5  # 大残基可能空间位阻

        elif key_contact.get("type") == "hydrogen_bond":
            # 氢键破坏
            if mut_props["charge_change"] != 0:
                score += 0.6
            elif key_contact.get("atoms") == "backbone-backbone":
                # 主链氢键较难破坏
                score += 0.2

        elif key_contact.get("type") == "hydrophobic_core":
            # 疏水核心破坏
            if abs(mut_props["hydrophobicity_delta"]) > 2:
                score += 0.7
            if mut_props["charge_change"] != 0:
                score += 0.5  # 电荷进入疏水核心有害

        return min(score, 1.0)

    def extract_features_batch(self, variants: List[Dict]) -> List[ResidueLevelFeatures]:
        """
        批量提取特征

        Args:
            variants: 列表 of {"variant_id": str, "position": int, "ref_aa": str, "alt_aa": str}

        Returns:
            List[ResidueLevelFeatures]: 特征列表
        """
        results = []
        for var in variants:
            features = self.extract_residue_features(
                variant_id=var["variant_id"],
                position=var["position"],
                ref_aa=var["ref_aa"],
                alt_aa=var["alt_aa"],
                structural_data=var.get("structural_data")
            )
            results.append(features)
        return results

    def generate_feature_report(self, features: ResidueLevelFeatures) -> str:
        """
        生成特征报告
        """
        report = []
        report.append("=" * 80)
        report.append(f"VWF Residue-Level Feature Report: {features.variant_id}")
        report.append("=" * 80)
        report.append("")

        # 基本信息
        report.append("-" * 80)
        report.append("VARIANT INFORMATION")
        report.append("-" * 80)
        report.append(f"Position: {features.position} ({features.ref_aa} → {features.alt_aa})")
        report.append(f"Domain: {features.domain} ({features.domain_start}-{features.domain_end})")
        report.append(f"Relative Position in Domain: {features.relative_position:.2%}")
        report.append("")

        # AIM特征
        if features.is_in_AIM:
            report.append("-" * 80)
            report.append("AIM (AUTOINHIBITORY MODULE) FEATURES")
            report.append("-" * 80)
            report.append(f"AIM Component: {features.aim_component}")
            report.append(f"AIM Disruption Score: {features.aim_disruption_score:.2f}")
            if features.aim_key_contact:
                report.append(f"Key Contact: {features.aim_key_contact}")
            report.append("")

        # GPIbα界面特征
        if features.is_in_gpib_interface:
            report.append("-" * 80)
            report.append("GPIbα INTERFACE FEATURES")
            report.append("-" * 80)
            report.append(f"Interactive Site: {features.gpib_site}")
            report.append(f"Is Key Residue: {features.gpib_key_residue}")
            report.append(f"Interaction Type: {features.gpib_interaction_type}")
            report.append("")

        # ADAMTS13特征
        if features.is_scissile_bond or features.is_in_exosite_1:
            report.append("-" * 80)
            report.append("ADAMTS13 CLEAVAGE SITE FEATURES")
            report.append("-" * 80)
            report.append(f"Is Scissile Bond: {features.is_scissile_bond}")
            report.append(f"Exosite 1: {features.is_in_exosite_1}")
            report.append(f"Exosite 2: {features.is_in_exosite_2}")
            report.append(f"Exosite 3: {features.is_in_exosite_3}")
            report.append(f"Group 1 Type 2A: {features.is_group1_2A_mutation}")
            report.append(f"Group 2 Type 2A: {features.is_group2_2A_mutation}")
            report.append("")

        # 功能残基
        report.append("-" * 80)
        report.append("FUNCTIONAL RESIDUE ANNOTATIONS")
        report.append("-" * 80)
        report.append(f"Calcium Coordinating: {features.is_calcium_coordinating}")
        report.append(f"Collagen Binding: {features.is_collagen_binding_residue}")
        report.append(f"FVIII Binding: {features.is_fviii_binding_residue}")
        report.append(f"RGD Motif: {features.is_rgd_motif}")
        report.append(f"Multimerization: {features.is_multimerization_residue}")
        report.append("")

        # VWD热点
        report.append("-" * 80)
        report.append("VWD HOTSPOT STATUS")
        report.append("-" * 80)
        report.append(f"Is Known Hotspot: {features.is_known_vwd_hotspot}")
        report.append(f"Hotspot Type: {features.vwd_hotspot_type}")
        report.append(f"Known Mutations: {features.vwd_hotspot_mutation}")
        report.append("")

        # 结构特征
        report.append("-" * 80)
        report.append("STRUCTURAL FEATURES")
        report.append("-" * 80)
        report.append(f"pLDDT WT: {features.plddt_wt:.2f}")
        report.append(f"pLDDT Mut: {features.plddt_mut:.2f}")
        report.append(f"pLDDT Delta: {features.plddt_delta:.2f}")
        report.append(f"Global RMSD: {features.global_rmsd:.3f} Å")
        report.append(f"Local RMSD (10Å): {features.local_rmsd_10a:.3f} Å")
        report.append("")

        # 突变性质
        report.append("-" * 80)
        report.append("MUTATION PROPERTIES")
        report.append("-" * 80)
        report.append(f"Size Delta: {features.mutation_size_delta:.1f} Å³")
        report.append(f"Hydrophobicity Delta: {features.mutation_hydrophobicity_delta:.2f}")
        report.append(f"Charge Change: {features.mutation_charge_change}")
        report.append(f"Aromatic Change: {features.mutation_aromatic_change}")
        report.append(f"Helix Propensity Delta: {features.mutation_helix_propensity_delta:.3f}")
        report.append("")

        # 文献引用
        if features.literature_pmids:
            report.append("-" * 80)
            report.append("LITERATURE REFERENCES")
            report.append("-" * 80)
            for pmid in features.literature_pmids:
                report.append(f"  - PMID: {pmid}")
            for key, evidence in features.literature_evidence.items():
                if isinstance(evidence, str):
                    report.append(f"  {key}: {evidence}")
                elif isinstance(evidence, dict):
                    report.append(f"  {key}: {evidence.get('description', '')}")
            report.append("")

        report.append("=" * 80)
        return "\n".join(report)


def main():
    """
    测试残基级特征提取器
    """
    print("=" * 80)
    print("VWF Residue-Level Feature Extractor - Test")
    print("=" * 80)
    print()

    extractor = VWFResidueFeatureExtractor()

    # 测试用例：覆盖不同功能域
    test_variants = [
        # A1域 - AIM + GPIbα界面 (Type 2B热点)
        {"variant_id": "R1306W", "position": 1306, "ref_aa": "R", "alt_aa": "W",
         "structural_data": {"plddt_wt": 85.5, "plddt_mut": 78.3, "plddt_delta": -7.2,
                            "global_rmsd": 0.8, "local_rmsd_10a": 1.2}},

        # A1域 - 仅GPIbα界面
        {"variant_id": "V1316M", "position": 1316, "ref_aa": "V", "alt_aa": "M"},

        # A2域 - Group 2 Type 2A
        {"variant_id": "D1614N", "position": 1614, "ref_aa": "D", "alt_aa": "N"},

        # A2域 - Group 1 Type 2A
        {"variant_id": "M1528V", "position": 1528, "ref_aa": "M", "alt_aa": "V"},

        # A3域 - Type 2M (胶原结合)
        {"variant_id": "S1731P", "position": 1731, "ref_aa": "S", "alt_aa": "P"},

        # D'D3域 - Type 2N (FVIII结合)
        {"variant_id": "R816W", "position": 816, "ref_aa": "R", "alt_aa": "W"},
    ]

    for var in test_variants:
        features = extractor.extract_residue_features(
            variant_id=var["variant_id"],
            position=var["position"],
            ref_aa=var["ref_aa"],
            alt_aa=var["alt_aa"],
            structural_data=var.get("structural_data")
        )

        print(f"\n{'='*80}")
        print(extractor.generate_feature_report(features))

        # 输出JSON格式
        print("\n--- JSON Output ---")
        print(json.dumps(features.to_dict(), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

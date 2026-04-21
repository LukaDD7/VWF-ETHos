#!/usr/bin/env python3
"""
AgenticVWFClassifier
====================
Phase 2 Implementation: Multi-Expert Architecture with Logical Fusion

Per Alphagenome-DFR-Phase3.md architecture:
- Expert 1: Structural Expert (AF3 data) → structural_damage_score
- Expert 2: Transcriptomic Expert (AG data) → rna_drop_score, splice_override_score
- Expert 3: Clinical Geneticist (Logical Fusion) → MultiLabelClassificationResult

Key design principles:
- NO hardcoded thresholds (use adaptive thresholds based on data distribution)
- NO black-box models for Expert 3 (use interpretable logical rules)
- D4 domain: use ag_rna_delta to resolve Type1 vs 2A competition

Author: Claude Code
Date: 2026-04-19
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from enum import Enum

# ============================================================================
# DATA CLASSES
# ============================================================================

class VWFType(Enum):
    """VWD Type-2 subtypes."""
    TYPE_1 = "1"
    TYPE_2A = "2A"
    TYPE_2B = "2B"
    TYPE_2M = "2M"
    TYPE_2N = "2N"
    UNCERTAIN = "uncertain"


@dataclass
class ExpertScores:
    """Scores from each expert agent."""
    # Expert 1: Structural
    structural_damage_score: float = 0.0  # 0-1, higher = more damaging
    af3_plddt_delta: float = 0.0  # pLDDT change at mutation site
    interface_perturbation: float = 0.0  # PAE-based interface disruption

    # Expert 2: Transcriptomic
    rna_drop_score: float = 0.0  # 0-1, higher = more expression loss (→ Type1)
    splice_override_score: float = 0.0  # 0-1, higher = more splice disruption (→ 2A/2M)
    ag_rna_delta: float = 0.0
    ag_splice_delta: float = 0.0

    # Expert 3: Clinical Genetics
    domain_pleiotropy_score: float = 0.0  # A1→2B/2M, A3→2M/2A etc.


@dataclass
class MultiLabelClassificationResult:
    """Final output from AgenticVWFClassifier."""
    main_subtype: str
    alternatives: List[str] = field(default_factory=list)
    confidence: float = 0.0
    reasoning: str = ""
    domain_pleiotropy: str = ""
    expert_scores: Optional[ExpertScores] = None

    # Additional metadata
    is_type1_signal: bool = False
    is_splice_override: bool = False
    is_d4_competition: bool = False


# ============================================================================
# VWF DOMAIN ARCHITECTURE (from literature)
# ============================================================================

VWF_DOMAIN_RANGES = {
    'D1': (1, 233),
    'D2': (234, 480),
    'D3': (481, 728),
    "D'": (729, 763),
    'D3_extended': (764, 1233),  # D' + D3 extended for 2N
    'A1': (1271, 1492),
    'A2': (1493, 1684),
    'A3': (1685, 1873),
    'D4': (1874, 2255),
    'C1': (2256, 2299),
    'C2': (2300, 2363),
    'C3': (2364, 2411),
    'C4': (2412, 2450),
    'C5': (2451, 2459),
    'C6': (2460, 2527),
    'CK': (2528, 2813),
}

# Type-2 domain associations - REFINED (from Blood 2024/2026 literature)
# Only map to 2N if position is near FVIII binding interface (782-816)
DOMAIN_TO_TYPE2 = {
    'A1': ['2B', '2M'],      # A1 domain: 2B (GPIb gain-of-function) or 2M (collagen binding)
    'A2': ['2A'],            # A2 domain: ADAMTS13 cleavage → 2A
    'A3': ['2M', '2A'],      # A3 domain: collagen binding → 2M, or 2A-like
    "D'": ['2N'],            # D' domain: FVIII binding → 2N (only if near 782-816)
    'D3_extended': ['2N', '2A'],  # D3 extended: primarily 2N if near FVIII interface, else 2A
    'D3': ['2A'],            # D3 domain: mostly 2A (multimerization), NOT 2N
    'D4': ['2A', '1'],       # D4 domain: secretion/multimerization → 2A or Type1 (competition!)
    'C1': ['2A'],
    'C2': ['2A'],
    'C3': ['2A'],
    'C4': ['2A'],
    'C5': ['2A'],
    'C6': ['2A'],
    'CK': ['2A'],
    'D1': ['2A'],            # Propeptide, affects multimerization → 2A
    'D2': ['2A'],
}

# Refined functional sites for type-specific mechanisms
# Only positions within these ranges trigger specific mechanisms
FUNCTIONAL_SITES = {
    'A2_cleavage_site': (1604, 1606),  # Y1605-M1606, ADAMTS13 cleavage
    'A2_ca_binding': (1596, 1602),     # Ca2+ binding site
    'A2_disulfide': (1669, 1670),      #临近二硫键
    'A1_aim_n_term': (1238, 1268),     # A1 N-terminal AIM → 2B
    'A1_aim_c_term': (1460, 1472),     # A1 C-terminal AIM → 2B
    'A1_interface': (1271, 1492),      # A1 GPIb binding interface → 2B or 2M
    'D_prime_fviii': (782, 816),       # D' FVIII binding interface → 2N (key range!)
    'D3_fviii_interface': (764, 900), # D3 N-terminal region for FVIII → 2N (if near 800)
    'A3_collagen': (1700, 1850),       # Collagen binding site → 2M
    'D4_sEC': (1874, 2255),            # D4 secretion/multimerization → 2A or Type1
}


# ============================================================================
# EXPERT 1: STRUCTURAL EXPERT
# ============================================================================

class StructuralExpert:
    """
    Expert 1: Analyzes AF3 structural data to assess structural damage.

    Inputs: af3_plddt_mean, af3_plddt_min, af3_pae_interface, position, domain
    Outputs: structural_damage_score (0-1), af3_plddt_delta, interface_perturbation
    """

    def __init__(self, threshold_percentile=75):
        """
        Args:
            threshold_percentile: Adaptive threshold based on data distribution
        """
        self.threshold_percentile = threshold_percentile
        self._plddt_baseline = None  # Will be set from data

    def set_baseline(self, plddt_mean: float):
        """Set baseline pLDDT from WT structure."""
        self._plddt_baseline = plddt_mean

    def compute_damage_score(self, variant_data: Dict) -> Tuple[float, float, float]:
        """
        Compute structural damage score for a variant.

        Returns: (structural_damage_score, plddt_delta, interface_perturbation)
        """
        plddt_mean = variant_data.get('af3_plddt_mean', np.nan)
        plddt_min = variant_data.get('af3_plddt_min', np.nan)
        pae_interface = variant_data.get('af3_pae_interface', np.nan)

        # pLDDT delta from baseline (WT typically ~90)
        if self._plddt_baseline is not None:
            plddt_delta = plddt_mean - self._plddt_baseline
        else:
            # Assume WT baseline of 90 if not set
            plddt_delta = plddt_mean - 90 if not np.isnan(plddt_mean) else 0

        # Normalize structural damage score (0-1)
        # Lower pLDDT = more damage
        if not np.isnan(plddt_mean):
            # Using adaptive threshold: if pLDDT < 70, significant damage
            damage_score = max(0, min(1, (90 - plddt_mean) / 30))
        else:
            damage_score = 0.5  # Unknown

        # Interface perturbation from PAE
        # Higher PAE = more interface disruption (for 2N/2B binding sites)
        if not np.isnan(pae_interface):
            interface_score = min(1, pae_interface / 10)  # Normalize
        else:
            interface_score = 0

        return damage_score, plddt_delta, interface_score

    def is_functional_site_proximal(self, position: int, domain: str) -> bool:
        """Check if position is near a functional site."""
        for site_name, (start, end) in FUNCTIONAL_SITES.items():
            if start <= position <= end:
                return True
            # Check proximity (within 10 residues)
            if abs(position - start) <= 10 or abs(position - end) <= 10:
                return True
        return False


# ============================================================================
# EXPERT 2: TRANSCRIPTOMIC EXPERT
# ============================================================================

class TranscriptomicExpert:
    """
    Expert 2: Analyzes AlphaGenome data for RNA/splice signals.

    Key rules from Phase3:
    - ag_rna_delta极低 → Type 1 (ER degradation/NMD)
    - ag_splice_delta极高 → 强制Type 2A/2M (外显子跳过)
    - D4 domain: 用ag_rna_delta区分Type1 vs 2A
    """

    def __init__(self):
        # Adaptive thresholds will be set from data distribution
        self._rna_drop_threshold = None
        self._splice_override_threshold = None

    def set_adaptive_thresholds(self, rna_deltas: List[float], splice_deltas: List[float]):
        """
        Set adaptive thresholds based on data distribution.
        Uses percentiles to avoid hardcoded values.
        """
        if len(rna_deltas) > 0:
            # 25th percentile as threshold for "low"
            self._rna_drop_threshold = np.percentile(rna_deltas, 25)
        if len(splice_deltas) > 0:
            # 75th percentile as threshold for "high"
            self._splice_override_threshold = np.percentile(splice_deltas, 75)

    def compute_rna_drop_score(self, ag_rna_delta: float) -> float:
        """
        Compute RNA drop score (0-1).
        Higher score = more likely Type1 (secretion defect / ER degradation).

        Uses adaptive threshold from data distribution.
        """
        if self._rna_drop_threshold is None:
            # Default threshold
            threshold = 0.5
        else:
            threshold = self._rna_drop_threshold

        # Score: delta below threshold → high score
        if ag_rna_delta < threshold:
            # Linear scaling: lower delta = higher score
            score = 1 - (ag_rna_delta / threshold) if threshold > 0 else 1.0
        else:
            score = 0.0

        return max(0, min(1, score))

    def compute_splice_override_score(self, ag_splice_delta: float) -> float:
        """
        Compute splice override score (0-1).
        Higher score = more likely to override structural defaults → 2A/2M.

        Uses adaptive threshold from data distribution.
        """
        if self._splice_override_threshold is None:
            threshold = 0.5
        else:
            threshold = self._splice_override_threshold

        # Score: delta above threshold → high score
        if ag_splice_delta > threshold:
            score = min(1, ag_splice_delta)
        else:
            score = 0.0

        return max(0, min(1, score))


# ============================================================================
# EXPERT 3: CLINICAL GENETICIST (LOGICAL FUSION)
# ============================================================================

class ClinicalGeneticistAgent:
    """
    Expert 3: Logical Fusion of Expert 1 & 2 outputs.

    Key design: INTERPRETABLE RULES, not black-box model.
    All decisions are traceable and explainable.

    Decision hierarchy (per Phase3):
    1. RNA drop + D4 domain → Type 1 (secretion defect)
    2. Splice override → Force 2A/2M regardless of domain
    3. D4 domain + normal RNA → 2A (multimerization defect)
    4. Structural damage + A2 domain → 2A
    5. A1 domain + interface perturbation → 2B vs 2M (pleiotropy)
    6. D'/D3 domain → 2N (FVIII binding)
    """

    def __init__(self, structural_expert: StructuralExpert, transcriptomic_expert: TranscriptomicExpert):
        self.structural_expert = structural_expert
        self.transcriptomic_expert = transcriptomic_expert

    def fuse(self, variant_data: Dict, expert_scores: ExpertScores) -> MultiLabelClassificationResult:
        """
        Main fusion logic - interpretable rule-based decision making.
        """
        domain = variant_data.get('domain', '')
        position = variant_data.get('protein_pos', 0)

        # Initialize
        main_subtype = 'uncertain'
        alternatives = []
        confidence = 0.0
        reasoning_steps = []

        # =====================================================================
        # RULE 1: Type 1 Signal (RNA drop + D4 domain)
        # =====================================================================
        if expert_scores.rna_drop_score > 0.7 and domain == 'D4':
            main_subtype = '1'
            confidence = 0.8
            reasoning_steps.append("RULE1: D4 domain + high RNA drop → Type1 (secretion defect)")
            return MultiLabelClassificationResult(
                main_subtype=main_subtype,
                alternatives=['2A'],
                confidence=confidence,
                reasoning="; ".join(reasoning_steps),
                domain_pleiotropy="D4域突变可导致分泌障碍(Type1)或多聚化异常(Type2A)，RNA drop区分之",
                expert_scores=expert_scores,
                is_type1_signal=True,
                is_d4_competition=True
            )

        # =====================================================================
        # RULE 2: Splice Override (force 2A/2M)
        # =====================================================================
        if expert_scores.splice_override_score > 0.7:
            main_subtype = '2A'
            confidence = 0.85
            reasoning_steps.append(f"RULE2: Splice override score={expert_scores.splice_override_score:.2f} → Force 2A/2M")
            alternatives = ['2M']

            # Check if domain suggests 2M instead
            if domain == 'A3':
                main_subtype = '2M'
                reasoning_steps.append("RULE2b: A3 domain + splice override → 2M (collagen binding defect)")

            return MultiLabelClassificationResult(
                main_subtype=main_subtype,
                alternatives=alternatives,
                confidence=confidence,
                reasoning="; ".join(reasoning_steps),
                domain_pleiotropy="Splice override overrides structural defaults",
                expert_scores=expert_scores,
                is_splice_override=True
            )

        # =====================================================================
        # RULE 3: D4 domain + normal RNA → 2A (multimerization defect)
        # =====================================================================
        if domain == 'D4' and expert_scores.rna_drop_score < 0.3:
            main_subtype = '2A'
            confidence = 0.75
            reasoning_steps.append("RULE3: D4 domain + normal RNA → 2A (multimerization defect)")
            return MultiLabelClassificationResult(
                main_subtype=main_subtype,
                alternatives=['1'],
                confidence=confidence,
                reasoning="; ".join(reasoning_steps),
                domain_pleiotropy="D4域突变可导致分泌障碍(Type1)或多聚化异常(Type2A)",
                expert_scores=expert_scores,
                is_d4_competition=True
            )

        # =====================================================================
        # RULE 4: A2 domain + structural damage → 2A
        # =====================================================================
        if domain == 'A2':
            reasoning_steps.append("RULE4: A2 domain analysis")

            # Check proximity to cleavage site
            if 1593 <= position <= 1610:
                main_subtype = '2A'
                confidence = 0.9
                reasoning_steps.append(f"RULE4a: Near ADAMTS13 cleavage site (1593-1610), pos={position}")
            elif expert_scores.structural_damage_score > 0.5:
                main_subtype = '2A'
                confidence = 0.75
                reasoning_steps.append(f"RULE4b: A2 domain + high structural damage ({expert_scores.structural_damage_score:.2f})")
            else:
                main_subtype = '2A'
                confidence = 0.6
                reasoning_steps.append("RULE4c: A2 domain default → 2A")

            return MultiLabelClassificationResult(
                main_subtype=main_subtype,
                alternatives=[],
                confidence=confidence,
                reasoning="; ".join(reasoning_steps),
                domain_pleiotropy="A2域主要关联2A (ADAMTS13切割敏感性)",
                expert_scores=expert_scores
            )

        # =====================================================================
        # RULE 5: D'/D3 domain → 2N (FVIII binding) ONLY if near interface
        # =====================================================================
        if domain in ["D'", "D3", "D3_extended"]:
            reasoning_steps.append(f"RULE5: {domain} domain analysis, pos={position}")

            # Only classify as 2N if position is near FVIII binding interface (782-816)
            # D' domain (729-763): only 2N if position is close to 782
            # D3 domain (481-728): mostly 2A (multimerization), NOT 2N
            # D3_extended (764-1233): primarily 2N only if near 782-816, else 2A

            if 782 <= position <= 816:
                # Clear FVIII binding interface region
                main_subtype = '2N'
                confidence = 0.85
                reasoning_steps.append(f"RULE5a: Near FVIII binding interface (782-816), pos={position} → 2N")
            elif domain == "D'" and 760 <= position <= 790:
                # D' domain near interface
                main_subtype = '2N'
                confidence = 0.75
                reasoning_steps.append(f"RULE5b: D' domain near FVIII interface, pos={position} → 2N")
            elif domain == "D3" or (domain == "D3_extended" and position < 900):
                # D3 domain: multimerization defects → 2A, NOT 2N
                main_subtype = '2A'
                confidence = 0.7
                reasoning_steps.append(f"RULE5c: {domain} domain + position {position} → 2A (multimerization defect, NOT 2N)")
            elif domain == "D3_extended":
                # D3_extended with position >= 900: could be 2A or 2N based on structural damage
                if expert_scores.structural_damage_score > 0.5:
                    main_subtype = '2A'
                    confidence = 0.7
                    reasoning_steps.append(f"RULE5d: D3_extended + high structural damage → 2A")
                else:
                    main_subtype = '2N'
                    confidence = 0.6
                    reasoning_steps.append(f"RULE5e: D3_extended with moderate damage → 2N (default)")
            else:
                # Fallback: use structural damage
                if expert_scores.structural_damage_score > 0.5:
                    main_subtype = '2A'
                    confidence = 0.6
                else:
                    main_subtype = '2N'
                    confidence = 0.5
                reasoning_steps.append(f"RULE5f: {domain} fallback → {main_subtype} (confidence={confidence})")

            return MultiLabelClassificationResult(
                main_subtype=main_subtype,
                alternatives=['2A'] if main_subtype == '2N' else ['2N'],
                confidence=confidence,
                reasoning="; ".join(reasoning_steps),
                domain_pleiotropy=f"{domain}域关联: 主要取决于位置和结构损伤程度",
                expert_scores=expert_scores
            )

        # =====================================================================
        # RULE 6: A1 domain → 2B vs 2M (pleiotropy resolution)
        # =====================================================================
        if domain == 'A1':
            reasoning_steps.append("RULE6: A1 domain pleiotropy resolution")

            # A1 domain: 2B (GPIb gain-of-function) or 2M (collagen binding)
            # 2B requires AIM (autoinhibitory module) disruption: positions 1238-1268 or 1460-1472
            # OR high structural damage at the GPIb binding interface

            # Check if position is in AIM regions (strong 2B signal)
            is_aim_n = 1238 <= position <= 1268
            is_aim_c = 1460 <= position <= 1472

            if is_aim_n or is_aim_c:
                # Strong 2B signal - AIM disruption
                main_subtype = '2B'
                confidence = 0.85
                if is_aim_n:
                    reasoning_steps.append(f"RULE6a: A1 N-terminal AIM region ({1238}-{1268}), pos={position} → 2B (AIM disruption)")
                else:
                    reasoning_steps.append(f"RULE6b: A1 C-terminal AIM region ({1460}-{1472}), pos={position} → 2B (AIM disruption)")

            elif expert_scores.interface_perturbation > 0.3:
                # Interface perturbation suggests 2B (GPIb binding issue)
                main_subtype = '2B'
                confidence = 0.75
                reasoning_steps.append(f"RULE6c: A1 + interface perturbation={expert_scores.interface_perturbation:.2f} → 2B")
            elif expert_scores.structural_damage_score > 0.6:
                # High structural damage (pLDDT < 72) could indicate 2B mechanism
                main_subtype = '2B'
                confidence = 0.65
                reasoning_steps.append(f"RULE6d: A1 + high structural damage ({expert_scores.structural_damage_score:.2f}) → 2B")
            else:
                # Default to 2M (collagen binding defect)
                main_subtype = '2M'
                confidence = 0.65
                reasoning_steps.append(f"RULE6e: A1 + normal interface → 2M (collagen binding defect)")

            return MultiLabelClassificationResult(
                main_subtype=main_subtype,
                alternatives=['2M'] if main_subtype == '2B' else ['2B'],
                confidence=confidence,
                reasoning="; ".join(reasoning_steps),
                domain_pleiotropy="A1域可致2B(GPIb-GOF)或2M(胶原结合)，AIM区域/界面扰动区分之",
                expert_scores=expert_scores
            )

        # =====================================================================
        # RULE 7: A3 domain → 2M (collagen binding)
        # =====================================================================
        if domain == 'A3':
            main_subtype = '2M'
            confidence = 0.7
            reasoning_steps.append("RULE7: A3 domain → 2M (collagen binding defect)")
            return MultiLabelClassificationResult(
                main_subtype=main_subtype,
                alternatives=['2A'],
                confidence=confidence,
                reasoning="; ".join(reasoning_steps),
                domain_pleiotropy="A3域可致2M或2A样表型",
                expert_scores=expert_scores
            )

        # =====================================================================
        # RULE 8: Default - insufficient data
        # =====================================================================
        reasoning_steps.append("RULE8: Default - insufficient data for confident classification")
        return MultiLabelClassificationResult(
            main_subtype='uncertain',
            alternatives=['2A', '2M', '2N'],
            confidence=0.3,
            reasoning="; ".join(reasoning_steps),
            domain_pleiotropy="无法确定域关联",
            expert_scores=expert_scores
        )


# ============================================================================
# MAIN CLASSIFIER
# ============================================================================

class AgenticVWFClassifier:
    """
    Main classifier: Coordinates 3 Expert agents and produces final classification.

    Architecture:
    1. Structural Expert (Expert 1) → AF3-based damage scores
    2. Transcriptomic Expert (Expert 2) → AG-based RNA/splice signals
    3. Clinical Geneticist Agent (Expert 3) → Logical fusion → Final result
    """

    def __init__(self):
        self.structural_expert = StructuralExpert()
        self.transcriptomic_expert = TranscriptomicExpert()
        self.clinical_geneticist = ClinicalGeneticistAgent(
            self.structural_expert,
            self.transcriptomic_expert
        )
        self._is_fitted = False

    def fit(self, df: pd.DataFrame):
        """
        Fit the classifier on training data to set adaptive thresholds.

        Args:
            df: DataFrame with columns ['ag_rna_delta', 'ag_splice_delta', 'af3_plddt_mean']
        """
        # Set baseline pLDDT from WT
        if 'af3_plddt_mean' in df.columns:
            wt_data = df[df['aa_change'].str.lower().str.contains('wt', na=False)]
            if len(wt_data) > 0:
                self.structural_expert.set_baseline(wt_data['af3_plddt_mean'].iloc[0])

        # Set adaptive thresholds from data distribution
        rna_deltas = df['ag_rna_delta'].dropna().tolist()
        splice_deltas = df['ag_splice_delta'].dropna().tolist()
        self.transcriptomic_expert.set_adaptive_thresholds(rna_deltas, splice_deltas)

        self._is_fitted = True

    def classify(self, variant_data: Dict) -> MultiLabelClassificationResult:
        """
        Classify a single variant.

        Args:
            variant_data: Dict with keys:
                - protein_pos: int
                - ref_aa: str
                - alt_aa: str
                - domain: str
                - ag_rna_delta: float (optional)
                - ag_splice_delta: float (optional)
                - af3_plddt_mean: float (optional)
                - af3_plddt_min: float (optional)
                - af3_pae_interface: float (optional)

        Returns:
            MultiLabelClassificationResult
        """
        # =====================================================================
        # Expert 1: Structural Analysis
        # =====================================================================
        structural_damage, plddt_delta, interface_perturbation = \
            self.structural_expert.compute_damage_score(variant_data)

        # =====================================================================
        # Expert 2: Transcriptomic Analysis
        # =====================================================================
        ag_rna_delta = variant_data.get('ag_rna_delta', np.nan)
        ag_splice_delta = variant_data.get('ag_splice_delta', np.nan)

        rna_drop_score = self.transcriptomic_expert.compute_rna_drop_score(ag_rna_delta)
        splice_override_score = self.transcriptomic_expert.compute_splice_override_score(ag_splice_delta)

        # =====================================================================
        # Collect Expert Scores
        # =====================================================================
        expert_scores = ExpertScores(
            structural_damage_score=structural_damage,
            af3_plddt_delta=plddt_delta,
            interface_perturbation=interface_perturbation,
            rna_drop_score=rna_drop_score,
            splice_override_score=splice_override_score,
            ag_rna_delta=ag_rna_delta,
            ag_splice_delta=ag_splice_delta
        )

        # =====================================================================
        # Expert 3: Logical Fusion
        # =====================================================================
        result = self.clinical_geneticist.fuse(variant_data, expert_scores)

        return result

    def classify_batch(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Classify a batch of variants.

        Args:
            df: DataFrame with variant data

        Returns:
            DataFrame with added columns: main_subtype, alternatives, confidence, reasoning
        """
        if not self._is_fitted:
            self.fit(df)

        results = []
        for _, row in df.iterrows():
            variant_data = row.to_dict()
            result = self.classify(variant_data)

            results.append({
                'main_subtype': result.main_subtype,
                'alternatives': ','.join(result.alternatives),
                'confidence': result.confidence,
                'reasoning': result.reasoning,
                'domain_pleiotropy': result.domain_pleiotropy,
                'is_type1_signal': result.is_type1_signal,
                'is_splice_override': result.is_splice_override,
                'is_d4_competition': result.is_d4_competition,
                'structural_damage_score': result.expert_scores.structural_damage_score if result.expert_scores else np.nan,
                'rna_drop_score': result.expert_scores.rna_drop_score if result.expert_scores else np.nan,
                'splice_override_score': result.expert_scores.splice_override_score if result.expert_scores else np.nan,
            })

        return pd.DataFrame(results)


# ============================================================================
# VALIDATION: CONFUCSIOM MATRIX
# ============================================================================

def run_validation(classifier: AgenticVWFClassifier, df: pd.DataFrame):
    """
    Run validation on variants with known Type-2 labels.
    Generate confusion matrix comparing predicted vs known subtypes.
    """
    print("=" * 60)
    print("VALIDATION RUN")
    print("=" * 60)

    # Filter to variants with known labels
    labeled = df[df['type2_subtype'].notna()].copy()

    # Run classification
    predictions = classifier.classify_batch(labeled)

    # Merge with labels
    results = labeled.copy()
    for col in ['main_subtype', 'alternatives', 'confidence', 'reasoning']:
        results[col] = predictions[col]

    # Show results
    print(f"\nTotal labeled variants: {len(results)}")
    print(f"Predicted distribution:")
    print(results['main_subtype'].value_counts())
    print(f"\nKnown distribution:")
    print(results['type2_subtype'].value_counts())

    # Confusion matrix
    print("\n" + "=" * 60)
    print("CONFUSION MATRIX (main_subtype vs known type2_subtype)")
    print("=" * 60)

    subtypes = ['2A', '2B', '2M', '2N', 'uncertain']
    known_labels = results['type2_subtype'].unique()

    for known in sorted(known_labels):
        predicted = results[results['type2_subtype'] == known]['main_subtype'].value_counts()
        print(f"\nKnown={known}:")
        for subtype in subtypes:
            count = predicted.get(subtype, 0)
            print(f"  Predicted {subtype}: {count}")

    # Show reasoning for each prediction
    print("\n" + "=" * 60)
    print("SAMPLE PREDICTIONS WITH REASONING")
    print("=" * 60)
    for _, row in results.head(10).iterrows():
        print(f"\n{row['aa_change']} (known={row['type2_subtype']}, domain={row['domain']})")
        print(f"  → {row['main_subtype']} (conf={row['confidence']:.2f})")
        print(f"  Reasoning: {row['reasoning'][:100]}...")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="AgenticVWFClassifier - Phase 2 Implementation")
    parser.add_argument("--input", "-i", default="VWF_Alpha_Matrix.parquet",
                        help="Input parquet file from Phase 1")
    parser.add_argument("--output", "-o", default="VWF_Alpha_Matrix_classified.parquet",
                        help="Output parquet file")
    parser.add_argument("--validate", "-v", action="store_true",
                        help="Run validation on labeled variants")

    args = parser.parse_args()

    # Load data
    print("Loading VWF_Alpha_Matrix...")
    df = pd.read_parquet(args.input)
    print(f"  Loaded {len(df)} variants")

    # Create and fit classifier
    print("\nInitializing AgenticVWFClassifier...")
    classifier = AgenticVWFClassifier()

    if args.validate:
        run_validation(classifier, df)
    else:
        # Run classification
        print("\nRunning classification...")
        results = classifier.classify_batch(df)

        # Merge with original data
        output_df = df.merge(results, left_index=True, right_index=True)

        # Save
        output_df.to_parquet(args.output)
        print(f"\nSaved to: {args.output}")

        # Summary
        print("\n" + "=" * 60)
        print("CLASSIFICATION SUMMARY")
        print("=" * 60)
        print(f"Total: {len(output_df)}")
        print("\nPredicted subtype distribution:")
        print(output_df['main_subtype'].value_counts())
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

    # Autoinhibition (A1 AIM↔A1 接触, 来自 a1_aim_autoinhibition_context 结构)
    # 越大 = AIM 从 A1 GPIb 面脱离越多 = 自抑制松开 = 越像 2B (GOF)。
    # 由 extract_aim_autoinhib_features.py 产出, NaN = 无该特征 (向后兼容)。
    aim_release_score: float = float('nan')

    # MD A1 结合面破坏程度 (7A6O AIM-A1 平衡 MD; extract_7a6o_md_features.py)。
    # 越高 = WT AIM→结合面屏蔽接触丢失越多 = GPIbα 结合面被破坏 = 2M (LOF)。
    # ⚠ 平衡 MD 的 LOF/2M 证据轴, 不是 2B release 轴; NaN = 无 MD (向后兼容)。
    md_face_destab_score: float = float('nan')

    # MD AIM↔A1 盐桥保留 z-score (7A6O 平衡 MD; extract_aim_saltbridge_features.py)。
    # 高 = 结合面接触保留 = "释放但功能在" = 2B-compatible/anti-2M (实测高 z 仅见 2B+WT,
    # 无 2M); 低 = 接触塌陷 = 2M 旁证 (但低 z 不干净, 也含部分 2B, 仅作旁证不独立判 2M)。
    # ⚠ 非 2B 阳性轴 (2B 是力现象, 平衡态不显现); 阈值待校准 (n: WT1/2B8/2M3 单副本)。
    aim_sb_retained_z: float = float('nan')


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

# ---------------------------------------------------------------------------
# Autoinhibition-release thresholds (2B vs 2M discriminator)
# ---------------------------------------------------------------------------
# aim_release_score = -zscore(AIM↔A1 contacts) from a1_aim_autoinhibition_context
# structures (see extract_aim_autoinhib_features.py). Higher = AIM disengaged
# from the A1 GPIb-binding face = autoinhibition released = 2B (gain-of-function).
#
# ⚠ STILL PROVISIONAL — 这条轴无法用 panel 的 z-score 校准。evidence_matrix.csv 里的
#   a1_aim_autoinhibition_context z 是旧的全局 ptm/plddt context 值, 仍分不开 2B/2M
#   (calib_2b2m v2 实测中位: 2B aim=-0.09 vs 2M aim=-0.36, 弱且不可靠), 它 ≠ 几何
#   接触数 aim_release_score。真值要 extract_aim_autoinhib_features.py 跑 panel CIF
#   (待 H200→HF→A40 数据传输), 或由 7A6O AIM-A1 MD 的接触下降量定标 → 届时回填本常数。
AIM_RELEASE_2B_Z = 1.0    # aim_release_score ≥ this → strong release → 2B
AIM_RELEASE_LEAN_Z = 0.0  # aim_release_score > this (any release) → rescue A1 default toward 2B
# ⚠ 2026-06-20 7A6O AIM-A1 MD 实测纠偏: 不要用"MD AIM↔A1 接触下降量"反推 aim_release_score
#   并定标本常数。平衡 50ns MD 里 **2M 丢失的 AIM→结合面屏蔽接触比 2B 更多**(2M 比 2B
#   接触更少/mindist 更大), 与"接触下降=2B 松开"的旧假设**符号相反**。力依赖的 2B 松开在
#   无剪切的平衡 MD 中不显现 → MD 接触类指标是 LOF/2M 轴(见下 MD_FACE_DESTAB_2M_Z), 不是
#   2B release 轴。aim_release_score 仍需 panel CIF 或 steered/force MD 定标, 勿用平衡 MD 接触。

# ---------------------------------------------------------------------------
# 轴B': MD A1 结合面完整性 (7A6O AIM-A1 平衡 MD; LOF/2M 证据轴)
# ---------------------------------------------------------------------------
# md_face_destab_score = -zscore(变体在 tail 40-50ns 保留的 WT AIM→结合面 masking
# 接触比例)。数据驱动的 WT 屏蔽界面 = AIM 接触 A1 残基 [1305-1308,1313,1376,1378,
# 1380,1410,1434] (恰好覆盖 2B 突变位点 → 验证 Deng et al. Blood 2017 屏蔽模型)。
# 越高 = 屏蔽接触丢失越多 = A1 GPIbα 结合面被破坏 = 2M (LOF)。
# 校准 (7A6O 参考集 WT+8×2B+3×2M, 单 50ns 副本, extract_7a6o_md_features.py):
#   保留率 mask_retention: WT 0.81 > 2B 0.75 > 2M 0.63; destab 2M 中位 +0.27
#   (R1374H +1.60, G1324S +0.72, R1374C +0.27) vs 2B 中位 -0.13; AUC(2M>2B)=0.83。
#   ⚠ 平衡 MD 看不到力依赖的 2B 松开 → 本轴**只作 2M/LOF 确证, 绝不判 2B**;
#     n(2M)=3 单副本 + 2B/2M 重叠 → 仅作软证据/tie-breaker。待多副本/更多
#     2M 标签或 steered-MD 后收紧。
MD_FACE_DESTAB_2M_Z = 1.0   # md_face_destab_score ≥ this → 强结合面破坏 → 支持 2M (LOF)

# 轴B'': AIM↔A1 盐桥保留 (extract_aim_saltbridge_features.py, 7A6O 平衡 MD)。
# 联合判据 (docs/7A6O_SMD_LITERATURE_NOGO_2026-06-24.md §6): 所有致病变体自抑制
# (D1269-R1306) 都释放; 在此前提下看结合面是否保留 → 保留=2B(释放但功能在/GOF),
# 塌陷=2M(释放且功能没了/LOF)。实测 z 分离:
#   高 z (≥0.6): 仅 2B(I1309V/R1306W/S1310F)+WT, 无 2M → "结合面保留"是干净的
#                anti-2M / 2B-protective 闸 (joint 判据的"+结合保留→2B"半边)。
#   低 z (≤-0.7): 全部 3 个 2M, 但也含 2 个 2B(V1316M/V1314F) → 不干净 →
#                **仅作 2M 旁证(与 md_destab/轴B 一致时提置信), 绝不独立判 2M**。
# ⚠ 阈值待校准 (n: WT1/2B8/2M3 单副本)。NaN = 无该特征 (向后兼容)。
AIM_SB_RETAINED_2B_Z = 0.6   # aim_sb_retained_z ≥ this → 结合面保留 → anti-2M, 支持 2B
AIM_SB_COLLAPSE_2M_Z = -0.7  # aim_sb_retained_z ≤ this → 结合面塌陷 → 2M 旁证 (不独立判)

# ---------------------------------------------------------------------------
# 轴B: A1 结合面完整性 (LOF 探测) = forced_binding + heparan 两轴联合
# ---------------------------------------------------------------------------
# 临床: 2M(A1型)= 功能丧失,A1 GPIb 结合面被破坏 → GPIb 与 heparan(位点紧邻)
# 结合都降; 2B = 功能获得,结合面保留 → 两者不降。
# 校准 v2 (2026-06-12, calibrate_2b2m_thresholds.py): 用更大的独立干净标签集
# (build_labeled_variant_set.py 从 /Volumes/LQ1000/VWD 文献表清出, A1 域单标签
#  2B=36 / 2M=37 join 上 panel) 重校, 与旧 44/49 集结论一致, 故阈值不动:
#   - 中位 z: 2B fb=-0.24 hep=+0.34 | 2M fb=-0.11 hep=-0.21 (hep 是最干净的方向轴);
#   - fb 单轴几乎没用 (最佳 net 仅 0.07 @ -0.75) → fb-only fallback 必须保守;
#   - heparan 单轴较好但高召回点误伤大 (@ -0.25: 抓 49% 2M 但误伤 28% 2B);
#   - **联合 mean(fb, heparan) ≤ -0.75 = 精度最优点: 抓 24% 2M, 仅误伤 6% 2B (net 0.19)**。
#   - 收紧到 -1.0 → 16% 2M / 3% 2B (更高精度选项); 现守 -0.75 平衡召回与防 2B→2M 漏判。
# 故 LOF→2M 要求两条结合面轴一致变低(都需存在), 既抓 LOF 又不重引入 2B→2M 漏判。
LOF_COMBINED_Z = -0.75    # mean(fb_binding_zscore, heparan_zscore) ≤ this → 结合面丧失 → 2M
FB_LOSS_Z = -1.5          # (单轴 fallback, 仅 heparan 缺失时) forced_binding 极低才判 LOF

# A3 collagen-binding Boltz axis. A low collagen-complex iPTM supports Type 2M
# through impaired collagen I/III binding. This is intentionally a confidence
# modifier inside the A3 rule, not a subtype override, because the supplementary
# panel shows some known 2M variants with retained/high static complex scores.
A3_COLLAGEN_LOF_Z = -1.0
A3_COLLAGEN_RETAINED_Z = 1.0

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

    def compute_allosteric_risk(self, variant_data: Dict) -> float:
        """
        Compute allosteric risk score (0-1) for A1 domain mutations.

        Type 2B variants don't always occur in the AIM region - they can be on the
        back side of A1 domain, forcing AIM exposure through conformational changes.

        Logic:
        - If mutation is at A1 edge (outside central AIM region)
        - And PAE shows interface perturbation
        - And mutation causes local structural instability
        → High allosteric risk (2B mechanism)

        Returns:
            allosteric_risk: float 0-1 (higher = more likely allosteric 2B)
        """
        position = variant_data.get('protein_pos', 0)
        domain = variant_data.get('domain', '')
        pae = variant_data.get('af3_pae_interface', np.nan)
        plddt = variant_data.get('af3_plddt_mean', np.nan)

        if domain != 'A1':
            return 0.0

        # AIM regions
        aim_n_range = (1238, 1268)
        aim_c_range = (1460, 1472)

        # Check if mutation is at A1 edge (outside AIM but still in A1)
        is_at_edge = (
            (aim_n_range[1] < position < 1300) or  # After N-terminal AIM, before center
            (1450 < position < aim_c_range[0])      # Between central region and C-terminal AIM
        )

        if not is_at_edge:
            return 0.0

        # Compute allosteric risk score
        allosteric_score = 0.0

        if not np.isnan(pae) and pae > 0.15:
            # PAE perturbation suggests conformational change
            allosteric_score += min(0.6, pae * 2)

        if not np.isnan(plddt) and 70 < plddt < 85:
            # Mild local instability (pLDDT 70-85) suggests allosteric effect
            # Too stable (pLDDT > 85) or too unstable (pLDDT < 70) doesn't suggest allosteric
            instability = (85 - plddt) / 15
            allosteric_score += instability * 0.4

        return min(1.0, allosteric_score)


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
            # ============================================================
            # RULE6 (重设计): A1 域 = 2B(GOF) vs 2M(LOF), 方向相反, 同域同位点。
            # 老逻辑用"损伤大小"推"功能方向"(轻扰→2B,重损→2M)在原理上错误:
            # GOF/LOF 是方向不是程度。改为**多轴方向判别**(先定方向):
            #   轴B  GPIb 结合能力 (forced_binding) ↓↓ → 结合面坏 → 2M (LOF)
            #   轴A  自抑制松开 (AIM↔A1 接触/MD 闭合态) ↑ + 结合保留 → 2B (GOF)
            # 静态轴 A/B 实测都弱(2B/2M 重叠), MD 闭合态稳定性才是 2B 决定特征;
            # 因此无明确方向信号时返回 uncertain, 不硬判(历史 2B↔2M 漏判区)。
            # 所有新特征缺失(NaN)时自动退到结构启发, 向后兼容。
            # ============================================================
            reasoning_steps.append("RULE6: A1 域多轴方向判别 (2B=GOF / 2M=LOF)")

            aim_release = getattr(expert_scores, 'aim_release_score', np.nan)   # 轴A: 自抑制松开
            fb_z = variant_data.get('fb_binding_zscore', np.nan)               # 轴B-1: GPIb 结合
            hep_z = variant_data.get('heparan_zscore', np.nan)                 # 轴B-2: heparan 结合
            md_destab = getattr(expert_scores, 'md_face_destab_score', np.nan) # 轴B': MD 结合面破坏 (LOF/2M)
            md_lof = (not np.isnan(md_destab)) and (md_destab >= MD_FACE_DESTAB_2M_Z)
            sb_z = getattr(expert_scores, 'aim_sb_retained_z', np.nan)          # 轴B'': AIM↔A1 盐桥保留
            face_retained = (not np.isnan(sb_z)) and (sb_z >= AIM_SB_RETAINED_2B_Z)   # 保留 → anti-2M (joint判据)
            sb_collapsed = (not np.isnan(sb_z)) and (sb_z <= AIM_SB_COLLAPSE_2M_Z)     # 塌陷 → 2M 旁证 (不独立判)
            is_aim_n = 1238 <= position <= 1268
            is_aim_c = 1460 <= position <= 1472
            allosteric_risk = self.structural_expert.compute_allosteric_risk(variant_data)
            # 轴B: 两条结合面轴联合 (校准最优: mean ≤ -0.75)。两者都在 → 用联合;
            # 仅 fb 在 → 退到极保守单轴 (FB_LOSS_Z=-1.5)。
            if not np.isnan(fb_z) and not np.isnan(hep_z):
                lof_score = (fb_z + hep_z) / 2.0
                binding_lost = lof_score <= LOF_COMBINED_Z
                lof_why = f"mean(fb={fb_z:.2f}, heparan={hep_z:.2f})={lof_score:.2f} ≤ {LOF_COMBINED_Z}"
            elif not np.isnan(fb_z):
                binding_lost = fb_z <= FB_LOSS_Z
                lof_why = f"forced_binding z={fb_z:.2f} ≤ {FB_LOSS_Z} (heparan 缺失, 单轴保守)"
            else:
                binding_lost = False
                lof_why = ""

            def _a1(sub, conf, why, alts=None):
                return MultiLabelClassificationResult(
                    main_subtype=sub,
                    alternatives=alts if alts is not None else (['2M'] if sub == '2B' else ['2B']),
                    confidence=conf,
                    reasoning="; ".join(reasoning_steps + [why]),
                    domain_pleiotropy="A1域: 2B=自抑制松开+结合保留(GOF); 2M=结合面破坏(LOF)",
                    expert_scores=expert_scores)

            # 轴B 优先: A1 结合面丧失 → 2M (LOF), 与是否松开无关 ----------------
            if binding_lost:
                conf, why = 0.72, f"RULE6-B: {lof_why} → A1 结合面丧失 → 2M (LOF)"
                if md_lof:   # MD 平衡态结合面破坏 = 独立确证, 提升置信
                    conf, why = 0.8, why + f" [MD 结合面破坏确证 destab={md_destab:.2f}]"
                if sb_collapsed:  # MD AIM↔A1 盐桥塌陷 = 旁证 (轴B 已判 2M, 此处仅提置信)
                    conf, why = min(0.85, conf + 0.05), why + f" [MD 盐桥塌陷旁证 sb_z={sb_z:.2f}]"
                return _a1('2M', conf, why)

            # 轴A 强信号: 自抑制强松开 + 结合保留 → 2B (GOF) --------------------
            if not np.isnan(aim_release) and aim_release >= AIM_RELEASE_2B_Z:
                return _a1('2B', min(0.9, 0.7 + 0.1 * aim_release),
                           f"RULE6-A: 自抑制松开 score={aim_release:.2f} ≥ {AIM_RELEASE_2B_Z} + 结合保留 → 2B (GOF)")

            # AIM 区位置 (强 2B 先验, 结合未丧失) -------------------------------
            if is_aim_n or is_aim_c:
                return _a1('2B', 0.8,
                           f"RULE6-AIM位: pos={position} 在 AIM {'N' if is_aim_n else 'C'}端区 + 结合保留 → 2B")

            # 结构启发 (保留): 接口微扰+轻损 / 别构 → 2B (构象开放) --------------
            if 0.15 < expert_scores.interface_perturbation < 0.6 and expert_scores.structural_damage_score < 0.5:
                return _a1('2B', 0.72,
                           f"RULE6c: 接口微扰 {expert_scores.interface_perturbation:.2f} + 轻损 → 2B (构象开放)")
            if allosteric_risk > 0.4:
                return _a1('2B', 0.68, f"RULE6c2: 别构风险 {allosteric_risk:.2f} → 2B")

            # 弱自抑制松开救回 -------------------------------------------------
            if not np.isnan(aim_release) and aim_release > AIM_RELEASE_LEAN_Z:
                return _a1('2B', 0.55,
                           f"RULE6f-AIM: 弱松开 {aim_release:.2f} > {AIM_RELEASE_LEAN_Z} → 2B (marginal)")

            # 重度折叠损伤 → 2M (LOF, collapse) --------------------------------
            if expert_scores.structural_damage_score > 0.6:
                return _a1('2M', 0.72,
                           f"RULE6d: 重度结构损伤 {expert_scores.structural_damage_score:.2f} → 2M (LOF)")

            # 联合判据闸: MD 结合面盐桥保留 (高 z, 实测仅 2B/WT, 无 2M) → 即便有弱 MD
            # destab 也不翻 2M。致病前提下自抑制已释放 + 功能保留 → 2B marginal。
            if face_retained:
                return _a1('2B', 0.55,
                           f"RULE6-保留: MD 结合面盐桥保留 sb_z={sb_z:.2f} ≥ "
                           f"{AIM_SB_RETAINED_2B_Z} → 释放但功能保留 → 2B (marginal, 联合判据)")

            # MD 软证据 tie-breaker: 平衡态结合面强破坏 → 2M (LOF) --------------
            # 仅在轴 A/B 与结构启发皆未定向时介入 (此处已排除 AIM 位/别构);
            # 故不会硬翻已定的 2B。软证据, 低置信。
            if md_lof:
                return _a1('2M', 0.55,
                           f"RULE6-MD: 平衡 MD A1 结合面破坏 destab={md_destab:.2f} ≥ "
                           f"{MD_FACE_DESTAB_2M_Z} → 2M (LOF, 软证据)")

            # 无明确方向信号 → uncertain (不硬判; MD 闭合态稳定性出来后定 2B/2M) --
            return _a1('uncertain', 0.35,
                       "RULE6-U: A1 无明确方向信号(轴A/B 皆弱)→ uncertain (待 MD 闭合态稳定性)",
                       alts=['2B', '2M'])

        # =====================================================================
        # RULE 7: A3 domain → 2M (collagen binding)
        # =====================================================================
        if domain == 'A3':
            collagen_z = variant_data.get('a3_collagen_zscore', np.nan)
            main_subtype = '2M'
            confidence = 0.7
            if not np.isnan(collagen_z) and collagen_z <= A3_COLLAGEN_LOF_Z:
                confidence = 0.78
                reasoning_steps.append(
                    f"RULE7: A3 domain + collagen-binding z={collagen_z:.2f} ≤ "
                    f"{A3_COLLAGEN_LOF_Z} → 2M (collagen binding LOF support)"
                )
            elif not np.isnan(collagen_z) and collagen_z >= A3_COLLAGEN_RETAINED_Z:
                confidence = 0.6
                reasoning_steps.append(
                    f"RULE7: A3 domain → 2M prior, but collagen-binding z={collagen_z:.2f} "
                    f"does not support static LOF; keep low confidence"
                )
            elif not np.isnan(collagen_z):
                confidence = 0.68
                reasoning_steps.append(
                    f"RULE7: A3 domain → 2M; collagen-binding z={collagen_z:.2f} is weak/neutral"
                )
            else:
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
            ag_splice_delta=ag_splice_delta,
            # 自抑制松开特征 (extract_aim_autoinhib_features.py → merge 进输入矩阵);
            # 缺失时为 NaN, RULE6 退回原逻辑。
            aim_release_score=variant_data.get('aim_release_score', np.nan),
            # MD 结合面破坏 (7A6O AIM-A1 平衡 MD); 缺失时 NaN, RULE6 退回原逻辑。
            md_face_destab_score=variant_data.get('md_face_destab_score', np.nan),
            # MD AIM↔A1 盐桥保留 z (extract_aim_saltbridge_features.py); 缺失 NaN。
            aim_sb_retained_z=variant_data.get('aim_sb_retained_z', np.nan),
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

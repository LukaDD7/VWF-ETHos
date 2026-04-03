#!/usr/bin/env python3
"""
VWF Type-2 Competitive Multi-Label Classifier
处理域多功能性的增强分类器

关键特性：
1. 竞争式决策 - 一个域可能对应多个分型时的解决机制
2. 多标签输出 - 主分型 + 次要分型 + 备选列表
3. 置信度校准 - 基于证据强度和竞争解决
4. 可解释推理 - 完整的决策路径记录

作者: Claude Code
日期: 2026-04-02
"""

import sys
sys.path.insert(0, '/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/domain_analysis')

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
from collections import defaultdict
import json

from vwf_type2_literature_based_classifier import (
    VWFType2Classifier,
    VWFVariant,
    VWF_DOMAIN_ARCHITECTURE,
    FUNCTIONAL_REGIONS,
    KNOWN_TYPE2_HOTSPOTS
)


@dataclass
class MultiLabelClassificationResult:
    """多标签分类结果"""
    variant_id: str
    position: int
    domain: str

    # 主分型（最主要）
    primary_type: str
    primary_confidence: float

    # 次要分型（可能共存）
    secondary_type: Optional[str] = None
    secondary_confidence: Optional[float] = None

    # 备选分型（按可能性排序）
    alternatives: List[Tuple[str, float]] = field(default_factory=list)

    # 推理详情
    reasoning: str = ""
    decision_path: List[str] = field(default_factory=list)

    # 竞争解决详情
    competition_resolved: bool = False
    competing_types: List[str] = field(default_factory=list)
    resolution_strategy: str = ""

    # 可靠性评估
    reliability: str = ""  # "high", "medium", "low"
    recommended_tests: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict:
        """转换为字典格式"""
        return {
            "variant_id": self.variant_id,
            "position": self.position,
            "domain": self.domain,
            "classification": {
                "primary": {
                    "type": self.primary_type,
                    "confidence": round(self.primary_confidence, 3)
                },
                "secondary": {
                    "type": self.secondary_type,
                    "confidence": round(self.secondary_confidence, 3)
                } if self.secondary_type else None,
                "alternatives": [(t, round(c, 3)) for t, c in self.alternatives]
            },
            "reasoning": self.reasoning,
            "decision_path": self.decision_path,
            "competition": {
                "resolved": self.competition_resolved,
                "competing_types": self.competing_types,
                "strategy": self.resolution_strategy
            },
            "reliability": self.reliability,
            "recommended_tests": self.recommended_tests
        }


class CompetitiveVWFClassifier(VWFType2Classifier):
    """
    竞争式VWF Type-2分类器
    继承并扩展基础分类器，添加多标签和竞争解决功能
    """

    # 域多功能性配置
    DOMAIN_PLEIOTROPY = {
        "A1": {
            "types": ["2B", "2M"],
            "resolution": "AIM_based",  # 基于AIM参与程度
            "description": "A1 domain: Type 2B (gain-of-function) vs Type 2M (loss-of-function)"
        },
        "A3": {
            "types": ["2M", "2A"],
            "resolution": "position_based",  # 基于位置
            "description": "A3 domain: Primarily Type 2M (collagen), some may affect proteolysis"
        },
        "C1": {
            "types": ["2A", "2M"],
            "resolution": "evidence_based",  # 证据不足时默认2A
            "description": "C1 domain: Multimerization (2A) vs collagen binding (2M)"
        },
        "C2": {
            "types": ["2A", "2M"],
            "resolution": "evidence_based",
            "description": "C2 domain: Multimerization (2A) vs collagen binding (2M)"
        },
        "D4": {
            "types": ["2A", "1"],
            "resolution": "multimer_based",  # 需要多聚体数据
            "description": "D4 domain: Multimerization defect (2A) vs secretion defect (Type 1)"
        }
    }

    # 置信度阈值
    CONFIDENCE_THRESHOLDS = {
        "high": 0.85,
        "medium": 0.60,
        "low": 0.40
    }

    def __init__(self):
        super().__init__()
        self.base_classifier = VWFType2Classifier()

    def competitive_classify(self, variant: VWFVariant,
                            multimer_data: Optional[Dict] = None) -> MultiLabelClassificationResult:
        """
        竞争式分类主函数

        Args:
            variant: VWF变异
            multimer_data: 可选的多聚体分析数据

        Returns:
            MultiLabelClassificationResult: 多标签分类结果
        """
        # 1. 获取基础分类得分
        base_scores = self._calculate_base_scores(variant)

        # 2. 确定域和可能的竞争
        domain = self.get_domain_for_position(variant.position)
        variant.domain = domain

        # 3. 检查是否需要竞争解决
        if domain in self.DOMAIN_PLEIOTROPY:
            # 需要竞争解决
            result = self._resolve_competition(
                variant, domain, base_scores, multimer_data
            )
        else:
            # 单一分型，直接输出
            result = self._single_type_classification(
                variant, domain, base_scores
            )

        # 4. 评估可靠性
        result.reliability = self._assess_reliability(result)
        result.recommended_tests = self._get_recommended_tests(result, domain)

        return result

    def _calculate_base_scores(self, variant: VWFVariant) -> Dict[str, float]:
        """计算所有Type-2分型的基础得分"""
        scores = {}

        # 计算特征
        features = self.calculate_structure_based_features(variant)

        # Type 2A scoring
        if variant.domain == "A2" or features.get("adamts13_sensitivity", 0) > 0:
            scores["2A"] = features.get("adamts13_sensitivity", 0) * 2.0
            if features.get("cleavage_site_proximity", 100) < 10:
                scores["2A"] = scores.get("2A", 0) + 1.0

        # Type 2B scoring
        if variant.domain == "A1":
            aim_score = features.get("aim_disruption_potential", 0)
            gpib_score = features.get("gpib_binding_impact", 0)

            if aim_score > 0.5 and gpib_score > 0.5:
                scores["2B"] = 3.0
            elif aim_score > 0.5:
                scores["2B"] = 1.5

        # Type 2M scoring
        if variant.domain == "A3":
            scores["2M"] = features.get("collagen_binding_impact", 0) * 2.0

        if variant.domain == "A1":
            gpib_score = features.get("gpib_binding_impact", 0)
            aim_score = features.get("aim_disruption_potential", 0)
            if gpib_score > 0.5 and aim_score < 0.3:
                scores["2M"] = scores.get("2M", 0) + 1.5

        # Type 2N scoring
        if variant.domain in ["D_prime", "D3"]:
            scores["2N"] = features.get("fviii_binding_impact", 0) * 2.0

        # Propeptide/D1-D2/D4 → Type 2A
        if variant.domain in ["propeptide_D1", "propeptide_D2", "D4", "CK"]:
            base_score = 1.5 if variant.domain == "D4" else 2.0
            if variant.domain == "propeptide_D1":
                base_score += 0.5
            scores["2A"] = base_score

        # C-terminal domains
        if variant.domain in ["C1", "C2"]:
            scores["2M"] = 0.8
            scores["2A"] = 0.5
        elif variant.domain in ["C3", "C4", "C5", "C6"]:
            scores["2A"] = 1.0
        elif variant.domain == "CK":
            scores["2A"] = 2.5

        return scores

    def _resolve_competition(self, variant: VWFVariant, domain: str,
                            scores: Dict[str, float],
                            multimer_data: Optional[Dict]) -> MultiLabelClassificationResult:
        """
        解决域内竞争
        """
        pleiotropy_info = self.DOMAIN_PLEIOTROPY[domain]
        competing_types = pleiotropy_info["types"]
        strategy = pleiotropy_info["resolution"]

        # 根据策略选择解决方法
        if strategy == "AIM_based":
            return self._resolve_A1_competition(variant, scores, competing_types)
        elif strategy == "position_based":
            return self._resolve_A3_competition(variant, scores, competing_types)
        elif strategy == "evidence_based":
            return self._resolve_C_domain_competition(variant, scores, competing_types)
        elif strategy == "multimer_based":
            return self._resolve_D4_competition(variant, scores, competing_types, multimer_data)
        else:
            # 默认：选择得分最高的
            return self._default_resolution(variant, domain, scores, competing_types)

    def _resolve_A1_competition(self, variant: VWFVariant,
                               scores: Dict[str, float],
                               competing_types: List[str]) -> MultiLabelClassificationResult:
        """
        解决A1域竞争：Type 2B vs Type 2M
        关键：是否破坏AIM（自抑制模块）
        """
        pos = variant.position

        # AIM区域检查
        in_AIM_N = 1238 <= pos <= 1268
        in_AIM_C = 1460 <= pos <= 1472
        in_AIM = in_AIM_N or in_AIM_C

        # GPIbα结合界面
        in_gpib_interface = 1296 <= pos <= 1350

        # 电荷改变检查
        charged_aas = {"R", "K", "D", "E", "H"}
        ref_charge = variant.ref_aa in charged_aas
        alt_charge = variant.alt_aa in charged_aas
        charge_change = (ref_charge != alt_charge)

        # 构建结果
        result = MultiLabelClassificationResult(
            variant_id=variant.variant_id,
            position=pos,
            domain="A1",
            primary_type="",  # 稍后填充
            primary_confidence=0.0,  # 稍后填充
            competing_types=competing_types,
            competition_resolved=True
        )

        if in_AIM:
            # 明确Type 2B
            result.primary_type = "2B"
            result.primary_confidence = 0.95
            result.reasoning = (
                f"Position {pos} is within the Autoinhibitory Module (AIM) "
                f"({ 'N-terminal' if in_AIM_N else 'C-terminal' } region). "
                f"AIM disruption causes gain-of-function → Type 2B."
            )
            result.decision_path = [
                f"1. Variant {variant.variant_id} located in A1 domain",
                f"2. Position {pos} within AIM region (1238-1268 or 1460-1472)",
                f"3. AIM disruption causes spontaneous GPIbα binding",
                f"4. Classification: Type 2B with high confidence"
            ]
            result.resolution_strategy = "AIM_involvement"

        elif in_gpib_interface and not in_AIM:
            # GPIb界面但不在AIM
            if charge_change:
                # 电荷改变 → 可能是2B
                result.primary_type = "2B"
                result.primary_confidence = 0.75
                result.secondary_type = "2M"
                result.secondary_confidence = 0.35
                result.reasoning = (
                    f"Position {pos} in GPIbα interface but not AIM. "
                    f"Charge change favors gain-of-function (2B), but "
                    f"loss-of-function (2M) cannot be excluded."
                )
                result.alternatives = [("2M", 0.35), ("unclassified", 0.05)]
                result.decision_path = [
                    f"1. Variant in GPIbα interface (1296-1350)",
                    f"2. No AIM involvement → excludes classic 2B",
                    f"3. Charge change detected → possible gain-of-function",
                    f"4. Primary: 2B (75%), Secondary: 2M (35%)"
                ]
                result.resolution_strategy = "charge_effect"
            else:
                # 无电荷改变 → Type 2M
                result.primary_type = "2M"
                result.primary_confidence = 0.85
                result.reasoning = (
                    f"Position {pos} in GPIbα interface without AIM disruption. "
                    f"No charge change → loss-of-function (Type 2M)."
                )
                result.decision_path = [
                    f"1. Variant in GPIbα interface",
                    f"2. No AIM involvement",
                    f"3. No charge change → reduced binding affinity",
                    f"4. Classification: Type 2M"
                ]
                result.resolution_strategy = "loss_of_function"
        else:
            # A1域其他位置
            result.primary_type = "2M"
            result.primary_confidence = 0.60
            result.reasoning = (
                f"Position {pos} in A1 domain but outside GPIb interface. "
                f"Likely affects overall domain stability → Type 2M."
            )
            result.alternatives = [("unclassified", 0.25), ("2B", 0.15)]
            result.resolution_strategy = "default_A1"

        return result

    def _resolve_A3_competition(self, variant: VWFVariant,
                               scores: Dict[str, float],
                               competing_types: List[str]) -> MultiLabelClassificationResult:
        """解决A3域竞争"""
        # A3主要是Type 2M，少数情况可能Type 2A

        result = MultiLabelClassificationResult(
            variant_id=variant.variant_id,
            position=variant.position,
            domain="A3",
            primary_type="2M",
            primary_confidence=0.85,
            competing_types=competing_types,
            competition_resolved=True,
            reasoning="A3 domain primarily associated with collagen binding (Type 2M).",
            decision_path=[
                f"1. Variant in A3 domain (collagen binding)",
                f"2. Primary function: collagen binding",
                f"3. Classification: Type 2M"
            ],
            resolution_strategy="primary_function"
        )

        # 如果靠近A2边界，可能有2A风险
        if variant.position < 1700:
            result.secondary_type = "2A"
            result.secondary_confidence = 0.20
            result.alternatives = [("2A", 0.20)]
            result.reasoning += " Proximity to A2 suggests possible ADAMTS13 sensitivity effect."

        return result

    def _resolve_C_domain_competition(self, variant: VWFVariant,
                                     scores: Dict[str, float],
                                     competing_types: List[str]) -> MultiLabelClassificationResult:
        """解决C1-C2域竞争"""

        # 计算氨基酸性质改变
        aa_size = {
            'A': 88.6, 'C': 108.5, 'D': 111.1, 'E': 138.4, 'F': 189.9,
            'G': 60.1, 'H': 153.2, 'I': 166.7, 'K': 168.6, 'L': 166.7,
            'M': 162.9, 'N': 114.1, 'P': 112.7, 'Q': 143.8, 'R': 173.4,
            'S': 89.0, 'T': 116.1, 'V': 140.0, 'W': 227.8, 'Y': 193.6
        }

        size_change = abs(aa_size.get(variant.alt_aa, 0) - aa_size.get(variant.ref_aa, 0))

        result = MultiLabelClassificationResult(
            variant_id=variant.variant_id,
            position=variant.position,
            domain=variant.domain,
            primary_type="",  # 稍后填充
            primary_confidence=0.0,  # 稍后填充
            competing_types=competing_types,
            competition_resolved=True
        )

        # 大残基改变 → 更可能影响多聚化 (2A)
        if size_change > 50:
            result.primary_type = "2A"
            result.primary_confidence = 0.70
            result.secondary_type = "2M"
            result.secondary_confidence = 0.30
            result.reasoning = (
                f"Large side chain change ({size_change:.1f} Å³) in {variant.domain} "
                f"suggests structural disruption → multimerization defect (2A). "
                f"Collagen binding may also be affected (2M)."
            )
            result.alternatives = [("2M", 0.30), ("1", 0.10)]
            result.resolution_strategy = "size_effect"
        else:
            # 默认2A（C域主要功能是多聚化支持）
            result.primary_type = "2A"
            result.primary_confidence = 0.55
            result.secondary_type = "2M"
            result.secondary_confidence = 0.25
            result.reasoning = (
                f"{variant.domain} domain primarily supports multimerization. "
                f"Default classification: Type 2A with secondary 2M possibility."
            )
            result.alternatives = [("2M", 0.25), ("unclassified", 0.20)]
            result.resolution_strategy = "default_C_domain"

        result.decision_path = [
            f"1. Variant in {variant.domain} domain",
            f"2. C-domains support multimerization (2A) and collagen binding (2M)",
            f"3. Size change: {size_change:.1f} Å³",
            f"4. Primary: {result.primary_type}, Secondary: {result.secondary_type}"
        ]

        return result

    def _resolve_D4_competition(self, variant: VWFVariant,
                               scores: Dict[str, float],
                               competing_types: List[str],
                               multimer_data: Optional[Dict]) -> MultiLabelClassificationResult:
        """解决D4域竞争：需要多聚体数据"""

        result = MultiLabelClassificationResult(
            variant_id=variant.variant_id,
            position=variant.position,
            domain="D4",
            primary_type="",  # 稍后填充
            primary_confidence=0.0,  # 稍后填充
            competing_types=competing_types,
            competition_resolved=True
        )

        if multimer_data:
            pattern = multimer_data.get("pattern", "unknown")

            if pattern == "loss_of_HMW":
                result.primary_type = "2A"
                result.primary_confidence = 0.90
                result.reasoning = "Loss of HMW multimers indicates multimerization defect (Type 2A)."
                result.resolution_strategy = "multimer_evidence"
            elif pattern == "normal":
                result.primary_type = "1"
                result.primary_confidence = 0.75
                result.reasoning = "Normal multimers with reduced VWF:Ag suggests Type 1 (secretion defect)."
                result.resolution_strategy = "multimer_evidence"
            else:
                result.primary_type = "2A"
                result.primary_confidence = 0.65
                result.reasoning = "Ambiguous multimer pattern; defaulting to Type 2A."
                result.alternatives = [("1", 0.25), ("unclassified", 0.10)]
                result.resolution_strategy = "ambiguous_multimer"
        else:
            # 无多聚体数据
            result.primary_type = "2A"
            result.primary_confidence = 0.60
            result.secondary_type = "1"
            result.secondary_confidence = 0.30
            result.reasoning = (
                "No multimer data available. D4 mutations most commonly "
                "cause Type 2A, but Type 1 possible if multimers normal."
            )
            result.alternatives = [("1", 0.30), ("unclassified", 0.10)]
            result.resolution_strategy = "no_multimer_data"

        result.decision_path = [
            f"1. Variant in D4 domain",
            f"2. D4 mutations can cause 2A (multimerization) or 1 (secretion)",
            f"3. Multimer data: {multimer_data.get('pattern', 'N/A') if multimer_data else 'Not available'}",
            f"4. Classification: {result.primary_type}"
        ]

        return result

    def _single_type_classification(self, variant: VWFVariant, domain: str,
                                   scores: Dict[str, float]) -> MultiLabelClassificationResult:
        """单一分型分类（无竞争）"""

        # 选择最高得分的类型
        if not scores:
            primary_type = "unclassified"
            confidence = 0.0
        else:
            primary_type = max(scores, key=scores.get)
            max_score = scores[primary_type]
            confidence = min(max_score / 3.0, 1.0)

        return MultiLabelClassificationResult(
            variant_id=variant.variant_id,
            position=variant.position,
            domain=domain,
            primary_type=primary_type,
            primary_confidence=confidence,
            competing_types=[],
            competition_resolved=False,
            reasoning=f"Single classification for {domain} domain.",
            decision_path=[
                f"1. Variant in {domain} domain",
                f"2. No competing Type-2 associations",
                f"3. Classification: Type {primary_type}"
            ]
        )

    def _default_resolution(self, variant: VWFVariant, domain: str,
                           scores: Dict[str, float],
                           competing_types: List[str]) -> MultiLabelClassificationResult:
        """默认竞争解决：选择得分最高的"""

        primary_type = max(scores, key=scores.get) if scores else "unclassified"

        # 构建备选列表
        alternatives = []
        for t, s in sorted(scores.items(), key=lambda x: x[1], reverse=True):
            if t != primary_type:
                alternatives.append((t, min(s / 3.0, 1.0)))

        return MultiLabelClassificationResult(
            variant_id=variant.variant_id,
            position=variant.position,
            domain=domain,
            primary_type=primary_type,
            primary_confidence=min(scores.get(primary_type, 0) / 3.0, 1.0) if scores else 0.0,
            alternatives=alternatives,
            competing_types=competing_types,
            competition_resolved=True,
            reasoning=f"Default resolution: highest scoring type selected.",
            resolution_strategy="default_highest_score"
        )

    def _assess_reliability(self, result: MultiLabelClassificationResult) -> str:
        """评估分类可靠性"""

        if result.primary_confidence >= self.CONFIDENCE_THRESHOLDS["high"]:
            return "high"
        elif result.primary_confidence >= self.CONFIDENCE_THRESHOLDS["medium"]:
            return "medium"
        else:
            return "low"

    def _get_recommended_tests(self, result: MultiLabelClassificationResult,
                              domain: str) -> List[str]:
        """获取推荐验证测试"""

        tests = []

        if result.reliability == "low":
            tests.append("Functional assays required")

        if domain == "A1":
            tests.append("RIPA (ristocetin-induced platelet aggregation)")
            if result.primary_type == "2B":
                tests.append("Platelet count (check for thrombocytopenia)")

        elif domain in ["propeptide_D1", "propeptide_D2", "D4", "CK"]:
            tests.append("VWF multimer analysis")
            if domain == "D4":
                tests.append("VWFpp/Ag ratio")

        elif domain == "A2":
            tests.append("VWF multimer analysis")
            tests.append("ADAMTS13 activity")

        elif domain in ["D_prime", "D3"]:
            tests.append("FVIII:C level")
            tests.append("VWF:FVIII binding assay")

        elif domain == "A3":
            tests.append("VWF:CB (collagen binding)")

        return tests

    def generate_report(self, result: MultiLabelClassificationResult) -> str:
        """生成人类可读的分类报告"""

        report = []
        report.append("=" * 80)
        report.append(f"Competitive Classification Report: {result.variant_id}")
        report.append("=" * 80)
        report.append("")

        # 基本信息
        report.append(f"Variant: {result.variant_id}")
        report.append(f"Position: {result.position}")
        report.append(f"Domain: {result.domain}")
        report.append("")

        # 分类结果
        report.append("-" * 80)
        report.append("CLASSIFICATION")
        report.append("-" * 80)
        report.append(f"Primary Type: Type {result.primary_type}")
        report.append(f"Confidence: {result.primary_confidence:.1%}")

        if result.secondary_type:
            report.append(f"Secondary Type: Type {result.secondary_type}")
            report.append(f"Secondary Confidence: {result.secondary_confidence:.1%}")

        if result.alternatives:
            report.append("\nAlternative Classifications:")
            for alt_type, alt_conf in result.alternatives:
                report.append(f"  - Type {alt_type}: {alt_conf:.1%}")

        report.append("")

        # 可靠性
        report.append("-" * 80)
        report.append("RELIABILITY ASSESSMENT")
        report.append("-" * 80)
        report.append(f"Reliability: {result.reliability.upper()}")
        report.append("")

        # 推理
        report.append("-" * 80)
        report.append("REASONING")
        report.append("-" * 80)
        report.append(result.reasoning)
        report.append("")

        # 决策路径
        report.append("-" * 80)
        report.append("DECISION PATH")
        report.append("-" * 80)
        for step in result.decision_path:
            report.append(step)
        report.append("")

        # 竞争解决
        if result.competition_resolved:
            report.append("-" * 80)
            report.append("COMPETITION RESOLUTION")
            report.append("-" * 80)
            report.append(f"Competing Types: {', '.join(result.competing_types)}")
            report.append(f"Resolution Strategy: {result.resolution_strategy}")
            report.append("")

        # 推荐测试
        if result.recommended_tests:
            report.append("-" * 80)
            report.append("RECOMMENDED VERIFICATION TESTS")
            report.append("-" * 80)
            for test in result.recommended_tests:
                report.append(f"  • {test}")
            report.append("")

        report.append("=" * 80)

        return "\n".join(report)


def main():
    """测试竞争式分类器"""

    print("=" * 80)
    print("Competitive VWF Type-2 Classifier - Test")
    print("=" * 80)
    print()

    classifier = CompetitiveVWFClassifier()

    # 测试用例：涵盖各种竞争场景
    test_cases = [
        # A1域：明确Type 2B（AIM区域）
        ("R1306W", 1306, "R", "W"),
        # A1域：明确Type 2M（GPIb界面，无AIM）
        ("I1327T", 1327, "I", "T"),
        # A1域：模糊（GPIb界面，电荷改变）
        ("G1324E", 1324, "G", "E"),
        # D4域：无多聚体数据
        ("P1888L", 1888, "P", "L"),
        # C1域：可能2A或2M
        ("V2465M", 2465, "V", "M"),
    ]

    for var_id, pos, ref, alt in test_cases:
        variant = VWFVariant(
            variant_id=var_id,
            protein_change=f"p.{ref}{pos}{alt}",
            position=pos,
            ref_aa=ref,
            alt_aa=alt,
            acmg_classification="VUS"
        )

        result = classifier.competitive_classify(variant)

        print(f"\n{'='*80}")
        print(classifier.generate_report(result))
        print()


if __name__ == "__main__":
    main()

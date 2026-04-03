#!/usr/bin/env python3
"""
集成示例：使用竞争式分类器替换原有Pipeline

演示如何在主Pipeline中集成新的竞争式多标签分类器
"""

import sys
sys.path.insert(0, '/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/domain_analysis')

import pandas as pd
from pathlib import Path
from typing import Dict, Optional

from vwf_competitive_classifier import (
    CompetitiveVWFClassifier,
    MultiLabelClassificationResult
)
from vwf_type2_literature_based_classifier import VWFVariant


class EnhancedVWFAnalysisPipeline:
    """
    增强版VWF分析Pipeline
    使用竞争式多标签分类器替代原有分类器
    """

    def __init__(self, wt_structure_path: Optional[str] = None):
        # 使用新的竞争式分类器
        self.classifier = CompetitiveVWFClassifier()
        self.results = []

    def process_variant_competitive(self, row: pd.Series) -> Dict:
        """使用竞争式分类器处理单个变异"""

        # 解析变异信息
        variant_id = row.get('Name', '')
        position = int(row.get('AA_Position', 0))
        ref_aa = row.get('REF_AA', '')
        alt_aa = row.get('ALT_AA', '')

        # 创建变异对象
        variant = VWFVariant(
            variant_id=variant_id,
            protein_change=f"p.{ref_aa}{position}{alt_aa}",
            position=position,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
            acmg_classification=row.get('INFO_acmg_classification_base', 'VUS'),
            type2_subtype=row.get('Type2_Subtype', None)
        )

        # 使用竞争式分类
        # 注意：可以传入多聚体数据（如果有）
        multimer_data = None  # 可以从row中获取如果有

        result = self.classifier.competitive_classify(variant, multimer_data)

        # 转换为字典
        result_dict = result.to_dict()

        # 添加额外信息
        result_dict['known_type'] = row.get('Type2_Subtype', 'N/A')
        result_dict['match'] = self._check_match(result, row)

        return result_dict

    def _check_match(self, result: MultiLabelClassificationResult, row: pd.Series) -> bool:
        """检查预测是否与已知类型匹配"""
        known = row.get('Type2_Subtype', '').upper().replace('TYPE', '').strip()
        predicted = result.primary_type.upper()

        # 考虑次要类型
        if result.secondary_type:
            return predicted == known or result.secondary_type.upper() == known
        return predicted == known

    def run_competitive_analysis(self, variants_csv: str) -> pd.DataFrame:
        """运行竞争式分析"""

        # 加载数据
        df = pd.read_csv(variants_csv)
        print(f"Loaded {len(df)} variants from {variants_csv}")

        # 筛选Type-2变异
        df = df[df['Type2_Subtype'] != 'WT_Control']

        print(f"\nProcessing {len(df)} Type-2 variants with competitive classifier...")
        print("=" * 80)

        results = []
        for idx, row in df.iterrows():
            if idx % 10 == 0:
                print(f"  Processed {idx}/{len(df)} variants...")

            result = self.process_variant_competitive(row)
            results.append(result)

        # 转换为DataFrame
        results_df = pd.json_normalize(results)

        # 生成报告
        self._generate_competitive_report(results_df)

        return results_df

    def _generate_competitive_report(self, results_df: pd.DataFrame):
        """生成竞争式分类报告"""

        report = []
        report.append("=" * 80)
        report.append("Competitive Classification Analysis Report")
        report.append("=" * 80)
        report.append("")

        # 统计
        total = len(results_df)
        report.append(f"Total variants analyzed: {total}")
        report.append("")

        # 按主分型统计
        report.append("Primary Classifications:")
        primary_counts = results_df['classification.primary.type'].value_counts()
        for ptype, count in primary_counts.items():
            report.append(f"  Type {ptype}: {count} ({count/total*100:.1f}%)")
        report.append("")

        # 有竞争解决的变异
        resolved = results_df[results_df['competition.resolved'] == True]
        report.append(f"Variants with competition resolution: {len(resolved)} ({len(resolved)/total*100:.1f}%)")
        report.append("")

        # 按域统计竞争类型
        report.append("Domain Pleiotropy Summary:")
        for domain in ['A1', 'A3', 'C1', 'C2', 'D4']:
            domain_variants = resolved[resolved['domain'] == domain]
            if len(domain_variants) > 0:
                competing = domain_variants['competition.competing_types'].iloc[0]
                report.append(f"  {domain}: Competing types {competing}")
        report.append("")

        # 可靠性分布
        report.append("Reliability Distribution:")
        reliability_counts = results_df['reliability'].value_counts()
        for rel, count in reliability_counts.items():
            report.append(f"  {rel.upper()}: {count} ({count/total*100:.1f}%)")
        report.append("")

        # 多标签输出示例
        with_secondary = results_df[results_df['classification.secondary.type'].notna()]
        report.append(f"Variants with secondary classification: {len(with_secondary)}")
        report.append("")

        if len(with_secondary) > 0:
            report.append("Examples of multi-label classification:")
            for _, row in with_secondary.head(5).iterrows():
                report.append(f"  {row['variant_id']}: "
                            f"Primary={row['classification.primary.type']} "
                            f"({row['classification.primary.confidence']:.0%}), "
                            f"Secondary={row['classification.secondary.type']} "
                            f"({row['classification.secondary.confidence']:.0%})")
        report.append("")

        # 与已知类型对比（如果有）
        if 'match' in results_df.columns:
            matches = results_df['match'].sum()
            report.append(f"Accuracy (vs known types): {matches}/{total} ({matches/total*100:.1f}%)")
            report.append("")

        report.append("=" * 80)

        # 保存报告
        report_text = "\n".join(report)
        print(report_text)

        # 保存到文件
        output_path = Path("competitive_classification_report.txt")
        with open(output_path, 'w') as f:
            f.write(report_text)
        print(f"\nReport saved to: {output_path}")


def demo_comparison():
    """演示：对比基础分类器和竞争式分类器"""

    from vwf_type2_literature_based_classifier import VWFType2Classifier

    print("=" * 80)
    print("Comparison: Basic vs Competitive Classifier")
    print("=" * 80)
    print()

    # 测试用例：涵盖竞争场景
    test_cases = [
        # (变异ID, 位置, ref, alt, 描述)
        ("R1306W", 1306, "R", "W", "A1-AIM region (classic 2B)"),
        ("V1316M", 1316, "V", "M", "A1-GPIb interface (ambiguous)"),
        ("I1628T", 1628, "I", "T", "A2-cleavage site (clear 2A)"),
        ("P1888L", 1888, "P", "L", "D4-multimerization/secretion"),
        ("V2465M", 2465, "V", "M", "C3-uncertain"),
    ]

    basic_classifier = VWFType2Classifier()
    competitive_classifier = CompetitiveVWFClassifier()

    print(f"{'Variant':<12} {'Description':<35} {'Basic':<15} {'Competitive':<30}")
    print("-" * 100)

    for var_id, pos, ref, alt, desc in test_cases:
        # 基础分类
        variant = VWFVariant(
            variant_id=var_id,
            protein_change=f"p.{ref}{pos}{alt}",
            position=pos,
            ref_aa=ref,
            alt_aa=alt,
            acmg_classification="VUS"
        )
        basic_type, basic_conf = basic_classifier.predict_subtype(variant)

        # 竞争式分类
        comp_result = competitive_classifier.competitive_classify(variant)
        comp_type = comp_result.primary_type
        comp_conf = comp_result.primary_confidence
        comp_secondary = comp_result.secondary_type if comp_result.secondary_type else "-"

        basic_str = f"{basic_type} ({basic_conf:.0%})"
        comp_str = f"{comp_type} ({comp_conf:.0%})"
        if comp_result.secondary_type:
            comp_str += f" + {comp_secondary}"

        print(f"{var_id:<12} {desc:<35} {basic_str:<15} {comp_str:<30}")

    print()
    print("Key Differences:")
    print("  - Competitive classifier identifies secondary types")
    print("  - Resolves domain pleiotropy (e.g., A1: 2B vs 2M)")
    print("  - Provides alternatives when uncertain")
    print("  - Recommends specific verification tests")


def main():
    """主函数"""

    print("=" * 80)
    print("Enhanced VWF Pipeline with Competitive Classifier")
    print("=" * 80)
    print()

    # 演示对比
    demo_comparison()

    print()
    print("=" * 80)
    print("Integration Instructions:")
    print("=" * 80)
    print()
    print("1. Replace in vwf_type2_domain_pipeline.py:")
    print()
    print("   OLD:")
    print("     from vwf_type2_literature_based_classifier import VWFType2Classifier")
    print("     self.classifier = VWFType2Classifier()")
    print("     predicted, confidence = self.classifier.predict_subtype(variant)")
    print()
    print("   NEW:")
    print("     from vwf_competitive_classifier import CompetitiveVWFClassifier")
    print("     self.classifier = CompetitiveVWFClassifier()")
    print("     result = self.classifier.competitive_classify(variant)")
    print("     predicted = result.primary_type")
    print("     confidence = result.primary_confidence")
    print("     secondary = result.secondary_type  # may be None")
    print()
    print("2. Access full result details:")
    print("     result.to_dict()  # Convert to dictionary")
    print("     classifier.generate_report(result)  # Human-readable report")
    print()
    print("3. Handle multi-label output in your downstream analysis")
    print()


if __name__ == "__main__":
    main()

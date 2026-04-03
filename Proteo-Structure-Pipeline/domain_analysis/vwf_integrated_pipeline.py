#!/usr/bin/env python3
"""
Integrated Pipeline for VWF Residue-Level Feature Extraction
整合残基级特征提取的Pipeline

整合内容：
1. 残基级文献特征 (vwf_residue_feature_extractor)
2. AF3结构特征 (af3_cif_parser)
3. 统一输出格式

Author: Claude Code
Date: 2026-04-03
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import asdict

# 添加domain_analysis目录到路径
sys.path.insert(0, str(Path(__file__).parent))

from vwf_residue_feature_extractor import (
    VWFResidueFeatureExtractor,
    ResidueLevelFeatures
)
from af3_cif_parser import AF3CIFParser, AF3StructureData


class IntegratedVWFPipeline:
    """
    整合Pipeline - 结合残基级特征和AF3结构特征
    """

    def __init__(self, af3_base_dir: Optional[Path] = None):
        """
        初始化Pipeline

        Args:
            af3_base_dir: AF3结果目录 (可选)
        """
        self.residue_extractor = VWFResidueFeatureExtractor()
        self.af3_parser = AF3CIFParser(af3_base_dir) if af3_base_dir else None

        # 加载WT数据（如果可用）
        self.wt_data = None
        if af3_base_dir:
            wt_dir = Path(af3_base_dir) / "fold_vwf_wt"
            if wt_dir.exists():
                self.wt_data = self.af3_parser.parse_variant_directory(wt_dir)
                print(f"Loaded WT data from {wt_dir}")

    def process_variant(self, variant_id: str, position: int,
                       ref_aa: str, alt_aa: str,
                       af3_variant_dir: Optional[Path] = None) -> Dict:
        """
        处理单个变异 - 提取完整特征

        Args:
            variant_id: 变异ID (如 "R1306W")
            position: 氨基酸位置
            ref_aa: 参考氨基酸
            alt_aa: 突变氨基酸
            af3_variant_dir: AF3变异体目录 (可选)

        Returns:
            完整特征字典
        """
        # 1. 提取残基级文献特征
        residue_features = self.residue_extractor.extract_residue_features(
            variant_id=variant_id,
            position=position,
            ref_aa=ref_aa,
            alt_aa=alt_aa
        )

        # 2. 提取AF3结构特征（如果可用）
        structural_data = {}
        if af3_variant_dir and self.af3_parser:
            af3_data = self.af3_parser.parse_variant_directory(af3_variant_dir)
            if af3_data:
                structural_data = self._extract_structural_features(af3_data, position)

        # 3. 合并特征
        combined_features = self._merge_features(residue_features, structural_data)

        return combined_features

    def _extract_structural_features(self, af3_data: AF3StructureData,
                                    position: int) -> Dict:
        """从AF3数据提取结构特征"""
        features = {
            "af3_available": True,
            "plddt_mean": af3_data.mean_plddt,
            "plddt_min": af3_data.min_plddt,
            "plddt_max": af3_data.max_plddt,
            "plddt_at_site": af3_data.get_plddt_at_position(position),
            "global_plddt_cif": af3_data.global_plddt,
            "ptm_score": af3_data.ptm_score,
            "iptm_score": af3_data.iptm_score,
            "ranking_score": af3_data.ranking_score,
            "has_clash": af3_data.has_clash,
            "fraction_disordered": af3_data.fraction_disordered,
        }

        # 计算与WT的差异
        if self.wt_data:
            wt_plddt = self.wt_data.get_plddt_at_position(position)
            var_plddt = af3_data.get_plddt_at_position(position)

            features["wt_plddt_at_site"] = wt_plddt
            features["plddt_delta"] = var_plddt - wt_plddt if var_plddt and wt_plddt else None
            features["plddt_mean_delta"] = af3_data.mean_plddt - self.wt_data.mean_plddt

        return features

    def _merge_features(self, residue_features: ResidueLevelFeatures,
                       structural_data: Dict) -> Dict:
        """合并残基特征和结构特征"""
        # 基础特征
        merged = {
            # 变异信息
            "variant_id": residue_features.variant_id,
            "position": residue_features.position,
            "ref_aa": residue_features.ref_aa,
            "alt_aa": residue_features.alt_aa,

            # 域信息
            "domain": residue_features.domain,
            "relative_position_in_domain": residue_features.relative_position,

            # AIM特征
            "is_in_aim": residue_features.is_in_AIM,
            "aim_component": residue_features.aim_component,
            "aim_disruption_score": residue_features.aim_disruption_score,

            # GPIbα界面特征
            "is_gpib_interface": residue_features.is_in_gpib_interface,
            "gpib_site": residue_features.gpib_site,
            "gpib_key_residue": residue_features.gpib_key_residue,
            "gpib_interaction_type": residue_features.gpib_interaction_type,

            # ADAMTS13特征
            "is_scissile_bond": residue_features.is_scissile_bond,
            "is_exosite_1": residue_features.is_in_exosite_1,
            "is_exosite_2": residue_features.is_in_exosite_2,
            "is_exosite_3": residue_features.is_in_exosite_3,
            "is_group1_2a": residue_features.is_group1_2A_mutation,
            "is_group2_2a": residue_features.is_group2_2A_mutation,

            # 功能残基
            "is_calcium_coordinating": residue_features.is_calcium_coordinating,
            "is_collagen_binding": residue_features.is_collagen_binding_residue,
            "is_fviii_binding": residue_features.is_fviii_binding_residue,
            "is_rgd_motif": residue_features.is_rgd_motif,
            "is_multimerization_residue": residue_features.is_multimerization_residue,

            # VWD热点
            "is_vwd_hotspot": residue_features.is_known_vwd_hotspot,
            "vwd_hotspot_type": residue_features.vwd_hotspot_type,
            "known_vwd_mutation": residue_features.vwd_hotspot_mutation,

            # 突变性质
            "mutation_size_delta": residue_features.mutation_size_delta,
            "mutation_charge_change": residue_features.mutation_charge_change,
            "mutation_hydrophobicity_delta": residue_features.mutation_hydrophobicity_delta,
            "mutation_aromatic_change": residue_features.mutation_aromatic_change,

            # 文献引用
            "literature_pmids": residue_features.literature_pmids,
        }

        # 添加结构特征（如果可用）
        if structural_data:
            merged.update({
                "af3_available": True,
                "plddt_mean": structural_data.get("plddt_mean"),
                "plddt_at_site": structural_data.get("plddt_at_site"),
                "plddt_delta": structural_data.get("plddt_delta"),
                "ptm_score": structural_data.get("ptm_score"),
                "ranking_score": structural_data.get("ranking_score"),
                "has_clash": structural_data.get("has_clash"),
            })
        else:
            merged["af3_available"] = False

        return merged

    def process_batch(self, variants_df: pd.DataFrame,
                     af3_base_dir: Optional[Path] = None) -> pd.DataFrame:
        """
        批量处理变异

        Args:
            variants_df: DataFrame with columns [variant_id, position, ref_aa, alt_aa]
            af3_base_dir: AF3结果目录

        Returns:
            特征DataFrame
        """
        records = []

        for idx, row in variants_df.iterrows():
            variant_id = row.get('variant_id', f"{row['ref_aa']}{row['position']}{row['alt_aa']}")

            # 查找AF3目录
            af3_dir = None
            if af3_base_dir:
                expected_dir = af3_base_dir / f"fold_vwf_{variant_id.lower()}"
                if expected_dir.exists():
                    af3_dir = expected_dir

            # 处理变异
            features = self.process_variant(
                variant_id=variant_id,
                position=int(row['position']),
                ref_aa=row['ref_aa'],
                alt_aa=row['alt_aa'],
                af3_variant_dir=af3_dir
            )

            records.append(features)

            if (idx + 1) % 10 == 0:
                print(f"Processed {idx + 1}/{len(variants_df)} variants")

        return pd.DataFrame(records)

    def generate_feature_summary(self, features_df: pd.DataFrame) -> str:
        """
        生成特征汇总报告
        """
        report = []
        report.append("=" * 80)
        report.append("VWF Residue-Level Feature Extraction Summary")
        report.append("=" * 80)
        report.append("")

        # 基本统计
        report.append(f"Total variants processed: {len(features_df)}")
        report.append(f"Variants with AF3 data: {features_df['af3_available'].sum()}")
        report.append("")

        # 域分布
        report.append("-" * 80)
        report.append("Domain Distribution")
        report.append("-" * 80)
        domain_counts = features_df['domain'].value_counts()
        for domain, count in domain_counts.items():
            report.append(f"  {domain}: {count}")
        report.append("")

        # 功能特征统计
        report.append("-" * 80)
        report.append("Functional Feature Statistics")
        report.append("-" * 80)
        report.append(f"In AIM: {features_df['is_in_aim'].sum()}")
        report.append(f"In GPIb interface: {features_df['is_gpib_interface'].sum()}")
        report.append(f"At scissile bond: {features_df['is_scissile_bond'].sum()}")
        report.append(f"VWD hotspots: {features_df['is_vwd_hotspot'].sum()}")
        report.append("")

        # 突变类型分布
        if 'vwd_hotspot_type' in features_df.columns:
            report.append("-" * 80)
            report.append("VWD Type Distribution")
            report.append("-" * 80)
            type_counts = features_df['vwd_hotspot_type'].value_counts()
            for vwd_type, count in type_counts.items():
                if vwd_type:
                    report.append(f"  {vwd_type}: {count}")
        report.append("")

        # AF3质量统计
        af3_data = features_df[features_df['af3_available'] == True]
        if len(af3_data) > 0:
            report.append("-" * 80)
            report.append("AF3 Structure Quality")
            report.append("-" * 80)
            report.append(f"Mean pLDDT: {af3_data['plddt_mean'].mean():.2f}")
            report.append(f"Mean pTM: {af3_data['ptm_score'].mean():.3f}")
            report.append(f"Mean ranking score: {af3_data['ranking_score'].mean():.3f}")
            report.append(f"Structures with clashes: {af3_data['has_clash'].sum()}")
        report.append("")

        report.append("=" * 80)
        return "\n".join(report)


def create_type2_variant_list() -> pd.DataFrame:
    """
    创建Type 2变异列表用于测试
    基于AF3目录中的变异
    """
    base_dir = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/output/af3_batches_type2/AF3_Results/extracted")

    variants = []
    for var_dir in sorted(base_dir.iterdir()):
        if not var_dir.is_dir() or not var_dir.name.startswith('fold_vwf_'):
            continue

        # 解析变异名
        name = var_dir.name.replace('fold_vwf_', '')

        # 尝试解析 (如 r1306w -> R1306W)
        if len(name) >= 3:
            import re
            match = re.match(r'([a-z]+)(\d+)([a-z]+)', name.lower())
            if match:
                ref_aa = match.group(1).upper()
                if len(ref_aa) > 1:
                    # 3字母转1字母
                    aa_map = {
                        'ala': 'A', 'arg': 'R', 'asn': 'N', 'asp': 'D', 'cys': 'C',
                        'gln': 'Q', 'glu': 'E', 'gly': 'G', 'his': 'H', 'ile': 'I',
                        'leu': 'L', 'lys': 'K', 'met': 'M', 'phe': 'F', 'pro': 'P',
                        'ser': 'S', 'thr': 'T', 'trp': 'W', 'tyr': 'Y', 'val': 'V'
                    }
                    ref_aa = aa_map.get(ref_aa, ref_aa[0])

                position = int(match.group(2))
                alt_aa = match.group(3).upper()
                if len(alt_aa) > 1:
                    alt_aa = aa_map.get(alt_aa, alt_aa[0])

                variants.append({
                    'variant_id': f"{ref_aa}{position}{alt_aa}",
                    'position': position,
                    'ref_aa': ref_aa,
                    'alt_aa': alt_aa,
                    'af3_dir': var_dir
                })

    return pd.DataFrame(variants)


def main():
    """
    测试整合Pipeline
    """
    print("=" * 80)
    print("Integrated VWF Pipeline - Test with Type 2 Variants")
    print("=" * 80)
    print()

    # 创建变异列表
    variants_df = create_type2_variant_list()
    print(f"Found {len(variants_df)} Type 2 variants")
    print()

    # 初始化Pipeline
    af3_base_dir = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/output/af3_batches_type2/AF3_Results/extracted")
    pipeline = IntegratedVWFPipeline(af3_base_dir)

    # 处理前10个变异
    print("Processing variants...")
    features_df = pipeline.process_batch(variants_df.head(10))

    # 生成报告
    print("\n" + pipeline.generate_feature_summary(features_df))

    # 保存结果
    output_file = Path("type2_residue_features.csv")
    features_df.to_csv(output_file, index=False)
    print(f"\nFeatures saved to {output_file}")

    # 显示详细特征
    print("\n" + "=" * 80)
    print("Detailed Features for First 3 Variants")
    print("=" * 80)
    for idx, row in features_df.head(3).iterrows():
        print(f"\n{row['variant_id']} (Position {row['position']}, Domain: {row['domain']}):")
        print(f"  - In AIM: {row['is_in_aim']}")
        print(f"  - GPIb interface: {row['is_gpib_interface']}")
        print(f"  - ADAMTS13 scissile: {row['is_scissile_bond']}")
        print(f"  - VWD hotspot: {row['is_vwd_hotspot']} ({row.get('vwd_hotspot_type', 'N/A')})")
        print(f"  - AF3 available: {row['af3_available']}")
        if row['af3_available']:
            print(f"  - pLDDT at site: {row.get('plddt_at_site', 'N/A')}")
            print(f"  - pLDDT delta: {row.get('plddt_delta', 'N/A')}")


if __name__ == "__main__":
    main()

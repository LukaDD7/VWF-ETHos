#!/usr/bin/env python3
"""
VWF Residue Feature Visualizer
残基级特征可视化模块

功能：
1. 特征热图 (domain-level heatmap)
2. 残基功能注释图
3. pLDDT分布图
4. 突变性质对比图

Author: Claude Code
Date: 2026-04-03
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import seaborn as sns

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class VWFResidueVisualizer:
    """
    VWF残基特征可视化器
    """

    def __init__(self, output_dir: Path = Path("./figures")):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # VWF域架构（用于绘图）
        self.domain_architecture = {
            "signal": (1, 22, "#808080"),
            "D1": (23, 386, "#E74C3C"),      # Propeptide D1
            "D2": (387, 763, "#C0392B"),     # Propeptide D2
            "D'": (764, 865, "#F39C12"),     # D' (FVIII binding)
            "D3": (866, 1233, "#E67E22"),    # D3 (FVIII binding)
            "A1": (1271, 1492, "#27AE60"),   # A1 (GPIb, Type 2B)
            "A2": (1493, 1684, "#2ECC71"),   # A2 (ADAMTS13, Type 2A)
            "A3": (1685, 1874, "#16A085"),   # A3 (Collagen, Type 2M)
            "D4": (1875, 2255, "#3498DB"),   # D4 (Multimerization)
            "C1": (2256, 2324, "#9B59B6"),   # C1
            "C2": (2325, 2392, "#8E44AD"),   # C2
            "C3": (2393, 2496, "#7D3C98"),   # C3
            "C4": (2497, 2577, "#6C3483"),   # C4 (RGD motif)
            "C5": (2578, 2658, "#5B2C6F"),   # C5
            "C6": (2659, 2722, "#4A235A"),   # C6
            "CK": (2723, 2813, "#1ABC9C"),   # CK (Dimerization)
        }

        # 功能位点颜色
        self.feature_colors = {
            "AIM": "#FF0000",
            "GPIb_interface": "#00FF00",
            "ADAMTS13_scissile": "#0000FF",
            "ADAMTS13_exosite": "#0099FF",
            "Calcium_site": "#FFFF00",
            "Collagen_binding": "#FF00FF",
            "FVIII_binding": "#FFA500",
            "RGD_motif": "#00FFFF",
            "VWD_hotspot": "#FF4500",
        }

    def plot_domain_landscape(self, features_df: pd.DataFrame,
                             title: str = "VWF Domain Landscape with Functional Features",
                             output_file: Optional[str] = None) -> None:
        """
        绘制VWF域景观图，显示功能特征分布
        """
        fig, ax = plt.subplots(figsize=(20, 8))

        # 绘制域背景
        y_pos = 0
        for domain, (start, end, color) in self.domain_architecture.items():
            width = end - start + 1
            rect = Rectangle((start, y_pos), width, 1,
                           facecolor=color, edgecolor='black', linewidth=1,
                           alpha=0.6, label=domain)
            ax.add_patch(rect)

            # 添加域标签
            if width > 100:
                ax.text(start + width/2, y_pos + 0.5, domain,
                       ha='center', va='center', fontsize=9, fontweight='bold')

        # 绘制功能特征标记
        feature_y = 1.2

        # AIM区域
        ax.axvspan(1238, 1268, ymin=0.9, ymax=0.95, color=self.feature_colors["AIM"],
                  label="AIM N-term", alpha=0.8)
        ax.axvspan(1460, 1472, ymin=0.9, ymax=0.95, color=self.feature_colors["AIM"],
                  label="AIM C-term", alpha=0.8)

        # GPIbα界面
        ax.axvspan(1296, 1350, ymin=0.85, ymax=0.9, color=self.feature_colors["GPIb_interface"],
                  label="GPIbα interface", alpha=0.8)

        # ADAMTS13切割位点
        ax.axvspan(1605, 1606, ymin=0.8, ymax=0.85, color=self.feature_colors["ADAMTS13_scissile"],
                  label="Scissile bond", alpha=0.9, linewidth=2)
        ax.axvspan(1594, 1622, ymin=0.75, ymax=0.8, color=self.feature_colors["ADAMTS13_exosite"],
                  label="Exosite 1", alpha=0.7)

        # 标记VWD热点
        if 'is_vwd_hotspot' in features_df.columns:
            hotspots = features_df[features_df['is_vwd_hotspot'] == True]
            for _, row in hotspots.iterrows():
                pos = row['position']
                vwd_type = row.get('vwd_hotspot_type', 'Unknown')
                ax.plot(pos, 1.1, 'v', markersize=8,
                       color=self.feature_colors.get("VWD_hotspot", "#FF4500"))

        # 设置图表
        ax.set_xlim(0, 2815)
        ax.set_ylim(-0.5, 1.5)
        ax.set_xlabel("Amino Acid Position", fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')

        # 添加图例
        handles = [mpatches.Patch(color=color, label=domain)
                  for domain, (_, _, color) in list(self.domain_architecture.items())[:8]]
        ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.15, 1))

        ax.set_yticks([])
        ax.grid(axis='x', alpha=0.3)

        plt.tight_layout()

        if output_file:
            plt.savefig(self.output_dir / output_file, dpi=300, bbox_inches='tight')
            print(f"Saved: {self.output_dir / output_file}")
        else:
            plt.savefig(self.output_dir / "vwf_domain_landscape.png", dpi=300, bbox_inches='tight')

        plt.close()

    def plot_feature_heatmap(self, features_df: pd.DataFrame,
                            output_file: Optional[str] = None) -> None:
        """
        绘制特征热图
        """
        # 选择要展示的特征
        feature_cols = [
            'is_in_aim', 'is_gpib_interface', 'is_scissile_bond',
            'is_exosite_1', 'is_exosite_2', 'is_exosite_3',
            'is_calcium_coordinating', 'is_collagen_binding',
            'is_fviii_binding', 'is_rgd_motif', 'is_vwd_hotspot'
        ]

        # 过滤存在的列
        available_cols = [col for col in feature_cols if col in features_df.columns]

        if not available_cols:
            print("No feature columns available for heatmap")
            return

        # 按域分组计算特征频率
        domain_features = features_df.groupby('domain')[available_cols].mean()

        # 绘制热图
        fig, ax = plt.subplots(figsize=(12, 8))

        sns.heatmap(domain_features.T, annot=True, fmt='.2f',
                   cmap='YlOrRd', cbar_kws={'label': 'Feature Frequency'},
                   ax=ax, linewidths=0.5)

        ax.set_xlabel('VWF Domain', fontsize=12)
        ax.set_ylabel('Functional Feature', fontsize=12)
        ax.set_title('VWF Functional Feature Distribution by Domain', fontsize=14, fontweight='bold')

        # 旋转y轴标签
        plt.setp(ax.get_yticklabels(), rotation=0)

        plt.tight_layout()

        if output_file:
            plt.savefig(self.output_dir / output_file, dpi=300, bbox_inches='tight')
        else:
            plt.savefig(self.output_dir / "vwf_feature_heatmap.png", dpi=300, bbox_inches='tight')

        plt.close()

    def plot_plddt_distribution(self, features_df: pd.DataFrame,
                               output_file: Optional[str] = None) -> None:
        """
        绘制pLDDT分布图
        """
        if 'plddt_mean' not in features_df.columns or features_df['plddt_mean'].isna().all():
            print("No pLDDT data available")
            return

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # 1. pLDDT按域分布
        ax1 = axes[0, 0]
        af3_data = features_df[features_df['af3_available'] == True]
        if len(af3_data) > 0:
            domain_plddt = af3_data.groupby('domain')['plddt_mean'].mean().sort_values(ascending=False)
            domain_plddt.plot(kind='bar', ax=ax1, color='steelblue')
            ax1.set_xlabel('Domain')
            ax1.set_ylabel('Mean pLDDT')
            ax1.set_title('Mean pLDDT by Domain')
            ax1.tick_params(axis='x', rotation=45)
            ax1.axhline(y=70, color='r', linestyle='--', label='Quality threshold (70)')
            ax1.legend()

        # 2. pLDDT分布直方图
        ax2 = axes[0, 1]
        if len(af3_data) > 0:
            ax2.hist(af3_data['plddt_mean'], bins=20, color='steelblue', edgecolor='black', alpha=0.7)
            ax2.set_xlabel('Mean pLDDT')
            ax2.set_ylabel('Count')
            ax2.set_title('Distribution of Mean pLDDT')
            ax2.axvline(x=70, color='r', linestyle='--', label='Threshold (70)')
            ax2.legend()

        # 3. pLDDT Delta分布（如果可用）
        ax3 = axes[1, 0]
        if 'plddt_delta' in af3_data.columns and not af3_data['plddt_delta'].isna().all():
            delta_data = af3_data['plddt_delta'].dropna()
            ax3.hist(delta_data, bins=20, color='coral', edgecolor='black', alpha=0.7)
            ax3.set_xlabel('pLDDT Delta (Mut - WT)')
            ax3.set_ylabel('Count')
            ax3.set_title('Distribution of pLDDT Changes')
            ax3.axvline(x=0, color='black', linestyle='-', linewidth=2)
            ax3.axvline(x=-5, color='r', linestyle='--', label='Significant destabilization')

        # 4. pTM分数分布
        ax4 = axes[1, 1]
        if 'ptm_score' in af3_data.columns:
            ax4.scatter(af3_data['plddt_mean'], af3_data['ptm_score'],
                       alpha=0.6, s=50, c='green')
            ax4.set_xlabel('Mean pLDDT')
            ax4.set_ylabel('pTM Score')
            ax4.set_title('pLDDT vs pTM Score')

        plt.tight_layout()

        if output_file:
            plt.savefig(self.output_dir / output_file, dpi=300, bbox_inches='tight')
        else:
            plt.savefig(self.output_dir / "vwf_plddt_analysis.png", dpi=300, bbox_inches='tight')

        plt.close()

    def plot_mutation_properties(self, features_df: pd.DataFrame,
                                output_file: Optional[str] = None) -> None:
        """
        绘制突变性质分布图
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # 1. 电荷变化分布
        ax1 = axes[0, 0]
        if 'mutation_charge_change' in features_df.columns:
            charge_counts = features_df['mutation_charge_change'].value_counts().sort_index()
            colors = ['red', 'lightcoral', 'gray', 'lightblue', 'blue']
            charge_counts.plot(kind='bar', ax=ax1, color=colors[:len(charge_counts)])
            ax1.set_xlabel('Charge Change')
            ax1.set_ylabel('Count')
            ax1.set_title('Distribution of Charge Changes')
            ax1.tick_params(axis='x', rotation=0)

        # 2. 大小变化分布
        ax2 = axes[0, 1]
        if 'mutation_size_delta' in features_df.columns:
            ax2.hist(features_df['mutation_size_delta'], bins=30,
                    color='steelblue', edgecolor='black', alpha=0.7)
            ax2.set_xlabel('Size Change (Å³)')
            ax2.set_ylabel('Count')
            ax2.set_title('Distribution of Size Changes')
            ax2.axvline(x=0, color='black', linestyle='-', linewidth=2)

        # 3. 疏水性变化
        ax3 = axes[1, 0]
        if 'mutation_hydrophobicity_delta' in features_df.columns:
            ax3.hist(features_df['mutation_hydrophobicity_delta'], bins=30,
                    color='coral', edgecolor='black', alpha=0.7)
            ax3.set_xlabel('Hydrophobicity Change (Kyte-Doolittle)')
            ax3.set_ylabel('Count')
            ax3.set_title('Distribution of Hydrophobicity Changes')
            ax3.axvline(x=0, color='black', linestyle='-', linewidth=2)

        # 4. 突变类型与域的关系
        ax4 = axes[1, 1]
        if 'is_vwd_hotspot' in features_df.columns:
            domain_hotspot = features_df.groupby('domain')['is_vwd_hotspot'].sum()
            domain_hotspot.plot(kind='barh', ax=ax4, color='orange')
            ax4.set_xlabel('Number of VWD Hotspots')
            ax4.set_title('VWD Hotspots by Domain')

        plt.tight_layout()

        if output_file:
            plt.savefig(self.output_dir / output_file, dpi=300, bbox_inches='tight')
        else:
            plt.savefig(self.output_dir / "vwf_mutation_properties.png", dpi=300, bbox_inches='tight')

        plt.close()

    def create_comprehensive_report(self, features_df: pd.DataFrame,
                                    output_prefix: str = "vwf_analysis") -> None:
        """
        创建综合分析报告
        """
        print("Generating comprehensive visualizations...")

        # 1. 域景观图
        self.plot_domain_landscape(features_df,
                                   output_file=f"{output_prefix}_domain_landscape.png")

        # 2. 特征热图
        self.plot_feature_heatmap(features_df,
                                 output_file=f"{output_prefix}_feature_heatmap.png")

        # 3. pLDDT分布
        self.plot_plddt_distribution(features_df,
                                    output_file=f"{output_prefix}_plddt_analysis.png")

        # 4. 突变性质
        self.plot_mutation_properties(features_df,
                                     output_file=f"{output_prefix}_mutation_properties.png")

        print(f"All visualizations saved to {self.output_dir}")


def main():
    """测试可视化功能"""
    import sys
    sys.path.insert(0, str(Path(__file__).parent))

    # 加载特征数据（如果存在）
    feature_file = Path("type2_residue_features.csv")
    if feature_file.exists():
        features_df = pd.read_csv(feature_file)
        print(f"Loaded {len(features_df)} variants from {feature_file}")

        # 创建可视化
        visualizer = VWFResidueVisualizer(output_dir=Path("./figures"))
        visualizer.create_comprehensive_report(features_df, output_prefix="type2_variants")
    else:
        print(f"Feature file not found: {feature_file}")
        print("Run vwf_integrated_pipeline.py first to generate features.")


if __name__ == "__main__":
    main()

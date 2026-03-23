#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本 4：分析与可视化（临床级高保真切片版）
目标：生成符合国际高水平期刊要求的多模态差值图谱

输入：../results/03_inference_results.pkl, ../results/03_inference_results.csv
输出：../results/04_analysis_summary.csv, figures/*.png
"""

import logging
import pickle
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ==================== InferenceResult 定义 ====================
@dataclass
class InferenceResult:
    chromosome: str
    position: int
    ref: str
    alt: str
    interval_start: int
    interval_end: int
    clinvar_classification: str
    consequences: str
    rna_seq_delta_max: float | None = None
    rna_seq_delta_sum: float | None = None
    splice_site_delta_max: float | None = None
    splice_site_delta_sum: float | None = None
    delta_score: float = 0.0
    raw_outputs: Any = field(default=None)
    success: bool = False
    error_message: str | None = None

# ==================== 配置区域 ====================
# 基础目录配置
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"

INPUT_PKL = RESULTS_DIR / "03_inference_results.pkl"
INPUT_CSV = RESULTS_DIR / "03_inference_results.csv"
OUTPUT_DIR = RESULTS_DIR / "figures"
TOP_N_VARIANTS = 20

FIGURE_DPI = 300
FIGURE_SIZE = (14, 12)  # 黄金比例，避免过度拉伸
FOCUS_REGION_WIDTH = 2**15  # 32kb

# 临床学术制图色彩
COLOR_REF = "#8C8C8C"   # 沉稳的灰色
COLOR_ALT = "#D9383A"   # 敏锐的警示红
COLOR_DELTA = "#533475" # 深邃的紫色差值

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def load_results(pkl_path: Path, csv_path: Path):
    with open(pkl_path, "rb") as f:
        full_results = pickle.load(f)
    summary_df = pd.read_csv(csv_path)
    return full_results, summary_df

def get_top_variants(summary_df: pd.DataFrame, n: int) -> pd.DataFrame:
    success_df = summary_df[summary_df["success"] == True]
    return success_df.nlargest(n, "delta_score")

def visualize_variant_dual_modal(result, output_dir: Path):
    if result.raw_outputs is None:
        return None

    outputs = result.raw_outputs
    has_rna = (outputs.reference.rna_seq is not None and outputs.reference.rna_seq.values.size > 0)
    has_splice = (outputs.reference.splice_sites is not None and outputs.reference.splice_sites.values.size > 0)

    if not has_rna and not has_splice:
        return None

    output_dir.mkdir(exist_ok=True)
    filename = f"{result.chromosome}_{result.position}_{result.ref}_{result.alt}.png"
    output_path = output_dir / filename

    # 物理坐标切片计算
    CENTER_IDX = 524288
    HALF_WINDOW = 16384
    start_idx = CENTER_IDX - HALF_WINDOW
    end_idx = CENTER_IDX + HALF_WINDOW
    x_axis = np.arange(result.position - HALF_WINDOW, result.position + HALF_WINDOW)

    fig, axes = plt.subplots(4, 1, figsize=FIGURE_SIZE, dpi=FIGURE_DPI, sharex=True)

    # ==================== 1. RNA-Seq (生物学重复求均值) ====================
    if has_rna:
        ref_rna_mean = np.mean(outputs.reference.rna_seq.values[start_idx:end_idx, :], axis=1)
        alt_rna_mean = np.mean(outputs.alternate.rna_seq.values[start_idx:end_idx, :], axis=1)

        # 粗灰细红策略
        axes[0].plot(x_axis, ref_rna_mean, color=COLOR_REF, linewidth=2.5, alpha=0.8, label='REF (Mean)')
        axes[0].plot(x_axis, alt_rna_mean, color=COLOR_ALT, linewidth=1.2, alpha=0.9, label='ALT (Mean)')
        axes[0].set_title("RNA-Seq Expression: Reference vs Mutated", fontsize=11, fontweight='bold', loc='left')
        axes[0].set_ylabel("Read Coverage")
        axes[0].legend(loc='upper right', frameon=False)

        # 差值图
        rna_diff = alt_rna_mean - ref_rna_mean
        axes[1].fill_between(x_axis, rna_diff, color=COLOR_DELTA, alpha=0.6)
        axes[1].axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        axes[1].set_title(f"RNA-Seq Difference (ALT - REF) | Max Δ: {result.rna_seq_delta_max:.3f}", fontsize=10, loc='left')
        axes[1].set_ylabel("Δ Expression")
    else:
        axes[0].text(0.5, 0.5, "No RNA-Seq Data", ha='center', va='center', transform=axes[0].transAxes)
        axes[1].text(0.5, 0.5, "No RNA-Seq Data", ha='center', va='center', transform=axes[1].transAxes)

    # ==================== 2. Splice Sites (保留 4 个独立 Track) ====================
    if has_splice:
        ref_splice = outputs.reference.splice_sites.values[start_idx:end_idx, :] # Shape: (32768, 4)
        alt_splice = outputs.alternate.splice_sites.values[start_idx:end_idx, :]

        # 绝不使用 fill_between，避免红色淹没。使用纯线条，且粗灰细红
        axes[2].plot(x_axis, ref_splice, color=COLOR_REF, linewidth=2.5, alpha=0.7)
        axes[2].plot(x_axis, alt_splice, color=COLOR_ALT, linewidth=1.0, alpha=0.9)

        # 伪造图例标识
        axes[2].plot([], [], color=COLOR_REF, linewidth=2.5, label='REF (4 Tracks)')
        axes[2].plot([], [], color=COLOR_ALT, linewidth=1.0, label='ALT (4 Tracks)')
        axes[2].set_title("Splice Sites Probabilities (All 4 Strands)", fontsize=11, fontweight='bold', loc='left')
        axes[2].set_ylabel("Probability")
        axes[2].legend(loc='upper right', frameon=False)

        # 差值图 (依然保留 4 个维度，找出哪个位点被破坏了)
        splice_diff = alt_splice - ref_splice
        axes[3].plot(x_axis, splice_diff, color=COLOR_DELTA, linewidth=1.2, alpha=0.8)
        axes[3].axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        axes[3].set_title(f"Splice Sites Difference (ALT - REF) | Max Δ: {result.splice_site_delta_max:.3f}", fontsize=10, loc='left')
        axes[3].set_ylabel("Δ Probability")
        axes[3].set_xlabel(f"Chromosome {result.chromosome} Position (hg38)", fontweight='bold')
    else:
        axes[2].text(0.5, 0.5, "No Splice Data", ha='center', va='center', transform=axes[2].transAxes)
        axes[3].text(0.5, 0.5, "No Splice Data", ha='center', va='center', transform=axes[3].transAxes)

    # ==================== 全局标注与防拥挤排版 ====================
    for ax in axes:
        ax.grid(True, alpha=0.2, linestyle='--')
        ax.axvline(x=result.position, color='black', linestyle=':', linewidth=1.5, alpha=0.6)
        ax.ticklabel_format(useOffset=False, style='plain', axis='x')

    # 优雅地拆分超长标题
    conseq_wrapped = textwrap.fill(result.consequences, width=100)
    title_text = (
        f"Variant: {result.chromosome}:{result.position} {result.ref} > {result.alt}  |  "
        f"ClinVar: {result.clinvar_classification}\n"
        f"Consequences: {conseq_wrapped}"
    )
    fig.suptitle(title_text, fontsize=13, y=0.97, fontweight='bold', family='sans-serif')

    plt.tight_layout(rect=[0, 0, 1, 0.94]) # 为主标题留出空间
    plt.savefig(output_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return output_path

def create_analysis_summary(summary_df: pd.DataFrame, output_path: Path):
    success_df = summary_df[summary_df["success"] == True]
    ranked_df = success_df.sort_values(by="delta_score", ascending=False).reset_index(drop=True)
    columns_to_keep = [
        "chromosome", "position", "ref", "alt", "clinvar_classification",
        "consequences", "delta_score", "rna_seq_delta_max", "splice_site_delta_max",
    ]
    existing_cols = [c for c in columns_to_keep if c in ranked_df.columns]
    ranked_df[existing_cols].to_csv(output_path, index=False)
    logger.info(f"成功导出医学级靶点清单，共 {len(ranked_df)} 个记录。")

def main():
    logger.info("=" * 50)
    logger.info("开始生成临床级数字病理图谱...")

    pkl_path = INPUT_PKL
    csv_path = INPUT_CSV
    full_results, summary_df = load_results(pkl_path, csv_path)
    top_variants = get_top_variants(summary_df, TOP_N_VARIANTS)

    output_dir = OUTPUT_DIR
    output_dir.mkdir(exist_ok=True)

    for _, row in top_variants.iterrows():
        result = next((r for r in full_results if r.chromosome == row["chromosome"] and r.position == row["position"] and r.ref == row["ref"] and r.alt == row["alt"]), None)
        if result:
            visualize_variant_dual_modal(result, output_dir)

    create_analysis_summary(summary_df, RESULTS_DIR / "04_analysis_summary.csv")
    logger.info("可视化分析全部完成！")

if __name__ == "__main__":
    main()

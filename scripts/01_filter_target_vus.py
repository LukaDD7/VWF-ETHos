#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本 1：过滤目标 VUS（意义未明变异）

目标：从原始 ClinVar+HGMD 融合数据中提取最有研究价值的靶点。
筛选策略：
  1. 临床分类为 'Uncertain_significance'（意义未明）
  2. 突变后果类型为 intron_variant、synonymous_variant 或 upstream_gene_variant
     （排除纯外显子错义突变，关注非编码区和同义突变）

输入：../data/Clinvar_HGMD_merge_annotated.xlsx 或 .csv
输出：../data/01_filtered_targets.csv
"""

import logging
import os
import sys
from pathlib import Path

import pandas as pd

# ==================== 配置区域 ====================
# 基础目录配置
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
DATA_DIR = BASE_DIR / "data"
LOGS_DIR = BASE_DIR / "logs"

# 支持的输入文件格式（在 data/ 目录中查找）
INPUT_FILES = [
    DATA_DIR / "Clinvar_HGMD_merge_annotated.xlsx - Sheet1.csv",  # CSV 格式
    DATA_DIR / "Clinvar_HGMD_merge_annotated.xlsx",               # Excel 格式
    DATA_DIR / "Clinvar_HGMD_merge_annotated.csv",                # 备选 CSV 名
]
OUTPUT_FILE = DATA_DIR / "01_filtered_targets.csv"

# 列名映射（根据实际数据调整）
COL_CLINVAR_CLASSIFICATION = "INFO_clinvar_classification_base"
COL_CONSEQUENCES = "INFO_consequences_base"

# 筛选条件
# 保留三类：VUS + 已知致病突变 (阳性对照)
TARGET_CLINICAL_CLASSIFICATIONS = [
    "Uncertain_significance",  # 意义未明 (主要研究目标)
    "Pathogenic",              # 致病突变 (阳性对照)
    "Likely_pathogenic",       # 可能致病 (阳性对照)
]
TARGET_CONSEQUENCES = [
    "intron_variant",
    "synonymous_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
]

# 确保日志目录存在
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# 日志配置
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS_DIR / "01_filter.log", mode="w"),
    ],
)
logger = logging.getLogger(__name__)


def find_input_file() -> Path | None:
    """查找存在的输入文件。"""
    for filepath in INPUT_FILES:
        if filepath.exists():
            return filepath
    return None


def load_data(filepath: Path) -> pd.DataFrame:
    """加载输入数据，支持 xlsx 和 csv 格式。"""
    logger.info(f"正在加载文件：{filepath}")

    if filepath.suffix == ".xlsx":
        df = pd.read_excel(filepath)
    elif filepath.suffix == ".csv":
        df = pd.read_csv(filepath)
    else:
        raise ValueError(f"不支持的文件格式：{filepath.suffix}")

    logger.info(f"加载完成：{len(df)} 行，{len(df.columns)} 列")
    return df


def filter_vus(df: pd.DataFrame) -> pd.DataFrame:
    """
    筛选目标变异。

    筛选逻辑：
    1. 临床分类为目标类型（Uncertain_significance, Pathogenic, Likely_pathogenic）
    2. 突变后果包含目标类型（如 intron_variant 等）
    """
    initial_count = len(df)
    logger.info(f"初始数据量：{initial_count} 个变异")

    # 检查必需列是否存在
    required_cols = [COL_CLINVAR_CLASSIFICATION, COL_CONSEQUENCES]
    for col in required_cols:
        if col not in df.columns:
            raise KeyError(f"必需列 '{col}' 不存在。可用列：{list(df.columns)}")

    # 筛选条件 A：临床分类（保留 VUS + 致病突变作为阳性对照）
    mask_classification = df[COL_CLINVAR_CLASSIFICATION].isin(TARGET_CLINICAL_CLASSIFICATIONS)
    logger.info(f"临床分类为 {TARGET_CLINICAL_CLASSIFICATIONS} 的变异数：{mask_classification.sum()}")

    # 输出临床分类分布
    logger.info("临床分类分布:")
    for cls in TARGET_CLINICAL_CLASSIFICATIONS:
        count = (df[COL_CLINVAR_CLASSIFICATION] == cls).sum()
        logger.info(f"  - {cls}: {count}")

    # 筛选条件 B：突变后果类型
    def matches_consequences(value) -> bool:
        """检查突变后果是否匹配目标类型。"""
        if pd.isna(value):
            return False
        value_str = str(value)
        return any(conseq in value_str for conseq in TARGET_CONSEQUENCES)

    mask_consequences = df[COL_CONSEQUENCES].apply(matches_consequences)
    logger.info(f"突变后果匹配目标类型的变异数：{mask_consequences.sum()}")

    # 组合筛选条件
    combined_mask = mask_classification & mask_consequences
    filtered_df = df[combined_mask].reset_index(drop=True)

    logger.info(f"筛选后数据量：{len(filtered_df)} 个变异")
    logger.info(f"筛选掉的变异数：{initial_count - len(filtered_df)}")

    # 确保至少有 1 个致病突变作为阳性对照
    pathogenic_count = (filtered_df[COL_CLINVAR_CLASSIFICATION] == "Pathogenic").sum()
    likely_pathogenic_count = (filtered_df[COL_CLINVAR_CLASSIFICATION] == "Likely_pathogenic").sum()
    logger.info(f"阳性对照 (Pathogenic + Likely_pathogenic): {pathogenic_count + likely_pathogenic_count} 个")

    if pathogenic_count + likely_pathogenic_count == 0:
        logger.warning("警告：没有筛选到致病突变作为阳性对照！")

    return filtered_df


def save_output(df: pd.DataFrame, output_path: Path) -> None:
    """保存筛选后的数据到 CSV。"""
    df.to_csv(output_path, index=False)
    logger.info(f"已保存筛选结果到：{output_path}")


def main():
    """主函数。"""
    logger.info("=" * 60)
    logger.info("开始执行脚本 1：过滤目标 VUS")
    logger.info("=" * 60)

    # 查找输入文件
    input_path = find_input_file()
    if input_path is None:
        logger.error("未找到输入文件。请确保以下文件之一存在：")
        for f in INPUT_FILES:
            logger.error(f"  - {f}")
        sys.exit(1)

    logger.info(f"找到输入文件：{input_path}")

    # 加载数据
    try:
        df = load_data(input_path)
    except Exception as e:
        logger.error(f"加载数据失败：{e}")
        sys.exit(1)

    # 筛选 VUS
    try:
        filtered_df = filter_vus(df)
    except Exception as e:
        logger.error(f"筛选数据失败：{e}")
        sys.exit(1)

    # 保存结果
    output_path = OUTPUT_FILE
    save_output(filtered_df, output_path)

    # 输出统计摘要
    logger.info("=" * 60)
    logger.info("筛选完成！")
    logger.info(f"输出文件：{output_path}")
    logger.info(f"筛选出的 VUS 数量：{len(filtered_df)}")

    if len(filtered_df) > 0:
        logger.info("临床分类分布：")
        logger.info(f"\n{filtered_df[COL_CLINVAR_CLASSIFICATION].value_counts()}")
        logger.info("突变后果类型分布（前 10）：")
        logger.info(f"\n{filtered_df[COL_CONSEQUENCES].value_counts().head(10)}")

    logger.info("=" * 60)


if __name__ == "__main__":
    main()

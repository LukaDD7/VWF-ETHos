#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本 6b：修复 Final_VWF_Target_List 中缺失的 AlphaGenome 数据

问题：
- 280个错义突变在Final_VWF_Target_List中缺失AlphaGenome数据
- 但07结果CSV中实际上有275个存在
- 原因是06脚本使用PKL文件，可能存在类定义不匹配问题

解决方案：
- 直接从07 CSV读取数据（绕过PKL）
- 基于hg19坐标重新匹配
- 补齐AlphaGenome_Max_Score等关键列
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path

# ==================== 配置 ====================
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"

INPUT_EXCEL = RESULTS_DIR / "Final_VWF_Target_List_with_AlphaGenome.xlsx"
AG_CSV = RESULTS_DIR / "07_VCF_AlphaGenome_Results.csv"
OUTPUT_EXCEL = RESULTS_DIR / "Final_VWF_Target_List_with_AlphaGenome_FIXED.xlsx"

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)


def main():
    logger.info("=" * 70)
    logger.info("修复 Final_VWF_Target_List 中缺失的 AlphaGenome 数据")
    logger.info("=" * 70)

    # 1. 读取输入文件
    logger.info(f"1. 读取医生表格: {INPUT_EXCEL}")
    df_doctor = pd.read_excel(INPUT_EXCEL)
    logger.info(f"   总行数: {len(df_doctor)}")

    logger.info(f"2. 读取07结果CSV: {AG_CSV}")
    df_ag = pd.read_csv(AG_CSV)
    logger.info(f"   总行数: {len(df_ag)}")

    # 2. 为两个数据集创建匹配键（基于hg19坐标）
    logger.info("3. 创建匹配键...")

    # 医生表格的匹配键（GRCh37Location是hg19坐标）
    df_doctor['match_key'] = (
        df_doctor['GRCh37Chromosome'].astype(str).str.replace('.0', '', regex=False) + '_' +
        df_doctor['GRCh37Location'].astype(str).str.replace('.0', '', regex=False)
    )

    # 07结果的匹配键
    df_ag['match_key'] = df_ag['hg19_chr'].astype(str) + '_' + df_ag['hg19_pos'].astype(str)

    # 去重：保留第一个出现的
    df_ag_dedup = df_ag.drop_duplicates(subset=['match_key'], keep='first')
    logger.info(f"   07结果去重后: {len(df_ag_dedup)} (去掉了 {len(df_ag) - len(df_ag_dedup)} 个重复)")

    # 3. 找出缺失AlphaGenome数据的错义突变
    missense_mask = df_doctor['Molecular.consequence'] == 'missense variant'
    no_ag_mask = df_doctor['AlphaGenome_Max_Score'].isna()
    missing_ag = df_doctor[missense_mask & no_ag_mask]

    logger.info(f"4. 无AlphaGenome数据的错义突变: {len(missing_ag)} 个")

    # 4. 基于match_key进行左连接，补齐数据
    logger.info("5. 合并数据...")

    # 选择需要合并的列
    merge_cols = ['match_key', 'Delta_Max_Core', 'AG_RNA_SEQ', 'AG_SPLICE_SITES',
                  'AG_CAGE', 'AG_PROCAP', 'AG_SPLICE_SITE_USAGE', 'AG_SPLICE_JUNCTIONS',
                  'AG_DNASE', 'AG_ATAC', 'AG_CHIP_HISTONE', 'AG_CHIP_TF', 'AG_CONTACT_MAPS']

    # 只保留存在的列
    merge_cols = [c for c in merge_cols if c in df_ag_dedup.columns]

    # 执行合并
    df_merged = df_doctor.merge(
        df_ag_dedup[merge_cols],
        on='match_key',
        how='left',
        suffixes=('', '_from07')
    )

    # 5. 补齐缺失的AlphaGenome_Max_Score
    before_fix = df_merged['AlphaGenome_Max_Score'].isna().sum()

    # 如果原始AlphaGenome_Max_Score为NaN，但合并后有Delta_Max_Core，则补齐
    # 列名在合并后会自动添加后缀 '_from07'
    delta_col = 'Delta_Max_Core_from07' if 'Delta_Max_Core_from07' in df_merged.columns else 'Delta_Max_Core'
    mask_missing = df_merged['AlphaGenome_Max_Score'].isna() & df_merged[delta_col].notna()
    df_merged.loc[mask_missing, 'AlphaGenome_Max_Score'] = df_merged.loc[mask_missing, delta_col]

    after_fix = df_merged['AlphaGenome_Max_Score'].isna().sum()
    fixed_count = before_fix - after_fix

    logger.info(f"6. 修复统计:")
    logger.info(f"   修复前缺失: {before_fix}")
    logger.info(f"   修复后缺失: {after_fix}")
    logger.info(f"   成功修复: {fixed_count}")

    # 6. 检查错义突变的修复情况
    missense_fixed = df_merged[missense_mask & df_merged['AlphaGenome_Max_Score'].notna()]
    missense_missing = df_merged[missense_mask & df_merged['AlphaGenome_Max_Score'].isna()]

    logger.info(f"\n7. 错义突变修复情况:")
    logger.info(f"   总数: {missense_mask.sum()}")
    logger.info(f"   修复后有数据: {len(missense_fixed)}")
    logger.info(f"   仍无数据: {len(missense_missing)}")

    # 8. 清理临时列
    df_merged = df_merged.drop(columns=[
        'match_key', 'Delta_Max_Core_from07', 'AG_RNA_SEQ_from07',
        'AG_SPLICE_SITES_from07', 'AG_CAGE_from07', 'AG_PROCAP_from07',
        'AG_SPLICE_SITE_USAGE_from07', 'AG_SPLICE_JUNCTIONS_from07',
        'AG_DNASE_from07', 'AG_ATAC_from07', 'AG_CHIP_HISTONE_from07',
        'AG_CHIP_TF_from07', 'AG_CONTACT_MAPS_from07'
    ], errors='ignore')

    # 9. 保存结果
    logger.info(f"\n8. 保存修复后的文件: {OUTPUT_EXCEL}")
    df_merged.to_excel(OUTPUT_EXCEL, index=False)

    logger.info("=" * 70)
    logger.info("修复完成！")
    logger.info(f"输出文件: {OUTPUT_EXCEL}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()

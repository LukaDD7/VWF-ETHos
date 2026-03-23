#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
08_merge_godmode_results.py - 精准补齐5种模态

功能：
1. 读取 07 结果（内皮细胞6模态金标准）
2. 读取 07e 结果（全组织5模态补充）
3. 根据 (chr, pos, ref, alt) 精准匹配
4. 将 07e 的5种模态（ATAC/Histone/CAGE/PROcap/Contact）补齐到 07 结果
5. 保留原始07结果作为备份
6. 生成最终合并结果

输出列说明：
- AG_ATAC -> 来自 07e (AllTissues)
- AG_CHIP_HISTONE -> 来自 07e (AllTissues)
- AG_CAGE -> 来自 07e (AllTissues)
- AG_PROCAP -> 来自 07e (AllTissues)
- AG_CONTACT_MAPS -> 来自 07e (AllTissues)
- 其他列保持 07 原值不变
"""

import logging
import sys
import shutil
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np

# ==================== 配置区域 ====================
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"
LOGS_DIR = BASE_DIR / "logs"

# 输入文件
RESULT_07_CSV = RESULTS_DIR / "07_VCF_AlphaGenome_Results.csv"
RESULT_07_PKL = RESULTS_DIR / "07_VCF_AlphaGenome_Results.pkl"
RESULT_07E_CSV = RESULTS_DIR / "07e_GodMode_Epigenome_Peaks.csv"

# 输出文件
MERGED_CSV = RESULTS_DIR / "08_VCF_AlphaGenome_Results_Merged_11_Modalities.csv"
MERGED_PKL = RESULTS_DIR / "08_VCF_AlphaGenome_Results_Merged_11_Modalities.pkl"
BACKUP_DIR = RESULTS_DIR / "backups"

# 确保目录存在
LOGS_DIR.mkdir(parents=True, exist_ok=True)
BACKUP_DIR.mkdir(parents=True, exist_ok=True)

# 日志配置
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS_DIR / "08_merge_godmode_log.txt", mode='a'),
    ]
)
logger = logging.getLogger(__name__)


def load_07_results():
    """加载 07 结果（优先使用pkl，回退到csv）"""
    if RESULT_07_PKL.exists():
        logger.info(f"📂 从 PKL 加载 07 结果: {RESULT_07_PKL}")
        try:
            import pickle
            with open(RESULT_07_PKL, 'rb') as f:
                data = pickle.load(f)
            if isinstance(data, list):
                df = pd.DataFrame(data)
            elif isinstance(data, pd.DataFrame):
                df = data
            else:
                raise ValueError(f"不支持的 07 数据类型: {type(data)}")
            logger.info(f"   ✓ 加载 {len(df)} 条记录")
            return df
        except Exception as e:
            logger.warning(f"   ⚠️ PKL 加载失败: {e}，回退到 CSV")

    logger.info(f"📂 从 CSV 加载 07 结果: {RESULT_07_CSV}")
    df = pd.read_csv(RESULT_07_CSV)
    logger.info(f"   ✓ 加载 {len(df)} 条记录")
    return df


def load_07e_results():
    """加载 07e 结果（只取成功的记录）"""
    logger.info(f"📂 加载 07e 结果: {RESULT_07E_CSV}")

    if not RESULT_07E_CSV.exists():
        logger.error(f"❌ 07e 结果文件不存在: {RESULT_07E_CSV}")
        sys.exit(1)

    df = pd.read_csv(RESULT_07E_CSV)
    logger.info(f"   ✓ 加载 {len(df)} 条记录")

    # 只保留成功的记录
    if 'God_Mode_Status' in df.columns:
        df_success = df[df['God_Mode_Status'] == 'Success'].copy()
        logger.info(f"   ✓ 其中成功: {len(df_success)} 条，失败: {len(df) - len(df_success)} 条")
        return df_success
    else:
        logger.warning("   ⚠️ 未找到 God_Mode_Status 列，使用全部记录")
        return df


def create_variant_key(df, chr_col, pos_col, ref_col, alt_col):
    """创建变异的唯一标识键"""
    return (
        df[chr_col].astype(str) + ':' +
        df[pos_col].astype(str) + ':' +
        df[ref_col].astype(str) + ':' +
        df[alt_col].astype(str)
    )


def backup_original_files():
    """备份原始 07 结果文件"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    for src_file in [RESULT_07_CSV, RESULT_07_PKL]:
        if src_file.exists():
            dst_file = BACKUP_DIR / f"{src_file.stem}_pre_merge_{timestamp}{src_file.suffix}"
            shutil.copy2(src_file, dst_file)
            logger.info(f"💾 备份: {src_file.name} -> {dst_file.name}")


def merge_results(df_07, df_07e):
    """
    精准合并 07 和 07e 结果

    策略：
    1. 为两个DataFrame创建variant_key
    2. 根据variant_key进行左连接（保留所有07记录）
    3. 将07e的5种模态填充到07的对应列
    """
    logger.info("\n🔧 开始合并...")

    # 创建variant_key
    df_07 = df_07.copy()
    df_07e = df_07e.copy()

    # 07的列名映射
    df_07['variant_key'] = create_variant_key(
        df_07, 'hg38_chr', 'hg38_pos', 'ref', 'alt'
    )

    # 07e的列名映射
    df_07e['variant_key'] = create_variant_key(
        df_07e, 'Chromosome', 'Position', 'Ref', 'Alt'
    )

    # 检查重复
    dup_07 = df_07['variant_key'].duplicated().sum()
    dup_07e = df_07e['variant_key'].duplicated().sum()
    if dup_07 > 0:
        logger.warning(f"   ⚠️ 07 结果中有 {dup_07} 个重复变异，将去重")
        df_07 = df_07.drop_duplicates(subset=['variant_key'], keep='first')
    if dup_07e > 0:
        logger.warning(f"   ⚠️ 07e 结果中有 {dup_07e} 个重复变异，将去重")
        df_07e = df_07e.drop_duplicates(subset=['variant_key'], keep='first')

    # 左连接（保留所有07记录）
    logger.info(f"   07 总记录: {len(df_07)}")
    logger.info(f"   07e 成功记录: {len(df_07e)}")

    df_merged = df_07.merge(
        df_07e[['variant_key', 'ATAC_Max_Delta', 'Histone_Max_Delta',
                'CAGE_Max_Delta', 'PROCAP_Max_Delta', 'Contact_Max_Delta']],
        on='variant_key',
        how='left',
        suffixes=('', '_07e')
    )

    # 统计匹配情况
    matched = df_merged['ATAC_Max_Delta'].notna().sum()
    logger.info(f"   ✓ 成功匹配: {matched} 条 ({matched/len(df_07)*100:.1f}%)")

    # 补齐5种模态（注意：列名可能有变化，需要确认）
    # 07e列名: ATAC_Max_Delta, Histone_Max_Delta, CAGE_Max_Delta, PROCAP_Max_Delta, Contact_Max_Delta
    # 07列名: AG_ATAC, AG_CHIP_HISTONE, AG_CAGE, AG_PROCAP, AG_CONTACT_MAPS

    # 记录原始值用于比较
    df_merged['AG_ATAC_original'] = df_merged.get('AG_ATAC', 0)
    df_merged['AG_CHIP_HISTONE_original'] = df_merged.get('AG_CHIP_HISTONE', 0)
    df_merged['AG_CAGE_original'] = df_merged.get('AG_CAGE', 0)
    df_merged['AG_PROCAP_original'] = df_merged.get('AG_PROCAP', 0)
    df_merged['AG_CONTACT_MAPS_original'] = df_merged.get('AG_CONTACT_MAPS', 0)

    # 填充新值（只填充非空值）
    df_merged['AG_ATAC'] = df_merged['ATAC_Max_Delta'].fillna(
        df_merged.get('AG_ATAC', 0)
    )
    df_merged['AG_CHIP_HISTONE'] = df_merged['Histone_Max_Delta'].fillna(
        df_merged.get('AG_CHIP_HISTONE', 0)
    )
    df_merged['AG_CAGE'] = df_merged['CAGE_Max_Delta'].fillna(
        df_merged.get('AG_CAGE', 0)
    )
    df_merged['AG_PROCAP'] = df_merged['PROCAP_Max_Delta'].fillna(
        df_merged.get('AG_PROCAP', 0)
    )
    df_merged['AG_CONTACT_MAPS'] = df_merged['Contact_Max_Delta'].fillna(
        df_merged.get('AG_CONTACT_MAPS', 0)
    )

    # 添加数据来源标记
    df_merged['GodMode_Filled'] = df_merged['ATAC_Max_Delta'].notna()

    # 清理临时列
    df_merged = df_merged.drop(columns=[
        'variant_key', 'ATAC_Max_Delta', 'Histone_Max_Delta',
        'CAGE_Max_Delta', 'PROCAP_Max_Delta', 'Contact_Max_Delta'
    ])

    return df_merged


def calculate_final_scores(df):
    """重新计算综合得分（包含所有11种模态）"""
    logger.info("\n📊 重新计算综合得分...")

    # 定义所有模态列
    modality_cols = [
        'AG_RNA_SEQ', 'AG_SPLICE_SITES', 'AG_CAGE', 'AG_PROCAP',
        'AG_SPLICE_SITE_USAGE', 'AG_SPLICE_JUNCTIONS', 'AG_DNASE',
        'AG_ATAC', 'AG_CHIP_HISTONE', 'AG_CHIP_TF', 'AG_CONTACT_MAPS'
    ]

    # 确保所有列存在（缺失的设为0）
    for col in modality_cols:
        if col not in df.columns:
            df[col] = 0.0
            logger.warning(f"   ⚠️ 列 {col} 不存在，已设为0")

    # 计算最大值（绝对值最大的）
    df['Delta_Max_All'] = df[modality_cols].abs().max(axis=1)

    # 计算平均值
    df['Delta_Mean_All'] = df[modality_cols].abs().mean(axis=1)

    # 计算非零模态数
    df['Modality_Count'] = (df[modality_cols].abs() > 0).sum(axis=1)

    logger.info(f"   ✓ 11模态最大值: {df['Delta_Max_All'].describe()}")
    logger.info(f"   ✓ 平均非零模态数: {df['Modality_Count'].mean():.1f}")

    return df


def generate_summary(df):
    """生成合并摘要"""
    logger.info("\n" + "=" * 70)
    logger.info("📋 合并摘要")
    logger.info("=" * 70)

    total = len(df)
    filled = df['GodMode_Filled'].sum() if 'GodMode_Filled' in df.columns else 0

    logger.info(f"总变异数: {total}")
    logger.info(f"已补齐: {filled} ({filled/total*100:.1f}%)")
    logger.info(f"未补齐: {total - filled} ({(total-filled)/total*100:.1f}%)")

    # 各模态统计
    logger.info("\n各模态统计（补齐后）:")
    modality_cols = [
        ('AG_ATAC', 'ATAC'),
        ('AG_CHIP_HISTONE', 'ChIP-Histone'),
        ('AG_CAGE', 'CAGE'),
        ('AG_PROCAP', 'PROcap'),
        ('AG_CONTACT_MAPS', 'ContactMaps')
    ]

    for col, name in modality_cols:
        if col in df.columns:
            non_zero = (df[col].abs() > 0).sum()
            max_val = df[col].abs().max()
            logger.info(f"  {name:15s}: 非零={non_zero:4d}/{total} ({non_zero/total*100:5.1f}%), 最大值={max_val:.2f}")

    logger.info("\nTop 10 变异（按11模态最大值排序）:")
    top10 = df.nlargest(10, 'Delta_Max_All')[['hg38_chr', 'hg38_pos', 'ref', 'alt', 'Delta_Max_All']]
    for _, row in top10.iterrows():
        logger.info(f"  {row['hg38_chr']}:{row['hg38_pos']} {row['ref']}>{row['alt']}: {row['Delta_Max_All']:.2f}")


def save_results(df):
    """保存合并结果"""
    logger.info("\n💾 保存结果...")

    # 保存为CSV
    df.to_csv(MERGED_CSV, index=False)
    logger.info(f"   ✓ CSV: {MERGED_CSV}")

    # 保存为PKL
    df.to_pickle(MERGED_PKL)
    logger.info(f"   ✓ PKL: {MERGED_PKL}")

    # 生成一个干净的版本（只保留关键列）
    key_cols = [
        'hg19_chr', 'hg19_pos', 'hg38_chr', 'hg38_pos', 'ref', 'alt',
        'Status', 'Delta_Max_Core', 'Delta_Max_All', 'Delta_Mean_All',
        'AG_RNA_SEQ', 'AG_SPLICE_SITES', 'AG_CAGE', 'AG_PROCAP',
        'AG_SPLICE_SITE_USAGE', 'AG_SPLICE_JUNCTIONS', 'AG_DNASE',
        'AG_ATAC', 'AG_CHIP_HISTONE', 'AG_CHIP_TF', 'AG_CONTACT_MAPS',
        'GodMode_Filled'
    ]

    # 只保留存在的列
    key_cols = [c for c in key_cols if c in df.columns]
    df_clean = df[key_cols].copy()

    clean_csv = RESULTS_DIR / "08_VCF_AlphaGenome_Results_Merged_Clean.csv"
    df_clean.to_csv(clean_csv, index=False)
    logger.info(f"   ✓ 精简版CSV: {clean_csv}")


def main():
    """主函数"""
    logger.info("=" * 70)
    logger.info("🚀 08 Merge GodMode Results - 精准补齐5种模态")
    logger.info("=" * 70)

    # 1. 备份原始文件
    logger.info("\n1️⃣ 备份原始文件...")
    backup_original_files()

    # 2. 加载数据
    logger.info("\n2️⃣ 加载数据...")
    df_07 = load_07_results()
    df_07e = load_07e_results()

    # 3. 合并
    logger.info("\n3️⃣ 合并数据...")
    df_merged = merge_results(df_07, df_07e)

    # 4. 重新计算得分
    logger.info("\n4️⃣ 重新计算得分...")
    df_merged = calculate_final_scores(df_merged)

    # 5. 生成摘要
    generate_summary(df_merged)

    # 6. 保存结果
    save_results(df_merged)

    logger.info("\n" + "=" * 70)
    logger.info("🎉 合并完成！")
    logger.info(f"   输出: {MERGED_CSV}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()

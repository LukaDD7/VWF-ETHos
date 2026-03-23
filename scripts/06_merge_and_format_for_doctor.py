#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本 6：终极临床报告生成器 (终极集大成版)
包含特性：
1. Peak-Finder 32kb 震中取值法
2. RNA-seq 与 Splice Sites 危险度双轨解耦
3. hg19 坐标桥接映射 (规避 hg38 缺失)
4. GRCh37Location / location / Begin 三重防线兜底匹配
5. 完美保留医生原版的多层复合表头

输入：../results/07_VCF_AlphaGenome_Results.pkl, ../results/07_VCF_AlphaGenome_Results.csv
      ../data/*Clinvar*.xlsx (医生原版 Excel)
输出：../results/Final_VWF_Target_List_with_AlphaGenome.xlsx
"""

import pandas as pd
import numpy as np
import pickle
import logging
import glob
import sys
import os
from pathlib import Path

# ==================== 配置 ====================
# 基础目录配置
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"

# 指向我们最新跑全量的结果
PKL_FILE = RESULTS_DIR / "07_VCF_AlphaGenome_Results.pkl"
CSV_MAP_FILE = RESULTS_DIR / "07_VCF_AlphaGenome_Results.csv"  # 用作 hg19 -> hg38 的桥梁
OUTPUT_EXCEL = RESULTS_DIR / "Final_VWF_Target_List_with_AlphaGenome.xlsx"
TARGET_SHEET = "Clinvar_HGMD_merge"  # 锁定医生要求的第二张大表

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

# ==================== 数据结构 (必须与 07 脚本一致) ====================
from dataclasses import dataclass, field
from typing import Any

@dataclass
class InferenceResult:
    chromosome: str
    position: int
    ref: str
    alt: str
    interval_start: int
    interval_end: int
    clinvar_classification: str = "N/A"
    consequences: str = "N/A"
    rna_seq_delta_max: float | None = None
    rna_seq_delta_sum: float | None = None
    splice_site_delta_max: float | None = None
    splice_site_delta_sum: float | None = None
    delta_score: float = 0.0
    raw_outputs: Any = field(default=None)
    success: bool = False
    error_message: str | None = None

def extract_peak_values(result):
    """Peak-Finder：在 32kb 核心区提取震荡最强烈的绝对物理数值"""
    vals = {
        "RNA最大变化点_REF": None, "RNA最大变化点_ALT": None,
        "Splice最大变化点_REF": None, "Splice最大变化点_ALT": None
    }
    if not result.success or result.raw_outputs is None: return vals

    outputs = result.raw_outputs
    CENTER_IDX = 524288
    HALF_WINDOW = 16384
    start_idx, end_idx = CENTER_IDX - HALF_WINDOW, CENTER_IDX + HALF_WINDOW

    if outputs.reference.rna_seq is not None:
        ref_rna_mean = np.mean(outputs.reference.rna_seq.values[start_idx:end_idx, :], axis=1)
        alt_rna_mean = np.mean(outputs.alternate.rna_seq.values[start_idx:end_idx, :], axis=1)
        diff = np.abs(alt_rna_mean - ref_rna_mean)
        peak_idx = np.argmax(diff)
        vals["RNA最大变化点_REF"] = round(float(ref_rna_mean[peak_idx]), 3)
        vals["RNA最大变化点_ALT"] = round(float(alt_rna_mean[peak_idx]), 3)

    if outputs.reference.splice_sites is not None:
        ref_sp = outputs.reference.splice_sites.values[start_idx:end_idx, :]
        alt_sp = outputs.alternate.splice_sites.values[start_idx:end_idx, :]
        diff_sp = np.abs(alt_sp - ref_sp)
        row_idx, col_idx = np.unravel_index(np.argmax(diff_sp), diff_sp.shape)
        vals["Splice最大变化点_REF"] = round(float(ref_sp[row_idx, col_idx]), 4)
        vals["Splice最大变化点_ALT"] = round(float(alt_sp[row_idx, col_idx]), 4)

    return vals

def main():
    # 1. 自动寻找新版 Excel 表格 (在 data/ 目录中)
    target_files = list(DATA_DIR.glob("*Clinvar*.xlsx"))
    if not target_files:
        logger.error("未找到包含 'Clinvar_HGMD_merge+' 的 Excel 文件！")
        sys.exit(1)
    doctor_table = sorted(target_files, key=os.path.getmtime, reverse=True)[0]
    logger.info(f"1. 成功锁定医生原版 Excel 表格: {doctor_table}")

    # 2. 读取指定的 Sheet 并保留双层表头
    logger.info(f"2. 正在提取 '{TARGET_SHEET}' Sheet (保留双层表头)...")
    try:
        df_doctor = pd.read_excel(doctor_table, sheet_name=TARGET_SHEET, header=[0, 1])
    except Exception as e:
        logger.error(f"读取 Excel 失败！报错信息: {e}")
        sys.exit(1)

    # 3. 读取桥接文件 (CSV) 和张量文件 (PKL)
    logger.info("3. 正在加载 AlphaGenome 桥接映射表与 PKL 张量数据...")
    if not PKL_FILE.exists() or not CSV_MAP_FILE.exists():
        logger.error(f"找不到 {PKL_FILE} 或 {CSV_MAP_FILE}，请确认 07 脚本全量预测已跑完并生成了这两个文件。")
        sys.exit(1)

    df_map = pd.read_csv(CSV_MAP_FILE)
    with open(PKL_FILE, "rb") as f:
        results = pickle.load(f)

    # 4. 提取解耦的 AI 分数与绝对值，建立基于 hg19 的匹配池
    logger.info("4. 正在提取双模态解耦分数及震中绝对值...")

    # 建立 hg38 pos -> hg19 pos 的字典映射 (桥梁)
    pos_38_to_19 = dict(zip(df_map['hg38_pos'], df_map['hg19_pos']))

    ai_data = []
    for r in results:
        if r.success and r.position in pos_38_to_19:
            row_data = {
                "Match_Pos": str(pos_38_to_19[r.position]), # 转换为 hg19 坐标用于匹配
                "AI_RNA_Delta_Score": round(r.rna_seq_delta_max, 4) if r.rna_seq_delta_max else 0.0,
                "AI_Splice_Delta_Score": round(r.splice_site_delta_max, 4) if r.splice_site_delta_max else 0.0,
            }
            row_data.update(extract_peak_values(r))
            ai_data.append(row_data)

    df_ai = pd.DataFrame(ai_data)

    # 为 AI 数据戴上双层表头
    new_cols = []
    for c in df_ai.columns:
        if c == "Match_Pos": new_cols.append(("Meta", c))
        else: new_cols.append(("AlphaGenome 独家预测", c))
    df_ai.columns = pd.MultiIndex.from_tuples(new_cols)

    # 5. ================= 核心修补逻辑：三重坐标防线兜底 =================
    logger.info("5. 正在执行三重兜底策略提取医生表内的 hg19 坐标...")

    pos_clinvar = None
    pos_hgmd = None
    pos_begin = None

    # 动态扫描寻找可能的坐标列
    for col in df_doctor.columns:
        col_str = str(col[1]).strip()
        if col_str == 'GRCh37Location': pos_clinvar = col
        if col_str == 'location': pos_hgmd = col
        if col_str == 'Begin': pos_begin = col

    # 第一道防线：ClinVar 标准列
    s_pos_clinvar = pd.Series(np.nan, index=df_doctor.index)
    if pos_clinvar:
        s_pos_clinvar = df_doctor[pos_clinvar].astype(str).str.split('.').str[0].replace(['nan', 'None'], np.nan)

    # 第二道防线：HGMD 复杂字符串提取
    s_pos_hgmd = pd.Series(np.nan, index=df_doctor.index)
    if pos_hgmd:
        s_pos_hgmd = df_doctor[pos_hgmd].astype(str).str.extract(r':(\d+)-?')[0]

    # 第三道防线：最底层的 Begin 区间列
    s_pos_begin = pd.Series(np.nan, index=df_doctor.index)
    if pos_begin:
        s_pos_begin = df_doctor[pos_begin].astype(str).str.split('.').str[0].replace(['nan', 'None'], np.nan)

    # 终极融合：层层兜底补齐坐标缺失
    s_final_pos = s_pos_clinvar.fillna(s_pos_hgmd).fillna(s_pos_begin)
    df_doctor[("Meta", "Match_Pos")] = s_final_pos
    # ===============================================================

    logger.info("6. 正在基于融合坐标进行绝对精准的合并...")
    df_ai = df_ai.drop_duplicates(subset=[("Meta", "Match_Pos")])

    df_master = pd.merge(
        df_doctor, df_ai,
        on=[("Meta", "Match_Pos")],
        how='left'
    )

    # 7. 双模态综合预警排序
    logger.info("7. 清理排版并按 RNA 或 Splice 最高危险度进行降序排序...")
    # 统一使用 tuple 格式添加临时列，彻底消除 PerformanceWarning 警告
    df_master[("Meta", "Temp_Sort_Score")] = df_master[[("AlphaGenome 独家预测", "AI_RNA_Delta_Score"),
                                                        ("AlphaGenome 独家预测", "AI_Splice_Delta_Score")]].max(axis=1)

    df_master = df_master.sort_values(by=("Meta", "Temp_Sort_Score"), ascending=False)
    df_master = df_master.drop(columns=[("Meta", "Temp_Sort_Score"), ("Meta", "Match_Pos")])

    # 清除 Pandas 自动生成的 'Unnamed:' 丑陋前缀
    clean_columns = []
    for col in df_master.columns:
        if str(col[0]).startswith('Unnamed:'): clean_columns.append(("", col[1]))
        else: clean_columns.append(col)
    df_master.columns = pd.MultiIndex.from_tuples(clean_columns)

    # 8. 终极导出：绕过 Pandas 限制的底层黑客写法
    logger.info("8. 正在利用 openpyxl 引擎绕过 Pandas 限制，生成完美 Excel...")

    with pd.ExcelWriter(OUTPUT_EXCEL, engine='openpyxl') as writer:
        # 第一步：顺从 Pandas 的脾气，允许 index=True 写入内存
        df_master.to_excel(writer, index=True, sheet_name="合并结果")

        # 第二步：直接操纵 Excel 底层，把强行生成的丑陋第一列（Index列）咔嚓掉！
        worksheet = writer.sheets['合并结果']
        worksheet.delete_cols(1)

    logger.info(f"🎉 大功告成！双层表头保留完美，报表已生成：{OUTPUT_EXCEL}")

if __name__ == "__main__":
    main()

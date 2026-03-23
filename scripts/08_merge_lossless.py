#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本 8：绝对无损临床报告合并器 (基于 hg19 搭桥匹配，解决 Liftover 漂移问题)

输入：../results/07_VCF_AlphaGenome_Results.csv (搭桥文件)
      ../results/07_VCF_AlphaGenome_Results.pkl
      ../data/*Clinvar*.xlsx (医生原版 Excel)
输出：../results/Final_VWF_Lossless_Merge.xlsx
"""

import pandas as pd
import numpy as np
import pickle
import logging
import glob
import sys
import os
import re
from pathlib import Path

# ==================== 配置 ====================
# 基础目录配置
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"

CSV_BRIDGE = RESULTS_DIR / "07_VCF_AlphaGenome_Results.csv"  # 必须用到它来搭桥
PKL_FILE = RESULTS_DIR / "07_VCF_AlphaGenome_Results.pkl"
OUTPUT_EXCEL = RESULTS_DIR / "Final_VWF_Lossless_Merge.xlsx"
TARGET_SHEET = "Clinvar_HGMD_merge"  # 锁定 2577 行的完整全集

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

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

def parse_consequence(info_str):
    parsed = {"AI_突变后果": "", "AI_转录本": "", "AI_HGVSc": "", "AI_HGVSp": ""}
    if pd.isna(info_str) or not isinstance(info_str, str): return parsed
    parts = str(info_str).split('|')
    for p in parts:
        if any(x in p for x in ["variant", "stop", "frameshift", "missense", "splice"]):
            parsed["AI_突变后果"] = p
            break
    for p in parts:
        if p.startswith("NM_") or p.startswith("ENST"):
            parsed["AI_转录本"] = p
            break
    for p in parts:
        if p.startswith("c.") or p.startswith("n."): parsed["AI_HGVSc"] = p
        elif p.startswith("p."): parsed["AI_HGVSp"] = p
    return parsed

def extract_peak_values(result):
    vals = {
        "RNA_REF基准": "", "RNA_ALT突变": "",
        "Splice_REF基准": "", "Splice_ALT突变": ""
    }
    if not result.success or result.raw_outputs is None: return vals
    outputs = result.raw_outputs
    CENTER_IDX = 524288
    HALF_WINDOW = 16384
    start_idx, end_idx = CENTER_IDX - HALF_WINDOW, CENTER_IDX + HALF_WINDOW

    if outputs.reference.rna_seq is not None:
        ref_rna = np.mean(outputs.reference.rna_seq.values[start_idx:end_idx, :], axis=1)
        alt_rna = np.mean(outputs.alternate.rna_seq.values[start_idx:end_idx, :], axis=1)
        peak_idx = np.argmax(np.abs(alt_rna - ref_rna))
        vals["RNA_REF基准"] = round(float(ref_rna[peak_idx]), 3)
        vals["RNA_ALT突变"] = round(float(alt_rna[peak_idx]), 3)

    if outputs.reference.splice_sites is not None:
        ref_sp = outputs.reference.splice_sites.values[start_idx:end_idx, :]
        alt_sp = outputs.alternate.splice_sites.values[start_idx:end_idx, :]
        diff_sp = np.abs(alt_sp - ref_sp)
        row_idx, col_idx = np.unravel_index(np.argmax(diff_sp), diff_sp.shape)
        vals["Splice_REF基准"] = round(float(ref_sp[row_idx, col_idx]), 4)
        vals["Splice_ALT突变"] = round(float(alt_sp[row_idx, col_idx]), 4)
    return vals

def clean_coord(coord_str):
    """提取坐标里的第一个纯数字"""
    match = re.search(r'\d+', str(coord_str))
    return match.group() if match else None

def main():
    # 1. 加载 07 脚本生成的 CSV，搭建 hg19 到 hg38 的桥梁
    logger.info(f"1. 加载搭桥文件: {CSV_BRIDGE}")
    try:
        df_csv = pd.read_csv(CSV_BRIDGE)
        hg19_to_hg38 = {}
        for _, row in df_csv.iterrows():
            hg19_pos = clean_coord(row['hg19_pos'])
            hg38_pos = clean_coord(row['hg38_pos'])
            if hg19_pos and hg38_pos:
                hg19_to_hg38[hg19_pos] = hg38_pos
    except Exception as e:
        logger.error(f"读取 {CSV_BRIDGE} 失败，请确认文件存在: {e}")
        sys.exit(1)

    # 2. 加载 PKL 数据
    logger.info(f"2. 加载 AI 预测 PKL: {PKL_FILE}")
    with open(PKL_FILE, "rb") as f:
        results = pickle.load(f)
    ai_dict = {str(r.position): r for r in results if r.success}

    # 3. 寻找 Excel (在 data/ 目录中)
    all_xlsx = list(DATA_DIR.glob("*Clinvar*.xlsx"))
    doctor_table = None
    for f in sorted(all_xlsx, key=os.path.getmtime, reverse=True):
        try:
            xl = pd.ExcelFile(f)
            if TARGET_SHEET in xl.sheet_names:
                doctor_table = f
                break
        except: continue

    if not doctor_table:
        logger.error(f"未找到包含 '{TARGET_SHEET}' Sheet 的 Excel 文件！")
        sys.exit(1)

    logger.info(f"3. 锁定原表: {doctor_table}")
    df_raw = pd.read_excel(doctor_table, sheet_name=TARGET_SHEET, header=None)

    # 寻找 hg19 的坐标列 (GRCh37Location)
    sub_headers = df_raw.iloc[1].fillna('').astype(str).str.lower()
    pos_idx, cons_idx = None, None
    for i, val in enumerate(sub_headers):
        if 'grch37location' in val: pos_idx = i
        if 'molecular.consequence' in val or 'consequence' in val: cons_idx = i

    if pos_idx is None:
        logger.error("在原表第二行找不到 GRCh37Location，无法进行 hg19 匹配！")
        sys.exit(1)

    new_cols = [
        "AI_hg38_Position (核对用)", "AlphaGenome_Score",
        "AI_突变后果", "AI_转录本", "AI_HGVSc",
        "RNA_REF基准", "RNA_ALT突变", "Splice_REF基准", "Splice_ALT突变"
    ]
    ai_matrix = [[""] * len(new_cols) for _ in range(len(df_raw))]
    ai_matrix[0] = ["AlphaGenome 预测 (100%原序)"] + [""] * (len(new_cols) - 1)
    ai_matrix[1] = new_cols

    logger.info("4. 开始逐行扫描原表，使用 hg19 搭桥匹配...")
    matched_count = 0
    for row_idx in range(2, len(df_raw)):
        raw_hg19 = str(df_raw.iat[row_idx, pos_idx])
        hg19_pos = clean_coord(raw_hg19)

        # 搭桥逻辑：hg19 -> hg38 -> PKL
        if hg19_pos in hg19_to_hg38:
            hg38_pos = hg19_to_hg38[hg19_pos]
            if hg38_pos in ai_dict:
                matched_count += 1
                r = ai_dict[hg38_pos]

                cons_text = str(df_raw.iat[row_idx, cons_idx]) if cons_idx is not None else ""
                parsed = parse_consequence(cons_text)
                peaks = extract_peak_values(r)

                ai_matrix[row_idx] = [
                    r.position, round(r.delta_score, 4),
                    parsed.get("AI_突变后果", ""), parsed.get("AI_转录本", ""), parsed.get("AI_HGVSc", ""),
                    peaks["RNA_REF基准"], peaks["RNA_ALT突变"], peaks["Splice_REF基准"], peaks["Splice_ALT突变"]
                ]

    logger.info(f"匹配完成！原表正文共有 {len(df_raw)-2} 行，成功匹配上了 {matched_count} 行。")

    df_ai = pd.DataFrame(ai_matrix)
    df_final = pd.concat([df_raw, df_ai], axis=1)
    df_final.to_excel(OUTPUT_EXCEL, index=False, header=False)
    logger.info(f"🎉 文件已保存为 {OUTPUT_EXCEL}。")

if __name__ == "__main__":
    main()

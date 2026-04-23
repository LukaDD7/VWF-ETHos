#!/usr/bin/env python3
"""
VWF Variant Table Normalizer
==========================

统一处理各种格式的 VWF 变异表格，自动检测格式并转换为标准输入格式。

支持表格格式:
1. 3-letter格式: "Tyr1258Cys", "Pro1266Gln"
2. 1-letter格式: "Y1258C", "P1266Q"
3. 分列格式: AA_Position | WT_AA | Mut_AA

输出格式:
- AA_Position: int
- WT_AA_1: str (1-letter)
- Mut_AA_1: str (1-letter)
- Domain: str (optional)
- Original: str (原始格式)

使用示例:
    normalizer = TableNormalizer("path/to/2B_variants.xlsx")
    normalized_df = normalizer.normalize()
    normalizer.validate()  # 检查错误
    normalizer.save_csv("output.csv")

    # 直接转换为 AF3 JSON
    normalizer.to_af3_json(output_dir="af3_batch")

Author: Claude Code
Date: 2026-04-23
"""

import re
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# 常量
# ============================================================================

THREE_TO_ONE = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Asx': 'B',
    'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Glx': 'Z', 'Gly': 'G',
    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M',
    'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
    'Tyr': 'Y', 'Val': 'V', 'Sec': 'U', 'Pyl': 'O'
}

ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

# 标准列名
STANDARD_COLUMNS = ['AA_Position', 'WT_AA_1', 'Mut_AA_1']
OPTIONAL_COLUMNS = ['Domain', 'Original', 'Type2_Subtype', 'Source']


# ============================================================================
# 表格格式检测与解析
# ============================================================================

class FormatDetector:
    """检测表格格式"""

    @staticmethod
    def detect_3letter_format(value: str) -> bool:
        """检测是否是3-letter格式: Tyr1258Cys"""
        if pd.isna(value):
            return False
        return bool(re.match(r'^[A-Za-z]{3}\d+[A-Za-z]{3}$', str(value)))

    @staticmethod
    def detect_1letter_format(value: str) -> bool:
        """检测是否是1-letter格式: Y1258C"""
        if pd.isna(value):
            return False
        return bool(re.match(r'^[A-Za-z]\d+[A-Za-z*]$', str(value)))

    @staticmethod
    def parse_3letter(aa_change: str) -> Optional[Tuple[str, int, str]]:
        """
        解析3-letter格式: Tyr1258Cys → (T, 1258, Y)
        返回 (wt_1letter, position, mut_1letter)
        """
        match = re.match(r'([A-Za-z]{3})(\d+)([A-Za-z]{3})', str(aa_change))
        if not match:
            return None

        wt_3letter = match.group(1)
        position = int(match.group(2))
        mut_3letter = match.group(3)

        wt_1letter = THREE_TO_ONE.get(wt_3letter.capitalize())
        mut_1letter = THREE_TO_ONE.get(mut_3letter.capitalize())

        if not wt_1letter or not mut_1letter:
            return None

        return wt_1letter, position, mut_1letter

    @staticmethod
    def parse_1letter(aa_change: str) -> Optional[Tuple[str, int, str]]:
        """
        解析1-letter格式: Y1258C → (Y, 1258, C)
        返回 (wt_1letter, position, mut_1letter)
        """
        match = re.match(r'([A-Za-z])(\d+)([A-Za-z*])', str(aa_change))
        if not match:
            return None

        return match.group(1).upper(), int(match.group(2)), match.group(3).upper()


# ============================================================================
# 主类
# ============================================================================

class TableNormalizer:
    """
    VWF 变异表格标准化器

    自动检测输入表格格式，统一转换为标准格式。
    """

    def __init__(self, excel_path: str, sheet: int = 0):
        """
        Args:
            excel_path: Excel 文件路径
            sheet: 工作表索引或名称
        """
        self.excel_path = Path(excel_path)
        self.df_raw = pd.read_excel(excel_path, sheet_name=sheet)
        self.df_normalized: Optional[pd.DataFrame] = None
        self.format_detected: Optional[str] = None
        self.errors: List[str] = []

        # 检测并跳过标题行
        self._skip_title_rows()

    def _skip_title_rows(self):
        """
        检测并跳过表格标题行
        有些 Excel 文件第一行是标题如 "Table S3: ...", 需要跳过
        可能还有列名行如 "Amino acid\nchange"
        """
        skipped = 0
        while len(self.df_raw) > 0:
            first_value = str(self.df_raw.iloc[0, 0]).lower()
            # 跳过标题行
            if 'table s' in first_value or 'vwd' in first_value or 'type 2' in first_value:
                self.df_raw = self.df_raw.iloc[1:].reset_index(drop=True)
                skipped += 1
                continue
            # 跳过列名行 (如 "Amino acid\nchange")
            if 'amino' in first_value or 'change' in first_value:
                self.df_raw = self.df_raw.iloc[1:].reset_index(drop=True)
                skipped += 1
                continue
            break

        if skipped > 0:
            logger.info(f"跳过 {skipped} 个标题/列名行")

    def detect_format(self) -> str:
        """
        自动检测表格格式

        Returns:
            format_name: '3letter', '1letter', 'separated', 'unknown'
        """
        # 检查第一行是否是header
        first_col = self.df_raw.columns[0]
        sample_values = self.df_raw.iloc[:5, 0].dropna()

        if len(sample_values) == 0:
            return 'unknown'

        # 尝试3-letter格式
        if all(FormatDetector.detect_3letter_format(v) for v in sample_values):
            return '3letter'

        # 尝试1-letter格式
        if all(FormatDetector.detect_1letter_format(v) for v in sample_values):
            return '1letter'

        # 检查是否已经是分列格式
        cols_lower = [str(c).lower() for c in self.df_raw.columns]
        if any('position' in c or 'pos' in c for c in cols_lower):
            return 'separated'

        return 'unknown'

    def normalize(self) -> pd.DataFrame:
        """
        标准化表格

        Returns:
            标准格式 DataFrame，包含列: AA_Position, WT_AA_1, Mut_AA_1, Domain, Original
        """
        self.format_detected = self.detect_format()
        logger.info(f"检测到表格格式: {self.format_detected}")

        if self.format_detected == '3letter':
            self.df_normalized = self._parse_3letter()
        elif self.format_detected == '1letter':
            self.df_normalized = self._parse_1letter()
        elif self.format_detected == 'separated':
            self.df_normalized = self._parse_separated()
        else:
            raise ValueError(f"无法识别的表格格式: {self.format_detected}")

        # 过滤同义突变
        initial_count = len(self.df_normalized)
        self.df_normalized = self.df_normalized[
            self.df_normalized['WT_AA_1'] != self.df_normalized['Mut_AA_1']
        ]
        filtered_count = initial_count - len(self.df_normalized)
        if filtered_count > 0:
            logger.warning(f"过滤了 {filtered_count} 个同义突变")

        return self.df_normalized

    def _parse_3letter(self) -> pd.DataFrame:
        """解析3-letter格式"""
        records = []
        first_col = self.df_raw.columns[0]  # AA change 列

        # 检测 Domain 列位置
        domain_col_idx = None
        for i, col in enumerate(self.df_raw.columns):
            if str(col).lower() in ['domain', 'd3', 'a1']:
                domain_col_idx = i
                break

        for _, row in self.df_raw.iterrows():
            aa_change = row[first_col]
            if pd.isna(aa_change):
                continue

            result = FormatDetector.parse_3letter(aa_change)
            if result:
                wt, pos, mut = result
                domain = row[domain_col_idx] if domain_col_idx is not None else None
                records.append({
                    'AA_Position': pos,
                    'WT_AA_1': wt,
                    'Mut_AA_1': mut,
                    'Domain': domain if pd.notna(domain) else None,
                    'Original': str(aa_change)
                })

        df = pd.DataFrame(records)
        logger.info(f"解析了 {len(df)} 条变异记录")
        return df

    def _parse_1letter(self) -> pd.DataFrame:
        """解析1-letter格式"""
        records = []
        first_col = self.df_raw.columns[0]

        domain_col_idx = None
        for i, col in enumerate(self.df_raw.columns):
            if str(col).lower() in ['domain', 'd3', 'a1']:
                domain_col_idx = i
                break

        for _, row in self.df_raw.iterrows():
            aa_change = row[first_col]
            if pd.isna(aa_change):
                continue

            result = FormatDetector.parse_1letter(aa_change)
            if result:
                wt, pos, mut = result
                domain = row[domain_col_idx] if domain_col_idx is not None else None
                records.append({
                    'AA_Position': pos,
                    'WT_AA_1': wt,
                    'Mut_AA_1': mut,
                    'Domain': domain if pd.notna(domain) else None,
                    'Original': str(aa_change)
                })

        df = pd.DataFrame(records)
        logger.info(f"解析了 {len(df)} 条变异记录")
        return df

    def _parse_separated(self) -> pd.DataFrame:
        """解析分列格式 (AA_Position, WT_AA, Mut_AA)"""
        # 查找标准列
        cols = {c: c.lower() for c in self.df_raw.columns}

        pos_col = next((c for c in cols if 'position' in cols[c] or 'pos' in cols[c]), None)
        wt_col = next((c for c in cols if 'wt' in cols[c] or 'wild' in cols[c]), None)
        mut_col = next((c for c in cols if 'mut' in cols[c] or 'alt' in cols[c]), None)

        if not pos_col:
            raise ValueError("缺少位置列 (Position/Pos)")

        # 如果没有WT/Mut列，尝试从name列解析
        if not wt_col or not mut_col:
            name_col = next((c for c in cols if 'name' in cols[c] or 'change' in cols[c]), self.df_raw.columns[0])
            records = []
            for _, row in self.df_raw.iterrows():
                pos = row[pos_col]
                if pd.isna(pos):
                    continue
                aa_change = row[name_col]
                if pd.isna(aa_change):
                    continue
                result = FormatDetector.parse_3letter(aa_change) or FormatDetector.parse_1letter(aa_change)
                if result:
                    wt, pos_num, mut = result
                    records.append({
                        'AA_Position': int(float(pos)),
                        'WT_AA_1': wt,
                        'Mut_AA_1': mut,
                        'Domain': None,
                        'Original': str(aa_change)
                    })
        else:
            # 直接使用分列
            records = []
            for _, row in self.df_raw.iterrows():
                if pd.isna(row[pos_col]):
                    continue
                records.append({
                    'AA_Position': int(float(row[pos_col])),
                    'WT_AA_1': str(row[wt_col]).strip().upper()[0] if pd.notna(row[wt_col]) else None,
                    'Mut_AA_1': str(row[mut_col]).strip().upper()[0] if pd.notna(row[mut_col]) else None,
                    'Domain': None,
                    'Original': f"{row[wt_col]}{row[pos_col]}{row[mut_col]}" if pd.notna(row[wt_col]) and pd.notna(row[mut_col]) else None
                })

        df = pd.DataFrame(records)
        logger.info(f"解析了 {len(df)} 条变异记录")
        return df

    def validate(self) -> List[str]:
        """
        验证标准化后的表格

        Returns:
            错误列表，空列表表示验证通过
        """
        self.errors = []

        if self.df_normalized is None:
            self.errors.append("请先调用 normalize()")
            return self.errors

        # 检查必需列
        for col in STANDARD_COLUMNS:
            if col not in self.df_normalized.columns:
                self.errors.append(f"缺少必需列: {col}")

        # 检查AA_Position范围
        positions = self.df_normalized['AA_Position']
        if positions.min() < 1 or positions.max() > 3000:
            self.errors.append(f"AA_Position 超出范围: {positions.min()}-{positions.max()}")

        # 检查氨基酸有效性
        valid_aa = set(THREE_TO_ONE.values())
        invalid_wt = ~self.df_normalized['WT_AA_1'].isin(valid_aa)
        invalid_mut = ~self.df_normalized['Mut_AA_1'].isin(valid_aa)

        if invalid_wt.any():
            bad = self.df_normalized[invalid_wt]['WT_AA_1'].unique()
            self.errors.append(f"无效的 WT_AA: {bad}")

        if invalid_mut.any():
            bad = self.df_normalized[invalid_mut]['Mut_AA_1'].unique()
            self.errors.append(f"无效的 Mut_AA: {bad}")

        # 检查重复
        dupes = self.df_normalized.duplicated(subset=['AA_Position', 'WT_AA_1', 'Mut_AA_1'])
        if dupes.any():
            n_dupes = dupes.sum()
            self.errors.append(f"发现 {n_dupes} 个重复变异")

        if self.errors:
            logger.warning(f"验证发现 {len(self.errors)} 个问题")
            for err in self.errors:
                logger.warning(f"  - {err}")
        else:
            logger.info("验证通过 ✓")

        return self.errors

    def save_csv(self, output_path: str):
        """保存标准化后的表格"""
        if self.df_normalized is None:
            raise ValueError("请先调用 normalize()")
        self.df_normalized.to_csv(output_path, index=False)
        logger.info(f"已保存到 {output_path}")

    def to_af3_json(
        self,
        wt_fasta: str,
        output_path: Optional[str] = None,
        chunk_size: int = 30,
        exclude_existing: Optional[List[str]] = None
    ) -> List[Dict]:
        """
        转换为 AlphaFold3 JSON 格式

        Args:
            wt_fasta: WT 序列 FASTA 文件路径
            output_path: 输出 JSON 文件路径 (可选)
            chunk_size: 每个 JSON 文件的任务数
            exclude_existing: 排除已存在的任务名 (如 L1460F, A1461V)

        Returns:
            AF3 batch jobs 列表
        """
        if self.df_normalized is None:
            raise ValueError("请先调用 normalize()")

        from Bio import SeqIO

        # 加载 WT 序列
        wt_record = SeqIO.read(wt_fasta, "fasta")
        wt_sequence = str(wt_record.seq)

        # 构建任务
        jobs = []
        for _, row in self.df_normalized.iterrows():
            pos = int(row['AA_Position'])
            wt_aa = row['WT_AA_1']
            mut_aa = row['Mut_AA_1']
            job_name = f"VWF_{wt_aa}{pos}{mut_aa}"

            # 排除已存在的任务
            if exclude_existing and job_name.replace("VWF_", "") in exclude_existing:
                continue

            # 生成突变序列
            mut_sequence = wt_sequence[:pos-1] + mut_aa + wt_sequence[pos:]

            jobs.append({
                "name": job_name,
                "sequences": [{
                    "proteinChain": {
                        "sequence": mut_sequence,
                        "count": 1
                    }
                }]
            })

        logger.info(f"生成了 {len(jobs)} 个 AF3 任务")

        # 保存到文件
        if output_path:
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(jobs, f, indent=2, ensure_ascii=False)
            logger.info(f"已保存到 {output_path}")

        return jobs


# ============================================================================
# 主函数
# ============================================================================

if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description="VWF 变异表格标准化工具")
    parser.add_argument("input", help="输入 Excel 文件路径")
    parser.add_argument("-o", "--output", help="输出 CSV 路径")
    parser.add_argument("--af3-json", help="同时生成 AF3 JSON 文件")
    parser.add_argument("--wt-fasta", help="WT 序列 FASTA 文件 (用于 AF3 JSON)")
    parser.add_argument("--exclude", nargs="+", help="排除已存在的任务名")
    parser.add_argument("--chunk-size", type=int, default=30, help="AF3 每文件任务数")
    parser.add_argument("--skip-validation", action="store_true", help="跳过验证")

    args = parser.parse_args()

    # 标准化
    normalizer = TableNormalizer(args.input)
    df = normalizer.normalize()

    # 验证
    if not args.skip_validation:
        errors = normalizer.validate()
        if errors:
            sys.exit(1)

    # 输出
    print(f"\n标准化结果 ({len(df)} 条记录):")
    print(df.head(10).to_string())

    if args.output:
        normalizer.save_csv(args.output)

    if args.af3_json:
        if not args.wt_fasta:
            print("错误: --af3-json 需要 --wt-fasta")
            sys.exit(1)
        jobs = normalizer.to_af3_json(
            wt_fasta=args.wt_fasta,
            output_path=args.af3_json,
            exclude_existing=args.exclude,
            chunk_size=args.chunk_size
        )
        print(f"\nAF3 JSON 已保存到 {args.af3_json}")

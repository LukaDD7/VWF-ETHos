#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1: Smart Filter & Parsing
Proteo-Structure Pipeline for VWF (Uniprot: P04275)

功能:
1. 从Excel过滤missense variants
2. 排除splice disruptors (|Splice_REF - Splice_ALT| > 0.3)
3. 提取氨基酸变化信息 (1-letter和3-letter)
4. 下载野生型FASTA
5. 生成突变型FASTA
"""

import pandas as pd
import numpy as np
import re
import requests
import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from datetime import datetime


@dataclass
class AminoAcidChange:
    """氨基酸变化数据类"""
    wt_aa_1: str          # 野生型1字母
    mut_aa_1: str         # 突变型1字母
    wt_aa_3: str          # 野生型3字母
    mut_aa_3: str         # 突变型3字母
    position: int         # 氨基酸位置
    raw_str: str          # 原始字符串


class Phase1Filter:
    """Phase 1 过滤和解析器"""

    # 氨基酸1字母 <-> 3字母转换表
    AA_1TO3 = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
        'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
        'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
        'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
        '*': 'Ter'  # 终止密码子
    }

    AA_3TO1 = {v: k for k, v in AA_1TO3.items()}

    # 支持的3字母缩写
    AA_3LETTER_PATTERN = '|'.join(AA_3TO1.keys())

    def __init__(self, excel_path: str, output_dir: str):
        self.excel_path = excel_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 创建子目录
        self.structures_dir = self.output_dir.parent / 'structures'
        self.wt_dir = self.structures_dir / 'wt'
        self.mutant_dir = self.structures_dir / 'mutant'
        self.wt_dir.mkdir(parents=True, exist_ok=True)
        self.mutant_dir.mkdir(parents=True, exist_ok=True)

        self.df = None
        self.filtered_df = None
        self.wt_sequence = None

    def load_excel(self, sheet_name: str = 'Sheet1') -> pd.DataFrame:
        """加载Excel文件"""
        print(f"[1/6] 正在加载Excel: {self.excel_path}")
        self.df = pd.read_excel(self.excel_path, sheet_name=sheet_name)
        print(f"      总计 {len(self.df)} 行, {len(self.df.columns)} 列")
        return self.df

    def filter_missense_variants(self) -> pd.DataFrame:
        """过滤missense variants"""
        print("\n[2/6] 过滤Missense Variants...")

        # 检查列是否存在
        mol_cons_col = 'Molecular.consequence'
        ai_cons_col = 'AI_突变后果'

        mask = pd.Series([False] * len(self.df))

        if mol_cons_col in self.df.columns:
            mol_mask = self.df[mol_cons_col].astype(str).str.contains('missense variant', na=False, case=False)
            print(f"      {mol_cons_col}: {mol_mask.sum()} 个missense")
            mask = mask | mol_mask

        if ai_cons_col in self.df.columns:
            ai_mask = self.df[ai_cons_col].astype(str).str.contains('missense variant', na=False, case=False)
            print(f"      {ai_cons_col}: {ai_mask.sum()} 个missense")
            mask = mask | ai_mask

        self.filtered_df = self.df[mask].copy()
        print(f"      => Missense过滤后: {len(self.filtered_df)} 个变异")
        return self.filtered_df

    def filter_splice_disruptors(self, threshold: float = 0.3) -> pd.DataFrame:
        """排除splice disruptors (|Splice_REF - Splice_ALT| > threshold)"""
        print(f"\n[3/6] 排除Splice Disruptors (阈值: |REF - ALT| > {threshold})...")

        splice_ref_col = 'Splice_REF正常概率'
        splice_alt_col = 'Splice_ALT突变概率'

        if splice_ref_col not in self.filtered_df.columns or splice_alt_col not in self.filtered_df.columns:
            print(f"      警告: 未找到splice列，跳过此过滤")
            return self.filtered_df

        # 转换为数值
        self.filtered_df[splice_ref_col] = pd.to_numeric(self.filtered_df[splice_ref_col], errors='coerce')
        self.filtered_df[splice_alt_col] = pd.to_numeric(self.filtered_df[splice_alt_col], errors='coerce')

        # 计算delta
        splice_delta = (self.filtered_df[splice_ref_col] - self.filtered_df[splice_alt_col]).abs()

        # 保留splice delta <= threshold的
        keep_mask = (splice_delta <= threshold) | (splice_delta.isna())
        excluded = (~keep_mask).sum()

        self.filtered_df = self.filtered_df[keep_mask].copy()
        self.filtered_df['Splice_Delta'] = splice_delta[keep_mask]

        print(f"      排除 {excluded} 个splice disruptors")
        print(f"      => Splice过滤后: {len(self.filtered_df)} 个变异")
        return self.filtered_df

    def parse_amino_acid_change(self, row) -> Optional[AminoAcidChange]:
        """解析氨基酸变化，尝试多种列"""
        # 优先尝试的列
        protein_change = str(row.get('Protein.change', ''))
        name = str(row.get('Name', ''))
        p_col = str(row.get('p.', ''))

        # 尝试1: Protein.change (如 G1531D)
        if protein_change and protein_change not in ['nan', 'None', '']:
            aa_change = self._parse_1letter(protein_change)
            if aa_change:
                return aa_change

        # 尝试2: Name列 (如 p.Gly1531Asp)
        if name and name not in ['nan', 'None', '']:
            aa_change = self._parse_name_column(name)
            if aa_change:
                return aa_change

        # 尝试3: p.列 (如 Val1409Phe)
        if p_col and p_col not in ['nan', 'None', '']:
            aa_change = self._parse_3letter(p_col)
            if aa_change:
                return aa_change

        return None

    def _parse_1letter(self, s: str) -> Optional[AminoAcidChange]:
        """解析1字母格式: G1531D, T1538K, P1551L"""
        # 格式: [A-Z][0-9]+[A-Z*] 或 [A-Z][0-9]+fs (frameshift)
        pattern = r'^([A-Z])(\d+)([A-Z*])$'
        match = re.match(pattern, s)

        if match:
            wt_aa_1 = match.group(1)
            position = int(match.group(2))
            mut_aa_1 = match.group(3)

            wt_aa_3 = self.AA_1TO3.get(wt_aa_1, wt_aa_1)
            mut_aa_3 = self.AA_1TO3.get(mut_aa_1, mut_aa_1)

            return AminoAcidChange(
                wt_aa_1=wt_aa_1,
                mut_aa_1=mut_aa_1,
                wt_aa_3=wt_aa_3,
                mut_aa_3=mut_aa_3,
                position=position,
                raw_str=s
            )
        return None

    def _parse_3letter(self, s: str) -> Optional[AminoAcidChange]:
        """解析3字母格式: Val1409Phe, Gly1531Asp"""
        # 清理字符串
        s = s.strip()
        s = re.sub(r'\s+', '', s)  # 移除空格

        # 格式: (Ala|Arg|...)(\d+)(Ala|Arg|...)
        pattern = rf'^({self.AA_3LETTER_PATTERN})(\d+)({self.AA_3LETTER_PATTERN})$'
        match = re.match(pattern, s, re.IGNORECASE)

        if match:
            wt_aa_3 = match.group(1).title()
            position = int(match.group(2))
            mut_aa_3 = match.group(3).title()

            wt_aa_1 = self.AA_3TO1.get(wt_aa_3, 'X')
            mut_aa_1 = self.AA_3TO1.get(mut_aa_3, 'X')

            return AminoAcidChange(
                wt_aa_1=wt_aa_1,
                mut_aa_1=mut_aa_1,
                wt_aa_3=wt_aa_3,
                mut_aa_3=mut_aa_3,
                position=position,
                raw_str=s
            )
        return None

    def _parse_name_column(self, s: str) -> Optional[AminoAcidChange]:
        """解析Name列中的p.信息"""
        # 格式: NM_000552.5(VWF):c.4592G>A (p.Gly1531Asp)
        # 或: (p.Gly1531Asp)
        pattern = r'\(p\.([A-Za-z]+\d+[A-Za-z]+)\)'
        match = re.search(pattern, s)

        if match:
            aa_str = match.group(1)
            return self._parse_3letter(aa_str)
        return None

    def extract_mutations(self) -> pd.DataFrame:
        """提取所有变异信息"""
        print("\n[4/6] 解析氨基酸变化...")

        mutations = []
        parsed_count = 0
        failed_count = 0

        for idx, row in self.filtered_df.iterrows():
            aa_change = self.parse_amino_acid_change(row)

            if aa_change:
                mut_data = {
                    'Original_Index': idx,
                    'Chromosome': row.get('Chromosome', ''),
                    'Position': row.get('Position', ''),
                    'Ref': row.get('Ref', ''),
                    'Alt': row.get('Alt', ''),
                    'WT_AA_1': aa_change.wt_aa_1,
                    'Mut_AA_1': aa_change.mut_aa_1,
                    'WT_AA_3': aa_change.wt_aa_3,
                    'Mut_AA_3': aa_change.mut_aa_3,
                    'AA_Position': aa_change.position,
                    'Raw_Change': aa_change.raw_str,
                    'Splice_Delta': row.get('Splice_Delta', np.nan),
                    'Source_Column': self._get_source_column(row),
                }
                # 添加原始行的其他列
                for k, v in row.items():
                    if k not in mut_data:
                        mut_data[k] = v
                mutations.append(mut_data)
                parsed_count += 1
            else:
                failed_count += 1
                if failed_count <= 3:
                    print(f"      警告: 无法解析行 {idx}: Protein.change={row.get('Protein.change')}, Name={row.get('Name')}, p.={row.get('p.')}")

        print(f"      成功解析: {parsed_count}, 失败: {failed_count}")

        self.filtered_df = pd.DataFrame(mutations)
        print(f"      => 最终可用变异: {len(self.filtered_df)} 个")
        return self.filtered_df

    def _get_source_column(self, row) -> str:
        """判断数据来自哪一列"""
        protein_change = str(row.get('Protein.change', ''))
        name = str(row.get('Name', ''))
        p_col = str(row.get('p.', ''))

        if protein_change and protein_change not in ['nan', 'None', '']:
            return 'Protein.change'
        elif 'p.' in name:
            return 'Name'
        elif p_col and p_col not in ['nan', 'None', '']:
            return 'p.'
        return 'unknown'

    def download_wt_fasta(self, uniprot_id: str = 'P04275') -> str:
        """从Uniprot下载野生型FASTA"""
        print(f"\n[5/6] 下载野生型FASTA (Uniprot: {uniprot_id})...")

        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        output_path = self.wt_dir / f"VWF_{uniprot_id}_WT.fasta"

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            with open(output_path, 'w') as f:
                f.write(response.text)

            # 解析序列
            lines = response.text.strip().split('\n')
            header = lines[0]
            self.wt_sequence = ''.join(lines[1:])

            print(f"      成功下载: {output_path}")
            print(f"      序列长度: {len(self.wt_sequence)} aa")
            print(f"      序列标题: {header[:80]}...")

            return str(output_path)

        except Exception as e:
            print(f"      错误: 下载失败 {e}")
            raise

    def generate_mutant_fasta(self) -> List[str]:
        """生成所有突变型FASTA文件"""
        print(f"\n[6/6] 生成突变型FASTA文件...")

        if self.wt_sequence is None:
            raise ValueError("请先下载野生型FASTA")

        fasta_files = []
        success_count = 0
        fail_count = 0

        for idx, row in self.filtered_df.iterrows():
            try:
                position = int(row['AA_Position']) - 1  # 0-based index
                wt_aa = row['WT_AA_1']
                mut_aa = row['Mut_AA_1']
                mut_pos = row['AA_Position']

                # 验证野生型氨基酸
                if position < 0 or position >= len(self.wt_sequence):
                    print(f"      警告: 位置 {position+1} 超出序列范围")
                    fail_count += 1
                    continue

                actual_wt = self.wt_sequence[position]
                if actual_wt != wt_aa:
                    print(f"      警告: 位置 {position+1} 期望WT={wt_aa}, 实际={actual_wt}")
                    fail_count += 1
                    continue

                # 生成突变序列
                mut_sequence = self.wt_sequence[:position] + mut_aa + self.wt_sequence[position+1:]

                # 创建FASTA
                chrom = row.get('Chromosome', 'chrNA')
                pos = row.get('Position', 'NA')
                filename = f"VWF_{wt_aa}{mut_pos}{mut_aa}_{chrom}_{pos}.fasta"
                filepath = self.mutant_dir / filename

                with open(filepath, 'w') as f:
                    f.write(f">VWF|{wt_aa}{mut_pos}{mut_aa}|{chrom}:{pos}|WT={wt_aa}|MUT={mut_aa}\n")
                    # 每行60个字符
                    for i in range(0, len(mut_sequence), 60):
                        f.write(mut_sequence[i:i+60] + '\n')

                fasta_files.append(str(filepath))
                success_count += 1

            except Exception as e:
                print(f"      错误: 处理行 {idx} 时出错: {e}")
                fail_count += 1

        print(f"      成功生成: {success_count} 个FASTA文件")
        print(f"      失败: {fail_count}")
        print(f"      输出目录: {self.mutant_dir}")

        return fasta_files

    def save_results(self) -> Dict[str, str]:
        """保存过滤结果"""
        output_csv = self.output_dir / 'phase1_filtered_mutations.csv'
        summary_txt = self.output_dir / 'phase1_summary.txt'

        # 保存CSV (只保留关键列)
        key_cols = ['Chromosome', 'Position', 'Ref', 'Alt',
                   'WT_AA_1', 'Mut_AA_1', 'WT_AA_3', 'Mut_AA_3',
                   'AA_Position', 'Splice_Delta', 'Raw_Change']
        available_cols = [c for c in key_cols if c in self.filtered_df.columns]
        self.filtered_df[available_cols].to_csv(output_csv, index=False)

        # 生成摘要
        with open(summary_txt, 'w') as f:
            f.write(f"Proteo-Structure Pipeline - Phase 1 Summary\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write(f"{'='*50}\n\n")
            f.write(f"Input File: {self.excel_path}\n")
            f.write(f"Total Input Variants: {len(self.df)}\n")
            f.write(f"After Missense Filter: {len(self.df[self.df['Molecular.consequence'].astype(str).str.contains('missense variant', na=False)])}\n")
            f.write(f"After Splice Filter: {len(self.filtered_df)}\n")
            f.write(f"Final Parsed Mutations: {len(self.filtered_df)}\n\n")

            f.write(f"Output Files:\n")
            f.write(f"  - CSV: {output_csv}\n")
            f.write(f"  - WT FASTA: {self.wt_dir}/VWF_P04275_WT.fasta\n")
            f.write(f"  - Mutant FASTAs: {self.mutant_dir}/\n")

        print(f"\n{'='*50}")
        print(f"Phase 1 完成!")
        print(f"结果保存在: {self.output_dir}")
        print(f"  - CSV: {output_csv}")
        print(f"  - 摘要: {summary_txt}")

        return {
            'csv': str(output_csv),
            'summary': str(summary_txt),
            'wt_fasta': str(self.wt_dir / 'VWF_P04275_WT.fasta'),
            'mutant_dir': str(self.mutant_dir)
        }

    def run(self) -> Dict[str, str]:
        """运行完整的Phase 1流程"""
        print("="*60)
        print("Proteo-Structure Pipeline - Phase 1: Smart Filter & Parsing")
        print("Target: VWF (Uniprot: P04275)")
        print("="*60)

        self.load_excel()
        self.filter_missense_variants()
        self.filter_splice_disruptors(threshold=0.3)
        self.extract_mutations()
        self.download_wt_fasta(uniprot_id='P04275')
        self.generate_mutant_fasta()
        results = self.save_results()

        return results


def main():
    parser = argparse.ArgumentParser(description='Phase 1: Smart Filter & Parsing for VWF')
    parser.add_argument('--excel', type=str,
                       default='../Final_VWF_Target_List_with_AlphaGenome.xlsx',
                       help='输入Excel文件路径')
    parser.add_argument('--output-dir', type=str, default='../output',
                       help='输出目录')
    parser.add_argument('--splice-threshold', type=float, default=0.3,
                       help='Splice filter阈值 (默认: 0.3)')

    args = parser.parse_args()

    # 如果路径是相对路径，转换为绝对路径
    excel_path = args.excel
    if not os.path.isabs(excel_path):
        excel_path = os.path.join(os.path.dirname(__file__), excel_path)

    filter_tool = Phase1Filter(excel_path, args.output_dir)
    filter_tool.run()


if __name__ == '__main__':
    main()

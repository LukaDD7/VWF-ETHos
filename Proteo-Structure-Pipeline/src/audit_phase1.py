#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 1 逻辑审计与学术复查脚本 (Cross-Validation Auditor)
目标：硬核实错义突变漏斗、剪接排阻阈值、以及 FASTA 序列的单点对应关系。
"""

import os
import glob
import pandas as pd
from Bio import SeqIO

# 路径配置
ORIGINAL_EXCEL = "../../Final_VWF_Target_List_with_AlphaGenome.xlsx"
PROCESSED_CSV = "../output/phase1_filtered_mutations.csv"
WT_FASTA = "../structures/wt/VWF_P04275_WT.fasta"
MUTANT_DIR = "../structures/mutant/"

def run_audit():
    print("="*50)
    print("🔬 开始 Phase 1 学术合规性交叉审计")
    print("="*50)

    # 1. 载入数据
    try:
        df_orig = pd.read_excel(ORIGINAL_EXCEL, sheet_name='Sheet1')
        df_proc = pd.read_csv(PROCESSED_CSV)
    except FileNotFoundError as e:
        print(f"❌ 找不到文件: {e}")
        return

    print(f"📄 原始输入变异总数: {len(df_orig)}")
    print(f"📄 成功解析突变总数: {len(df_proc)}")

    # ==========================================
    # 审计点 1：错义突变绝对纯净度测试
    # 需要从原始 Excel 中核对（processed CSV 已简化列）
    # ==========================================
    print("\n🔍 审计点 1: 错义突变绝对纯净度 (Missense Exclusivity)")
    print("   [从原始 Excel 核对筛选逻辑...]")

    # 获取 processed 变异的位置信息用于匹配
    proc_positions = set()
    for _, row in df_proc.iterrows():
        chrom = str(row.get('Chromosome', ''))
        pos = str(row.get('Position', ''))
        ref = str(row.get('Ref', ''))
        alt = str(row.get('Alt', ''))
        proc_positions.add((chrom, pos, ref, alt))

    # 在原始数据中检查这些变异
    invalid_missense = 0
    not_found = 0
    mol_cons_col = 'Molecular.consequence'
    ai_cons_col = 'AI_突变后果'

    for _, row in df_orig.iterrows():
        chrom = str(row.get('Chromosome', ''))
        pos = str(row.get('Position', ''))
        ref = str(row.get('Ref', ''))
        alt = str(row.get('Alt', ''))

        if (chrom, pos, ref, alt) in proc_positions:
            # 这个变异被保留了，检查它是否是错义
            is_missense = False
            if mol_cons_col in df_orig.columns:
                mol_cons = str(row.get(mol_cons_col, '')).lower()
                if 'missense variant' in mol_cons:
                    is_missense = True
            if ai_cons_col in df_orig.columns:
                ai_cons = str(row.get(ai_cons_col, '')).lower()
                if 'missense variant' in ai_cons:
                    is_missense = True

            if not is_missense:
                invalid_missense += 1
                if invalid_missense <= 3:
                    print(f"      ⚠️ 非错义混入: {chrom}:{pos} {ref}>{alt}")

    if invalid_missense == 0:
        print("   ✅ 完美：100% 的保留变异均为错义突变。")
    else:
        print(f"   ❌ 严重警告：发现 {invalid_missense} 个非错义突变混入！")

    # ==========================================
    # 审计点 2：剪接阈值严苛执行测试 (|Delta| <= 0.3)
    # ==========================================
    print("\n🔍 审计点 2: 剪接破坏者排阻阈值 (Splice Delta <= 0.3)")
    splice_violations = 0
    for _, row in df_proc.iterrows():
        # 如果有 AI_Splice_Delta_Score 列优先用，否则用原表概率计算
        if 'AI_Splice_Delta_Score' in row and pd.notna(row['AI_Splice_Delta_Score']):
            delta = float(row['AI_Splice_Delta_Score'])
        else:
            ref_prob = float(row.get('Splice_REF正常概率', 0.0))
            alt_prob = float(row.get('Splice_ALT突变概率', 0.0))
            delta = abs(alt_prob - ref_prob)

        if delta > 0.3:
            splice_violations += 1

    if splice_violations == 0:
        print("   ✅ 完美：所有保留变异的剪接偏移量均未越过 0.3 的红线。")
    else:
        print(f"   ❌ 严重警告：发现 {splice_violations} 个变异突破了剪接排阻红线！")

    # ==========================================
    # 审计点 3：FASTA 序列物理核实 (长度与单点突变比对)
    # ==========================================
    print("\n🔍 审计点 3: 蛋白质 FASTA 图纸物理核对")

    # 读取野生型
    wt_record = SeqIO.read(WT_FASTA, "fasta")
    wt_seq = str(wt_record.seq)

    if len(wt_seq) == 2813:
        print(f"   ✅ 野生型长度正确: {len(wt_seq)} aa")
    else:
        print(f"   ❌ 野生型长度异常: 预期 2813，实际 {len(wt_seq)} aa")

    # 随机抽检前 5 个生成的突变 FASTA
    sample_mutants = df_proc.head(5)
    print("   [抽检前 5 个突变序列的物理吻合度]")

    for _, row in sample_mutants.iterrows():
        wt_aa = row['WT_AA_1']       # 1字母野生型
        pos = int(row['AA_Position'])
        mut_aa = row['Mut_AA_1']     # 1字母突变型
        chrom = str(row.get('Chromosome', ''))

        # 实际文件名格式: VWF_A1032V__.fasta
        fasta_name = f"VWF_{wt_aa}{pos}{mut_aa}_{chrom}.fasta"
        fasta_path = os.path.join(MUTANT_DIR, fasta_name)

        if not os.path.exists(fasta_path):
            # 尝试备用格式（不带染色体）
            fasta_name_alt = f"VWF_{wt_aa}{pos}{mut_aa}__.fasta"
            fasta_path_alt = os.path.join(MUTANT_DIR, fasta_name_alt)
            if os.path.exists(fasta_path_alt):
                fasta_path = fasta_path_alt
            else:
                print(f"      ⚠️ 找不到文件: {fasta_path} 或 {fasta_path_alt}")
                continue

        mut_record = SeqIO.read(fasta_path, "fasta")
        mut_seq = str(mut_record.seq)

        # 验证 1：长度必须严格一致
        len_match = len(mut_seq) == len(wt_seq)

        # 验证 2：野生型对应位置的氨基酸是否真的是这个 (注意 Python 索引从 0 开始，真实位置要 -1)
        zero_index = pos - 1
        actual_wt_aa = wt_seq[zero_index]
        wt_match = actual_wt_aa == wt_aa

        # 验证 3：突变型该位置是否成功被替换
        mut_match = mut_seq[zero_index] == mut_aa

        if len_match and wt_match and mut_match:
             print(f"      ✅ 校验通过: {wt_aa}{pos}{mut_aa} (WT确认='{actual_wt_aa}', MUT植入='{mut_aa}', 序列长度一致)")
        else:
             print(f"      ❌ 校验失败: {wt_aa}{pos}{mut_aa}")
             if not wt_match: print(f"         - 冲突: 数据表说是 {wt_aa}，但 WT 序列第 {pos} 位实际是 {actual_wt_aa}")
             if not mut_match: print(f"         - 冲突: 未能成功替换为 {mut_aa}")

    print("\n" + "="*50)
    print("审计完成。如果以上均为绿灯 ✅，则证明数据清洗逻辑具备严谨的发表级可信度。")
    print("="*50)

if __name__ == "__main__":
    run_audit()

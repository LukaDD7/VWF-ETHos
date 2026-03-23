#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
重新生成AlphaFold3 Batch JSON（基于修复后的数据）

基于: Final_VWF_Target_List_with_AlphaGenome_FIXED.xlsx
筛选条件:
1. Molecular.consequence = 'missense variant'
2. AlphaGenome_Max_Score 不为空（有AlphaGenome数据）
3. |Splice_REF - Splice_ALT| <= 0.3（剪接位点正常）
4. 排除已完成的变异（batch_01_completed目录）
"""

import json
import logging
import pandas as pd
import requests
from pathlib import Path
from typing import Dict, List

# ==================== 配置 ====================
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/results")
OUTPUT_DIR = BASE_DIR / "output" / "af3_batches"
COMPLETED_DIR = BASE_DIR / "AF3_Results" / "batch_01_completed_20250318"

INPUT_EXCEL = RESULTS_DIR / "Final_VWF_Target_List_with_AlphaGenome_FIXED.xlsx"
VWF_UNIPROT_ID = "P04275"
BATCH_SIZE = 30

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)


# 氨基酸映射
AA_3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*'
}


def get_wt_sequence() -> str:
    """从Uniprot获取VWT WT序列"""
    cache_file = OUTPUT_DIR / "VWF_WT.fasta"
    if cache_file.exists():
        with open(cache_file, 'r') as f:
            lines = f.readlines()
            return ''.join(line.strip() for line in lines if not line.startswith('>'))

    logger.info(f"从Uniprot下载WT序列: {VWF_UNIPROT_ID}")
    url = f"https://rest.uniprot.org/uniprotkb/{VWF_UNIPROT_ID}.fasta"
    response = requests.get(url)
    response.raise_for_status()

    sequence = ''.join(line.strip() for line in response.text.split('\n') if not line.startswith('>'))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache_file, 'w') as f:
        f.write(response.text)

    return sequence


def parse_protein_change(change_str: str) -> tuple:
    """
    解析蛋白质变化，支持:
    - 1字母: G1531D
    - 3字母: Gly1531Asp (从p.列)
    """
    if pd.isna(change_str) or change_str == '':
        return None, None, None

    change_str = str(change_str).strip()

    # 1字母格式 (如 G1531D)
    if len(change_str) >= 3 and change_str[0].isalpha() and change_str[-1].isalpha():
        wt_aa = change_str[0]
        mut_aa = change_str[-1]
        try:
            pos = int(change_str[1:-1])
            return wt_aa, pos, mut_aa
        except ValueError:
            pass

    return None, None, None


def generate_mutant_sequence(wt_seq: str, wt_aa: str, pos: int, mut_aa: str) -> str:
    """生成突变体序列"""
    if pos < 1 or pos > len(wt_seq):
        raise ValueError(f"位置{pos}超出序列范围(1-{len(wt_seq)})")

    actual_aa = wt_seq[pos - 1]  # 序列是0-based
    if actual_aa != wt_aa:
        logger.warning(f"WT氨基酸不匹配: 预期{wt_aa}@{pos}, 实际{actual_aa}")
        # 继续处理，但记录警告

    mutant_seq = wt_seq[:pos - 1] + mut_aa + wt_seq[pos:]
    return mutant_seq


def load_completed_variants() -> set:
    """加载已完成的变异"""
    completed = set()
    if not COMPLETED_DIR.exists():
        return completed

    for zip_file in COMPLETED_DIR.glob("*.zip"):
        # 文件名格式: fold_vwf_g1531d.zip -> VWF_G1531D
        name = zip_file.stem.replace("fold_vwf_", "")
        name_upper = name.upper()
        completed.add(f"VWF_{name_upper}")

    logger.info(f"已完成的变异: {len(completed)} 个")
    return completed


def create_batch_job(name: str, sequence: str) -> dict:
    """创建单个batch任务"""
    return {
        "name": name,
        "sequences": [{
            "proteinChain": {
                "sequence": sequence,
                "count": 1
            }
        }]
    }


def main():
    logger.info("=" * 70)
    logger.info("重新生成AlphaFold3 Batch JSON（基于修复后的数据）")
    logger.info("=" * 70)

    # 1. 读取修复后的Excel
    logger.info(f"读取: {INPUT_EXCEL}")
    df = pd.read_excel(INPUT_EXCEL)
    logger.info(f"总变异数: {len(df)}")

    # 2. 筛选条件
    logger.info("\n应用筛选条件...")

    # 2.1 错义突变
    missense = df[df['Molecular.consequence'] == 'missense variant'].copy()
    logger.info(f"1. 错义突变: {len(missense)}")

    # 2.2 有AlphaGenome数据
    missense_ag = missense[missense['AlphaGenome_Max_Score'].notna()].copy()
    logger.info(f"2. 有AlphaGenome数据: {len(missense_ag)}")

    # 2.3 剪接位点正常 (|Splice_REF - Splice_ALT| <= 0.3)
    missense_ag['splice_delta'] = abs(missense_ag['Splice_REF正常概率'] - missense_ag['Splice_ALT突变概率'])
    missense_good = missense_ag[missense_ag['splice_delta'] <= 0.3].copy()
    excluded_splice = missense_ag[missense_ag['splice_delta'] > 0.3]
    logger.info(f"3. 剪接位点正常(|Δ|<=0.3): {len(missense_good)} (排除剪接异常: {len(excluded_splice)})")

    # 2.3.5 去重：同一个Protein.change只保留第一个
    missense_good = missense_good.drop_duplicates(subset=['Protein.change'], keep='first')
    logger.info(f"3.5 去重后: {len(missense_good)} (同一氨基酸变化只保留第一个)")

    # 2.4 排除已完成的
    completed = load_completed_variants()
    if completed:
        # 从Protein.change生成任务名
        def get_task_name(change):
            wt_aa, pos, mut_aa = parse_protein_change(change)
            if wt_aa and pos and mut_aa:
                return f"VWF_{wt_aa}{pos}{mut_aa}"
            return None

        missense_good['task_name'] = missense_good['Protein.change'].apply(get_task_name)
        missense_todo = missense_good[~missense_good['task_name'].isin(completed)].copy()
        logger.info(f"4. 排除已完成({len(completed)}): {len(missense_todo)} 个待处理")
    else:
        missense_todo = missense_good.copy()
        logger.info(f"4. 无已完成记录，全部待处理: {len(missense_todo)}")

    # 3. 获取WT序列
    logger.info("\n获取VWT WT序列...")
    wt_seq = get_wt_sequence()
    logger.info(f"WT序列长度: {len(wt_seq)} aa")

    # 4. 生成任务列表
    logger.info("\n生成任务...")
    jobs = []

    # 先生成所有WT任务（只生成一个）
    jobs.append(create_batch_job("VWF_WT", wt_seq))

    # 生成变异任务
    skipped = []
    for _, row in missense_todo.iterrows():
        change = row['Protein.change']
        wt_aa, pos, mut_aa = parse_protein_change(change)

        if not wt_aa:
            skipped.append((change, "无法解析Protein.change"))
            continue

        try:
            mut_seq = generate_mutant_sequence(wt_seq, wt_aa, pos, mut_aa)
            task_name = f"VWF_{wt_aa}{pos}{mut_aa}"
            jobs.append(create_batch_job(task_name, mut_seq))
        except Exception as e:
            skipped.append((change, str(e)))

    logger.info(f"总任务数: {len(jobs)} (WT: 1 + 变异: {len(jobs)-1})")
    if skipped:
        logger.warning(f"跳过 {len(skipped)} 个变异: {skipped[:5]}")

    # 5. 拆分为batch文件
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 清理旧的batch文件（保留_REMAINING备份）
    for old_file in OUTPUT_DIR.glob("af3_batch_*.json"):
        if '_REMAINING' not in old_file.name:
            old_file.unlink()

    # 生成新的batch文件
    batch_count = 0
    for i in range(0, len(jobs), BATCH_SIZE):
        batch = jobs[i:i + BATCH_SIZE]
        batch_count += 1
        batch_file = OUTPUT_DIR / f"af3_batch_{batch_count:02d}.json"

        with open(batch_file, 'w') as f:
            json.dump(batch, f, indent=2)

    logger.info(f"\n生成 {batch_count} 个batch文件，每个最多 {BATCH_SIZE} 任务")

    # 6. 生成manifest
    manifest = {
        "total_jobs": len(jobs),
        "wt_jobs": 1,
        "variant_jobs": len(jobs) - 1,
        "batch_count": batch_count,
        "batch_size": BATCH_SIZE,
        "batches": [f"af3_batch_{i:02d}.json" for i in range(1, batch_count + 1)],
        "excluded_splice_disruptors": len(excluded_splice),
        "excluded_completed": len(completed) if completed else 0,
        "skipped_parsing": len(skipped)
    }

    manifest_file = OUTPUT_DIR / "af3_jobs_manifest.json"
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)

    logger.info(f"\n任务清单: {manifest_file}")
    logger.info("=" * 70)
    logger.info("重新生成完成！")
    logger.info(f"输出目录: {OUTPUT_DIR}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()

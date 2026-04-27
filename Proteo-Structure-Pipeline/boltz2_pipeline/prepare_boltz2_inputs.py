#!/usr/bin/env python3
"""
VWF Boltz-2 Input Generator (Blind-Scan Mode)
==============================================
表格预处理 + Boltz-2 JSON 批次生成

设计原则（无上帝视角，诊疗智能体模式）：
  - 不预设分型，对每个突变位点自动关联其所在结构域的全部配体
  - 每个突变会生成 N 个独立 job（N = 该位点相关配体数量）
  - 供 Boltz-2 并行计算亲和力，再由下游分析脚本差异比对

输入：
  - 标准化突变表（FIXED Excel 或 CSV，使用 VWF 统一位置编号）
  - VWF WT FASTA

输出：
  - output/boltz2_blind_scan/batch_XX.json（每批 ≤ batch_size 个 job）
  - output/boltz2_blind_scan/job_manifest.csv（全部 job 的索引表）

用法：
  python prepare_boltz2_inputs.py \\
    --excel ../../results/Final_VWF_Target_List_with_AlphaGenome_FIXED.xlsx \\
    --wt-fasta ../structures/wt/VWF_P04275_WT.fasta \\
    --output-dir ../../output/boltz2_blind_scan \\
    --batch-size 30
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

# 引入本地配体数据库（避免任何外部 API 调用）
sys.path.insert(0, str(Path(__file__).parent))
from vwf_ligand_database import (
    VWF_DOMAIN_MAP,
    VWF_LIGAND_DATABASE,
    get_domain_for_position,
    get_ligands_for_position,
)


# =============================================================================
# 工具函数
# =============================================================================

def read_fasta(path: str) -> str:
    """读取 FASTA 文件，返回纯序列字符串。"""
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def parse_mutation(s: str) -> Optional[Tuple[str, int, str]]:
    """
    解析多种格式的氨基酸变化字符串。
    支持: G1531D, R1306Trp, p.Arg1306Trp, NM_*..(VWF):..(p.Gly1531Asp)
    返回 (ref_aa_1letter, position, alt_aa_1letter) 或 None
    """
    AA3 = {
        "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C",
        "Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
        "Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
        "Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
    }
    if not s or pd.isna(s):
        return None
    s = str(s).strip()

    # 优先从 Name 列的 (p.XXXNNNyyy) 里提取
    match_nm = re.search(r'\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})\)', s)
    if match_nm:
        ref3, pos, alt3 = match_nm.groups()
        ref = AA3.get(ref3.capitalize())
        alt = AA3.get(alt3.capitalize())
        if ref and alt:
            return ref, int(pos), alt

    # 3-letter: Pro1266Leu
    match3 = re.match(r'^([A-Za-z]{3})(\d+)([A-Za-z]{3})$', s)
    if match3:
        ref3, pos, alt3 = match3.groups()
        ref = AA3.get(ref3.capitalize())
        alt = AA3.get(alt3.capitalize())
        if ref and alt:
            return ref, int(pos), alt

    # 1-letter: G1531D
    match1 = re.match(r'^([A-Z])(\d+)([A-Z])$', s)
    if match1:
        ref, pos, alt = match1.groups()
        return ref, int(pos), alt

    return None


def load_variants_from_excel(path: str) -> pd.DataFrame:
    """
    从 FIXED Excel 中提取 missense 突变。
    自动适配列名，返回标准 DataFrame:
      [Variant_ID, WT_AA, Position, Mut_AA, raw_change]
    """
    df = pd.read_excel(path)
    
    # 找突变列（Protein.change 优先，其次 p., 再 Name）
    priority_cols = ["Protein.change", "p.", "Protein_change", "AA_Change"]
    name_col = "Name"
    change_col = next(
        (c for c in priority_cols if c in df.columns), None
    )
    
    rows = []
    for _, row in df.iterrows():
        # 先试主突变列，再试 Name 列里的注释
        parsed = None
        if change_col:
            parsed = parse_mutation(row.get(change_col, ""))
        if parsed is None and name_col in df.columns:
            parsed = parse_mutation(row.get(name_col, ""))
        if parsed is None:
            continue
        
        ref, pos, alt = parsed
        # 过滤非 missense（包含 fs / * / splice）
        if alt in ("*", "X") or "fs" in str(row.get(change_col, "")):
            continue
        
        rows.append({
            "Variant_ID": f"VWF_{ref}{pos}{alt}",
            "WT_AA": ref,
            "Position": pos,
            "Mut_AA": alt,
        })
    
    df_out = pd.DataFrame(rows).drop_duplicates(
        subset=["WT_AA", "Position", "Mut_AA"]
    ).reset_index(drop=True)
    return df_out


def apply_mutation_and_slice(
    wt_seq: str,
    ref: str,
    pos: int,
    alt: str,
    domain_start: int,
    domain_end: int,
) -> Optional[str]:
    """
    1. 验证 WT 氨基酸匹配
    2. 在全长 WT 序列上应用点突变
    3. 返回指定结构域的切片（domain_start, domain_end 均为 1-indexed）
    """
    if pos < 1 or pos > len(wt_seq):
        return None
    actual = wt_seq[pos - 1]
    if actual != ref:
        print(f"    [WARN] WT mismatch at {pos}: expected {ref}, found {actual}")
        # 不强制跳过，但标注警告
    mutated = wt_seq[:pos - 1] + alt + wt_seq[pos:]
    return mutated[domain_start - 1: domain_end]


# =============================================================================
# Boltz-2 JSON 生成核心
# =============================================================================

def build_boltz2_job(
    variant_id: str,
    vwf_domain_seq: str,
    ligand_key: str,
    ligand_seq: str,
    n_chains: int,
) -> dict:
    """
    构建单个 Boltz-2 job 对象。
    Chain A = 突变 VWF 结构域
    Chain B, C, D ... = 配体链（支持三股螺旋）
    """
    sequences = [
        {
            "protein": {
                "id": "A",
                "sequence": vwf_domain_seq,
            }
        }
    ]
    for i in range(n_chains):
        sequences.append({
            "protein": {
                "id": chr(ord("B") + i),
                "sequence": ligand_seq,
            }
        })
    
    return {
        "name": f"{variant_id}_vs_{ligand_key}",
        "sequences": sequences,
        "properties": [
            {"affinity": {"binder": "A"}}   # Boltz-2 亲和力评分，A链=VWF
        ],
    }


def generate_batches(
    variants_df: pd.DataFrame,
    wt_seq: str,
    output_dir: Path,
    batch_size: int = 30,
) -> pd.DataFrame:
    """
    主生成函数：遍历所有突变 × 所有相关配体，生成 Boltz-2 batch JSON 文件。
    返回 job_manifest DataFrame。
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    all_jobs = []       # 所有 job dict（用于分批写 JSON）
    manifest_rows = []  # 索引表

    skipped = 0
    total = len(variants_df)

    for i, row in variants_df.iterrows():
        ref      = row["WT_AA"]
        pos      = int(row["Position"])
        alt      = row["Mut_AA"]
        var_id   = row["Variant_ID"]

        # 1. 确定结构域
        domain = get_domain_for_position(pos)
        if domain is None:
            print(f"  [SKIP] {var_id}: pos {pos} not in any known domain")
            skipped += 1
            continue

        d_start, d_end = VWF_DOMAIN_MAP[domain]

        # 2. 获取该位点所有相关配体（盲扫模式）
        ligands = get_ligands_for_position(pos)
        if not ligands:
            print(f"  [SKIP] {var_id} in domain {domain}: no ligands defined")
            skipped += 1
            continue

        # 3. 生成突变结构域序列
        domain_seq = apply_mutation_and_slice(wt_seq, ref, pos, alt, d_start, d_end)
        if domain_seq is None:
            print(f"  [SKIP] {var_id}: sequence slicing failed")
            skipped += 1
            continue

        # 4. 每个配体生成一个独立 job
        for lig in ligands:
            job = build_boltz2_job(
                variant_id=var_id,
                vwf_domain_seq=domain_seq,
                ligand_key=lig.key,
                ligand_seq=lig.sequence,
                n_chains=lig.n_chains,
            )
            all_jobs.append(job)
            manifest_rows.append({
                "job_name":       job["name"],
                "variant_id":     var_id,
                "position":       pos,
                "wt_aa":          ref,
                "mut_aa":         alt,
                "domain":         domain,
                "domain_range":   f"{d_start}-{d_end}",
                "domain_seq_len": len(domain_seq),
                "ligand_key":     lig.key,
                "ligand_uniprot": lig.uniprot_id,
                "n_ligand_chains":lig.n_chains,
                "vwf_subtypes":   "|".join(lig.vwf_subtypes),
                "interaction_type":lig.interaction_type,
            })

        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{total} variants...")

    # 5. 分批写出 JSON
    n_batches = (len(all_jobs) + batch_size - 1) // batch_size
    print(f"\nTotal jobs: {len(all_jobs)} | Batches: {n_batches} | Skipped: {skipped}")

    for b_idx in range(n_batches):
        batch_jobs = all_jobs[b_idx * batch_size: (b_idx + 1) * batch_size]
        batch_data = {
            "name": f"vwf_blind_scan_batch_{b_idx+1:03d}",
            "jobs": batch_jobs,
        }
        out_path = output_dir / f"batch_{b_idx+1:03d}.json"
        with open(out_path, "w") as f:
            json.dump(batch_data, f, indent=2)
        print(f"  Saved: {out_path.name}  [{len(batch_jobs)} jobs]")

    # 6. 保存 manifest
    manifest_df = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / "job_manifest.csv"
    manifest_df.to_csv(manifest_path, index=False)
    print(f"\nJob manifest saved: {manifest_path}")
    return manifest_df


# =============================================================================
# CLI 入口
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="VWF Boltz-2 Input Generator (Blind-Scan, No Subtype Assumption)"
    )
    parser.add_argument(
        "--excel", type=str,
        default="../../results/Final_VWF_Target_List_with_AlphaGenome_FIXED.xlsx",
        help="Input variant Excel file"
    )
    parser.add_argument(
        "--wt-fasta", type=str,
        default="../structures/wt/VWF_P04275_WT.fasta",
        help="VWF WT FASTA file"
    )
    parser.add_argument(
        "--output-dir", type=str,
        default="../../output/boltz2_blind_scan",
        help="Output directory for JSON batches"
    )
    parser.add_argument(
        "--batch-size", type=int, default=30,
        help="Max jobs per batch JSON (default: 30, matches Boltz-2 server quota)"
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Limit to first N variants (for testing)"
    )
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    excel_path = (script_dir / args.excel).resolve()
    fasta_path = (script_dir / args.wt_fasta).resolve()
    output_dir = (script_dir / args.output_dir).resolve()

    print("=" * 70)
    print("VWF Boltz-2 Blind-Scan Input Generator")
    print("=" * 70)
    print(f"Excel   : {excel_path}")
    print(f"WT FASTA: {fasta_path}")
    print(f"Output  : {output_dir}")
    print()

    # 加载 WT 序列
    wt_seq = read_fasta(str(fasta_path))
    print(f"WT sequence loaded: {len(wt_seq)} aa")

    # 加载变异列表
    variants_df = load_variants_from_excel(str(excel_path))
    print(f"Variants parsed: {len(variants_df)}")
    if args.limit:
        variants_df = variants_df.head(args.limit)
        print(f"  (Limited to first {args.limit})")
    print()

    # 生成批次
    manifest = generate_batches(variants_df, wt_seq, output_dir, args.batch_size)

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(manifest.groupby(["domain", "ligand_key"])["job_name"].count().to_string())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本 7：VCF 全量变异打分引擎 (终极版：11大模态全景提取 + PKL/CSV 双轨存档)

输入：../data/Clinvar_HGMD_merge_Hg37_fixed.vcf
      ../data/hg19ToHg38.over.chain.gz
输出：../results/07_VCF_AlphaGenome_Results.csv
      ../results/07_VCF_AlphaGenome_Results.pkl
"""

import logging
import sys
import time
import pickle
import json
import pandas as pd
import numpy as np
from pathlib import Path
from pyliftover import LiftOver
from dataclasses import dataclass, field, asdict
from typing import Any

# ==================== 配置 ====================
# 基础目录配置
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"

INPUT_VCF = DATA_DIR / "Clinvar_HGMD_merge_Hg37_fixed.vcf"
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
OUTPUT_CSV = RESULTS_DIR / "07_VCF_AlphaGenome_Results.csv"
OUTPUT_PKL = RESULTS_DIR / "07_VCF_AlphaGenome_Results.pkl"
CHECKPOINT_JSONL = RESULTS_DIR / "07_checkpoint.jsonl"  # 断点续传文件

API_KEY = "AIzaSyC25qItDVDM6jJJFqNDSvEyhR4P56E0CtY"
ONTOLOGY_TERMS = ['CL:0000115']  # 内皮细胞
CHECKPOINT_INTERVAL = 50

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)

# ==================== 数据结构 (扩充至 11 模态) ====================
@dataclass
class InferenceResult:
    chromosome: str
    position: int
    ref: str
    alt: str
    interval_start: int
    interval_end: int

    # 11 模态的最高破坏极值
    rna_seq_delta: float = 0.0
    splice_sites_delta: float = 0.0
    cage_delta: float = 0.0
    procap_delta: float = 0.0
    splice_site_usage_delta: float = 0.0
    splice_junctions_delta: float = 0.0
    dnase_delta: float = 0.0
    atac_delta: float = 0.0
    chip_histone_delta: float = 0.0
    chip_tf_delta: float = 0.0
    contact_maps_delta: float = 0.0

    delta_score: float = 0.0  # 核心预警总分 (通常取 RNA 和 Splice 的最大值)
    raw_outputs: Any = field(default=None)
    success: bool = False
    error_message: str | None = None

def parse_vcf(vcf_path):
    logger.info(f"正在解析 VCF 文件: {vcf_path}")
    data = []
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                chrom = parts[0].replace('chr', '')
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                if '<' in alt or '[' in alt or ']' in alt: continue
                data.append({"hg19_chr": chrom, "hg19_pos": pos, "ref": ref, "alt": alt})
    df = pd.DataFrame(data)
    return df

def liftover_to_hg38(df, chain_path):
    lo = LiftOver(str(chain_path))
    hg38_chroms, hg38_positions = [], []
    for _, row in df.iterrows():
        res = lo.convert_coordinate(f"chr{row['hg19_chr']}", row['hg19_pos'] - 1)
        if res:
            hg38_chroms.append(res[0][0])
            hg38_positions.append(res[0][1] + 1)
        else:
            hg38_chroms.append(None)
            hg38_positions.append(None)
    df['hg38_chr'] = hg38_chroms
    df['hg38_pos'] = hg38_positions
    df_valid = df.dropna(subset=['hg38_pos']).copy()
    df_valid['hg38_pos'] = df_valid['hg38_pos'].astype(int)
    return df_valid

def init_alphagenome():
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    model = dna_client.create(API_KEY)
    return dna_client, genome, model

def extract_all_scores(outputs):
    """自动提取 11 个模态的最大绝对差值"""
    scores = {}
    def get_max_delta(ref_track, alt_track):
        if ref_track is not None and alt_track is not None:
            ref_val = ref_track.values
            alt_val = alt_track.values
            if ref_val.size > 0 and alt_val.size > 0:
                return float(np.max(np.abs(alt_val - ref_val)))
        return 0.0

    scores['rna_seq'] = get_max_delta(outputs.reference.rna_seq, outputs.alternate.rna_seq)
    scores['splice_sites'] = get_max_delta(outputs.reference.splice_sites, outputs.alternate.splice_sites)
    scores['cage'] = get_max_delta(outputs.reference.cage, outputs.alternate.cage)
    scores['procap'] = get_max_delta(outputs.reference.procap, outputs.alternate.procap)
    scores['splice_site_usage'] = get_max_delta(outputs.reference.splice_site_usage, outputs.alternate.splice_site_usage)
    scores['splice_junctions'] = get_max_delta(outputs.reference.splice_junctions, outputs.alternate.splice_junctions)
    scores['dnase'] = get_max_delta(outputs.reference.dnase, outputs.alternate.dnase)
    scores['atac'] = get_max_delta(outputs.reference.atac, outputs.alternate.atac)
    scores['chip_histone'] = get_max_delta(outputs.reference.chip_histone, outputs.alternate.chip_histone)
    scores['chip_tf'] = get_max_delta(outputs.reference.chip_tf, outputs.alternate.chip_tf)
    scores['contact_maps'] = get_max_delta(outputs.reference.contact_maps, outputs.alternate.contact_maps)

    return scores


# ==================== 断点续传机制 ====================

def get_variant_key_from_row(row) -> str:
    """生成变异的唯一标识键。"""
    chrom = row['hg38_chr'] if 'hg38_chr' in row else row.get('chrom', '')
    pos = row['hg38_pos'] if 'hg38_pos' in row else row.get('pos', 0)
    ref = row['ref']
    alt = row['alt']
    return f"{chrom}:{pos}:{ref}:{alt}"


def load_completed_variants(checkpoint_path: Path) -> set[str]:
    """从 checkpoint 文件加载已完成的变异标识。"""
    completed = set()
    if not checkpoint_path.exists():
        return completed

    logger.info(f"检测到已有 checkpoint 文件：{checkpoint_path}")
    try:
        with open(checkpoint_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                    key = f"{data['hg38_chr']}:{data['hg38_pos']}:{data['ref']}:{data['alt']}"
                    completed.add(key)
                except (json.JSONDecodeError, KeyError) as e:
                    logger.warning(f"跳过损坏的 checkpoint 行: {e}")
        logger.info(f"已加载 {len(completed)} 个已完成的变异记录")
    except Exception as e:
        logger.error(f"读取 checkpoint 失败: {e}")
    return completed


def append_result_to_checkpoint(result: InferenceResult, checkpoint_path: Path):
    """立即追加单个结果到 checkpoint 文件（JSONL 格式）。"""
    record = {
        "hg19_chr": getattr(result, 'hg19_chr', ''),
        "hg19_pos": getattr(result, 'hg19_pos', 0),
        "hg38_chr": result.chromosome,
        "hg38_pos": result.position,
        "ref": result.ref,
        "alt": result.alt,
        "interval_start": result.interval_start,
        "interval_end": result.interval_end,
        "rna_seq_delta": result.rna_seq_delta,
        "splice_sites_delta": result.splice_sites_delta,
        "cage_delta": result.cage_delta,
        "procap_delta": result.procap_delta,
        "splice_site_usage_delta": result.splice_site_usage_delta,
        "splice_junctions_delta": result.splice_junctions_delta,
        "dnase_delta": result.dnase_delta,
        "atac_delta": result.atac_delta,
        "chip_histone_delta": result.chip_histone_delta,
        "chip_tf_delta": result.chip_tf_delta,
        "contact_maps_delta": result.contact_maps_delta,
        "delta_score": result.delta_score,
        "success": result.success,
        "error_message": result.error_message,
    }
    with open(checkpoint_path, 'a', encoding='utf-8') as f:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')


def load_existing_results_from_checkpoint(checkpoint_path: Path) -> tuple[list, list]:
    """从 checkpoint 文件加载所有已有结果。"""
    pkl_results = []
    csv_results = []
    if not checkpoint_path.exists():
        return pkl_results, csv_results

    try:
        with open(checkpoint_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                    # 重建 InferenceResult 对象
                    result = InferenceResult(
                        chromosome=data['hg38_chr'],
                        position=data['hg38_pos'],
                        ref=data['ref'],
                        alt=data['alt'],
                        interval_start=data['interval_start'],
                        interval_end=data['interval_end'],
                        rna_seq_delta=data.get('rna_seq_delta', 0.0),
                        splice_sites_delta=data.get('splice_sites_delta', 0.0),
                        cage_delta=data.get('cage_delta', 0.0),
                        procap_delta=data.get('procap_delta', 0.0),
                        splice_site_usage_delta=data.get('splice_site_usage_delta', 0.0),
                        splice_junctions_delta=data.get('splice_junctions_delta', 0.0),
                        dnase_delta=data.get('dnase_delta', 0.0),
                        atac_delta=data.get('atac_delta', 0.0),
                        chip_histone_delta=data.get('chip_histone_delta', 0.0),
                        chip_tf_delta=data.get('chip_tf_delta', 0.0),
                        contact_maps_delta=data.get('contact_maps_delta', 0.0),
                        delta_score=data.get('delta_score', 0.0),
                        success=data.get('success', False),
                        error_message=data.get('error_message'),
                    )
                    pkl_results.append(result)

                    # 重建 CSV 行
                    csv_row = {
                        "hg19_chr": data.get('hg19_chr', ''),
                        "hg19_pos": data.get('hg19_pos', 0),
                        "hg38_chr": data['hg38_chr'],
                        "hg38_pos": data['hg38_pos'],
                        "ref": data['ref'],
                        "alt": data['alt'],
                        "Status": "Success" if data.get('success', False) else "Failed",
                        "Delta_Max_Core": data.get('delta_score', 0.0),
                    }
                    if data.get('success', False):
                        for key in ['rna_seq', 'splice_sites', 'cage', 'procap',
                                   'splice_site_usage', 'splice_junctions', 'dnase',
                                   'atac', 'chip_histone', 'chip_tf', 'contact_maps']:
                            csv_row[f"AG_{key.upper()}"] = data.get(f'{key}_delta', 0.0)
                    csv_results.append(csv_row)
                except (json.JSONDecodeError, KeyError) as e:
                    logger.warning(f"跳过损坏的 checkpoint 记录: {e}")
        logger.info(f"从 checkpoint 恢复了 {len(pkl_results)} 个结果")
    except Exception as e:
        logger.error(f"恢复 checkpoint 失败: {e}")
    return pkl_results, csv_results


def main():
    df_vcf = parse_vcf(INPUT_VCF)
    df_hg38 = liftover_to_hg38(df_vcf, CHAIN_FILE)

    # ========== 断点续传：加载已完成的变异 ==========
    completed_variants = load_completed_variants(CHECKPOINT_JSONL)

    # 筛选出未完成的变异
    df_hg38['variant_key'] = df_hg38.apply(get_variant_key_from_row, axis=1)
    df_pending = df_hg38[~df_hg38['variant_key'].isin(completed_variants)].copy()
    df_pending = df_pending.drop(columns=['variant_key'])

    total = len(df_hg38)
    pending_count = len(df_pending)
    completed_count = len(completed_variants)

    logger.info(f"总变异: {total}, 已完成: {completed_count}, 待处理: {pending_count}")

    if pending_count == 0:
        logger.info("所有变异都已处理完成！正在生成最终输出文件...")
        pkl_results, csv_results = load_existing_results_from_checkpoint(CHECKPOINT_JSONL)
        with open(OUTPUT_PKL, "wb") as f:
            pickle.dump(pkl_results, f, protocol=pickle.HIGHEST_PROTOCOL)
        pd.DataFrame(csv_results).to_csv(OUTPUT_CSV, index=False)
        logger.info(f"============== 全部完成！==============")
        return

    dna_client, genome, model = init_alphagenome()
    # 开启全模态请求 (11 模态)
    req_outputs = [
        dna_client.OutputType.RNA_SEQ, dna_client.OutputType.SPLICE_SITES,
        dna_client.OutputType.CAGE, dna_client.OutputType.PROCAP,
        dna_client.OutputType.SPLICE_SITE_USAGE, dna_client.OutputType.SPLICE_JUNCTIONS,
        dna_client.OutputType.DNASE, dna_client.OutputType.ATAC,
        dna_client.OutputType.CHIP_HISTONE, dna_client.OutputType.CHIP_TF,
        dna_client.OutputType.CONTACT_MAPS
    ]

    # 加载已有结果
    pkl_results, csv_results = load_existing_results_from_checkpoint(CHECKPOINT_JSONL)

    logger.info("============== 开始批量请求 AlphaGenome (11 模态全开) ==============")
    for i, row in df_pending.iterrows():
        chrom = f"chr{row['hg38_chr']}" if not str(row['hg38_chr']).startswith('chr') else str(row['hg38_chr'])
        pos = row['hg38_pos']
        ref = row['ref']
        alt = row['alt']
        hg19_chr = row.get('hg19_chr', '')
        hg19_pos = row.get('hg19_pos', 0)

        half_len = 524288
        interval_start = max(1, pos - half_len)
        interval_end = pos + half_len

        result_obj = InferenceResult(
            chromosome=chrom, position=pos, ref=ref, alt=alt,
            interval_start=interval_start, interval_end=interval_end
        )
        # 保存原始坐标用于 checkpoint
        result_obj.hg19_chr = hg19_chr
        result_obj.hg19_pos = hg19_pos

        try:
            interval = genome.Interval(chromosome=chrom, start=interval_start, end=interval_end)
            variant = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)

            outputs = model.predict_variant(
                interval=interval, variant=variant,
                ontology_terms=ONTOLOGY_TERMS, requested_outputs=req_outputs
            )

            result_obj.raw_outputs = outputs
            scores = extract_all_scores(outputs)

            # 将字典赋值给对象属性
            result_obj.rna_seq_delta = scores['rna_seq']
            result_obj.splice_sites_delta = scores['splice_sites']
            result_obj.cage_delta = scores['cage']
            result_obj.procap_delta = scores['procap']
            result_obj.splice_site_usage_delta = scores['splice_site_usage']
            result_obj.splice_junctions_delta = scores['splice_junctions']
            result_obj.dnase_delta = scores['dnase']
            result_obj.atac_delta = scores['atac']
            result_obj.chip_histone_delta = scores['chip_histone']
            result_obj.chip_tf_delta = scores['chip_tf']
            result_obj.contact_maps_delta = scores['contact_maps']

            result_obj.delta_score = max(scores['rna_seq'], scores['splice_sites'])
            result_obj.success = True

            actual_count = completed_count + len(pkl_results) + 1
            logger.info(f"[{actual_count}/{total}] 成功: {chrom}:{pos} {ref}>{alt} | 核心预警: {result_obj.delta_score:.3f}")

        except Exception as e:
            actual_count = completed_count + len(pkl_results) + 1
            logger.warning(f"[{actual_count}/{total}] 失败: {chrom}:{pos} - {e}")
            result_obj.error_message = str(e)
            result_obj.success = False

        pkl_results.append(result_obj)

        # 写入 11 模态到 CSV
        csv_row = {
            "hg19_chr": hg19_chr, "hg19_pos": hg19_pos,
            "hg38_chr": chrom, "hg38_pos": pos,
            "ref": ref, "alt": alt, "Status": "Success" if result_obj.success else "Failed",
            "Delta_Max_Core": result_obj.delta_score
        }
        if result_obj.success:
            csv_row.update({f"AG_{k.upper()}": v for k, v in scores.items()})
        csv_results.append(csv_row)

        # ========== 断点续传：立即追加保存 ==========
        append_result_to_checkpoint(result_obj, CHECKPOINT_JSONL)

        time.sleep(1)

        # 定期保存（额外保险）
        if len(pkl_results) % CHECKPOINT_INTERVAL == 0:
            with open(OUTPUT_PKL, "wb") as f:
                pickle.dump(pkl_results, f, protocol=pickle.HIGHEST_PROTOCOL)
            pd.DataFrame(csv_results).to_csv(OUTPUT_CSV, index=False)
            logger.info(f"🔄 已自动备份前 {completed_count + len(pkl_results)} 个变异的原始数据 (PKL + CSV)...")

    # 最终保存
    with open(OUTPUT_PKL, "wb") as f:
        pickle.dump(pkl_results, f, protocol=pickle.HIGHEST_PROTOCOL)
    pd.DataFrame(csv_results).to_csv(OUTPUT_CSV, index=False)
    logger.info(f"============== 全部完成！==============")

if __name__ == "__main__":
    main()

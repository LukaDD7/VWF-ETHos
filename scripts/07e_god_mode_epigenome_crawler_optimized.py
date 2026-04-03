#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07e_god_mode_epigenome_crawler_optimized.py - 优化版

优化内容：
1. 从CSV读取任务清单（313KB）替代184GB pickle，秒级启动
2. 流式写入CSV（csv.writer替代pandas），低内存占用
3. 实时进度显示

输入：../results/07_VCF_AlphaGenome_Results_Backup.csv（轻量级任务清单）
输出：../results/07e_GodMode_Epigenome_Peaks.csv
"""

import os
import logging
import sys
import time
import json
import csv
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Optional, Set, Dict, List
import random

import numpy as np

# 导入 AlphaGenome
from alphagenome.data import genome
from alphagenome.models import dna_client

# ==================== 配置区域 ====================
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"
LOGS_DIR = BASE_DIR / "logs"

# 输入文件 - 使用轻量级CSV替代184GB pickle
GOLD_STANDARD_CSV = RESULTS_DIR / "07_VCF_AlphaGenome_Results_Backup.csv"
TASK_LIST_JSON = RESULTS_DIR / "07e_task_list.json"

# 输出文件
OUTPUT_CSV = RESULTS_DIR / "07e_GodMode_Epigenome_Peaks.csv"
CHECKPOINT_JSONL = RESULTS_DIR / "07e_godmode_checkpoint.jsonl"

# API 配置
API_KEY = "AIzaSyC25qItDVDM6jJJFqNDSvEyhR4P56E0CtY"
ONTOLOGY_TERMS = None  # 所有组织/细胞系

# 网络容错配置
MAX_RETRIES = 5
INITIAL_BACKOFF = 5
BACKOFF_FACTOR = 2
MAX_BACKOFF = 300
REQUEST_DELAY = 2

# 确保目录存在
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# 日志配置
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS_DIR / "07e_godmode_log.txt", mode='a'),
    ]
)
logger = logging.getLogger(__name__)


# ==================== 数据结构 ====================
@dataclass
class GodModeResult:
    """God Mode 单个变异结果"""
    chromosome: str
    position: int
    ref: str
    alt: str

    # 核心预测值
    rna_delta: float = 0.0
    splice_delta: float = 0.0

    # 补充的5种模态
    atac_max_delta: float = 0.0
    atac_track_idx: int = -1
    histone_max_delta: float = 0.0
    histone_track_idx: int = -1
    cage_max_delta: float = 0.0
    cage_track_idx: int = -1
    procap_max_delta: float = 0.0
    procap_track_idx: int = -1
    contact_max_delta: float = 0.0
    contact_track_idx: int = -1

    # 数据来源标记
    source: str = "AllTissues"
    status: str = "Pending"  # Success, Failed, Skipped
    error_message: str = ""


# ==================== 轻量级任务清单 ====================

def generate_task_list_from_csv():
    """从CSV生成轻量级任务清单（一次性）"""
    if TASK_LIST_JSON.exists():
        logger.info(f"✅ 任务清单已存在: {TASK_LIST_JSON}")
        return

    if not GOLD_STANDARD_CSV.exists():
        logger.error(f"❌ 金标准CSV不存在: {GOLD_STANDARD_CSV}")
        sys.exit(1)

    logger.info(f"📂 从CSV生成任务清单: {GOLD_STANDARD_CSV}")

    tasks = []
    with open(GOLD_STANDARD_CSV, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            tasks.append({
                'chromosome': row.get('hg38_chr', ''),
                'position': int(row.get('hg38_pos', 0)),
                'ref': row.get('ref', ''),
                'alt': row.get('alt', ''),
            })

    with open(TASK_LIST_JSON, 'w', encoding='utf-8') as f:
        json.dump(tasks, f, ensure_ascii=False, indent=2)

    logger.info(f"✅ 任务清单生成完成: {len(tasks)} 个变异")


def load_task_list() -> List[Dict]:
    """加载轻量级任务清单（秒级）"""
    # 如果不存在，先生成
    if not TASK_LIST_JSON.exists():
        generate_task_list_from_csv()

    with open(TASK_LIST_JSON, 'r', encoding='utf-8') as f:
        tasks = json.load(f)

    logger.info(f"✅ 加载任务清单: {len(tasks)} 个变异")
    return tasks


# ==================== CSV 流式操作 ====================

CSV_HEADER = [
    'Chromosome', 'Position', 'Ref', 'Alt',
    'RNA_Delta_HUVEC', 'Splice_Delta_HUVEC',
    'ATAC_Max_Delta', 'ATAC_Track_Idx',
    'Histone_Max_Delta', 'Histone_Track_Idx',
    'CAGE_Max_Delta', 'CAGE_Track_Idx',
    'PROCAP_Max_Delta', 'PROCAP_Track_Idx',
    'Contact_Max_Delta', 'Contact_Track_Idx',
    'Data_Source', 'God_Mode_Status'
]


def load_completed_variants() -> Set[str]:
    """从CSV加载已完成的变异（轻量级）"""
    completed = set()
    if not OUTPUT_CSV.exists():
        logger.info("📂 CSV不存在，从头开始")
        return completed

    try:
        with open(OUTPUT_CSV, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # 只读取简单状态
                status = row.get('God_Mode_Status', '')
                if status in ['Success', 'Failed', 'Pending']:
                    key = f"{row.get('Chromosome', '')}:{row.get('Position', '')}"
                    completed.add(key)
        logger.info(f"✅ 从CSV加载了 {len(completed)} 个已完成变异")
    except Exception as e:
        logger.warning(f"⚠️ 读取CSV失败: {e}，将从头开始")

    return completed


def init_csv_file():
    """初始化CSV文件（如果不存在）"""
    if not OUTPUT_CSV.exists():
        with open(OUTPUT_CSV, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(CSV_HEADER)
        logger.info(f"✅ 创建新CSV: {OUTPUT_CSV}")


def append_result_to_csv_stream(result: GodModeResult):
    """流式追加单条记录到CSV（低内存）"""
    row = [
        result.chromosome,
        result.position,
        result.ref,
        result.alt,
        result.rna_delta,
        result.splice_delta,
        result.atac_max_delta,
        result.atac_track_idx,
        result.histone_max_delta,
        result.histone_track_idx,
        result.cage_max_delta,
        result.cage_track_idx,
        result.procap_max_delta,
        result.procap_track_idx,
        result.contact_max_delta,
        result.contact_track_idx,
        result.source,
        result.status if result.status in ['Success', 'Failed', 'Pending'] else 'Failed'
    ]

    with open(OUTPUT_CSV, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(row)


def append_result_to_checkpoint(result: GodModeResult):
    """追加单条记录到 checkpoint JSONL"""
    record = {
        'chromosome': result.chromosome,
        'position': result.position,
        'ref': result.ref,
        'alt': result.alt,
        'rna_delta': result.rna_delta,
        'splice_delta': result.splice_delta,
        'atac_max_delta': result.atac_max_delta,
        'histone_max_delta': result.histone_max_delta,
        'cage_max_delta': result.cage_max_delta,
        'procap_max_delta': result.procap_max_delta,
        'contact_max_delta': result.contact_max_delta,
        'source': result.source,
        'status': result.status,
        'error_message': str(result.error_message)[:200] if result.error_message else ""
    }
    with open(CHECKPOINT_JSONL, 'a', encoding='utf-8') as f:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')


# ==================== 网络容错机制 ====================

def exponential_backoff_retry(func, max_retries=MAX_RETRIES,
                               initial_backoff=INITIAL_BACKOFF,
                               backoff_factor=BACKOFF_FACTOR,
                               max_backoff=MAX_BACKOFF):
    """指数退避重试装饰器"""
    for attempt in range(max_retries):
        try:
            return func()
        except Exception as e:
            error_str = str(e)

            is_network_error = any([
                "UNAVAILABLE" in error_str,
                "Timeout occurred" in error_str,
                "FD Shutdown" in error_str,
                "failed to connect" in error_str,
                "Connection timed out" in error_str,
                "Socket closed" in error_str,
                "Broken pipe" in error_str,
                "Connection reset by peer" in error_str,
                "tcp handshaker shutdown" in error_str,
                "RST_STREAM" in error_str,
                "recvmsg" in error_str,
            ])

            if not is_network_error:
                raise

            if attempt < max_retries - 1:
                backoff = min(initial_backoff * (backoff_factor ** attempt), max_backoff)
                jitter = random.uniform(0, 0.1 * backoff)
                sleep_time = backoff + jitter

                logger.warning(f"🌐 网络错误 (尝试 {attempt + 1}/{max_retries}): {e}")
                logger.info(f"⏳ 指数退避: 等待 {sleep_time:.1f} 秒后重试...")
                time.sleep(sleep_time)
            else:
                logger.error(f"❌ 最大重试次数用尽 ({max_retries})")
                raise

    return None


# ==================== God Mode 提取逻辑 ====================

def get_max_delta_with_idx(ref_track, alt_track, expected_tracks=None):
    """获取最大差值及其索引（带索引补偿）"""
    if ref_track is not None and alt_track is not None:
        ref_val = ref_track.values
        alt_val = alt_track.values
        if ref_val.size > 0 and alt_val.size > 0:
            delta = np.abs(alt_val - ref_val)
            actual_tracks = delta.shape[-1] if delta.ndim > 1 else 1

            if expected_tracks and actual_tracks < expected_tracks:
                effective_tracks = min(expected_tracks, actual_tracks)
                delta_effective = delta[..., :effective_tracks] if delta.ndim > 1 else delta
                max_idx = int(np.argmax(delta_effective))
                max_val = float(delta_effective.flat[max_idx])
                return max_val, max_idx
            else:
                max_idx = int(np.argmax(delta))
                max_val = float(delta.flat[max_idx] if delta.ndim == 1 else delta[..., max_idx].max())
                return max_val, max_idx
    return 0.0, -1


def extract_all_max_deltas(outputs) -> dict:
    """从 AlphaGenome 输出中提取所有模态的最大差值"""
    scores = {}

    # RNA & Splice (核心参考)
    scores['rna'], _ = get_max_delta_with_idx(
        outputs.reference.rna_seq, outputs.alternate.rna_seq, expected_tracks=1)
    scores['splice'], _ = get_max_delta_with_idx(
        outputs.reference.splice_sites, outputs.alternate.splice_sites, expected_tracks=1)

    # 补充的5种模态（全组织数据）
    scores['atac'], scores['atac_idx'] = get_max_delta_with_idx(
        outputs.reference.atac, outputs.alternate.atac, expected_tracks=200)
    scores['histone'], scores['histone_idx'] = get_max_delta_with_idx(
        outputs.reference.chip_histone, outputs.alternate.chip_histone, expected_tracks=800)
    scores['cage'], scores['cage_idx'] = get_max_delta_with_idx(
        outputs.reference.cage, outputs.alternate.cage, expected_tracks=300)
    scores['procap'], scores['procap_idx'] = get_max_delta_with_idx(
        outputs.reference.procap, outputs.alternate.procap, expected_tracks=20)
    scores['contact'], scores['contact_idx'] = get_max_delta_with_idx(
        outputs.reference.contact_maps, outputs.alternate.contact_maps, expected_tracks=20)

    return scores


def run_prediction_with_retry(model, dna_client, genome, task: dict, half_len=524288):
    """带重试机制的预测函数"""
    chrom = task['chromosome']
    if not str(chrom).startswith('chr'):
        chrom = f"chr{chrom}"

    pos = int(task['position'])
    ref = task['ref']
    alt = task['alt']

    result = GodModeResult(
        chromosome=chrom,
        position=pos,
        ref=ref,
        alt=alt
    )

    def do_prediction():
        interval_start = max(1, pos - half_len)
        interval_end = pos + half_len

        interval = genome.Interval(chromosome=chrom, start=interval_start, end=interval_end)
        variant = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)

        requested_outputs = [
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.SPLICE_SITES,
            dna_client.OutputType.ATAC,
            dna_client.OutputType.CHIP_HISTONE,
            dna_client.OutputType.CAGE,
            dna_client.OutputType.PROCAP,
            dna_client.OutputType.CONTACT_MAPS
        ]

        outputs = model.predict_variant(
            interval=interval,
            variant=variant,
            ontology_terms=ONTOLOGY_TERMS,
            requested_outputs=requested_outputs
        )
        return outputs

    try:
        outputs = exponential_backoff_retry(do_prediction)
        scores = extract_all_max_deltas(outputs)

        result.rna_delta = scores['rna']
        result.splice_delta = scores['splice']
        result.atac_max_delta = scores['atac']
        result.atac_track_idx = scores['atac_idx']
        result.histone_max_delta = scores['histone']
        result.histone_track_idx = scores['histone_idx']
        result.cage_max_delta = scores['cage']
        result.cage_track_idx = scores['cage_idx']
        result.procap_max_delta = scores['procap']
        result.procap_track_idx = scores['procap_idx']
        result.contact_max_delta = scores['contact']
        result.contact_track_idx = scores['contact_idx']
        result.status = "Success"

    except Exception as e:
        result.status = "Failed"
        result.error_message = str(e)

    return result


# ==================== 主函数 ====================

def main():
    """主函数 - 优化版"""
    logger.info("=" * 70)
    logger.info("🚀 07e God Mode Epigenome Crawler - 优化版")
    logger.info("   优化: CSV任务清单 + 流式写入")
    logger.info("=" * 70)

    # 1. 加载轻量级任务清单（秒级）
    tasks = load_task_list()
    total_tasks = len(tasks)

    # 2. 初始化CSV
    init_csv_file()

    # 3. 加载已完成变异
    completed_variants = load_completed_variants()
    logger.info(f"✅ 断点续传: 已跳过 {len(completed_variants)} 个已完成变异")

    # 4. 初始化 AlphaGenome 客户端
    logger.info("🔌 连接 AlphaGenome God Mode...")
    try:
        model = dna_client.create(API_KEY)
        logger.info("✅ AlphaGenome API 连接成功!")
    except Exception as e:
        logger.error(f"❌ API 连接失败: {e}")
        sys.exit(1)

    # 5. 主循环 - 流式处理
    success_count = len(completed_variants)
    fail_count = 0
    skip_count = 0

    for idx, task in enumerate(tasks, start=1):
        chrom = task['chromosome']
        pos = task['position']
        variant_key = f"{chrom}:{pos}"

        # 断点续传
        if variant_key in completed_variants:
            skip_count += 1
            if skip_count <= 5 or skip_count % 100 == 0:
                logger.info(f"⏭️  [{idx}/{total_tasks}] 跳过已完成: {chrom}:{pos}")
            continue

        current_num = success_count + fail_count + 1
        logger.info(f"🧬 [{current_num}/{total_tasks}] 正在爬取 {chrom}:{pos} 的表观暗物质...")

        try:
            # 运行预测
            result = run_prediction_with_retry(model, dna_client, genome, task)

            # 流式写入（低内存）
            append_result_to_csv_stream(result)
            append_result_to_checkpoint(result)

            if result.status == "Success":
                success_count += 1
                logger.info(f"✅ [{current_num}/{total_tasks}] 成功: {chrom}:{pos} | "
                           f"RNA={result.rna_delta:.2f}, Splice={result.splice_delta:.4f}, "
                           f"ATAC={result.atac_max_delta:.1f}, Histone={result.histone_max_delta:.1f}")
            else:
                fail_count += 1
                logger.warning(f"❌ [{current_num}/{total_tasks}] 失败: {chrom}:{pos}")

        except Exception as e:
            fail_count += 1
            logger.error(f"❌ [{current_num}/{total_tasks}] 提取失败: {e}")

            # 记录失败
            failed_result = GodModeResult(
                chromosome=str(chrom),
                position=int(pos),
                ref=str(task.get('ref', '')),
                alt=str(task.get('alt', '')),
                status="Failed",
                error_message=str(e)
            )
            append_result_to_csv_stream(failed_result)
            append_result_to_checkpoint(failed_result)

        # 请求间隔
        time.sleep(REQUEST_DELAY)

    # 6. 完成统计
    logger.info("=" * 70)
    logger.info("🎉 God Mode 爬取完成!")
    logger.info(f"   总计: {total_tasks}")
    logger.info(f"   成功: {success_count}")
    logger.info(f"   失败: {fail_count}")
    logger.info(f"   跳过: {skip_count}")
    logger.info(f"   输出: {OUTPUT_CSV}")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()

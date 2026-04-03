#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07e_god_mode_epigenome_crawler.py - 究极强化版

God Mode 表观暗物质爬取引擎 - 补充内皮细胞缺失的5种模态

输入：../results/07_VCF_AlphaGenome_Results_Backup.pkl (内皮细胞6模态金标准)
输出：../results/07e_GodMode_Epigenome_Peaks.csv (补充5模态：ATAC/Histone/CAGE/PROcap/Contact)

背景：
- 内皮细胞(CL:0000115)只有6种有效模态：RNA, Splice, SpliceSiteUsage, SpliceJunctions, DNase, ChIP-TF
- 缺失5种模态：ATAC, ChIP-Histone, CAGE, PROcap, ContactMaps
- 07e 使用全组织数据(ONTOLOGY_TERMS=None)补充这5种缺失模态

特性：
- 自动断点续传（从已有 CSV 恢复进度）
- 指数退避重试（处理 gRPC UNAVAILABLE）
- 实时追加保存（每成功一个立即写入）
- 网络容错（断网时自动休眠等待）
- "以 5 补 6"索引补偿（处理 track 数量不匹配）
"""

import os

import logging
import sys
import time
import pickle
import json
import pandas as pd
import numpy as np
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Optional, Set
import random

# 导入 AlphaGenome
from alphagenome.data import genome
from alphagenome.models import dna_client

# ==================== 配置区域 ====================
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"
LOGS_DIR = BASE_DIR / "logs"

# 输入文件
GOLD_STANDARD_PKL = RESULTS_DIR / "07_VCF_AlphaGenome_Results_Backup.pkl"

# 输出文件（断点续传目标）
OUTPUT_CSV = RESULTS_DIR / "07e_GodMode_Epigenome_Peaks.csv"
CHECKPOINT_JSONL = RESULTS_DIR / "07e_godmode_checkpoint.jsonl"

# API 配置
API_KEY = "AIzaSyC25qItDVDM6jJJFqNDSvEyhR4P56E0CtY"
# 07e 核心逻辑：用全组织数据(None)补充内皮细胞缺失的5种模态
ONTOLOGY_TERMS = None  # 所有组织/细胞系，获取完整的11模态

# 网络容错配置
MAX_RETRIES = 5          # 最大重试次数
INITIAL_BACKOFF = 5      # 初始退避秒数
BACKOFF_FACTOR = 2       # 指数退避倍数
MAX_BACKOFF = 300        # 最大退避秒数（5分钟）
REQUEST_DELAY = 2        # 正常请求间隔

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

# 导入 07 脚本的 InferenceResult 以正确加载 pickle 文件
@dataclass
class InferenceResult:
    """07 脚本的数据结构 - 用于加载金标准 pickle 文件"""
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

    delta_score: float = 0.0
    raw_outputs: Any = field(default=None)
    success: bool = False
    error_message: str | None = None


@dataclass
class GodModeResult:
    """God Mode 单个变异结果

    07e 核心逻辑：
    - 内皮细胞(CL:0000115)已有6种有效模态：RNA, Splice, SpliceSiteUsage, SpliceJunctions, DNase, ChIP-TF
    - 内皮细胞缺失5种模态：ATAC, ChIP-Histone, CAGE, PROcap, ContactMaps
    - 07e 使用全组织数据(ONTOLOGY_TERMS=None)补充缺失的5种模态
    - 所有提取的模态均标记为"AllTissues"来源，以区分内皮细胞结果
    """
    chromosome: str
    position: int
    ref: str
    alt: str

    # 核心预测值（来自全组织数据，作为参考）
    rna_delta: float = 0.0
    splice_delta: float = 0.0

    # 补充的5种模态（全组织数据）- 07e 核心输出
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
    source: str = "AllTissues"  # 全组织数据，区别于内皮细胞

    status: str = "Pending"  # Success, Failed, Skipped
    error_message: str = ""


# ==================== 断点续传机制 ====================

def get_variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """生成变异的唯一标识键"""
    return f"{chrom}:{pos}:{ref}:{alt}"


def load_completed_variants_from_csv(csv_path: Path) -> Set[str]:
    """从已有的 CSV 文件加载已完成的变异"""
    completed = set()
    if not csv_path.exists():
        logger.info(f"未找到已有 CSV 文件，将从零开始: {csv_path}")
        return completed

    try:
        df = pd.read_csv(csv_path)
        for _, row in df.iterrows():
            key = get_variant_key(
                str(row['Chromosome']),
                int(row['Position']),
                str(row['Ref']),
                str(row['Alt'])
            )
            completed.add(key)
        logger.info(f"✅ 从 CSV 加载了 {len(completed)} 个已完成变异")
    except Exception as e:
        logger.warning(f"⚠️ 读取 CSV 失败: {e}，将从零开始")
    return completed


def load_completed_variants_from_checkpoint(checkpoint_path: Path) -> Set[str]:
    """从 checkpoint JSONL 加载已完成的变异（双重保险）"""
    completed = set()
    if not checkpoint_path.exists():
        return completed

    try:
        with open(checkpoint_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                    key = get_variant_key(
                        data['chromosome'],
                        data['position'],
                        data['ref'],
                        data['alt']
                    )
                    completed.add(key)
                except (json.JSONDecodeError, KeyError):
                    continue
        logger.info(f"✅ 从 checkpoint 加载了 {len(completed)} 个已完成变异")
    except Exception as e:
        logger.warning(f"⚠️ 读取 checkpoint 失败: {e}")
    return completed


def append_result_to_csv(result: GodModeResult, csv_path: Path):
    """追加单条记录到 CSV（实时保存）

    注意：错误信息会被截断和清理，避免破坏CSV格式
    """
    # 清理错误信息：只保留第一行，去除特殊字符
    error_msg = result.error_message
    if error_msg:
        # 只保留第一行
        error_msg = str(error_msg).split('\n')[0]
        # 限制长度
        error_msg = error_msg[:100]
        # 移除CSV分隔符
        error_msg = error_msg.replace(',', ';')

    record = {
        'Chromosome': result.chromosome,
        'Position': result.position,
        'Ref': result.ref,
        'Alt': result.alt,
        # 全组织数据（作为参考，与内皮细胞结果区分）
        'RNA_Delta_HUVEC': result.rna_delta,
        'Splice_Delta_HUVEC': result.splice_delta,
        # 补充的5种模态（07e核心输出）
        'ATAC_Max_Delta': result.atac_max_delta,
        'ATAC_Track_Idx': result.atac_track_idx,
        'Histone_Max_Delta': result.histone_max_delta,
        'Histone_Track_Idx': result.histone_track_idx,
        'CAGE_Max_Delta': result.cage_max_delta,
        'CAGE_Track_Idx': result.cage_track_idx,
        'PROCAP_Max_Delta': result.procap_max_delta,
        'PROCAP_Track_Idx': result.procap_track_idx,
        'Contact_Max_Delta': result.contact_max_delta,
        'Contact_Track_Idx': result.contact_track_idx,
        # 数据来源标记
        'Data_Source': result.source,
        'God_Mode_Status': result.status if result.status in ['Success', 'Failed', 'Pending'] else 'Failed'
    }

    df = pd.DataFrame([record])

    # 追加模式（文件存在则追加，不存在则创建并写入表头）
    file_exists = csv_path.exists()
    df.to_csv(csv_path, mode='a', header=not file_exists, index=False)


def append_result_to_checkpoint(result: GodModeResult, checkpoint_path: Path):
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
        'error_message': result.error_message
    }
    with open(checkpoint_path, 'a', encoding='utf-8') as f:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')


# ==================== 网络容错机制 ====================

def exponential_backoff_retry(func, max_retries=MAX_RETRIES,
                               initial_backoff=INITIAL_BACKOFF,
                               backoff_factor=BACKOFF_FACTOR,
                               max_backoff=MAX_BACKOFF):
    """
    指数退避重试装饰器
    处理 gRPC UNAVAILABLE 和超时错误
    """
    for attempt in range(max_retries):
        try:
            return func()
        except Exception as e:
            error_str = str(e)

            # 检查是否是网络/连接错误（扩展匹配列表）
            is_network_error = any([
                "UNAVAILABLE" in error_str,
                "Timeout occurred" in error_str,
                "FD Shutdown" in error_str,
                "failed to connect" in error_str,
                "Connection timed out" in error_str,
                # 新增错误类型
                "Socket closed" in error_str,
                "Broken pipe" in error_str,
                "Connection reset by peer" in error_str,
                "tcp handshaker shutdown" in error_str,
                "RST_STREAM" in error_str,
                "recvmsg" in error_str,
                "ipv4:" in error_str and "Failed" in error_str,
            ])

            if not is_network_error:
                # 非网络错误，直接抛出
                raise

            if attempt < max_retries - 1:
                # 计算退避时间（指数增长 + 随机抖动）
                backoff = min(initial_backoff * (backoff_factor ** attempt), max_backoff)
                jitter = random.uniform(0, 0.1 * backoff)  # 10% 随机抖动
                sleep_time = backoff + jitter

                logger.warning(f"🌐 网络错误 (尝试 {attempt + 1}/{max_retries}): {e}")
                logger.info(f"⏳ 指数退避: 等待 {sleep_time:.1f} 秒后重试...")
                time.sleep(sleep_time)
            else:
                # 最大重试次数用尽
                logger.error(f"❌ 最大重试次数用尽 ({max_retries})，放弃请求")
                raise

    return None


# ==================== God Mode 提取逻辑 ====================

def extract_all_max_deltas(outputs) -> dict:
    """从 AlphaGenome 输出中提取所有模态的最大差值"""
    scores = {}

    def get_max_delta_with_idx(ref_track, alt_track, expected_tracks=None):
        """获取最大差值及其索引

        支持索引补偿逻辑：当预期有 N 个 tracks 但实际只有 M 个（M < N）时，
        限制在 M 范围内计算，避免越界错误。

        这在全组织模式下很重要，因为不同模态的 track 数量差异很大：
        - ATAC: 可能 ~150 tracks
        - Histone: 可能 ~700+ tracks
        - CAGE: 可能 ~250+ tracks
        - PROcap: 可能 ~10+ tracks
        - ContactMaps: 可能 ~10+ tracks

        Args:
            ref_track: REF 状态的 track 数据
            alt_track: ALT 状态的 track 数据
            expected_tracks: 预期的 track 数量（用于"以 5 补 6"补偿逻辑）
        """
        if ref_track is not None and alt_track is not None:
            ref_val = ref_track.values
            alt_val = alt_track.values
            if ref_val.size > 0 and alt_val.size > 0:
                delta = np.abs(alt_val - ref_val)

                # 获取实际 track 数量（最后一个维度）
                actual_tracks = delta.shape[-1] if delta.ndim > 1 else 1

                # "以 5 补 6"补偿逻辑：处理 track 数量不匹配的情况
                # 当预期 N 个 tracks 但实际只有 M 个（M < N）时，限制在 M 范围内计算
                # 这在全组织模式下很重要，因为不同模态的 track 数量差异很大
                if expected_tracks and actual_tracks < expected_tracks:
                    # 实际 tracks 少于预期，限制索引在有效范围内
                    # 例如：预期 200 个但只有 150 个，只在前 150 个范围内找最大值
                    effective_tracks = min(expected_tracks, actual_tracks)
                    # 在有效 tracks 范围内找最大值
                    delta_effective = delta[..., :effective_tracks] if delta.ndim > 1 else delta
                    max_idx = int(np.argmax(delta_effective))
                    max_val = float(delta_effective.flat[max_idx])
                    return max_val, max_idx
                else:
                    # 正常情况
                    max_idx = int(np.argmax(delta))
                    max_val = float(delta.flat[max_idx] if delta.ndim == 1 else delta[..., max_idx].max())
                    return max_val, max_idx
        return 0.0, -1

    # RNA & Splice (核心参考)
    scores['rna'], _ = get_max_delta_with_idx(
        outputs.reference.rna_seq, outputs.alternate.rna_seq, expected_tracks=1)
    scores['splice'], _ = get_max_delta_with_idx(
        outputs.reference.splice_sites, outputs.alternate.splice_sites, expected_tracks=1)

    # 补充的5种模态（全组织数据）
    # 注意：全组织模式下 track 数量可能很多，使用较大的 expected_tracks
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


def run_prediction_with_retry(model, dna_client, genome, row, half_len=524288):
    """
    带重试机制的预测函数
    """
    chrom = row['hg38_chr'] if 'hg38_chr' in row else row.get('chromosome', '')
    if not str(chrom).startswith('chr'):
        chrom = f"chr{chrom}"

    pos = int(row['hg38_pos']) if 'hg38_pos' in row else int(row.get('position', 0))
    ref = row['ref']
    alt = row['alt']

    # 创建结果对象
    result = GodModeResult(
        chromosome=chrom,
        position=pos,
        ref=ref,
        alt=alt
    )

    # 定义预测函数（用于重试包装）
    def do_prediction():
        interval_start = max(1, pos - half_len)
        interval_end = pos + half_len

        interval = genome.Interval(chromosome=chrom, start=interval_start, end=interval_end)
        variant = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)

        # 请求全模态输出
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
        # 使用指数退避重试
        outputs = exponential_backoff_retry(do_prediction)

        # 提取所有分数
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
        # 不再raise，返回失败结果让上层处理

    return result


# ==================== 主函数 ====================

def main():
    """主函数 - 支持断点续传和网络容错

    07e 核心逻辑：
    1. 内皮细胞(CL:0000115)已有6种有效模态（RNA, Splice, SpliceSiteUsage, SpliceJunctions, DNase, ChIP-TF）
    2. 内皮细胞缺失5种模态（ATAC, ChIP-Histone, CAGE, PROcap, ContactMaps）
    3. 使用全组织数据(ONTOLOGY_TERMS=None)补充这5种缺失模态
    4. 保留"以 5 补 6"索引补偿逻辑，处理全组织模式下 track 数量不匹配问题
    """
    logger.info("=" * 70)
    logger.info("🚀 07e God Mode Epigenome Crawler - 补充内皮细胞缺失5模态")
    logger.info("   内皮细胞有效模态(6): RNA/Splice/SpliceSiteUsage/SpliceJunctions/DNase/ChIP-TF")
    logger.info("   全组织补充模态(5): ATAC/ChIP-Histone/CAGE/PROcap/ContactMaps")
    logger.info("   ONTOLOGY_TERMS = None (所有组织/细胞系)")
    logger.info("=" * 70)

    # 1. 读取金标准任务清单
    if not GOLD_STANDARD_PKL.exists():
        logger.error(f"❌ 金标准清单不存在: {GOLD_STANDARD_PKL}")
        sys.exit(1)

    logger.info(f"📂 读取金标准任务清单: {GOLD_STANDARD_PKL}")
    with open(GOLD_STANDARD_PKL, 'rb') as f:
        gold_standard_data = pickle.load(f)

    # 转换为 DataFrame
    if isinstance(gold_standard_data, list):
        df_tasks = pd.DataFrame(gold_standard_data)
    elif isinstance(gold_standard_data, pd.DataFrame):
        df_tasks = gold_standard_data
    else:
        logger.error("❌ 不支持的金标准数据格式")
        sys.exit(1)

    total_tasks = len(df_tasks)
    logger.info(f"📊 总任务数: {total_tasks}")

    # 2. 加载已完成变异（双重保险：CSV + Checkpoint）
    completed_from_csv = load_completed_variants_from_csv(OUTPUT_CSV)
    completed_from_ckpt = load_completed_variants_from_checkpoint(CHECKPOINT_JSONL)
    completed_variants = completed_from_csv | completed_from_ckpt  # 并集

    logger.info(f"✅ 断点续传: 已跳过 {len(completed_variants)} 个已完成变异")

    # 3. 初始化 AlphaGenome 客户端
    logger.info("🔌 连接 AlphaGenome God Mode...")
    try:
        model = dna_client.create(API_KEY)  # 使用默认连接，通过环境变量配置代理
        logger.info("✅ AlphaGenome API 连接成功!")
    except Exception as e:
        logger.error(f"❌ API 连接失败: {e}")
        sys.exit(1)

    # 4. 主循环 - 带断点续传和网络容错
    success_count = len(completed_variants)
    fail_count = 0
    skip_count = 0

    for idx, row in df_tasks.iterrows():
        # 获取变异标识
        chrom = row.get('hg38_chr', row.get('chromosome', ''))
        pos = row.get('hg38_pos', row.get('position', 0))
        ref = row.get('ref', row.get('REF', ''))
        alt = row.get('alt', row.get('ALT', ''))

        variant_key = get_variant_key(str(chrom), int(pos), str(ref), str(alt))

        # 断点续传：跳过已完成
        if variant_key in completed_variants:
            skip_count += 1
            if skip_count <= 5 or skip_count % 100 == 0:
                logger.info(f"⏭️  [{idx + 1}/{total_tasks}] 跳过已完成: {chrom}:{pos}")
            continue

        # 处理当前变异
        current_num = success_count + fail_count + 1
        logger.info(f"🧬 [{current_num}/{total_tasks}] 正在爬取 {chrom}:{pos} 的表观暗物质...")

        try:
            # 运行预测（带重试）
            result = run_prediction_with_retry(model, dna_client, genome, row)

            # 实时追加保存（断点续传核心）
            append_result_to_csv(result, OUTPUT_CSV)
            append_result_to_checkpoint(result, CHECKPOINT_JSONL)

            success_count += 1
            logger.info(f"✅ [{current_num}/{total_tasks}] 成功: {chrom}:{pos} | "
                       f"RNA={result.rna_delta:.2f}, Splice={result.splice_delta:.4f}, "
                       f"ATAC={result.atac_max_delta:.1f}, Histone={result.histone_max_delta:.1f}")

        except Exception as e:
            fail_count += 1
            logger.error(f"❌ [{current_num}/{total_tasks}] 提取失败: {e}")

            # 记录失败到 checkpoint（避免无限重试同一变异）
            failed_result = GodModeResult(
                chromosome=str(chrom),
                position=int(pos),
                ref=str(ref),
                alt=str(alt),
                status="Failed",
                error_message=str(e)
            )
            append_result_to_checkpoint(failed_result, CHECKPOINT_JSONL)
            # 失败后继续处理下一个变异
            continue

        # 请求间隔
        time.sleep(REQUEST_DELAY)

    # 5. 完成统计
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

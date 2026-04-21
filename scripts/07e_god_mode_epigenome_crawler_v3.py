#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07e_god_mode_epigenome_crawler_v3.py - 代理流量保护版

修复内容：
1. 【代理流量检测】通过错误信息识别代理限制 (ipv4:10.5.36.27:7897)
2. 【智能暂停】检测到代理问题时，自动暂停等待代理恢复（5-30分钟）
3. 【渐进式间隔】根据代理状态动态调整请求间隔
4. 【健康检查】暂停后先验证代理是否恢复，再继续请求
5. 【批量控制】限制每批请求数量，避免流量暴增

使用方法：
    # 正常启动（自动检测代理状态）
    python3 scripts/07e_god_mode_epigenome_crawler_v3.py

    # 强制从头开始（跳过checkpoint）
    rm results/07e_godmode_checkpoint.jsonl
    rm results/07e_GodMode_Epigenome_Peaks.csv
    python3 scripts/07e_god_mode_epigenome_crawler_v3.py

输入：../results/07_VCF_AlphaGenome_Results_Backup.csv
输出：../results/07e_GodMode_Epigenome_Peaks.csv
"""

import os
import logging
import sys
import time
import json
import csv
import subprocess
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Optional, Set, Dict, List
from datetime import datetime, timedelta

import numpy as np

# ==================== 设置代理环境变量 ====================
os.environ['HTTP_PROXY'] = 'http://10.5.36.27:7897'
os.environ['HTTPS_PROXY'] = 'http://10.5.36.27:7897'
os.environ['http_proxy'] = 'http://10.5.36.27:7897'
os.environ['https_proxy'] = 'http://10.5.36.27:7897'
logger.info(f"🌐 代理已设置: http://10.5.36.27:7897")

# 导入 AlphaGenome
from alphagenome.data import genome
from alphagenome.models import dna_client

# ==================== 配置区域 ====================
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"
LOGS_DIR = BASE_DIR / "logs"

# 输入文件 - 使用主结果CSV作为任务清单
GOLD_STANDARD_CSV = RESULTS_DIR / "07_VCF_AlphaGenome_Results.csv"
TASK_LIST_JSON = RESULTS_DIR / "07e_task_list.json"

# 输出文件
OUTPUT_CSV = RESULTS_DIR / "07e_GodMode_Epigenome_Peaks.csv"
CHECKPOINT_JSONL = RESULTS_DIR / "07e_godmode_checkpoint.jsonl"
STATUS_FILE = RESULTS_DIR / "07e_proxy_status.json"  # 代理状态文件

# API 配置
API_KEY = "AIzaSyC25qItDVDM6jJJFqNDSvEyhR4P56E0CtY"
ONTOLOGY_TERMS = None  # 所有组织/细胞系
API_ADDRESS = "gdmscience.googleapis.com:443"  # AlphaGenome 实际端点

# ==================== 代理流量保护配置 ====================

# 代理地址（用于识别代理相关错误）
PROXY_ADDRESS = "10.5.36.27:7897"

# 智能暂停配置
INITIAL_REQUEST_DELAY = 2    # 初始请求间隔（秒）
PROXY_COOLDOWN_INITIAL = 300  # 初始冷却时间 5 分钟
PROXY_COOLDOWN_MAX = 1800    # 最大冷却时间 30 分钟
PROXY_COOLDOWN_MULTIPLIER = 2  # 连续触发冷却时的指数增长

# 连续失败阈值
CONSECUTIVE_PROXY_ERRORS_THRESHOLD = 2  # 连续代理错误2次就触发冷却
MAX_CONSECUTIVE_SUCCESS = 50  # 成功50次后可以重置冷却时间

# 批量控制
BATCH_SIZE = 100  # 每批处理100个后检查代理状态

# 当前代理冷却状态（运行时）
current_cooldown_until: Optional[float] = None
consecutive_proxy_errors: int = 0
current_request_delay: float = INITIAL_REQUEST_DELAY
success_since_last_cooldown: int = 0


# ==================== 日志配置 ====================

LOG_FILE = LOGS_DIR / "07e_godmode_v3_log.txt"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOG_FILE, mode='a'),
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

    rna_delta: float = 0.0
    splice_delta: float = 0.0
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

    source: str = "AllTissues"
    status: str = "Pending"
    error_message: str = ""


# ==================== 代理状态管理 ====================

def load_proxy_status() -> dict:
    """从文件加载代理状态"""
    if not STATUS_FILE.exists():
        return {
            "cooldown_until": None,
            "last_cooldown_time": None,
            "total_cooldowns": 0,
            "last_proxy_error": None
        }

    try:
        with open(STATUS_FILE, 'r') as f:
            return json.load(f)
    except Exception:
        return {}


def save_proxy_status(status: dict):
    """保存代理状态到文件"""
    try:
        with open(STATUS_FILE, 'w') as f:
            json.dump(status, f, indent=2, default=str)
    except Exception as e:
        logger.warning(f"无法保存代理状态: {e}")


def is_proxy_error(error_str: str) -> bool:
    """判断错误是否是代理相关错误"""
    proxy_indicators = [
        PROXY_ADDRESS,           # 包含代理地址
        "RST_STREAM",             # RST_STREAM 通常是代理重置
        "Socket closed",          # 代理关闭连接
        "Timeout occurred",       # 代理超时
        "Failed to connect to remote host",  # 代理连接失败
        "grpc_status:13",         # gRPC 内部错误（通常是代理问题）
        "grpc_status:14",         # gRPC UNAVAILABLE（代理不可用）
    ]
    return any(indicator in error_str for indicator in proxy_indicators)


def check_proxy_cooldown():
    """检查是否处于代理冷却期"""
    global current_cooldown_until

    if current_cooldown_until is None:
        return False

    now = time.time()
    if now < current_cooldown_until:
        remaining = current_cooldown_until - now
        logger.info(f"⏸️  代理冷却中... 剩余 {remaining/60:.1f} 分钟")
        return True

    # 冷却结束
    logger.info(f"✅ 代理冷却结束，开始恢复请求")
    current_cooldown_until = None
    return False


def trigger_proxy_cooldown(reason: str):
    """触发代理冷却"""
    global current_cooldown_until, consecutive_proxy_errors, current_request_delay

    consecutive_proxy_errors += 1

    # 读取历史冷却记录
    status = load_proxy_status()

    # 计算冷却时间（指数退避）
    cooldown_time = min(
        PROXY_COOLDOWN_INITIAL * (PROXY_COOLDOWN_MULTIPLIER ** (status.get("total_cooldowns", 0))),
        PROXY_COOLDOWN_MAX
    )

    current_cooldown_until = time.time() + cooldown_time

    # 更新状态
    status["cooldown_until"] = datetime.fromtimestamp(current_cooldown_until).isoformat()
    status["last_cooldown_time"] = datetime.now().isoformat()
    status["last_proxy_error"] = reason
    status["total_cooldowns"] = status.get("total_cooldowns", 0) + 1
    save_proxy_status(status)

    # 增加请求间隔
    current_request_delay = min(current_request_delay * 1.5, 30)

    logger.warning(f"🚨 代理流量保护触发！")
    logger.warning(f"   原因: {reason}")
    logger.warning(f"   冷却时间: {cooldown_time/60:.1f} 分钟")
    logger.warning(f"   累计触发次数: {status['total_cooldowns']}")
    logger.warning(f"   新请求间隔: {current_request_delay:.1f} 秒")

    # 发送桌面通知
    send_desktop_notification(
        "🚨 07e 代理保护触发",
        f"原因: {reason}\n冷却时间: {cooldown_time/60:.1f} 分钟\n累计触发: {status['total_cooldowns']} 次",
        urgency="critical"
    )


def reset_proxy_success():
    """代理请求成功后重置计数器"""
    global consecutive_proxy_errors, success_since_last_cooldown, current_request_delay

    success_since_last_cooldown += 1
    consecutive_proxy_errors = 0

    # 渐进恢复请求间隔
    if success_since_last_cooldown > MAX_CONSECUTIVE_SUCCESS:
        current_request_delay = max(current_request_delay * 0.9, INITIAL_REQUEST_DELAY)
        success_since_last_cooldown = 0


# ==================== 桌面通知 ====================

def send_desktop_notification(title, message, urgency="critical"):
    """发送桌面通知"""
    try:
        subprocess.run(
            ["notify-send", "-u", urgency, "-t", "5000", title, message],
            check=True, capture_output=True, timeout=5
        )
        return True
    except Exception:
        return False


# ==================== 任务清单 ====================

def load_task_list() -> List[Dict]:
    """加载任务清单"""
    if not TASK_LIST_JSON.exists():
        logger.info(f"📂 从CSV生成任务清单...")
        if not GOLD_STANDARD_CSV.exists():
            logger.error(f"❌ 金标准CSV不存在: {GOLD_STANDARD_CSV}")
            sys.exit(1)

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
    else:
        with open(TASK_LIST_JSON, 'r', encoding='utf-8') as f:
            tasks = json.load(f)
        logger.info(f"✅ 加载任务清单: {len(tasks)} 个变异")

    return tasks


# ==================== CSV 操作 ====================

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
    """从CSV加载已完成的变异"""
    completed = set()
    if not OUTPUT_CSV.exists():
        return completed

    try:
        with open(OUTPUT_CSV, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = f"{row.get('Chromosome', '')}:{row.get('Position', '')}"
                completed.add(key)
        logger.info(f"✅ 从CSV加载了 {len(completed)} 个已完成变异")
    except Exception as e:
        logger.warning(f"⚠️ 读取CSV失败: {e}")
    return completed


def init_csv_file():
    """初始化CSV文件"""
    if not OUTPUT_CSV.exists():
        with open(OUTPUT_CSV, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(CSV_HEADER)


def append_result_to_csv(result: GodModeResult):
    """追加单条记录到CSV"""
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
        result.status
    ]
    with open(OUTPUT_CSV, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(row)


def append_result_to_checkpoint(result: GodModeResult):
    """追加单条记录到 checkpoint"""
    record = {
        'chromosome': result.chromosome,
        'position': result.position,
        'ref': result.ref,
        'alt': result.alt,
        'status': result.status,
        'error_message': str(result.error_message)[:200] if result.error_message else "",
        'timestamp': datetime.now().isoformat()
    }
    with open(CHECKPOINT_JSONL, 'a', encoding='utf-8') as f:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')


# ==================== 预测逻辑 ====================

def get_max_delta_with_idx(ref_track, alt_track, expected_tracks=None):
    """获取最大差值及其索引"""
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

    scores['rna'], _ = get_max_delta_with_idx(
        outputs.reference.rna_seq, outputs.alternate.rna_seq, expected_tracks=1)
    scores['splice'], _ = get_max_delta_with_idx(
        outputs.reference.splice_sites, outputs.alternate.splice_sites, expected_tracks=1)
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


def run_prediction(model, dna_client, genome, task: dict, half_len=524288) -> GodModeResult:
    """执行单次预测"""
    global consecutive_proxy_errors

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

    try:
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

        # 成功后重置代理错误计数
        reset_proxy_success()

    except Exception as e:
        error_str = str(e)
        result.status = "Failed"
        result.error_message = error_str

        # 检查是否是代理错误
        if is_proxy_error(error_str):
            consecutive_proxy_errors += 1
            logger.warning(f"🌐 代理错误 (连续 {consecutive_proxy_errors}/{CONSECUTIVE_PROXY_ERRORS_THRESHOLD}): {error_str[:100]}")

            # 达到阈值触发冷却
            if consecutive_proxy_errors >= CONSECUTIVE_PROXY_ERRORS_THRESHOLD:
                trigger_proxy_cooldown(error_str[:150])
        else:
            # 非代理错误，不计入连续代理错误
            consecutive_proxy_errors = 0

    return result


# ==================== 代理健康检查 ====================

def wait_for_proxy_recovery(model, dna_client, genome, max_attempts=3):
    """等待代理恢复并进行健康检查"""
    global current_cooldown_until

    logger.info("🔍 开始代理健康检查...")

    for attempt in range(max_attempts):
        # 等待冷却时间
        if current_cooldown_until:
            now = time.time()
            if now < current_cooldown_until:
                wait_time = current_cooldown_until - now
                logger.info(f"⏳ 等待冷却结束... ({wait_time/60:.1f} 分钟剩余)")
                time.sleep(min(wait_time, 60))  # 最多等1分钟检查一次
                current_cooldown_until = None

        # 尝试一个简单的请求
        logger.info(f"🧪 健康检查尝试 {attempt + 1}/{max_attempts}...")
        try:
            test_interval = genome.Interval(chromosome="chr12", start=6000000, end=6100000)
            test_variant = genome.Variant(chromosome="chr12", position=6050000, reference_bases="A", alternate_bases="G")

            model.predict_variant(
                interval=test_interval,
                variant=test_variant,
                ontology_terms=ONTOLOGY_TERMS,
                requested_outputs=[dna_client.OutputType.RNA_SEQ]
            )
            logger.info("✅ 健康检查成功！代理已恢复")
            return True

        except Exception as e:
            logger.warning(f"❌ 健康检查失败: {str(e)[:80]}")
            if attempt < max_attempts - 1:
                wait_time = 60 * (attempt + 1)  # 递增等待
                logger.info(f"⏳ 等待 {wait_time} 秒后重试...")
                time.sleep(wait_time)

    logger.error("❌ 健康检查全部失败，代理可能仍不可用")
    return False


# ==================== 主函数 ====================

def main():
    global current_request_delay, current_cooldown_until

    logger.info("=" * 70)
    logger.info("🚀 07e God Mode Epigenome Crawler - 代理流量保护版 v3")
    logger.info(f"   代理地址: {PROXY_ADDRESS}")
    logger.info(f"   初始请求间隔: {INITIAL_REQUEST_DELAY}s")
    logger.info(f"   代理冷却初始时间: {PROXY_COOLDOWN_INITIAL/60}分钟")
    logger.info("=" * 70)

    # 加载代理状态
    status = load_proxy_status()
    if status.get("cooldown_until"):
        try:
            cooldown_ts = datetime.fromisoformat(status["cooldown_until"]).timestamp()
            if time.time() < cooldown_ts:
                current_cooldown_until = cooldown_ts
                logger.info(f"📍 检测到未完成的冷却期: {status['cooldown_until']}")
        except Exception:
            pass

    # 1. 加载任务清单
    tasks = load_task_list()
    total_tasks = len(tasks)

    # 2. 初始化CSV
    init_csv_file()

    # 3. 加载已完成变异
    completed_variants = load_completed_variants()
    pending_tasks = [t for t in tasks if f"{t['chromosome']}:{t['position']}" not in completed_variants]
    logger.info(f"📊 总任务: {total_tasks}, 已完成: {len(completed_variants)}, 待处理: {len(pending_tasks)}")

    # 4. 初始化 AlphaGenome 客户端
    logger.info("🔌 连接 AlphaGenome API...")
    try:
        model = dna_client.create(API_KEY)
        logger.info("✅ API 连接成功!")
    except Exception as e:
        logger.error(f"❌ API 连接失败: {e}")
        sys.exit(1)

    # 5. 主循环
    success_count = len(completed_variants)
    fail_count = 0
    batch_count = 0

    for idx, task in enumerate(pending_tasks, start=1):
        chrom = task['chromosome']
        pos = task['position']
        variant_key = f"{chrom}:{pos}"

        current_num = success_count + fail_count + 1
        total_processed = len(completed_variants) + success_count + fail_count

        # 检查代理冷却
        if check_proxy_cooldown():
            # 执行健康检查
            if not wait_for_proxy_recovery(model, dna_client, genome):
                logger.warning("⚠️ 代理仍不可用，继续等待...")

        logger.info(f"🧬 [{current_num}/{len(pending_tasks)}] 处理 {chrom}:{pos}...")

        # 执行预测
        result = run_prediction(model, dna_client, genome, task)

        # 保存结果
        append_result_to_csv(result)
        append_result_to_checkpoint(result)

        if result.status == "Success":
            success_count += 1
            logger.info(f"✅ 成功 | RNA={result.rna_delta:.2f}, ATAC={result.atac_max_delta:.1f}, Histone={result.histone_max_delta:.1f}")
        else:
            fail_count += 1
            logger.warning(f"❌ 失败 | {result.error_message[:80] if result.error_message else 'Unknown'}")

        batch_count += 1

        # 每批次检查
        if batch_count >= BATCH_SIZE:
            batch_count = 0
            logger.info(f"📊 批次完成: {success_count} 成功, {fail_count} 失败, {total_processed}/{total_tasks} 总进度")

            # 检查是否需要冷却
            if consecutive_proxy_errors >= CONSECUTIVE_PROXY_ERRORS_THRESHOLD:
                trigger_proxy_cooldown("批次检查触发冷却")

        # 请求间隔
        time.sleep(current_request_delay)

    # 6. 完成统计
    final_status = load_proxy_status()
    final_status["total_success"] = success_count
    final_status["total_failed"] = fail_count
    final_status["completed_at"] = datetime.now().isoformat()
    save_proxy_status(final_status)

    logger.info("=" * 70)
    logger.info("🎉 07e God Mode 爬取完成!")
    logger.info(f"   总任务: {total_tasks}")
    logger.info(f"   成功: {success_count}")
    logger.info(f"   失败: {fail_count}")
    logger.info(f"   代理冷却触发: {final_status.get('total_cooldowns', 0)} 次")
    logger.info(f"   输出: {OUTPUT_CSV}")
    logger.info("=" * 70)

    send_desktop_notification(
        "✅ 07e 完成",
        f"成功: {success_count}, 失败: {fail_count}\n代理冷却: {final_status.get('total_cooldowns', 0)} 次",
        urgency="normal"
    )


if __name__ == "__main__":
    main()
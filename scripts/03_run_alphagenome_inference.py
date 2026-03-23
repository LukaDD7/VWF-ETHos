#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本 3：AlphaGenome API 推理

目标：调用 DeepMind AlphaGenome API 进行变异效应预测。

预测内容：
  - RNA_SEQ：基因表达预测
  - SPLICE_SITES：剪接位点预测

输入：../data/02_ready_for_inference.csv
输出：
  - ../results/03_inference_results.pkl (完整预测结果，包含 raw_outputs)
  - ../results/03_inference_results.csv (摘要统计)

依赖：alphagenome (已安装在当前虚拟环境)
"""

import logging
import pickle
import sys
import time
import json
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# ==================== 配置区域 ====================
# 基础目录配置
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
LOGS_DIR = BASE_DIR / "logs"

INPUT_FILE = DATA_DIR / "02_ready_for_inference.csv"
OUTPUT_PKL = RESULTS_DIR / "03_inference_results.pkl"
OUTPUT_CSV = RESULTS_DIR / "03_inference_results.csv"
# 断点续传：临时保存文件（每次成功立即追加）
CHECKPOINT_JSONL = RESULTS_DIR / "03_inference_checkpoint.jsonl"

API_KEY = "AIzaSyC25qItDVDM6jJJFqNDSvEyhR4P56E0CtY"

# 请求的输出类型
REQUESTED_OUTPUTS = ["RNA_SEQ", "SPLICE_SITES"]

# 本体论术语：CL:0000115 = 内皮细胞 (VWF 的主要表达细胞)
# 注意：UBERON:0002363（血管内皮）没有对应的 RNA-seq tracks，故改用 CL:0000115
ONTOLOGY_TERMS = ['CL:0000115']

# API 调用配置
REQUEST_DELAY = 2  # 秒

# 限制处理的变异数量（None 表示全部）
MAX_VARIANTS = None

# 确保输出和日志目录存在
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# 日志配置
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS_DIR / "03_inference.log", mode="w"),
    ],
)
logger = logging.getLogger(__name__)


@dataclass
class InferenceResult:
    """存储单个变异的推理结果。"""
    chromosome: str
    position: int
    ref: str
    alt: str
    interval_start: int
    interval_end: int
    clinvar_classification: str
    consequences: str

    # 预测结果摘要
    rna_seq_delta_max: float | None = None
    rna_seq_delta_sum: float | None = None
    splice_site_delta_max: float | None = None
    splice_site_delta_sum: float | None = None
    delta_score: float = 0.0

    # 完整输出对象（用于后续可视化）
    raw_outputs: Any = field(default=None)  # 保存完整的 outputs 对象
    success: bool = False
    error_message: str | None = None


def initialize_alphagenome_client():
    """初始化 AlphaGenome 客户端。"""
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    logger.info("正在连接 AlphaGenome API...")
    model = dna_client.create(API_KEY)
    logger.info("AlphaGenome API 连接成功！")
    return dna_client, genome, model


def get_output_type_enum(dna_client, output_name: str):
    """获取 OutputType 枚举值。"""
    output_map = {
        "RNA_SEQ": dna_client.OutputType.RNA_SEQ,
        "SPLICE_SITES": dna_client.OutputType.SPLICE_SITES,
        "DNASE": dna_client.OutputType.DNASE,
        "ATAC": dna_client.OutputType.ATAC,
    }
    return output_map.get(output_name)


def run_prediction(dna_client, genome, model, row: pd.Series) -> InferenceResult:
    """对单个变异运行 AlphaGenome 预测，保存完整输出。"""
    result = InferenceResult(
        chromosome=row["chromosome"],
        position=row["hg38_POS"],
        ref=row["REF"],
        alt=row["ALT"],
        interval_start=int(row["interval_start"]),
        interval_end=int(row["interval_end"]),
        clinvar_classification=row.get("INFO_clinvar_classification_base", "N/A"),
        consequences=row.get("INFO_consequences_base", "N/A"),
        raw_outputs=None,
    )

    try:
        interval = genome.Interval(
            chromosome=result.chromosome,
            start=result.interval_start,
            end=result.interval_end,
        )
        variant = genome.Variant(
            chromosome=result.chromosome,
            position=result.position,
            reference_bases=result.ref,
            alternate_bases=result.alt,
        )

        requested_outputs = [
            get_output_type_enum(dna_client, out) for out in REQUESTED_OUTPUTS
        ]
        requested_outputs = [o for o in requested_outputs if o is not None]

        # 调用 API
        outputs = model.predict_variant(
            interval=interval,
            variant=variant,
            ontology_terms=ONTOLOGY_TERMS,
            requested_outputs=requested_outputs,
        )

        # 保存完整输出对象
        result.raw_outputs = outputs

        # 提取 RNA-seq 结果
        if outputs.reference.rna_seq is not None and outputs.reference.rna_seq.values.size > 0:
            ref_rna = outputs.reference.rna_seq.values
            alt_rna = outputs.alternate.rna_seq.values
            result.rna_seq_delta_max = float(np.max(np.abs(alt_rna - ref_rna)))
            result.rna_seq_delta_sum = float(np.sum(np.abs(alt_rna - ref_rna)))

        # 提取剪接位点结果
        if outputs.reference.splice_sites is not None and outputs.reference.splice_sites.values.size > 0:
            ref_splice = outputs.reference.splice_sites.values
            alt_splice = outputs.alternate.splice_sites.values
            result.splice_site_delta_max = float(np.max(np.abs(alt_splice - ref_splice)))
            result.splice_site_delta_sum = float(np.sum(np.abs(alt_splice - ref_splice)))

        # 计算综合 Delta Score
        scores = [s for s in [result.rna_seq_delta_max, result.splice_site_delta_max] if s is not None]
        result.delta_score = max(scores) if scores else 0.0

        result.success = True
        logger.info(f"✓ 预测成功：{result.chromosome}:{result.position} {result.ref}>{result.alt}")

    except Exception as e:
        result.error_message = str(e)
        logger.warning(f"✗ 预测失败：{result.chromosome}:{result.position} {result.ref}>{result.alt} - {e}")

    return result


def save_results(results: list[InferenceResult], output_pkl: Path, output_csv: Path) -> pd.DataFrame:
    """保存推理结果。"""
    # 保存完整结果（pickle）- 包含 raw_outputs
    with open(output_pkl, "wb") as f:
        pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)
    logger.info(f"已保存完整结果到：{output_pkl}")

    # 保存摘要（CSV）
    summary_data = []
    for r in results:
        summary_data.append({
            "chromosome": r.chromosome,
            "position": r.position,
            "ref": r.ref,
            "alt": r.alt,
            "interval_start": r.interval_start,
            "interval_end": r.interval_end,
            "clinvar_classification": r.clinvar_classification,
            "consequences": r.consequences,
            "rna_seq_delta_max": r.rna_seq_delta_max,
            "rna_seq_delta_sum": r.rna_seq_delta_sum,
            "splice_site_delta_max": r.splice_site_delta_max,
            "splice_site_delta_sum": r.splice_site_delta_sum,
            "delta_score": r.delta_score,
            "success": r.success,
            "error_message": r.error_message,
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_csv, index=False)
    logger.info(f"已保存摘要结果到：{output_csv}")
    return summary_df


# ==================== 断点续传（Checkpoint）机制 ====================

def get_variant_key(row: pd.Series) -> str:
    """生成变异的唯一标识键。"""
    return f"{row['chromosome']}:{row['hg38_POS']}:{row['REF']}:{row['ALT']}"


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
                    # 从 checkpoint 记录中提取变异标识
                    key = f"{data['chromosome']}:{data['position']}:{data['ref']}:{data['alt']}"
                    completed.add(key)
                except (json.JSONDecodeError, KeyError) as e:
                    logger.warning(f"跳过损坏的 checkpoint 行: {e}")
        logger.info(f"已加载 {len(completed)} 个已完成的变异记录")
    except Exception as e:
        logger.error(f"读取 checkpoint 失败: {e}")
    return completed


def append_result_to_checkpoint(result: InferenceResult, checkpoint_path: Path):
    """立即追加单个结果到 checkpoint 文件（JSONL 格式）。"""
    # 将 result 转换为可 JSON 序列化的字典
    record = {
        "chromosome": result.chromosome,
        "position": result.position,
        "ref": result.ref,
        "alt": result.alt,
        "interval_start": result.interval_start,
        "interval_end": result.interval_end,
        "clinvar_classification": result.clinvar_classification,
        "consequences": result.consequences,
        "rna_seq_delta_max": result.rna_seq_delta_max,
        "rna_seq_delta_sum": result.rna_seq_delta_sum,
        "splice_site_delta_max": result.splice_site_delta_max,
        "splice_site_delta_sum": result.splice_site_delta_sum,
        "delta_score": result.delta_score,
        "success": result.success,
        "error_message": result.error_message,
    }
    # 追加写入 JSONL（每行一个 JSON 对象）
    with open(checkpoint_path, 'a', encoding='utf-8') as f:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')


def load_existing_results_from_checkpoint(checkpoint_path: Path) -> list[InferenceResult]:
    """从 checkpoint 文件加载所有已有结果（用于最终 PKL 保存）。"""
    results = []
    if not checkpoint_path.exists():
        return results

    try:
        with open(checkpoint_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                    result = InferenceResult(
                        chromosome=data['chromosome'],
                        position=data['position'],
                        ref=data['ref'],
                        alt=data['alt'],
                        interval_start=data['interval_start'],
                        interval_end=data['interval_end'],
                        clinvar_classification=data['clinvar_classification'],
                        consequences=data['consequences'],
                        rna_seq_delta_max=data.get('rna_seq_delta_max'),
                        rna_seq_delta_sum=data.get('rna_seq_delta_sum'),
                        splice_site_delta_max=data.get('splice_site_delta_max'),
                        splice_site_delta_sum=data.get('splice_site_delta_sum'),
                        delta_score=data.get('delta_score', 0.0),
                        success=data.get('success', False),
                        error_message=data.get('error_message'),
                    )
                    results.append(result)
                except (json.JSONDecodeError, KeyError) as e:
                    logger.warning(f"跳过损坏的 checkpoint 记录: {e}")
        logger.info(f"从 checkpoint 恢复了 {len(results)} 个结果")
    except Exception as e:
        logger.error(f"恢复 checkpoint 失败: {e}")
    return results


def main():
    """主函数 - 支持断点续传。"""
    logger.info("=" * 60)
    logger.info("开始执行脚本 3：AlphaGenome API 推理")
    logger.info("=" * 60)

    input_path = INPUT_FILE
    if not input_path.exists():
        logger.error(f"输入文件不存在：{input_path}")
        sys.exit(1)

    df = pd.read_csv(input_path)
    logger.info(f"加载完成：{len(df)} 个变异待预测")

    # 优先确保包含致病变异（阳性对照）
    if "INFO_clinvar_classification_base" in df.columns:
        pathogenic = df[df["INFO_clinvar_classification_base"] == "Pathogenic"]
        if len(pathogenic) > 0:
            logger.info(f"检测到 {len(pathogenic)} 个 Pathogenic 变异，将优先处理")
            remaining = df[~df.index.isin(pathogenic.index)]
            df = pd.concat([pathogenic, remaining]).reset_index(drop=True)

    # 限制处理数量
    if MAX_VARIANTS:
        df = df.head(MAX_VARIANTS)
        logger.info(f"限制处理前 {MAX_VARIANTS} 个变异")

    # ========== 断点续传：加载已完成的变异 ==========
    completed_variants = load_completed_variants(CHECKPOINT_JSONL)
    logger.info(f"已跳过 {len(completed_variants)} 个已完成的变异")

    # 筛选出未完成的变异
    df['variant_key'] = df.apply(get_variant_key, axis=1)
    df_pending = df[~df['variant_key'].isin(completed_variants)].copy()
    df_pending = df_pending.drop(columns=['variant_key'])

    pending_count = len(df_pending)
    completed_count = len(completed_variants)
    total_count = len(df)

    logger.info(f"待处理: {pending_count}/{total_count} (已完成: {completed_count})")

    if pending_count == 0:
        logger.info("所有变异都已处理完成！正在生成最终输出文件...")
        existing_results = load_existing_results_from_checkpoint(CHECKPOINT_JSONL)
        save_results(existing_results, OUTPUT_PKL, OUTPUT_CSV)
        logger.info("=" * 60)
        logger.info("断点续传：全部数据已导出！")
        logger.info("=" * 60)
        return

    # 初始化客户端
    try:
        dna_client, genome, model = initialize_alphagenome_client()
    except Exception as e:
        logger.error(f"初始化失败：{e}")
        sys.exit(1)

    # 加载已有结果（用于最终合并）
    results = load_existing_results_from_checkpoint(CHECKPOINT_JSONL)

    # 运行预测（只处理未完成的变异）
    for idx, row in df_pending.iterrows():
        actual_idx = completed_count + idx + 1  # 实际进度计数
        logger.info(f"处理变异 {actual_idx}/{total_count} (当前批次: {idx + 1}/{pending_count})")
        result = run_prediction(dna_client, genome, model, row)
        results.append(result)

        # ========== 断点续传：立即追加保存 ==========
        append_result_to_checkpoint(result, CHECKPOINT_JSONL)

        if idx < len(df_pending) - 1:
            time.sleep(REQUEST_DELAY)

    # 最终保存：合并所有结果生成 PKL 和 CSV
    logger.info("生成最终输出文件...")
    summary_df = save_results(results, OUTPUT_PKL, OUTPUT_CSV)

    # 可选：删除 checkpoint 文件（如果全部成功）
    if CHECKPOINT_JSONL.exists():
        logger.info(f"保留 checkpoint 文件用于审计: {CHECKPOINT_JSONL}")

    # 输出统计
    success_count = sum(1 for r in results if r.success)
    logger.info("=" * 60)
    logger.info("推理完成！")
    logger.info(f"成功：{success_count}/{len(results)}")

    if success_count > 0:
        top = summary_df.nlargest(5, "delta_score")
        logger.info("\nDelta Score 最高的 5 个变异:")
        for _, row in top.iterrows():
            logger.info(
                f"  {row['chromosome']}:{row['position']} {row['ref']}>{row['alt']} - "
                f"Delta: {row['delta_score']:.4f} | ClinVar: {row['clinvar_classification']}"
            )

    logger.info("=" * 60)


if __name__ == "__main__":
    main()

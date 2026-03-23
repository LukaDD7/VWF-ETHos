#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本 2：预处理与坐标转换（hg19 -> hg38）

目标：
  1. 使用 pyliftover 将 hg19/GRCh37 坐标精确转换为 hg38/GRCh38
  2. 构建 1Mb 上下文区间（以变异位置为中心，左右各扩展 500kb）

注意：
  - AlphaGenome 使用的染色体格式为 'chr1', 'chr2', ..., 'chrX', 'chrY', 'chrM'
  - 使用本地 hg19ToHg38.over.chain.gz 文件进行精确转换
  - 转换失败的变异将被严格丢弃

输入：../data/01_filtered_targets.csv
输出：../data/02_ready_for_inference.csv

依赖安装：
  conda install -c bioconda pyliftover
"""

import logging
import sys
from pathlib import Path

import pandas as pd

# ==================== 配置区域 ====================
# 基础目录配置
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
DATA_DIR = BASE_DIR / "data"
LOGS_DIR = BASE_DIR / "logs"

INPUT_FILE = DATA_DIR / "01_filtered_targets.csv"
OUTPUT_FILE = DATA_DIR / "02_ready_for_inference.csv"

# AlphaGenome 支持的序列长度
SEQUENCE_LENGTH_1MB = 2**20  # 1,048,576 bp
HALF_LENGTH = SEQUENCE_LENGTH_1MB // 2  # 524,288 bp

# 染色体长度参考（hg38，用于边界检查）
CHROMOSOME_LENGTHS_HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415, "chrM": 16569,
}

# Chain 文件配置 - 使用本地已下载的文件（在 data/ 目录中）
CHAIN_FILE_PATH = DATA_DIR / "hg19ToHg38.over.chain.gz"

# 确保日志目录存在
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# 日志配置
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS_DIR / "02_preprocess.log", mode="w"),
    ],
)
logger = logging.getLogger(__name__)


def get_chain_file() -> Path:
    """获取 hg19->hg38 chain 文件（使用本地文件）。"""
    if not CHAIN_FILE_PATH.exists():
        logger.error(f"Chain 文件不存在：{CHAIN_FILE_PATH}")
        logger.error("请确保 hg19ToHg38.over.chain.gz 文件在 data/ 目录下")
        raise FileNotFoundError(f"Chain file not found: {CHAIN_FILE_PATH}")

    logger.info(f"使用本地 Chain 文件：{CHAIN_FILE_PATH}")
    return CHAIN_FILE_PATH


def load_filtered_data(filepath: Path) -> pd.DataFrame:
    """加载筛选后的数据。"""
    logger.info(f"正在加载文件：{filepath}")
    df = pd.read_csv(filepath)
    logger.info(f"加载完成：{len(df)} 行变异")
    return df


def convert_chromosome_name(chrom) -> str:
    """转换染色体名称为 AlphaGenome 格式。"""
    chrom_str = str(chrom)
    return chrom_str if chrom_str.startswith("chr") else f"chr{chrom_str}"


def liftover_hg19_to_hg38(df: pd.DataFrame, chain_path: Path) -> pd.DataFrame:
    """
    使用 pyliftover 将坐标从 hg19 精确转换到 hg38。
    转换失败的变异将被严格丢弃。
    """
    from pyliftover import LiftOver

    logger.info(f"初始化 LiftOver 工具（hg19 -> hg38）...")
    # 直接使用本地 chain 文件路径
    lo = LiftOver(str(chain_path))

    hg38_chroms, hg38_positions, failures, reasons = [], [], [], []

    for idx, row in df.iterrows():
        chrom_hg19, pos_hg19 = str(row["CHROM"]), int(row["POS"])
        try:
            result = lo.convert_coordinate(f"chr{chrom_hg19}", pos_hg19 - 1)
            if result:
                new_chrom = result[0][0].replace("chr", "")
                new_pos = result[0][1] + 1  # 0-based -> 1-based
                hg38_chroms.append(new_chrom)
                hg38_positions.append(new_pos)
                failures.append(False)
                reasons.append("")
            else:
                logger.warning(f"行 {idx}: 坐标无 hg38 映射 - {chrom_hg19}:{pos_hg19}")
                hg38_chroms.append(None)
                hg38_positions.append(None)
                failures.append(True)
                reasons.append("no mapping")
        except Exception as e:
            logger.warning(f"行 {idx}: 转换异常 - {chrom_hg19}:{pos_hg19}: {e}")
            hg38_chroms.append(None)
            hg38_positions.append(None)
            failures.append(True)
            reasons.append(str(e))

    df["hg38_CHROM"], df["hg38_POS"] = hg38_chroms, hg38_positions
    failed = sum(failures)

    if failed > 0:
        logger.info(f"严格剔除 {failed} 个转换失败的变异")
        df = df[~df["hg38_CHROM"].isna()].reset_index(drop=True)

    logger.info(f"坐标转换完成：{len(df)} 个有效变异")
    return df


def build_interval(df: pd.DataFrame) -> pd.DataFrame:
    """构建 AlphaGenome 所需的 1Mb 区间。"""
    logger.info("正在构建 1Mb 上下文区间...")
    df["chromosome"] = df["hg38_CHROM"].apply(convert_chromosome_name)
    df["interval_start"] = df["hg38_POS"] - HALF_LENGTH
    df["interval_end"] = df["hg38_POS"] + HALF_LENGTH

    # 边界调整
    for idx, row in df.iterrows():
        chrom, max_len = row["chromosome"], CHROMOSOME_LENGTHS_HG38.get(row["chromosome"], 0)
        if max_len and row["interval_start"] < 1:
            df.loc[idx, "interval_start"] = 1
            df.loc[idx, "interval_end"] = min(SEQUENCE_LENGTH_1MB, max_len)
        elif max_len and row["interval_end"] > max_len:
            df.loc[idx, "interval_end"] = max_len
            df.loc[idx, "interval_start"] = max(1, max_len - SEQUENCE_LENGTH_1MB + 1)

    df["interval_width"] = df["interval_end"] - df["interval_start"] + 1
    logger.info(f"区间构建完成。平均长度：{df['interval_width'].mean():.0f} bp")
    return df


def save_ready_data(df: pd.DataFrame, output_path: Path) -> None:
    """保存准备好的数据。"""
    cols = ["chromosome", "hg38_POS", "interval_start", "interval_end", "interval_width",
            "CHROM", "POS", "REF", "ALT", "INFO_clinvar_classification_base", "INFO_consequences_base"]
    existing = [c for c in cols if c in df.columns]
    df.to_csv(output_path, index=False, columns=existing)
    logger.info(f"已保存数据到：{output_path}")


def main():
    """主函数。"""
    logger.info("=" * 60)
    logger.info("开始执行脚本 2：预处理与坐标转换")
    logger.info("=" * 60)

    # 获取本地 chain 文件
    try:
        chain_path = get_chain_file()
    except FileNotFoundError as e:
        logger.error(f"无法获取 chain 文件：{e}")
        sys.exit(1)

    # 加载数据
    input_path = INPUT_FILE
    if not input_path.exists():
        logger.error(f"输入文件不存在：{input_path}")
        sys.exit(1)

    df = load_filtered_data(input_path)

    # 坐标转换
    try:
        df = liftover_hg19_to_hg38(df, chain_path)
    except ImportError:
        logger.error("未安装 pyliftover。请运行：conda install -c bioconda pyliftover")
        sys.exit(1)
    except Exception as e:
        logger.error(f"坐标转换失败：{e}")
        sys.exit(1)

    # 构建区间
    df = build_interval(df)

    # 保存
    save_ready_data(df, OUTPUT_FILE)

    logger.info("=" * 60)
    logger.info("预处理完成！")
    logger.info(f"输出文件：{OUTPUT_FILE}")
    logger.info(f"准备的变异数量：{len(df)}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()

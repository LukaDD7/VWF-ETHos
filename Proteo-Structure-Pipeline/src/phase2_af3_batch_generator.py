#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 2: AlphaFold3 Batch JSON 生成器 (官方 Batch JSON 格式)

严格遵循 AlphaFold Server 官方 Batch JSON 语法：
1. 最外层必须是列表（Array [ ]）
2. 每个任务对象包含 name 和 sequences（含 proteinChain）
3. 输出纯 JSON 文件（非 ZIP），每天限额 30 Jobs

AlphaFold Server 官方 Batch JSON 格式：
[
  {
    "name": "VWF_WT",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "MIPARF...",
          "count": 1
        }
      }
    ]
  },
  {
    "name": "VWF_G1531D",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "MIPARF...D...",
          "count": 1
        }
      }
    ]
  }
]

输出文件：
- af3_batch_01.json       # 每文件 30 个任务（默认）
- af3_batch_02.json
- ...
- af3_jobs_manifest.json  # 任务清单
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import pandas as pd
from Bio import SeqIO

# 配置日志
log_dir = Path('../logs')
log_dir.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_dir / 'phase2_batch_generator.log')
    ]
)
logger = logging.getLogger(__name__)


class AF3BatchJob:
    """AlphaFold3 单个任务配置"""

    def __init__(self, name: str, sequence: str, count: int = 1):
        self.name = name
        self.sequence = sequence
        self.count = count

    def to_dict(self) -> Dict:
        """
        转换为 AlphaFold Server 官方 Batch JSON 格式
        注意：不包含 modelSeeds，官方示例中没有此字段
        """
        return {
            "name": self.name,
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": self.sequence,
                        "count": self.count
                    }
                }
            ]
        }


class AF3BatchGenerator:
    """
    AlphaFold3 批量 JSON 生成器（官方 Batch 格式）

    输出 JSON 数组格式，每个文件包含多个任务
    适配 AF3 官网每日限额（默认 30 个任务/文件）
    """

    def __init__(
        self,
        csv_path: str,
        wt_fasta: str,
        output_dir: str,
        limit: Optional[int] = None,
        chunk_size: int = 30,
        completed_jobs: Optional[List[str]] = None
    ):
        self.csv_path = Path(csv_path)
        self.wt_fasta = Path(wt_fasta)
        self.output_dir = Path(output_dir)
        self.limit = limit
        self.chunk_size = chunk_size
        self.completed_jobs = set(completed_jobs) if completed_jobs else set()

        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.wt_sequence = None
        self.mutations_df = None
        self.batch_jobs: List[AF3BatchJob] = []

    def load_data(self):
        """加载变异列表和 WT 序列"""
        logger.info("="*60)
        logger.info("Phase 2: AlphaFold3 Batch JSON 生成器 (官方格式)")
        logger.info("="*60)

        # 读取 CSV
        logger.info(f"[1/5] 加载变异列表: {self.csv_path}")
        self.mutations_df = pd.read_csv(self.csv_path)

        original_count = len(self.mutations_df)

        if self.limit:
            self.mutations_df = self.mutations_df.head(self.limit)
            logger.info(f"      限制模式: 只取前 {self.limit} 个变异")

        logger.info(f"      原始变异数: {original_count}")

        # 强制去重：基于物理结构（位置+野生型氨基酸+突变型氨基酸）
        self.mutations_df = self.mutations_df.drop_duplicates(
            subset=['AA_Position', 'WT_AA_1', 'Mut_AA_1'],
            keep='first'  # 保留第一次出现的记录
        )

        duplicates_removed = original_count - len(self.mutations_df)
        unique_count = len(self.mutations_df)

        if duplicates_removed > 0:
            logger.info(f"      ⚠️ 发现并移除 {duplicates_removed} 个重复物理结构")
            logger.info(f"      ✅ 剩余 {unique_count} 个唯一突变任务")
        else:
            logger.info(f"      ✅ 无重复结构，共 {unique_count} 个唯一突变任务")

        # 读取 WT 序列
        logger.info(f"[2/5] 加载 WT 序列: {self.wt_fasta}")
        wt_record = SeqIO.read(self.wt_fasta, "fasta")
        self.wt_sequence = str(wt_record.seq)
        logger.info(f"      序列长度: {len(self.wt_sequence)} aa")

        # 验证序列长度（AlphaFold3 限制）
        if len(self.wt_sequence) > 5000:
            logger.warning(f"      ⚠️ 警告: 序列长度超过 5000 aa，AlphaFold3 可能不支持")

    def create_batch_jobs(self):
        """创建所有批量任务"""
        logger.info("[3/5] 创建 AlphaFold3 Batch Jobs...")

        # 1. 添加 WT 任务（只需一次）
        wt_job = AF3BatchJob(
            name="VWF_WT",
            sequence=self.wt_sequence,
            count=1
        )
        self.batch_jobs.append(wt_job)
        logger.info(f"      [1] VWF_WT (WT序列, {len(self.wt_sequence)} aa)")

        # 2. 为每个变异创建任务
        for idx, row in self.mutations_df.iterrows():
            pos = int(row['AA_Position'])
            wt_aa = row['WT_AA_1']
            mut_aa = row['Mut_AA_1']

            # 创建突变序列
            mut_sequence = (
                self.wt_sequence[:pos-1] +
                mut_aa +
                self.wt_sequence[pos:]
            )

            # 任务名称（简洁直观）
            job_name = f"VWF_{wt_aa}{pos}{mut_aa}"

            # 创建任务
            mut_job = AF3BatchJob(
                name=job_name,
                sequence=mut_sequence,
                count=1
            )
            self.batch_jobs.append(mut_job)

            if (idx + 1) % 100 == 0 or (idx + 1) == len(self.mutations_df):
                logger.info(f"      [{idx+1}/{len(self.mutations_df)}] {job_name}")

        logger.info(f"      => 共创建 {len(self.batch_jobs)} 个 Batch Jobs (含 WT)")

    def save_batch_jsons(self) -> List[Path]:
        """
        按 chunk_size 分组，生成 JSON 数组文件
        每个文件是一个 JSON 数组，包含多个任务对象
        """
        logger.info("[4/5] 生成 Batch JSON 文件...")

        batch_files = []
        total_jobs = len(self.batch_jobs)
        num_batches = (total_jobs + self.chunk_size - 1) // self.chunk_size

        logger.info(f"      分组: {total_jobs} 个任务, 每文件 {self.chunk_size} 个")
        logger.info(f"      预计生成 {num_batches} 个 JSON 文件")

        for batch_idx in range(num_batches):
            start = batch_idx * self.chunk_size
            end = min(start + self.chunk_size, total_jobs)
            batch_jobs = self.batch_jobs[start:end]

            # 构建 JSON 数组（官方格式：最外层是列表）
            batch_data = [job.to_dict() for job in batch_jobs]

            # 生成文件名: af3_batch_01.json, af3_batch_02.json...
            batch_filename = f"af3_batch_{batch_idx+1:02d}.json"
            batch_path = self.output_dir / batch_filename

            # 写入 JSON 数组
            with open(batch_path, 'w', encoding='utf-8') as f:
                json.dump(batch_data, f, indent=2, ensure_ascii=False)

            batch_files.append(batch_path)
            logger.info(f"      ✓ [{batch_idx+1}/{num_batches}] {batch_filename} ({len(batch_jobs)} 个任务)")

        logger.info(f"      => 共生成 {len(batch_files)} 个 Batch JSON 文件")
        return batch_files

    def save_manifest(self, batch_files: List[Path]):
        """保存任务清单文件，并标记已完成的任务"""
        total_size_mb = sum(f.stat().st_size for f in batch_files) / (1024 * 1024)

        # 计算已完成和待运行的任务
        completed_in_batch = []
        pending_in_batch = []
        for job in self.batch_jobs:
            if job.name in self.completed_jobs:
                completed_in_batch.append(job.name)
            else:
                pending_in_batch.append(job.name)

        # 构建每个 Batch 包含的任务列表（带状态）
        batch_details = []
        for i, batch_path in enumerate(batch_files):
            start = i * self.chunk_size
            end = min(start + self.chunk_size, len(self.batch_jobs))
            jobs_in_batch = self.batch_jobs[start:end]

            batch_details.append({
                "batch_number": i + 1,
                "filename": batch_path.name,
                "size_mb": round(batch_path.stat().st_size / (1024 * 1024), 2),
                "jobs_count": len(jobs_in_batch),
                "jobs": [
                    {
                        "name": job.name,
                        "status": "COMPLETED" if job.name in self.completed_jobs else "PENDING"
                    }
                    for job in jobs_in_batch
                ]
            })

        manifest = {
            "generated_at": datetime.now().isoformat(),
            "total_jobs": len(self.batch_jobs),
            "completed_jobs_count": len(completed_in_batch),
            "pending_jobs_count": len(pending_in_batch),
            "wt_sequence_length": len(self.wt_sequence),
            "chunk_size": self.chunk_size,
            "total_batch_files": len(batch_files),
            "total_size_mb": round(total_size_mb, 2),
            "batch_files": batch_details,
            "all_jobs": [
                {
                    "name": job.name,
                    "sequence_length": len(job.sequence),
                    "status": "COMPLETED" if job.name in self.completed_jobs else "PENDING"
                }
                for job in self.batch_jobs
            ],
            "completed_jobs_list": completed_in_batch,
            "pending_jobs_list": pending_in_batch
        }

        manifest_path = self.output_dir / "af3_jobs_manifest.json"
        with open(manifest_path, 'w', encoding='utf-8') as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False)

        logger.info(f"      ✓ 清单文件: {manifest_path.name}")

        # 如果有已完成的任务，生成剩余任务文件
        if self.completed_jobs:
            self._save_remaining_batch(pending_in_batch)

    def _save_remaining_batch(self, pending_job_names: List[str]):
        """生成只包含待运行任务的 batch 文件"""
        if not pending_job_names:
            logger.info("      ✅ 所有任务都已完成！")
            return

        # 创建待运行任务的 job 对象列表
        pending_jobs = [job for job in self.batch_jobs if job.name in pending_job_names]

        logger.info(f"      📋 发现 {len(pending_jobs)} 个待运行任务")

        # 按 chunk_size 切分待运行任务
        chunk_size = self.chunk_size
        num_batches = (len(pending_jobs) + chunk_size - 1) // chunk_size

        for batch_idx in range(num_batches):
            start = batch_idx * chunk_size
            end = min(start + chunk_size, len(pending_jobs))
            batch_jobs = pending_jobs[start:end]

            # 构建 JSON 数组
            batch_data = [job.to_dict() for job in batch_jobs]

            # 生成文件名: af3_batch_01_remaining.json, af3_batch_02_remaining.json...
            batch_filename = f"af3_batch_{batch_idx+1:02d}_REMAINING.json"
            batch_path = self.output_dir / batch_filename

            # 写入 JSON
            with open(batch_path, 'w', encoding='utf-8') as f:
                json.dump(batch_data, f, indent=2, ensure_ascii=False)

            completed_count = len([j for j in batch_jobs if j.name in self.completed_jobs])
            pending_count = len(batch_jobs) - completed_count

            logger.info(f"      ✓ [{batch_idx+1}/{num_batches}] {batch_filename} ({len(batch_jobs)} 个任务, {completed_count} 已完成, {pending_count} 待运行)")

        # 生成第一个剩余 batch 的摘要
        if pending_jobs:
            first_batch_jobs = pending_jobs[:chunk_size]
            logger.info(f"      📌 af3_batch_01_REMAINING.json 包含 {len(first_batch_jobs)} 个任务:")
            for job in first_batch_jobs[:10]:  # 只显示前 10 个
                status = "✓ 已完成" if job.name in self.completed_jobs else "⏳ 待运行"
                logger.info(f"         - {job.name} ({status})")
            if len(first_batch_jobs) > 10:
                logger.info(f"         ... 还有 {len(first_batch_jobs) - 10} 个任务")

    def run(self):
        """运行完整的 Phase 2 流程"""
        self.load_data()
        self.create_batch_jobs()
        batch_files = self.save_batch_jsons()
        self.save_manifest(batch_files)

        # 计算已完成和待运行的任务数
        completed_count = len([job for job in self.batch_jobs if job.name in self.completed_jobs])
        pending_count = len(self.batch_jobs) - completed_count

        # 打印摘要
        logger.info("="*60)
        logger.info("Phase 2 完成!")
        logger.info("="*60)
        logger.info(f"输出目录: {self.output_dir}")
        logger.info(f"总任务数: {len(self.batch_jobs)} (WT + {len(self.batch_jobs)-1} 个突变)")

        if self.completed_jobs:
            logger.info(f"  ✅ 已完成: {completed_count} 个")
            logger.info(f"  ⏳ 待运行: {pending_count} 个")

        logger.info(f"Batch JSON: {len(batch_files)} 个文件")
        for i, bf in enumerate(batch_files[:5]):  # 只显示前5个
            logger.info(f"  - {bf.name}")
        if len(batch_files) > 5:
            logger.info(f"  ... 还有 {len(batch_files) - 5} 个")

        # 如果有已完成任务，提示 REMAINING 文件
        if self.completed_jobs and pending_count > 0:
            logger.info("")
            logger.info("【🎯 增量上传模式】")
            logger.info(f"已生成 *_REMAINING.json 文件，只包含 {pending_count} 个待运行任务")
            logger.info(f"请上传: af3_batch_01_REMAINING.json (包含第一批待运行的 {min(self.chunk_size, pending_count)} 个任务)")

        logger.info("")
        logger.info("【上传策略建议】")
        logger.info(f"AF3 官网每日限额约 30 个任务")
        logger.info(f"当前每个 JSON 文件包含 {self.chunk_size} 个任务")
        logger.info(f"建议每天上传 1 个 JSON 文件")

        if self.completed_jobs:
            remaining_batches = (pending_count + self.chunk_size - 1) // self.chunk_size
            logger.info(f"预计还需要 {remaining_batches} 天完成剩余预测")
        else:
            logger.info(f"预计需要 {len(batch_files)} 天完成全部预测")

        logger.info("")
        logger.info("【下一步 - 人工操作】")
        logger.info("1. 访问 https://alphafoldserver.com/")
        logger.info("2. 登录后点击 'Batch Submission'")

        if self.completed_jobs and pending_count > 0:
            logger.info(f"3. 上传剩余任务: af3_batch_01_REMAINING.json → af3_batch_{remaining_batches:02d}_REMAINING.json")
        else:
            logger.info(f"3. 每天上传 1 个 JSON: af3_batch_01.json → af3_batch_{len(batch_files):02d}.json")

        logger.info("4. 等待服务器计算完成（通常 1-3 天/批次）")
        logger.info("5. 下载结果 CIF 文件，解压到 structures/predictions/")
        logger.info("6. 运行 Phase 3 进行结构评分")
        logger.info("="*60)


def main():
    parser = argparse.ArgumentParser(
        description='Phase 2: AlphaFold3 Batch JSON 生成器 (官方格式)'
    )
    parser.add_argument(
        '--csv',
        type=str,
        default='../output/phase1_filtered_mutations.csv',
        help='Phase 1 输出的 CSV 文件路径 (默认: ../output/phase1_filtered_mutations.csv)'
    )
    parser.add_argument(
        '--wt-fasta',
        type=str,
        default='../structures/wt/VWF_P04275_WT.fasta',
        help='WT FASTA 文件路径 (默认: ../structures/wt/VWF_P04275_WT.fasta)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='../output/af3_batches',
        help='输出目录 (默认: ../output/af3_batches)'
    )
    parser.add_argument(
        '--limit',
        type=int,
        default=None,
        help='限制只处理前 N 个变异（测试用）'
    )
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=30,
        help='每个 JSON 文件包含的任务数（默认: 30，匹配官方每日限额）'
    )
    parser.add_argument(
        '--completed-jobs',
        type=str,
        default=None,
        help='已完成的作业名称列表，逗号分隔（例如：VWF_WT,VWF_G1531D,VWF_T1538K）'
    )

    args = parser.parse_args()

    # 解析已完成的作业列表
    completed_jobs = None
    if args.completed_jobs:
        completed_jobs = [name.strip() for name in args.completed_jobs.split(',')]

    generator = AF3BatchGenerator(
        csv_path=args.csv,
        wt_fasta=args.wt_fasta,
        output_dir=args.output_dir,
        limit=args.limit,
        chunk_size=args.chunk_size,
        completed_jobs=completed_jobs
    )

    generator.run()


if __name__ == '__main__':
    main()

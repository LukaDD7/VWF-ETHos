#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
安全运行 07e God Mode 的监控脚本
- 预先检查文件状态
- 实时显示进度
- 自动备份保护
- 网络断点续传
"""

import subprocess
import sys
import time
import os
from pathlib import Path
from datetime import datetime

# 路径配置
BASE_DIR = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan")
RESULTS_DIR = BASE_DIR / "results"
LOGS_DIR = BASE_DIR / "logs"
SCRIPT_PATH = BASE_DIR / "scripts" / "07e_god_mode_epigenome_crawler.py"

# 关键文件
GOLD_PKL = RESULTS_DIR / "07_VCF_AlphaGenome_Results_Backup.pkl"
OUTPUT_CSV = RESULTS_DIR / "07e_GodMode_Epigenome_Peaks.csv"
CHECKPOINT_JSONL = RESULTS_DIR / "07e_godmode_checkpoint.jsonl"
LOG_FILE = LOGS_DIR / "07e_godmode_log.txt"

def print_header(msg):
    print(f"\n{'='*70}")
    print(f"  {msg}")
    print(f"{'='*70}")

def print_status(label, value, status="info"):
    icons = {"info": "ℹ️", "ok": "✅", "warn": "⚠️", "error": "❌"}
    icon = icons.get(status, "ℹ️")
    print(f"{icon} {label}: {value}")

def check_files():
    """预先检查文件状态"""
    print_header("步骤 1: 文件状态检查")

    # 检查金标准文件
    if GOLD_PKL.exists():
        size_gb = GOLD_PKL.stat().st_size / (1024**3)
        print_status("金标准文件", f"{GOLD_PKL.name} ({size_gb:.1f} GB)", "ok")
    else:
        print_status("金标准文件", "不存在! 需要 07_VCF_AlphaGenome_Results_Backup.pkl", "error")
        return False

    # 检查已完成的07e结果
    if OUTPUT_CSV.exists():
        # 计算行数
        with open(OUTPUT_CSV) as f:
            lines = sum(1 for _ in f)
        completed = lines - 1 if lines > 0 else 0  # 减去表头
        print_status("07e结果文件", f"{OUTPUT_CSV.name} ({completed} 条已完成)", "ok")
    else:
        print_status("07e结果文件", "不存在 (将新建)", "warn")

    # 检查checkpoint
    if CHECKPOINT_JSONL.exists():
        with open(CHECKPOINT_JSONL) as f:
            lines = sum(1 for _ in f)
        print_status("Checkpoint文件", f"{CHECKPOINT_JSONL.name} ({lines} 条记录)", "ok")
    else:
        print_status("Checkpoint文件", "不存在 (将新建)", "warn")

    # 备份保护检查
    backups = list(RESULTS_DIR.glob("backups/*07e*"))
    if backups:
        print_status("备份文件", f"找到 {len(backups)} 个备份", "ok")
    else:
        print_status("备份文件", "暂无", "info")

    return True

def count_remaining_tasks():
    """计算剩余任务数"""
    import pickle
    import pandas as pd

    # 加载总任务
    with open(GOLD_PKL, 'rb') as f:
        gold_data = pickle.load(f)

    if isinstance(gold_data, list):
        total = len(gold_data)
    elif isinstance(gold_data, pd.DataFrame):
        total = len(gold_data)
    else:
        total = "未知"

    # 加载已完成
    completed = 0
    if OUTPUT_CSV.exists():
        try:
            df = pd.read_csv(OUTPUT_CSV)
            completed = len(df)
        except:
            pass

    if isinstance(total, int):
        remaining = total - completed
        return total, completed, remaining
    return total, completed, "未知"

def run_07e_with_monitoring():
    """运行07e并监控"""
    print_header("步骤 2: 启动 07e God Mode")

    total, completed, remaining = count_remaining_tasks()
    print_status("总任务数", total)
    print_status("已完成", completed)
    print_status("待处理", remaining)

    print(f"\n⏳ 正在启动 07e_god_mode_epigenome_crawler.py...")
    print(f"📊 将自动断点续传，跳过已完成的 {completed} 个变异")
    print(f"💾 实时保存到: {OUTPUT_CSV}")
    print(f"📝 日志写入: {LOG_FILE}")
    print(f"\n⚠️  注意: 请确保VPN已连接，代理服务器 10.5.36.27:7897 可访问")
    print(f"⏹️  按 Ctrl+C 可随时中断，进度会自动保存\n")

    # 启动子进程
    env = os.environ.copy()
    # 确保Python路径包含alphagenome
    env["PYTHONPATH"] = str(Path("/media/luzhenyang/project/alphagenome/alphagenome"))

    process = subprocess.Popen(
        [sys.executable, str(SCRIPT_PATH)],
        cwd=str(BASE_DIR / "scripts"),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        bufsize=1,
        env=env
    )

    # 实时监控输出
    try:
        for line in process.stdout:
            print(line, end='')
            sys.stdout.flush()

            # 检测关键状态
            if "API 连接失败" in line:
                print("\n" + "="*70)
                print("❌ 检测到API连接失败!")
                print("   可能原因: VPN未连接 或 代理服务器不可达")
                print("   请检查:")
                print("   1. VPN是否已连接")
                print("   2. 能否访问 10.5.36.27:7897")
                print("="*70 + "\n")

            elif "成功:" in line and "ATAC=" in line:
                # 提取进度信息
                pass

    except KeyboardInterrupt:
        print("\n\n⚠️ 收到中断信号，正在停止...")
        process.terminate()
        time.sleep(2)
        process.kill()
        print("✅ 已安全停止。进度已保存，可稍后重新运行断点续传。")
        return

    process.wait()

    if process.returncode == 0:
        print_header("步骤 3: 运行完成")
        print_status("退出码", "0 (成功)", "ok")
    else:
        print_header("步骤 3: 运行异常")
        print_status("退出码", process.returncode, "error")

def main():
    print_header("07e God Mode 安全运行监控器")
    print(f"启动时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # 检查文件
    if not check_files():
        print("\n❌ 文件检查失败，请确保金标准文件存在后重试")
        sys.exit(1)

    # 确认运行
    print(f"\n{'='*70}")
    response = input("是否继续运行 07e? (yes/no): ").strip().lower()
    if response not in ['yes', 'y', '']:
        print("已取消")
        sys.exit(0)

    # 运行
    run_07e_with_monitoring()

if __name__ == "__main__":
    main()

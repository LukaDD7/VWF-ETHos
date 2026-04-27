#!/usr/bin/env python3
"""
Boltz-2 Progress Monitor
========================
在 GPU 服务器上独立运行，实时显示预测进度。
与 run_boltz2.sh 完全解耦，只读文件系统，不干扰预测进程。

用法：
  # 实时刷新（每 30 秒）
  python monitor_progress.py --watch

  # 单次快照
  python monitor_progress.py
"""

import argparse
import json
import os
import time
from pathlib import Path


def count_jobs(results_dir: Path):
    """扫描结果目录，统计 job 完成情况。"""
    done_jobs    = list(results_dir.glob("*/.done"))
    all_job_dirs = [d for d in results_dir.iterdir()
                    if d.is_dir() and not d.name.startswith("_")]
    total     = len(all_job_dirs)
    completed = len(done_jobs)
    failed    = []

    for d in all_job_dirs:
        pred_dir = d / "predictions"
        done_f   = d / ".done"
        # 有 predictions 目录但没有 .done → 可能是上次崩溃的脏目录
        if pred_dir.exists() and not done_f.exists():
            # 检查是否有 cif 文件（表示部分完成）
            cifs = list(pred_dir.glob("*.cif"))
            if cifs:
                failed.append(d.name)

    return total, completed, failed


def read_progress_json(results_dir: Path):
    pf = results_dir / "progress.json"
    if not pf.exists():
        return None
    try:
        return json.load(open(pf))
    except Exception:
        return None


def show_status(results_dir: Path):
    total, completed, failed = count_jobs(results_dir)
    pj = read_progress_json(results_dir)

    remaining = total - completed
    pct = completed / total * 100 if total > 0 else 0

    # 进度条
    bar_len = 40
    filled  = int(bar_len * completed / total) if total > 0 else 0
    bar     = "█" * filled + "░" * (bar_len - filled)

    print("\n" + "=" * 60)
    print("  VWF Boltz-2 Progress Monitor")
    print("=" * 60)
    print(f"  [{bar}]  {pct:.1f}%")
    print(f"  Completed : {completed:>5} / {total}")
    print(f"  Remaining : {remaining:>5}")
    print(f"  Dirty (crash remnants): {len(failed)}")
    if failed:
        print(f"  ↳ To clean: rm -rf {results_dir}/<name> && ./run_boltz2.sh")
        for f in failed[:5]:
            print(f"      - {f}")
        if len(failed) > 5:
            print(f"      ... and {len(failed)-5} more")

    if pj:
        print(f"\n  Current job    : {pj.get('current_job', 'N/A')}")
        print(f"  Last update    : {pj.get('last_update', 'N/A')}")

    # 速度估算（基于最近完成的 .done 文件时间戳）
    done_files = sorted(
        Path(results_dir).glob("*/.done"),
        key=lambda f: f.stat().st_mtime,
        reverse=True
    )
    if len(done_files) >= 2:
        newest = done_files[0].stat().st_mtime
        oldest = done_files[-1].stat().st_mtime
        elapsed_h = (newest - oldest) / 3600
        if elapsed_h > 0:
            rate = len(done_files) / elapsed_h
            eta_h = remaining / rate if rate > 0 else float("inf")
            print(f"\n  Speed          : {rate:.1f} jobs/hour")
            print(f"  ETA            : {eta_h:.1f} hours")

    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Monitor Boltz-2 prediction progress")
    parser.add_argument(
        "--results-dir", type=str,
        default=str(Path(__file__).parent / "../../output/boltz2_results"),
        help="Path to boltz2_results directory"
    )
    parser.add_argument(
        "--watch", action="store_true",
        help="Continuously refresh every 30 seconds"
    )
    parser.add_argument(
        "--interval", type=int, default=30,
        help="Refresh interval in seconds (default: 30)"
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir).resolve()

    if not results_dir.exists():
        print(f"[ERROR] Results dir not found: {results_dir}")
        print("  Run run_boltz2.sh first.")
        return

    if args.watch:
        print(f"Watching {results_dir} (Ctrl+C to stop, refresh every {args.interval}s)")
        try:
            while True:
                # 清屏
                os.system("clear" if os.name == "posix" else "cls")
                show_status(results_dir)
                print(f"\n  [Auto-refresh every {args.interval}s — Ctrl+C to exit]")
                time.sleep(args.interval)
        except KeyboardInterrupt:
            print("\nMonitor stopped.")
    else:
        show_status(results_dir)


if __name__ == "__main__":
    main()

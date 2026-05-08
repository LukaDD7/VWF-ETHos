#!/usr/bin/env python3
"""
upload_boltz2_results_huggingface.py
======================================
把 Boltz-2 A1+GPIbα 预测结果上传到 HuggingFace Dataset。

用法：
    # 先登录（只需一次）
    pip install huggingface_hub
    huggingface-cli login

    # 上传
    python scripts/pipeline/upload_boltz2_results_huggingface.py

上传内容：
    - 74 个 job 的 .cif 结构文件（5 samples/job = 370 个 .cif）
    - 74 个 job 的 confidence_*.json（370 个文件）
    - iptm_results.csv 汇总

HuggingFace 公开数据集限额：无限（免费）
"""

import argparse
import json
import shutil
from pathlib import Path

try:
    from huggingface_hub import HfApi
    HAS_HUGGINGFACE = True
except ImportError:
    HAS_HUGGINGFACE = False
    print("[ERROR] huggingface_hub not installed.")
    print("  Run: pip install huggingface_hub")
    import sys; sys.exit(1)


REPO_ID = "LukaDD7/VWF-A1-GPIb-alpha-structures"
RESULTS_DIR = Path(__file__).parent.parent.parent / "output" / "boltz2_a1_gpiba_results"
ANALYSIS_DIR = Path(__file__).parent.parent.parent / "output" / "boltz2_a1_gpiba_analysis"


def upload_results():
    api = HfApi()

    # 确认仓库存在
    try:
        api.repo_info(repo_id=REPO_ID, repo_type="dataset")
        print(f"[OK] Repo exists: {REPO_ID}")
    except Exception:
        print(f"[INFO] Repo not found, creating: {REPO_ID}")
        api.create_repo(repo_id=REPO_ID, repo_type="dataset", exist_ok=True)

    # 扫描要上传的文件
    cif_files = list(RESULTS_DIR.glob("*/predictions/*/*.cif"))
    conf_files = list(RESULTS_DIR.glob("*/predictions/*/confidence_*.json"))
    csv_files = list(ANALYSIS_DIR.glob("*.csv"))

    print(f"Found {len(cif_files)} .cif files")
    print(f"Found {len(conf_files)} confidence .json files")
    print(f"Found {len(csv_files)} CSV files")

    total = len(cif_files) + len(conf_files) + len(csv_files)
    print(f"\nTotal files to upload: {total}")

    # 上传文件夹（递归）
    folder_path = RESULTS_DIR
    print(f"\nUploading {folder_path} ...")
    api.upload_folder(
        folder_path=str(folder_path),
        repo_id=REPO_ID,
        repo_type="dataset",
        path_in_repo="boltz2_a1_gpiba_results",
    )
    print("[OK] boltz2_a1_gpiba_results/ uploaded")

    # 上传分析结果
    for csv in csv_files:
        print(f"Uploading {csv.name} ...")
        api.upload_file(
            path_or_fileobj=str(csv),
            path_in_repo=f"analysis/{csv.name}",
            repo_id=REPO_ID,
            repo_type="dataset",
        )
        print(f"[OK] analysis/{csv.name}")

    print(f"\n[OK] All done! Browse at: https://huggingface.co/datasets/{REPO_ID}")


def download_results():
    """在另一台设备上下载结果。"""
    print("Download instructions:")
    print(f"  from datasets import load_dataset")
    print(f"  ds = load_dataset('{REPO_ID}', split='train')")
    print()
    print(f"  # Or just the CSV:")
    print(f"  import pandas as pd")
    print(f"  df = pd.read_csv('https://huggingface.co/datasets/{REPO_ID}/raw/main/analysis/iptm_results.csv')")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--download", action="store_true", help="Show download instructions")
    args = parser.parse_args()

    if args.download:
        download_results()
    else:
        upload_results()
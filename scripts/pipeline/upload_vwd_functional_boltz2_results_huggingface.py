#!/usr/bin/env python3
"""
upload_vwd_functional_boltz2_results_huggingface.py
===================================================
上传 VWD/VWF functional Boltz-2 panel 结果到 HuggingFace Dataset。

修复 v2：使用 upload_folder 分批上传，每次批量算一次 commit。
速率限制：128 commits/hour，单文件上传太慢。

用法：
    python scripts/pipeline/upload_vwd_functional_boltz2_results_huggingface.py
"""

import os
import time
import shutil
from pathlib import Path

try:
    from huggingface_hub import HfApi, login
except ImportError:
    print("[ERROR] huggingface_hub not installed.")
    import sys; sys.exit(1)

REPO_ID = "lucachangretta/VWF"
VWD_RESULTS_DIR = Path(__file__).parent.parent.parent / "output" / "boltz2_vwd_functional_panel" / "boltz_results"
VWD_ANALYSIS_DIR = Path(__file__).parent.parent.parent / "output" / "boltz2_vwd_functional_panel"
HF_SUBFOLDER = "vwd_functional_panel"
BATCH_SIZE = 200  # 每个子目录多少个文件


def get_token():
    token = os.environ.get("HF_TOKEN")
    if not token:
        raise ValueError("HF_TOKEN environment variable not set")
    return token


def upload_analysis_files(api):
    """上传分析文件。"""
    print("\n=== Step 1: Upload analysis files ===")
    csv_files = [
        VWD_ANALYSIS_DIR / "diagnostic_panel.csv",
        VWD_ANALYSIS_DIR / "job_manifest.csv",
        VWD_ANALYSIS_DIR / "summary.json",
    ]
    for f in csv_files:
        if f.exists():
            hf_path = f"vwd_analysis/{f.name}"
            print(f"Uploading {f.name} ...")
            for attempt in range(3):
                try:
                    api.upload_file(
                        path_or_fileobj=str(f),
                        path_in_repo=hf_path,
                        repo_id=REPO_ID,
                        repo_type="dataset",
                    )
                    print(f"  [OK] vwd_analysis/{f.name}")
                    break
                except Exception as e:
                    if attempt < 2:
                        print(f"  [WARN] Retry {attempt+1}: {e}")
                        time.sleep(5)
                    else:
                        print(f"  [ERROR] Failed: {e}")
    time.sleep(5)  # 分析文件只需3次commit，休息一下


def batch_upload_folder(api, local_dir, hf_path_prefix, batch_name, max_per_batch=200):
    """
    把大量文件分批组织到临时目录，然后用 upload_folder 上传。
    每次上传只产生一个 commit。
    """
    files = sorted(local_dir.glob("*"))
    if not files:
        return

    total = len(files)
    n_batches = (total + max_per_batch - 1) // max_per_batch

    for i in range(0, total, max_per_batch):
        batch_files = files[i:i+max_per_batch]
        batch_num = i // max_per_batch + 1

        # 创建临时批量目录
        batch_dir = Path("/tmp/hf_batch_upload")
        batch_subdir = batch_dir / f"{batch_name}_batch{batch_num}"
        if batch_subdir.exists():
            shutil.rmtree(batch_subdir)
        batch_subdir.mkdir(parents=True)

        # 复制文件到临时目录（保持文件名）
        for f in batch_files:
            shutil.copy2(f, batch_subdir / f.name)

        print(f"\n  [{batch_name}] Batch {batch_num}/{n_batches}: {len(batch_files)} files")
        hf_path = f"{hf_path_prefix}/{batch_name}_batch{batch_num}"

        for attempt in range(3):
            try:
                api.upload_folder(
                    folder_path=str(batch_subdir),
                    repo_id=REPO_ID,
                    repo_type="dataset",
                    path_in_repo=hf_path,
                )
                print(f"    [OK] {hf_path}")
                break
            except Exception as e:
                if attempt < 2:
                    print(f"    [WARN] Retry {attempt+1}: {e}")
                    time.sleep(10)
                else:
                    print(f"    [ERROR] Failed: {e}")

        # 清理临时目录
        shutil.rmtree(batch_subdir)

        # 每批次之间休息一下，减少 rate limit 压力
        if i + max_per_batch < total:
            time.sleep(3)


def find_job_dirs():
    """查找所有 job 子目录。"""
    job_dirs = []
    for d in VWD_RESULTS_DIR.iterdir():
        if d.is_dir():
            pred_dir = d / "predictions"
            if pred_dir.exists():
                job_dirs.append(d)
    return sorted(job_dirs)


def upload_results():
    token = get_token()
    login(token)
    api = HfApi()

    # 确认仓库
    try:
        api.repo_info(repo_id=REPO_ID, repo_type="dataset")
        print(f"[OK] Repo exists: {REPO_ID}")
    except Exception:
        print(f"[INFO] Creating repo: {REPO_ID}")
        api.create_repo(repo_id=REPO_ID, repo_type="dataset", exist_ok=True)

    # Step 1: 分析文件
    upload_analysis_files(api)

    # Step 2: 收集所有 job 的 predictions 目录
    job_dirs = find_job_dirs()
    print(f"\n=== Step 2: Found {len(job_dirs)} job directories")

    # Step 3: 对每个 job 目录，分批上传 CIF 和 confidence 文件
    print("\n=== Step 3: Upload CIF + CONF files by job directory")

    for j, job_dir in enumerate(job_dirs):
        job_name = job_dir.name
        print(f"\n[{j+1}/{len(job_dirs)}] {job_name}")

        pred_dir = job_dir / "predictions"
        if not pred_dir.exists():
            print(f"  [SKIP] No predictions/")
            continue

        # 找 CIF 和 confidence 文件
        cif_files = sorted(pred_dir.glob("**/*.cif"))
        conf_files = sorted(pred_dir.glob("**/confidence_*.json"))

        print(f"  CIF: {len(cif_files)}, CONF: {len(conf_files)}")

        # 临时目录存放要上传的文件
        temp_upload = Path("/tmp/hf_upload_temp")
        if temp_upload.exists():
            shutil.rmtree(temp_upload)
        temp_upload.mkdir()

        # 复制 CIF 到 temp
        cif_temp = temp_upload / "cif"
        cif_temp.mkdir()
        for f in cif_files:
            shutil.copy2(f, cif_temp / f.name)

        # 复制 CONF 到 temp
        conf_temp = temp_upload / "conf"
        conf_temp.mkdir()
        for f in conf_files:
            shutil.copy2(f, conf_temp / f.name)

        # 上传 CIF
        if cif_files:
            print(f"  Uploading CIF...")
            for attempt in range(3):
                try:
                    api.upload_folder(
                        folder_path=str(cif_temp),
                        repo_id=REPO_ID,
                        repo_type="dataset",
                        path_in_repo=f"{HF_SUBFOLDER}/{job_name}/cif",
                    )
                    print(f"    [OK] CIF uploaded")
                    break
                except Exception as e:
                    if attempt < 2:
                        print(f"    [WARN] CIF retry {attempt+1}: {e}")
                        time.sleep(10)
                    else:
                        print(f"    [ERROR] CIF failed: {e}")

        # 上传 CONF
        if conf_files:
            print(f"  Uploading CONF...")
            for attempt in range(3):
                try:
                    api.upload_folder(
                        folder_path=str(conf_temp),
                        repo_id=REPO_ID,
                        repo_type="dataset",
                        path_in_repo=f"{HF_SUBFOLDER}/{job_name}/confidence",
                    )
                    print(f"    [OK] CONF uploaded")
                    break
                except Exception as e:
                    if attempt < 2:
                        print(f"    [WARN] CONF retry {attempt+1}: {e}")
                        time.sleep(10)
                    else:
                        print(f"    [ERROR] CONF failed: {e}")

        # 清理临时目录
        shutil.rmtree(temp_upload)

        # 每个 job 之休息一下
        time.sleep(2)

    print(f"\n[OK] All done!")
    print(f"  Browse: https://huggingface.co/datasets/{REPO_ID}/tree/main/{HF_SUBFOLDER}")


if __name__ == "__main__":
    upload_results()
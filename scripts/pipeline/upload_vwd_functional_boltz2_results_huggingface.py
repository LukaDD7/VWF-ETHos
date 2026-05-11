#!/usr/bin/env python3
"""Upload VWD/VWF functional Boltz-2 panel results to HuggingFace.

The first implementation uploaded each job directory separately. That creates
roughly 2 commits per job and hits HuggingFace's commit-rate limit. This version
uses a small number of commits:

  * analysis files are uploaded as one folder commit
  * raw results are either uploaded as one browsable folder commit, or as a
    small set of tar.gz shards in one archive-folder commit

Recommended robust mode for the 2.9GB / ~15k-file functional panel:

    HF_TOKEN=... python scripts/pipeline/upload_vwd_functional_boltz2_results_huggingface.py --mode archive

Browsable tree mode, if the Hub accepts a very large single folder commit:

    HF_TOKEN=... python scripts/pipeline/upload_vwd_functional_boltz2_results_huggingface.py --mode folder
"""

import argparse
import os
import shutil
import tarfile
import tempfile
from pathlib import Path

try:
    from huggingface_hub import HfApi
except ImportError:
    print("[ERROR] huggingface_hub not installed.")
    import sys; sys.exit(1)

REPO_ID = "lucachangretta/VWF"
VWD_RESULTS_DIR = Path(__file__).parent.parent.parent / "output" / "boltz2_vwd_functional_panel" / "boltz_results"
VWD_ANALYSIS_DIR = Path(__file__).parent.parent.parent / "output" / "boltz2_vwd_functional_panel"
HF_SUBFOLDER = "vwd_functional_panel"
HF_ANALYSIS_SUBFOLDER = "vwd_analysis"
HF_ARCHIVE_SUBFOLDER = "vwd_functional_panel_archives"
DEFAULT_JOBS_PER_ARCHIVE = 75


def get_token():
    token = os.environ.get("HF_TOKEN")
    if not token:
        raise ValueError("HF_TOKEN environment variable not set")
    return token


def ensure_repo(api, repo_id):
    try:
        api.repo_info(repo_id=repo_id, repo_type="dataset")
        print(f"[OK] Repo exists: {repo_id}")
    except Exception:
        print(f"[INFO] Creating repo: {repo_id}")
        api.create_repo(repo_id=repo_id, repo_type="dataset", exist_ok=True)


def upload_analysis_files(api, repo_id):
    """Upload analysis files in one commit."""
    print("\n=== Step 1: Upload analysis files ===")
    analysis_files = [
        VWD_ANALYSIS_DIR / "diagnostic_panel.csv",
        VWD_ANALYSIS_DIR / "job_manifest.csv",
        VWD_ANALYSIS_DIR / "summary.json",
        VWD_ANALYSIS_DIR / "boltz_results_summary.csv",
    ]
    existing_files = [f for f in analysis_files if f.exists()]
    if not existing_files:
        print("  [SKIP] No analysis files found.")
        return

    with tempfile.TemporaryDirectory(prefix="vwd_hf_analysis_") as tmp:
        stage = Path(tmp)
        for f in existing_files:
            shutil.copy2(f, stage / f.name)

        print(f"  Uploading {len(existing_files)} files as one commit -> {HF_ANALYSIS_SUBFOLDER}/")
        api.upload_folder(
            folder_path=str(stage),
            repo_id=repo_id,
            repo_type="dataset",
            path_in_repo=HF_ANALYSIS_SUBFOLDER,
            commit_message="Upload VWD functional panel analysis files",
        )
        print("  [OK] analysis upload complete")


def find_job_dirs():
    """查找所有 job 子目录。"""
    job_dirs = []
    if not VWD_RESULTS_DIR.exists():
        return job_dirs
    for d in VWD_RESULTS_DIR.iterdir():
        if d.is_dir():
            pred_dir = d / "predictions"
            if pred_dir.exists():
                job_dirs.append(d)
    return sorted(job_dirs)


def upload_results_folder(api, repo_id):
    """Upload the raw results tree in one folder commit."""
    if not VWD_RESULTS_DIR.exists():
        raise FileNotFoundError(f"Results directory not found: {VWD_RESULTS_DIR}")

    print("\n=== Step 2: Upload raw results as one folder commit ===")
    print(f"  Local : {VWD_RESULTS_DIR}")
    print(f"  Remote: {HF_SUBFOLDER}/boltz_results")
    api.upload_folder(
        folder_path=str(VWD_RESULTS_DIR),
        repo_id=repo_id,
        repo_type="dataset",
        path_in_repo=f"{HF_SUBFOLDER}/boltz_results",
        commit_message="Upload VWD functional panel Boltz-2 raw results",
        ignore_patterns=["**/lightning_logs/**", "**/.DS_Store"],
    )
    print("  [OK] raw folder upload complete")


def create_result_archives(jobs_per_archive, archive_dir):
    """Create tar.gz shards by job directory to avoid thousands of Hub file ops."""
    job_dirs = find_job_dirs()
    if not job_dirs:
        raise FileNotFoundError(f"No job directories with predictions found in {VWD_RESULTS_DIR}")

    archive_dir.mkdir(parents=True, exist_ok=True)
    n_archives = (len(job_dirs) + jobs_per_archive - 1) // jobs_per_archive
    archive_paths = []

    print("\n=== Step 2: Create tar.gz result shards ===")
    print(f"  Jobs          : {len(job_dirs)}")
    print(f"  Jobs/archive  : {jobs_per_archive}")
    print(f"  Archives      : {n_archives}")
    print(f"  Archive dir   : {archive_dir}")

    for batch_idx, start in enumerate(range(0, len(job_dirs), jobs_per_archive), start=1):
        batch = job_dirs[start:start + jobs_per_archive]
        archive_path = archive_dir / f"vwd_functional_panel_boltz_results_{batch_idx:03d}_of_{n_archives:03d}.tar.gz"
        if archive_path.exists():
            print(f"  [SKIP] exists: {archive_path.name}")
            archive_paths.append(archive_path)
            continue

        print(f"  Creating {archive_path.name}: {len(batch)} jobs")
        with tarfile.open(archive_path, "w:gz") as tar:
            for job_dir in batch:
                tar.add(job_dir, arcname=f"boltz_results/{job_dir.name}")
        archive_paths.append(archive_path)

    return archive_paths


def upload_results_archives(api, repo_id, jobs_per_archive, archive_dir):
    """Upload archive shards in one commit."""
    archive_paths = create_result_archives(jobs_per_archive, archive_dir)

    with tempfile.TemporaryDirectory(prefix="vwd_hf_archives_") as tmp:
        stage = Path(tmp)
        for archive_path in archive_paths:
            # Hardlink first to avoid another multi-GB copy; fall back to copy.
            dest = stage / archive_path.name
            try:
                os.link(archive_path, dest)
            except OSError:
                shutil.copy2(archive_path, dest)

        print("\n=== Step 3: Upload result archives as one commit ===")
        print(f"  Files : {len(archive_paths)}")
        print(f"  Remote: {HF_ARCHIVE_SUBFOLDER}/")
        api.upload_folder(
            folder_path=str(stage),
            repo_id=repo_id,
            repo_type="dataset",
            path_in_repo=HF_ARCHIVE_SUBFOLDER,
            commit_message="Upload VWD functional panel Boltz-2 result archives",
        )
        print("  [OK] archive upload complete")


def upload_results(mode, repo_id, jobs_per_archive, archive_dir, analysis_only, results_only):
    token = get_token()
    api = HfApi(token=token)

    ensure_repo(api, repo_id)

    if not results_only:
        upload_analysis_files(api, repo_id)

    if not analysis_only:
        if mode == "folder":
            upload_results_folder(api, repo_id)
        elif mode == "archive":
            upload_results_archives(api, repo_id, jobs_per_archive, archive_dir)
        else:
            raise ValueError(f"Unsupported mode: {mode}")

    print(f"\n[OK] All done!")
    print(f"  Browse: https://huggingface.co/datasets/{repo_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-id", default=REPO_ID)
    parser.add_argument(
        "--mode",
        choices=["archive", "folder"],
        default="archive",
        help="archive is robust and commit-light; folder keeps individual files browsable.",
    )
    parser.add_argument("--jobs-per-archive", type=int, default=DEFAULT_JOBS_PER_ARCHIVE)
    parser.add_argument(
        "--archive-dir",
        type=Path,
        default=VWD_ANALYSIS_DIR / "hf_archives",
        help="Where to create/reuse tar.gz result shards.",
    )
    parser.add_argument("--analysis-only", action="store_true")
    parser.add_argument("--results-only", action="store_true")
    args = parser.parse_args()

    if args.analysis_only and args.results_only:
        raise SystemExit("--analysis-only and --results-only are mutually exclusive")

    upload_results(
        mode=args.mode,
        repo_id=args.repo_id,
        jobs_per_archive=args.jobs_per_archive,
        archive_dir=args.archive_dir,
        analysis_only=args.analysis_only,
        results_only=args.results_only,
    )

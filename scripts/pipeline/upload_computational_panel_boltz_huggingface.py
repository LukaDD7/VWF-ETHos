#!/usr/bin/env python3
"""Upload the 16-case computational panel Boltz raw results to HuggingFace.

Source: outputs/computational_panel_20260829/boltz/raw/  (32 jobs, ~231 MB)
Destination: lucachangretta/VWF -> computational_panel_20260829/boltz_raw/

Each job directory (boltz_results_VWF_<variant>__<assay>) contains the full
Boltz-2 output tree: predictions/ (5 model CIFs + confidence JSONs +
pae/pde/plddt NPZs), processed/, lightning_logs/, and a .done marker.

One tar.gz per job preserves the on-disk tree and keeps the repo clean of
thousands of small files; a manifest.json records per-job file counts,
sizes, and sha256 of every archive for integrity checks.

Usage:
    HF_TOKEN=... python scripts/pipeline/upload_computational_panel_boltz_huggingface.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tarfile
from pathlib import Path

try:
    from huggingface_hub import HfApi
except ImportError:
    raise SystemExit("[ERROR] huggingface_hub not installed")

REPO_ID = "lucachangretta/VWF"
REMOTE_SUBFOLDER = "computational_panel_20260829/boltz_raw"

ROOT = Path(__file__).parent.parent.parent
RAW_DIR = ROOT / "outputs" / "computational_panel_20260829" / "boltz" / "raw"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_job(job_dir: Path) -> dict:
    pred_dir = job_dir / "predictions" / job_dir.name.replace("boltz_results_", "")
    return {
        "cif_models": len(list(pred_dir.glob("*model_*.cif"))),
        "confidence_jsons": len(list(pred_dir.glob("confidence_*model_*.json"))),
        "pae_npz": len(list(pred_dir.glob("pae_*model_*.npz"))),
        "pde_npz": len(list(pred_dir.glob("pde_*model_*.npz"))),
        "plddt_npz": len(list(pred_dir.glob("plddt_*model_*.npz"))),
        "done_marker": (job_dir / ".done").exists(),
    }


def build_archives(stage_dir: Path) -> tuple[list[Path], dict]:
    stage_dir.mkdir(parents=True, exist_ok=True)
    jobs = sorted([d for d in RAW_DIR.iterdir() if d.is_dir()])
    archive_paths: list[Path] = []
    manifest: dict[str, dict] = {}

    for job in jobs:
        counts = verify_job(job)
        ok = (
            counts["cif_models"] == 5
            and counts["confidence_jsons"] == 5
            and counts["pae_npz"] == 5
            and counts["pde_npz"] == 5
            and counts["plddt_npz"] == 5
            and counts["done_marker"]
        )
        if not ok:
            print(f"  [SKIP] incomplete job {job.name}: {counts}")
            continue

        arc = stage_dir / f"{job.name}.tar.gz"
        if not arc.exists():
            print(f"  tar {arc.name}")
            subprocess.run(
                f"tar -czf {arc} -C {RAW_DIR} {job.name}",
                shell=True,
                check=True,
            )
        manifest[job.name] = {
            **counts,
            "archive": arc.name,
            "size_bytes": arc.stat().st_size,
            "sha256": sha256(arc),
        }
        archive_paths.append(arc)

    return archive_paths, manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", default="/tmp/vwf_boltz_raw_stage")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    token = os.environ.get("HF_TOKEN")
    if not token and not args.dry_run:
        raise SystemExit("HF_TOKEN environment variable not set")

    stage = Path(args.stage)
    archives, manifest = build_archives(stage)
    print(f"\n{len(archives)} verified job archives, "
          f"total {sum(a.stat().st_size for a in archives) / 1e6:.1f} MB")

    (stage / "manifest.json").write_text(
        json.dumps(
            {
                "description": "16-case computational panel Boltz-2 raw results "
                               "(experimental-context assays, 5 models per job)",
                "source": "outputs/computational_panel_20260829/boltz/raw/",
                "job_count": len(archives),
                "jobs": manifest,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    if args.dry_run:
        print("[DRY-RUN] archives built, manifest written, upload skipped")
        return 0

    api = HfApi(token=token)
    try:
        api.repo_info(repo_id=REPO_ID, repo_type="dataset")
        print(f"[OK] Repo exists: {REPO_ID}")
    except Exception:
        api.create_repo(repo_id=REPO_ID, repo_type="dataset", exist_ok=True)

    print(f"=== Upload {len(archives)} archives + manifest -> {REMOTE_SUBFOLDER}/ ===")
    api.upload_folder(
        folder_path=str(stage),
        repo_id=REPO_ID,
        repo_type="dataset",
        path_in_repo=REMOTE_SUBFOLDER,
        commit_message="Upload 16-case computational panel Boltz raw results (32 jobs)",
    )
    print(f"[OK] https://huggingface.co/datasets/{REPO_ID}/tree/main/{REMOTE_SUBFOLDER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

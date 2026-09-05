#!/usr/bin/env python3
"""Upload the 7A6O AIM-A1 50 ns MD trajectories to HuggingFace (optional layer).

Source: outputs/computational_panel_20260829/md/trajectories/ (7 XTC, ~62 MB
each: PANEL_WT_MATCHED + A1461D, P1413L, R1308C, R1341W, S1310F, V1316M)
Destination: lucachangretta/VWF -> computational_panel_20260829/md/trajectories/

Uploaded as-is (XTC is already a compressed binary format; re-tarring saves
<1%). A manifest.json records per-trajectory size + sha256. Backup files
matching '#...xtc.1#' are excluded.

Usage:
    HF_TOKEN=... python scripts/pipeline/upload_7a6o_trajectories_huggingface.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

try:
    from huggingface_hub import HfApi
except ImportError:
    raise SystemExit("[ERROR] huggingface_hub not installed")

REPO_ID = "lucachangretta/VWF"
REMOTE_SUBFOLDER = "computational_panel_20260829/md/trajectories"

ROOT = Path(__file__).parent.parent.parent
SRC_DIR = ROOT / "outputs" / "computational_panel_20260829" / "md" / "trajectories"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    token = os.environ.get("HF_TOKEN")
    if not token and not args.dry_run:
        raise SystemExit("HF_TOKEN environment variable not set")

    xtcs = sorted(SRC_DIR.glob("*_prod_concat.xtc"))
    if not xtcs:
        raise SystemExit(f"[ERROR] no XTC files found in {SRC_DIR}")

    manifest = {
        p.name: {"size_bytes": p.stat().st_size, "sha256": sha256(p)}
        for p in xtcs
    }
    print(f"{len(xtcs)} XTC trajectories, "
          f"total {sum(p.stat().st_size for p in xtcs) / 1e6:.1f} MB")

    manifest_path = SRC_DIR / "upload_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "description": "7A6O AIM-A1 experimental-structure 50 ns MD "
                               "trajectories (P0 panel, GROMACS 2025.4 CHARMM27)",
                "source": "outputs/computational_panel_20260829/md/trajectories/",
                "trajectories": manifest,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    if args.dry_run:
        print("[DRY-RUN] manifest written, upload skipped")
        return 0

    api = HfApi(token=token)
    try:
        api.repo_info(repo_id=REPO_ID, repo_type="dataset")
        print(f"[OK] Repo exists: {REPO_ID}")
    except Exception:
        api.create_repo(repo_id=REPO_ID, repo_type="dataset", exist_ok=True)

    print(f"=== Upload {len(xtcs)} XTC + manifest -> {REMOTE_SUBFOLDER}/ ===")
    for xtc in xtcs:
        print(f"  {xtc.name}")
        api.upload_file(
            path_or_fileobj=str(xtc),
            path_in_repo=f"{REMOTE_SUBFOLDER}/{xtc.name}",
            repo_id=REPO_ID,
            repo_type="dataset",
            commit_message=f"Add {xtc.name} (7A6O 50ns MD trajectory)",
        )
    api.upload_file(
        path_or_fileobj=str(manifest_path),
        path_in_repo=f"{REMOTE_SUBFOLDER}/upload_manifest.json",
        repo_id=REPO_ID,
        repo_type="dataset",
        commit_message="Add 7A6O trajectory upload manifest (sha256)",
    )
    print(f"[OK] https://huggingface.co/datasets/{REPO_ID}/tree/main/{REMOTE_SUBFOLDER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

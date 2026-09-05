#!/usr/bin/env python3
"""Package and upload VWF P0 experimental-structure MD results to HuggingFace.

Three-axis P0 MD (20 ns production, experimental-structure backbone):
    a1_gpiba   (1SQ0) : WT, R1308C, S1310F, A1461D, R1341W
    a2_folded  (3GXB) : WT, A1500V
    dprime_d3  (6N29) : WT, R1205H (md_20ns_p0_retry_safe)

One tar.gz per (axis, variant), preserving the on-disk tree:
    <archive>  ->  gromacs_md_<axis>/<variant>/...

Uploads archives + a manifest to lucachangretta/VWF -> md_p0_experimental_archives/.

Usage:
    HF_TOKEN=... python scripts/pipeline/upload_p0_md_huggingface.py
"""

import argparse
import json
import os
import shutil
import subprocess
import tarfile
from pathlib import Path

try:
    from huggingface_hub import HfApi
except ImportError:
    print("[ERROR] huggingface_hub not installed.")
    import sys; sys.exit(1)

REPO_ID = "lucachangretta/VWF"
REMOTE_SUBFOLDER = "md_p0_experimental_archives"

ROOT = Path(__file__).parent.parent.parent
OUT = ROOT / "output"

# axis -> {variant -> production md_tag dir under output/gromacs_md_<axis>/<variant>/}
RUNS = {
    "a1_gpiba": {
        "WT": "md_20ns_p0",
        "R1308C": "md_20ns_p0",
        "S1310F": "md_20ns_p0",
        "A1461D": "md_20ns_p0",
        "R1341W": "md_20ns_p0",
    },
    "a2_folded": {
        "WT": "md_20ns_p0",
        "A1500V": "md_20ns_p0",
    },
    "dprime_d3": {
        "WT": "md_20ns_p0",
        "R1205H": "md_20ns_p0_retry_safe",
    },
}

# Merge runs we already joined via trjcat (full 20 ns trajectory present).
MERGED_FULL = {"dprime_d3"}


def get_token():
    token = os.environ.get("HF_TOKEN")
    if not token:
        raise ValueError("HF_TOKEN environment variable not set")
    return token


def make_archives(stage_dir: Path) -> tuple[list[Path], dict]:
    """Create one tar.gz per (axis, variant). Return (archive paths, manifest)."""
    stage_dir.mkdir(parents=True, exist_ok=True)
    archive_paths = []
    manifest = {}

    for axis, variants in RUNS.items():
        axis_root = OUT / f"gromacs_md_{axis}"
        manifest[axis] = {}
        for variant, md_tag in variants.items():
            variant_dir = axis_root / variant
            md_dir = variant_dir / md_tag
            if not md_dir.exists():
                print(f"  [WARN] missing {md_dir}; skipping {axis}/{variant}")
                continue

            name = f"{axis}_{variant}.tar.gz"
            arc = stage_dir / name
            if arc.exists():
                print(f"  [SKIP] exists: {name}")
            else:
                # Exclude the archived failed CUDA-700 run dir (R1205H only).
                exclude_cuda700 = ""
                if (variant_dir / "md_20ns_p0_cuda700_v1").exists():
                    exclude_cuda700 = "--exclude=md_20ns_p0_cuda700_v1"
                rel = variant_dir.relative_to(OUT)
                print(f"  tar {name}  <- {rel}")
                comp = "-I pigz" if shutil.which("pigz") else "-z"
                cmd = (
                    f"tar {exclude_cuda700} -c {comp} -f {arc} -C {OUT} "
                    f"{rel}"
                )
                subprocess.run(cmd, shell=True, check=True)

            size = arc.stat().st_size if arc.exists() else 0

            # Inventory the production run.
            full_xtc = None
            if axis in MERGED_FULL and (md_dir / "md_prod_full_20ns.xtc").exists():
                full_xtc = "md_prod_full_20ns.xtc"
            elif (md_dir / "md_prod.xtc").exists() and not list(md_dir.glob("md_prod.part*.xtc")):
                full_xtc = "md_prod.xtc"
            frames = "-"
            if full_xtc:
                # frame count via gmx check is expensive; trust trjcat/WT merge.
                pass

            manifest[axis][variant] = {
                "md_tag": md_tag,
                "archive": name,
                "size_bytes": size,
                "merged_full_xtc": full_xtc is not None,
                "has_final_gro": (md_dir / "md_prod.gro").exists()
                or bool(list(md_dir.glob("md_prod.part*.gro"))),
            }

            archive_paths.append(arc)

    return archive_paths, manifest


def upload(stage_dir: Path):
    token = get_token()
    api = HfApi(token=token)

    try:
        api.repo_info(repo_id=REPO_ID, repo_type="dataset")
        print(f"[OK] Repo exists: {REPO_ID}")
    except Exception:
        api.create_repo(repo_id=REPO_ID, repo_type="dataset", exist_ok=True)

    # Build archives into stage, then upload the whole staged folder in one commit.
    archive_paths, manifest = make_archives(stage_dir)

    with open(stage_dir / "manifest.json", "w") as f:
        json.dump(
            {
                "description": "VWF P0 experimental-structure MD (20 ns production)",
                "axes": manifest,
                "note": "xtc trajectories + final gro + topology/mdp per variant; "
                        "dprime_d3 includes md_prod_full_20ns.xtc (merged split run).",
            },
            f,
            indent=2,
        )

    # Remove any other files from a previous staging run.
    print(f"\n=== Upload {len(archive_paths)} archives + manifest -> {REMOTE_SUBFOLDER}/ ===")
    api.upload_folder(
        folder_path=str(stage_dir),
        repo_id=REPO_ID,
        repo_type="dataset",
        path_in_repo=REMOTE_SUBFOLDER,
        commit_message="Upload VWF P0 experimental-structure MD results",
    )
    print("[OK] upload complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        default="/tmp/vwf_p0_md_stage",
        help="Local staging directory for archives (default: %(default)s)",
    )
    args = parser.parse_args()

    stage = Path(args.stage)
    upload(stage)
    print(f"\n[OK] Done. Browse: https://huggingface.co/datasets/{REPO_ID}/tree/main/{REMOTE_SUBFOLDER}")
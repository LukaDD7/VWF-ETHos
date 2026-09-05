#!/usr/bin/env python3
"""Upload the 16-case AlphaGenome paired-track raw results to HuggingFace.

Source: outputs/computational_evidence_16case_20260905/alphagenome/
  - raw/*.npz + raw/*.json sidecars (16 cases x 9 output types = 144 pairs,
    ~2.8 GB compressed; ATAC/CONTACT_MAPS genuinely not returned by the model)
  - derived/paired_track_summaries.csv
  - ledger_flat.csv, run_ledger.jsonl, run_manifest.json, ontology_plan.json
Destination: lucachangretta/VWF -> computational_evidence_16case_20260905/alphagenome/

Grouped by case (one tar.gz per case, ~176 MB each, keeps files count low and
mirrors the per-case access pattern). A manifest.json records per-case file
counts, archive sizes and sha256 of every archive. The API key is NOT part of
any artifact (runner read it only from the environment).

Usage:
    HF_TOKEN=... python scripts/pipeline/upload_16case_alphagenome_huggingface.py
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
ROOT = Path(__file__).parent.parent.parent
SRC_DIR = ROOT / "outputs" / "computational_evidence_16case_20260905" / "alphagenome"

REMOTE_SUBFOLDER = "computational_evidence_16case_20260905/alphagenome"

N_OUTPUT_TYPES = 9  # per case: 9 persisted + 2 not_supported_by_model


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def case_ids(raw_dir: Path) -> list[str]:
    ids = sorted({p.name.split("__")[0] for p in raw_dir.glob("*.npz")})
    return ids


def verify_case(raw_dir: Path, case_id: str) -> dict:
    npzs = sorted(raw_dir.glob(f"{case_id}__*.npz"))
    sidecars = sorted(raw_dir.glob(f"{case_id}__*.json"))
    return {
        "npz_count": len(npzs),
        "sidecar_count": len(sidecars),
        "npz_total_bytes": sum(p.stat().st_size for p in npzs),
        "sidecar_total_bytes": sum(p.stat().st_size for p in sidecars),
    }


def build_archives(stage_dir: Path) -> tuple[list[Path], dict]:
    stage_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = SRC_DIR / "raw"
    ids = case_ids(raw_dir)
    archive_paths: list[Path] = []
    manifest: dict[str, dict] = {}

    for case_id in ids:
        counts = verify_case(raw_dir, case_id)
        ok = (
            counts["npz_count"] == N_OUTPUT_TYPES
            and counts["sidecar_count"] == N_OUTPUT_TYPES
        )
        if not ok:
            print(f"  [SKIP] incomplete case {case_id}: {counts}")
            continue

        arc = stage_dir / f"{case_id}.tar.gz"
        if not arc.exists():
            print(f"  tar {arc.name} ({counts['npz_total_bytes'] / 1e6:.1f} MB npz)")
            members = sorted(raw_dir.glob(f"{case_id}__*"))
            with tarfile.open(arc, "w:gz") as tar:
                for member in members:
                    tar.add(member, arcname=member.name)
        manifest[case_id] = {
            **counts,
            "archive": arc.name,
            "size_bytes": arc.stat().st_size,
            "sha256": sha256(arc),
        }
        archive_paths.append(arc)

    return archive_paths, manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", default="/tmp/vwf_ag16_stage")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    token = os.environ.get("HF_TOKEN")
    if not token and not args.dry_run:
        raise SystemExit("HF_TOKEN environment variable not set")

    stage = Path(args.stage)
    archives, manifest = build_archives(stage)
    print(f"\n{len(archives)} verified case archives, "
          f"total {sum(a.stat().st_size for a in archives) / 1e9:.2f} GB")

    # Shared small files go up as-is (not tarred).
    small_files = [
        SRC_DIR / "derived" / "paired_track_summaries.csv",
        SRC_DIR / "ledger_flat.csv",
        SRC_DIR / "run_ledger.jsonl",
        SRC_DIR / "run_manifest.json",
        SRC_DIR / "ontology_plan.json",
    ]
    small_files = [p for p in small_files if p.exists()]

    (stage / "manifest.json").write_text(
        json.dumps(
            {
                "description": "16-case AlphaGenome paired REF/ALT raw arrays "
                               "(9 output types per case; ATAC and CONTACT_MAPS "
                               "not supported by model for VWF biosamples, "
                               "recorded as not_supported_by_model in ledger)",
                "source": "outputs/computational_evidence_16case_20260905/alphagenome/",
                "case_count": len(archives),
                "output_types_persisted": [
                    "RNA_SEQ", "CAGE", "DNASE", "CHIP_HISTONE", "CHIP_TF",
                    "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS",
                    "PROCAP",
                ],
                "output_types_not_supported": ["ATAC", "CONTACT_MAPS"],
                "array_layout": "npz keys: reference, alternate; "
                                "float32, shape (1048576, n_tracks); "
                                "track order matches JSON sidecar track_metadata",
                "cases": manifest,
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

    print(f"=== Upload {len(archives)} case archives + "
          f"{len(small_files)} shared files + manifest -> {REMOTE_SUBFOLDER}/ ===")
    api.upload_folder(
        folder_path=str(stage),
        repo_id=REPO_ID,
        repo_type="dataset",
        path_in_repo=REMOTE_SUBFOLDER,
        commit_message="Upload 16-case AlphaGenome paired raw arrays (16 cases, 9 output types)",
    )
    for extra in small_files:
        api.upload_file(
            path_or_fileobj=str(extra),
            path_in_repo=f"{REMOTE_SUBFOLDER}/{extra.relative_to(SRC_DIR)}",
            repo_id=REPO_ID,
            repo_type="dataset",
            commit_message=f"Add {extra.name} (AlphaGenome 16-case run)",
        )
    print(f"[OK] https://huggingface.co/datasets/{REPO_ID}/tree/main/{REMOTE_SUBFOLDER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

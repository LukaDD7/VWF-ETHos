#!/usr/bin/env python3
"""Download legacy Boltz A1-GPIb result files from the VWF Hugging Face dataset."""

from __future__ import annotations

import argparse
import time
import urllib.request
from pathlib import Path


HF_BASE = "https://huggingface.co/datasets/lucachangretta/VWF/resolve/main"
HF_PREFIX = "boltz2_a1_gpiba_results"


def download(url: str, dest: Path, retries: int = 3) -> bool:
    if dest.exists() and dest.stat().st_size > 100:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                data = resp.read()
            if len(data) < 100:
                raise ValueError(f"too small: {len(data)} bytes")
            dest.write_bytes(data)
            return True
        except Exception as exc:
            if attempt == retries - 1:
                print(f"[WARN] failed: {url} ({exc})")
                return False
            time.sleep(2**attempt)
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--yaml-dir", default="output/boltz2_a1_gpiba")
    parser.add_argument("--out-dir", default="output/boltz2_a1_gpiba_results")
    parser.add_argument("--models", default="0")
    parser.add_argument("--file-kinds", default="cif,pae,plddt,confidence")
    parser.add_argument(
        "--jobs",
        default="",
        help="Optional comma-separated job names or variant IDs, e.g. VWF_R1306W,VWF_V1316M_vs_GPIb_alpha",
    )
    args = parser.parse_args()

    yaml_dir = Path(args.yaml_dir)
    out_dir = Path(args.out_dir)
    models = [int(x) for x in args.models.split(",") if x.strip()]
    kinds = {x.strip() for x in args.file_kinds.split(",") if x.strip()}
    jobs = sorted(p.stem for p in yaml_dir.glob("*.yaml"))
    if args.jobs.strip():
        wanted = []
        available = set(jobs)
        for item in args.jobs.split(","):
            item = item.strip()
            if not item:
                continue
            job = item if item.endswith("_vs_GPIb_alpha") else f"{item}_vs_GPIb_alpha"
            if job not in available:
                print(f"[WARN] requested job not in {yaml_dir}: {job}")
            else:
                wanted.append(job)
        jobs = wanted
    if not jobs:
        raise SystemExit(f"No YAML jobs found in {yaml_dir}")

    total = 0
    ok = 0
    for job in jobs:
        pred_rel = f"{HF_PREFIX}/boltz_results_{job}/predictions/{job}"
        pred_dir = out_dir / f"boltz_results_{job}" / "predictions" / job
        for mid in models:
            files = []
            if "cif" in kinds:
                files.append((f"{job}_model_{mid}.cif", pred_dir / f"{job}_model_{mid}.cif"))
            if "pae" in kinds:
                files.append((f"pae_{job}_model_{mid}.npz", pred_dir / f"pae_{job}_model_{mid}.npz"))
            if "plddt" in kinds:
                files.append((f"plddt_{job}_model_{mid}.npz", pred_dir / f"plddt_{job}_model_{mid}.npz"))
            if "confidence" in kinds:
                files.append((f"confidence_{job}_model_{mid}.json", pred_dir / f"confidence_{job}_model_{mid}.json"))
            for filename, dest in files:
                total += 1
                url = f"{HF_BASE}/{pred_rel}/{filename}"
                if download(url, dest):
                    ok += 1
        print(f"{job}: downloaded/available")
    print(f"Done: {ok}/{total} files available under {out_dir}")


if __name__ == "__main__":
    main()

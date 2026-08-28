#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time

import httpx


ROOT = Path(__file__).resolve().parents[1]


def write_snapshot(path: Path, payload: dict) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(payload, indent=2, ensure_ascii=False)
    path.write_text(text + "\n", encoding="utf-8")
    return sha256(text.encode("utf-8")).hexdigest()


def download_clingen(output_dir: Path) -> dict:
    url = "https://erepo.genome.network/evrepo/api/classifications"
    params = {"gene": "VWF", "matchMode": "exact", "matchLimit": 1000}
    response = httpx.get(url, params=params, timeout=60)
    response.raise_for_status()
    payload = response.json()
    path = output_dir / "clingen_vwf_erepo.json"
    digest = write_snapshot(path, payload)
    return {
        "tool": "clingen_erepo",
        "url": str(response.url),
        "path": str(path),
        "sha256": digest,
        "records": len(payload.get("variantInterpretations", [])),
    }


def download_clinvar(output_dir: Path, retmax: int) -> dict:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    search = httpx.get(
        f"{base}/esearch.fcgi",
        params={"db": "clinvar", "term": "VWF[Gene]", "retmode": "json", "retmax": retmax},
        timeout=60,
    )
    search.raise_for_status()
    ids = search.json().get("esearchresult", {}).get("idlist", [])
    time.sleep(0.5)
    result: dict = {"result": {"uids": []}}
    for offset in range(0, len(ids), 200):
        batch = ids[offset : offset + 200]
        summary = httpx.post(
            f"{base}/esummary.fcgi",
            data={"db": "clinvar", "id": ",".join(batch), "retmode": "json"},
            timeout=60,
        )
        summary.raise_for_status()
        payload_batch = summary.json()
        result["result"]["uids"].extend(payload_batch.get("result", {}).get("uids", []))
        for uid, value in payload_batch.get("result", {}).items():
            if uid not in {"uids"}:
                result["result"][uid] = value
        time.sleep(0.5)
    payload = result
    path = output_dir / "clinvar_vwf_esummary.json"
    digest = write_snapshot(path, payload)
    return {
        "tool": "clinvar_eutils",
        "url": str(summary.url),
        "path": str(path),
        "sha256": digest,
        "records": len(ids),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Download small VWF-specific public database snapshots.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/external/vwf_biomedical_snapshots")
    parser.add_argument("--tools", nargs="+", choices=["clingen", "clinvar"], default=["clingen", "clinvar"])
    parser.add_argument("--clinvar-retmax", type=int, default=1000)
    args = parser.parse_args()
    records = []
    if "clingen" in args.tools:
        records.append(download_clingen(args.output_dir))
    if "clinvar" in args.tools:
        records.append(download_clinvar(args.output_dir, args.clinvar_retmax))
    manifest = {"downloaded_at": datetime.now(timezone.utc).isoformat(), "snapshots": records}
    path = args.output_dir / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

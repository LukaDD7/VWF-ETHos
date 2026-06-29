#!/usr/bin/env python3
"""Summarize completed Type 2M LOF MD runs without touching active jobs."""

from __future__ import annotations

import argparse
import csv
import re
from datetime import datetime
from pathlib import Path


TIME_RE = re.compile(r"\[(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})\]")
PERF_RE = re.compile(r"Performance:\s+([0-9.]+)\s+ns/day")


def read_variants(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]


def parse_runner_times(prod_dir: Path, variant: str) -> dict[str, str | float | None]:
    starts: list[datetime] = []
    completes: list[datetime] = []
    for log in sorted(prod_dir.glob("*.log")) + sorted(prod_dir.glob("*.stdout")):
        text = log.read_text(errors="ignore")
        for line in text.splitlines():
            if variant in line and ("start" in line.lower() or "resume" in line.lower() or "complete" in line.lower()):
                m = TIME_RE.search(line)
                if not m:
                    continue
                ts = datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
                if "complete" in line.lower():
                    completes.append(ts)
                else:
                    starts.append(ts)
    start = min(starts) if starts else None
    end = max(completes) if completes else None
    hours = round((end - start).total_seconds() / 3600, 2) if start and end else None
    return {
        "runner_start": start.isoformat(sep=" ") if start else None,
        "runner_complete": end.isoformat(sep=" ") if end else None,
        "wall_hours_from_runner": hours,
    }


def parse_performance(prod_dir: Path) -> float | None:
    values: list[float] = []
    for log in sorted(prod_dir.glob("md_prod*.log")) + sorted(prod_dir.glob("md_prod*.stdout")):
        text = log.read_text(errors="ignore")
        for m in PERF_RE.finditer(text):
            values.append(float(m.group(1)))
    return round(values[-1], 3) if values else None


def file_time(path: Path) -> str | None:
    if not path.exists():
        return None
    return datetime.fromtimestamp(path.stat().st_mtime).isoformat(sep=" ", timespec="seconds")


def summarize_group(root: Path, variants: list[str], md_subdir: str, system: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for variant in variants:
        prod_dir = root / variant / md_subdir
        complete = (prod_dir / "md_prod.gro").exists() or bool(list(prod_dir.glob("md_prod.part*.gro")))
        started = (prod_dir / "md_prod.cpt").exists() or bool(list(prod_dir.glob("md_prod.part*.cpt")))
        xtc = prod_dir / "md_prod.xtc"
        if not xtc.exists():
            parts = sorted(prod_dir.glob("md_prod.part*.xtc"))
            xtc = parts[-1] if parts else xtc
        gro = prod_dir / "md_prod.gro"
        if not gro.exists():
            parts = sorted(prod_dir.glob("md_prod.part*.gro"))
            gro = parts[-1] if parts else gro
        perf = parse_performance(prod_dir)
        rows.append({
            "system": system,
            "variant": variant,
            "status": "complete" if complete else ("running_or_checkpointed" if started else "pending"),
            "prod_dir": str(prod_dir),
            "final_gro": str(gro) if gro.exists() else "",
            "trajectory": str(xtc) if xtc.exists() else "",
            "tpr_exists": (prod_dir / "md_prod.tpr").exists(),
            "xtc_size_mb": round(xtc.stat().st_size / 1024 / 1024, 2) if xtc.exists() else None,
            "gro_mtime": file_time(gro),
            "cpt_mtime": file_time(prod_dir / "md_prod.cpt"),
            "ns_per_day": perf,
            **parse_runner_times(prod_dir, variant),
        })
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def estimate_remaining(rows: list[dict[str, object]]) -> dict[str, object]:
    completed = [r for r in rows if r["status"] == "complete"]
    pending = [r for r in rows if r["status"] != "complete"]
    durations = [float(r["wall_hours_from_runner"]) for r in completed if r.get("wall_hours_from_runner")]
    per_run = sorted(durations)[len(durations) // 2] if durations else None
    return {
        "total": len(rows),
        "complete": len(completed),
        "remaining": len(pending),
        "median_wall_hours_per_run": round(per_run, 2) if per_run else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--inputs", type=Path, default=Path("output/type2m_lof_server_inputs_2026-06-28"))
    parser.add_argument("--output", type=Path, default=Path("output/type2m_lof_md_analysis_2026-06-29"))
    args = parser.parse_args()

    root = args.root.resolve()
    inputs = (root / args.inputs).resolve()
    out = (root / args.output).resolve()

    a1_variants = read_variants(inputs / "md_a1_gpiba" / "a1_gpiba_all_recommended_variants.txt")
    a1_p0 = set(read_variants(inputs / "md_a1_gpiba" / "a1_gpiba_p0_plus_anchor_variants.txt"))
    closed_variants = read_variants(inputs / "md_7a6o_closed_state" / "7a6o_2m_closed_state_new_variants.txt")

    a1_rows = summarize_group(root / "output/gromacs_md_a1_gpiba", a1_variants, "md_a1_gpiba", "a1_gpiba")
    for row in a1_rows:
        row["queue"] = "P0" if row["variant"] in a1_p0 else "all_recommended"
    closed_rows = summarize_group(root / "output/gromacs_md_autoinhib", closed_variants, "md_7a6o", "7a6o_closed_state")

    write_csv(out / "a1_gpiba_completed_and_running_summary.csv", a1_rows)
    write_csv(out / "7a6o_closed_state_completed_and_running_summary.csv", closed_rows)
    write_csv(out / "combined_completed_and_running_summary.csv", a1_rows + closed_rows)

    summary_rows = [
        {"system": "a1_gpiba", **estimate_remaining(a1_rows)},
        {"system": "7a6o_closed_state", **estimate_remaining(closed_rows)},
    ]
    write_csv(out / "queue_progress_summary.csv", summary_rows)

    readme = out / "README.md"
    readme.write_text(
        "# Type 2M LOF MD partial analysis\n\n"
        "Generated from completed or checkpointed production MD runs only. Large trajectory files are not copied here.\n\n"
        "Files:\n"
        "- `queue_progress_summary.csv`: per-queue completion counts and median completed wall time.\n"
        "- `a1_gpiba_completed_and_running_summary.csv`: A1-GPIb queue status and completed-run metadata.\n"
        "- `7a6o_closed_state_completed_and_running_summary.csv`: 7A6O closed-state queue status and completed-run metadata.\n"
        "- `combined_completed_and_running_summary.csv`: combined view.\n"
        "- `7a6o_completed_qc/`: GROMACS RMSD/contact QC for completed 7A6O closed-state runs, when generated.\n",
        encoding="utf-8",
    )
    print(f"Wrote {out}")
    for row in summary_rows:
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

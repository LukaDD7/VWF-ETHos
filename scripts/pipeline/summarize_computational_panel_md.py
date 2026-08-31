#!/usr/bin/env python3
"""Build lightweight derived summaries for the 2026-08-29 targeted MD panel."""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path
from statistics import mean


DEFAULT_MD_DIR = Path("outputs/computational_panel_20260829/md")
SELECTIONS = ("AIM_all", "N_AIM", "C_AIM")


def parse_xvg(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "@")):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            try:
                rows.append((float(fields[0]), float(fields[1])))
            except ValueError:
                continue
    return rows


def summarize(rows: list[tuple[float, float]]) -> dict[str, float | int | None]:
    if not rows:
        return {key: None for key in (
            "frames", "mean", "first0_5", "mid20_30", "tail40_50",
            "final", "min", "max", "tail_minus_first",
        )}
    values = [value for _, value in rows]
    first = [value for time, value in rows if time <= 5]
    middle = [value for time, value in rows if 20 <= time <= 30]
    tail = [value for time, value in rows if time >= 40]
    first_mean = mean(first) if first else None
    tail_mean = mean(tail) if tail else None
    return {
        "frames": len(rows),
        "mean": round(mean(values), 6),
        "first0_5": round(first_mean, 6) if first_mean is not None else None,
        "mid20_30": round(mean(middle), 6) if middle else None,
        "tail40_50": round(tail_mean, 6) if tail else None,
        "final": round(rows[-1][1], 6),
        "min": round(min(values), 6),
        "max": round(max(values), 6),
        "tail_minus_first": round(tail_mean - first_mean, 6)
        if first_mean is not None and tail_mean is not None else None,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_labels(qc_path: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    if not qc_path.exists():
        return labels
    with qc_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            labels[row["variant"]] = row.get("label") or "?"
    return labels


def collect_rmsd(md_dir: Path, labels: dict[str, str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    summary: list[dict[str, object]] = []
    timeseries: list[dict[str, object]] = []
    for path in sorted((md_dir / "rmsd").glob("*_backbone_rmsd.xvg")):
        variant = path.name.removesuffix("_backbone_rmsd.xvg")
        rows = parse_xvg(path)
        stats = summarize(rows)
        summary.append({"variant": variant, "label": labels.get(variant, "?"), **stats})
        timeseries.extend(
            {"variant": variant, "label": labels.get(variant, "?"), "time_ns": time, "rmsd_nm": value}
            for time, value in rows
        )
    return summary, timeseries


def collect_distances(md_dir: Path, labels: dict[str, str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    summary: list[dict[str, object]] = []
    timeseries: list[dict[str, object]] = []
    pattern = re.compile(r"^(?P<variant>.+)_(?P<selection>AIM_all|N_AIM|C_AIM)_a1_nonlocal_mindist\.xvg$")
    for path in sorted((md_dir / "contacts").glob("*_a1_nonlocal_mindist.xvg")):
        match = pattern.match(path.name)
        if not match:
            continue
        variant = match.group("variant")
        selection = match.group("selection")
        rows = parse_xvg(path)
        stats = summarize(rows)
        summary.append({
            "variant": variant,
            "label": labels.get(variant, "?"),
            "selection": selection,
            **stats,
        })
        timeseries.extend(
            {
                "variant": variant,
                "label": labels.get(variant, "?"),
                "selection": selection,
                "time_ns": time,
                "min_distance_nm": value,
            }
            for time, value in rows
        )
    return summary, timeseries


def add_wt_deltas(rows: list[dict[str, object]], value_key: str, wt_variant: str) -> None:
    grouped: dict[object, dict[str, dict[str, object]]] = {}
    for row in rows:
        group_key = row.get("selection", "all")
        grouped.setdefault(group_key, {})[str(row.get("variant"))] = row
    for variants in grouped.values():
        wt = variants.get(wt_variant)
        if wt is None:
            continue
        for row in variants.values():
            row[f"{value_key}_delta_tail_vs_WT"] = (
                round(float(row[value_key]) - float(wt[value_key]), 6)
                if row.get(value_key) is not None and wt.get(value_key) is not None
                else None
            )


def build_feature_matrix(
    qc_path: Path,
    rmsd_summary: list[dict[str, object]],
    distance_summary: list[dict[str, object]],
    wt_variant: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    rmsd_by_variant = {row["variant"]: row for row in rmsd_summary}
    distance_by_key = {(row["variant"], row["selection"]): row for row in distance_summary}
    with qc_path.open(newline="", encoding="utf-8") as handle:
        for qc in csv.DictReader(handle):
            variant = qc["variant"]
            rmsd = rmsd_by_variant.get(variant, {})
            row: dict[str, object] = {
                "variant": variant,
                "label": qc.get("label") or "?",
                "frames": qc.get("frames"),
                "timestep_ps": qc.get("timestep_ps"),
                "atoms": qc.get("atoms"),
                "rmsd_mean_nm": rmsd.get("mean"),
                "rmsd_tail40_50_mean_nm": rmsd.get("tail40_50"),
                "rmsd_final_nm": rmsd.get("final"),
            }
            for selection in SELECTIONS:
                distance = distance_by_key.get((variant, selection), {})
                row[f"{selection}_min_distance_tail40_50_nm"] = distance.get("tail40_50")
                row[f"{selection}_min_distance_mean_nm"] = distance.get("mean")
                row[f"{selection}_min_distance_final_nm"] = distance.get("final")
            rows.append(row)

    wt = next((row for row in rows if row["variant"] == wt_variant), None)
    if wt is not None:
        metric_keys = [key for key in rows[0] if key not in {"variant", "label", "frames", "timestep_ps", "atoms"}]
        for row in rows:
            for key in metric_keys:
                row[f"{key}_delta_vs_WT"] = (
                    round(float(row[key]) - float(wt[key]), 6)
                    if row.get(key) not in (None, "") and wt.get(key) not in (None, "")
                    else None
                )
    return rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(md_dir: Path, derived_paths: list[Path]) -> None:
    rows: list[dict[str, object]] = []
    for path in derived_paths + [
        *(md_dir / "contacts").glob("*.xvg"),
        *(md_dir / "index").glob("*.ndx"),
        *(md_dir / "rmsd").glob("*.xvg"),
        *(md_dir / "trajectories").glob("*.xtc"),
    ]:
        if not path.is_file() or path.name.startswith("#"):
            continue
        rows.append({
            "relative_path": str(path.relative_to(md_dir)),
            "artifact_type": "derived_summary" if path in derived_paths else path.suffix.lstrip(".").lower(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "git_policy": "commit" if path.suffix in {".csv", ".xvg", ".ndx"} else "local_only",
        })
    write_csv(md_dir / "artifact_manifest.csv", sorted(rows, key=lambda row: str(row["relative_path"])))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--md-dir", type=Path, default=DEFAULT_MD_DIR)
    parser.add_argument("--wt-variant", default="PANEL_WT_MATCHED")
    args = parser.parse_args()
    md_dir = args.md_dir.resolve()
    labels = read_labels(md_dir / "qc_summary.csv")

    rmsd_summary, rmsd_timeseries = collect_rmsd(md_dir, labels)
    distance_summary, distance_timeseries = collect_distances(md_dir, labels)
    add_wt_deltas(rmsd_summary, "tail40_50", args.wt_variant)
    add_wt_deltas(distance_summary, "tail40_50", args.wt_variant)
    feature_matrix = build_feature_matrix(md_dir / "qc_summary.csv", rmsd_summary, distance_summary, args.wt_variant)

    outputs = [
        md_dir / "rmsd_timeseries.csv",
        md_dir / "aim_a1_min_distance_summary.csv",
        md_dir / "aim_a1_min_distance_timeseries.csv",
        md_dir / "md_feature_matrix.csv",
    ]
    write_csv(outputs[0], rmsd_timeseries)
    write_csv(outputs[1], distance_summary)
    write_csv(outputs[2], distance_timeseries)
    write_csv(outputs[3], feature_matrix)
    write_manifest(md_dir, outputs)
    for path in outputs + [md_dir / "artifact_manifest.csv"]:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

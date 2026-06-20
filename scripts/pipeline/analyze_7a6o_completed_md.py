#!/usr/bin/env python3
"""Analyze completed 7A6O AIM-A1 GROMACS production runs.

This script is intentionally independent of gemmi. It uses GROMACS tools to:
  1. discover completed WT/mutant production runs,
  2. concatenate split/noappend md_prod trajectories into QC copies,
  3. verify frame counts with gmx check,
  4. compute backbone RMSD,
  5. compute AIM-A1 contact counts with gmx select + gmx mindist.

It never edits original production files. Outputs are written under:
  output/gromacs_md_autoinhib/analysis_completed_7a6o/
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Iterable

DEFAULT_VARIANTS = [
    "WT", "WT_r1", "WT_r2", "WT_r3",
    "R1306W", "R1306Q", "R1308C", "I1309V", "S1310F", "W1313C", "V1314F", "V1316M",
    "R1374C", "R1374H", "G1324S",
    "P1337L", "R1341Q", "R1341W",
    "C1458R", "C1458R_r1", "C1458R_r2", "C1458R_r3",
]

LABELS = {
    "WT": "WT", "WT_r1": "WT", "WT_r2": "WT", "WT_r3": "WT",
    "R1306W": "2B", "R1306Q": "2B", "R1308C": "2B", "I1309V": "2B", "S1310F": "2B",
    "W1313C": "2B", "V1314F": "2B", "V1316M": "2B",
    "R1374C": "2M", "R1374H": "2M", "G1324S": "2M",
    "P1337L": "?", "R1341Q": "?", "R1341W": "?",
    "C1458R": "patient", "C1458R_r1": "patient", "C1458R_r2": "patient", "C1458R_r3": "patient",
}

SELECTIONS = {
    "AIM_all": ('"AIM" resnr 1262 to 1271 or resnr 1460 to 1466', "AIM"),
    "N_AIM": ('"N_AIM" resnr 1262 to 1271', "N_AIM"),
    "C_AIM": ('"C_AIM" resnr 1460 to 1466', "C_AIM"),
}
A1_NONLOCAL = '"A1_nonlocal" resnr 1280 to 1451'


@dataclass
class Run:
    variant: str
    label: str
    prod_dir: Path
    tpr: Path
    final_gro: Path
    xtc_parts: list[Path]


def run_cmd(cmd: list[str], cwd: Path, stdin: str | None = None, log: Path | None = None) -> str:
    env = os.environ.copy()
    env["GMX_RESET_COUNTERS"] = "0"
    result = subprocess.run(cmd, cwd=cwd, input=stdin, text=True, capture_output=True, env=env)
    output = result.stdout + result.stderr
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(output)
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(cmd)}\n{output[-4000:]}")
    return output


def parse_xvg(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    if not path.exists():
        return rows
    for line in path.read_text(errors="ignore").splitlines():
        if not line.strip() or line.startswith(("#", "@")):
            continue
        parts = line.split()
        if len(parts) >= 2:
            rows.append((float(parts[0]), float(parts[1])))
    return rows


def summarize_series(rows: list[tuple[float, float]]) -> dict[str, float | int | None]:
    if not rows:
        return {
            "frames": 0, "mean": None, "first0_5": None, "mid20_30": None,
            "tail40_50": None, "final": None, "min": None, "max": None,
            "tail_minus_first": None,
        }
    vals = [v for _, v in rows]
    first = [v for t, v in rows if t <= 5]
    mid = [v for t, v in rows if 20 <= t <= 30]
    tail = [v for t, v in rows if t >= 40]
    first_mean = mean(first) if first else None
    tail_mean = mean(tail) if tail else None
    return {
        "frames": len(rows),
        "mean": round(mean(vals), 4),
        "first0_5": round(first_mean, 4) if first_mean is not None else None,
        "mid20_30": round(mean(mid), 4) if mid else None,
        "tail40_50": round(tail_mean, 4) if tail_mean is not None else None,
        "final": round(rows[-1][1], 4),
        "min": round(min(vals), 4),
        "max": round(max(vals), 4),
        "tail_minus_first": round(tail_mean - first_mean, 4) if first_mean is not None and tail_mean is not None else None,
    }


def discover_run(root: Path, variant: str) -> Run | None:
    if variant == "WT":
        prod = root / "7A6O_WT" / "md_7a6o"
        tpr = prod / "md_prod.tpr"
        final = prod / "md_prod_run.part0002.gro"
        parts = [p for p in [prod / "md_prod.xtc", prod / "md_prod_run.part0002.xtc"] if p.exists() and p.stat().st_size > 0]
    else:
        prod = root / variant / "md_7a6o"
        tpr = prod / "md_prod.tpr"
        final_candidates = sorted(list(prod.glob("md_prod.part*.gro")) + list(prod.glob("md_prod.gro")))
        # Exclude probe/test files. Only official md_prod and md_prod.partNNNN files are accepted.
        final_candidates = [p for p in final_candidates if re.match(r"md_prod(\.part\d+)?\.gro$", p.name)]
        final = final_candidates[-1] if final_candidates else prod / "md_prod.gro"
        part_candidates = [prod / "md_prod.xtc"] + sorted(prod.glob("md_prod.part*.xtc"))
        parts = [p for p in part_candidates if p.exists() and p.stat().st_size > 0 and re.match(r"md_prod(\.part\d+)?\.xtc$", p.name)]
    if not (prod.exists() and tpr.exists() and final.exists() and parts):
        return None
    return Run(variant=variant, label=LABELS.get(variant, "?"), prod_dir=prod, tpr=tpr, final_gro=final, xtc_parts=parts)


def parse_gmx_check(log_text: str) -> dict[str, int | float | None]:
    frames = None
    dt = None
    atoms = None
    for line in log_text.splitlines():
        m = re.search(r"# Atoms\s+(\d+)", line)
        if m:
            atoms = int(m.group(1))
        m = re.match(r"Coords\s+(\d+)\s+([0-9.]+)", line.strip())
        if m:
            frames = int(m.group(1))
            dt = float(m.group(2))
    return {"frames": frames, "timestep_ps": dt, "atoms": atoms}


def concat_and_check(gmx: str, run: Run, out_dir: Path, force: bool) -> tuple[Path, dict[str, int | float | None]]:
    concat = out_dir / "trajectories" / f"{run.variant}_prod_concat.xtc"
    concat.parent.mkdir(parents=True, exist_ok=True)
    if force or not concat.exists():
        cmd = [gmx, "trjcat", "-f", *map(str, run.xtc_parts), "-o", str(concat), "-cat"]
        run_cmd(cmd, cwd=Path.cwd(), log=out_dir / "logs" / f"{run.variant}_trjcat.log")
    check_out = run_cmd([gmx, "check", "-f", str(concat)], cwd=Path.cwd(), log=out_dir / "logs" / f"{run.variant}_check.log")
    return concat, parse_gmx_check(check_out)


def compute_rmsd(gmx: str, run: Run, traj: Path, out_dir: Path, force: bool) -> dict[str, float | int | None]:
    xvg = out_dir / "rmsd" / f"{run.variant}_backbone_rmsd.xvg"
    xvg.parent.mkdir(parents=True, exist_ok=True)
    if force or not xvg.exists():
        run_cmd(
            [gmx, "rms", "-s", str(run.tpr), "-f", str(traj), "-o", str(xvg), "-tu", "ns"],
            cwd=Path.cwd(), stdin="Backbone\nBackbone\n", log=out_dir / "logs" / f"{run.variant}_rmsd.log",
        )
    rows = parse_xvg(xvg)
    summary = summarize_series(rows)
    return {
        "rmsd_frames": summary["frames"],
        "rmsd_mean_nm": summary["mean"],
        "rmsd_max_nm": summary["max"],
        "rmsd_tail40_50_mean_nm": summary["tail40_50"],
        "rmsd_final_nm": summary["final"],
    }


def compute_contacts(gmx: str, run: Run, traj: Path, out_dir: Path, force: bool) -> tuple[dict[str, float | int | None], list[dict[str, float | str]]]:
    summary: dict[str, float | int | None] = {}
    timeseries_rows: list[dict[str, float | str]] = []
    for key, (aim_sel, group_name) in SELECTIONS.items():
        ndx = out_dir / "index" / f"{run.variant}_{key}_a1_nonlocal.ndx"
        numcont = out_dir / "contacts" / f"{run.variant}_{key}_a1_nonlocal_numcont.xvg"
        mindist = out_dir / "contacts" / f"{run.variant}_{key}_a1_nonlocal_mindist.xvg"
        ndx.parent.mkdir(parents=True, exist_ok=True)
        numcont.parent.mkdir(parents=True, exist_ok=True)
        if force or not ndx.exists():
            run_cmd(
                [gmx, "select", "-s", str(run.tpr), "-on", str(ndx), "-select", f"{aim_sel}; {A1_NONLOCAL}"],
                cwd=Path.cwd(), log=out_dir / "logs" / f"{run.variant}_{key}_select.log",
            )
        if force or not numcont.exists():
            run_cmd(
                [gmx, "mindist", "-s", str(run.tpr), "-f", str(traj), "-n", str(ndx),
                 "-od", str(mindist), "-on", str(numcont), "-d", "0.45", "-group", "-tu", "ns", "-xvg", "none"],
                cwd=Path.cwd(), stdin=f"{group_name}\nA1_nonlocal\n", log=out_dir / "logs" / f"{run.variant}_{key}_mindist.log",
            )
        rows = parse_xvg(numcont)
        stats = summarize_series(rows)
        prefix = f"{key}_contacts"
        for stat_key, value in stats.items():
            summary[f"{prefix}_{stat_key}"] = value
        for t, value in rows:
            timeseries_rows.append({"variant": run.variant, "label": run.label, "selection": key, "time_ns": t, "contacts": value})
    return summary, timeseries_rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fields: list[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_variants(raw: str | None) -> list[str]:
    if not raw:
        return DEFAULT_VARIANTS
    return [v.strip() for v in raw.split(",") if v.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", default="output/gromacs_md_autoinhib", help="7A6O GROMACS output root")
    parser.add_argument("--output", default=None, help="analysis output dir")
    parser.add_argument("--gmx", default="envs/gromacs/bin.AVX2_256/gmx", help="GROMACS executable")
    parser.add_argument("--variants", default=None, help="comma-separated variants; default WT + known 14")
    parser.add_argument("--force", action="store_true", help="rebuild intermediate trajectory/contact files")
    args = parser.parse_args()

    root = Path(args.input).resolve()
    out_dir = Path(args.output).resolve() if args.output else root / "analysis_completed_7a6o"
    gmx = str(Path(args.gmx).resolve()) if Path(args.gmx).exists() else args.gmx
    variants = parse_variants(args.variants)

    runs = []
    skipped = []
    for variant in variants:
        run = discover_run(root, variant)
        if run:
            runs.append(run)
        else:
            skipped.append(variant)

    print(f"Input: {root}")
    print(f"Output: {out_dir}")
    print(f"Completed runs found: {len(runs)}")
    if skipped:
        print("Skipped incomplete/not found:", ",".join(skipped))

    qc_rows: list[dict] = []
    contact_summary_rows: list[dict] = []
    contact_timeseries_rows: list[dict] = []

    for run in runs:
        print(f"Analyzing {run.variant} ({run.label})")
        traj, check = concat_and_check(gmx, run, out_dir, args.force)
        rmsd = compute_rmsd(gmx, run, traj, out_dir, args.force)
        contacts, ts = compute_contacts(gmx, run, traj, out_dir, args.force)
        qc_rows.append({
            "variant": run.variant,
            "label": run.label,
            "prod_dir": str(run.prod_dir),
            "final_gro": str(run.final_gro),
            "concat_xtc": str(traj),
            "xtc_parts": ";".join(str(p) for p in run.xtc_parts),
            **check,
            **rmsd,
        })
        contact_summary_rows.append({"variant": run.variant, "label": run.label, **contacts})
        contact_timeseries_rows.extend(ts)

    if contact_summary_rows:
        wt = next((row for row in contact_summary_rows if row["variant"] == "WT"), None)
        if wt:
            for row in contact_summary_rows:
                for key in SELECTIONS:
                    tail_key = f"{key}_contacts_tail40_50"
                    delta_key = f"{key}_contacts_delta_tail_vs_WT"
                    if row.get(tail_key) is not None and wt.get(tail_key) is not None:
                        row[delta_key] = round(float(row[tail_key]) - float(wt[tail_key]), 4)

    write_csv(out_dir / "qc_summary.csv", qc_rows)
    write_csv(out_dir / "aim_a1_contacts_summary.csv", contact_summary_rows)
    write_csv(out_dir / "aim_a1_contacts_timeseries.csv", contact_timeseries_rows)
    print(f"Written: {out_dir / 'qc_summary.csv'}")
    print(f"Written: {out_dir / 'aim_a1_contacts_summary.csv'}")
    print(f"Written: {out_dir / 'aim_a1_contacts_timeseries.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

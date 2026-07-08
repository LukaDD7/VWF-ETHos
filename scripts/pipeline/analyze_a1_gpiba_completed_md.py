#!/usr/bin/env python3
"""Analyze completed A1-GPIbalpha complex MD runs for classifier-facing features.

This script is meant to run on the A40/H200 server after A1-GPIbalpha production
MD has completed. It extracts directional interface features from trajectories,
instead of using backbone RMSD as a classifier proxy.

Outputs:
  output/type2m_lof_md_analysis_2026-06-29/a1_gpiba_interface_qc/
    - a1_gpiba_interface_summary.csv
    - a1_gpiba_interface_timeseries.csv
    - a1_gpiba_classifier_features.csv
    - logs/*.log

Feature direction:
  higher contact retention / retained_z = interface is retained
  higher contact loss / loss_z          = interface is lost, 2M/LOF-compatible

The script never edits production files. It creates only analysis xvg/csv/log
files under the output directory.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[2]

DEFAULT_STATUS = ROOT / "output/type2m_lof_md_analysis_2026-06-29/a1_gpiba_completed_and_running_summary.csv"
DEFAULT_OUT = ROOT / "output/type2m_lof_md_analysis_2026-06-29/a1_gpiba_interface_qc"
DEFAULT_LABELS = ROOT / "output/labeled_variants_all.csv"
DEFAULT_ROOT = ROOT / "output/gromacs_md_a1_gpiba"

# 1SQ0-derived A1-GPIbalpha MD. FoldX uses VWF native position = PDB resnr + 763.
DEFAULT_A1_SELECTION = 'group "Protein" and resnr 505 to 703'
DEFAULT_GPIBA_SELECTION = 'group "Protein" and resnr 1 to 267'


@dataclass
class Run:
    variant: str
    label: str
    prod_dir: Path
    tpr: Path
    final_gro: Path
    xtc_parts: list[Path]
    status: str = "complete"


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fields: list[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_labels(labels_path: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    for row in read_csv(labels_path):
        aa = row.get("aa_change") or ""
        lab = row.get("labels") or row.get("true_label") or ""
        conflict = str(row.get("conflict", "")).lower() == "true"
        if aa and lab and not conflict:
            labels[aa] = lab
    return labels


def parse_variants(raw: str | None) -> set[str] | None:
    if not raw:
        return None
    vals = {v.strip() for v in raw.split(",") if v.strip()}
    return vals or None


def file_candidates(prod_dir: Path, stem: str, suffix: str) -> list[Path]:
    direct = prod_dir / f"{stem}{suffix}"
    parts = sorted(prod_dir.glob(f"{stem}.part*{suffix}"))
    out = []
    if direct.exists() and direct.stat().st_size > 0:
        out.append(direct)
    out.extend(p for p in parts if p.exists() and p.stat().st_size > 0)
    return out


def discover_from_status(status_path: Path, labels: dict[str, str], variants: set[str] | None) -> tuple[list[Run], list[dict]]:
    runs: list[Run] = []
    skipped: list[dict] = []
    for row in read_csv(status_path):
        variant = row.get("variant", "")
        if not variant or (variants and variant not in variants):
            continue
        status = row.get("status", "")
        prod_dir = Path(row.get("prod_dir", ""))
        if status != "complete":
            skipped.append({"variant": variant, "reason": f"status={status or 'unknown'}", "prod_dir": str(prod_dir)})
            continue
        tpr = prod_dir / "md_prod.tpr"
        final = Path(row.get("final_gro") or prod_dir / "md_prod.gro")
        if not final.exists():
            gro_parts = file_candidates(prod_dir, "md_prod", ".gro")
            final = gro_parts[-1] if gro_parts else final
        xtc_parts = file_candidates(prod_dir, "md_prod", ".xtc")
        if not (prod_dir.exists() and tpr.exists() and final.exists() and xtc_parts):
            skipped.append({"variant": variant, "reason": "missing prod_dir/tpr/final_gro/xtc on this machine", "prod_dir": str(prod_dir)})
            continue
        runs.append(Run(variant=variant, label=labels.get(variant, "?"), prod_dir=prod_dir, tpr=tpr, final_gro=final, xtc_parts=xtc_parts, status=status))
    return runs, skipped


def discover_from_root(root: Path, labels: dict[str, str], variants: set[str] | None) -> tuple[list[Run], list[dict]]:
    runs: list[Run] = []
    skipped: list[dict] = []
    candidates = sorted(p.parent for p in root.glob("*/md_a1_gpiba/md_prod.tpr"))
    for prod_dir in candidates:
        variant = prod_dir.parent.name
        if variants and variant not in variants:
            continue
        tpr = prod_dir / "md_prod.tpr"
        gro_parts = file_candidates(prod_dir, "md_prod", ".gro")
        xtc_parts = file_candidates(prod_dir, "md_prod", ".xtc")
        if not (tpr.exists() and gro_parts and xtc_parts):
            skipped.append({"variant": variant, "reason": "incomplete files", "prod_dir": str(prod_dir)})
            continue
        runs.append(Run(variant=variant, label=labels.get(variant, "?"), prod_dir=prod_dir, tpr=tpr, final_gro=gro_parts[-1], xtc_parts=xtc_parts))
    return runs, skipped


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
            "frames": 0,
            "mean": None,
            "first0_5": None,
            "mid20_30": None,
            "tail40_50": None,
            "final": None,
            "min": None,
            "max": None,
            "tail_minus_first": None,
        }
    vals = [v for _, v in rows]
    first = [v for t, v in rows if 0 <= t <= 5]
    mid = [v for t, v in rows if 20 <= t <= 30]
    tail = [v for t, v in rows if 40 <= t <= 50]
    if not tail:
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
        run_cmd(cmd, cwd=ROOT, log=out_dir / "logs" / f"{run.variant}_trjcat.log")
    check_out = run_cmd([gmx, "check", "-f", str(concat)], cwd=ROOT, log=out_dir / "logs" / f"{run.variant}_check.log")
    return concat, parse_gmx_check(check_out)


def make_index(
    gmx: str,
    run: Run,
    out_dir: Path,
    a1_selection: str,
    gpiba_selection: str,
    force: bool,
) -> Path:
    ndx = out_dir / "index" / f"{run.variant}_a1_gpiba.ndx"
    ndx.parent.mkdir(parents=True, exist_ok=True)
    if force or not ndx.exists():
        select = f'"VWF_A1" {a1_selection}; "GPIBA" {gpiba_selection}'
        run_cmd(
            [gmx, "select", "-s", str(run.tpr), "-on", str(ndx), "-select", select],
            cwd=ROOT,
            log=out_dir / "logs" / f"{run.variant}_select.log",
        )
    return ndx


def compute_interface(
    gmx: str,
    run: Run,
    traj: Path,
    out_dir: Path,
    ndx: Path,
    cutoff_nm: float,
    force: bool,
    dt_ns: float | None,
) -> tuple[dict[str, float | int | None], list[dict]]:
    prefix = out_dir / "interface" / run.variant
    prefix.parent.mkdir(parents=True, exist_ok=True)
    mindist = prefix.with_name(f"{run.variant}_a1_gpiba_mindist.xvg")
    numcont = prefix.with_name(f"{run.variant}_a1_gpiba_numcont.xvg")
    if force or not (mindist.exists() and numcont.exists()):
        cmd = [
            gmx,
            "mindist",
            "-s",
            str(run.tpr),
            "-f",
            str(traj),
            "-n",
            str(ndx),
            "-od",
            str(mindist),
            "-on",
            str(numcont),
            "-d",
            str(cutoff_nm),
            "-group",
            "-tu",
            "ns",
            "-xvg",
            "none",
        ]
        if dt_ns is not None and dt_ns > 0:
            cmd.extend(["-dt", str(dt_ns)])
        run_cmd(
            cmd,
            cwd=ROOT,
            stdin="VWF_A1\nGPIBA\n",
            log=out_dir / "logs" / f"{run.variant}_mindist.log",
        )

    contact_rows = parse_xvg(numcont)
    dist_rows = parse_xvg(mindist)
    cstats = summarize_series(contact_rows)
    dstats = summarize_series(dist_rows)

    first = cstats["first0_5"]
    tail = cstats["tail40_50"]
    loss_abs = (first - tail) if isinstance(first, (int, float)) and isinstance(tail, (int, float)) else None
    retention = (tail / first) if isinstance(first, (int, float)) and first > 0 and isinstance(tail, (int, float)) else None
    loss_frac = (loss_abs / first) if isinstance(first, (int, float)) and first > 0 and isinstance(loss_abs, (int, float)) else None

    summary = {
        "a1_gpiba_contacts_frames": cstats["frames"],
        "a1_gpiba_contacts_mean": cstats["mean"],
        "a1_gpiba_contacts_first0_5": cstats["first0_5"],
        "a1_gpiba_contacts_mid20_30": cstats["mid20_30"],
        "a1_gpiba_contacts_tail40_50": cstats["tail40_50"],
        "a1_gpiba_contacts_final": cstats["final"],
        "a1_gpiba_contacts_min": cstats["min"],
        "a1_gpiba_contacts_max": cstats["max"],
        "a1_gpiba_contacts_tail_minus_first": cstats["tail_minus_first"],
        "a1_gpiba_contact_loss_abs": round(loss_abs, 4) if loss_abs is not None else None,
        "a1_gpiba_contact_retention": round(retention, 4) if retention is not None else None,
        "a1_gpiba_contact_loss_frac": round(loss_frac, 4) if loss_frac is not None else None,
        "a1_gpiba_mindist_frames": dstats["frames"],
        "a1_gpiba_mindist_mean_nm": dstats["mean"],
        "a1_gpiba_mindist_first0_5_nm": dstats["first0_5"],
        "a1_gpiba_mindist_mid20_30_nm": dstats["mid20_30"],
        "a1_gpiba_mindist_tail40_50_nm": dstats["tail40_50"],
        "a1_gpiba_mindist_final_nm": dstats["final"],
        "a1_gpiba_mindist_min_nm": dstats["min"],
        "a1_gpiba_mindist_max_nm": dstats["max"],
        "a1_gpiba_mindist_tail_minus_first_nm": dstats["tail_minus_first"],
    }

    timeseries = []
    by_time = {t: {"time_ns": t, "contacts": v} for t, v in contact_rows}
    for t, v in dist_rows:
        by_time.setdefault(t, {"time_ns": t})["mindist_nm"] = v
    for t in sorted(by_time):
        row = {"variant": run.variant, "label": run.label, **by_time[t]}
        timeseries.append(row)
    return summary, timeseries


def zscore(values: list[float | None], value: float | None) -> float | None:
    vals = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    if value is None or not vals:
        return None
    mu = mean(vals)
    sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / len(vals))
    if sd == 0:
        return None
    return round((float(value) - mu) / sd, 4)


def add_classifier_features(rows: list[dict], loss_scale_contacts: float) -> list[dict]:
    retention_values = [r.get("a1_gpiba_contact_retention") for r in rows]
    tail_values = [r.get("a1_gpiba_contacts_tail40_50") for r in rows]
    mindist_values = [r.get("a1_gpiba_mindist_tail40_50_nm") for r in rows]

    out = []
    for row in rows:
        loss_abs = row.get("a1_gpiba_contact_loss_abs")
        loss_score = (float(loss_abs) / loss_scale_contacts) if isinstance(loss_abs, (int, float)) else None
        retained_parts = []
        retention_z = zscore(retention_values, row.get("a1_gpiba_contact_retention"))
        tail_z = zscore(tail_values, row.get("a1_gpiba_contacts_tail40_50"))
        # Larger min distance is less retained, so invert its z-score.
        mindist_z = zscore(mindist_values, row.get("a1_gpiba_mindist_tail40_50_nm"))
        if retention_z is not None:
            retained_parts.append(retention_z)
        if tail_z is not None:
            retained_parts.append(tail_z)
        if mindist_z is not None:
            retained_parts.append(-mindist_z)
        retained_z = round(mean(retained_parts), 4) if retained_parts else None
        loss_z = round(-retained_z, 4) if retained_z is not None else None

        rec = {
            "variant": row["variant"],
            "label": row.get("label", "?"),
            "md_gpiba_interface_loss_score": round(loss_score, 4) if loss_score is not None else None,
            "md_gpiba_interface_retained_z": retained_z,
            "md_gpiba_interface_loss_z": loss_z,
            "md_gpiba_contact_retention": row.get("a1_gpiba_contact_retention"),
            "md_gpiba_contact_loss_frac": row.get("a1_gpiba_contact_loss_frac"),
            "md_gpiba_contact_loss_abs": row.get("a1_gpiba_contact_loss_abs"),
            "md_gpiba_contacts_first0_5": row.get("a1_gpiba_contacts_first0_5"),
            "md_gpiba_contacts_tail40_50": row.get("a1_gpiba_contacts_tail40_50"),
            "md_gpiba_mindist_tail40_50_nm": row.get("a1_gpiba_mindist_tail40_50_nm"),
            "md_feature_source": "a1_gpiba_equilibrium_md_interface_retention",
        }
        out.append(rec)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help="A1-GPIbalpha MD root, used if --status-csv is absent")
    parser.add_argument("--status-csv", type=Path, default=DEFAULT_STATUS, help="completed/running summary CSV")
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS, help="optional labels table for readout only")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--gmx", default=str(ROOT / "envs/gromacs/bin.AVX2_256/gmx"))
    parser.add_argument("--variants", default=None, help="comma-separated variant subset")
    parser.add_argument("--a1-selection", default=DEFAULT_A1_SELECTION)
    parser.add_argument("--gpiba-selection", default=DEFAULT_GPIBA_SELECTION)
    parser.add_argument("--contact-cutoff-nm", type=float, default=0.45)
    parser.add_argument("--dt-ns", type=float, default=None, help="optional trajectory sampling interval for gmx mindist, in ns")
    parser.add_argument(
        "--loss-scale-contacts",
        type=float,
        default=20.0,
        help="Contact loss that maps to md_gpiba_interface_loss_score=1.0; provisional, for downstream sweeps.",
    )
    parser.add_argument("--force", action="store_true", help="rebuild intermediate xvg/index files")
    parser.add_argument("--dry-run", action="store_true", help="discover runs and write skipped table, but do not call GROMACS")
    args = parser.parse_args()

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    labels = load_labels(args.labels.resolve())
    variants = parse_variants(args.variants)

    runs, skipped = discover_from_status(args.status_csv.resolve(), labels, variants)
    if not runs and not args.status_csv.exists():
        runs, skipped = discover_from_root(args.root.resolve(), labels, variants)

    write_csv(out_dir / "a1_gpiba_interface_skipped.csv", skipped)
    print(f"Runs ready for A1-GPIbalpha interface analysis: {len(runs)}")
    if skipped:
        print(f"Skipped: {len(skipped)} (see {out_dir / 'a1_gpiba_interface_skipped.csv'})")
    if args.dry_run:
        for run in runs:
            print(f"DRY-RUN {run.variant}: {run.prod_dir}")
        return 0
    if not runs:
        print("No analyzable completed runs found on this machine.")
        return 0

    gmx = str(Path(args.gmx).resolve()) if Path(args.gmx).exists() else args.gmx
    summary_rows: list[dict] = []
    timeseries_rows: list[dict] = []
    for run in runs:
        print(f"Analyzing {run.variant} ({run.label})")
        traj, check = concat_and_check(gmx, run, out_dir, args.force)
        ndx = make_index(gmx, run, out_dir, args.a1_selection, args.gpiba_selection, args.force)
        iface, ts = compute_interface(gmx, run, traj, out_dir, ndx, args.contact_cutoff_nm, args.force, args.dt_ns)
        summary_rows.append({
            "variant": run.variant,
            "label": run.label,
            "prod_dir": str(run.prod_dir),
            "final_gro": str(run.final_gro),
            "concat_xtc": str(traj),
            "xtc_parts": ";".join(str(p) for p in run.xtc_parts),
            "a1_selection": args.a1_selection,
            "gpiba_selection": args.gpiba_selection,
            "contact_cutoff_nm": args.contact_cutoff_nm,
            **check,
            **iface,
        })
        timeseries_rows.extend(ts)

    classifier_rows = add_classifier_features(summary_rows, args.loss_scale_contacts)
    write_csv(out_dir / "a1_gpiba_interface_summary.csv", summary_rows)
    write_csv(out_dir / "a1_gpiba_interface_timeseries.csv", timeseries_rows)
    write_csv(out_dir / "a1_gpiba_classifier_features.csv", classifier_rows)

    print(f"Wrote {out_dir / 'a1_gpiba_interface_summary.csv'}")
    print(f"Wrote {out_dir / 'a1_gpiba_classifier_features.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

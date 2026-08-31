#!/usr/bin/env python3
"""Consolidate completed VWF MD/SMD result layers into the computational panel.

This script copies only lightweight CSV/JSON artifacts. It does not copy raw
trajectories, TPR files, GRO files, or machine-specific absolute paths.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "outputs/computational_panel_20260829/md"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def copy_csv(source: Path, destination: Path) -> tuple[int, int]:
    copy_file(source, destination)
    with destination.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return len(rows), len(reader.fieldnames or [])


def copy_json(source: Path, destination: Path) -> None:
    copy_file(source, destination)


def sanitize_a1_gpiba_summary(source: Path, destination: Path) -> tuple[int, int]:
    keep_columns = [
        "variant",
        "label",
        "contact_cutoff_nm",
        "frames",
        "timestep_ps",
        "atoms",
        "a1_gpiba_contacts_frames",
        "a1_gpiba_contacts_mean",
        "a1_gpiba_contacts_first0_5",
        "a1_gpiba_contacts_mid20_30",
        "a1_gpiba_contacts_tail40_50",
        "a1_gpiba_contacts_final",
        "a1_gpiba_contacts_min",
        "a1_gpiba_contacts_max",
        "a1_gpiba_contacts_tail_minus_first",
        "a1_gpiba_contact_loss_abs",
        "a1_gpiba_contact_retention",
        "a1_gpiba_contact_loss_frac",
        "a1_gpiba_mindist_frames",
        "a1_gpiba_mindist_mean_nm",
        "a1_gpiba_mindist_first0_5_nm",
        "a1_gpiba_mindist_mid20_30_nm",
        "a1_gpiba_mindist_tail40_50_nm",
        "a1_gpiba_mindist_final_nm",
        "a1_gpiba_mindist_min_nm",
        "a1_gpiba_mindist_max_nm",
        "a1_gpiba_mindist_tail_minus_first_nm",
    ]

    with source.open(newline="") as handle:
        reader = csv.DictReader(handle)
        input_columns = reader.fieldnames or []
        missing = [column for column in keep_columns if column not in input_columns]
        if missing:
            raise ValueError(f"Missing expected columns in {source}: {missing}")
        rows = [{column: row[column] for column in keep_columns} for row in reader]

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keep_columns)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows), len(keep_columns)


def write_layer_manifest(
    output_dir: Path,
    layers: list[dict],
) -> None:
    manifest_path = output_dir / "md_result_layers.csv"
    fieldnames = [
        "layer",
        "output_path",
        "source_path",
        "n_rows",
        "n_columns",
        "sha256",
        "status",
        "interpretation",
    ]
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(layers)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination MD result directory.",
    )
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()

    sources = {
        "aim_masking_features": REPO_ROOT / "output/md_7a6o_features.csv",
        "aim_masking_metadata": REPO_ROOT / "output/md_7a6o_masking_interface.json",
        "aim_saltbridge_features": REPO_ROOT / "output/md_7a6o_saltbridge_features.csv",
        "smd_summary": REPO_ROOT / "output/md_7a6o_smd_slow025_features.csv",
        "smd_per_replicate": REPO_ROOT / "output/md_7a6o_smd_slow025_features_perrep.csv",
        "a1_gpiba_classifier_features": REPO_ROOT
        / "output/type2m_lof_md_analysis_2026-06-29/a1_gpiba_interface_qc_dt1/a1_gpiba_classifier_features.csv",
        "a1_gpiba_interface_summary": REPO_ROOT
        / "output/type2m_lof_md_analysis_2026-06-29/a1_gpiba_interface_qc_dt1/a1_gpiba_interface_summary.csv",
        "a1_gpiba_interface_timeseries": REPO_ROOT
        / "output/type2m_lof_md_analysis_2026-06-29/a1_gpiba_interface_qc_dt1/a1_gpiba_interface_timeseries.csv",
    }

    for name, path in sources.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing source file for {name}: {path}")

    layers: list[dict] = []

    def add_layer(
        layer: str,
        destination: Path,
        source: Path,
        n_rows: int | None,
        n_columns: int | None,
        status: str,
        interpretation: str,
    ) -> None:
        layers.append(
            {
                "layer": layer,
                "output_path": str(destination.relative_to(REPO_ROOT)),
                "source_path": str(source.relative_to(REPO_ROOT)),
                "n_rows": n_rows if n_rows is not None else "",
                "n_columns": n_columns if n_columns is not None else "",
                "sha256": sha256(destination),
                "status": status,
                "interpretation": interpretation,
            }
        )

    # A1/AIM masking interface.
    destination = output_dir / "a1_aim_masking/a1_aim_masking_features.csv"
    rows, columns = copy_csv(sources["aim_masking_features"], destination)
    add_layer(
        "A1_AIM_masking_equilibrium_MD",
        destination,
        sources["aim_masking_features"],
        rows,
        columns,
        "complete",
        "AIM-A1 contact retention and binding-face stability features for the 7A6O reference panel",
    )

    destination = output_dir / "a1_aim_masking/a1_aim_masking_interface.json"
    copy_json(sources["aim_masking_metadata"], destination)
    add_layer(
        "A1_AIM_masking_metadata",
        destination,
        sources["aim_masking_metadata"],
        None,
        None,
        "complete",
        "Masking-interface residue definitions, cutoff, and z-score calibration metadata",
    )

    # A1/AIM salt-bridge equilibrium MD.
    destination = output_dir / "a1_aim_saltbridge/a1_aim_saltbridge_features.csv"
    rows, columns = copy_csv(sources["aim_saltbridge_features"], destination)
    add_layer(
        "A1_AIM_saltbridge_equilibrium_MD",
        destination,
        sources["aim_saltbridge_features"],
        rows,
        columns,
        "complete",
        "AIM-A1 salt-bridge occupancy and retained-z features across the completed 7A6O cohort",
    )

    # A1-GPIbα interface equilibrium MD.
    destination = output_dir / "a1_gpiba_interface/a1_gpiba_classifier_features.csv"
    rows, columns = copy_csv(sources["a1_gpiba_classifier_features"], destination)
    add_layer(
        "A1_GPIb_alpha_interface_equilibrium_MD",
        destination,
        sources["a1_gpiba_classifier_features"],
        rows,
        columns,
        "complete",
        "A1-GPIbα interface retention, contact loss, and minimum-distance features",
    )

    destination = output_dir / "a1_gpiba_interface/a1_gpiba_interface_summary.csv"
    rows, columns = sanitize_a1_gpiba_summary(
        sources["a1_gpiba_interface_summary"], destination
    )
    add_layer(
        "A1_GPIb_alpha_interface_summary",
        destination,
        sources["a1_gpiba_interface_summary"],
        rows,
        columns,
        "complete_sanitized",
        "Portable A1-GPIbα interface summary with trajectory and machine-specific paths removed",
    )

    destination = output_dir / "a1_gpiba_interface/a1_gpiba_interface_timeseries.csv"
    rows, columns = copy_csv(sources["a1_gpiba_interface_timeseries"], destination)
    add_layer(
        "A1_GPIb_alpha_interface_timeseries",
        destination,
        sources["a1_gpiba_interface_timeseries"],
        rows,
        columns,
        "complete",
        "Per-frame A1-GPIbα contacts and minimum distance",
    )

    # 7A6O slow025 SMD calibration.
    destination = output_dir / "smd/smd_slow025_features.csv"
    rows, columns = copy_csv(sources["smd_summary"], destination)
    add_layer(
        "A1_AIM_slow025_SMD_summary",
        destination,
        sources["smd_summary"],
        rows,
        columns,
        "complete_no_go",
        "Three-control slow025 SMD calibration; direction is reversed and must not enter the classifier",
    )

    destination = output_dir / "smd/smd_slow025_features_perrep.csv"
    rows, columns = copy_csv(sources["smd_per_replicate"], destination)
    add_layer(
        "A1_AIM_slow025_SMD_per_replicate",
        destination,
        sources["smd_per_replicate"],
        rows,
        columns,
        "complete_no_go",
        "Per-replicate force and work curves used for the slow025 no-go diagnosis",
    )

    write_layer_manifest(output_dir, layers)
    print(f"Wrote {len(layers)} MD/SMD result layers to {output_dir}")


if __name__ == "__main__":
    main()

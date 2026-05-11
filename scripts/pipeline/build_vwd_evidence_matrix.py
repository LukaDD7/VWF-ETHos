#!/usr/bin/env python3
"""
Build the VWD/VWF functional evidence matrix.

This script joins:
  - variants_master.csv
  - diagnostic_panel.csv
  - job_manifest.csv
  - boltz_results_summary.csv

It writes both a long evidence table (one row per variant x assay) and a wide
matrix (one row per variant). This is the central input table for downstream
subtype calibration and the diagnostic Agent.

Important metric rule:
  - complex/interface assays use iPTM as the primary structural metric
  - monomer/context assays use pTM first, then complex pLDDT as a fallback
  - monomer iPTM=0 is expected and is not used as evidence
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PANEL_DIR = PROJECT_ROOT / "output" / "boltz2_vwd_functional_panel"

POSITIVE_LABELS = ("Type1", "Type2A", "Type2B", "Type2M", "Type2N", "Type3")
TYPE2_LABELS = ("Type2A", "Type2B", "Type2M", "Type2N")
NEGATIVE_LABELS = ("Negative_GeneBe", "Negative_ClinVar")


def read_csv(path: Path, required: bool = True) -> pd.DataFrame:
    if not path.exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    return pd.read_csv(path)


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def split_labels(value: object) -> set[str]:
    text = clean_text(value)
    if not text:
        return set()
    return {part for part in text.split("|") if part}


def add_label_policy_columns(variants: pd.DataFrame) -> pd.DataFrame:
    variants = variants.copy()
    label_sets = variants["source_labels"].apply(split_labels)

    for label in POSITIVE_LABELS + NEGATIVE_LABELS:
        variants[f"is_{label}"] = label_sets.apply(lambda labels, lab=label: lab in labels)

    variants["positive_any"] = label_sets.apply(lambda labels: any(label in labels for label in POSITIVE_LABELS))
    variants["type2_any"] = label_sets.apply(lambda labels: any(label in labels for label in TYPE2_LABELS))
    variants["negative_any"] = label_sets.apply(lambda labels: any(label in labels for label in NEGATIVE_LABELS))
    variants["negative_only"] = variants["negative_any"] & ~variants["positive_any"]
    variants["positive_only"] = variants["positive_any"] & ~variants["negative_any"]
    variants["negative_positive_overlap"] = variants["negative_any"] & variants["positive_any"]
    variants["n_source_labels"] = label_sets.apply(len)
    variants["label_policy"] = variants.apply(classify_label_policy, axis=1)
    return variants


def classify_label_policy(row: pd.Series) -> str:
    if row.get("negative_only", False):
        return "negative_only"
    if row.get("negative_positive_overlap", False):
        return "negative_positive_overlap"
    if row.get("positive_only", False):
        return "positive_only"
    return "unlabeled_or_other"


def is_complex_assay(ligand_keys: object) -> bool:
    text = clean_text(ligand_keys)
    return text.lower() not in ("", "nan", "none")


def normalize_summary(summary: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    summary = summary.copy()

    for column in [
        "avg_iptm",
        "best_iptm",
        "avg_ptm",
        "best_ptm",
        "avg_complex_plddt",
        "best_complex_plddt",
    ]:
        if column not in summary.columns:
            summary[column] = pd.NA

    manifest_cols = [
        "job_name",
        "variant_id",
        "assay_key",
        "clinical_axis",
        "ligand_keys",
        "n_chains",
        "run_decision",
    ]
    manifest_for_join = manifest[[c for c in manifest_cols if c in manifest.columns]].drop_duplicates("job_name")
    summary = summary.merge(
        manifest_for_join,
        on="job_name",
        how="left",
        suffixes=("", "_manifest"),
    )

    for column in ("variant_id", "assay_key", "clinical_axis"):
        fallback = f"{column}_manifest"
        if fallback in summary.columns:
            summary[column] = summary[column].where(summary[column].notna(), summary[fallback])
            summary = summary.drop(columns=[fallback])

    if "primary_metric" not in summary.columns:
        summary["primary_metric"] = summary["ligand_keys"].apply(
            lambda ligands: "iptm" if is_complex_assay(ligands) else "ptm_or_plddt"
        )
    if "primary_metric_reason" not in summary.columns:
        summary["primary_metric_reason"] = summary["primary_metric"].map(
            {
                "iptm": "complex/interface assay",
                "ptm_or_plddt": "monomer assay; iPTM is not meaningful without an interface",
            }
        )

    summary["primary_value"] = summary.apply(select_primary_value, axis=1)
    summary["primary_value_source"] = summary.apply(select_primary_value_source, axis=1)
    summary["metric_status"] = summary.apply(metric_status, axis=1)
    return summary


def select_primary_value(row: pd.Series):
    if row.get("primary_metric") == "iptm":
        return row.get("avg_iptm")
    for column in ("avg_ptm", "avg_complex_plddt"):
        value = row.get(column)
        if pd.notna(value):
            return value
    return pd.NA


def select_primary_value_source(row: pd.Series) -> str:
    if row.get("primary_metric") == "iptm":
        return "avg_iptm"
    if pd.notna(row.get("avg_ptm")):
        return "avg_ptm"
    if pd.notna(row.get("avg_complex_plddt")):
        return "avg_complex_plddt"
    return "missing"


def metric_status(row: pd.Series) -> str:
    if pd.notna(row.get("primary_value")):
        return "ok"
    if row.get("primary_metric") == "ptm_or_plddt":
        return "missing_monomer_ptm_plddt_rerun_parser"
    return "missing_primary_metric"


def add_wt_delta_features(long_df: pd.DataFrame) -> pd.DataFrame:
    long_df = long_df.copy()
    for column in ("primary_value", "avg_iptm", "avg_ptm", "avg_complex_plddt"):
        if column in long_df.columns:
            long_df[column] = pd.to_numeric(long_df[column], errors="coerce")
    wt = (
        long_df[long_df["variant_id"] == "VWF_WT"][
            [
                "assay_key",
                "primary_value",
                "avg_iptm",
                "avg_ptm",
                "avg_complex_plddt",
            ]
        ]
        .rename(
            columns={
                "primary_value": "wt_primary_value",
                "avg_iptm": "wt_avg_iptm",
                "avg_ptm": "wt_avg_ptm",
                "avg_complex_plddt": "wt_avg_complex_plddt",
            }
        )
        .drop_duplicates("assay_key")
    )
    long_df = long_df.merge(wt, on="assay_key", how="left")

    for metric in ("primary_value", "avg_iptm", "avg_ptm", "avg_complex_plddt"):
        wt_col = f"wt_{metric}" if metric != "primary_value" else "wt_primary_value"
        delta_col = f"{metric}_delta_vs_wt"
        if wt_col in long_df.columns:
            long_df[delta_col] = long_df[metric] - long_df[wt_col]

    non_wt = long_df["variant_id"] != "VWF_WT"
    long_df["primary_delta_direction"] = "not_available"
    long_df.loc[non_wt & (long_df["primary_value_delta_vs_wt"] > 0), "primary_delta_direction"] = "higher_than_wt"
    long_df.loc[non_wt & (long_df["primary_value_delta_vs_wt"] < 0), "primary_delta_direction"] = "lower_than_wt"
    long_df.loc[non_wt & (long_df["primary_value_delta_vs_wt"] == 0), "primary_delta_direction"] = "same_as_wt"
    return long_df


def add_distribution_features(long_df: pd.DataFrame) -> pd.DataFrame:
    long_df = long_df.copy()
    long_df["primary_value"] = pd.to_numeric(long_df["primary_value"], errors="coerce")
    non_wt_mask = long_df["variant_id"] != "VWF_WT"
    values = long_df.loc[non_wt_mask].copy()

    assay_stats = values.groupby("assay_key")["primary_value"].agg(["mean", "std"]).rename(
        columns={"mean": "assay_primary_mean", "std": "assay_primary_std"}
    )
    long_df = long_df.merge(assay_stats, on="assay_key", how="left")
    long_df["primary_zscore_within_assay"] = (
        (long_df["primary_value"] - long_df["assay_primary_mean"]) / long_df["assay_primary_std"]
    )
    long_df.loc[long_df["assay_primary_std"].fillna(0) == 0, "primary_zscore_within_assay"] = pd.NA

    long_df["primary_percentile_within_assay"] = pd.NA
    for assay_key, idx in long_df.loc[non_wt_mask].groupby("assay_key").groups.items():
        ranks = long_df.loc[idx, "primary_value"].rank(pct=True, method="average")
        long_df.loc[idx, "primary_percentile_within_assay"] = ranks

    domain_stats = values.groupby(["assay_key", "inferred_domain"])["primary_value"].agg(["mean", "std"]).rename(
        columns={"mean": "domain_primary_mean", "std": "domain_primary_std"}
    )
    long_df = long_df.merge(domain_stats, on=["assay_key", "inferred_domain"], how="left")
    long_df["primary_zscore_within_assay_domain"] = (
        (long_df["primary_value"] - long_df["domain_primary_mean"]) / long_df["domain_primary_std"]
    )
    long_df.loc[long_df["domain_primary_std"].fillna(0) == 0, "primary_zscore_within_assay_domain"] = pd.NA
    return long_df


def sanitize_column(value: str) -> str:
    value = re.sub(r"[^0-9A-Za-z]+", "_", str(value)).strip("_").lower()
    return value or "unknown"


def make_wide_matrix(long_df: pd.DataFrame, variants: pd.DataFrame) -> pd.DataFrame:
    assay_feature_cols = [
        "run_decision",
        "n_samples",
        "primary_metric",
        "primary_value",
        "primary_value_source",
        "metric_status",
        "wt_primary_value",
        "primary_value_delta_vs_wt",
        "primary_delta_direction",
        "primary_zscore_within_assay",
        "primary_percentile_within_assay",
        "avg_iptm",
        "best_iptm",
        "avg_ptm",
        "best_ptm",
        "avg_complex_plddt",
        "best_complex_plddt",
    ]

    rows = []
    for variant_id, sub in long_df[long_df["variant_id"] != "VWF_WT"].groupby("variant_id", sort=True):
        row = {"variant_id": variant_id}
        for _, assay_row in sub.iterrows():
            prefix = sanitize_column(assay_row["assay_key"])
            for col in assay_feature_cols:
                row[f"{prefix}__{col}"] = assay_row.get(col)
        rows.append(row)

    wide = pd.DataFrame(rows)
    metadata_cols = [
        c for c in variants.columns
        if c not in wide.columns or c == "variant_id"
    ]
    wide = variants[metadata_cols].merge(wide, on="variant_id", how="left")
    return wide


def build_evidence_matrix(panel_dir: Path, output_long: Path, output_wide: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    variants = read_csv(panel_dir / "variants_master.csv")
    diagnostic = read_csv(panel_dir / "diagnostic_panel.csv")
    manifest = read_csv(panel_dir / "job_manifest.csv")
    summary = read_csv(panel_dir / "boltz_results_summary.csv")

    variants = add_label_policy_columns(variants)
    summary = normalize_summary(summary, manifest)

    diagnostic_keys = [
        "variant_id",
        "assay_key",
        "aa_change",
        "position",
        "source_labels",
        "inferred_domain",
        "clinical_axis",
        "construct_role",
        "vwf_range",
        "ligand_keys",
        "run_decision",
        "expected_signal",
        "interpretability",
        "notes",
    ]
    diagnostic = diagnostic[[c for c in diagnostic_keys if c in diagnostic.columns]].copy()
    wt_rows = manifest[manifest["variant_id"] == "VWF_WT"][
        [
            "variant_id",
            "assay_key",
            "aa_change",
            "position",
            "source_labels",
            "inferred_domain",
            "clinical_axis",
            "construct_role",
            "vwf_range",
            "ligand_keys",
            "run_decision",
            "expected_signal",
            "interpretability",
            "notes",
        ]
    ].copy()
    assay_panel = pd.concat([diagnostic, wt_rows], ignore_index=True)

    summary_cols = [
        "job_name",
        "variant_id",
        "assay_key",
        "n_samples",
        "avg_iptm",
        "best_iptm",
        "avg_ptm",
        "best_ptm",
        "avg_complex_plddt",
        "best_complex_plddt",
        "primary_metric",
        "primary_metric_reason",
        "primary_value",
        "primary_value_source",
        "metric_status",
    ]
    long_df = assay_panel.merge(
        summary[[c for c in summary_cols if c in summary.columns]],
        on=["variant_id", "assay_key"],
        how="left",
    )

    long_df["is_assay_run"] = long_df["run_decision"].isin(["RUN", "WT_BASELINE"])
    long_df["has_boltz_result"] = long_df["n_samples"].notna()
    long_df = long_df.merge(
        variants[
            [
                "variant_id",
                "wt_aa",
                "mut_aa",
                "aa_change_3letter",
                "source_files",
                "source_domains",
                "wt_validation_status",
                "label_policy",
                "positive_any",
                "type2_any",
                "negative_any",
                "negative_only",
                "positive_only",
                "negative_positive_overlap",
                *[f"is_{label}" for label in POSITIVE_LABELS + NEGATIVE_LABELS],
            ]
        ],
        on="variant_id",
        how="left",
    )

    long_df = add_wt_delta_features(long_df)
    long_df = add_distribution_features(long_df)

    output_long.parent.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(output_long, index=False)

    wide = make_wide_matrix(long_df, variants)
    output_wide.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(output_wide, index=False)
    return long_df, wide


def print_summary(long_df: pd.DataFrame, wide: pd.DataFrame, output_long: Path, output_wide: Path) -> None:
    non_wt = long_df[long_df["variant_id"] != "VWF_WT"]
    print("=" * 72)
    print("VWD/VWF Evidence Matrix")
    print("=" * 72)
    print(f"Variants             : {len(wide)}")
    print(f"Long rows            : {len(non_wt)} variant-assay rows")
    print(f"Assays               : {non_wt['assay_key'].nunique()}")
    print(f"Rows with Boltz data : {int(non_wt['has_boltz_result'].sum())}")
    print(f"Rows with primary    : {int(non_wt['primary_value'].notna().sum())}")
    print()
    print("Metric status:")
    print(non_wt["metric_status"].fillna("no_result").value_counts().to_string())
    print()
    print("Label policy:")
    label_cols = ["variant_id", "label_policy"]
    label_policy = wide[label_cols].drop_duplicates()["label_policy"].value_counts()
    print(label_policy.to_string())
    print()
    print(f"Long evidence : {output_long}")
    print(f"Wide matrix   : {output_wide}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", type=Path, default=DEFAULT_PANEL_DIR)
    parser.add_argument(
        "--output-long",
        type=Path,
        default=DEFAULT_PANEL_DIR / "evidence_long.csv",
    )
    parser.add_argument(
        "--output-wide",
        type=Path,
        default=DEFAULT_PANEL_DIR / "evidence_matrix.csv",
    )
    args = parser.parse_args()

    panel_dir = args.panel_dir if args.panel_dir.is_absolute() else PROJECT_ROOT / args.panel_dir
    output_long = args.output_long if args.output_long.is_absolute() else PROJECT_ROOT / args.output_long
    output_wide = args.output_wide if args.output_wide.is_absolute() else PROJECT_ROOT / args.output_wide

    long_df, wide = build_evidence_matrix(panel_dir, output_long, output_wide)
    print_summary(long_df, wide, output_long, output_wide)


if __name__ == "__main__":
    main()

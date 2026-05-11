#!/usr/bin/env python3
"""
Generate first-pass analysis tables from the VWD evidence matrix.

Inputs:
  - output/boltz2_vwd_functional_panel/evidence_long.csv
  - output/boltz2_vwd_functional_panel/evidence_matrix.csv

Outputs:
  - analysis/assay_metric_summary.csv
  - analysis/label_distribution.csv
  - analysis/negative_vs_positive_by_assay.csv
  - analysis/subtype_one_vs_rest_by_assay.csv
  - analysis/top_primary_delta_outliers.csv
  - analysis/README.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PANEL_DIR = PROJECT_ROOT / "output" / "boltz2_vwd_functional_panel"
POSITIVE_LABELS = ("Type1", "Type2A", "Type2B", "Type2M", "Type2N", "Type3")
A1_ASSAYS = (
    "a1_aim_autoinhibition_context",
    "a1_gpiba_forced_binding",
    "a1_heparan_sulfate_binding",
)


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, low_memory=False)


def cohen_d(a: pd.Series, b: pd.Series) -> float | None:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return None
    var_a = a.var(ddof=1)
    var_b = b.var(ddof=1)
    pooled = (((len(a) - 1) * var_a + (len(b) - 1) * var_b) / (len(a) + len(b) - 2)) ** 0.5
    if pooled == 0:
        return None
    return (a.mean() - b.mean()) / pooled


def bool_mask(series: pd.Series) -> pd.Series:
    return series.map(lambda value: bool(value) if pd.notna(value) else False)


def assay_metric_summary(long_df: pd.DataFrame) -> pd.DataFrame:
    non_wt = long_df[long_df["variant_id"] != "VWF_WT"].copy()
    grouped = non_wt.groupby("assay_key", dropna=False)
    rows = []
    for assay_key, sub in grouped:
        primary = pd.to_numeric(sub["primary_value"], errors="coerce")
        delta = pd.to_numeric(sub["primary_value_delta_vs_wt"], errors="coerce")
        rows.append(
            {
                "assay_key": assay_key,
                "clinical_axis": first_non_null(sub, "clinical_axis"),
                "primary_metric": first_non_null(sub, "primary_metric"),
                "n_rows": len(sub),
                "n_run": int(sub["is_assay_run"].fillna(False).sum()),
                "n_not_run": int((~bool_mask(sub["is_assay_run"])).sum()),
                "n_with_boltz_result": int(sub["has_boltz_result"].fillna(False).sum()),
                "n_with_primary_value": int(primary.notna().sum()),
                "mean_primary": primary.mean(),
                "std_primary": primary.std(),
                "mean_delta_vs_wt": delta.mean(),
                "std_delta_vs_wt": delta.std(),
                "min_delta_vs_wt": delta.min(),
                "max_delta_vs_wt": delta.max(),
                "n_metric_ok": int((sub["metric_status"] == "ok").sum()),
                "n_run_missing_metric": int(
                    (bool_mask(sub["is_assay_run"]) & (sub["metric_status"].fillna("no_result") != "ok")).sum()
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(["assay_key"]).reset_index(drop=True)


def first_non_null(df: pd.DataFrame, column: str):
    if column not in df.columns:
        return None
    values = df[column].dropna()
    if values.empty:
        return None
    return values.iloc[0]


def label_distribution(matrix: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.append({"label": "all_variants", "n": len(matrix)})
    for col in [
        "positive_any",
        "type2_any",
        "negative_any",
        "negative_only",
        "positive_only",
        "negative_positive_overlap",
    ]:
        if col in matrix.columns:
            rows.append({"label": col, "n": int(matrix[col].fillna(False).sum())})
    for label in POSITIVE_LABELS:
        col = f"is_{label}"
        if col in matrix.columns:
            rows.append({"label": label, "n": int(matrix[col].fillna(False).sum())})
    return pd.DataFrame(rows)


def negative_vs_positive_by_assay(long_df: pd.DataFrame) -> pd.DataFrame:
    non_wt = long_df[long_df["variant_id"] != "VWF_WT"].copy()
    rows = []
    for assay_key, sub in non_wt.groupby("assay_key", dropna=False):
        negative_mask = bool_mask(sub["negative_only"])
        positive_mask = bool_mask(sub["positive_only"])
        neg = sub[negative_mask]
        pos = sub[positive_mask]
        neg_values = pd.to_numeric(neg["primary_value_delta_vs_wt"], errors="coerce")
        pos_values = pd.to_numeric(pos["primary_value_delta_vs_wt"], errors="coerce")
        rows.append(
            {
                "assay_key": assay_key,
                "clinical_axis": first_non_null(sub, "clinical_axis"),
                "primary_metric": first_non_null(sub, "primary_metric"),
                "n_negative_only": int(neg_values.notna().sum()),
                "n_positive_only": int(pos_values.notna().sum()),
                "negative_mean_delta": neg_values.mean(),
                "positive_mean_delta": pos_values.mean(),
                "positive_minus_negative_delta": pos_values.mean() - neg_values.mean(),
                "cohen_d_positive_vs_negative": cohen_d(pos_values, neg_values),
            }
        )
    return pd.DataFrame(rows).sort_values("assay_key").reset_index(drop=True)


def subtype_one_vs_rest_by_assay(long_df: pd.DataFrame) -> pd.DataFrame:
    non_wt = long_df[long_df["variant_id"] != "VWF_WT"].copy()
    rows = []
    for label in POSITIVE_LABELS:
        label_col = f"is_{label}"
        if label_col not in non_wt.columns:
            continue
        for assay_key, sub in non_wt.groupby("assay_key", dropna=False):
            target_mask = bool_mask(sub[label_col])
            overlap_mask = bool_mask(sub["negative_positive_overlap"])
            target = sub[target_mask]
            rest = sub[(~target_mask) & (~overlap_mask)]
            target_values = pd.to_numeric(target["primary_value_delta_vs_wt"], errors="coerce")
            rest_values = pd.to_numeric(rest["primary_value_delta_vs_wt"], errors="coerce")
            rows.append(
                {
                    "target_label": label,
                    "assay_key": assay_key,
                    "clinical_axis": first_non_null(sub, "clinical_axis"),
                    "primary_metric": first_non_null(sub, "primary_metric"),
                    "n_target": int(target_values.notna().sum()),
                    "n_rest": int(rest_values.notna().sum()),
                    "target_mean_delta": target_values.mean(),
                    "rest_mean_delta": rest_values.mean(),
                    "target_minus_rest_delta": target_values.mean() - rest_values.mean(),
                    "cohen_d_target_vs_rest": cohen_d(target_values, rest_values),
                }
            )
    out = pd.DataFrame(rows)
    return out.sort_values(["target_label", "assay_key"]).reset_index(drop=True)


def type2b_vs_type2m_a1_assays(long_df: pd.DataFrame) -> pd.DataFrame:
    non_wt = long_df[
        (long_df["variant_id"] != "VWF_WT")
        & (long_df["assay_key"].isin(A1_ASSAYS))
    ].copy()

    for label in ("Type2A", "Type2B", "Type2M", "Type2N"):
        col = f"is_{label}"
        if col not in non_wt.columns:
            non_wt[col] = False
        non_wt[col] = bool_mask(non_wt[col])

    non_wt["type2b_only"] = (
        non_wt["is_Type2B"]
        & ~non_wt["is_Type2M"]
        & ~non_wt["is_Type2A"]
        & ~non_wt["is_Type2N"]
    )
    non_wt["type2m_only"] = (
        non_wt["is_Type2M"]
        & ~non_wt["is_Type2B"]
        & ~non_wt["is_Type2A"]
        & ~non_wt["is_Type2N"]
    )

    rows = []
    for assay_key, sub in non_wt.groupby("assay_key", dropna=False):
        type2b = pd.to_numeric(
            sub.loc[sub["type2b_only"], "primary_value_delta_vs_wt"],
            errors="coerce",
        ).dropna()
        type2m = pd.to_numeric(
            sub.loc[sub["type2m_only"], "primary_value_delta_vs_wt"],
            errors="coerce",
        ).dropna()
        rows.append(
            {
                "assay_key": assay_key,
                "clinical_axis": first_non_null(sub, "clinical_axis"),
                "primary_metric": first_non_null(sub, "primary_metric"),
                "n_type2b_only": len(type2b),
                "n_type2m_only": len(type2m),
                "type2b_mean_delta": type2b.mean(),
                "type2m_mean_delta": type2m.mean(),
                "type2b_minus_type2m_delta": type2b.mean() - type2m.mean(),
                "cohen_d_type2b_vs_type2m": cohen_d(type2b, type2m),
            }
        )
    return pd.DataFrame(rows).sort_values("assay_key").reset_index(drop=True)


def top_primary_delta_outliers(long_df: pd.DataFrame, n_per_assay: int = 20) -> pd.DataFrame:
    non_wt = long_df[long_df["variant_id"] != "VWF_WT"].copy()
    non_wt["abs_primary_zscore_within_assay"] = pd.to_numeric(
        non_wt["primary_zscore_within_assay"], errors="coerce"
    ).abs()
    out = (
        non_wt.sort_values(["assay_key", "abs_primary_zscore_within_assay"], ascending=[True, False])
        .groupby("assay_key", as_index=False)
        .head(n_per_assay)
    )
    cols = [
        "assay_key",
        "variant_id",
        "aa_change",
        "source_labels",
        "label_policy",
        "inferred_domain",
        "primary_metric",
        "primary_value",
        "wt_primary_value",
        "primary_value_delta_vs_wt",
        "primary_zscore_within_assay",
        "primary_percentile_within_assay",
        "metric_status",
    ]
    return out[[c for c in cols if c in out.columns]].reset_index(drop=True)


def write_readme(output_dir: Path) -> None:
    readme = """# VWD Functional Panel Analysis

Generated by `scripts/pipeline/analyze_vwd_functional_panel.py`.

Files:

- `assay_metric_summary.csv`: per-assay metric coverage and distributions.
- `label_distribution.csv`: label-policy counts.
- `negative_vs_positive_by_assay.csv`: negative-only vs positive-only delta comparison.
- `subtype_one_vs_rest_by_assay.csv`: weakly supervised subtype-vs-rest assay effects.
- `type2b_vs_type2m_a1_assays.csv`: direct A1-axis Type 2B-only vs Type 2M-only comparison.
- `top_primary_delta_outliers.csv`: strongest per-assay structural outliers.

Interpretation:

- Complex/interface assays use iPTM as the primary value.
- Monomer/context assays use pTM first, then complex pLDDT if pTM is unavailable.
- `n_not_run` means a variant was intentionally outside that assay/domain, not
  a failed Boltz job.
"""
    (output_dir / "README.md").write_text(readme)


def run_analysis(panel_dir: Path, output_dir: Path) -> None:
    long_df = read_csv(panel_dir / "evidence_long.csv")
    matrix = read_csv(panel_dir / "evidence_matrix.csv")
    output_dir.mkdir(parents=True, exist_ok=True)

    tables = {
        "assay_metric_summary.csv": assay_metric_summary(long_df),
        "label_distribution.csv": label_distribution(matrix),
        "negative_vs_positive_by_assay.csv": negative_vs_positive_by_assay(long_df),
        "subtype_one_vs_rest_by_assay.csv": subtype_one_vs_rest_by_assay(long_df),
        "type2b_vs_type2m_a1_assays.csv": type2b_vs_type2m_a1_assays(long_df),
        "top_primary_delta_outliers.csv": top_primary_delta_outliers(long_df),
    }
    for name, df in tables.items():
        df.to_csv(output_dir / name, index=False)

    write_readme(output_dir)

    print("=" * 72)
    print("VWD Functional Panel Analysis")
    print("=" * 72)
    print(f"Panel dir : {panel_dir}")
    print(f"Output dir: {output_dir}")
    for name, df in tables.items():
        print(f"  {name:<40} {df.shape[0]:>6} rows")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-dir", type=Path, default=DEFAULT_PANEL_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_PANEL_DIR / "analysis")
    args = parser.parse_args()

    panel_dir = args.panel_dir if args.panel_dir.is_absolute() else PROJECT_ROOT / args.panel_dir
    output_dir = args.output_dir if args.output_dir.is_absolute() else PROJECT_ROOT / args.output_dir
    run_analysis(panel_dir, output_dir)


if __name__ == "__main__":
    main()

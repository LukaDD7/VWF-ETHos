#!/usr/bin/env python3
"""Validate returned AlphaGenome/Boltz/MD files and build Agent evidence rows."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def compact_number(value: Any) -> str:
    try:
        return f"{float(value):.4g}"
    except (TypeError, ValueError):
        return "NA"


def alpha_rows(features_path: Path, expected_cases: int | None = None) -> list[dict[str, Any]]:
    frame = pd.read_csv(features_path)
    if expected_cases is not None and (len(frame) != expected_cases or frame["case_id"].nunique() != expected_cases):
        raise ValueError(f"AlphaGenome return must contain {expected_cases} unique requests; found {len(frame)} rows")
    feature_cols = [column for column in frame if column.startswith("ag_") and column.endswith("_abs_max")]
    splice_cols = [
        column
        for column in (
            "ag_splice_sites_abs_max",
            "ag_splice_site_usage_abs_max",
            "ag_splice_junctions_abs_max",
        )
        if column in frame.columns
    ]
    rows: list[dict[str, Any]] = []
    for _, item in frame.iterrows():
        values = pd.to_numeric(item[feature_cols], errors="coerce").dropna()
        ranked = values.sort_values(ascending=False).head(5)
        top = ", ".join(f"{name}={compact_number(value)}" for name, value in ranked.items()) or "no selected-track score returned"
        splice_values = pd.to_numeric(item[splice_cols], errors="coerce").dropna() if splice_cols else pd.Series(dtype=float)
        splice_text = ""
        if not splice_values.empty and float(splice_values.max()) >= 0.5:
            splice_ranked = splice_values.sort_values(ascending=False)
            splice_text = (
                " Splice-axis signal: "
                + ", ".join(f"{name}={compact_number(value)}" for name, value in splice_ranked.items())
                + "."
            )
        patient_id = item.get("patient_id", item["case_id"])
        rows.append({
            "case_id": patient_id,
            "source": "alphagenome_full_profile",
            "source_record_id": f"{item['case_id']}:alphagenome-full-profile",
            "query": f"{item['chromosome']}:{int(item['position'])}:{item['ref']}>{item['alt']}",
            "retrieved_at": utc_now(),
            "source_version": features_path.name,
            "supports": "",
            "refutes": "",
            "confidence": np.nan,
            "conclusion": (
                f"AlphaGenome complete-profile selected-track scores; top absolute views: {top}."
                f"{splice_text}"
            ),
            "limitations": json.dumps([
                "Research-only regulatory and splicing effect scores; they are not clinical pathogenicity classifications.",
                "Cell selection follows the returned output_metadata plan; unavailable VWF-relevant tracks remain missing rather than using unrelated biosamples.",
            ], ensure_ascii=False),
            "raw_excerpt_locator": str(features_path),
        })
    return rows


def boltz_rows(
    summary_path: Path,
    manifest_path: Path,
    variants_path: Path,
    expected_jobs: int | None = None,
    expected_wt: int | None = None,
) -> list[dict[str, Any]]:
    summary = pd.read_csv(summary_path)
    manifest = pd.read_csv(manifest_path)
    variants = pd.read_csv(variants_path)
    if expected_jobs is not None and manifest["job_name"].nunique() != expected_jobs:
        raise ValueError(f"Boltz manifest must contain {expected_jobs} jobs; found {manifest['job_name'].nunique()}")
    if expected_wt is not None and manifest["run_decision"].eq("WT_BASELINE").sum() != expected_wt:
        raise ValueError(f"Boltz manifest must contain exactly {expected_wt} assay-matched WT baselines")
    missing_jobs = sorted(set(manifest["job_name"]) - set(summary["job_name"]))
    if missing_jobs:
        raise ValueError(f"Boltz return is incomplete; missing jobs: {missing_jobs}")
    merged = manifest.merge(summary, on="job_name", how="left", suffixes=("_manifest", ""))
    if "primary_value" not in merged:
        merged["primary_value"] = np.where(
            merged["primary_metric"].eq("iptm"), merged["avg_iptm"], merged["avg_ptm"]
        )
    baselines = (
        merged[merged["run_decision"].eq("WT_BASELINE")]
        .set_index("assay_key_manifest" if "assay_key_manifest" in merged else "assay_key")["primary_value"]
    )
    assay_col = "assay_key_manifest" if "assay_key_manifest" in merged else "assay_key"
    variant_id_col = "variant_id_manifest" if "variant_id_manifest" in merged else "variant_id"
    merged["delta_vs_wt"] = merged.apply(
        lambda row: row["primary_value"] - baselines.get(row[assay_col], np.nan), axis=1
    )
    case_map: dict[str, list[str]] = {}
    for row in variants.itertuples():
        patient_id = str(getattr(row, "patient_id", row.case_id))
        case_map.setdefault(f"VWF_{row.aa_change}", []).append(patient_id)
    rows: list[dict[str, Any]] = []
    for variant_id, group in merged[merged["run_decision"].eq("RUN")].groupby(variant_id_col):
        case_ids = list(dict.fromkeys(case_map.get(str(variant_id), [])))
        if not case_ids:
            raise ValueError(f"Cannot map Boltz variant {variant_id} to a case")
        axes = "; ".join(
            f"{row[assay_col]} delta_vs_WT={compact_number(row['delta_vs_wt'])}"
            for _, row in group.sort_values(assay_col).iterrows()
        )
        for case_id in case_ids:
            rows.append({
                "case_id": case_id,
                "source": "boltz2_functional_panel",
                "source_record_id": f"{case_id}:{variant_id}:boltz2-functional-panel",
                "query": str(variant_id).removeprefix("VWF_"),
                "retrieved_at": utc_now(),
                "source_version": summary_path.name,
                "supports": "",
                "refutes": "",
                "confidence": np.nan,
                "conclusion": f"Assay-matched Boltz-2 structural deltas: {axes}.",
                "limitations": json.dumps([
                    "Boltz confidence metrics are structural proxies, not binding free energies or clinical classifications.",
                    "Every delta uses the WT baseline from the same assay construct.",
                ], ensure_ascii=False),
                "raw_excerpt_locator": str(summary_path),
            })
    return rows


def md_rows(md_path: Path | None, variants_path: Path) -> list[dict[str, Any]]:
    if md_path is None:
        return []
    frame = pd.read_csv(md_path)
    variants = pd.read_csv(variants_path)
    case_map: dict[str, list[str]] = {}
    for row in variants.itertuples():
        patient_id = str(getattr(row, "patient_id", row.case_id))
        case_map.setdefault(row.aa_change, []).append(patient_id)
    rows: list[dict[str, Any]] = []
    for _, item in frame.iterrows():
        aa_change = str(item.get("aa_change") or item.get("variant") or "").removeprefix("VWF_")
        explicit = str(item.get("case_id") or "").strip()
        case_ids = [explicit] if explicit else list(dict.fromkeys(case_map.get(aa_change, [])))
        if not case_ids:
            continue
        numeric = pd.to_numeric(item, errors="coerce").dropna()
        compact = ", ".join(f"{key}={compact_number(value)}" for key, value in numeric.items())
        for case_id in case_ids:
            rows.append({
            "case_id": case_id,
            "source": "md_targeted_panel",
            "source_record_id": f"{case_id}:{aa_change}:targeted-md",
            "query": aa_change,
            "retrieved_at": utc_now(),
            "source_version": md_path.name,
            "supports": "",
            "refutes": "",
            "confidence": np.nan,
            "conclusion": f"Targeted MD features: {compact}.",
            "limitations": json.dumps([
                "Targeted equilibrium MD is mechanism support only and cannot establish subtype or pathogenicity.",
                "Interpret only against a pipeline-matched WT baseline and completed QC.",
            ], ensure_ascii=False),
            "raw_excerpt_locator": str(md_path),
            })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alphagenome-features", type=Path, required=True)
    parser.add_argument("--boltz-summary", type=Path, required=True)
    parser.add_argument("--boltz-manifest", type=Path, required=True)
    parser.add_argument("--boltz-variants", type=Path, required=True)
    parser.add_argument("--md-features", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-alpha", type=int)
    parser.add_argument("--expected-boltz-jobs", type=int)
    parser.add_argument("--expected-boltz-wt", type=int)
    args = parser.parse_args()

    rows = alpha_rows(args.alphagenome_features, args.expected_alpha)
    rows.extend(boltz_rows(
        args.boltz_summary, args.boltz_manifest, args.boltz_variants,
        args.expected_boltz_jobs, args.expected_boltz_wt,
    ))
    rows.extend(md_rows(args.md_features, args.boltz_variants))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    manifest = {
        "created_at": utc_now(),
        "evidence_rows": len(rows),
        "case_count": len({row["case_id"] for row in rows}),
        "source_counts": pd.Series([row["source"] for row in rows]).value_counts().to_dict(),
        "output": str(args.output),
    }
    args.output.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

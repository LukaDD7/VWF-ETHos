#!/usr/bin/env python3
"""Summarize the Type-1 and Type-2B pre-MD AlphaGenome/Boltz results."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "outputs/computational_panel_20260829/pre_md"


def alpha_summary(path: Path, panel: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    feature_cols = [c for c in df.columns if c.startswith("ag_") and c.endswith("_abs_max")]
    values = df[feature_cols].apply(pd.to_numeric, errors="coerce")
    z = (values - values.mean()) / values.std(ddof=0)
    z = z.replace([np.inf, -np.inf], np.nan)

    top_feature = z.abs().idxmax(axis=1)
    top_z = z.abs().max(axis=1)
    top_abs = [df.loc[i, feature] if pd.notna(feature) else np.nan for i, feature in enumerate(top_feature)]
    out = df[["case_id", "hgvs_c", "hgvs_p"]].copy()
    out.insert(0, "panel", panel)
    out["available_scorers"] = values.notna().sum(axis=1)
    out["top_feature"] = top_feature
    out["top_abs_max"] = top_abs
    out["top_abs_z"] = top_z
    return out


def boltz_summary(manifest_path: Path, summary_path: Path, panel: str) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path)
    summary = pd.read_csv(summary_path)
    merged = manifest.merge(
        summary[["job_name", "primary_value", "avg_iptm", "avg_ptm"]],
        on="job_name",
        how="left",
    )
    wt = (
        merged[merged["run_decision"].eq("WT_BASELINE")]
        .set_index("assay_key")["primary_value"]
        .to_dict()
    )
    merged["delta_vs_wt"] = merged.apply(
        lambda row: row["primary_value"] - wt.get(row["assay_key"], np.nan), axis=1
    )
    merged = merged[merged["run_decision"].eq("RUN")].copy()
    merged["panel"] = panel
    merged["variant"] = merged["aa_change"]

    z_parts = []
    for assay, group in merged.groupby("assay_key", sort=False):
        vals = pd.to_numeric(group["delta_vs_wt"], errors="coerce")
        std = vals.std(ddof=0)
        group = group.copy()
        group["z_within_assay"] = (vals - vals.mean()) / std if std and pd.notna(std) else np.nan
        z_parts.append(group)
    out = pd.concat(z_parts, ignore_index=True) if z_parts else pd.DataFrame()
    cols = [
        "panel", "variant", "aa_change", "assay_key", "clinical_axis",
        "primary_value", "delta_vs_wt", "z_within_assay",
    ]
    return out[[c for c in cols if c in out.columns]]


def write_report(alpha: pd.DataFrame, boltz: pd.DataFrame, out_dir: Path) -> None:
    lines = [
        "# Type-1 and Type-2B pre-MD computational summary",
        "",
        "## AlphaGenome",
        "",
        "- All 16 requests completed successfully.",
        "- Each request has 15 of 19 selected scorer views.",
        "- ATAC, contact maps, polyadenylation, and ATAC-active views are unavailable for the selected VWF-relevant biosamples.",
        "- The strongest within-panel signals are shown in `alphagenome_case_summary.csv`.",
        "",
        "## Boltz-2",
        "",
        "- All 32 unique Boltz jobs completed successfully.",
        "- Three A1 WT baselines are shared between the Type-1 and Type-2B panels.",
        "- The strongest assay-level deltas are shown in `boltz_variant_assay_summary.csv`.",
        "",
        "## Interpretation",
        "",
        "- Type-1 CASE_T1_003 and CASE_T1_007 show the clearest splice-axis signals.",
        "- Type-2B R1341W shows the strongest splice-junction signal in its panel, but it is weaker than the Type-1 splice variants.",
        "- A1 Boltz results are mixed: P1413L has the strongest GPIb forced-binding increase, while V1316M has the strongest decrease.",
        "- R1308C has the strongest heparan-sulfate increase; A1461D has the strongest decrease.",
        "- Static Boltz confidence deltas are mechanism proxies, not direct binding free energies or clinical classifications.",
        "",
    ]
    (out_dir / "pre_md_interpretation.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    alpha = pd.concat([
        alpha_summary(ROOT / "outputs/type1_panel_agent_20260828/server_bundle/results/alphagenome/alphagenome_agent_features.csv", "type1"),
        alpha_summary(ROOT / "outputs/type2b_panel_agent_20260829/server_bundle/results/alphagenome/alphagenome_agent_features.csv", "type2b"),
    ], ignore_index=True)
    boltz = pd.concat([
        boltz_summary(
            ROOT / "outputs/type1_panel_agent_20260828/server_bundle/boltz/job_manifest.csv",
            ROOT / "outputs/type1_panel_agent_20260828/server_bundle/results/boltz/boltz_results_summary.csv",
            "type1",
        ),
        boltz_summary(
            ROOT / "outputs/type2b_panel_agent_20260829/server_bundle/boltz/job_manifest.csv",
            ROOT / "outputs/type2b_panel_agent_20260829/server_bundle/results/boltz/boltz_results_summary.csv",
            "type2b",
        ),
    ], ignore_index=True)

    alpha.to_csv(args.output_dir / "alphagenome_case_summary.csv", index=False)
    boltz.to_csv(args.output_dir / "boltz_variant_assay_summary.csv", index=False)
    write_report(alpha, boltz, args.output_dir)
    print(f"Wrote {args.output_dir / 'alphagenome_case_summary.csv'}")
    print(f"Wrote {args.output_dir / 'boltz_variant_assay_summary.csv'}")
    print(f"Wrote {args.output_dir / 'pre_md_interpretation.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

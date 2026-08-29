from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd

from scripts.pipeline.run_type1_10case_alphagenome import (
    RAW_OUTPUTS,
    build_ontology_plan,
    summarize_raw_delta,
)
from scripts.pipeline.ingest_type1_10case_results import boltz_rows


def test_metadata_plan_selects_relevant_tracks_and_preserves_global_outputs() -> None:
    metadata = SimpleNamespace()
    for output_name in RAW_OUTPUTS:
        if output_name == "SPLICE_SITES":
            frame = pd.DataFrame({"name": ["donor+", "acceptor+"], "ontology_curie": [None, None]})
        elif output_name == "PROCAP":
            frame = pd.DataFrame({
                "name": ["unrelated PROCAP"],
                "ontology_curie": ["CL:9999999"],
                "biosample_name": ["unrelated cell"],
            })
        else:
            frame = pd.DataFrame({
                "name": ["CL:0000115 assay", "CL:9999999 assay"],
                "ontology_curie": ["CL:0000115", "CL:9999999"],
                "biosample_name": ["endothelial cell", "unrelated cell"],
            })
        setattr(metadata, output_name.casefold(), frame)

    plan, inventory = build_ontology_plan(metadata, ["CL:0000115"])
    assert len(plan) == 11
    assert plan["RNA_SEQ"]["ontology_terms"] == ["CL:0000115"]
    assert plan["SPLICE_SITES"]["strategy"] == "global_non_biosample"
    assert plan["PROCAP"]["strategy"] == "unavailable"
    assert inventory["selected_for_type1_panel"].fillna(False).sum() == 11


def test_raw_delta_summary_retains_alt_minus_ref_direction() -> None:
    ref = SimpleNamespace(values=np.array([[1.0, 3.0], [2.0, 0.0]]))
    alt = SimpleNamespace(values=np.array([[0.0, 5.0], [2.5, -1.0]]))
    row = summarize_raw_delta("CASE_T1_001", "RNA_SEQ", ref, alt)
    assert row["status"] == "success"
    assert row["signed_min"] == -1.0
    assert row["signed_max"] == 2.0
    assert row["abs_max"] == 2.0


def test_boltz_duplicate_variant_maps_back_to_both_patients(tmp_path) -> None:
    manifest = pd.DataFrame([
        {"job_name": "wt", "run_decision": "WT_BASELINE", "assay_key": "a1", "variant_id": "VWF_WT"},
        {"job_name": "mut", "run_decision": "RUN", "assay_key": "a1", "variant_id": "VWF_V1316M"},
    ])
    summary = pd.DataFrame([
        {"job_name": "wt", "primary_metric": "iptm", "primary_value": 0.7},
        {"job_name": "mut", "primary_metric": "iptm", "primary_value": 0.8},
    ])
    variants = pd.DataFrame([
        {"case_id": "CASE_T2B_001", "patient_id": "CASE_T2B_001", "aa_change": "V1316M"},
        {"case_id": "CASE_T2B_002", "patient_id": "CASE_T2B_002", "aa_change": "V1316M"},
    ])
    manifest_path, summary_path, variants_path = (
        tmp_path / "manifest.csv", tmp_path / "summary.csv", tmp_path / "variants.csv"
    )
    manifest.to_csv(manifest_path, index=False)
    summary.to_csv(summary_path, index=False)
    variants.to_csv(variants_path, index=False)
    rows = boltz_rows(summary_path, manifest_path, variants_path, expected_jobs=2, expected_wt=1)
    assert {row["case_id"] for row in rows} == {"CASE_T2B_001", "CASE_T2B_002"}
    assert all("delta_vs_WT=0.1" in row["conclusion"] for row in rows)

from pathlib import Path

import pandas as pd

from src.vwd_clinical_agent.schemas import VariantContext
from src.vwd_clinical_agent.tools.computational_panel import LocalComputationalPanelProvider


ROOT = Path(__file__).resolve().parents[2]


def test_existing_alphagenome_and_boltz_evidence_are_embedded() -> None:
    provider = LocalComputationalPanelProvider(ROOT)
    variant = VariantContext(
        source_row_id="CASE_T1_001",
        variant_index=1,
        gene="VWF",
        hgvs_c="c.4382C>A",
        hgvs_p="p.Ala1461Asp",
        genome_build="GRCh38",
        hg38_chromosome="chr12",
        hg38_position=6019036,
        hg38_ref="G",
        hg38_alt="T",
        alphagenome_request_status="READY",
        boltz_request_status="READY",
    )
    items, status = provider.collect(patient_id="CASE_T1_001", variant=variant)
    assert {item.source for item in items} == {
        "alphagenome_full_profile",
        "boltz2_functional_panel",
    }
    assert status["alphagenome"] == "returned_full_profile_embedded"
    assert status["boltz"] == "returned_panel_embedded"
    alpha = next(item for item in items if item.source == "alphagenome_full_profile")
    assert any("research-only" in limitation.casefold() for limitation in alpha.limitations)


def test_case_t1_007_vwf_splice_variant_is_model_ready_despite_f8_comorbidity() -> None:
    provider = LocalComputationalPanelProvider(ROOT)
    variant = VariantContext(
        source_row_id="CASE_T1_007",
        variant_index=1,
        gene="VWF",
        hgvs_c="c.3379+1G>A",
        genome_build="GRCh38",
        hg38_chromosome="chr12",
        hg38_position=6023630,
        hg38_ref="C",
        hg38_alt="T",
        alphagenome_request_status="READY",
        boltz_request_status="NOT_APPLICABLE_NON_MISSENSE",
    )
    items, status = provider.collect(patient_id="CASE_T1_007", variant=variant)
    assert [item.source for item in items] == ["alphagenome_full_profile"]
    assert status["alphagenome"] == "returned_full_profile_embedded"
    assert status["boltz"] == "NOT_APPLICABLE_NON_MISSENSE"


def test_returned_multi_variant_evidence_is_variant_specific(tmp_path, monkeypatch) -> None:
    returned = tmp_path / "returned.csv"
    pd.DataFrame([
        {
            "case_id": "CASE_T2B_005", "source": "alphagenome_full_profile",
            "source_record_id": "CASE_T2B_005_VARIANT_1:alphagenome-full-profile",
            "query": "chr12:6019036:G>T", "conclusion": "A1461D alpha",
        },
        {
            "case_id": "CASE_T2B_005", "source": "boltz2_functional_panel",
            "source_record_id": "CASE_T2B_005:VWF_D2449N:boltz2-functional-panel",
            "query": "D2449N", "conclusion": "D2449N boltz",
        },
    ]).to_csv(returned, index=False)
    monkeypatch.setenv("VWF_COMPUTATIONAL_EVIDENCE", str(returned))
    provider = LocalComputationalPanelProvider(ROOT)
    first = VariantContext(
        source_row_id="CASE_T2B_005_VARIANT_1", variant_index=1, gene="VWF",
        hgvs_p="p.Ala1461Asp", hg38_chromosome="chr12", hg38_position=6019036,
        hg38_ref="G", hg38_alt="T", alphagenome_request_status="READY", boltz_request_status="READY",
    )
    second = VariantContext(
        source_row_id="CASE_T2B_005_VARIANT_2_BENIGN", variant_index=2, gene="VWF",
        hgvs_p="p.Asp2449Asn", alphagenome_request_status="READY", boltz_request_status="READY",
    )
    first_items, _ = provider.collect(patient_id="CASE_T2B_005", variant=first)
    second_items, _ = provider.collect(patient_id="CASE_T2B_005", variant=second)
    assert [item.conclusion for item in first_items if item.source == "alphagenome_full_profile"] == ["A1461D alpha"]
    assert all(item.conclusion != "D2449N boltz" for item in first_items)
    assert [item.conclusion for item in second_items if item.source == "boltz2_functional_panel"] == ["D2449N boltz"]
    assert all(item.conclusion != "A1461D alpha" for item in second_items)

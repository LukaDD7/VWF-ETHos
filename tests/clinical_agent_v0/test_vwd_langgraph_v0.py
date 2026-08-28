from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.azure import DeterministicLLMProvider
from src.vwd_clinical_agent.graph import build_workflow
from src.vwd_clinical_agent.workbook import LocalWorkbookProvider, observed_value, parse_case_ref


WORKBOOK = ROOT / "data/clinical_agent_pilot/vwd_agentic_workflow_deidentified_v3.xlsx"


def test_workbook_patient_and_variant_counts() -> None:
    cases, audit = LocalWorkbookProvider(WORKBOOK).load_cases()
    assert audit["unique_patients"] == 47
    assert audit["first_level_rows"] == 59
    assert audit["variant_rows"] == 59
    assert len(cases) == 47


def test_benign_suffix_is_not_a_new_patient() -> None:
    assert parse_case_ref("CASE_021_VARIANT_2_BENIGN") == ("CASE_021", 2, True)
    cases, _ = LocalWorkbookProvider(WORKBOOK).load_cases()
    case = next(case for case in cases if case.patient_id == "CASE_021")
    assert len(case.variants) == 2
    assert case.variants[1].benign_reported is True


def test_zero_is_observed_and_na_is_missing() -> None:
    zero = observed_value(0)
    assert zero.observed is True and zero.value == 0
    missing = observed_value("NA")
    assert missing.observed is False and missing.missing_reason == "not_available"


def test_graph_smoke_and_safety_gate() -> None:
    cases, _ = LocalWorkbookProvider(WORKBOOK).load_cases()
    case = next(item for item in cases if item.patient_id == "CASE_024")
    graph = build_workflow(DeterministicLLMProvider())
    result = graph.invoke(
        {
            "run_id": "test-run",
            "case": case,
            "mode": "retrospective",
            "trace": [],
            "provider_calls": [],
        }
    )
    assert result["status"] == "completed"
    assert result["final_opinion"].abstention is True
    assert result["final_opinion"].expert_review_required is True
    assert all(variant.phase_status == "unknown" for variant in result["variants"])
    assert result["second_level_status"] == "not_observed"
    assert all(action.action_code != "UNKNOWN_ACTION" for action in result["recommended_actions"])

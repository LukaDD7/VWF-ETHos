from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from src.vwd_clinical_agent.mechanistic_tasks import (
    DecisionContext,
    MechanismAssessment,
    MechanismTaskResult,
    ProtocolRegistry,
    QCResult,
    TaskProposal,
    VariantSpec,
    create_task_request,
    result_template,
    result_to_fhir_bundle,
    task_to_fhir,
    validate_task_result,
)


ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = ROOT / "protocols/vwd_mechanistic_v1/registry.json"


def proposal() -> TaskProposal:
    return TaskProposal(
        patient_id="EVAL-TEST",
        variant=VariantSpec(
            hgvs_c="c.3946G>A",
            hgvs_p="p.Val1316Met",
            zygosity="heterozygous",
            domain="A1",
            variant_class="missense",
        ),
        mechanism_id="M6_2B_A1_GOF",
        protocol_id="VWF_A1_2B_GOF_LEGAN2023_V1",
        hypothesis_rationale=(
            "The current evidence leaves a competitive A1 gain-versus-loss-of-function question."
        ),
        competing_hypotheses=["M7_2M_A1_DISORDER", "M8_2M_A1_HYPERSTABLE", "M_UNKNOWN"],
        decision_context=DecisionContext(
            evidence_retrieval_complete=True,
            evidence_sufficiency="ambiguous",
            candidate_hypotheses=[
                "M1_SECRETION",
                "M6_2B_A1_GOF",
                "M7_2M_A1_DISORDER",
                "M8_2M_A1_HYPERSTABLE",
                "M_UNKNOWN",
            ],
            current_subtype_probabilities={
                "type_1": 0.02,
                "type_2A": 0.04,
                "type_2B": 0.46,
                "type_2M": 0.43,
                "type_2N": 0.02,
                "type_3": 0.01,
                "unresolved": 0.02,
            },
            current_evidence=["ClinGen/ClinVar/PubMed retrieval complete"],
            missing_evidence=["low-dose RIPA", "multimer analysis"],
            why_this_task_changes_decision=(
                "A calibrated activated-A1 pattern would separate the leading 2B and 2M hypotheses."
            ),
            expected_information_gain=0.42,
        ),
    )


def completed_result(request, registry: ProtocolRegistry) -> MechanismTaskResult:
    protocol = registry.protocols[request.acquisition.protocol_id]
    measurements = []
    for item in protocol.measurements:
        if item.comparison == "case_vs_matched_wt":
            case_value = 2.0
            reference_value = 1.0
            delta = 1.0
        else:
            case_value = {"activated": 0.7}
            reference_value = {"activated": 0.2}
            delta = {"activated": 0.5}
        measurements.append(
            {
                "measurement_id": item.measurement_id,
                "definition_id": item.definition_id,
                "status": "available",
                "unit": registry.measurements[item.definition_id].default_unit,
                "case_value": case_value,
                "reference_value": reference_value,
                "delta_case_minus_reference": delta,
                "comparison_kind": item.comparison,
                "replicate_values": [case_value],
                "provenance": {"roi": item.roi},
            }
        )
    return MechanismTaskResult(
        task_id=request.task_id,
        request_digest=request.request_digest,
        protocol_id=protocol.protocol_id,
        protocol_version=protocol.version,
        protocol_digest=protocol.digest,
        status="completed",
        completed_at="2026-09-10T00:00:00+00:00",
        software_versions={"gromacs": "2025.4"},
        runtime={"wall_seconds": 1200},
        qc=QCResult(
            passed=True,
            failures=[],
            replicate_consistency=0.9,
            equilibration_passed=True,
            structure_integrity_passed=True,
            out_of_distribution=False,
        ),
        measurements=measurements,
        benchmark={
            "similarity_to_positive_controls": 0.82,
            "similarity_to_negative_controls": 0.24,
            "calibration_set_id": "AIM-A1-v1",
            "out_of_distribution": False,
        },
        mechanism_assessment=MechanismAssessment(
            supports=["M6_2B_A1_GOF"],
            contradicts=[],
            indeterminate=[],
            confidence="moderate",
            reasons=["Composite measurement pattern matches the positive-control set."],
        ),
        artifacts=[],
        limitations=["Equilibrium MD does not measure platelet binding under flow."],
    )


def test_registry_integrity_and_open_set_candidates() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    assert "M_UNKNOWN" in registry.candidate_mechanisms("A1", "missense")
    assert len(registry.protocols) == 3
    assert all(protocol.digest.startswith("sha256:") for protocol in registry.protocols.values())


def test_task_is_protocol_locked_and_fhir_requested() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    request = create_task_request(
        proposal(),
        registry,
        source_commit="a" * 40,
        task_id="mtask-" + "1" * 32,
    )
    protocol = registry.protocols[request.acquisition.protocol_id]
    assert request.acquisition.requested_measurements == protocol.required_measurement_ids
    task = task_to_fhir(request, registry)
    assert task["resourceType"] == "Task"
    assert task["status"] == "requested"
    assert protocol.protocol_id in task["instantiatesCanonical"]


def test_submission_requires_retrieval_unknown_and_information_gain() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    raw = proposal().model_dump(mode="json")
    raw["decision_context"]["evidence_retrieval_complete"] = False
    with pytest.raises(ValueError, match="retrieval must complete"):
        create_task_request(TaskProposal.model_validate(raw), registry, source_commit="a" * 40)

    raw = proposal().model_dump(mode="json")
    raw["decision_context"]["candidate_hypotheses"].remove("M_UNKNOWN")
    with pytest.raises(ValueError, match="M_UNKNOWN"):
        create_task_request(TaskProposal.model_validate(raw), registry, source_commit="a" * 40)

    raw = proposal().model_dump(mode="json")
    raw["decision_context"]["expected_information_gain"] = 0.01
    with pytest.raises(ValueError, match="information gain"):
        create_task_request(TaskProposal.model_validate(raw), registry, source_commit="a" * 40)


def test_noncompetitive_target_is_rejected() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    raw = proposal().model_dump(mode="json")
    raw["decision_context"]["current_subtype_probabilities"] = {
        "type_1": 0.01,
        "type_2A": 0.01,
        "type_2B": 0.05,
        "type_2M": 0.90,
        "type_2N": 0.01,
        "type_3": 0.01,
        "unresolved": 0.01,
    }
    with pytest.raises(ValueError, match="noncompetitive"):
        create_task_request(TaskProposal.model_validate(raw), registry, source_commit="a" * 40)


def test_completed_result_requires_exact_measurements_and_becomes_fhir_observation() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    request = create_task_request(
        proposal(),
        registry,
        source_commit="a" * 40,
        task_id="mtask-" + "2" * 32,
    )
    result = completed_result(request, registry)
    validate_task_result(result, request, registry)
    bundle = result_to_fhir_bundle(result, request, registry)
    resources = [entry["resource"] for entry in bundle["entry"]]
    assert [item["resourceType"] for item in resources] == ["Task", "Observation"]
    assert resources[0]["status"] == "completed"
    assert resources[1]["basedOn"] == [{"reference": f"Task/{request.task_id}"}]

    result.measurements = result.measurements[1:]
    with pytest.raises(ValueError, match="missing required measurements"):
        validate_task_result(result, request, registry)


def test_result_rejects_wrong_protocol_digest() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    request = create_task_request(proposal(), registry, source_commit="a" * 40)
    result = completed_result(request, registry)
    result.protocol_digest = "sha256:" + "0" * 64
    with pytest.raises(ValueError, match="exact task/protocol snapshot"):
        validate_task_result(result, request, registry)


def test_untouched_result_template_is_intentionally_invalid() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    request = create_task_request(proposal(), registry, source_commit="a" * 40)
    with pytest.raises(ValidationError, match="template_only"):
        MechanismTaskResult.model_validate(result_template(request, registry))

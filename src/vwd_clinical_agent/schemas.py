from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum
from typing import Annotated, Any, Literal, Optional, TypedDict
import operator

from pydantic import BaseModel, Field


ActionCode = Literal[
    "COMPLETE_FIRST_LEVEL",
    "RIPA",
    "VWF_MULTIMER",
    "VWF_CB",
    "VWF_FVIIIB",
    "VWF_PP",
    "DDAVP_0_1_4H",
    "WAIT_SECOND_LEVEL_RESULTS",
    "WAIT_GENETICS",
    "EXPERT_REVIEW",
]


class MissingReason(str, Enum):
    NOT_RECORDED = "not_recorded"
    NOT_AVAILABLE = "not_available"
    NOT_DONE = "not_done"
    PENDING = "pending"


class ObservedValue(BaseModel):
    value: Optional[float] = None
    raw: str = ""
    observed: bool = False
    unit: Optional[str] = None
    reference_range: Optional[str] = None
    missing_reason: Optional[MissingReason] = None


class FirstLevelLabs(BaseModel):
    vwf_ag: ObservedValue
    vwf_act: ObservedValue
    fviii_c: ObservedValue
    platelet_count: ObservedValue


class VariantContext(BaseModel):
    source_row_id: str
    variant_index: int
    hgvs_c: Optional[str] = None
    hgvs_p: Optional[str] = None
    chromosomal_position: Optional[str] = None
    benign_reported: bool = False
    phase_status: Literal["single", "unknown", "not_applicable"] = "single"


class PatientCase(BaseModel):
    patient_id: str
    episode_id: str
    source_row_ids: list[str]
    first_level: FirstLevelLabs
    variants: list[VariantContext]


class ClinicalAction(BaseModel):
    action_code: ActionCode
    rank: int
    status: Literal["recommended", "not_available", "deferred", "not_applicable"]
    availability: Literal["unknown", "available", "not_available"]
    clinical_hypothesis: list[str]
    expected_discriminator: str
    rationale: str
    requires_human_confirmation: bool = True
    provenance: list[str]


class EvidenceItem(BaseModel):
    source: str
    source_record_id: str
    query: str
    retrieved_at: str
    source_version: str
    supports: list[str] = Field(default_factory=list)
    refutes: list[str] = Field(default_factory=list)
    conclusion: str
    confidence: Optional[float] = None
    evidence_level: Literal["research_input", "clinical_database", "guideline", "mechanism_model"] = "research_input"
    provenance_url: Optional[str] = None
    citation: Optional[str] = None
    raw_payload_hash: Optional[str] = None
    limitations: list[str] = Field(default_factory=list)
    raw_excerpt_locator: Optional[str] = None


class SafetyFlag(BaseModel):
    code: str
    severity: Literal["info", "warning", "critical"]
    message: str
    requires_expert_review: bool = True


class FinalClinicalPackage(BaseModel):
    patient_id: str
    candidate_subtypes: list[str]
    confidence: Literal["low", "moderate", "high"]
    opinion: str
    abstention: bool
    expert_review_required: bool
    recommended_actions: list[ClinicalAction]
    supporting_evidence: list[str]
    contradicting_evidence: list[str]
    missing_information: list[str]
    limitations: list[str]
    provenance: list[str]


class EvidenceConflict(BaseModel):
    conflict_id: str
    conflict_type: Literal[
        "classification_vs_population_frequency",
        "classification_vs_prediction",
        "database_vs_database",
        "missing_discriminator",
        "insufficient_evidence",
    ]
    severity: Literal["info", "warning", "critical"]
    evidence_refs: list[str]
    description: str
    explanation: str
    recommended_action: Literal["expert_review", "additional_evidence", "accept_uncertainty"]


class TraceEvent(BaseModel):
    run_id: str
    case_id: str
    node: str
    timestamp: str
    detail: dict[str, Any] = Field(default_factory=dict)


class VWDWorkflowState(TypedDict, total=False):
    run_id: str
    case: PatientCase
    mode: Literal["retrospective", "interactive"]
    decision_time: Optional[str]
    first_level: FirstLevelLabs
    variants: list[VariantContext]
    vwf_act_ag_ratio: Optional[float]
    missing_critical_fields: list[str]
    pre_genetic_hypotheses: list[str]
    pre_genetic_route: Literal["second_level", "wait_genetics"]
    recommended_actions: list[ClinicalAction]
    second_level_environment: Literal["provided", "auto_unavailable"]
    second_level_status: Literal["not_observed", "waiting", "observed", "not_available"]
    evidence_plan: list[dict[str, Any]]
    evidence_items: list[EvidenceItem]
    provider_calls: Annotated[list[dict[str, Any]], operator.add]
    second_level_bundle: Any | None
    fhir_evidence_bundle: Any | None
    evidence_conflicts: list[EvidenceConflict]
    evidence_missing: list[str]
    evidence_conflict_summary: str
    acmg_evidence_hints: list[str]
    candidate_subtypes: list[str]
    safety_flags: list[SafetyFlag]
    llm_summary: Optional[str]
    llm_provider: Optional[str]
    final_opinion: Optional[FinalClinicalPackage]
    status: Literal[
        "running",
        "waiting_second_level",
        "waiting_genetics",
        "expert_review",
        "completed",
        "failed",
    ]
    trace: Annotated[list[TraceEvent], operator.add]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()

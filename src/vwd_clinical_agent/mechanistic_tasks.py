"""Versioned, mechanism-guided computational task handoff.

The local clinical agent decides *what* mechanism should be tested.  A remote
runner decides none of the scientific parameters: it executes the referenced
protocol and returns a result conforming to the same registry snapshot.

Git is deliberately treated as a transport/audit envelope, not as a scheduler.
The models in this module are independent of a particular Git host and can be
used by both the networked local agent and an offline GPU worker.
"""

from __future__ import annotations

import base64
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Literal
from uuid import uuid4

from pydantic import BaseModel, ConfigDict, Field, model_validator


SCHEMA_VERSION = "1.0.0"
DEFAULT_REGISTRY_PATH = (
    Path(__file__).resolve().parents[2]
    / "protocols"
    / "vwd_mechanistic_v1"
    / "registry.json"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonical_json(value: Any) -> str:
    """Serialize a contract deterministically for content-addressed handoff."""

    if isinstance(value, BaseModel):
        value = value.model_dump(mode="json")
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=lambda item: item.model_dump(mode="json") if isinstance(item, BaseModel) else str(item),
    )


def content_digest(value: Any) -> str:
    return "sha256:" + sha256(canonical_json(value).encode("utf-8")).hexdigest()


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", protected_namespaces=())


class MechanismDefinition(StrictModel):
    mechanism_id: str
    subtype_targets: list[str]
    applicable_domains: list[str]
    applicable_variant_classes: list[str]
    causal_chain: list[str]
    experimental_endpoint: str
    computational_observability: Literal["high", "medium", "low", "none"]
    allowed_protocols: list[str] = Field(default_factory=list)
    limitations: list[str] = Field(default_factory=list)


class MeasurementDefinition(StrictModel):
    definition_id: str
    layer: Literal["generic", "mechanism", "comparison", "qc"]
    source: str
    semantic: str
    value_shape: Literal["scalar", "distribution", "mapping"]
    default_unit: str | None = None
    limitations: list[str] = Field(default_factory=list)


class ProtocolMeasurement(StrictModel):
    measurement_id: str
    definition_id: str
    required: bool = True
    roi: dict[str, Any] = Field(default_factory=dict)
    comparison: Literal["case_vs_matched_wt", "within_case", "calibration_set"]
    expected_pattern: str


class ProtocolControl(StrictModel):
    control_id: str
    role: Literal["matched_wt", "positive_reference", "negative_reference", "boundary_reference"]
    description: str


class ProtocolDefinition(StrictModel):
    protocol_id: str
    version: str
    title: str
    status: Literal["active", "draft", "reserved"]
    submission_enabled: bool
    research_only: bool = True
    mechanisms: list[str]
    eligible_domains: list[str]
    eligible_variant_classes: list[str]
    tool_family: list[str]
    starting_model: dict[str, Any]
    variant_construction: dict[str, Any]
    execution_tiers: dict[str, dict[str, Any]]
    fixed_parameters: dict[str, Any]
    controls: list[ProtocolControl]
    measurements: list[ProtocolMeasurement]
    qc_gates: list[str]
    provenance: list[dict[str, str]]
    interpretation_boundary: list[str]

    @property
    def digest(self) -> str:
        return content_digest(self)

    @property
    def required_measurement_ids(self) -> list[str]:
        return [item.measurement_id for item in self.measurements if item.required]


class PlanningPolicy(StrictModel):
    minimum_expected_information_gain: float = Field(ge=0, le=1)
    dominant_posterior_threshold: float = Field(ge=0, le=1)
    minimum_competing_hypotheses: int = Field(ge=1)
    require_unknown_hypothesis: bool = True


class RegistryDocument(StrictModel):
    schema_version: str
    registry_id: str
    version: str
    canonical_base: str
    research_only: bool = True
    planning_policy: PlanningPolicy
    mechanisms: list[MechanismDefinition]
    measurements: list[MeasurementDefinition]
    protocols: list[ProtocolDefinition]


class ProtocolRegistry:
    def __init__(self, document: RegistryDocument, source_path: Path | None = None):
        self.document = document
        self.source_path = source_path
        self.mechanisms = {item.mechanism_id: item for item in document.mechanisms}
        self.measurements = {item.definition_id: item for item in document.measurements}
        self.protocols = {item.protocol_id: item for item in document.protocols}
        self._validate_integrity()

    @classmethod
    def load(cls, path: str | Path = DEFAULT_REGISTRY_PATH) -> "ProtocolRegistry":
        source = Path(path)
        return cls(
            RegistryDocument.model_validate_json(source.read_text(encoding="utf-8")),
            source,
        )

    @property
    def digest(self) -> str:
        return content_digest(self.document)

    def _validate_integrity(self) -> None:
        if len(self.mechanisms) != len(self.document.mechanisms):
            raise ValueError("Duplicate mechanism_id in registry")
        if len(self.measurements) != len(self.document.measurements):
            raise ValueError("Duplicate definition_id in measurement registry")
        if len(self.protocols) != len(self.document.protocols):
            raise ValueError("Duplicate protocol_id in registry")
        if "M_UNKNOWN" not in self.mechanisms:
            raise ValueError("Registry must retain M_UNKNOWN to prevent closed-world classification")

        for mechanism in self.document.mechanisms:
            unknown = set(mechanism.allowed_protocols) - set(self.protocols)
            if unknown:
                raise ValueError(f"{mechanism.mechanism_id} references unknown protocols: {sorted(unknown)}")
        for protocol in self.document.protocols:
            unknown_mechanisms = set(protocol.mechanisms) - set(self.mechanisms)
            if unknown_mechanisms:
                raise ValueError(f"{protocol.protocol_id} references unknown mechanisms: {sorted(unknown_mechanisms)}")
            measurement_ids = [item.measurement_id for item in protocol.measurements]
            if len(measurement_ids) != len(set(measurement_ids)):
                raise ValueError(f"Duplicate measurement_id in {protocol.protocol_id}")
            for measurement in protocol.measurements:
                if measurement.definition_id not in self.measurements:
                    raise ValueError(
                        f"{protocol.protocol_id}/{measurement.measurement_id} references unknown "
                        f"definition {measurement.definition_id}"
                    )
            roles = {control.role for control in protocol.controls}
            if "matched_wt" not in roles or "positive_reference" not in roles:
                raise ValueError(f"{protocol.protocol_id} requires matched WT and positive-reference controls")

    def candidate_mechanisms(self, domain: str, variant_class: str) -> list[str]:
        candidates = []
        normalized_domain = domain.casefold()
        normalized_class = variant_class.casefold()
        for item in self.document.mechanisms:
            domain_match = "*" in item.applicable_domains or any(
                normalized_domain == candidate.casefold() for candidate in item.applicable_domains
            )
            class_match = "*" in item.applicable_variant_classes or any(
                normalized_class == candidate.casefold() for candidate in item.applicable_variant_classes
            )
            if domain_match and class_match:
                candidates.append(item.mechanism_id)
        if "M_UNKNOWN" not in candidates:
            candidates.append("M_UNKNOWN")
        return candidates

    def fhir_plan_definition(self) -> dict[str, Any]:
        actions = []
        for protocol in self.document.protocols:
            actions.append(
                {
                    "id": protocol.protocol_id,
                    "title": protocol.title,
                    "description": "Mechanism-specific computational acquisition; research use only.",
                    "code": {"text": ", ".join(protocol.mechanisms)},
                    "definitionCanonical": (
                        f"{self.document.canonical_base}/ActivityDefinition/{protocol.protocol_id}"
                        f"|{protocol.version}"
                    ),
                }
            )
        return {
            "resourceType": "PlanDefinition",
            "id": self.document.registry_id,
            "url": f"{self.document.canonical_base}/PlanDefinition/{self.document.registry_id}",
            "version": self.document.version,
            "name": "VWDMechanismGuidedAcquisition",
            "title": "VWD mechanism-guided computational evidence acquisition",
            "status": "active",
            "experimental": True,
            "type": {"text": "clinical-protocol"},
            "description": (
                "Bounded action space for selecting a versioned computational experiment after "
                "clinical and low-cost genetic/literature evidence retrieval."
            ),
            "action": actions,
        }

    def fhir_activity_definition(self, protocol_id: str) -> dict[str, Any]:
        protocol = self.protocols[protocol_id]
        return {
            "resourceType": "ActivityDefinition",
            "id": protocol.protocol_id,
            "url": f"{self.document.canonical_base}/ActivityDefinition/{protocol.protocol_id}",
            "version": protocol.version,
            "name": protocol.protocol_id.replace("_", ""),
            "title": protocol.title,
            "status": "active" if protocol.status == "active" else "draft",
            "experimental": True,
            "kind": "Task",
            "code": {"text": "VWD mechanism-specific computational experiment"},
            "description": (
                f"Execute {protocol.protocol_id}@{protocol.version}; registry digest {self.digest}; "
                f"protocol digest {protocol.digest}."
            ),
            "extension": [
                {
                    "url": f"{self.document.canonical_base}/StructureDefinition/protocol-digest",
                    "valueString": protocol.digest,
                },
                {
                    "url": f"{self.document.canonical_base}/StructureDefinition/research-only",
                    "valueBoolean": True,
                },
            ],
        }


class VariantSpec(StrictModel):
    gene: Literal["VWF"] = "VWF"
    hgvs_c: str | None = None
    hgvs_p: str
    zygosity: str
    domain: str
    variant_class: str = "missense"


class DecisionContext(StrictModel):
    evidence_retrieval_complete: bool
    evidence_sufficiency: Literal["insufficient", "ambiguous", "adequate"]
    candidate_hypotheses: list[str]
    current_subtype_probabilities: dict[str, float]
    current_evidence: list[str] = Field(default_factory=list)
    missing_evidence: list[str]
    why_this_task_changes_decision: str = Field(min_length=20)
    expected_information_gain: float = Field(ge=0, le=1)

    @model_validator(mode="after")
    def validate_probabilities(self) -> "DecisionContext":
        if not self.current_subtype_probabilities:
            raise ValueError("current_subtype_probabilities must not be empty")
        if any(value < 0 or value > 1 for value in self.current_subtype_probabilities.values()):
            raise ValueError("Subtype probabilities must be in [0, 1]")
        total = sum(self.current_subtype_probabilities.values())
        if not math.isclose(total, 1.0, abs_tol=0.02):
            raise ValueError(f"Subtype probabilities must sum to 1 (received {total:.4f})")
        return self


class TaskProposal(StrictModel):
    patient_id: str = Field(pattern=r"^EVAL-[A-Za-z0-9-]+$")
    variant: VariantSpec
    mechanism_id: str
    protocol_id: str
    hypothesis_rationale: str = Field(min_length=20)
    competing_hypotheses: list[str]
    decision_context: DecisionContext
    execution_tier: Literal["pilot", "production"] = "pilot"
    priority: Literal["routine", "urgent", "asap", "stat"] = "routine"
    resource_budget: dict[str, Any] = Field(default_factory=dict)


class TaskHypothesis(StrictModel):
    mechanism_id: str
    subtype_targets: list[str]
    rationale: str
    competing_hypotheses: list[str]


class TaskAcquisition(StrictModel):
    protocol_id: str
    protocol_version: str
    protocol_digest: str
    registry_digest: str
    tool_family: list[str]
    controls: list[str]
    requested_measurements: list[str]


class TaskExecution(StrictModel):
    tier: Literal["pilot", "production"]
    priority: Literal["routine", "urgent", "asap", "stat"]
    resource_budget: dict[str, Any]
    status: Literal["requested"] = "requested"


class MechanismTaskRequest(StrictModel):
    schema_version: Literal[SCHEMA_VERSION] = SCHEMA_VERSION
    task_id: str = Field(pattern=r"^mtask-[0-9a-f]{32}$")
    authored_at: str
    source_commit: str = Field(pattern=r"^[0-9a-f]{7,64}$")
    patient_id: str = Field(pattern=r"^EVAL-[A-Za-z0-9-]+$")
    variant: VariantSpec
    hypothesis: TaskHypothesis
    acquisition: TaskAcquisition
    decision_context: DecisionContext
    execution: TaskExecution
    request_digest: str


class QCResult(StrictModel):
    passed: bool
    failures: list[str] = Field(default_factory=list)
    replicate_consistency: float | None = Field(default=None, ge=0, le=1)
    equilibration_passed: bool | None = None
    structure_integrity_passed: bool | None = None
    out_of_distribution: bool | None = None


class MeasurementResult(StrictModel):
    measurement_id: str
    definition_id: str
    status: Literal["available", "unavailable", "failed_qc"]
    unit: str | None = None
    case_value: Any = None
    reference_value: Any = None
    delta_case_minus_reference: Any = None
    comparison_kind: Literal["case_vs_matched_wt", "within_case", "calibration_set"]
    replicate_values: list[Any] = Field(default_factory=list)
    provenance: dict[str, Any] = Field(default_factory=dict)


class BenchmarkResult(StrictModel):
    similarity_to_positive_controls: float | None = Field(default=None, ge=0, le=1)
    similarity_to_negative_controls: float | None = Field(default=None, ge=0, le=1)
    calibration_set_id: str | None = None
    out_of_distribution: bool | None = None


class MechanismAssessment(StrictModel):
    supports: list[str] = Field(default_factory=list)
    contradicts: list[str] = Field(default_factory=list)
    indeterminate: list[str] = Field(default_factory=list)
    confidence: Literal["low", "moderate", "high"]
    reasons: list[str]


class ArtifactResult(StrictModel):
    path: str
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    media_type: str
    role: str


class MechanismTaskResult(StrictModel):
    schema_version: Literal[SCHEMA_VERSION] = SCHEMA_VERSION
    task_id: str
    request_digest: str
    protocol_id: str
    protocol_version: str
    protocol_digest: str
    status: Literal["completed", "failed", "inconclusive"]
    completed_at: str
    software_versions: dict[str, str]
    model_versions: dict[str, str] = Field(default_factory=dict)
    runtime: dict[str, Any]
    qc: QCResult
    measurements: list[MeasurementResult]
    benchmark: BenchmarkResult
    mechanism_assessment: MechanismAssessment
    artifacts: list[ArtifactResult] = Field(default_factory=list)
    limitations: list[str] = Field(min_length=1)


def request_digest_for(request: MechanismTaskRequest) -> str:
    payload = request.model_dump(mode="json", exclude={"request_digest"})
    return content_digest(payload)


def validate_task_request(request: MechanismTaskRequest, registry: ProtocolRegistry) -> None:
    protocol = registry.protocols.get(request.acquisition.protocol_id)
    if protocol is None:
        raise ValueError(f"Unknown protocol: {request.acquisition.protocol_id}")
    if not protocol.submission_enabled or protocol.status != "active":
        raise ValueError(f"Protocol is not enabled for submission: {protocol.protocol_id}")
    if request.acquisition.protocol_version != protocol.version:
        raise ValueError("Protocol version does not match registry")
    if request.acquisition.protocol_digest != protocol.digest:
        raise ValueError("Protocol digest does not match registry")
    if request.acquisition.registry_digest != registry.digest:
        raise ValueError("Registry digest does not match registry")
    if request.request_digest != request_digest_for(request):
        raise ValueError("Request digest does not match request contents")
    if request.hypothesis.mechanism_id not in protocol.mechanisms:
        raise ValueError("Selected mechanism is not testable by selected protocol")
    if request.variant.domain.casefold() not in {item.casefold() for item in protocol.eligible_domains}:
        raise ValueError("Variant domain is not eligible for selected protocol")
    if request.variant.variant_class.casefold() not in {
        item.casefold() for item in protocol.eligible_variant_classes
    }:
        raise ValueError("Variant class is not eligible for selected protocol")

    context = request.decision_context
    policy = registry.document.planning_policy
    if not context.evidence_retrieval_complete:
        raise ValueError("Low-cost evidence retrieval must complete before computational submission")
    if context.evidence_sufficiency == "adequate":
        raise ValueError("Computational task rejected because current evidence is already adequate")
    if context.expected_information_gain < policy.minimum_expected_information_gain:
        raise ValueError("Expected information gain is below the submission threshold")
    if request.hypothesis.mechanism_id not in context.candidate_hypotheses:
        raise ValueError("Selected mechanism is absent from candidate_hypotheses")
    if policy.require_unknown_hypothesis and "M_UNKNOWN" not in context.candidate_hypotheses:
        raise ValueError("M_UNKNOWN must remain in the hypothesis space")
    if len(set(request.hypothesis.competing_hypotheses)) < policy.minimum_competing_hypotheses:
        raise ValueError("At least one competing hypothesis is required")
    if request.hypothesis.mechanism_id in request.hypothesis.competing_hypotheses:
        raise ValueError("Selected mechanism cannot compete with itself")
    unknown_competitors = set(request.hypothesis.competing_hypotheses) - set(context.candidate_hypotheses)
    if unknown_competitors:
        raise ValueError(f"Competing hypotheses are not in the candidate set: {sorted(unknown_competitors)}")

    probabilities = context.current_subtype_probabilities
    top_subtype, top_probability = max(probabilities.items(), key=lambda item: item[1])
    target_probabilities = [probabilities.get(item, 0.0) for item in request.hypothesis.subtype_targets]
    if (
        top_probability >= policy.dominant_posterior_threshold
        and top_subtype not in request.hypothesis.subtype_targets
        and max(target_probabilities, default=0.0) < 0.10
    ):
        raise ValueError("Task target is noncompetitive under the current subtype posterior")

    if request.acquisition.requested_measurements != protocol.required_measurement_ids:
        raise ValueError("Task must request the protocol's exact required measurement set")
    expected_controls = [item.control_id for item in protocol.controls]
    if request.acquisition.controls != expected_controls:
        raise ValueError("Task controls do not match the locked protocol")


def create_task_request(
    proposal: TaskProposal,
    registry: ProtocolRegistry,
    *,
    source_commit: str,
    task_id: str | None = None,
) -> MechanismTaskRequest:
    protocol = registry.protocols.get(proposal.protocol_id)
    mechanism = registry.mechanisms.get(proposal.mechanism_id)
    if protocol is None:
        raise ValueError(f"Unknown protocol: {proposal.protocol_id}")
    if mechanism is None:
        raise ValueError(f"Unknown mechanism: {proposal.mechanism_id}")
    request = MechanismTaskRequest(
        task_id=task_id or f"mtask-{uuid4().hex}",
        authored_at=utc_now(),
        source_commit=source_commit,
        patient_id=proposal.patient_id,
        variant=proposal.variant,
        hypothesis=TaskHypothesis(
            mechanism_id=proposal.mechanism_id,
            subtype_targets=mechanism.subtype_targets,
            rationale=proposal.hypothesis_rationale,
            competing_hypotheses=proposal.competing_hypotheses,
        ),
        acquisition=TaskAcquisition(
            protocol_id=protocol.protocol_id,
            protocol_version=protocol.version,
            protocol_digest=protocol.digest,
            registry_digest=registry.digest,
            tool_family=protocol.tool_family,
            controls=[item.control_id for item in protocol.controls],
            requested_measurements=protocol.required_measurement_ids,
        ),
        decision_context=proposal.decision_context,
        execution=TaskExecution(
            tier=proposal.execution_tier,
            priority=proposal.priority,
            resource_budget=proposal.resource_budget,
        ),
        request_digest="sha256:" + "0" * 64,
    )
    request.request_digest = request_digest_for(request)
    validate_task_request(request, registry)
    return request


def validate_task_result(
    result: MechanismTaskResult,
    request: MechanismTaskRequest,
    registry: ProtocolRegistry,
) -> None:
    validate_task_request(request, registry)
    protocol = registry.protocols[request.acquisition.protocol_id]
    identity = (
        result.task_id == request.task_id
        and result.request_digest == request.request_digest
        and result.protocol_id == protocol.protocol_id
        and result.protocol_version == protocol.version
        and result.protocol_digest == protocol.digest
    )
    if not identity:
        raise ValueError("Result does not match the exact task/protocol snapshot")

    result_ids = [item.measurement_id for item in result.measurements]
    if len(result_ids) != len(set(result_ids)):
        raise ValueError("Result contains duplicate measurement_id values")
    allowed = {item.measurement_id: item for item in protocol.measurements}
    unknown = set(result_ids) - set(allowed)
    if unknown:
        raise ValueError(f"Result contains measurements not requested by protocol: {sorted(unknown)}")
    for item in result.measurements:
        expected = allowed[item.measurement_id]
        if item.definition_id != expected.definition_id:
            raise ValueError(f"Definition mismatch for {item.measurement_id}")
        if item.comparison_kind != expected.comparison:
            raise ValueError(f"Comparison-kind mismatch for {item.measurement_id}")
        definition = registry.measurements[item.definition_id]
        if item.status == "available" and item.case_value is None:
            raise ValueError(f"{item.measurement_id} requires a case value")
        if item.status == "available" and item.comparison_kind in {
            "case_vs_matched_wt",
            "calibration_set",
        } and item.reference_value is None:
            raise ValueError(f"{item.measurement_id} requires a reference value")
        if item.status == "available" and definition.default_unit is not None:
            if item.unit != definition.default_unit:
                raise ValueError(f"Unit mismatch for {item.measurement_id}")
        if item.status == "available" and item.comparison_kind == "case_vs_matched_wt":
            if item.delta_case_minus_reference is None:
                raise ValueError(f"{item.measurement_id} requires case-minus-reference delta")
            if all(
                isinstance(value, (int, float)) and not isinstance(value, bool)
                for value in (item.case_value, item.reference_value, item.delta_case_minus_reference)
            ):
                expected_delta = float(item.case_value) - float(item.reference_value)
                if not math.isclose(expected_delta, float(item.delta_case_minus_reference), rel_tol=1e-6, abs_tol=1e-8):
                    raise ValueError(f"Incorrect case-minus-reference delta for {item.measurement_id}")

    if result.status == "completed" and not result.qc.passed:
        raise ValueError("Completed result must pass QC; otherwise return inconclusive or failed")
    if result.status != "completed" and result.qc.passed:
        raise ValueError("QC-passing result must use completed status")
    if result.status == "completed":
        available = {item.measurement_id for item in result.measurements if item.status == "available"}
        missing = set(protocol.required_measurement_ids) - available
        if missing:
            raise ValueError(f"QC-passing completed result is missing required measurements: {sorted(missing)}")
        if not result.software_versions or not result.runtime:
            raise ValueError("Completed result requires software versions and runtime metadata")
        tier = protocol.execution_tiers[request.execution.tier]
        expected_replicates = int(tier.get("replicates", 1))
        if expected_replicates > 1:
            if result.qc.replicate_consistency is None:
                raise ValueError("Production result requires replicate consistency")
            missing_replicates = [
                item.measurement_id
                for item in result.measurements
                if item.measurement_id in protocol.required_measurement_ids
                and len(item.replicate_values) < expected_replicates
            ]
            if missing_replicates:
                raise ValueError(
                    "Production result lacks protocol-required replicate values: "
                    f"{sorted(missing_replicates)}"
                )
    if result.status in {"failed", "inconclusive"} and not result.limitations:
        raise ValueError("Failed or inconclusive results must explain their limitations")

    buckets = {
        "supports": set(result.mechanism_assessment.supports),
        "contradicts": set(result.mechanism_assessment.contradicts),
        "indeterminate": set(result.mechanism_assessment.indeterminate),
    }
    if any(left & right for index, left in enumerate(buckets.values()) for right in list(buckets.values())[index + 1 :]):
        raise ValueError("Mechanism-assessment buckets must be disjoint")
    selected = request.hypothesis.mechanism_id
    if sum(selected in bucket for bucket in buckets.values()) != 1:
        raise ValueError("Selected hypothesis must appear in exactly one assessment bucket")
    known = set(registry.mechanisms)
    reported = set().union(*buckets.values())
    if reported - known:
        raise ValueError(f"Assessment contains unknown mechanisms: {sorted(reported - known)}")


def task_to_fhir(request: MechanismTaskRequest, registry: ProtocolRegistry) -> dict[str, Any]:
    protocol = registry.protocols[request.acquisition.protocol_id]
    base = registry.document.canonical_base
    task = {
        "resourceType": "Task",
        "id": request.task_id,
        "status": "requested",
        "intent": "order",
        "priority": request.execution.priority,
        "code": {"text": "Execute VWD mechanism-specific computational experiment"},
        "subject": {"reference": f"Patient/{request.patient_id}"},
        "authoredOn": request.authored_at,
        "instantiatesCanonical": (
            f"{base}/ActivityDefinition/{protocol.protocol_id}|{protocol.version}"
        ),
        "input": [
            {"type": {"text": "mechanism-task-request"}, "valueString": canonical_json(request)},
            {"type": {"text": "request-digest"}, "valueString": request.request_digest},
            {"type": {"text": "source-commit"}, "valueString": request.source_commit},
        ],
        "extension": [
            {
                "url": f"{base}/StructureDefinition/originating-plan",
                "valueCanonical": f"{base}/PlanDefinition/{registry.document.registry_id}|{registry.document.version}",
            },
            {
                "url": f"{base}/StructureDefinition/research-only",
                "valueBoolean": True,
            },
        ],
    }
    return task


def result_to_fhir_bundle(
    result: MechanismTaskResult,
    request: MechanismTaskRequest,
    registry: ProtocolRegistry,
) -> dict[str, Any]:
    validate_task_result(result, request, registry)
    base = registry.document.canonical_base
    observation_id = f"mechanism-result-{request.task_id.removeprefix('mtask-')}"
    document_ids = [f"artifact-{request.task_id.removeprefix('mtask-')}-{index}" for index, _ in enumerate(result.artifacts, 1)]
    observation = {
        "resourceType": "Observation",
        "id": observation_id,
        "status": "final",
        "code": {
            "coding": [
                {
                    "system": f"{base}/CodeSystem/computational-mechanism",
                    "code": request.hypothesis.mechanism_id,
                }
            ],
            "text": "Protocol-conformant VWD computational mechanism result",
        },
        "subject": {"reference": f"Patient/{request.patient_id}"},
        "basedOn": [{"reference": f"Task/{request.task_id}"}],
        "effectiveDateTime": result.completed_at,
        "component": [
            {
                "code": {
                    "coding": [
                        {
                            "system": f"{base}/CodeSystem/computational-measurement",
                            "code": item.measurement_id,
                        }
                    ],
                    "text": item.definition_id,
                },
                "valueString": canonical_json(item),
            }
            for item in result.measurements
        ],
        "note": [
            {
                "text": canonical_json(
                    {
                        "qc": result.qc,
                        "benchmark": result.benchmark,
                        "mechanism_assessment": result.mechanism_assessment,
                        "limitations": result.limitations,
                    }
                )
            }
        ],
    }
    documents = []
    for document_id, artifact in zip(document_ids, result.artifacts):
        documents.append(
            {
                "resourceType": "DocumentReference",
                "id": document_id,
                "status": "current",
                "subject": {"reference": f"Patient/{request.patient_id}"},
                "content": [
                    {
                        "attachment": {
                            "contentType": artifact.media_type,
                            "url": artifact.path,
                            "hash": base64.b64encode(bytes.fromhex(artifact.sha256)).decode("ascii"),
                            "title": artifact.role,
                        }
                    }
                ],
            }
        )
    task = task_to_fhir(request, registry)
    task["status"] = "failed" if result.status == "failed" else "completed"
    task["lastModified"] = result.completed_at
    task["output"] = [
        {
            "type": {"text": "mechanism-measurement"},
            "valueReference": {"reference": f"Observation/{observation_id}"},
        },
        *[
            {
                "type": {"text": "computational-artifact"},
                "valueReference": {"reference": f"DocumentReference/{document_id}"},
            }
            for document_id in document_ids
        ],
    ]
    resources = [task, observation, *documents]
    return {
        "resourceType": "Bundle",
        "type": "collection",
        "timestamp": utc_now(),
        "entry": [
            {"fullUrl": f"urn:uuid:{resource['id']}", "resource": resource}
            for resource in resources
        ],
    }


def result_template(request: MechanismTaskRequest, registry: ProtocolRegistry) -> dict[str, Any]:
    protocol = registry.protocols[request.acquisition.protocol_id]
    return {
        "template_only": True,
        "schema_version": SCHEMA_VERSION,
        "task_id": request.task_id,
        "request_digest": request.request_digest,
        "protocol_id": protocol.protocol_id,
        "protocol_version": protocol.version,
        "protocol_digest": protocol.digest,
        "status": "inconclusive",
        "completed_at": utc_now(),
        "software_versions": {},
        "model_versions": {},
        "runtime": {},
        "qc": {
            "passed": False,
            "failures": ["replace with executed QC outcome"],
            "replicate_consistency": None,
            "equilibration_passed": None,
            "structure_integrity_passed": None,
            "out_of_distribution": None,
        },
        "measurements": [
            {
                "measurement_id": item.measurement_id,
                "definition_id": item.definition_id,
                "status": "unavailable",
                "unit": registry.measurements[item.definition_id].default_unit,
                "case_value": None,
                "reference_value": None,
                "delta_case_minus_reference": None,
                "comparison_kind": item.comparison,
                "replicate_values": [],
                "provenance": {"roi": item.roi},
            }
            for item in protocol.measurements
        ],
        "benchmark": {
            "similarity_to_positive_controls": None,
            "similarity_to_negative_controls": None,
            "calibration_set_id": None,
            "out_of_distribution": None,
        },
        "mechanism_assessment": {
            "supports": [],
            "contradicts": [],
            "indeterminate": [request.hypothesis.mechanism_id],
            "confidence": "low",
            "reasons": ["replace with protocol-derived assessment"],
        },
        "artifacts": [],
        "limitations": ["Template only; replace with the executed result."],
    }


class MechanismTaskStore:
    """Append-only file layout designed for one Git branch per task."""

    def __init__(self, root: str | Path):
        self.root = Path(root)

    def request_dir(self, task_id: str) -> Path:
        return self.root / "requests" / task_id

    def result_dir(self, task_id: str) -> Path:
        return self.root / "results" / task_id

    def ingested_dir(self, task_id: str) -> Path:
        return self.root / "ingested" / task_id

    @staticmethod
    def _write_once(path: Path, payload: Any) -> None:
        serialized = json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
        if path.exists():
            if path.read_text(encoding="utf-8") != serialized:
                raise FileExistsError(f"Refusing to overwrite immutable contract: {path}")
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(serialized, encoding="utf-8")

    def submit(self, request: MechanismTaskRequest, registry: ProtocolRegistry) -> list[Path]:
        validate_task_request(request, registry)
        request_path = self.request_dir(request.task_id) / "request.json"
        fhir_path = self.request_dir(request.task_id) / "task.fhir.json"
        self._write_once(request_path, request.model_dump(mode="json"))
        self._write_once(fhir_path, task_to_fhir(request, registry))
        return [request_path, fhir_path]

    def load_request(self, task_id: str, registry: ProtocolRegistry) -> MechanismTaskRequest:
        request = MechanismTaskRequest.model_validate_json(
            (self.request_dir(task_id) / "request.json").read_text(encoding="utf-8")
        )
        validate_task_request(request, registry)
        return request

    def accept_result(
        self,
        result: MechanismTaskResult,
        request: MechanismTaskRequest,
        registry: ProtocolRegistry,
    ) -> Path:
        validate_task_result(result, request, registry)
        path = self.result_dir(request.task_id) / "result.json"
        self._write_once(path, result.model_dump(mode="json"))
        return path

    def ingest_result(
        self,
        result: MechanismTaskResult,
        request: MechanismTaskRequest,
        registry: ProtocolRegistry,
    ) -> Path:
        bundle = result_to_fhir_bundle(result, request, registry)
        path = self.ingested_dir(request.task_id) / "bundle.fhir.json"
        self._write_once(path, bundle)
        return path

from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
from typing import Any, Literal
from uuid import uuid4

from pydantic import BaseModel, Field


SECOND_LEVEL_CODE_SYSTEM = "https://vwf-ethos.org/CodeSystem/vwd-second-level-action"
TOOL_CODE_SYSTEM = "https://vwf-ethos.org/CodeSystem/biomedical-evidence-tool"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def stable_id(prefix: str, *parts: Any) -> str:
    digest = sha256("\x1f".join(map(str, parts)).encode("utf-8")).hexdigest()[:20]
    return f"{prefix}-{digest}"


class Coding(BaseModel):
    system: str
    code: str
    display: str | None = None


class CodeableConcept(BaseModel):
    coding: list[Coding] = Field(default_factory=list)
    text: str | None = None


class Reference(BaseModel):
    reference: str | None = None
    display: str | None = None


class Quantity(BaseModel):
    value: float | int | None
    unit: str | None = None
    system: str | None = None
    code: str | None = None


class FHIRMeta(BaseModel):
    profile: list[str] = Field(default_factory=list)
    source: str | None = None
    versionId: str | None = None
    lastUpdated: str | None = None


class FHIRResource(BaseModel):
    resourceType: str
    id: str
    meta: FHIRMeta | None = None
    status: str | None = None
    intent: str | None = None
    category: list[CodeableConcept] = Field(default_factory=list)
    code: CodeableConcept | None = None
    subject: Reference | None = None
    basedOn: list[Reference] = Field(default_factory=list)
    effectiveDateTime: str | None = None
    issued: str | None = None
    valueCodeableConcept: CodeableConcept | None = None
    valueQuantity: Quantity | None = None
    valueString: str | None = None
    valueBoolean: bool | None = None
    component: list[dict[str, Any]] = Field(default_factory=list)
    result: list[Reference] = Field(default_factory=list)
    contained: list["FHIRResource"] = Field(default_factory=list)
    extension: list[dict[str, Any]] = Field(default_factory=list)
    target: list[Reference] = Field(default_factory=list)
    occurredDateTime: str | None = None
    recorded: str | None = None
    agent: list[dict[str, Any]] = Field(default_factory=list)
    type: str | None = None
    entry: list[dict[str, Any]] = Field(default_factory=list)
    total: int | None = None

    class Config:
        extra = "allow"


class FHIRBundle(BaseModel):
    resourceType: Literal["Bundle"] = "Bundle"
    id: str = Field(default_factory=lambda: str(uuid4()))
    type: Literal["collection", "transaction", "searchset", "document"] = "collection"
    timestamp: str = Field(default_factory=utc_now)
    entry: list[dict[str, Any]] = Field(default_factory=list)
    extension: list[dict[str, Any]] = Field(default_factory=list)

    @classmethod
    def of(cls, resources: list[FHIRResource], bundle_type: str = "collection") -> "FHIRBundle":
        return cls(
            type=bundle_type,  # type: ignore[arg-type]
            entry=[{"fullUrl": f"urn:uuid:{resource.id}", "resource": resource.model_dump()} for resource in resources],
        )

    def resources(self, resource_type: str | None = None) -> list[dict[str, Any]]:
        resources = [entry.get("resource", {}) for entry in self.entry]
        if resource_type is None:
            return resources
        return [resource for resource in resources if resource.get("resourceType") == resource_type]

    def add(self, resource: FHIRResource) -> None:
        self.entry.append({"fullUrl": f"urn:uuid:{resource.id}", "resource": resource.model_dump()})


FHIRResource.model_rebuild()


def action_code(action: str) -> CodeableConcept:
    return CodeableConcept(
        coding=[Coding(system=SECOND_LEVEL_CODE_SYSTEM, code=action, display=action.replace("_", " ").title())],
        text=action,
    )


def service_request(
    *,
    action: str,
    patient_id: str,
    status: Literal["draft", "active", "on-hold", "revoked", "completed", "entered-in-error", "unknown"] = "active",
    reason: str | None = None,
) -> FHIRResource:
    request = FHIRResource(
        resourceType="ServiceRequest",
        id=stable_id("sr", patient_id, action),
        status=status,
        intent="plan",
        category=[CodeableConcept(text="VWD second-level diagnostic workup")],
        code=action_code(action),
        subject=Reference(reference=f"Patient/{patient_id}"),
    )
    if reason:
        request.extension.append({"url": "https://vwf-ethos.org/StructureDefinition/action-rationale", "valueString": reason})
    return request


def observation(
    *,
    observation_id: str,
    patient_id: str,
    display: str,
    value: Any,
    unit: str | None = None,
    system: str = TOOL_CODE_SYSTEM,
    based_on: list[str] | None = None,
    issued: str | None = None,
    components: list[dict[str, Any]] | None = None,
) -> FHIRResource:
    resource = FHIRResource(
        resourceType="Observation",
        id=observation_id,
        status="final",
        code=CodeableConcept(coding=[Coding(system=system, code=stable_id("code", display), display=display)], text=display),
        subject=Reference(reference=f"Patient/{patient_id}"),
        effectiveDateTime=issued or utc_now(),
        issued=issued or utc_now(),
        basedOn=[Reference(reference=item) for item in based_on or []],
        component=components or [],
    )
    if isinstance(value, bool):
        resource.valueBoolean = value
    elif isinstance(value, (int, float)):
        resource.valueQuantity = Quantity(value=value, unit=unit or "proportion", code=unit)
    else:
        resource.valueString = str(value)
    return resource


def diagnostic_report(
    *,
    report_id: str,
    patient_id: str,
    display: str,
    result_ids: list[str],
    conclusion: str,
) -> FHIRResource:
    return FHIRResource(
        resourceType="DiagnosticReport",
        id=report_id,
        status="final",
        code=CodeableConcept(text=display),
        subject=Reference(reference=f"Patient/{patient_id}"),
        effectiveDateTime=utc_now(),
        issued=utc_now(),
        result=[Reference(reference=item) for item in result_ids],
        conclusion=conclusion,
    )


def operation_outcome(severity: str, code: str, diagnostics: str) -> FHIRResource:
    return FHIRResource(
        resourceType="OperationOutcome",
        id=stable_id("oo", severity, code, diagnostics),
        issue=[{"severity": severity, "code": code, "diagnostics": diagnostics}],
    )


def provenance(
    *,
    resource_id: str,
    tool_name: str,
    endpoint: str,
    request_payload: dict[str, Any] | None,
    response_payload: dict[str, Any] | None,
    occurred: str,
) -> FHIRResource:
    request_hash = sha256(str(request_payload or {}).encode("utf-8")).hexdigest()
    response_hash = sha256(str(response_payload or {}).encode("utf-8")).hexdigest()
    return FHIRResource(
        resourceType="Provenance",
        id=stable_id("prov", resource_id, tool_name, request_hash, response_hash),
        target=[Reference(reference=resource_id)],
        occurredDateTime=occurred,
        recorded=utc_now(),
        agent=[
            {
                "type": {"coding": [{"system": TOOL_CODE_SYSTEM, "code": tool_name}]},
                "who": {"display": tool_name},
                "onBehalfOf": {"display": endpoint},
            }
        ],
        extension=[
            {"url": "https://vwf-ethos.org/StructureDefinition/request-hash", "valueString": request_hash},
            {"url": "https://vwf-ethos.org/StructureDefinition/response-hash", "valueString": response_hash},
        ],
    )

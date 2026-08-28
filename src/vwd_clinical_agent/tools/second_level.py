from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pydantic import BaseModel

from .fhir import FHIRBundle, FHIRResource, diagnostic_report, observation, operation_outcome, service_request, stable_id, utc_now


SECOND_LEVEL_ACTIONS = {
    "COMPLETE_FIRST_LEVEL",
    "RIPA",
    "VWF_MULTIMER",
    "VWF_CB",
    "VWF_FVIIIB",
    "VWF_PP",
    "DDAVP_0_1_4H",
}


class SecondLevelRecord(BaseModel):
    action: str
    value: float | str | bool
    unit: str | None = None
    observed_at: str | None = None
    operator: str | None = None
    method: str | None = None


class SecondLevelFHIRStore:
    """Validated FHIR interaction store for second-level workflow actions."""

    def __init__(self, patient_id: str, bundle: FHIRBundle | None = None):
        self.patient_id = patient_id
        self.bundle = bundle or FHIRBundle(type="document")

    @classmethod
    def from_actions(cls, patient_id: str, actions: list[str], rationale: str | None = None) -> "SecondLevelFHIRStore":
        store = cls(patient_id)
        for action in actions:
            if action not in SECOND_LEVEL_ACTIONS:
                store.bundle.add(operation_outcome("error", "business-rule", f"Action outside second-level FHIR space: {action}"))
                continue
            store.bundle.add(service_request(action=action, patient_id=patient_id, reason=rationale))
        return store

    def save(self, path: str | Path) -> None:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(self.bundle.model_dump_json(indent=2), encoding="utf-8")

    @classmethod
    def load(cls, patient_id: str, path: str | Path) -> "SecondLevelFHIRStore":
        bundle = FHIRBundle.model_validate_json(Path(path).read_text(encoding="utf-8"))
        return cls(patient_id, bundle)

    def record_observation(self, record: SecondLevelRecord) -> FHIRResource:
        if record.action not in SECOND_LEVEL_ACTIONS:
            raise ValueError(f"Unsupported second-level action: {record.action}")
        request_id = self._request_id(record.action)
        if not request_id:
            raise ValueError(f"No ServiceRequest exists for {record.action}; create the request before recording a result")
        resource = observation(
            observation_id=stable_id("sl-obs", self.patient_id, record.action, record.value, record.observed_at or utc_now()),
            patient_id=self.patient_id,
            display=record.action,
            value=record.value,
            unit=record.unit,
            based_on=[f"ServiceRequest/{request_id}"],
            issued=record.observed_at,
            components=[
                {"code": {"text": "operator"}, "valueString": record.operator or "not_recorded"},
                {"code": {"text": "method"}, "valueString": record.method or "not_recorded"},
            ],
        )
        self.bundle.add(resource)
        self._set_request_status(request_id, "completed")
        return resource

    def mark_unavailable(self, action: str, reason: str) -> FHIRResource:
        if action not in SECOND_LEVEL_ACTIONS:
            raise ValueError(f"Unsupported second-level action: {action}")
        request_id = self._request_id(action)
        if not request_id:
            raise ValueError(f"No ServiceRequest exists for {action}")
        task = FHIRResource(
            resourceType="Task",
            id=stable_id("sl-task", self.patient_id, action, reason),
            status="rejected",
            intent="plan",
            code={"text": action},
            focus={"reference": f"ServiceRequest/{request_id}"},
            for_={"reference": f"Patient/{self.patient_id}"},
            statusReason={"text": reason},
        )
        self.bundle.add(task)
        self._set_request_status(request_id, "revoked")
        return task

    def finalize(self, conclusion: str) -> FHIRResource:
        requests = self.bundle.resources("ServiceRequest")
        if not requests:
            raise ValueError("No second-level ServiceRequests exist")
        unresolved = [request for request in requests if request.get("status") == "active"]
        if unresolved:
            actions = [request.get("code", {}).get("text") for request in unresolved]
            raise ValueError(f"Cannot finalize while requests remain active: {actions}")
        result_ids = [f"Observation/{resource['id']}" for resource in self.bundle.resources("Observation")]
        report = diagnostic_report(
            report_id=stable_id("sl-report", self.patient_id, conclusion),
            patient_id=self.patient_id,
            display="VWD second-level workup",
            result_ids=result_ids,
            conclusion=conclusion,
        )
        self.bundle.add(report)
        return report

    def is_final(self) -> bool:
        return any(report.get("status") == "final" for report in self.bundle.resources("DiagnosticReport"))

    def _request_id(self, action: str) -> str | None:
        expected = stable_id("sr", self.patient_id, action)
        return expected if any(request.get("id") == expected for request in self.bundle.resources("ServiceRequest")) else None

    def _set_request_status(self, request_id: str, status: str) -> None:
        for entry in self.bundle.entry:
            resource = entry.get("resource", {})
            if resource.get("id") == request_id:
                resource["status"] = status

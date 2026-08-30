from __future__ import annotations

from typing import Any

from ..schemas import ObservedValue, PatientCase
from .fhir import FHIRBundle, FHIRResource, observation


class PatientContextFHIRProvider:
    """Convert patient-level clinical context into FHIR Observations."""

    name = "patient_context_fhir"
    version = "v1"

    def collect(self, case: PatientCase) -> FHIRBundle:
        resources: list[FHIRResource] = [self._patient_resource(case)]
        labs = case.first_level
        resources.extend(
            [
                self._observed_value_observation(
                    case.patient_id,
                    "VWF antigen",
                    labs.vwf_ag,
                ),
                self._observed_value_observation(
                    case.patient_id,
                    "VWF activity",
                    labs.vwf_act,
                ),
                self._observed_value_observation(
                    case.patient_id,
                    "Factor VIII activity",
                    labs.fviii_c,
                ),
                self._observed_value_observation(
                    case.patient_id,
                    "Platelet count",
                    labs.platelet_count,
                ),
            ]
        )

        context = case.clinical_context
        resources.extend(
            [
                self._observed_value_observation(
                    case.patient_id,
                    "ISTH bleeding assessment tool score",
                    context.isth_bat,
                ),
                self._observed_value_observation(
                    case.patient_id,
                    "High-dose ristocetin platelet aggregation",
                    context.high_dose_ristocetin,
                ),
                self._text_observation(case.patient_id, "Reported symptoms", context.symptoms),
                self._text_observation(case.patient_id, "Family history", context.family_history),
                self._text_observation(case.patient_id, "Prior treatment", context.prior_treatment),
                self._text_observation(case.patient_id, "Comorbidity", context.comorbidity),
                self._text_observation(
                    case.patient_id,
                    "Interpretation constraints",
                    " | ".join(context.interpretation_constraints),
                ),
                self._text_observation(
                    case.patient_id,
                    "DDAVP reported response",
                    context.ddavp.reported_response,
                ),
            ]
        )
        return FHIRBundle.of(resources)

    @staticmethod
    def _patient_resource(case: PatientCase) -> FHIRResource:
        gender = None
        if case.clinical_context.sex:
            gender = "male" if case.clinical_context.sex == "男" else "female"
        return FHIRResource(
            resourceType="Patient",
            id=case.patient_id,
            gender=gender,
        )

    @staticmethod
    def _observed_value_observation(
        patient_id: str,
        display: str,
        value: ObservedValue,
    ) -> FHIRResource:
        components: list[dict[str, Any]] = []
        if value.unit:
            components.append({"code": {"text": "unit"}, "valueString": value.unit})
        if value.reference_range:
            components.append({"code": {"text": "reference_range"}, "valueString": value.reference_range})
        if value.missing_reason:
            components.append({"code": {"text": "missing_reason"}, "valueString": value.missing_reason.value})

        if value.observed and value.value is not None:
            return observation(
                observation_id=f"context-{display.lower().replace(' ', '-')}",
                patient_id=patient_id,
                display=display,
                value=value.value,
                unit=value.unit,
                components=components,
            )
        return observation(
            observation_id=f"context-{display.lower().replace(' ', '-')}",
            patient_id=patient_id,
            display=display,
            value="not_available" if value.missing_reason else "not_recorded",
            components=components,
        )

    @staticmethod
    def _text_observation(
        patient_id: str,
        display: str,
        value: str | None,
    ) -> FHIRResource:
        return observation(
            observation_id=f"context-{display.lower().replace(' ', '-')}",
            patient_id=patient_id,
            display=display,
            value=value or "not_recorded",
        )

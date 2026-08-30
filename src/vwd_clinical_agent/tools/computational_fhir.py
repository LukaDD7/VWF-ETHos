from __future__ import annotations

from typing import Any, Iterable

from ..schemas import EvidenceItem, VariantContext
from .fhir import (
    CodeableConcept,
    FHIRBundle,
    FHIRResource,
    Reference,
    observation,
    operation_outcome,
    provenance,
    stable_id,
    utc_now,
)


SOURCE_CATEGORIES = {
    "alphagenome": {
        "alphagenome_full_profile",
        "alphagenome_existing_panel",
    },
    "boltz": {
        "boltz2_type1_panel",
        "boltz2_functional_panel",
        "boltz_mechanism_classifier",
    },
    "md": {
        "md_type1_panel",
        "md_targeted_panel",
    },
}


class ComputationalPanelFHIRProvider:
    """Wrap local computational results in FHIR request/response resources."""

    name = "computational_panel_fhir"
    version = "v1"

    def collect(
        self,
        *,
        patient_id: str,
        variant: VariantContext,
        items: Iterable[EvidenceItem],
        statuses: dict[str, Any],
    ) -> FHIRBundle:
        item_list = list(items)
        resources: list[FHIRResource] = []
        for category, allowed_sources in SOURCE_CATEGORIES.items():
            request = self._service_request(patient_id, variant, category)
            resources.append(request)
            category_items = [item for item in item_list if item.source in allowed_sources]
            if category_items:
                for item in category_items:
                    result = self._evidence_observation(
                        patient_id,
                        variant,
                        category,
                        item,
                        request_id=request.id,
                    )
                    resources.append(result)
                    resources.append(
                        provenance(
                            resource_id=f"Observation/{result.id}",
                            tool_name=item.source,
                            endpoint=f"local://{category}",
                            request_payload={
                                "patient_id": patient_id,
                                "variant_id": variant.source_row_id,
                                "hgvs_c": variant.hgvs_c,
                                "hgvs_p": variant.hgvs_p,
                            },
                            response_payload=item.model_dump(mode="json"),
                            occurred=utc_now(),
                        )
                    )
            else:
                status = str(statuses.get(category, "not_available"))
                resources.append(
                    operation_outcome(
                        "warning",
                        "not-available",
                        f"{category} result unavailable for {variant.source_row_id}: {status}",
                    )
                )
        return FHIRBundle.of(resources)

    @staticmethod
    def _service_request(
        patient_id: str,
        variant: VariantContext,
        category: str,
    ) -> FHIRResource:
        return FHIRResource(
            resourceType="ServiceRequest",
            id=stable_id("computational-request", patient_id, variant.source_row_id, category),
            status="completed",
            intent="order",
            code=CodeableConcept(text=f"{category} computational analysis"),
            subject=Reference(reference=f"Patient/{patient_id}"),
            extension=[
                {
                    "url": "https://vwf-ethos.org/StructureDefinition/variant-context",
                    "valueString": variant.hgvs_c or variant.hgvs_p or variant.source_row_id,
                }
            ],
        )

    @staticmethod
    def _evidence_observation(
        patient_id: str,
        variant: VariantContext,
        category: str,
        item: EvidenceItem,
        *,
        request_id: str,
    ) -> FHIRResource:
        components = [
            {"code": {"text": "query"}, "valueString": item.query},
            {"code": {"text": "source_version"}, "valueString": item.source_version},
            {"code": {"text": "evidence_level"}, "valueString": item.evidence_level},
            {"code": {"text": "supports"}, "valueString": "|".join(item.supports)},
            {"code": {"text": "refutes"}, "valueString": "|".join(item.refutes)},
            {"code": {"text": "limitations"}, "valueString": "|".join(item.limitations)},
            {"code": {"text": "raw_excerpt_locator"}, "valueString": item.raw_excerpt_locator or ""},
        ]
        if item.confidence is not None:
            components.append(
                {"code": {"text": "confidence"}, "valueQuantity": {"value": item.confidence}}
            )
        return observation(
            observation_id=stable_id(
                "computational-result",
                patient_id,
                variant.source_row_id,
                category,
                item.source,
            ),
            patient_id=patient_id,
            display=f"{category} computational result",
            value=item.conclusion,
            based_on=[f"ServiceRequest/{request_id}"],
            components=components,
        )

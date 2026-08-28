from __future__ import annotations

import re
from typing import Any

from .base import BaseBiomedicalTool, ToolRequest
from .fhir import FHIRResource, observation, operation_outcome


# Mature VWF domain boundaries are approximate mature-protein intervals used only
# as a display/explanation aid, not as a structural predictor.
DOMAINS: list[tuple[int, int, str]] = [
    (1, 763, "D1-D2 propeptide"),
    (764, 1244, "D'-D3"),
    (1245, 1671, "A1"),
    (1672, 1875, "A2"),
    (1876, 2255, "A3"),
    (2256, 2571, "D4"),
    (2572, 2813, "C1-C6"),
    (2814, 2820, "CK"),
]


class VWFDomainAnnotator(BaseBiomedicalTool):
    name = "vwf_domain_annotator"
    version = "mature-protein-map-v1"
    endpoint = "local://vwf-domain-map"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        hgvs_p = str(request.parameters.get("hgvs_p", "")).strip()
        if not hgvs_p or hgvs_p.upper() in {"NA", "N/A"}:
            return [operation_outcome("information", "not-found", "No protein change was provided")], "not_found"
        match = re.search(r"p\.([A-Za-z]+)(\d+)([A-Za-z*+=]+)", hgvs_p)
        if not match:
            return [operation_outcome("information", "not-found", f"Could not parse protein position from {hgvs_p}")], "not_found"
        raw_ref, position_text, raw_alt = match.groups()
        position = int(position_text)
        domain = next((name for start, end, name in DOMAINS if start <= position <= end), "outside_annotated_domains")
        if "*" in raw_alt or raw_alt.lower().startswith("x"):
            variant_type = "nonsense"
        elif "=" in raw_alt or raw_alt == "=":
            variant_type = "synonymous"
        elif "fs" in raw_alt.lower() or "dup" in hgvs_p.lower() or "del" in hgvs_p.lower():
            variant_type = "frameshift_or_indel"
        elif raw_ref and raw_alt and len(raw_ref) == len(raw_alt) and raw_ref != raw_alt:
            variant_type = "missense"
        else:
            variant_type = "coding_effect_unknown"
        resource = observation(
            observation_id=f"vwf-domain-{position}-{raw_ref}-{raw_alt}",
            patient_id=request.patient_id,
            display="VWF mature-protein domain and variant class",
            value=f"p.{raw_ref}{position}{raw_alt}; {domain}; {variant_type}",
            based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
            components=[
                {"code": {"text": "protein_position"}, "valueQuantity": {"value": position}},
                {"code": {"text": "domain"}, "valueString": domain},
                {"code": {"text": "variant_class"}, "valueString": variant_type},
                {"code": {"text": "boundary_precision"}, "valueString": "approximate_display_only"},
            ],
        )
        return [resource, self.provenance_for(f"Observation/{resource.id}", request, {"hgvs_p": hgvs_p, "domain": domain, "variant_type": variant_type})], "success"

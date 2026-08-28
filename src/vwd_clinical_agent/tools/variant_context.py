from __future__ import annotations

import re
from typing import Any

from .base import BaseBiomedicalTool, ToolRequest
from .fhir import FHIRResource, observation, operation_outcome


# Canonical UniProt P04275 feature coordinates. HGVS p. positions and these
# coordinates both use prepro-VWF numbering; no mature-protein conversion is applied.
DOMAINS: list[tuple[int, int, str]] = [
    (1, 33, "VWFD1"),
    (295, 348, "TIL1"),
    (386, 560, "VWFD2"),
    (652, 707, "TIL2"),
    (776, 827, "TIL3"),
    (865, 1032, "VWFD3"),
    (1146, 1196, "TIL4"),
    (1277, 1453, "VWFA1"),
    (1498, 1665, "VWFA2"),
    (1691, 1871, "VWFA3"),
    (1948, 2124, "VWFD4"),
    (2255, 2328, "VWFC1"),
    (2429, 2495, "VWFC2"),
    (2580, 2645, "VWFC3"),
    (2724, 2812, "CTCK"),
]

DOMAIN_ORDER = [name for _, _, name in DOMAINS]


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
        domain = next((name for start, end, name in DOMAINS if start <= position <= end), "")
        if not domain:
            domain = _interdomain_label(position)
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
                {"code": {"text": "boundary_precision"}, "valueString": "canonical_feature_boundary"},
                {"code": {"text": "coordinate_system"}, "valueString": "UniProt_P04275_prepro_numbering"},
                {"code": {"text": "source"}, "valueString": "UniProt P04275 canonical features"},
            ],
        )
        return [resource, self.provenance_for(f"Observation/{resource.id}", request, {"hgvs_p": hgvs_p, "domain": domain, "variant_type": variant_type})], "success"


def _interdomain_label(position: int) -> str:
    previous = next((name for start, end, name in reversed(DOMAINS) if start <= position), None)
    following = next((name for start, end, name in DOMAINS if position <= start), None)
    if previous and following:
        return f"{previous}-{following} linker"
    return "outside_annotated_domains"

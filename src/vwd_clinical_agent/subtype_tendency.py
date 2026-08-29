from __future__ import annotations

from collections import defaultdict

from .schemas import EvidenceItem, SubtypeTendency
from .tools.fhir import FHIRBundle


def infer_subtype_tendency(
    *,
    fhir_bundle: FHIRBundle,
    ratio: float | None,
    candidate_subtypes: list[str],
    evidence_items: list[EvidenceItem] | None = None,
) -> list[SubtypeTendency]:
    scores: dict[str, float] = defaultdict(float)
    refs: dict[str, set[str]] = defaultdict(set)
    reasons: dict[str, list[str]] = defaultdict(list)
    seen_literature: set[str] = set()

    def add(label: str, amount: float, reason: str, ref: str = "") -> None:
        scores[label] += amount
        reasons[label].append(reason)
        if ref:
            refs[label].add(ref)

    candidates = set(candidate_subtypes)
    if "non_vwf_case_out_of_scope" in candidates:
        return []
    if "type_3_candidate" in candidates:
        add("type_3", 3.0, "Severe quantitative first-level deficit")
    if ratio is not None and ratio < 0.7:
        add("type_2_vwd", 2.0, "VWF:Act/VWF:Ag ratio below 0.70")
    elif "type_1_candidate_provisional" in candidates:
        add("type_1_or_low_vwf", 1.0, "Concordantly reduced VWF:Ag and VWF:Act")

    for item in evidence_items or []:
        if item.source != "boltz_mechanism_classifier":
            continue
        for label in item.supports:
            subtype = {
                "1": "type_1_or_low_vwf",
                "2A": "type_2A",
                "2B": "type_2B",
                "2M": "type_2M",
                "2N": "type_2N",
            }.get(label)
            if subtype:
                add(
                    subtype,
                    1.0,
                    "Boltz/structure mechanism classifier provides research-only directional support",
                    item.source_record_id,
                )

    for observation in fhir_bundle.resources("Observation"):
        name = (observation.get("code") or {}).get("text", "")
        ref = f"Observation/{observation['id']}"
        components = {
            (component.get("code") or {}).get("text"): _component_value(component)
            for component in observation.get("component", [])
        }
        value = observation.get("valueString") or str((observation.get("valueQuantity") or {}).get("value") or "")
        value_lower = value.casefold()
        literature_specific = str(components.get("variant_specific", "")).casefold() == "true"
        if name.startswith("Literature ") and name in seen_literature:
            continue
        if name.startswith("Literature "):
            seen_literature.add(name)

        if name == "VWF mature-protein domain and variant class":
            domain = str(components.get("domain", ""))
            if domain == "VWFA1":
                add("type_2B", 0.6, "Variant lies in the VWFA1/GPIb-binding domain", ref)
                add("type_2M", 0.4, "VWFA1 variants can also produce type 2M functional loss", ref)
            elif domain == "VWFA2":
                add("type_2A", 0.8, "Variant lies in VWFA2, where multimer processing defects are typical", ref)
            elif domain in {"VWFD3", "TIL4", "TIL3"}:
                add("type_2N", 0.4, "Variant lies in the D/D3 region, where FVIII-binding defects can occur", ref)

        if name in {"ClinVar clinical significance", "ClinGen expert variant classification"}:
            condition = str(components.get("condition", "")).casefold()
            for subtype in ("2a", "2b", "2m", "2n"):
                if f"type {subtype}" in condition or f"type{subtype}" in condition:
                    add(f"type_{subtype.upper()}", 2.0, f"Curated condition mentions type {subtype.upper()}", ref)
            if "pathogenic" in value_lower:
                add("pathogenic_variant", 1.0, "Curated database classification is pathogenic or likely pathogenic", ref)

        if name == "Literature RIPA evidence" and literature_specific:
            if any(term in value_lower for term in ("enhanced", "increased")):
                add("type_2B", 1.5, "Literature reports enhanced ristocetin-induced platelet aggregation", ref)
            elif any(term in value_lower for term in ("reduced", "poor", "abnormal")):
                add("type_2A", 0.6, "Literature reports reduced or abnormal RIPA", ref)
                add("type_2M", 0.6, "Reduced RIPA can also accompany type 2M functional loss", ref)

        if name == "Literature VWF multimer evidence" and literature_specific:
            if "normal" in value_lower:
                add(
                    "type_1_or_low_vwf",
                    0.4,
                    "Literature contains normal-multimer language; verify whether it refers to baseline or treatment response",
                    ref,
                )
            elif any(term in value_lower for term in ("abnormal", "reduced", "loss")):
                add("type_2A", 1.2, "Literature reports reduced or abnormal high-molecular-weight multimers", ref)

        if name == "Literature Secretion evidence" and literature_specific and any(
            term in value_lower for term in ("normal", "unaffected", "no significant")
        ):
            add("type_1_or_low_vwf", 0.4, "Literature reports normal secretion or no significant functional effect", ref)

        if name == "Literature Dominant-negative effect evidence" and literature_specific and any(
            term in value_lower for term in ("absent", "not observed")
        ):
            add("type_1_or_low_vwf", 0.3, "Literature reports no dominant-negative effect", ref)

    ranked = sorted(scores.items(), key=lambda item: item[1], reverse=True)[:3]
    return [
        SubtypeTendency(
            subtype_label=label,
            score=round(score, 3),
                confidence="moderate" if score >= 2.0 else "low",
            rationale="; ".join(dict.fromkeys(reasons[label])),
            evidence_refs=sorted(refs[label]),
        )
        for label, score in ranked
    ]


def _component_value(component: dict) -> object:
    for key in ("valueString", "valueBoolean", "valueQuantity", "valueInteger"):
        if key in component:
            return component[key]
    return None

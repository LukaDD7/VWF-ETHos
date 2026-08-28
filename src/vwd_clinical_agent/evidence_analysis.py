from __future__ import annotations

from typing import Any

from .schemas import EvidenceConflict
from .tools.fhir import FHIRBundle


def analyze_evidence_conflicts(
    *,
    fhir_bundle: FHIRBundle,
    second_level_status: str,
    candidate_subtypes: list[str],
) -> dict[str, Any]:
    conflicts: list[EvidenceConflict] = []
    missing: list[str] = []

    clinvar = _observation_by_name(fhir_bundle, "ClinVar clinical significance")
    clingen = _observation_by_name(fhir_bundle, "ClinGen expert variant classification")
    gnomad_exome = _observation_by_name(fhir_bundle, "gnomAD exome allele frequency")
    gnomad_genome = _observation_by_name(fhir_bundle, "gnomAD genome allele frequency")
    gnomad_aggregate = _observation_by_name(fhir_bundle, "gnomAD allele frequency")
    revel = _observation_by_name(fhir_bundle, "REVEL score")
    cadd = _observation_by_name(fhir_bundle, "CADD PHRED")
    alphamissense = _observation_by_name(fhir_bundle, "AlphaMissense pathogenicity")

    classification = _value_string(clinvar) or _value_string(clingen)
    classification_ref = _ref(clinvar) or _ref(clingen)
    frequencies = [
        value
        for value in (_value_float(gnomad_exome), _value_float(gnomad_genome), _value_float(gnomad_aggregate))
        if value is not None
    ]
    if classification and _is_pathogenic(classification) and frequencies and max(frequencies) >= 0.01:
        conflicts.append(
            EvidenceConflict(
                conflict_id="classification_vs_population_frequency",
                conflict_type="classification_vs_population_frequency",
                severity="critical",
                evidence_refs=[classification_ref, _ref(gnomad_exome), _ref(gnomad_genome), _ref(gnomad_aggregate)],
                description=f"Pathogenic classification conflicts with high population AF ({max(frequencies):.4g}).",
                explanation="A common allele is unlikely to be a highly penetrant dominant disease allele without unusual context.",
                recommended_action="expert_review",
            )
        )

    predictor_refs: list[str] = []
    predictor_votes: list[str] = []
    if _value_float(revel) is not None:
        predictor_refs.append(_ref(revel))
        predictor_votes.append("pathogenic" if _value_float(revel) >= 0.5 else "benign")
    if _value_float(cadd) is not None:
        predictor_refs.append(_ref(cadd))
        predictor_votes.append("pathogenic" if _value_float(cadd) >= 20 else "benign")
    if _value_float(alphamissense) is not None:
        predictor_refs.append(_ref(alphamissense))
        predictor_votes.append("pathogenic" if _value_float(alphamissense) >= 0.5 else "benign")
    if predictor_votes and len(set(predictor_votes)) > 1:
        conflicts.append(
            EvidenceConflict(
                conflict_id="internal_predictor_disagreement",
                conflict_type="classification_vs_prediction",
                severity="warning",
                evidence_refs=predictor_refs,
                description="Pathogenicity predictors disagree with each other.",
                explanation="REVEL, CADD, and AlphaMissense capture different signals; disagreement must remain explicit.",
                recommended_action="additional_evidence",
            )
        )

    clinvar_text = _value_string(clinvar)
    clingen_text = _value_string(clingen)
    if clinvar_text and clingen_text and _is_pathogenic(clinvar_text) != _is_pathogenic(clingen_text):
        conflicts.append(
            EvidenceConflict(
                conflict_id="clinvar_vs_clingen",
                conflict_type="database_vs_database",
                severity="critical",
                evidence_refs=[_ref(clinvar), _ref(clingen)],
                description=f"ClinVar ({clinvar_text}) and ClinGen ({clingen_text}) disagree.",
                explanation="Expert-panel curation and aggregate ClinVar submissions can diverge; expert review is required.",
                recommended_action="expert_review",
            )
        )

    if second_level_status == "not_available" and "type_2_candidate" in candidate_subtypes:
        missing.extend(["RIPA", "VWF_MULTIMER", "VWF_CB", "VWF_FVIIIB", "VWF_PP", "DDAVP_0_1_4H"])
        conflicts.append(
            EvidenceConflict(
                conflict_id="type2_discriminator_unavailable",
                conflict_type="missing_discriminator",
                severity="critical",
                evidence_refs=[],
                description="Type 2 candidate exists, but all second-level discriminators are unavailable.",
                explanation="The evidence cannot distinguish type 2 subtypes or exclude qualitative VWF defects.",
                recommended_action="accept_uncertainty",
            )
        )

    if not classification:
        missing.append("expert_variant_classification")
    if not frequencies:
        missing.append("population_frequency")
    if not predictor_refs:
        missing.append("pathogenicity_prediction")
    if not fhir_bundle.resources("DocumentReference"):
        missing.append("variant_specific_literature")

    acmg_hints = _acmg_hints(
        fhir_bundle=fhir_bundle,
        frequencies=frequencies,
    )
    if any("type_1" in subtype for subtype in candidate_subtypes):
        missing.extend(["ABO_blood_group_genotype", "VWF_clearance_context"])

    if not conflicts and missing:
        conflicts.append(
            EvidenceConflict(
                conflict_id="insufficient_evidence",
                conflict_type="insufficient_evidence",
                severity="warning",
                evidence_refs=[],
                description="No hard conflict, but required evidence is incomplete.",
                explanation="The remaining uncertainty is due to missing evidence rather than contradictory evidence.",
                recommended_action="additional_evidence",
            )
        )

    summary = f"{len(conflicts)} structured evidence issue(s); {len(missing)} missing evidence item(s)."
    return {
        "evidence_conflicts": conflicts,
        "evidence_missing": missing,
        "evidence_conflict_summary": summary,
        "acmg_evidence_hints": acmg_hints,
    }


def _observation_by_name(bundle: FHIRBundle, name: str) -> dict[str, Any] | None:
    for resource in bundle.resources("Observation"):
        if (resource.get("code") or {}).get("text") == name:
            return resource
    return None


def _value_string(resource: dict[str, Any] | None) -> str | None:
    return resource.get("valueString") if resource else None


def _value_float(resource: dict[str, Any] | None) -> float | None:
    if not resource:
        return None
    value = (resource.get("valueQuantity") or {}).get("value")
    return float(value) if value is not None else None


def _ref(resource: dict[str, Any] | None) -> str:
    return f"Observation/{resource['id']}" if resource else ""


def _is_pathogenic(value: str) -> bool:
    text = value.casefold()
    return "pathogenic" in text and "benign" not in text


def _acmg_hints(*, fhir_bundle: FHIRBundle, frequencies: list[float]) -> list[str]:
    hints: list[str] = []
    if frequencies and frequencies and max(frequencies) < 0.0001:
        hints.append("PM2_supporting_candidate: allele is absent or very rare in available gnomAD annotation")
    predictor_support = 0
    revel = _value_float(_observation_by_name(fhir_bundle, "REVEL score"))
    cadd = _value_float(_observation_by_name(fhir_bundle, "CADD PHRED"))
    alphamissense = _value_float(_observation_by_name(fhir_bundle, "AlphaMissense pathogenicity"))
    if revel is not None and revel >= 0.5:
        predictor_support += 1
    if cadd is not None and cadd >= 20:
        predictor_support += 1
    if alphamissense is not None and alphamissense >= 0.5:
        predictor_support += 1
    if predictor_support >= 2:
        hints.append("PP3_supporting_candidate: at least two computational predictors indicate deleterious missense effect")
    elif predictor_support == 1:
        hints.append("PP3_insufficient_alone: only one predictor currently supports deleterious effect")
    functional_normal = any(
        (_observation_by_name(fhir_bundle, name) or {}).get("valueString", "").casefold()
        in {"normal", "unaffected", "absent", "not observed"}
        for name in ["Secretion", "Dominant-negative effect", "VWF multimer"]
    )
    if functional_normal:
        hints.append("BS3_candidate_requires_review: full-text literature reports a normal functional or multimer phenotype")
    return hints

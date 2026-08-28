from __future__ import annotations

from pathlib import Path
import json
from typing import Any
from uuid import uuid4

from langgraph.graph import END, START, StateGraph

from .azure import LLMProvider
from .tools.fhir import FHIRBundle
from .tools.matrix import EvidenceToolMatrix
from .tools.second_level import SECOND_LEVEL_ACTIONS, SecondLevelFHIRStore
from .evidence_analysis import analyze_evidence_conflicts
from .schemas import (
    ClinicalAction,
    EvidenceItem,
    FinalClinicalPackage,
    PatientCase,
    SafetyFlag,
    TraceEvent,
    VWDWorkflowState,
    utc_now,
)


POLICY_DIR = Path(__file__).parent / "policies"


def _load_policy(name: str) -> dict[str, Any]:
    return json.loads((POLICY_DIR / name).read_text(encoding="utf-8"))


def _trace(state: VWDWorkflowState, node: str, **detail: Any) -> TraceEvent:
    case: PatientCase = state["case"]
    return TraceEvent(
        run_id=state["run_id"],
        case_id=case.patient_id,
        node=node,
        timestamp=utc_now(),
        detail=detail,
    )


def _lab_summary(case: PatientCase) -> dict[str, Any]:
    labs = case.first_level
    return {
        "VWF:Ag": labs.vwf_ag.model_dump(),
        "VWF:Act": labs.vwf_act.model_dump(),
        "FVIII:C": labs.fviii_c.model_dump(),
        "platelet_count": labs.platelet_count.model_dump(),
    }


def _fhir_evidence_summary(state: VWDWorkflowState) -> list[dict[str, Any]]:
    bundle = state.get("fhir_evidence_bundle")
    if not isinstance(bundle, FHIRBundle):
        return []
    evidence: list[dict[str, Any]] = []
    for resource in bundle.resources("Observation") + bundle.resources("DocumentReference"):
        code = resource.get("code") or {}
        value = resource.get("valueQuantity", {}).get("value") if resource.get("valueQuantity") else None
        value = value or resource.get("valueString") or resource.get("description")
        extensions = [
            extension
            for extension in resource.get("extension", [])
            if "pub-full-text" not in extension.get("url", "")
        ]
        full_text_extension = next(
            (
                extension
                for extension in resource.get("extension", [])
                if "pub-full-text" in extension.get("url", "")
            ),
            None,
        )
        evidence.append(
            {
                "resource_id": f"{resource['resourceType']}/{resource['id']}",
                "type": resource["resourceType"],
                "name": code.get("text") or resource.get("description"),
                "value": value,
                "extension": extensions,
                "full_text_available": bool(full_text_extension and full_text_extension.get("valueString")),
            }
        )
    return evidence


def load_case(state: VWDWorkflowState) -> dict[str, Any]:
    return {
        "run_id": state.get("run_id") or str(uuid4()),
        "status": "running",
        "trace": [_trace(state, "load_case", source_rows=state["case"].source_row_ids)],
    }


def validate_first_level(state: VWDWorkflowState) -> dict[str, Any]:
    case = state["case"]
    labs = case.first_level
    missing: list[str] = []
    if not labs.vwf_ag.observed:
        missing.append("VWF:Ag")
    if not labs.vwf_act.observed:
        missing.append("VWF:Act")
    if not labs.fviii_c.observed:
        missing.append("FVIII:C")
    ratio: float | None = None
    if labs.vwf_ag.observed and labs.vwf_act.observed and labs.vwf_ag.value not in (None, 0):
        ratio = labs.vwf_act.value / labs.vwf_ag.value
    elif labs.vwf_ag.value == 0 and labs.vwf_act.observed:
        missing.append("VWF:Act/VWF:Ag:denominator_zero")
    else:
        missing.append("VWF:Act/VWF:Ag:not_computable")
    return {
        "first_level": labs,
        "vwf_act_ag_ratio": ratio,
        "missing_critical_fields": missing,
        "trace": [_trace(state, "validate_first_level", missing=missing, ratio=ratio)],
    }


def pre_genetic_triage(state: VWDWorkflowState) -> dict[str, Any]:
    policy = _load_policy("first_level_v0.json")
    labs = state["case"].first_level
    ratio = state.get("vwf_act_ag_ratio")
    ag = labs.vwf_ag.value if labs.vwf_ag.observed else None
    act = labs.vwf_act.value if labs.vwf_act.observed else None
    fviii = labs.fviii_c.value if labs.fviii_c.observed else None
    platelets = labs.platelet_count.value if labs.platelet_count.observed else None

    hypotheses: list[str]
    route = "wait_genetics"
    if (
        ag is not None
        and act is not None
        and fviii is not None
        and ag <= policy["type3_candidate_max_vwf"]
        and act <= policy["type3_candidate_max_vwf"]
        and fviii <= policy["type3_candidate_max_fviii"]
    ):
        hypotheses = ["type_3_candidate"]
    elif ratio is not None and ratio < policy["act_ag_ratio_threshold"]:
        hypotheses = ["type_2_candidate"]
        if platelets is not None and platelets < policy["normal_platelet_min"]:
            hypotheses.append("platelet_type_vwd_candidate")
    elif (
        ratio is not None
        and ratio >= policy["act_ag_ratio_threshold"]
        and ag is not None
        and act is not None
        and fviii is not None
        and max(ag, act, fviii) <= policy["type1_candidate_max_lab"]
        and (platelets is None or platelets >= policy["normal_platelet_min"])
    ):
        hypotheses = ["type_1_candidate_provisional"]
    else:
        hypotheses = ["unresolved"]
    return {
        "pre_genetic_hypotheses": hypotheses,
        "pre_genetic_route": route,
        "trace": [_trace(state, "pre_genetic_triage", hypotheses=hypotheses, route=route)],
    }


def recommend_second_level(state: VWDWorkflowState) -> dict[str, Any]:
    policy = _load_policy("second_level_actions_v0.json")
    labs = state["case"].first_level
    missing = state.get("missing_critical_fields", [])
    candidates = set(state.get("candidate_subtypes", []))
    hypotheses = state.get("candidate_subtypes", state.get("pre_genetic_hypotheses", ["unresolved"]))
    codes: list[str] = []
    if "platelet_type_vwd_candidate" in candidates or (
        labs.platelet_count.observed and labs.platelet_count.value is not None and labs.platelet_count.value < 150
    ):
        codes.extend(["RIPA", "VWF_MULTIMER"])
    elif "type_2_candidate" in candidates:
        codes.extend(["VWF_MULTIMER", "VWF_CB"])
        if labs.fviii_c.observed and labs.fviii_c.value is not None and labs.fviii_c.value < 40:
            codes.append("VWF_FVIIIB")
    elif "type_1_candidate_provisional" in candidates:
        codes.extend(["VWF_PP", "DDAVP_0_1_4H"])
    elif "type_3_candidate" in candidates:
        codes.extend(["VWF_PP"])
    else:
        codes.extend(["VWF_MULTIMER", "VWF_CB"])
    codes = list(dict.fromkeys(codes))[:3]
    if any("VWF:Ag" in field or "VWF:Act" in field for field in missing):
        codes.insert(0, "COMPLETE_FIRST_LEVEL")
        codes = codes[:3]

    actions: list[ClinicalAction] = []
    for rank, code in enumerate(dict.fromkeys(codes), start=1):
        rationale = (
            f"Case {state['case'].patient_id}: first-level Ag={labs.vwf_ag.raw}, "
            f"Act={labs.vwf_act.raw}, ratio={state.get('vwf_act_ag_ratio')}; "
            f"missing={missing}; hypotheses={hypotheses}."
        )
        evidence_note = state.get("evidence_conflict_summary", "no structured conflict summary")
        actions.append(
            ClinicalAction(
                action_code=code,  # type: ignore[arg-type]
                rank=rank,
                status="recommended",
                availability="unknown",
                clinical_hypothesis=hypotheses,
                expected_discriminator=policy["discriminators"][code],
                rationale=f"{rationale} Evidence summary: {evidence_note}.",
                provenance=[f"policy:{policy['policy_id']}@{policy['version']}"],
            )
        )
    return {
        "recommended_actions": actions,
        "status": "waiting_second_level",
        "trace": [_trace(state, "recommend_second_level", actions=[a.action_code for a in actions])],
    }


def check_lab_availability(state: VWDWorkflowState) -> dict[str, Any]:
    provided_bundle = state.get("second_level_bundle")
    if isinstance(provided_bundle, FHIRBundle) and any(
        report.get("status") == "final" for report in provided_bundle.resources("DiagnosticReport")
    ):
        return {
            "second_level_status": "observed",
            "trace": [_trace(state, "check_lab_availability", environment="provided", observed_results=True)],
        }

    if state.get("second_level_environment") == "auto_unavailable":
        actions = [
            action.action_code
            for action in state.get("recommended_actions", [])
            if action.action_code in SECOND_LEVEL_ACTIONS
        ]
        store = SecondLevelFHIRStore.from_actions(
            state["case"].patient_id,
            actions,
            rationale="Automated retrospective environment feedback: requested second-level tests are unavailable.",
        )
        for action in actions:
            store.mark_unavailable(action, "environment feedback: not available")
        store.finalize(
            "All requested second-level tests were returned as unavailable by the automated retrospective environment; no values were imputed."
        )
        return {
            "second_level_bundle": store.bundle,
            "second_level_status": "not_available",
            "trace": [
                _trace(
                    state,
                    "check_lab_availability",
                    environment="auto_unavailable",
                    unavailable_actions=actions,
                    observed_results=0,
                    imputed_results=0,
                )
            ],
        }

    return {
        "second_level_status": "not_observed",
        "trace": [_trace(state, "check_lab_availability", availability="unknown", observed_results=False)],
    }


def ingest_second_level_results(state: VWDWorkflowState) -> dict[str, Any]:
    if state.get("second_level_status") == "not_available":
        return {
            "second_level_status": "not_available",
            "trace": [_trace(state, "ingest_second_level_results", fabricated_results=0, unavailable_tests=True)],
        }
    return {
        "second_level_status": "not_observed",
        "trace": [_trace(state, "ingest_second_level_results", fabricated_results=0)],
    }


def update_pre_genetic_assessment(state: VWDWorkflowState) -> dict[str, Any]:
    return {
        "status": "running",
        "trace": [_trace(state, "update_post_evidence_assessment", second_level=state.get("second_level_status"))],
    }


def wait_genetics(state: VWDWorkflowState) -> dict[str, Any]:
    has_variants = bool(state["case"].variants)
    return {
        "status": "waiting_genetics" if has_variants else "completed",
        "trace": [_trace(state, "wait_genetics", variants=len(state["case"].variants))],
    }


def normalize_variants(state: VWDWorkflowState) -> dict[str, Any]:
    variants = state["case"].variants
    return {
        "variants": variants,
        "trace": [
            _trace(
                state,
                "normalize_variants",
                variants=[variant.model_dump() for variant in variants],
                phase_inferred=False,
            )
        ],
    }


def plan_evidence_calls(state: VWDWorkflowState) -> dict[str, Any]:
    plan = [
        {
            "provider": "local_workbook",
            "query": {
                "patient_id": state["case"].patient_id,
                "variant": variant.hgvs_c or variant.source_row_id,
            },
        }
        for variant in state.get("variants", [])
    ]
    return {
        "evidence_plan": plan,
        "trace": [_trace(state, "plan_evidence_calls", planned_calls=len(plan))],
    }


def run_evidence_providers(state: VWDWorkflowState) -> dict[str, Any]:
    case = state["case"]
    items: list[EvidenceItem] = [
        EvidenceItem(
            source="local_workbook",
            source_record_id=variant.source_row_id,
            query=variant.hgvs_c or variant.source_row_id,
            retrieved_at=utc_now(),
            source_version="deidentified_v3",
            conclusion="Genetic report row was recorded; V0 performs no external variant interpretation.",
            confidence=None,
            limitations=["No gnomAD, ClinVar, HGMD, predictor, or guideline lookup was run in the offline profile."],
            raw_excerpt_locator=f"2.基因后!{variant.source_row_id}",
        )
        for variant in state.get("variants", [])
    ]
    provider_calls = [
        {
            "provider": "local_workbook",
            "version": "deidentified_v3",
            "case_id": case.patient_id,
            "status": "ok",
            "items_returned": len(items),
        }
    ]
    return {
        "evidence_items": items,
        "provider_calls": provider_calls,
        "trace": [_trace(state, "run_evidence_providers", evidence_items=len(items))],
    }


def run_fhir_evidence_tools(state: VWDWorkflowState, matrix: EvidenceToolMatrix) -> dict[str, Any]:
    second_level = state.get("second_level_bundle")
    if not isinstance(second_level, FHIRBundle):
        second_level = FHIRBundle()
    merged = FHIRBundle(type="collection")
    calls: list[dict[str, Any]] = []
    flags: list[SafetyFlag] = []
    for variant in state.get("variants", []):
        result = matrix.run(
            patient_id=state["case"].patient_id,
            variant_id=variant.source_row_id,
            hgvs_c=variant.hgvs_c or "",
            hgvs_p=variant.hgvs_p,
            second_level_bundle=second_level,
            allow_pre_second_level=True,
        )
        merged.entry.extend(result.bundle.entry)
        calls.extend(
            [
                {
                    "provider": call.tool,
                    "operation": call.operation,
                    "case_id": state["case"].patient_id,
                    "variant_id": variant.source_row_id,
                    "status": call.status,
                    "cache_hit": call.cache_hit,
                    "diagnostics": call.diagnostics,
                }
                for call in result.calls
            ]
        )
        for call in result.calls:
            if call.status in {"error", "not_found"}:
                flags.append(
                    SafetyFlag(
                        code=f"tool_{call.tool}_{call.status}",
                        severity="critical" if call.status == "error" else "warning",
                        message=f"{call.tool} returned {call.status}; this is not benign evidence.",
                    )
                )
    return {
        "fhir_evidence_bundle": merged,
        "provider_calls": calls,
        "safety_flags": flags,
        "trace": [_trace(state, "run_fhir_evidence_tools", authorized=True, calls=len(calls))],
    }


def integrate_patient_variant_evidence(state: VWDWorkflowState) -> dict[str, Any]:
    candidates = list(state.get("pre_genetic_hypotheses", ["unresolved"]))
    if len(state.get("variants", [])) > 1:
        candidates.append("multi_variant_unresolved")
    return {
        "candidate_subtypes": candidates,
        "trace": [_trace(state, "integrate_patient_variant_evidence", candidates=candidates)],
    }


def analyze_evidence(state: VWDWorkflowState) -> dict[str, Any]:
    bundle = state.get("fhir_evidence_bundle")
    if not isinstance(bundle, FHIRBundle):
        bundle = FHIRBundle()
    analysis = analyze_evidence_conflicts(
        fhir_bundle=bundle,
        second_level_status="not_available"
        if state.get("second_level_environment") == "auto_unavailable"
        else state.get("second_level_status", "not_observed"),
        candidate_subtypes=state.get("candidate_subtypes", ["unresolved"]),
    )
    return {
        **analysis,
        "acmg_evidence_hints": analysis.get("acmg_evidence_hints", []),
        "trace": [
            _trace(
                state,
                "analyze_evidence",
                conflict_count=len(analysis["evidence_conflicts"]),
                missing_count=len(analysis["evidence_missing"]),
            )
        ],
    }


def apply_type2_clingen_acmg(state: VWDWorkflowState) -> dict[str, Any]:
    flags = list(state.get("safety_flags", []))
    if "type_2_candidate" in state.get("candidate_subtypes", []):
        flags.append(
            SafetyFlag(
                code="type2_guideline_provider_not_enabled",
                severity="critical",
                message="ClinGen/ACMG type-2 rule provider is not enabled in V0 offline mode.",
            )
        )
    return {
        "safety_flags": flags,
        "trace": [_trace(state, "apply_type2_clingen_acmg", enabled=False)],
    }


def synthesize_opinion(state: VWDWorkflowState, llm_provider: LLMProvider) -> dict[str, Any]:
    system_prompt = (
        "You are a retrospective VWD research summarizer. Return only JSON with keys "
        "summary, abstention, expert_review_required, candidate_subtypes. Do not diagnose, "
        "do not invent laboratory results or evidence, and do not give treatment advice. "
        "Every factual statement must cite a FHIR resource ID from the provided evidence list; "
        "if no evidence supports a statement, omit it."
    )
    user_prompt = json.dumps(
        {
            "patient_id": state["case"].patient_id,
            "first_level": _lab_summary(state["case"]),
            "act_ag_ratio": state.get("vwf_act_ag_ratio"),
            "pre_genetic_hypotheses": state.get("pre_genetic_hypotheses", []),
            "recommended_actions": [a.action_code for a in state.get("recommended_actions", [])],
            "variants": [v.model_dump() for v in state.get("variants", [])],
            "second_level_status": state.get("second_level_status"),
            "evidence": _fhir_evidence_summary(state),
            "evidence_conflicts": [conflict.model_dump() for conflict in state.get("evidence_conflicts", [])],
            "evidence_missing": state.get("evidence_missing", []),
            "acmg_evidence_hints": state.get("acmg_evidence_hints", []),
        },
        ensure_ascii=False,
    )
    result = llm_provider.complete_json(system_prompt, user_prompt)
    summary = str(result.get("summary", ""))
    abstention = result.get("abstention", True)
    expert_review = result.get("expert_review_required", True)
    call = {
        "provider": llm_provider.name,
        "version": llm_provider.version,
        "deployment": getattr(llm_provider, "deployment", None),
        "case_id": state["case"].patient_id,
        "status": "ok",
        "json_mode_valid": isinstance(result, dict),
    }
    return {
        "llm_summary": summary,
        "llm_provider": llm_provider.name,
        "provider_calls": [call],
        "trace": [_trace(state, "synthesize_opinion", provider=llm_provider.name, abstention=abstention, expert_review=expert_review)],
    }


def safety_conflict_gate(state: VWDWorkflowState) -> dict[str, Any]:
    flags = list(state.get("safety_flags", []))
    if state.get("missing_critical_fields"):
        flags.append(
            SafetyFlag(
                code="critical_fields_missing",
                severity="critical",
                message=f"Missing critical fields: {', '.join(state['missing_critical_fields'])}.",
            )
        )
    flags.append(
        SafetyFlag(
            code="first_level_metadata_missing",
            severity="warning",
            message="Units, reference ranges, assay methods, and collection times are absent.",
        )
    )
    if state.get("second_level_status") == "not_observed":
        flags.append(
            SafetyFlag(
                code="second_level_not_observed",
                severity="critical",
                message="No real second-level result is present; retrospective mode must not synthesize one.",
            )
        )
    elif state.get("second_level_status") == "not_available":
        flags.append(
            SafetyFlag(
                code="second_level_unavailable",
                severity="warning",
                message="Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.",
            )
        )
    if len(state.get("variants", [])) > 1:
        flags.append(
            SafetyFlag(
                code="multi_variant_phase_unknown",
                severity="critical",
                message="Multiple variants are present and phase is unknown; compound mechanism cannot be inferred.",
            )
        )
    return {
        "safety_flags": flags,
        "trace": [_trace(state, "safety_conflict_gate", flags=[flag.code for flag in flags])],
    }


def expert_review(state: VWDWorkflowState) -> dict[str, Any]:
    return {
        "status": "expert_review",
        "trace": [_trace(state, "expert_review", reason="safety gate triggered")],
    }


def package_final_opinion(state: VWDWorkflowState) -> dict[str, Any]:
    flags = state.get("safety_flags", [])
    actions = list(state.get("recommended_actions", []))
    if not actions:
        actions = [
            ClinicalAction(
                action_code="EXPERT_REVIEW",
                rank=1,
                status="recommended",
                availability="unknown",
                clinical_hypothesis=state.get("candidate_subtypes", ["unresolved"]),
                expected_discriminator="Expert adjudication of incomplete or conflicting evidence.",
                rationale="No second-level recommendation was generated because first-level triage was not available.",
                provenance=["policy:vwd_safety_v0@2026-08-25"],
            )
        ]
    fhir_bundle = state.get("fhir_evidence_bundle")
    fhir_resource_ids: list[str] = []
    if isinstance(fhir_bundle, FHIRBundle):
        fhir_resource_ids = [
            f"{resource['resourceType']}/{resource['id']}"
            for resource in fhir_bundle.resources("Observation") + fhir_bundle.resources("DocumentReference")
        ]
    final = FinalClinicalPackage(
        patient_id=state["case"].patient_id,
        candidate_subtypes=state.get("candidate_subtypes", ["unresolved"]),
        confidence="low",
        opinion=state.get("llm_summary") or "Retrospective review abstains pending expert review.",
        abstention=True,
        expert_review_required=True,
        recommended_actions=actions,
        supporting_evidence=[item.source_record_id for item in state.get("evidence_items", [])] + fhir_resource_ids,
        contradicting_evidence=[
            f"{conflict.description} Evidence: {', '.join(conflict.evidence_refs)}"
            for conflict in state.get("evidence_conflicts", [])
        ],
        missing_information=state.get("missing_critical_fields", []) + state.get("evidence_missing", []),
        limitations=[flag.message for flag in flags]
        + [conflict.explanation for conflict in state.get("evidence_conflicts", [])]
        + [f"ACMG hint: {hint}" for hint in state.get("acmg_evidence_hints", [])],
        provenance=[
            "policy:vwd_first_level_v0@2026-08-25",
            "policy:vwd_second_level_actions_v0@2026-08-25",
            f"llm:{state.get('llm_provider', 'deterministic_policy')}",
        ],
    )
    return {
        "final_opinion": final,
        "status": "completed",
        "trace": [_trace(state, "package_final_opinion", abstention=True, expert_review_required=True)],
    }


def terminal_waiting(state: VWDWorkflowState) -> dict[str, Any]:
    return {
        "status": "waiting_genetics",
        "final_opinion": None,
        "trace": [_trace(state, "terminal_waiting", reason="genetic report not available")],
    }


def build_workflow(
    llm_provider: LLMProvider,
    evidence_tool_matrix: EvidenceToolMatrix | None = None,
    checkpointer: Any | None = None,
) -> Any:
    graph = StateGraph(VWDWorkflowState)
    graph.add_node("load_case", load_case)
    graph.add_node("validate_first_level", validate_first_level)
    graph.add_node("pre_genetic_triage", pre_genetic_triage)
    graph.add_node("recommend_second_level", recommend_second_level)
    graph.add_node("check_lab_availability", check_lab_availability)
    graph.add_node("ingest_second_level_results", ingest_second_level_results)
    graph.add_node("update_pre_genetic_assessment", update_pre_genetic_assessment)
    graph.add_node("wait_genetics", wait_genetics)
    graph.add_node("terminal_waiting", terminal_waiting)
    graph.add_node("normalize_variants", normalize_variants)
    graph.add_node("plan_evidence_calls", plan_evidence_calls)
    graph.add_node("run_evidence_providers", run_evidence_providers)
    if evidence_tool_matrix is not None:
        graph.add_node("run_fhir_evidence_tools", lambda state: run_fhir_evidence_tools(state, evidence_tool_matrix))
    graph.add_node("integrate_patient_variant_evidence", integrate_patient_variant_evidence)
    graph.add_node("analyze_evidence", analyze_evidence)
    graph.add_node("apply_type2_clingen_acmg", apply_type2_clingen_acmg)
    graph.add_node("synthesize_opinion", lambda state: synthesize_opinion(state, llm_provider))
    graph.add_node("safety_conflict_gate", safety_conflict_gate)
    graph.add_node("expert_review", expert_review)
    graph.add_node("package_final_opinion", package_final_opinion)

    graph.add_edge(START, "load_case")
    graph.add_edge("load_case", "validate_first_level")
    graph.add_edge("validate_first_level", "pre_genetic_triage")
    graph.add_edge("pre_genetic_triage", "wait_genetics")
    graph.add_conditional_edges(
        "wait_genetics",
        lambda state: "normalize_variants" if state["case"].variants else "terminal_waiting",
    )
    graph.add_edge("normalize_variants", "plan_evidence_calls")
    graph.add_edge("plan_evidence_calls", "run_evidence_providers")
    if evidence_tool_matrix is not None:
        graph.add_edge("run_evidence_providers", "run_fhir_evidence_tools")
        graph.add_edge("run_fhir_evidence_tools", "integrate_patient_variant_evidence")
    else:
        graph.add_edge("run_evidence_providers", "integrate_patient_variant_evidence")
    graph.add_edge("integrate_patient_variant_evidence", "analyze_evidence")
    graph.add_edge("analyze_evidence", "apply_type2_clingen_acmg")
    graph.add_edge("apply_type2_clingen_acmg", "recommend_second_level")
    graph.add_edge("recommend_second_level", "check_lab_availability")
    graph.add_edge("check_lab_availability", "ingest_second_level_results")
    graph.add_edge("ingest_second_level_results", "update_pre_genetic_assessment")
    graph.add_edge("update_pre_genetic_assessment", "synthesize_opinion")
    graph.add_edge("synthesize_opinion", "safety_conflict_gate")
    graph.add_conditional_edges(
        "safety_conflict_gate",
        lambda state: "expert_review" if any(flag.requires_expert_review for flag in state.get("safety_flags", [])) else "package_final_opinion",
    )
    graph.add_edge("expert_review", "package_final_opinion")
    graph.add_edge("package_final_opinion", END)
    graph.add_edge("terminal_waiting", END)
    return graph.compile(checkpointer=checkpointer)

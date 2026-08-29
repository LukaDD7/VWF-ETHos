from __future__ import annotations

from pathlib import Path
import json
import re
from typing import Any
from uuid import uuid4

from langgraph.graph import END, START, StateGraph

from .azure import LLMProvider
from .tools.fhir import FHIRBundle
from .tools.matrix import EvidenceToolMatrix
from .tools.computational_panel import LocalComputationalPanelProvider
from .tools.variant_context import DOMAINS
from .tools.second_level import SECOND_LEVEL_ACTIONS, SecondLevelFHIRStore
from .evidence_analysis import analyze_evidence_conflicts
from .subtype_tendency import infer_subtype_tendency
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


def _clinical_context_summary(case: PatientCase) -> dict[str, Any]:
    return case.clinical_context.model_dump(mode="json")


def _variant_domain(variant: Any) -> str:
    """Return the canonical VWF domain for a variant, if it can be parsed."""
    if not variant.hgvs_p:
        return ""
    match = re.search(r"p\.[A-Za-z]+(\d+)", variant.hgvs_p or "")
    if not match:
        return ""
    position = int(match.group(1))
    for start, end, name in DOMAINS:
        if start <= position <= end:
            return name
    nearest_name, nearest_distance = min(
        ((name, min(abs(position - start), abs(position - end))) for start, end, name in DOMAINS),
        key=lambda item: item[1],
    )
    if nearest_distance <= 20:
        return f"near {nearest_name}"
    return ""


def _mechanism_analysis(state: VWDWorkflowState) -> str:
    """Build a doctor-facing mechanism chain from labs, variants, and AI evidence."""
    case: PatientCase = state["case"]
    variants = state.get("variants", [])
    labs = case.first_level
    ratio = state.get("vwf_act_ag_ratio")
    context = case.clinical_context
    evidence_items = state.get("evidence_items", [])
    tendencies = state.get("subtype_tendencies", [])
    actions = state.get("recommended_actions", [])
    missing = state.get("evidence_missing", [])

    lines: list[str] = []
    lines.append("### 1. 临床与实验室表型")
    ag = labs.vwf_ag.value if labs.vwf_ag.observed else None
    act = labs.vwf_act.value if labs.vwf_act.observed else None
    fviii = labs.fviii_c.value if labs.fviii_c.observed else None
    platelets = labs.platelet_count.value if labs.platelet_count.observed else None
    lines.append(
        f"VWF:Ag={ag if ag is not None else '未提供'}, "
        f"VWF:Act={act if act is not None else '未提供'}, "
        f"FVIII:C={fviii if fviii is not None else '未提供'}, "
        f"血小板计数={platelets if platelets is not None else '未提供'}。"
    )
    if ratio is not None and ratio < 0.7:
        lines.append(
            f"VWF:Act/VWF:Ag 比值为 {ratio:.3f}，低于 0.70，提示存在不成比例的功能性 VWF 缺陷。"
        )
    else:
        lines.append("VWF:Act/VWF:Ag 比值未显示明显不成比例的功能下降。")
    clinical_text = " ".join(
        [
            context.symptoms or "",
            context.disease_course or "",
            " ".join(context.interpretation_constraints or []),
        ]
    )
    if platelets is not None and platelets < 150 or "血小板减少" in clinical_text:
        lines.append("临床背景中存在血小板减少线索，需与 2B/血小板型 VWD 鉴别。")

    lines.append("")
    lines.append("### 2. 变异与功能域")
    if not variants:
        lines.append("未提供 VWF 变异。")
    for variant in variants:
        domain = _variant_domain(variant)
        lines.append(
            f"- {variant.hgvs_c or variant.source_row_id} / {variant.hgvs_p or '无蛋白改变'}"
            f"（{domain or '未能映射到已知功能域'}）。"
        )

    lines.append("")
    lines.append("### 3. AI 机制证据")
    computational_sources = {
        "alphagenome_existing_panel",
        "alphagenome_full_profile",
        "boltz_mechanism_classifier",
        "boltz2_type1_panel",
        "boltz2_functional_panel",
        "md_type1_panel",
        "md_targeted_panel",
    }
    computational_items = [item for item in evidence_items if item.source in computational_sources]
    if not computational_items:
        lines.append("当前未获得可用的 AlphaGenome、Boltz 或 MD 机制证据。")
    for item in computational_items:
        lines.append(f"- {item.source}: {item.conclusion}")

    lines.append("")
    lines.append("### 4. 机制解释")
    domains = {_variant_domain(variant) for variant in variants}
    if any("VWFA1" in domain for domain in domains):
        lines.append(
            "该变异位于 A1 功能域，A1 参与 GPIb 结合并受 AIM 自抑制调控；"
            "该区域的错义变异可能通过改变自抑制界面、A1 暴露或表面电荷，影响血小板结合功能。"
        )
        if any(item.source == "md_targeted_panel" for item in computational_items):
            lines.append(
                "MD 结果提示 AIM-A1 接触动力学发生改变，可作为自抑制释放或 A1 功能面暴露的动态证据。"
            )
        if any(item.source == "boltz2_functional_panel" for item in computational_items):
            lines.append(
                "静态 Boltz 结果提示相关结构轴发生扰动，但静态置信度不等于结合自由能，需与 MD 和功能实验联合解释。"
            )
    if any("VWFA2" in domain for domain in domains):
        lines.append(
            "该变异位于 A2 功能域，A2 与多聚体加工和 ADAMTS13 易感性相关；"
            "该区域变异可能通过影响折叠稳定性或切割敏感性，导致高危分子量多聚体缺失。"
        )
    if any(domain in {"VWFD3", "TIL4", "TIL3"} or "VWFD3" in domain or "TIL4" in domain or "TIL3" in domain for domain in domains):
        lines.append(
            "该变异位于 D'/D3 相关区域，该区域参与 FVIII 结合并影响 2N 轴；"
            "若 FVIII:C 与 VWF:Ag 不成比例下降，应进一步评估 FVIII 结合功能。"
        )
    if any("splice" in item.conclusion.lower() for item in computational_items):
        lines.append(
            "AlphaGenome 提示剪接相关信号，需考虑变异通过转录/剪接异常导致 VWF 表达或分泌下降，"
            "而不一定直接改变蛋白结构。"
        )
    if ratio is not None and ratio < 0.7:
        lines.append(
            "实验室表型与功能性 VWF 缺陷方向一致；若同时存在高危分子量多聚体缺失，更支持 2A/2B 样机制。"
        )

    lines.append("")
    lines.append("### 5. 分型鉴别与不确定性")
    lines.append(
        "HGMD/ClinVar 的致病性评级本身不能直接给出 VWD 亚型方向；"
        "即使标注为致病或 uncertain，也必须结合实验室检查、出血表型和机制证据综合判断。"
    )
    if tendencies:
        tendency_text = "; ".join(
            f"{item.subtype_label}（{item.confidence}）" for item in tendencies
        )
        lines.append(f"当前倾向：{tendency_text}。")
    if "VWF_MULTIMER" in missing or "VWF_CB" in missing:
        lines.append(
            "VWF 多聚体分析和 VWF:CB/Ag 比值是区分 2A、2B 和 2M-A1 轴的关键检查；"
            "当前结果缺失，不能仅凭 AI 模型确定亚型。"
        )
    if actions:
        action_text = ", ".join(action.action_code for action in actions)
        lines.append(f"建议优先补充：{action_text}。")

    return "\n".join(lines)


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
    genes = {variant.gene.upper() for variant in state["case"].variants}
    if genes and "VWF" not in genes:
        hypotheses = ["non_vwf_case_out_of_scope"]
        return {
            "pre_genetic_hypotheses": hypotheses,
            "pre_genetic_route": "wait_genetics",
            "trace": [_trace(state, "pre_genetic_triage", hypotheses=hypotheses, route="wait_genetics")],
        }
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
        and max(ag, act) <= policy["type1_candidate_max_lab"]
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
    if "non_vwf_case_out_of_scope" in candidates:
        codes.append("EXPERT_REVIEW")
    elif "platelet_type_vwd_candidate" in candidates or (
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
        if not actions:
            return {
                "second_level_status": "not_available",
                "trace": [
                    _trace(
                        state,
                        "check_lab_availability",
                        environment="auto_unavailable",
                        unavailable_actions=[],
                        observed_results=0,
                        imputed_results=0,
                        reason="no in-scope second-level VWF action",
                    )
                ],
            }
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
    plan = []
    for variant in state.get("variants", []):
        plan.append(
            {
                "provider": "local_workbook",
                "query": {
                    "patient_id": state["case"].patient_id,
                    "variant": variant.hgvs_c or variant.source_row_id,
                },
            }
        )
        if variant.gene.upper() == "VWF":
            plan.extend(
                [
                    {"provider": "alphagenome_existing_panel", "query": variant.source_row_id},
                    {"provider": "boltz_mechanism_classifier", "query": variant.hgvs_p or variant.source_row_id},
                ]
            )
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


def run_computational_panels(
    state: VWDWorkflowState,
    provider: LocalComputationalPanelProvider,
) -> dict[str, Any]:
    items = list(state.get("evidence_items", []))
    statuses: list[dict[str, Any]] = []
    added = 0
    for variant in state.get("variants", []):
        variant_items, status = provider.collect(patient_id=state["case"].patient_id, variant=variant)
        items.extend(variant_items)
        statuses.append(status)
        added += len(variant_items)
    return {
        "evidence_items": items,
        "provider_calls": [
            {
                "provider": provider.name,
                "version": provider.version,
                "case_id": state["case"].patient_id,
                "status": "ok",
                "items_returned": added,
                "variant_statuses": statuses,
            }
        ],
        "trace": [_trace(state, "run_computational_panels", evidence_items_added=added, statuses=statuses)],
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
    if "non_vwf_case_out_of_scope" in candidates:
        return {
            "candidate_subtypes": candidates,
            "trace": [_trace(state, "integrate_patient_variant_evidence", candidates=candidates, out_of_scope=True)],
        }
    if len(state.get("variants", [])) > 1:
        candidates.append("multi_variant_unresolved")
    for item in state.get("evidence_items", []):
        if item.source != "boltz_mechanism_classifier":
            continue
        for label in item.supports:
            normalized = {
                "1": "type_1_candidate_provisional",
                "2A": "type_2A_candidate",
                "2B": "type_2B_candidate",
                "2M": "type_2M_candidate",
                "2N": "type_2N_candidate",
            }.get(label)
            if normalized and normalized not in candidates:
                candidates.append(normalized)
            if label.startswith("2") and "type_2_candidate" not in candidates:
                candidates.append("type_2_candidate")
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


def infer_subtype_tendency_node(state: VWDWorkflowState) -> dict[str, Any]:
    bundle = state.get("fhir_evidence_bundle")
    if not isinstance(bundle, FHIRBundle):
        bundle = FHIRBundle()
    tendencies = infer_subtype_tendency(
        fhir_bundle=bundle,
        ratio=state.get("vwf_act_ag_ratio"),
        candidate_subtypes=state.get("candidate_subtypes", ["unresolved"]),
        evidence_items=state.get("evidence_items", []),
    )
    return {
        "subtype_tendencies": tendencies,
        "trace": [_trace(state, "infer_subtype_tendency", tendencies=[item.model_dump() for item in tendencies])],
    }


def synthesize_opinion(state: VWDWorkflowState, llm_provider: LLMProvider) -> dict[str, Any]:
    mechanism_analysis = _mechanism_analysis(state)
    system_prompt = (
        "You are a retrospective VWD research summarizer. Return only JSON with keys "
        "summary, abstention, expert_review_required, candidate_subtypes. "
        "The summary must be a doctor-facing mechanism analysis that integrates clinical phenotype, "
        "variant domain, AI model evidence, and uncertainty; do not reduce it to a short checklist. "
        "Do not diagnose, do not invent laboratory results or evidence, and do not give treatment advice. "
        "Every factual statement must cite a FHIR resource ID from the provided evidence list; "
        "if no evidence supports a statement, omit it."
    )
    user_prompt = json.dumps(
        {
            "patient_id": state["case"].patient_id,
            "first_level": _lab_summary(state["case"]),
            "clinical_context": _clinical_context_summary(state["case"]),
            "act_ag_ratio": state.get("vwf_act_ag_ratio"),
            "pre_genetic_hypotheses": state.get("pre_genetic_hypotheses", []),
            "recommended_actions": [a.action_code for a in state.get("recommended_actions", [])],
            "variants": [v.model_dump() for v in state.get("variants", [])],
            "second_level_status": state.get("second_level_status"),
            "evidence": _fhir_evidence_summary(state),
            "computational_evidence": [
                item.model_dump(mode="json")
                for item in state.get("evidence_items", [])
                if item.source in {
                    "alphagenome_existing_panel",
                    "alphagenome_full_profile",
                    "boltz_mechanism_classifier",
                    "boltz2_type1_panel",
                    "md_type1_panel",
                    "boltz2_functional_panel",
                    "md_targeted_panel",
                }
            ],
            "evidence_conflicts": [conflict.model_dump() for conflict in state.get("evidence_conflicts", [])],
            "evidence_missing": state.get("evidence_missing", []),
            "acmg_evidence_hints": state.get("acmg_evidence_hints", []),
            "subtype_tendencies": [item.model_dump() for item in state.get("subtype_tendencies", [])],
            "mechanism_analysis": mechanism_analysis,
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
    labs = state["case"].first_level
    if any(not item.unit or not item.reference_range for item in (labs.vwf_ag, labs.vwf_act, labs.fviii_c)):
        flags.append(
            SafetyFlag(
                code="first_level_metadata_missing",
                severity="warning",
                message="Some units, reference ranges, assay methods, or collection times are absent.",
            )
        )
    non_vwf_genes = sorted({variant.gene for variant in state.get("variants", []) if variant.gene.upper() != "VWF"})
    if non_vwf_genes:
        flags.append(
            SafetyFlag(
                code="non_vwf_gene_out_of_scope",
                severity="critical",
                message=f"Reported gene(s) {', '.join(non_vwf_genes)} are outside the VWF mechanism model; no VWF-panel inference is allowed.",
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
    if len(state.get("variants", [])) == 1 and any(
        "复合" in (variant.zygosity or "") for variant in state.get("variants", [])
    ):
        flags.append(
            SafetyFlag(
                code="compound_heterozygous_report_incomplete",
                severity="critical",
                message="The report says compound heterozygous, but only one allele is represented; phase and the second allele must not be inferred.",
            )
        )
    if any(
        variant.gene.upper() == "VWF" and variant.alphagenome_request_status not in {None, "READY"}
        for variant in state.get("variants", [])
    ):
        flags.append(
            SafetyFlag(
                code="computational_request_not_modelable",
                severity="warning",
                message="At least one VWF variant lacks a model-ready normalized request; missing computational output is not negative evidence.",
            )
        )
    context = state["case"].clinical_context
    if (
        "hemophilia a" in (context.comorbidity or "").casefold()
        or any("fviii" in item.casefold() and "confound" in item.casefold() for item in context.interpretation_constraints)
    ):
        flags.append(
            SafetyFlag(
                code="fviii_confounded_by_hemophilia_a",
                severity="critical",
                message="Coexisting hemophilia A confounds FVIII:C; FVIII:C must not be used to support or refute VWD type 2N.",
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
        subtype_tendencies=state.get("subtype_tendencies", []),
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
    computational_panel_provider: LocalComputationalPanelProvider | None = None,
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
    if computational_panel_provider is not None:
        graph.add_node(
            "run_computational_panels",
            lambda state: run_computational_panels(state, computational_panel_provider),
        )
    if evidence_tool_matrix is not None:
        graph.add_node("run_fhir_evidence_tools", lambda state: run_fhir_evidence_tools(state, evidence_tool_matrix))
    graph.add_node("integrate_patient_variant_evidence", integrate_patient_variant_evidence)
    graph.add_node("analyze_evidence", analyze_evidence)
    graph.add_node("infer_subtype_tendency", infer_subtype_tendency_node)
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
    evidence_start = "run_evidence_providers"
    if computational_panel_provider is not None:
        graph.add_edge("run_evidence_providers", "run_computational_panels")
        evidence_start = "run_computational_panels"
    if evidence_tool_matrix is not None:
        graph.add_edge(evidence_start, "run_fhir_evidence_tools")
        graph.add_edge("run_fhir_evidence_tools", "integrate_patient_variant_evidence")
    else:
        graph.add_edge(evidence_start, "integrate_patient_variant_evidence")
    graph.add_edge("integrate_patient_variant_evidence", "analyze_evidence")
    graph.add_edge("analyze_evidence", "infer_subtype_tendency")
    graph.add_edge("infer_subtype_tendency", "apply_type2_clingen_acmg")
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

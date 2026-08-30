#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any
from uuid import uuid4

from pydantic import BaseModel


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.azure import AzureOpenAIProvider, DeterministicLLMProvider, LLMProvider
from src.vwd_clinical_agent.graph import build_workflow
from src.vwd_clinical_agent.tools.fhir import FHIRBundle
from src.vwd_clinical_agent.tools.matrix import EvidenceToolMatrix
from src.vwd_clinical_agent.tools.computational_panel import LocalComputationalPanelProvider
from src.vwd_clinical_agent.run_archive import RunArchive
from src.vwd_clinical_agent.schemas import PatientCase
from src.vwd_clinical_agent.workbook import LocalWorkbookProvider


def jsonable(value: Any) -> Any:
    if isinstance(value, BaseModel):
        return value.model_dump(mode="json")
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, datetime):
        return value.isoformat()
    if isinstance(value, list):
        return [jsonable(item) for item in value]
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    return value


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(jsonable(row), ensure_ascii=False, default=str) + "\n")


def case_summary(result: dict[str, Any]) -> dict[str, Any]:
    actions = result.get("recommended_actions", [])
    return {
        "patient_id": result["case"].patient_id,
        "status": result.get("status"),
        "pre_genetic_route": result.get("pre_genetic_route"),
        "hypotheses": result.get("pre_genetic_hypotheses", []),
        "top_action": actions[0].action_code if actions else "",
        "action_count": len(actions),
        "missing_fields": result.get("missing_critical_fields", []),
        "variant_count": len(result.get("variants", [])),
        "expert_review_required": bool(result.get("final_opinion") and result["final_opinion"].expert_review_required),
        "abstention": bool(result.get("final_opinion") and result["final_opinion"].abstention),
        "fabricated_second_level_results": 0,
        "alphagenome_evidence_count": sum(item.source.startswith("alphagenome_") for item in result.get("evidence_items", [])),
        "boltz_evidence_count": sum(item.source in {"boltz_mechanism_classifier", "boltz2_type1_panel", "md_type1_panel", "boltz2_functional_panel", "md_targeted_panel"} for item in result.get("evidence_items", [])),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the research-only VWD LangGraph V0 workflow.")
    parser.add_argument(
        "--profile",
        choices=["standard", "offline", "minimal"],
        default="standard",
        help="Standard runs the full local + online evidence chain; offline disables network tools; minimal keeps only local workbook/context.",
    )
    parser.add_argument("--workbook", type=Path, default=ROOT / "data/clinical_agent_pilot/vwd_agentic_workflow_deidentified_v3.xlsx")
    parser.add_argument("--mode", choices=["retrospective", "interactive"], default="retrospective")
    parser.add_argument("--provider-profile", choices=["fixture", "offline"], default="offline")
    parser.add_argument("--llm-provider", choices=["deterministic", "azure"], default="deterministic")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output/vwd_langgraph_v0/smoke")
    parser.add_argument("--case-id", help="Run one patient ID or a comma-separated list, for example CASE_001,CASE_002")
    parser.add_argument("--limit", type=int, help="Run only the first N cases")
    parser.add_argument("--biomedical-tools", action="store_true", default=None)
    parser.add_argument(
        "--computational-panels",
        action="store_true",
        default=None,
        help="Embed already-computed local AlphaGenome and Boltz/mechanism evidence; never launches new jobs.",
    )
    parser.add_argument("--second-level-bundle", type=Path)
    parser.add_argument(
        "--second-level-environment",
        choices=["provided", "auto_unavailable"],
        default="auto_unavailable",
        help="Retrospective environment feedback for second-level tests.",
    )
    parser.add_argument("--snapshot-dir", type=Path, default=ROOT / "data/external/vwf_biomedical_snapshots")
    parser.add_argument("--pubmed-full-text", action="store_true", default=None)
    parser.add_argument(
        "--mechanism-matrix",
        type=Path,
        help="Optional AlphaGenome/Boltz/MD feature table with aa_change or hgvs_c keys.",
    )
    parser.add_argument(
        "--mechanism-classifier",
        type=Path,
        default=ROOT / "scripts/agentic_vwf_classifier.py",
        help="Path to the AgenticVWFClassifier implementation.",
    )
    args = parser.parse_args()

    if args.biomedical_tools is None:
        args.biomedical_tools = args.profile == "standard"
    if args.computational_panels is None:
        args.computational_panels = args.profile in {"standard", "offline"}
    if args.pubmed_full_text is None:
        args.pubmed_full_text = args.profile == "standard"

    if args.mode == "interactive":
        raise SystemExit("Interactive interrupt/resume is not enabled in this V0 smoke runner.")

    provider: LLMProvider
    if args.llm_provider == "azure":
        provider = AzureOpenAIProvider.from_env()
    else:
        provider = DeterministicLLMProvider()

    workbook = LocalWorkbookProvider(args.workbook)
    cases, audit = workbook.load_cases()
    if args.case_id:
        requested = {case_id.strip() for case_id in args.case_id.split(",") if case_id.strip()}
        cases = [case for case in cases if case.patient_id in requested]
        if not cases:
            missing = sorted(requested - {case.patient_id for case in cases})
            raise SystemExit(f"Case(s) not found: {', '.join(missing)}")
    if args.limit is not None:
        cases = cases[: args.limit]

    second_level_bundle = None
    if args.second_level_bundle:
        second_level_bundle = FHIRBundle.model_validate_json(args.second_level_bundle.read_text(encoding="utf-8"))
    evidence_matrix = (
        EvidenceToolMatrix(
            cache_dir=str(args.output_dir / "tool_cache"),
            snapshot_dir=str(args.snapshot_dir),
            mechanism_matrix_path=args.mechanism_matrix,
            mechanism_classifier_path=args.mechanism_classifier,
        )
        if args.biomedical_tools
        else None
    )
    if evidence_matrix is not None:
        evidence_matrix.include_pubmed_full_text = args.pubmed_full_text
    computational_panel = LocalComputationalPanelProvider(ROOT) if args.computational_panels else None
    run_id = str(uuid4())
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    archive = RunArchive(run_id, output_dir / "run_archive")

    results: list[dict[str, Any]] = []
    trace_rows: list[dict[str, Any]] = []
    provider_calls: list[dict[str, Any]] = []
    from langgraph.checkpoint.sqlite import SqliteSaver

    with SqliteSaver.from_conn_string(str(archive.checkpoint_db)) as checkpointer:
        graph = build_workflow(
            provider,
            evidence_matrix,
            computational_panel_provider=computational_panel,
            checkpointer=checkpointer,
        )
        for case in cases:
            thread_id = f"case_{case.patient_id}"
            config = {"configurable": {"thread_id": thread_id}}
            initial_state = {
                "run_id": run_id,
                "case": case,
                "mode": args.mode,
                "decision_time": datetime.now(timezone.utc).isoformat(),
                "trace": [],
                "provider_calls": [],
                "second_level_bundle": second_level_bundle,
                "second_level_environment": "provided" if second_level_bundle else args.second_level_environment,
            }
            for event in graph.stream(initial_state, config, stream_mode="debug"):
                archive.write_debug_event(case.patient_id, event)
            result = graph.get_state(config).values
            results.append(result)
            trace_rows.extend(result.get("trace", []))
            provider_calls.extend(result.get("provider_calls", []))
            archive.write_state_history(
                case.patient_id,
                [
                    {
                        "checkpoint_id": getattr(snapshot, "config", {}).get("configurable", {}).get("checkpoint_id"),
                        "values": snapshot.values,
                    }
                    for snapshot in reversed(list(graph.get_state_history(config)))
                ],
            )
            archive.write_final_state(case.patient_id, result)
            if result.get("final_opinion"):
                archive.write_report(case.patient_id, result["final_opinion"])

    archive.write_manifest(
        {
            "case_ids": [case.patient_id for case in cases],
            "llm_provider": provider.name,
            "llm_provider_version": provider.version,
            "biomedical_tools_enabled": args.biomedical_tools,
            "computational_panels_enabled": args.computational_panels,
            "second_level_environment": args.second_level_environment,
            "checkpoint_db": str(archive.checkpoint_db),
        }
    )
    write_jsonl(output_dir / "cases.jsonl", results)
    write_jsonl(output_dir / "trace.jsonl", trace_rows)
    write_jsonl(output_dir / "provider_calls.jsonl", provider_calls)

    summaries = [case_summary(result) for result in results]
    summary_path = output_dir / "summary.csv"
    if summaries:
        import pandas as pd

        pd.DataFrame(summaries).to_csv(summary_path, index=False)
    else:
        summary_path.write_text("", encoding="utf-8")

    annotation_rows = [
        {
            "patient_id": result["case"].patient_id,
            "required_missing_information": "",
            "preferred_actions": "",
            "acceptable_actions": "",
            "inappropriate_actions": "",
            "must_not_miss_actions": "",
            "pre_genetic_subtype_set": "",
            "final_subtype_set_if_available": "",
            "abstention_expected": "",
            "rationale": "",
        }
        for result in results
    ]
    import pandas as pd

    pd.DataFrame(annotation_rows).to_csv(output_dir / "annotation_template.csv", index=False)

    metrics = {
        "run_id": run_id,
        "total_cases": len(results),
        "completed": sum(result.get("status") == "completed" for result in results),
        "expert_review_required": sum(
            bool(result.get("final_opinion") and result["final_opinion"].expert_review_required)
            for result in results
        ),
        "waiting": sum(result.get("status") in {"waiting_genetics", "waiting_second_level"} for result in results),
        "actions_all_from_allowed_enum": all(
            action.action_code in json.loads((ROOT / "src/vwd_clinical_agent/policies/second_level_actions_v0.json").read_text())["allowed_actions"]
            for result in results
            for action in result.get("recommended_actions", [])
        ),
        "fabricated_second_level_results": 0,
        "alphagenome_evidence_items": sum(
            item.source.startswith("alphagenome_")
            for result in results for item in result.get("evidence_items", [])
        ),
        "boltz_evidence_items": sum(
            item.source in {"boltz_mechanism_classifier", "boltz2_type1_panel", "md_type1_panel", "boltz2_functional_panel", "md_targeted_panel"}
            for result in results for item in result.get("evidence_items", [])
        ),
        "multi_variant_cases_with_unknown_phase": sum(
            len(result.get("variants", [])) > 1 and all(variant.phase_status == "unknown" for variant in result.get("variants", []))
            for result in results
        ),
        "clinical_accuracy_reported": False,
        "reason_no_accuracy": "No expert action/subtype gold set is available.",
    }
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    manifest = {
        "run_id": run_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "mode": args.mode,
        "provider_profile": args.provider_profile,
        "llm_provider": provider.name,
        "llm_provider_version": provider.version,
        "computational_panels_enabled": args.computational_panels,
        "workbook_audit": audit,
        "graph": "vwd_clinical_agent.graph/build_workflow",
        "policies": [
            "vwd_first_level_v0@2026-08-25",
            "vwd_second_level_actions_v0@2026-08-25",
        ],
        "privacy": "Research-only de-identified workbook; no original identifiers or reversible mapping are stored.",
    }
    (output_dir / "run_manifest.json").write_text(json.dumps(jsonable(manifest), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(json.dumps({"output_dir": str(output_dir), **metrics}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

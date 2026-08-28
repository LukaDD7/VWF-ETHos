#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.tools.fhir import FHIRBundle, operation_outcome
from src.vwd_clinical_agent.tools.matrix import EvidenceToolMatrix
from src.vwd_clinical_agent.workbook import LocalWorkbookProvider


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the gated FHIR-native VWF evidence-tool matrix.")
    parser.add_argument("--workbook", type=Path, default=ROOT / "data/clinical_agent_pilot/vwd_agentic_workflow_deidentified_v3.xlsx")
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--variant-index", type=int, default=1)
    parser.add_argument("--second-level-bundle", type=Path)
    parser.add_argument("--research-bypass", action="store_true")
    parser.add_argument("--snapshot-dir", type=Path, default=ROOT / "data/external/vwf_biomedical_snapshots")
    parser.add_argument("--cache-dir", type=Path, default=ROOT / "output/vwd_langgraph_v0/tool_cache")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    cases, _ = LocalWorkbookProvider(args.workbook).load_cases()
    case = next((item for item in cases if item.patient_id == args.case_id), None)
    if case is None:
        raise SystemExit(f"Case not found: {args.case_id}")
    variant = next((item for item in case.variants if item.variant_index == args.variant_index), None)
    if variant is None:
        raise SystemExit(f"Variant not found: {args.case_id} index {args.variant_index}")
    second_level = FHIRBundle()
    if args.second_level_bundle:
        second_level = FHIRBundle.model_validate_json(args.second_level_bundle.read_text(encoding="utf-8"))

    matrix = EvidenceToolMatrix(cache_dir=str(args.cache_dir))
    result = matrix.run(
        patient_id=case.patient_id,
        variant_id=variant.source_row_id,
        hgvs_c=variant.hgvs_c or "",
        hgvs_p=variant.hgvs_p,
        second_level_bundle=second_level,
        allow_research_bypass=args.research_bypass,
        local_parameters={"snapshot_dir": str(args.snapshot_dir)},
    )
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "evidence_bundle.json").write_text(result.bundle.model_dump_json(indent=2), encoding="utf-8")
    (output_dir / "tool_calls.json").write_text(json.dumps([call.model_dump() for call in result.calls], indent=2), encoding="utf-8")
    manifest = {
        "run_at": datetime.now(timezone.utc).isoformat(),
        "patient_id": case.patient_id,
        "variant": variant.model_dump(mode="json"),
        "authorized": result.authorized,
        "research_bypass": args.research_bypass,
        "normalized": result.normalized,
        "call_count": len(result.calls),
        "cache_dir": str(args.cache_dir),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    if not result.authorized:
        return 2
    errors = [call for call in result.calls if call.status == "error"]
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())

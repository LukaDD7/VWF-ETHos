#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def read_cases(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def split_actions(value: Any) -> set[str]:
    if pd.isna(value) or value is None:
        return set()
    return {part.strip() for part in str(value).split(";") if part.strip()}


def engineering_metrics(cases: list[dict[str, Any]]) -> dict[str, Any]:
    allowed = set(
        json.loads((ROOT / "src/vwd_clinical_agent/policies/second_level_actions_v0.json").read_text(encoding="utf-8"))["allowed_actions"]
    )
    action_pairs = [
        (case["case"]["patient_id"], action["action_code"])
        for case in cases
        for action in case.get("recommended_actions", [])
    ]
    invalid_actions = [pair for pair in action_pairs if pair[1] not in allowed]
    multi_variant = [case for case in cases if len(case.get("variants", [])) > 1]
    return {
        "total_cases": len(cases),
        "total_variant_rows": sum(len(case.get("variants", [])) for case in cases),
        "completed": sum(case.get("status") == "completed" for case in cases),
        "expert_review_final": sum(
            bool(case.get("final_opinion", {}) and case["final_opinion"].get("expert_review_required"))
            for case in cases
        ),
        "abstention_final": sum(
            bool(case.get("final_opinion", {}) and case["final_opinion"].get("abstention"))
            for case in cases
        ),
        "total_actions": len(action_pairs),
        "invalid_action_count": len(invalid_actions),
        "invalid_actions": invalid_actions,
        "fabricated_second_level_results": sum(case.get("second_level_status") == "observed" for case in cases),
        "multi_variant_cases": len(multi_variant),
        "multi_variant_phase_unknown": sum(
            all(variant.get("phase_status") == "unknown" for variant in case.get("variants", []))
            for case in multi_variant
        ),
        "clinical_accuracy_reported": False,
        "reason_no_accuracy": "No expert gold labels were supplied.",
    }


def gold_metrics(cases: list[dict[str, Any]], gold: pd.DataFrame) -> dict[str, Any]:
    gold_by_id = {str(row["patient_id"]).strip(): row for _, row in gold.iterrows()}
    matched = [case for case in cases if case["case"]["patient_id"] in gold_by_id]
    top1 = 0
    acceptable_recall_hits = 0
    acceptable_recall_total = 0
    inappropriate = 0
    must_miss_hits = 0
    must_miss_total = 0
    abstention_correct = 0
    for case in matched:
        row = gold_by_id[case["case"]["patient_id"]]
        actions = [action["action_code"] for action in case.get("recommended_actions", [])]
        preferred = split_actions(row.get("preferred_actions"))
        acceptable = split_actions(row.get("acceptable_actions"))
        inappropriate_set = split_actions(row.get("inappropriate_actions"))
        must_miss = split_actions(row.get("must_not_miss_actions"))
        if actions and actions[0] in preferred:
            top1 += 1
        acceptable_recall_total += len(acceptable)
        acceptable_recall_hits += len(set(actions) & acceptable)
        inappropriate += len(set(actions) & inappropriate_set)
        must_miss_total += len(must_miss)
        must_miss_hits += len(set(actions) & must_miss)
        expected_abstain = str(row.get("abstention_expected", "")).strip().casefold() in {"1", "true", "yes"}
        final = case.get("final_opinion") or {}
        if expected_abstain == bool(final.get("abstention", False)):
            abstention_correct += 1
    denominator = max(len(matched), 1)
    return {
        "gold_matched_cases": len(matched),
        "top1_preferred_action_accuracy": top1 / denominator,
        "topk_acceptable_action_recall": acceptable_recall_hits / max(acceptable_recall_total, 1),
        "inappropriate_action_rate": inappropriate / max(sum(len(case.get("recommended_actions", [])) for case in matched), 1),
        "must_not_miss_action_recall": must_miss_hits / max(must_miss_total, 1),
        "abstention_accuracy": abstention_correct / denominator,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate VWD LangGraph V0 outputs.")
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--gold", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = read_cases(args.predictions)
    metrics: dict[str, Any] = engineering_metrics(cases)
    if args.gold:
        gold = pd.read_csv(args.gold)
        required = {"patient_id", "preferred_actions", "acceptable_actions", "inappropriate_actions", "must_not_miss_actions"}
        missing = required - set(gold.columns)
        if missing:
            raise SystemExit(f"Gold file is missing columns: {sorted(missing)}")
        metrics.update(gold_metrics(cases, gold))
        metrics["clinical_accuracy_reported"] = True
        metrics.pop("reason_no_accuracy", None)
    output = json.dumps(metrics, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output, encoding="utf-8")
    else:
        sys.stdout.write(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

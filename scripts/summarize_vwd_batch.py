#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


OBSERVATION_NAMES = {
    "VWF mature-protein domain and variant class": "domain",
    "gnomAD exome allele frequency": "gnomad_exome_af",
    "gnomAD genome allele frequency": "gnomad_genome_af",
    "gnomAD allele frequency": "gnomad_aggregate_af",
    "ClinVar clinical significance": "clinvar",
    "ClinGen expert variant classification": "clingen",
    "REVEL score": "revel",
    "CADD PHRED": "cadd",
    "AlphaMissense pathogenicity": "alphamissense",
}


def load_states(archives: list[Path], wanted: set[str]) -> dict[str, dict[str, Any]]:
    states: dict[str, dict[str, Any]] = {}
    for archive in archives:
        if not archive.exists():
            continue
        paths = sorted(archive.glob("*/cases/*/final_state.json"), key=lambda path: path.stat().st_mtime)
        for path in paths:
            state = json.loads(path.read_text(encoding="utf-8"))
            case_id = state["case"]["patient_id"]
            if case_id in wanted:
                states[case_id] = state
    return states


def observation_values(state: dict[str, Any]) -> dict[str, list[str]]:
    values: dict[str, list[str]] = {name: [] for name in OBSERVATION_NAMES.values()}
    for entry in state.get("fhir_evidence_bundle", {}).get("entry", []):
        resource = entry.get("resource") or {}
        if resource.get("resourceType") != "Observation":
            continue
        label = (resource.get("code") or {}).get("text")
        if label not in OBSERVATION_NAMES:
            continue
        value = (resource.get("valueQuantity") or {}).get("value")
        value = value if value is not None else resource.get("valueString", "")
        values[OBSERVATION_NAMES[label]].append(str(value))
    return values


def summary_row(state: dict[str, Any]) -> dict[str, Any]:
    values = observation_values(state)
    variants = state.get("variants", [])
    return {
        "patient_id": state["case"]["patient_id"],
        "variants": "; ".join(
            f"{variant['hgvs_c']} ({variant['hgvs_p']}{' benign' if variant['benign_reported'] else ''})"
            for variant in variants
        ),
        "act_ag_ratio": state.get("vwf_act_ag_ratio"),
        "recommended_second_level": "; ".join(action["action_code"] for action in state.get("recommended_actions", [])),
        "subtype_tendencies": "; ".join(
            f"{item['subtype_label']}({item['score']},{item['confidence']})"
            for item in state.get("subtype_tendencies", [])
        ),
        "conflicts": "; ".join(conflict["conflict_id"] for conflict in state.get("evidence_conflicts", [])),
        "acmg_hints": "; ".join(state.get("acmg_evidence_hints", [])),
        **{key: "; ".join(items) for key, items in values.items()},
        "expert_review_required": bool((state.get("final_opinion") or {}).get("expert_review_required")),
    }


def markdown_table(rows: list[dict[str, Any]]) -> str:
    keys = [
        "patient_id", "variants", "act_ag_ratio", "domain", "subtype_tendencies", "gnomad_exome_af", "gnomad_genome_af",
        "clinvar", "clingen", "revel", "cadd", "alphamissense", "recommended_second_level",
        "conflicts", "acmg_hints", "expert_review_required",
    ]
    lines = ["| " + " | ".join(keys) + " |", "|" + "---|" * len(keys)]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(key, "")).replace("|", "\\|") for key in keys) + " |")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize VWD LangGraph final-state archives.")
    parser.add_argument("--run-archive", action="append", type=Path, required=True)
    parser.add_argument("--case-id", required=True, help="Comma-separated case IDs")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    wanted = {item.strip() for item in args.case_id.split(",") if item.strip()}
    states = load_states(args.run_archive, wanted)
    missing = wanted - set(states)
    if missing:
        raise SystemExit(f"Missing final states: {', '.join(sorted(missing))}")
    rows = [summary_row(states[case_id]) for case_id in sorted(wanted)]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "batch_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "batch_summary.md").write_text(markdown_table(rows), encoding="utf-8")
    with (args.output_dir / "batch_cases.jsonl").open("w", encoding="utf-8") as handle:
        for case_id in sorted(wanted):
            handle.write(json.dumps(states[case_id], ensure_ascii=False, default=str) + "\n")
    print(args.output_dir / "batch_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

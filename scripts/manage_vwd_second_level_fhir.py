#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.tools.second_level import SecondLevelFHIRStore, SecondLevelRecord


def main() -> int:
    parser = argparse.ArgumentParser(description="Manage FHIR-native VWD second-level orders and observations.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    create = subparsers.add_parser("create")
    create.add_argument("--patient-id", required=True)
    create.add_argument("--actions", nargs="+", required=True)
    create.add_argument("--rationale")
    create.add_argument("--output", type=Path, required=True)

    record = subparsers.add_parser("record")
    record.add_argument("--patient-id", required=True)
    record.add_argument("--bundle", type=Path, required=True)
    record.add_argument("--action", required=True)
    record.add_argument("--value", required=True)
    record.add_argument("--unit")
    record.add_argument("--observed-at")
    record.add_argument("--operator")
    record.add_argument("--method")

    unavailable = subparsers.add_parser("unavailable")
    unavailable.add_argument("--patient-id", required=True)
    unavailable.add_argument("--bundle", type=Path, required=True)
    unavailable.add_argument("--action", required=True)
    unavailable.add_argument("--reason", required=True)

    finalize = subparsers.add_parser("finalize")
    finalize.add_argument("--patient-id", required=True)
    finalize.add_argument("--bundle", type=Path, required=True)
    finalize.add_argument("--conclusion", required=True)

    args = parser.parse_args()
    if args.command == "create":
        store = SecondLevelFHIRStore.from_actions(args.patient_id, args.actions, args.rationale)
        store.save(args.output)
        print(args.output)
        return 0

    store = SecondLevelFHIRStore.load(args.patient_id, args.bundle)
    if args.command == "record":
        value: float | str | bool
        try:
            parsed = float(args.value)
            value = int(parsed) if parsed.is_integer() else parsed
        except ValueError:
            value = args.value.lower() in {"true", "false"} if args.value.lower() in {"true", "false"} else args.value
        store.record_observation(
            SecondLevelRecord(
                action=args.action,
                value=value,
                unit=args.unit,
                observed_at=args.observed_at,
                operator=args.operator,
                method=args.method,
            )
        )
    elif args.command == "unavailable":
        store.mark_unavailable(args.action, args.reason)
    elif args.command == "finalize":
        store.finalize(args.conclusion)
    store.save(args.bundle)
    print(args.bundle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

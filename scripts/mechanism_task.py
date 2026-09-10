#!/usr/bin/env python3
"""Create, validate, and ingest versioned VWD computational tasks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.mechanistic_tasks import (  # noqa: E402
    DEFAULT_REGISTRY_PATH,
    MechanismTaskResult,
    MechanismTaskStore,
    ProtocolRegistry,
    TaskProposal,
    create_task_request,
    result_template,
    validate_task_result,
)


def _source_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _write(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY_PATH)
    parser.add_argument("--store", type=Path, default=ROOT / "mechanism_tasks")
    sub = parser.add_subparsers(dest="command", required=True)

    list_parser = sub.add_parser("list", help="List mechanisms and submission-ready protocols")
    list_parser.add_argument("--domain")
    list_parser.add_argument("--variant-class", default="missense")

    submit_parser = sub.add_parser("submit", help="Validate a proposal and write an immutable request")
    submit_parser.add_argument("proposal", type=Path)
    submit_parser.add_argument("--task-id")
    submit_parser.add_argument("--source-commit", default=None)

    request_parser = sub.add_parser("validate-request", help="Validate one immutable request")
    request_parser.add_argument("task_id")

    template_parser = sub.add_parser("result-template", help="Create a server-side result template")
    template_parser.add_argument("task_id")
    template_parser.add_argument("--output", type=Path)

    validate_parser = sub.add_parser("validate-result", help="Validate a returned server result")
    validate_parser.add_argument("task_id")
    validate_parser.add_argument("result", type=Path)

    ingest_parser = sub.add_parser("ingest-result", help="Validate, store, and convert a result to FHIR")
    ingest_parser.add_argument("task_id")
    ingest_parser.add_argument("result", type=Path)

    export_parser = sub.add_parser("export-fhir", help="Export PlanDefinition and ActivityDefinitions")
    export_parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "protocols" / "vwd_mechanistic_v1" / "fhir",
    )

    args = parser.parse_args()
    registry = ProtocolRegistry.load(args.registry)
    store = MechanismTaskStore(args.store)

    if args.command == "list":
        candidates = (
            registry.candidate_mechanisms(args.domain, args.variant_class)
            if args.domain
            else list(registry.mechanisms)
        )
        print(
            json.dumps(
                {
                    "registry": registry.document.registry_id,
                    "version": registry.document.version,
                    "digest": registry.digest,
                    "candidate_mechanisms": candidates,
                    "submission_protocols": [
                        {
                            "protocol_id": item.protocol_id,
                            "version": item.version,
                            "digest": item.digest,
                            "mechanisms": item.mechanisms,
                        }
                        for item in registry.document.protocols
                        if item.submission_enabled and item.status == "active"
                    ],
                },
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0

    if args.command == "submit":
        proposal = TaskProposal.model_validate_json(args.proposal.read_text(encoding="utf-8"))
        request = create_task_request(
            proposal,
            registry,
            source_commit=args.source_commit or _source_commit(),
            task_id=args.task_id,
        )
        paths = store.submit(request, registry)
        print(json.dumps({"task_id": request.task_id, "request_digest": request.request_digest, "paths": [str(path) for path in paths]}, indent=2))
        return 0

    if args.command == "export-fhir":
        _write(args.output_dir / "PlanDefinition-vwd-mechanism-guided-acquisition-v1.json", registry.fhir_plan_definition())
        for protocol in registry.document.protocols:
            _write(
                args.output_dir / f"ActivityDefinition-{protocol.protocol_id}.json",
                registry.fhir_activity_definition(protocol.protocol_id),
            )
        print(json.dumps({"output_dir": str(args.output_dir), "registry_digest": registry.digest}, indent=2))
        return 0

    request = store.load_request(args.task_id, registry)
    if args.command == "validate-request":
        print(
            json.dumps(
                {
                    "valid": True,
                    "task_id": request.task_id,
                    "request_digest": request.request_digest,
                    "protocol_digest": request.acquisition.protocol_digest,
                },
                indent=2,
            )
        )
        return 0
    if args.command == "result-template":
        payload = result_template(request, registry)
        output = args.output or (store.result_dir(args.task_id) / "result.template.json")
        _write(output, payload)
        print(json.dumps({"output": str(output)}, indent=2))
        return 0

    result = MechanismTaskResult.model_validate_json(args.result.read_text(encoding="utf-8"))
    validate_task_result(result, request, registry)
    if args.command == "validate-result":
        print(json.dumps({"valid": True, "task_id": request.task_id, "status": result.status}, indent=2))
        return 0
    result_path = store.accept_result(result, request, registry)
    bundle_path = store.ingest_result(result, request, registry)
    print(json.dumps({"valid": True, "result": str(result_path), "fhir_bundle": str(bundle_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

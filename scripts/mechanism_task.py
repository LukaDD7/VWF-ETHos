#!/usr/bin/env python3
"""Create, validate, inspect, and ingest versioned VWD compute tasks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

from pydantic import ValidationError


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.mechanistic_tasks import (  # noqa: E402
    ArtifactIntegrityError,
    DEFAULT_REGISTRY_PATH,
    GitProvenanceError,
    MechanismTaskResult,
    MechanismTaskStore,
    ProtocolRegistry,
    TaskProposal,
    create_task_request,
    result_template,
    validate_request_git_provenance,
    validate_task_result,
)


TASK_ID_PATTERN = re.compile(r"^mtask-[0-9a-f]{32}$")


def _source_commit(repo_root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _write(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )


def _emit(payload: Any, *, stream: Any = sys.stdout) -> None:
    print(json.dumps(payload, ensure_ascii=False, indent=2), file=stream)


def _task_id(store: MechanismTaskStore, reference: str) -> str:
    if TASK_ID_PATTERN.fullmatch(reference):
        return reference
    path = Path(reference)
    if path.is_dir():
        path = path / "request.json"
    if not path.is_file():
        raise FileNotFoundError(
            f"Task reference must be a task ID, request.json, or request directory: {reference}"
        )
    payload = json.loads(path.read_text(encoding="utf-8"))
    task_id = str(payload.get("task_id", ""))
    if not TASK_ID_PATTERN.fullmatch(task_id):
        raise ValueError(f"Task reference does not contain a valid task_id: {path}")
    expected = store.request_dir(task_id) / "request.json"
    if expected.resolve() != path.resolve():
        raise ValueError(
            f"Request path is outside the configured store; expected {expected}, received {path}"
        )
    return task_id


def _verified_request(
    store: MechanismTaskStore,
    task_reference: str,
    *,
    repo_root: Path,
    registry_path: Path,
) -> tuple[Any, ProtocolRegistry]:
    task_id = _task_id(store, task_reference)
    request = store.read_request(task_id)
    registry = validate_request_git_provenance(
        request,
        repo_root,
        registry_path,
    )
    return request, registry


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY_PATH)
    parser.add_argument("--store", type=Path, default=ROOT / "mechanism_tasks")
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=ROOT,
        help="Approved root for repository-relative artifact paths.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    list_parser = sub.add_parser("list", help="List mechanisms and submission-ready protocols")
    list_parser.add_argument("--domain")
    list_parser.add_argument("--variant-class", default="missense")

    submit_parser = sub.add_parser("submit", help="Validate a proposal and write an immutable request")
    submit_parser.add_argument("proposal", type=Path)
    submit_parser.add_argument("--task-id")
    submit_parser.add_argument("--source-commit", default=None)

    request_parser = sub.add_parser("validate-request", help="Validate one immutable request")
    request_parser.add_argument("task")

    template_parser = sub.add_parser("result-template", help="Create a server-side result template")
    template_parser.add_argument("task")
    template_parser.add_argument("--output", type=Path)

    validate_parser = sub.add_parser("validate-result", help="Validate a returned server result")
    validate_parser.add_argument("task")
    validate_parser.add_argument("result", type=Path)

    ingest_parser = sub.add_parser("ingest-result", help="Validate, store, and convert a result to FHIR")
    ingest_parser.add_argument("task")
    ingest_parser.add_argument("result", type=Path)

    status_parser = sub.add_parser("status", help="Read and validate the local task lifecycle state")
    status_parser.add_argument("task")

    export_parser = sub.add_parser("export-fhir", help="Export PlanDefinition and ActivityDefinitions")
    export_parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "protocols" / "vwd_mechanistic_v1" / "fhir",
    )
    return parser


def run(args: argparse.Namespace) -> int:
    store = MechanismTaskStore(args.store)

    if args.command == "list":
        registry = ProtocolRegistry.load(args.registry)
        candidates = (
            registry.candidate_mechanisms(args.domain, args.variant_class)
            if args.domain
            else list(registry.mechanisms)
        )
        _emit(
            {
                "ok": True,
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
            }
        )
        return 0

    if args.command == "submit":
        registry = ProtocolRegistry.load(args.registry)
        proposal = TaskProposal.model_validate_json(args.proposal.read_text(encoding="utf-8"))
        request = create_task_request(
            proposal,
            registry,
            source_commit=args.source_commit or _source_commit(args.repo_root),
            task_id=args.task_id,
        )
        validate_request_git_provenance(request, args.repo_root, args.registry)
        paths = store.submit(request, registry)
        _emit(
            {
                "ok": True,
                "state": "awaiting_compute",
                "task_id": request.task_id,
                "request_digest": request.request_digest,
                "paths": [str(path) for path in paths],
            }
        )
        return 0

    if args.command == "export-fhir":
        registry = ProtocolRegistry.load(args.registry)
        _write(
            args.output_dir / "PlanDefinition-vwd-mechanism-guided-acquisition-v1.json",
            registry.fhir_plan_definition(),
        )
        for protocol in registry.document.protocols:
            _write(
                args.output_dir / f"ActivityDefinition-{protocol.protocol_id}.json",
                registry.fhir_activity_definition(protocol.protocol_id),
            )
        _emit(
            {
                "ok": True,
                "output_dir": str(args.output_dir),
                "registry_digest": registry.digest,
            }
        )
        return 0

    request, registry = _verified_request(
        store,
        args.task,
        repo_root=args.repo_root,
        registry_path=args.registry,
    )
    if args.command == "validate-request":
        _emit(
            {
                "ok": True,
                "valid": True,
                "task_id": request.task_id,
                "source_commit": request.source_commit,
                "request_digest": request.request_digest,
                "protocol_digest": request.acquisition.protocol_digest,
            }
        )
        return 0
    if args.command == "result-template":
        payload = result_template(request, registry)
        output = args.output or (store.result_dir(request.task_id) / "result.template.json")
        _write(output, payload)
        _emit({"ok": True, "output": str(output)})
        return 0
    if args.command == "status":
        result_path = store.result_dir(request.task_id) / "result.json"
        bundle_path = store.ingested_dir(request.task_id) / "bundle.fhir.json"
        state = "awaiting_compute"
        result_status = None
        if result_path.is_file():
            result = MechanismTaskResult.model_validate_json(result_path.read_text(encoding="utf-8"))
            validate_task_result(result, request, registry, artifact_root=args.artifact_root)
            state = "result_ready"
            result_status = result.status
        if bundle_path.is_file():
            json.loads(bundle_path.read_text(encoding="utf-8"))
            state = "ingested"
        _emit(
            {
                "ok": True,
                "task_id": request.task_id,
                "state": state,
                "result_status": result_status,
            }
        )
        return 0

    result = MechanismTaskResult.model_validate_json(args.result.read_text(encoding="utf-8"))
    validate_task_result(result, request, registry, artifact_root=args.artifact_root)
    if args.command == "validate-result":
        _emit(
            {
                "ok": True,
                "valid": True,
                "task_id": request.task_id,
                "status": result.status,
                "artifacts_verified": len(result.artifacts),
            }
        )
        return 0
    result_path = store.accept_result(
        result,
        request,
        registry,
        artifact_root=args.artifact_root,
    )
    bundle_path = store.ingest_result(
        result,
        request,
        registry,
        artifact_root=args.artifact_root,
    )
    _emit(
        {
            "ok": True,
            "state": "ingested",
            "result": str(result_path),
            "fhir_bundle": str(bundle_path),
        }
    )
    return 0


def main() -> int:
    args = _parser().parse_args()
    try:
        return run(args)
    except (
        ValidationError,
        GitProvenanceError,
        ArtifactIntegrityError,
        ValueError,
        FileNotFoundError,
        subprocess.SubprocessError,
    ) as exc:
        if isinstance(exc, GitProvenanceError):
            exit_code = 3
        elif isinstance(exc, ArtifactIntegrityError):
            exit_code = 4
        elif isinstance(exc, FileNotFoundError):
            exit_code = 5
        else:
            exit_code = 2
        _emit(
            {
                "ok": False,
                "error": {
                    "type": type(exc).__name__,
                    "message": str(exc),
                },
            },
            stream=sys.stderr,
        )
        return exit_code


if __name__ == "__main__":
    raise SystemExit(main())

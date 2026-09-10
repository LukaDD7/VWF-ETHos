#!/usr/bin/env python3
"""Publish, claim, run, and collect VWD mechanism tasks through Git branches."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

from pydantic import ValidationError


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.mechanistic_git import (  # noqa: E402
    GitMechanismTaskTransport,
    GitTransportError,
    ProtocolRunnerError,
    runner_command_from_config,
)
from src.vwd_clinical_agent.mechanistic_tasks import (  # noqa: E402
    ArtifactIntegrityError,
    GitProvenanceError,
    MechanismTaskStore,
    ProtocolRegistry,
    TaskProposal,
    create_task_request,
)


def _emit(payload: Any, *, stream: Any = sys.stdout) -> None:
    print(json.dumps(payload, ensure_ascii=False, indent=2), file=stream)


def _head(repo_root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--remote", default="origin")
    parser.add_argument(
        "--registry-relative",
        default="protocols/vwd_mechanistic_v1/registry.json",
    )
    parser.add_argument("--store-relative", default="mechanism_tasks")
    sub = parser.add_subparsers(dest="command", required=True)

    publish = sub.add_parser("publish", help="Create and push an immutable request branch")
    publish.add_argument("proposal", type=Path)
    publish.add_argument("--task-id")
    publish.add_argument("--source-commit")
    publish.add_argument("--no-push", action="store_true")

    claim = sub.add_parser("claim", help="Claim a request into an isolated server worktree")
    claim.add_argument("task_id")
    claim.add_argument("--worktree-root", type=Path, required=True)
    claim.add_argument("--worker-id", required=True)
    claim.add_argument("--no-push", action="store_true")

    run = sub.add_parser("run", help="Run the allowlisted protocol command and publish its result")
    run.add_argument("task_id")
    run.add_argument("--worktree", type=Path, required=True)
    run.add_argument("--runner-config", type=Path, required=True)
    run.add_argument("--no-push", action="store_true")

    return_result = sub.add_parser(
        "return",
        help="Validate and push a result written directly by the server agent",
    )
    return_result.add_argument("task_id")
    return_result.add_argument("--worktree", type=Path, required=True)
    return_result.add_argument("--no-push", action="store_true")

    collect = sub.add_parser("collect", help="Fetch, verify, and ingest a result branch")
    collect.add_argument("task_id")
    return parser


def run(args: argparse.Namespace) -> int:
    transport = GitMechanismTaskTransport(
        args.repo_root,
        remote=args.remote,
        registry_relative=args.registry_relative,
        store_relative=args.store_relative,
    )
    if args.command == "publish":
        registry_path = args.repo_root / args.registry_relative
        registry = ProtocolRegistry.load(registry_path)
        proposal = TaskProposal.model_validate_json(args.proposal.read_text(encoding="utf-8"))
        request = create_task_request(
            proposal,
            registry,
            source_commit=args.source_commit or _head(args.repo_root),
            task_id=args.task_id,
        )
        state = transport.publish_request(request, registry, push=not args.no_push)
    elif args.command == "claim":
        state = transport.claim_request(
            args.task_id,
            args.worktree_root,
            worker_id=args.worker_id,
            push=not args.no_push,
        )
    elif args.command == "run":
        store = MechanismTaskStore(args.worktree / args.store_relative)
        request = store.read_request(args.task_id)
        command = runner_command_from_config(args.runner_config, request.acquisition.protocol_id)
        state = transport.run_claimed(
            args.task_id,
            args.worktree,
            command,
            push=not args.no_push,
        )
    elif args.command == "return":
        state = transport.publish_claimed_result(
            args.task_id,
            args.worktree,
            push=not args.no_push,
        )
    else:
        state = transport.collect_result(args.task_id)
    _emit({"ok": True, **state.as_dict()})
    return 0


def main() -> int:
    args = _parser().parse_args()
    try:
        return run(args)
    except (
        ValidationError,
        GitProvenanceError,
        ArtifactIntegrityError,
        GitTransportError,
        ProtocolRunnerError,
        FileNotFoundError,
        FileExistsError,
        ValueError,
        subprocess.SubprocessError,
    ) as exc:
        _emit(
            {
                "ok": False,
                "error": {"type": type(exc).__name__, "message": str(exc)},
            },
            stream=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

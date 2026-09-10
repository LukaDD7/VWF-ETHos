#!/usr/bin/env python3
"""Run or resume the checkpointed asynchronous mechanism-compute subworkflow."""

from __future__ import annotations

import argparse
from contextlib import closing
import json
from pathlib import Path
import sqlite3
import sys
from typing import Any

from langgraph.checkpoint.sqlite import SqliteSaver
from langgraph.types import Command


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.mechanistic_workflow import (  # noqa: E402
    build_mechanism_compute_workflow,
)


def _emit(payload: Any, *, stream: Any = sys.stdout) -> None:
    print(json.dumps(payload, ensure_ascii=False, indent=2, default=str), file=stream)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument(
        "--registry",
        type=Path,
        default=ROOT / "protocols/vwd_mechanistic_v1/registry.json",
    )
    parser.add_argument("--store", type=Path, default=ROOT / "mechanism_tasks")
    parser.add_argument(
        "--checkpoint-db",
        type=Path,
        default=ROOT / "mechanism_tasks/checkpoints.sqlite",
    )
    parser.add_argument("--thread-id", required=True)
    parser.add_argument("--publish-to-git", action="store_true")
    parser.add_argument("--remote", default="origin")
    sub = parser.add_subparsers(dest="command", required=True)

    submit = sub.add_parser("submit", help="Submit and stop at the compute interrupt")
    submit.add_argument("proposal", type=Path)
    submit.add_argument("--task-id")
    submit.add_argument("--source-commit")

    sub.add_parser("resume", help="Resume after mechanism_git.py collect has validated the result")
    return parser


def _summary(state: dict[str, Any]) -> dict[str, Any]:
    interrupts = state.get("__interrupt__") or []
    return {
        "ok": True,
        "task_id": state.get("task_id"),
        "state": state.get("mechanism_compute_status"),
        "task_branch": state.get("task_branch"),
        "request_digest": state.get("request_digest"),
        "result_status": state.get("result_status"),
        "interrupts": [getattr(item, "value", item) for item in interrupts],
        "fhir_bundle_ready": bool(state.get("fhir_bundle")),
    }


def run(args: argparse.Namespace) -> int:
    args.checkpoint_db.parent.mkdir(parents=True, exist_ok=True)
    with closing(sqlite3.connect(str(args.checkpoint_db), check_same_thread=False)) as connection:
        checkpointer = SqliteSaver(connection)
        graph = build_mechanism_compute_workflow(
            repo_root=args.repo_root,
            registry_path=args.registry,
            store_root=args.store,
            publish_to_git=args.publish_to_git,
            remote=args.remote,
            checkpointer=checkpointer,
        )
        config = {"configurable": {"thread_id": args.thread_id}}
        if args.command == "submit":
            proposal = json.loads(args.proposal.read_text(encoding="utf-8"))
            initial: dict[str, Any] = {
                "proposal": proposal,
                "mechanism_compute_status": "proposed",
            }
            if args.task_id:
                initial["task_id"] = args.task_id
            if args.source_commit:
                initial["source_commit"] = args.source_commit
            state = graph.invoke(initial, config)
        else:
            saved = graph.get_state(config)
            if not saved.values or not saved.next:
                raise ValueError(
                    f"No interrupted mechanism task is available for thread {args.thread_id}"
                )
            state = graph.invoke(Command(resume={"result_collected": True}), config)
    _emit(_summary(state))
    return 0


def main() -> int:
    args = _parser().parse_args()
    try:
        return run(args)
    except Exception as exc:
        _emit(
            {"ok": False, "error": {"type": type(exc).__name__, "message": str(exc)}},
            stream=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

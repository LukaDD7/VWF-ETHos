"""Checkpointable LangGraph subworkflow for asynchronous mechanism compute."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
from typing import Any, Literal, TypedDict

from langgraph.graph import END, START, StateGraph
from langgraph.types import interrupt

from .mechanistic_git import GitMechanismTaskTransport
from .mechanistic_tasks import (
    DEFAULT_REGISTRY_PATH,
    MechanismTaskResult,
    MechanismTaskStore,
    ProtocolRegistry,
    TaskProposal,
    create_task_request,
    validate_request_git_provenance,
    validate_task_result,
)


class MechanismComputeState(TypedDict, total=False):
    proposal: dict[str, Any]
    source_commit: str
    task_id: str
    request: dict[str, Any]
    request_digest: str
    task_branch: str | None
    mechanism_compute_status: Literal["proposed", "awaiting_compute", "ingested"]
    result_status: Literal["completed", "failed", "inconclusive"]
    result: dict[str, Any]
    fhir_bundle: dict[str, Any]
    events: list[dict[str, Any]]


def _head(repo_root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def build_mechanism_compute_workflow(
    *,
    repo_root: str | Path,
    registry_path: str | Path = DEFAULT_REGISTRY_PATH,
    store_root: str | Path | None = None,
    publish_to_git: bool = False,
    remote: str = "origin",
    checkpointer: Any | None = None,
) -> Any:
    """Build a submit -> interrupt -> verified-ingest subgraph.

    The caller supplies a strictly validated ``TaskProposal`` payload.  When
    ``publish_to_git`` is enabled the submit node pushes a request branch.  The
    wait node interrupts until ``mechanism_git.py collect`` has materialized a
    validated result in the local append-only store, after which a resume
    ingests the FHIR bundle and returns it to the parent workflow state.
    """

    repo = Path(repo_root).resolve()
    registry_file = Path(registry_path).resolve()
    store_path = Path(store_root).resolve() if store_root is not None else repo / "mechanism_tasks"
    store = MechanismTaskStore(store_path)
    if publish_to_git and (
        registry_file != repo / "protocols/vwd_mechanistic_v1/registry.json"
        or store_path != repo / "mechanism_tasks"
    ):
        raise ValueError(
            "Git publication requires the canonical repository registry and mechanism_tasks store"
        )
    transport = GitMechanismTaskTransport(repo, remote=remote) if publish_to_git else None

    def submit(state: MechanismComputeState) -> dict[str, Any]:
        proposal = TaskProposal.model_validate(state["proposal"])
        registry = ProtocolRegistry.load(registry_file)
        request = create_task_request(
            proposal,
            registry,
            source_commit=state.get("source_commit") or _head(repo),
            task_id=state.get("task_id"),
        )
        validate_request_git_provenance(request, repo, registry_file)
        branch = None
        if transport is None:
            store.submit(request, registry)
        else:
            published = transport.publish_request(request, registry)
            branch = published.branch
        return {
            "task_id": request.task_id,
            "request": request.model_dump(mode="json"),
            "request_digest": request.request_digest,
            "task_branch": branch,
            "mechanism_compute_status": "awaiting_compute",
            "events": [
                {
                    "event": "mechanism_task_submitted",
                    "task_id": request.task_id,
                    "request_digest": request.request_digest,
                    "task_branch": branch,
                }
            ],
        }

    def await_and_ingest(state: MechanismComputeState) -> dict[str, Any]:
        task_id = state["task_id"]
        result_path = store.result_dir(task_id) / "result.json"
        if not result_path.is_file():
            interrupt(
                {
                    "kind": "awaiting_mechanism_compute",
                    "task_id": task_id,
                    "request_digest": state["request_digest"],
                    "task_branch": state.get("task_branch"),
                    "expected_result": str(result_path),
                    "resume_condition": (
                        "Collect and validate the mechanism-result branch, then resume this checkpoint."
                    ),
                }
            )
        if not result_path.is_file():
            raise FileNotFoundError(
                f"Task {task_id} was resumed before a validated result was collected"
            )
        request = store.read_request(task_id)
        registry = validate_request_git_provenance(request, repo, registry_file)
        result = MechanismTaskResult.model_validate_json(result_path.read_text(encoding="utf-8"))
        validate_task_result(result, request, registry, artifact_root=repo)
        bundle_path = store.ingest_result(
            result,
            request,
            registry,
            artifact_root=repo,
        )
        bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
        return {
            "mechanism_compute_status": "ingested",
            "result_status": result.status,
            "result": result.model_dump(mode="json"),
            "fhir_bundle": bundle,
            "events": [
                *state.get("events", []),
                {
                    "event": "mechanism_result_ingested",
                    "task_id": task_id,
                    "result_status": result.status,
                    "fhir_bundle": str(bundle_path),
                },
            ],
        }

    graph = StateGraph(MechanismComputeState)
    graph.add_node("submit_mechanism_task", submit)
    graph.add_node("await_mechanism_result", await_and_ingest)
    graph.add_edge(START, "submit_mechanism_task")
    graph.add_edge("submit_mechanism_task", "await_mechanism_result")
    graph.add_edge("await_mechanism_result", END)
    return graph.compile(checkpointer=checkpointer)

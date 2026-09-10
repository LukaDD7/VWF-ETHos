from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import pytest
from langgraph.checkpoint.memory import MemorySaver
from langgraph.types import Command
from pydantic import ValidationError

from src.vwd_clinical_agent.mechanistic_git import GitMechanismTaskTransport, GitTransportError
from src.vwd_clinical_agent.mechanistic_tasks import (
    ArtifactResult,
    ArtifactIntegrityError,
    DecisionContext,
    GitProvenanceError,
    MechanismAssessment,
    MechanismTaskResult,
    MechanismTaskStore,
    ProtocolRegistry,
    QCResult,
    TaskProposal,
    VariantSpec,
    create_task_request,
    result_template,
    result_to_fhir_bundle,
    task_to_fhir,
    validate_request_git_provenance,
    validate_task_result,
)
from src.vwd_clinical_agent.mechanistic_workflow import build_mechanism_compute_workflow


ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = ROOT / "protocols/vwd_mechanistic_v1/registry.json"


def git(repo: Path, *args: str) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def committed_registry_repo(tmp_path: Path, *, with_remote: bool = False) -> tuple[Path, ProtocolRegistry, str]:
    repo = tmp_path / "client"
    repo.mkdir()
    git(repo, "init", "-b", "main")
    registry_path = repo / "protocols/vwd_mechanistic_v1/registry.json"
    registry_path.parent.mkdir(parents=True)
    registry_path.write_text(REGISTRY_PATH.read_text(encoding="utf-8"), encoding="utf-8")
    git(repo, "add", "protocols/vwd_mechanistic_v1/registry.json")
    git(
        repo,
        "-c",
        "user.name=Test",
        "-c",
        "user.email=test@example.invalid",
        "commit",
        "-m",
        "registry",
    )
    if with_remote:
        remote = tmp_path / "remote.git"
        remote.mkdir()
        git(remote, "init", "--bare", "--initial-branch=main")
        git(repo, "remote", "add", "origin", str(remote))
        git(repo, "push", "-u", "origin", "main")
    return repo, ProtocolRegistry.load(registry_path), git(repo, "rev-parse", "HEAD")


def proposal() -> TaskProposal:
    return TaskProposal(
        patient_id="EVAL-TEST",
        variant=VariantSpec(
            hgvs_c="c.3946G>A",
            hgvs_p="p.Val1316Met",
            zygosity="heterozygous",
            domain="A1",
            variant_class="missense",
        ),
        mechanism_id="M6_2B_A1_GOF",
        protocol_id="VWF_A1_2B_GOF_LEGAN2023_V1",
        hypothesis_rationale=(
            "The current evidence leaves a competitive A1 gain-versus-loss-of-function question."
        ),
        competing_hypotheses=["M7_2M_A1_DISORDER", "M8_2M_A1_HYPERSTABLE", "M_UNKNOWN"],
        decision_context=DecisionContext(
            evidence_retrieval_complete=True,
            evidence_sufficiency="ambiguous",
            candidate_hypotheses=[
                "M1_SECRETION",
                "M6_2B_A1_GOF",
                "M7_2M_A1_DISORDER",
                "M8_2M_A1_HYPERSTABLE",
                "M_UNKNOWN",
            ],
            current_subtype_probabilities={
                "type_1": 0.02,
                "type_2A": 0.04,
                "type_2B": 0.46,
                "type_2M": 0.43,
                "type_2N": 0.02,
                "type_3": 0.01,
                "unresolved": 0.02,
            },
            current_evidence=["ClinGen/ClinVar/PubMed retrieval complete"],
            missing_evidence=["low-dose RIPA", "multimer analysis"],
            why_this_task_changes_decision=(
                "A calibrated activated-A1 pattern would separate the leading 2B and 2M hypotheses."
            ),
            expected_information_gain=0.42,
        ),
    )


def completed_result(request, registry: ProtocolRegistry) -> MechanismTaskResult:
    protocol = registry.protocols[request.acquisition.protocol_id]
    measurements = []
    for item in protocol.measurements:
        if item.comparison == "case_vs_matched_wt":
            case_value = 2.0
            reference_value = 1.0
            delta = 1.0
        else:
            case_value = {"activated": 0.7}
            reference_value = {"activated": 0.2}
            delta = {"activated": 0.5}
        measurements.append(
            {
                "measurement_id": item.measurement_id,
                "definition_id": item.definition_id,
                "status": "available",
                "unit": registry.measurements[item.definition_id].default_unit,
                "case_value": case_value,
                "reference_value": reference_value,
                "delta_case_minus_reference": delta,
                "comparison_kind": item.comparison,
                "replicate_values": [case_value],
                "provenance": {"roi": item.roi},
            }
        )
    return MechanismTaskResult(
        task_id=request.task_id,
        request_digest=request.request_digest,
        protocol_id=protocol.protocol_id,
        protocol_version=protocol.version,
        protocol_digest=protocol.digest,
        status="completed",
        completed_at="2026-09-10T00:00:00+00:00",
        software_versions={"gromacs": "2025.4"},
        runtime={"wall_seconds": 1200},
        qc=QCResult(
            passed=True,
            failures=[],
            replicate_consistency=0.9,
            equilibration_passed=True,
            structure_integrity_passed=True,
            out_of_distribution=False,
        ),
        measurements=measurements,
        benchmark={
            "similarity_to_positive_controls": 0.82,
            "similarity_to_negative_controls": 0.24,
            "calibration_set_id": "AIM-A1-v1",
            "out_of_distribution": False,
        },
        mechanism_assessment=MechanismAssessment(
            supports=["M6_2B_A1_GOF"],
            contradicts=[],
            indeterminate=[],
            confidence="moderate",
            reasons=["Composite measurement pattern matches the positive-control set."],
        ),
        artifacts=[],
        limitations=["Equilibrium MD does not measure platelet binding under flow."],
    )


def test_registry_integrity_and_open_set_candidates() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    candidates = registry.candidate_mechanisms("A1", "missense")
    assert "M_UNKNOWN" in candidates
    assert "M4_2A_ASSEMBLY" not in candidates
    assert len(registry.protocols) == 3
    assert all(protocol.digest.startswith("sha256:") for protocol in registry.protocols.values())


def test_checked_in_fhir_definitions_match_registry_export() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    fhir_dir = REGISTRY_PATH.parent / "fhir"
    plan = json.loads(
        (fhir_dir / "PlanDefinition-vwd-mechanism-guided-acquisition-v1.json").read_text(
            encoding="utf-8"
        )
    )
    assert plan == registry.fhir_plan_definition()
    for protocol_id in registry.protocols:
        activity = json.loads(
            (fhir_dir / f"ActivityDefinition-{protocol_id}.json").read_text(encoding="utf-8")
        )
        assert activity == registry.fhir_activity_definition(protocol_id)


def test_task_is_protocol_locked_and_fhir_requested() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    request = create_task_request(
        proposal(),
        registry,
        source_commit="a" * 40,
        task_id="mtask-" + "1" * 32,
    )
    protocol = registry.protocols[request.acquisition.protocol_id]
    assert request.acquisition.requested_measurements == protocol.required_measurement_ids
    task = task_to_fhir(request, registry)
    assert task["resourceType"] == "Task"
    assert task["status"] == "requested"
    assert protocol.protocol_id in task["instantiatesCanonical"]


def test_submission_requires_retrieval_unknown_and_information_gain() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    raw = proposal().model_dump(mode="json")
    raw["decision_context"]["evidence_retrieval_complete"] = False
    with pytest.raises(ValueError, match="retrieval must complete"):
        create_task_request(TaskProposal.model_validate(raw), registry, source_commit="a" * 40)

    raw = proposal().model_dump(mode="json")
    raw["decision_context"]["candidate_hypotheses"].remove("M_UNKNOWN")
    with pytest.raises(ValueError, match="M_UNKNOWN"):
        create_task_request(TaskProposal.model_validate(raw), registry, source_commit="a" * 40)

    raw = proposal().model_dump(mode="json")
    raw["decision_context"]["expected_information_gain"] = 0.01
    with pytest.raises(ValueError, match="information gain"):
        create_task_request(TaskProposal.model_validate(raw), registry, source_commit="a" * 40)


def test_noncompetitive_target_is_rejected() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    raw = proposal().model_dump(mode="json")
    raw["decision_context"]["current_subtype_probabilities"] = {
        "type_1": 0.01,
        "type_2A": 0.01,
        "type_2B": 0.05,
        "type_2M": 0.90,
        "type_2N": 0.01,
        "type_3": 0.01,
        "unresolved": 0.01,
    }
    with pytest.raises(ValueError, match="noncompetitive"):
        create_task_request(TaskProposal.model_validate(raw), registry, source_commit="a" * 40)


def test_completed_result_requires_exact_measurements_and_becomes_fhir_observation() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    request = create_task_request(
        proposal(),
        registry,
        source_commit="a" * 40,
        task_id="mtask-" + "2" * 32,
    )
    result = completed_result(request, registry)
    validate_task_result(result, request, registry)
    bundle = result_to_fhir_bundle(result, request, registry)
    resources = [entry["resource"] for entry in bundle["entry"]]
    assert [item["resourceType"] for item in resources] == ["Task", "Observation"]
    assert resources[0]["status"] == "completed"
    assert resources[1]["basedOn"] == [{"reference": f"Task/{request.task_id}"}]

    result.measurements = result.measurements[1:]
    with pytest.raises(ValueError, match="missing required measurements"):
        validate_task_result(result, request, registry)


def test_result_rejects_wrong_protocol_digest() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    request = create_task_request(proposal(), registry, source_commit="a" * 40)
    result = completed_result(request, registry)
    result.protocol_digest = "sha256:" + "0" * 64
    with pytest.raises(ValueError, match="exact task/protocol snapshot"):
        validate_task_result(result, request, registry)


def test_untouched_result_template_is_intentionally_invalid() -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    request = create_task_request(proposal(), registry, source_commit="a" * 40)
    with pytest.raises(ValidationError, match="template_only"):
        MechanismTaskResult.model_validate(result_template(request, registry))


def test_git_provenance_requires_existing_commit_and_exact_registry(tmp_path: Path) -> None:
    repo, registry, commit = committed_registry_repo(tmp_path)
    request = create_task_request(
        proposal(),
        registry,
        source_commit=commit,
        task_id="mtask-" + "3" * 32,
    )
    snapshot = validate_request_git_provenance(
        request,
        repo,
        repo / "protocols/vwd_mechanistic_v1/registry.json",
    )
    assert snapshot.digest == registry.digest

    request.source_commit = "aaaaaaa"
    request.request_digest = "sha256:" + "0" * 64
    from src.vwd_clinical_agent.mechanistic_tasks import request_digest_for

    request.request_digest = request_digest_for(request)
    with pytest.raises(GitProvenanceError, match="cat-file"):
        validate_request_git_provenance(
            request,
            repo,
            repo / "protocols/vwd_mechanistic_v1/registry.json",
        )

    registry_path = repo / "protocols/vwd_mechanistic_v1/registry.json"
    dirty_payload = json.loads(registry_path.read_text(encoding="utf-8"))
    dirty_payload["version"] += ".dirty"
    registry_path.write_text(json.dumps(dirty_payload), encoding="utf-8")
    dirty_registry = ProtocolRegistry.load(registry_path)
    dirty_request = create_task_request(
        proposal(),
        dirty_registry,
        source_commit=commit,
        task_id="mtask-" + "7" * 32,
    )
    with pytest.raises(GitProvenanceError, match="registry digest"):
        validate_request_git_provenance(dirty_request, repo, registry_path)


def test_artifacts_are_confined_and_rehashed(tmp_path: Path) -> None:
    registry = ProtocolRegistry.load(REGISTRY_PATH)
    request = create_task_request(proposal(), registry, source_commit="a" * 40)
    result = completed_result(request, registry)
    artifact = tmp_path / "mechanism_tasks/artifacts" / request.task_id / "summary.json"
    artifact.parent.mkdir(parents=True)
    artifact.write_bytes(b'{"verified":true}\n')
    result.artifacts = [
        ArtifactResult(
            path=artifact.relative_to(tmp_path).as_posix(),
            sha256=sha256(artifact.read_bytes()).hexdigest(),
            media_type="application/json",
            role="summary",
        )
    ]
    validate_task_result(result, request, registry, artifact_root=tmp_path)

    artifact.write_bytes(b"tampered\n")
    with pytest.raises(ArtifactIntegrityError, match="SHA256 mismatch"):
        validate_task_result(result, request, registry, artifact_root=tmp_path)

    result.artifacts[0].path = "../outside.json"
    with pytest.raises(ArtifactIntegrityError, match="relative POSIX"):
        validate_task_result(result, request, registry, artifact_root=tmp_path)

    outside = tmp_path / "outside.json"
    outside.write_bytes(b"outside\n")
    link = tmp_path / "mechanism_tasks/artifacts" / request.task_id / "linked.json"
    link.symlink_to(outside)
    result.artifacts[0].path = link.relative_to(tmp_path).as_posix()
    result.artifacts[0].sha256 = sha256(outside.read_bytes()).hexdigest()
    with pytest.raises(ArtifactIntegrityError, match="symlinks"):
        validate_task_result(result, request, registry, artifact_root=tmp_path)


def test_cli_returns_structured_contract_error_without_traceback(tmp_path: Path) -> None:
    repo, registry, commit = committed_registry_repo(tmp_path)
    task_id = "mtask-" + "4" * 32
    request = create_task_request(proposal(), registry, source_commit=commit, task_id=task_id)
    store = MechanismTaskStore(repo / "mechanism_tasks")
    store.submit(request, registry)
    template_path = repo / "result.template.json"
    template_path.write_text(json.dumps(result_template(request, registry)), encoding="utf-8")
    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/mechanism_task.py"),
            "--repo-root",
            str(repo),
            "--registry",
            str(repo / "protocols/vwd_mechanistic_v1/registry.json"),
            "--store",
            str(store.root),
            "validate-result",
            str(store.request_dir(task_id) / "request.json"),
            str(template_path),
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 2
    error = json.loads(completed.stderr)
    assert error["ok"] is False
    assert error["error"]["type"] == "ValidationError"
    assert "Traceback" not in completed.stderr


def test_git_publish_claim_run_collect_roundtrip(tmp_path: Path) -> None:
    client, registry, commit = committed_registry_repo(tmp_path, with_remote=True)
    task_id = "mtask-" + "5" * 32
    request = create_task_request(proposal(), registry, source_commit=commit, task_id=task_id)
    client_transport = GitMechanismTaskTransport(client)
    published = client_transport.publish_request(request, registry)
    assert published.state == "awaiting_compute"

    server = tmp_path / "server"
    subprocess.run(
        ["git", "clone", str(tmp_path / "remote.git"), str(server)],
        check=True,
        capture_output=True,
        text=True,
    )
    server_transport = GitMechanismTaskTransport(server)
    claimed = server_transport.claim_request(
        task_id,
        tmp_path / "server-worktrees",
        worker_id="gpu-test-01",
    )
    assert claimed.state == "claimed"
    with pytest.raises(GitTransportError, match="already claimed"):
        server_transport.claim_request(
            task_id,
            tmp_path / "second-worker",
            worker_id="gpu-test-02",
        )

    runner = tmp_path / "fixture_runner.py"
    runner.write_text(
        "import json, sys\n"
        "template, result = map(__import__('pathlib').Path, sys.argv[1:3])\n"
        "payload = json.loads(template.read_text())\n"
        "payload.pop('template_only')\n"
        "payload['status'] = 'completed'\n"
        "payload['software_versions'] = {'runner': 'fixture-v1'}\n"
        "payload['runtime'] = {'wall_seconds': 1}\n"
        "payload['qc'] = {'passed': True, 'failures': [], 'replicate_consistency': 1.0, 'equilibration_passed': True, 'structure_integrity_passed': True, 'out_of_distribution': False}\n"
        "payload['benchmark'] = {'similarity_to_positive_controls': 0.8, 'similarity_to_negative_controls': 0.2, 'calibration_set_id': 'fixture', 'out_of_distribution': False}\n"
        "for item in payload['measurements']:\n"
        "    item.update(status='available', case_value=2.0, reference_value=1.0, delta_case_minus_reference=1.0, replicate_values=[2.0])\n"
        "payload['mechanism_assessment'] = {'supports': ['M6_2B_A1_GOF'], 'contradicts': [], 'indeterminate': [], 'confidence': 'moderate', 'reasons': ['fixture contract pattern']}\n"
        "payload['limitations'] = ['Fixture result; not scientific evidence.']\n"
        "result.parent.mkdir(parents=True, exist_ok=True)\n"
        "result.write_text(json.dumps(payload, sort_keys=True, indent=2) + '\\n')\n",
        encoding="utf-8",
    )
    returned = server_transport.run_claimed(
        task_id,
        claimed.worktree or "",
        [sys.executable, str(runner), "{template}", "{result}"],
    )
    assert returned.state == "result_ready"
    assert returned.result_status == "completed"

    collected = client_transport.collect_result(task_id)
    assert collected.state == "ingested"
    assert (client / f"mechanism_tasks/results/{task_id}/result.json").is_file()
    assert (client / f"mechanism_tasks/ingested/{task_id}/bundle.fhir.json").is_file()


def test_git_transport_rejects_task_id_path_escape(tmp_path: Path) -> None:
    client, _, _ = committed_registry_repo(tmp_path, with_remote=True)
    transport = GitMechanismTaskTransport(client)
    with pytest.raises(GitTransportError, match="Invalid mechanism task ID"):
        transport.claim_request(
            "../outside",
            tmp_path / "server-worktrees",
            worker_id="gpu-test-01",
        )


def test_langgraph_subworkflow_interrupts_and_resumes_after_result(tmp_path: Path) -> None:
    repo, registry, commit = committed_registry_repo(tmp_path)
    store = MechanismTaskStore(repo / "mechanism_tasks")
    task_id = "mtask-" + "6" * 32
    graph = build_mechanism_compute_workflow(
        repo_root=repo,
        registry_path=repo / "protocols/vwd_mechanistic_v1/registry.json",
        store_root=store.root,
        checkpointer=MemorySaver(),
    )
    config = {"configurable": {"thread_id": "mechanism-test"}}
    first = graph.invoke(
        {
            "proposal": proposal().model_dump(mode="json"),
            "source_commit": commit,
            "task_id": task_id,
            "mechanism_compute_status": "proposed",
        },
        config,
    )
    assert first["mechanism_compute_status"] == "awaiting_compute"
    assert first["__interrupt__"]

    request = store.load_request(task_id, registry)
    result = completed_result(request, registry)
    store.accept_result(result, request, registry)
    resumed = graph.invoke(Command(resume={"result_available": True}), config)
    assert resumed["mechanism_compute_status"] == "ingested"
    assert resumed["result_status"] == "completed"
    assert resumed["fhir_bundle"]["resourceType"] == "Bundle"

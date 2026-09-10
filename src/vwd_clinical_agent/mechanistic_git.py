"""Auditable Git transport for asynchronous mechanism-compute tasks.

The clinical side publishes an immutable request branch.  A server claims that
branch into an isolated worktree, invokes a server-owned protocol runner, and
pushes a result branch.  No executable command is accepted from the request.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tempfile
from typing import Any, Sequence

from .mechanistic_tasks import (
    ArtifactIntegrityError,
    GitProvenanceError,
    MechanismTaskRequest,
    MechanismTaskResult,
    MechanismTaskStore,
    ProtocolRegistry,
    TASK_ID_PATTERN,
    result_template,
    utc_now,
    validate_request_git_provenance,
    validate_task_result,
)


TASK_BRANCH_PREFIX = "mechanism-task"
RESULT_BRANCH_PREFIX = "mechanism-result"


class GitTransportError(RuntimeError):
    """Raised when a Git task cannot move safely to its next state."""


class ProtocolRunnerError(RuntimeError):
    """Raised when a server-owned runner fails to produce a valid result."""


@dataclass(frozen=True)
class GitTaskState:
    task_id: str
    state: str
    branch: str
    commit: str
    worktree: str | None = None
    result_status: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "task_id": self.task_id,
            "state": self.state,
            "branch": self.branch,
            "commit": self.commit,
            "worktree": self.worktree,
            "result_status": self.result_status,
        }


class GitMechanismTaskTransport:
    """One-request-branch/one-result-branch transport with isolated worktrees."""

    def __init__(
        self,
        repo_root: str | Path,
        *,
        remote: str = "origin",
        registry_relative: str = "protocols/vwd_mechanistic_v1/registry.json",
        store_relative: str = "mechanism_tasks",
    ):
        self.repo_root = Path(repo_root).resolve()
        self.remote = remote
        self.registry_relative = self._safe_relative(registry_relative, "registry")
        self.store_relative = self._safe_relative(store_relative, "store")
        top = Path(self._git("rev-parse", "--show-toplevel")).resolve()
        if top != self.repo_root:
            raise GitTransportError(
                f"repo_root must be the Git top level: expected {top}, received {self.repo_root}"
            )

    @staticmethod
    def _safe_relative(value: str, label: str) -> PurePosixPath:
        path = PurePosixPath(value)
        if not value or path.is_absolute() or ".." in path.parts or "." in path.parts:
            raise GitTransportError(f"{label} path must be normalized and repository-relative: {value}")
        return path

    def _git(
        self,
        *args: str,
        cwd: str | Path | None = None,
        check: bool = True,
    ) -> subprocess.CompletedProcess[str] | str:
        completed = subprocess.run(
            ["git", *args],
            cwd=Path(cwd) if cwd is not None else self.repo_root,
            check=False,
            capture_output=True,
            text=True,
        )
        if check and completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip() or f"exit {completed.returncode}"
            raise GitTransportError(f"git {' '.join(args)} failed: {detail}")
        return completed.stdout.strip() if check else completed

    @staticmethod
    def task_branch(task_id: str) -> str:
        if not TASK_ID_PATTERN.fullmatch(task_id):
            raise GitTransportError(f"Invalid mechanism task ID: {task_id!r}")
        return f"{TASK_BRANCH_PREFIX}/{task_id}"

    @staticmethod
    def result_branch(task_id: str) -> str:
        if not TASK_ID_PATTERN.fullmatch(task_id):
            raise GitTransportError(f"Invalid mechanism task ID: {task_id!r}")
        return f"{RESULT_BRANCH_PREFIX}/{task_id}"

    def _local_branch_exists(self, branch: str) -> bool:
        completed = self._git(
            "show-ref",
            "--verify",
            "--quiet",
            f"refs/heads/{branch}",
            check=False,
        )
        assert isinstance(completed, subprocess.CompletedProcess)
        return completed.returncode == 0

    def _remote_branch_exists(self, branch: str) -> bool:
        output = self._git("ls-remote", "--heads", self.remote, f"refs/heads/{branch}")
        assert isinstance(output, str)
        return bool(output)

    def _require_ancestor(self, commit: str, worktree: Path) -> None:
        completed = self._git(
            "merge-base",
            "--is-ancestor",
            commit,
            "HEAD",
            cwd=worktree,
            check=False,
        )
        assert isinstance(completed, subprocess.CompletedProcess)
        if completed.returncode != 0:
            raise GitProvenanceError(
                f"Declared source_commit {commit} is not an ancestor of the task/result branch"
            )

    def _commit(self, worktree: Path, message: str, paths: Sequence[str]) -> str:
        self._git("add", "--", *paths, cwd=worktree)
        self._git(
            "-c",
            "user.name=VWF Mechanism Agent",
            "-c",
            "user.email=mechanism-agent@vwf-ethos.invalid",
            "commit",
            "-m",
            message,
            cwd=worktree,
        )
        commit = self._git("rev-parse", "HEAD", cwd=worktree)
        assert isinstance(commit, str)
        return commit

    def _temporary_worktree(self, ref: str, *, branch: str | None = None) -> tuple[Path, Path]:
        parent = Path(tempfile.mkdtemp(prefix="vwd-mechanism-git-"))
        worktree = parent / "worktree"
        if branch is None:
            self._git("worktree", "add", "--detach", str(worktree), ref)
        else:
            self._git("worktree", "add", "-b", branch, str(worktree), ref)
        return parent, worktree

    def _remove_worktree(self, parent: Path, worktree: Path) -> None:
        try:
            self._git("worktree", "remove", "--force", str(worktree))
        finally:
            shutil.rmtree(parent, ignore_errors=True)

    def publish_request(
        self,
        request: MechanismTaskRequest,
        registry: ProtocolRegistry,
        *,
        push: bool = True,
    ) -> GitTaskState:
        """Commit an immutable request on a branch rooted at source_commit."""

        validate_request_git_provenance(
            request,
            self.repo_root,
            self.repo_root.joinpath(*self.registry_relative.parts),
        )
        if registry.digest != request.acquisition.registry_digest:
            raise GitProvenanceError("Publisher registry does not match the request registry digest")
        branch = self.task_branch(request.task_id)
        if self._local_branch_exists(branch) or (push and self._remote_branch_exists(branch)):
            raise GitTransportError(f"Task branch already exists: {branch}")
        parent, worktree = self._temporary_worktree(request.source_commit, branch=branch)
        try:
            worktree_registry = ProtocolRegistry.load(
                worktree.joinpath(*self.registry_relative.parts)
            )
            store = MechanismTaskStore(worktree.joinpath(*self.store_relative.parts))
            paths = store.submit(request, worktree_registry)
            relative_paths = [path.relative_to(worktree).as_posix() for path in paths]
            commit = self._commit(
                worktree,
                f"task: submit {request.task_id}",
                relative_paths,
            )
            if push:
                self._git("push", self.remote, f"HEAD:refs/heads/{branch}", cwd=worktree)
            return GitTaskState(
                task_id=request.task_id,
                state="awaiting_compute",
                branch=branch,
                commit=commit,
            )
        finally:
            self._remove_worktree(parent, worktree)

    def claim_request(
        self,
        task_id: str,
        worktree_root: str | Path,
        *,
        worker_id: str,
        push: bool = True,
    ) -> GitTaskState:
        """Claim a remote task into a persistent, isolated server worktree."""

        task_branch = self.task_branch(task_id)
        result_branch = self.result_branch(task_id)
        if self._local_branch_exists(result_branch) or (push and self._remote_branch_exists(result_branch)):
            raise GitTransportError(f"Result branch already exists; task is already claimed: {result_branch}")
        remote_ref = f"refs/remotes/{self.remote}/{task_branch}"
        self._git(
            "fetch",
            self.remote,
            f"refs/heads/{task_branch}:{remote_ref}",
        )
        worktree = Path(worktree_root).resolve() / task_id
        if worktree.exists():
            raise GitTransportError(f"Claim worktree already exists: {worktree}")
        worktree.parent.mkdir(parents=True, exist_ok=True)
        self._git("worktree", "add", "-b", result_branch, str(worktree), remote_ref)
        try:
            store = MechanismTaskStore(worktree.joinpath(*self.store_relative.parts))
            request = store.read_request(task_id)
            validate_request_git_provenance(
                request,
                worktree,
                worktree.joinpath(*self.registry_relative.parts),
            )
            self._require_ancestor(request.source_commit, worktree)
            task_commit = self._git("rev-parse", "HEAD", cwd=worktree)
            assert isinstance(task_commit, str)
            claim_path = store.claim_dir(task_id) / "claim.json"
            store._write_once(
                claim_path,
                {
                    "schema_version": "1.0.0",
                    "task_id": task_id,
                    "request_digest": request.request_digest,
                    "task_branch": task_branch,
                    "task_commit": task_commit,
                    "result_branch": result_branch,
                    "worker_id": worker_id,
                    "claimed_at": utc_now(),
                },
            )
            commit = self._commit(
                worktree,
                f"task: claim {task_id}",
                [claim_path.relative_to(worktree).as_posix()],
            )
            if push:
                self._git("push", self.remote, f"HEAD:refs/heads/{result_branch}", cwd=worktree)
            return GitTaskState(
                task_id=task_id,
                state="claimed",
                branch=result_branch,
                commit=commit,
                worktree=str(worktree),
            )
        except Exception:
            self._git("worktree", "remove", "--force", str(worktree), check=False)
            raise

    def _validate_git_artifact_locations(
        self,
        result: MechanismTaskResult,
        task_id: str,
    ) -> None:
        prefix = self.store_relative / "artifacts" / task_id
        for artifact in result.artifacts:
            path = PurePosixPath(artifact.path)
            try:
                path.relative_to(prefix)
            except ValueError as exc:
                raise ArtifactIntegrityError(
                    f"Git-delivered artifact must be below {prefix.as_posix()}: {artifact.path}"
                ) from exc

    def run_claimed(
        self,
        task_id: str,
        worktree: str | Path,
        runner_command: Sequence[str],
        *,
        push: bool = True,
    ) -> GitTaskState:
        """Execute a server-owned runner and publish its validated result."""

        worktree = Path(worktree).resolve()
        expected = Path(worktree).parent / task_id
        if worktree != expected:
            raise GitTransportError(
                f"Worker path must end in the exact task ID: expected {expected}, received {worktree}"
            )
        branch = self._git("branch", "--show-current", cwd=worktree)
        assert isinstance(branch, str)
        if branch != self.result_branch(task_id):
            raise GitTransportError(
                f"Worker must run on {self.result_branch(task_id)}, received {branch or 'detached HEAD'}"
            )
        store = MechanismTaskStore(worktree.joinpath(*self.store_relative.parts))
        request = store.read_request(task_id)
        registry = validate_request_git_provenance(
            request,
            worktree,
            worktree.joinpath(*self.registry_relative.parts),
        )
        self._require_ancestor(request.source_commit, worktree)
        result_path = store.result_dir(task_id) / "result.json"
        template_path = store.result_dir(task_id) / "result.template.json"
        store._write_once(template_path, result_template(request, registry))
        replacements = {
            "task_id": task_id,
            "request": str(store.request_dir(task_id) / "request.json"),
            "result": str(result_path),
            "template": str(template_path),
            "repo_root": str(worktree),
        }
        command = [part.format_map(replacements) for part in runner_command]
        if not command:
            raise ProtocolRunnerError("Server runner command must not be empty")
        environment = os.environ.copy()
        environment.update(
            {
                "VWF_MECHANISM_TASK_ID": task_id,
                "VWF_MECHANISM_REQUEST": replacements["request"],
                "VWF_MECHANISM_RESULT": replacements["result"],
                "VWF_MECHANISM_RESULT_TEMPLATE": replacements["template"],
                "VWF_MECHANISM_REPO_ROOT": replacements["repo_root"],
            }
        )
        completed = subprocess.run(
            command,
            cwd=worktree,
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )
        if not result_path.is_file():
            raise ProtocolRunnerError(
                f"Runner exited {completed.returncode} without writing {result_path}; "
                f"stderr={completed.stderr.strip()[:1000]}"
            )
        result = MechanismTaskResult.model_validate_json(result_path.read_text(encoding="utf-8"))
        if completed.returncode != 0 and result.status == "completed":
            raise ProtocolRunnerError(
                "A non-zero runner exit cannot publish a completed result; return failed or inconclusive"
            )
        return self.publish_claimed_result(task_id, worktree, push=push)

    def publish_claimed_result(
        self,
        task_id: str,
        worktree: str | Path,
        *,
        push: bool = True,
    ) -> GitTaskState:
        """Validate and publish a result written by a human- or AI-operated worker."""

        worktree = Path(worktree).resolve()
        expected = Path(worktree).parent / task_id
        if worktree != expected:
            raise GitTransportError(
                f"Worker path must end in the exact task ID: expected {expected}, received {worktree}"
            )
        branch = self._git("branch", "--show-current", cwd=worktree)
        assert isinstance(branch, str)
        if branch != self.result_branch(task_id):
            raise GitTransportError(
                f"Worker must publish from {self.result_branch(task_id)}, received {branch or 'detached HEAD'}"
            )
        store = MechanismTaskStore(worktree.joinpath(*self.store_relative.parts))
        request = store.read_request(task_id)
        registry = validate_request_git_provenance(
            request,
            worktree,
            worktree.joinpath(*self.registry_relative.parts),
        )
        self._require_ancestor(request.source_commit, worktree)
        result_path = store.result_dir(task_id) / "result.json"
        if not result_path.is_file():
            raise FileNotFoundError(f"Worker result does not exist: {result_path}")
        result = MechanismTaskResult.model_validate_json(result_path.read_text(encoding="utf-8"))
        self._validate_git_artifact_locations(result, task_id)
        validate_task_result(result, request, registry, artifact_root=worktree)
        store.accept_result(result, request, registry, artifact_root=worktree)
        bundle_path = store.ingest_result(result, request, registry, artifact_root=worktree)
        commit_paths = [
            result_path.relative_to(worktree).as_posix(),
            bundle_path.relative_to(worktree).as_posix(),
            *[artifact.path for artifact in result.artifacts],
        ]
        commit = self._commit(
            worktree,
            f"task: return {task_id} ({result.status})",
            commit_paths,
        )
        if push:
            self._git("push", self.remote, f"HEAD:refs/heads/{branch}", cwd=worktree)
        return GitTaskState(
            task_id=task_id,
            state="result_ready",
            branch=branch,
            commit=commit,
            worktree=str(worktree),
            result_status=result.status,
        )

    @staticmethod
    def _copy_once(source: Path, destination: Path) -> None:
        if destination.exists():
            if not destination.is_file() or source.read_bytes() != destination.read_bytes():
                raise FileExistsError(f"Refusing to overwrite immutable collected file: {destination}")
            return
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    def collect_result(self, task_id: str) -> GitTaskState:
        """Fetch, verify, and ingest a server result without merging its branch."""

        branch = self.result_branch(task_id)
        remote_ref = f"refs/remotes/{self.remote}/{branch}"
        self._git("fetch", self.remote, f"refs/heads/{branch}:{remote_ref}")
        parent, worktree = self._temporary_worktree(remote_ref)
        try:
            remote_store = MechanismTaskStore(worktree.joinpath(*self.store_relative.parts))
            request = remote_store.read_request(task_id)
            registry = validate_request_git_provenance(
                request,
                worktree,
                worktree.joinpath(*self.registry_relative.parts),
            )
            self._require_ancestor(request.source_commit, worktree)
            result_path = remote_store.result_dir(task_id) / "result.json"
            result = MechanismTaskResult.model_validate_json(result_path.read_text(encoding="utf-8"))
            self._validate_git_artifact_locations(result, task_id)
            validate_task_result(result, request, registry, artifact_root=worktree)

            local_store = MechanismTaskStore(self.repo_root.joinpath(*self.store_relative.parts))
            local_store.submit(request, registry)
            for artifact in result.artifacts:
                self._copy_once(
                    worktree.joinpath(*PurePosixPath(artifact.path).parts),
                    self.repo_root.joinpath(*PurePosixPath(artifact.path).parts),
                )
            local_result = local_store.accept_result(
                result,
                request,
                registry,
                artifact_root=self.repo_root,
            )
            local_store.ingest_result(
                result,
                request,
                registry,
                artifact_root=self.repo_root,
            )
            commit = self._git("rev-parse", "HEAD", cwd=worktree)
            assert isinstance(commit, str)
            return GitTaskState(
                task_id=task_id,
                state="ingested",
                branch=branch,
                commit=commit,
                result_status=result.status,
                worktree=str(local_result.parent),
            )
        finally:
            self._remove_worktree(parent, worktree)


def runner_command_from_config(
    config_path: str | Path,
    protocol_id: str,
) -> list[str]:
    """Resolve a command only from a server-owned allowlist keyed by protocol."""

    payload = json.loads(Path(config_path).read_text(encoding="utf-8"))
    runners = payload.get("runners")
    if not isinstance(runners, dict):
        raise ProtocolRunnerError("Runner config requires a 'runners' mapping")
    command = runners.get(protocol_id)
    if not isinstance(command, list) or not command or not all(isinstance(item, str) for item in command):
        raise ProtocolRunnerError(f"No server-owned runner configured for {protocol_id}")
    return command

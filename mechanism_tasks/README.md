# Mechanism task exchange

This directory is an append-only Git handoff surface between the networked
clinical agent and an offline compute agent.

```text
requests/<task_id>/request.json       immutable machine contract
requests/<task_id>/task.fhir.json     FHIR R5 Task representation
claims/<task_id>/claim.json           server claim and task-commit provenance
results/<task_id>/result.json         server result after contract validation
ingested/<task_id>/bundle.fhir.json   local Task + Observation + artifacts
artifacts/<task_id>/*                  compact, checksummed Git artifacts only
```

The transport uses one `mechanism-task/<task_id>` request branch and one
`mechanism-result/<task_id>` result branch. The result branch is also the
atomic claim: a second worker cannot push a competing claim without a
non-fast-forward failure. Do not put raw trajectories or large binary
checkpoints in Git.

Client lifecycle:

```bash
python scripts/mechanism_git.py publish examples/mechanism_tasks/a1_v1316m_proposal.json
python scripts/mechanism_git.py collect <task_id>
```

Server lifecycle:

```bash
python scripts/mechanism_git.py claim <task_id> \
  --worktree-root /srv/vwf-mechanism-worktrees --worker-id gpu-01
# server agent executes the locked protocol and writes result.json
python scripts/mechanism_git.py return <task_id> \
  --worktree /srv/vwf-mechanism-worktrees/<task_id>
```

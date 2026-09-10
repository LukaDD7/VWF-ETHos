# Mechanism task exchange

This directory is an append-only Git handoff surface between the networked
clinical agent and an offline compute agent.

```text
requests/<task_id>/request.json       immutable machine contract
requests/<task_id>/task.fhir.json     FHIR R5 Task representation
results/<task_id>/result.json         server result after contract validation
ingested/<task_id>/bundle.fhir.json   local Task + Observation + artifacts
```

Use one `compute/<task_id>` branch per task. Do not put raw trajectories or
large binary checkpoints in Git.

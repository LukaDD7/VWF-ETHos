# Type 2M LOF MD partial analysis

Generated on 2026-06-29 from completed or checkpointed production MD runs only. Large trajectory files are intentionally ignored by git.

## Current progress

- A1-GPIb queue: 9/39 complete, 7 running/checkpointed, 23 not started by the scheduler yet.
- 7A6O closed-state queue: 10/32 complete, 4 running/checkpointed, 18 not started by the scheduler yet.
- 7A6O completed-run QC was generated for the first 10 complete closed-state runs.

## Remaining time estimate

These are queue-level estimates from the tmux scheduler timestamps, not guarantees. They assume current GPU sharing and no repeated failures.

- A1-GPIb: roughly 36-54 hours remaining. The queue uses GPUs 0-6, with 7 production runs in flight. Two variants (`L1276R`, `R1315G`) had an early `rc=134` retry event and may add delay if they fail again.
- 7A6O closed-state: roughly 40-55 hours remaining. The queue uses GPUs 0-3, with 4 production runs in flight.
- Both queues should likely finish in about 2 days if the current throughput holds.

## Files

- `queue_progress_summary.csv`: per-queue completion counts.
- `a1_gpiba_completed_and_running_summary.csv`: A1-GPIb queue status and completed-run metadata.
- `7a6o_closed_state_completed_and_running_summary.csv`: 7A6O closed-state queue status and completed-run metadata.
- `combined_completed_and_running_summary.csv`: combined view.
- `7a6o_completed_qc/qc_summary.csv`: frame counts and backbone RMSD for completed 7A6O closed-state runs.
- `7a6o_completed_qc/aim_a1_contacts_summary.csv`: AIM-A1 contact summary for completed 7A6O closed-state runs.
- `7a6o_completed_qc/aim_a1_contacts_timeseries.csv`: AIM-A1 contact time series.
- `scheduler_logs/`: tmux scheduler snapshots used for progress and estimates.

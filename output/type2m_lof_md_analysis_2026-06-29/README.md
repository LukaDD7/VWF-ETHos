# Type 2M LOF MD partial analysis

Generated on 2026-06-30 from completed or checkpointed production MD runs only. Large trajectory files are intentionally ignored by git.

## Current progress

- A1-GPIb queue: 25/39 complete, 7 running/checkpointed, 7 pending or waiting for retry.
- 7A6O closed-state queue: 22/32 complete, 3 running/checkpointed, 7 pending or waiting for retry.
- A1-GPIb completed-run QC was generated for the first 25 complete runs.
- 7A6O completed-run QC and AIM-A1 contact features were generated for the first 22 complete closed-state runs.

## Remaining time estimate

These are queue-level estimates from the tmux scheduler timestamps, not guarantees. They assume current GPU sharing and no repeated failures.

- A1-GPIb: roughly 12-24 hours remaining. Retry candidates include `L1276R`, `R1315G`, and `E1290L`.
- 7A6O closed-state: roughly 12-24 hours remaining. Retry candidate observed: `R1374S`.

## Files

- `queue_progress_summary.csv`: per-queue completion counts.
- `a1_gpiba_completed_and_running_summary.csv`: A1-GPIb queue status and completed-run metadata.
- `7a6o_closed_state_completed_and_running_summary.csv`: 7A6O closed-state queue status and completed-run metadata.
- `combined_completed_and_running_summary.csv`: combined view.
- `a1_gpiba_completed_qc/qc_summary.csv`: frame counts and backbone RMSD for completed A1-GPIb runs.
- `7a6o_completed_qc/qc_summary.csv`: frame counts and backbone RMSD for completed 7A6O closed-state runs.
- `7a6o_completed_qc/aim_a1_contacts_summary.csv`: AIM-A1 contact summary for completed 7A6O closed-state runs.
- `7a6o_completed_qc/aim_a1_contacts_timeseries.csv`: AIM-A1 contact time series.
- `scheduler_logs/`: tmux scheduler snapshots used for progress and estimates.

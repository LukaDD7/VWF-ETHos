# Type 2M LOF MD analysis

Updated: 2026-07-08

Generated from completed production MD runs only. Large trajectory files are not committed here.

Current usable completed set:
- A1-GPIb: 37 / 39 complete.
- 7A6O closed state: 31 / 32 complete.
- Overall: 68 / 71 usable completed MD runs.

The remaining three jobs were stopped as failed or withheld after repeated retry/hang states:
- A1-GPIb: `L1276R`, `R1315G`.
- 7A6O closed state: `R1374S`.

See `FAILED_VARIANTS.md` and `failed_variants.csv` for details. These three variants should be excluded from the current classifier and feature-signal validation.

Files:
- `queue_progress_summary.csv`: per-queue completion counts and median completed wall time.
- `a1_gpiba_completed_and_running_summary.csv`: A1-GPIb queue status and completed-run metadata.
- `7a6o_closed_state_completed_and_running_summary.csv`: 7A6O closed-state queue status and completed-run metadata.
- `combined_completed_and_running_summary.csv`: combined view.
- `a1_gpiba_completed_qc/`: GROMACS RMSD QC for 37 completed A1-GPIb runs.
- `7a6o_completed_qc/`: GROMACS RMSD and AIM-A1 contact QC for 31 completed 7A6O closed-state runs.

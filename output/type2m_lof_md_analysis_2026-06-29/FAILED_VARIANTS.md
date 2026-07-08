# Failed or withheld MD variants

Updated: 2026-07-08

The long-running MD queue was not simply slow. Three remaining jobs were in repeated retry or hang states and were stopped to prevent an unbounded loop. The usable completed MD set is 68 / 71 runs.

Summary:
- A1-GPIb: 37 / 39 complete, 2 failed or withheld.
- 7A6O closed state: 31 / 32 complete, 1 failed or withheld.

Failed or withheld variants:
- `a1_gpiba` `L1276R`: stuck in NVT/retry state without a completed `md_prod.gro`.
- `a1_gpiba` `R1315G`: repeated GROMACS CUDA illegal memory access failure during NVT, without a completed `md_prod.gro`.
- `7a6o_closed_state` `R1374S`: repeated NVT failures/hangs, without a completed `md_prod.gro`.

Use recommendation:
- Exclude these three variants from the current classifier and feature-signal validation.
- Treat the 68 completed runs as the current usable MD-derived feature set.
- The resilient MD schedulers now have `--max-attempts` with default `5`, so future failures stop with a `.failed` marker instead of retrying indefinitely.

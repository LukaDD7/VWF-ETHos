# Type-1 10-case server bundle

- AlphaGenome: 9 normalized variants x 11 official raw OutputTypes; metadata-selected VWF-relevant tracks; all official recommended scorer views.
- Boltz-2: 5 missense variants, 8 variant-assay jobs, 7 matched WT baselines = 15 jobs total.
- MD: one targeted 7A6O matched-WT/P1413L pair, 50 ns each; not a blanket MD run for all missense variants.
- Excluded from direct AlphaGenome/Boltz input: CASE_T1_004 large deletion/duplication only.
- CASE_T1_007: VWF c.3379+1G>A is included; hemophilia A is retained as a comorbidity and FVIII:C is blocked from type-2N inference.

Run in order:

1. `bash run_alphagenome.sh` (API; needs `ALPHAGENOME_API_KEY`)
2. `GPUS=4 bash run_boltz.sh` (GPU)
3. `FOLDX=/path/to/foldx bash run_md_p1413l.sh` (2 GPUs recommended)
4. Return the `results/` directory and run `bash ingest_and_test_agent.sh`

Expected inventory: AlphaGenome 9/9 successful cases, Boltz 15/15 completed jobs, and MD QC for `TYPE1_WT_MATCHED` plus `P1413L`.

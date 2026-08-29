# Type-2B six-patient server bundle

- Deidentified input: six patients, seven variant records; patient 5 has A1461D plus benign-control D2449N.
- AlphaGenome: 7 requests x 11 official raw outputs, metadata-selected VWF-relevant tracks and all recommended scorer views.
- Boltz-2: 6 unique missense variants; 16 variant-assay jobs + 4 assay-matched WT baselines = 20 jobs.
- MD: five unique A1 variants plus one pipeline-matched WT, 50 ns each. D2449N (C3 benign control) is not sent to A1 MD.
- V1316M is computed once and mapped back to both patients.

Run `run_alphagenome.sh`, `run_boltz.sh`, and `run_md_a1_panel.sh`; then run `ingest_and_test_agent.sh` after all result inventories pass.

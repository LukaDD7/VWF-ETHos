# VWD/VWF Functional Boltz-2 Panel

This directory contains pre-generated Boltz-2 YAML inputs for the VWD/VWF
functional evidence panel.

## Contents

- `yamls/`: one Boltz-2 YAML per job.
- `job_manifest.csv`: index of all runnable jobs and their biological axis.
- `variants_master.csv`: de-duplicated source variant table.
- `diagnostic_panel.csv`: variant-by-assay matrix, including skipped axes.
- `json_batches/`: optional batch JSONs for older local runner conventions.
- `summary.json`: generated counts.

## Run On GPU Server

From the repository root:

```bash
bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh
```

Common overrides:

```bash
DEVICES=4 \
RECYCLING_STEPS=3 \
DIFFUSION_SAMPLES=5 \
NUM_WORKERS=8 \
bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh
```

If the server can access the MSA server:

```bash
EXTRA_BOLTZ_ARGS="--use_msa_server" \
bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh
```

## Panel Size

- 597 unique missense variants.
- 15 functional assay axes.
- 990 YAML jobs total: 975 variant jobs plus 15 WT baselines.

WT baselines are for per-assay delta metrics. Negative/control variants are
included separately from GeneBe and ClinVar sources.

# New-sample mechanism-model run plan

Date: 2026-08-29

## Required input contract

For each new sample, provide at minimum:

```text
patient_id,source_row_id,hgvs_c,hgvs_p,aa_change,
chromosome,position,ref,alt
```

`aa_change` must use the compact one-letter form, for example `R1306W`. The
current `AgenticVWFClassifier` is designed for missense variants; nonsense,
frameshift, splice, synonymous, and structural variants need a separate rule
path rather than being forced through the missense feature matrix.

## AlphaGenome

1. Primary biosample: `CL:0002618` (endothelial cell of umbilical vein).
2. Request all 11 output types where available.
3. For ATAC and contact maps, issue a second all-track request with
   `ontology_terms=None`; label those features `source_scope=all_tracks`.
4. Preserve signed `alternate - reference` values, peak positions, track names,
   and ontology terms. Missing modalities must remain missing with a reason.

## Boltz-2

For missense variants:

```bash
python scripts/pipeline/generate_vwd_functional_boltz2_yamls.py \
  --variants-csv /path/to/new_samples.csv \
  --output-dir output/boltz2_new_samples

bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh \
  --input-dir output/boltz2_new_samples/yamls \
  --out-dir output/boltz2_new_samples/boltz_results \
  --gpus 8
```

The generator supports a clean `aa_change/wt_aa/position/mut_aa` CSV and writes
mechanism-aware assay YAMLs plus a job manifest.

## Agent integration

Merge AlphaGenome, Boltz, AF3/FoldX, and any MD features into one matrix keyed
by `aa_change` or `hgvs_c`, then run:

```bash
python scripts/run_vwd_langgraph_v0.py \
  --mode retrospective \
  --provider-profile offline \
  --llm-provider deterministic \
  --biomedical-tools \
  --mechanism-matrix output/new_samples_mechanism_features.csv \
  --mechanism-classifier scripts/agentic_vwf_classifier.py \
  --output-dir output/vwd_langgraph_v0/new_samples
```

The mechanism adapter is read-only, ignores label columns at inference time,
and does not bypass expert review.

## MD planning

Preferred starting structure for the AIM-A1/2B axis is experimental `7A6O`, not
de novo Boltz D'D3-A1. Mutants should be built on the same WT backbone and
passed through staged EM/NVT/NPT before production.

On the shared H200/NFS setup, the existing GROMACS environment can be reused.
On an independent A40 server, rebuild a CUDA GROMACS environment, install
`gemmi/numpy/pandas`, set `GMXLIB` to the patched force-field directory, and
pass `gpu_smoke_test.sh` before any batch run.

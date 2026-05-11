# VWD Functional Panel: Next Analysis Plan

> Purpose: turn Boltz-2 / FoldX / AF3 / AlphaGenome outputs into an evidence
> matrix for an agent that simulates the VWD diagnostic reasoning path.

## 1. Research Objective

The project is not trying to classify VWD from one structure score. The goal is
to build a diagnostic agent that receives patient-level molecular evidence and
asks the same functional questions a clinician would ask:

1. Is this a quantitative VWF problem? Type 1 / Type 3 direction.
2. Is this a qualitative VWF problem? Type 2 direction.
3. If Type 2, which functional axis is impaired or overactive?
4. Is the model evidence consistent with the reported phenotype and lab assay
   pattern?

The structural panel should therefore be interpreted as an evidence generator,
not as the final classifier.

## 2. Data Now Available

Main generated files:

- `output/boltz2_vwd_functional_panel/job_manifest.csv`
- `output/boltz2_vwd_functional_panel/diagnostic_panel.csv`
- `output/boltz2_vwd_functional_panel/boltz_results_summary.csv`
- `output/boltz2_vwd_functional_panel/variants_master.csv`

Approximate panel scale:

- 597 unique missense variants.
- 989 finished Boltz-2 jobs in the current parsed summary.
- 8 diffusion samples per job.
- Positive VWD labels include Type 1, Type 2A, Type 2B, Type 2M, Type 2N, Type 3.
- Negative/control variants include GeneBe and ClinVar benign/non-synonymous
  controls.

## 3. Important Metric Rule

Do not use iPTM for every assay.

Use iPTM for complex/interface assays:

- `dprime_d3_fviii_binding`
- `a1_gpiba_forced_binding`
- `a1_heparan_sulfate_binding`
- `a2_adamts13_folded_complex`
- `vwf73_adamts13_substrate`
- `a3_collagen_binding`
- `c1_collagen_binding`
- `c2_collagen_binding`
- `c4_integrin_binding`

Use pTM / pLDDT / local structural features for monomer or conformational
context assays:

- `a1_aim_autoinhibition_context`
- `a2_folded_stability`
- `d1d2_propeptide_context`
- `d4_assembly_context`
- `c_domain_assembly_context`
- `ck_dimerization_context`

For monomer assays, iPTM equal to zero is expected because there is no
inter-chain interface.

## 4. First Analysis Layer: Feature Table

Build a single row per variant with one feature block per assay.

Recommended feature families:

- Completion features:
  - whether the assay was run
  - number of successful samples
  - missing/error state

- Confidence features:
  - average and best iPTM for complex assays
  - average and best pTM for monomer assays
  - average and best complex pLDDT
  - sample dispersion across 8 diffusion samples

- Delta-to-WT features:
  - mutant metric minus assay-specific WT metric
  - mutant rank among all variants for that assay
  - z-score within assay and within domain

- Label policy features:
  - `positive_any`
  - `negative_only`
  - `negative_positive_overlap`
  - non-exclusive source labels

The immediate deliverable should be:

- `output/boltz2_vwd_functional_panel/evidence_matrix.csv`

## 5. Second Analysis Layer: Functional Interpretation

Map structural signals to clinical axes.

| Axis | Primary evidence | Biological interpretation |
|---|---|---|
| Type 2N | D′D3-FVIII delta iPTM | FVIII carrier binding loss |
| Type 2B | A1-AIM exposure proxy + A1-GPIbα context | A1 opens too easily / gain of platelet binding |
| Type 2M-A1 | A1-GPIbα loss without 2B-like exposure | platelet binding loss |
| Type 2M-A3 | A3-collagen loss | collagen binding loss |
| Type 2A | A2 stability + ADAMTS13 substrate/docking | HMW multimer loss / cleavage susceptibility |
| Type 1/3 | D1D2, D4, C-domain, CK stability plus AlphaGenome | secretion, dimerization, multimerization, expression |

For 2B vs 2M, do not rely on the forced A1-GPIbα complex alone. It answers
"can this local complex be modeled as bound?", not "does A1 open at low shear?".
The AIM-flanked A1 construct is exploratory but closer to the mechanistic
question.

## 6. Third Analysis Layer: Model Calibration

Use labels as weak supervision, not ground truth for every mechanism.

Recommended evaluation:

1. Negative-only controls vs VWD positives:
   - test whether the feature panel separates benign variants from curated VWD
     variants.

2. One-vs-rest subtype calibration:
   - Type 2N vs others using D′D3-FVIII features.
   - Type 2M-A3 vs others using A3/collagen features.
   - Type 2A vs others using A2 stability and ADAMTS13 features.
   - Type 2B vs Type 2M-A1 using A1-GPIbα plus AIM context features.

3. Within-domain controls:
   - compare only variants in the same domain to reduce domain-location leakage.

4. Report uncertainty:
   - high confidence when multiple independent axes agree.
   - uncertain when signal exists only in exploratory monomer/conformational
     features.

## 7. Fourth Analysis Layer: Agent Inputs

The diagnostic agent should receive structured evidence like this:

```json
{
  "variant_id": "VWF_R1306W",
  "clinical_labs": {
    "VWF_Ag": null,
    "VWF_activity": null,
    "FVIII_C": null,
    "multimers": null
  },
  "genomic_evidence": {
    "alphagenome_delta": null,
    "splice_delta": null
  },
  "structural_evidence": {
    "A1_GPIb": {"delta_iptm": null, "direction": "unknown"},
    "A1_AIM": {"delta_ptm": null, "exposure_proxy": null},
    "Dprime_D3_FVIII": {"delta_iptm": null},
    "A2_ADAMTS13": {"delta_iptm": null},
    "A3_collagen": {"delta_iptm": null}
  }
}
```

The final output should be a reasoning trace in clinical order, not just a
single class label.

## 8. Immediate Next Coding Tasks

1. Build `build_vwd_evidence_matrix.py`. Done in
   `scripts/pipeline/build_vwd_evidence_matrix.py`.
   - Join `boltz_results_summary.csv`, `job_manifest.csv`, `variants_master.csv`,
     and WT baselines.
   - Compute per-assay delta, z-score, rank, and primary metric.

2. Build `analyze_vwd_functional_panel.py`. Done in
   `scripts/pipeline/analyze_vwd_functional_panel.py`.
   - Generate assay-level distributions.
   - Generate negative vs positive comparisons.
   - Generate subtype one-vs-rest summaries.

3. Build figures for the presentation.
   - Heatmap: variants x functional axes.
   - Boxplots: subtype vs assay metric delta.
   - Focus panels: 2B vs 2M-A1, 2N, 2A, 2M-A3.

4. Update Agent v3.
   - Ingest evidence matrix instead of hard-coded subtype rules.
   - Keep rule-like clinical ordering, but score evidence probabilistically.

## 9. What Would Count As Success

A successful next phase is not necessarily high multiclass accuracy. It is:

- negative controls mostly look structurally quiet;
- 2N variants show D′D3/FVIII signal enrichment;
- 2M-A3 variants show collagen-axis signal enrichment;
- 2A variants show A2/ADAMTS13 or stability-axis signal enrichment;
- 2B and 2M-A1 are explicitly marked as difficult unless AIM/exposure features
  add signal beyond forced A1-GPIbα binding;
- the agent can explain why a subtype is favored or why the evidence is
  insufficient.

## 10. Server Run Order

After pulling the latest repository on the GPU/result server, run:

```bash
python scripts/pipeline/parse_vwd_functional_boltz2_results.py \
  --results-dir output/boltz2_vwd_functional_panel/boltz_results \
  --manifest output/boltz2_vwd_functional_panel/job_manifest.csv \
  --output output/boltz2_vwd_functional_panel/boltz_results_summary.csv

python scripts/pipeline/build_vwd_evidence_matrix.py

python scripts/pipeline/analyze_vwd_functional_panel.py
```

Expected new outputs:

- `output/boltz2_vwd_functional_panel/evidence_long.csv`
- `output/boltz2_vwd_functional_panel/evidence_matrix.csv`
- `output/boltz2_vwd_functional_panel/analysis/*.csv`

Then commit and push those generated analysis files for review.

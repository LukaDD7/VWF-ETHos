# Boltz-2 VWD Functional Panel Initial Readout

Date: 2026-05-11

This note summarizes the first complete server-side Boltz-2 functional panel
run after parsing into the VWD evidence matrix.

## 1. Data Integrity

The pulled server results contain:

- 597 unique non-WT variants.
- 990 Boltz-2 jobs in `boltz_results_summary.csv`.
- 15 functional assay axes.
- 975 non-WT variant-assay runs plus 15 WT assay baselines.
- 0 planned runs missing a primary metric.

The large number of `no_result` cells in `evidence_long.csv` is expected: the
panel is domain-aware, so each variant is only sent to biologically relevant
assays instead of all 15 assays.

## 2. Metric Rule Used

The evidence matrix uses assay-specific primary metrics:

- Complex/interface assays: iPTM.
- Monomer/context assays: pTM first, then complex pLDDT fallback.

This distinction matters because monomer iPTM is not a binding or structure
quality signal.

## 3. First Biological Signals

Strongest first-pass separation appears in:

- A2 folded stability:
  - negative-only mean delta: +0.066
  - positive-only mean delta: -0.035
  - Cohen d: -1.27
  - Interpretation: this is currently the clearest structural signal in the
    panel, but it should be interpreted as stability/confidence evidence, not
    direct ADAMTS13 cleavage kinetics.

- C-domain assembly context:
  - negative-only mean delta: -0.002
  - positive-only mean delta: -0.018
  - Cohen d: -0.80
  - Interpretation: potentially useful for quantitative VWD axes, especially
    secretion/multimerization/dimerization mechanisms.

- CK dimerization context:
  - negative-only mean delta: +0.008
  - positive-only mean delta: -0.011
  - Cohen d: -0.72
  - Interpretation: plausible quantitative VWD signal, but sample size is small
    because only CK-domain variants are eligible.

## 4. Type 2B vs Type 2M-A1

The direct A1 comparison does not strongly separate Type 2B-only from
Type 2M-only variants:

| Assay | Type 2B-only n | Type 2M-only n | Mean delta difference | Cohen d |
|---|---:|---:|---:|---:|
| A1 AIM autoinhibition context | 38 | 40 | +0.006 | +0.13 |
| A1 GPIbα forced binding | 37 | 40 | +0.005 | +0.08 |
| A1 heparan-sulfate binding | 37 | 40 | +0.050 | +0.48 |

Interpretation:

- The forced A1-GPIbα complex is not enough to distinguish 2B from 2M-A1.
- The AIM-context task also remains weak in the current feature form.
- The heparan-sulfate/charge-modulation axis shows a moderate difference, but
  it is not a canonical clinical discriminator by itself.

This supports the earlier hypothesis: Type 2B vs Type 2M-A1 is a conformational
activation problem, not just a static forced-binding problem. The next feature
layer should look for exposure/contact geometry, AIM packing, A1 opening, and
sample-level conformational heterogeneity rather than only average confidence
metrics.

## 5. Type 2N / FVIII Axis

The D′D3-FVIII assay shows only modest subtype signal:

- Type 2N one-vs-rest Cohen d: +0.41.
- Several Type 2N variants are high outliers, for example G785E and C788Y.

However, the sign is not automatically interpretable as stronger binding. iPTM
is a model confidence/interface confidence score, not binding affinity. For
Type 2N, the next pass should add interface-level features:

- interface PAE or contact confidence;
- FVIII-contact residue distance changes;
- buried surface / contact count proxies;
- comparison with known Type 2N residues in D′D3.

## 6. ADAMTS13 Axis

ADAMTS13 should not be treated as a simple ligand-affinity problem. The relevant
clinical biology is closer to:

- A2 folded stability;
- exposure of the VWF73 cleavage segment;
- accessibility of the Tyr1605-Met1606 cleavage site;
- force/shear-dependent unfolding susceptibility.

The current panel gives useful first evidence from A2 folded stability, but the
ADAMTS13 complex iPTM should be considered supportive, not definitive.

## 7. Implication For The Diagnostic Agent

The next agent should not take a single Boltz metric and output a subtype.
Instead, each assay should become an evidence item with:

- axis name;
- primary metric;
- delta vs WT;
- percentile / z-score;
- whether the axis is mechanistically direct or exploratory;
- caveat text for the clinician-style reasoning trace.

Recommended high-level agent logic:

1. Use VWF:Ag, VWF activity, FVIII:C, multimers, and sequencing evidence as the
   clinical scaffold when available.
2. Use AlphaGenome/splicing for quantitative or transcript-disruptive evidence.
3. Use Boltz/FoldX/AF3 as mechanism-specific evidence generators.
4. Report uncertainty explicitly, especially for 2B vs 2M-A1.

## 8. Next Coding Tasks

1. Add structure-derived local feature extractors:
   - interface contact counts;
   - interface PAE summaries if PAE is available;
   - residue-distance changes around known functional residues;
   - per-sample dispersion across 8 Boltz samples.

2. Build presentation figures from `evidence_matrix.csv`:
   - functional-axis heatmap;
   - subtype-vs-assay boxplots;
   - 2B vs 2M-A1 focused panel;
   - Type 2N D′D3/FVIII panel;
   - A2/ADAMTS13 panel.

3. Add an Agent v3 evidence schema:
   - one row per patient/variant;
   - multiple independent evidence channels;
   - clinician-style conclusion with caveats.

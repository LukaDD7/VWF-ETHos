# Hotspot prior leakage risk and removal decision

Date: 2026-06-29

## Decision

Remove the A1 Type 2B hotspot/position prior from the production classifier.

The previous rule treated recurrent Type 2B-associated A1 residues as a direct
2B prior. Even though the biological motivation is literature-based clustering
of Type 2B mutations in the VWF A1 / exon 28 / GPIb-binding region, this axis is
too close to memorizing subtype labels by residue position. It is likely to be
challenged in review as circular prior knowledge unless the hotspot set is
externally frozen, prospectively validated, and never derived from the Eval set.

For the current research prototype, accepting a lower metric is preferable to
retaining a position-memory axis.

## Latest Eval impact

Latest Eval artifact:

`output/type2m_lof_md_fast_validation_2026-06-29/eval_with_fast_md/eval_v2_predictions.csv`

Hotspot-rule calls in the with-MD classifier:

| Metric | Value |
|---|---:|
| Calls explicitly made by hotspot rule | 29 |
| Correct hotspot-rule calls | 22 |
| Incorrect hotspot-rule calls | 7 |
| Correct Type 2 predictions overall | 134 |
| Correct Type 2 predictions due to hotspot rule | 22 |
| Share of correct Type 2 predictions due to hotspot rule | 16.4% |

True labels among the 29 hotspot-rule calls:

| True label | n | Predicted label |
|---|---:|---|
| 2B | 22 | 2B |
| 2M | 6 | 2B |
| 1 | 1 | 2B |

Temporary hotspot-ablation result:

| Label | Current recall | No-hotspot recall | Delta correct |
|---|---:|---:|---:|
| 2A | 70/118 = 59.3% | 70/118 = 59.3% | 0 |
| 2B | 29/38 = 76.3% | 15/38 = 39.5% | -14 |
| 2M | 23/53 = 43.4% | 25/53 = 47.2% | +2 |
| 2N | 12/16 = 75.0% | 12/16 = 75.0% | 0 |
| Type 2 total | 134/225 = 59.6% | 122/225 = 54.2% | -12 |

The hotspot prior mainly boosts 2B recall but also masks A1 Type 2M LOF cases
that should be decided by independent functional and dynamics evidence.

## Replacement principle

After removal, A1 2B/2M calls should be driven by independent mechanism axes:

- A1-GPIb / heparan static complex evidence from Boltz-like predictors.
- A1 binding-face destabilization from equilibrium MD, used as soft Type 2M LOF
  evidence.
- Salt-bridge/contact retention from equilibrium MD, used only as anti-2M or
  face-retained evidence, not as a standalone 2B label.
- Future SMD / force-MD only for high-value mechanistic validation of force
  release, not as a broad position prior.

Position can still be reported for human review, but it should not determine a
subtype by itself.

# Targeted 7A6O AIM-A1 MD results

This directory contains the completed equilibrium-MD analysis for seven matched
7A6O AIM-flanked A1 systems:

- `PANEL_WT_MATCHED`
- `P1413L`
- `R1308C`
- `S1310F`
- `V1316M`
- `R1341W`
- `A1461D`

Each production trajectory has 501 frames at 100 ps spacing (50 ns). The
matched WT baseline is `PANEL_WT_MATCHED`; it should not be replaced by an
unmatched generic WT run when comparing this panel.

## Result layers

- `qc_summary.csv`: trajectory QC and backbone-RMSD summary.
- `aim_a1_contacts_summary.csv`: contact-count summary for AIM_all, N-AIM, and
  C-AIM versus A1 non-local residues.
- `aim_a1_contacts_timeseries.csv`: integrated per-frame contact counts.
- `rmsd_timeseries.csv`: integrated per-frame backbone RMSD.
- `aim_a1_min_distance_summary.csv`: minimum-distance summary for the three AIM
  selections.
- `aim_a1_min_distance_timeseries.csv`: integrated per-frame minimum distances.
- `md_feature_matrix.csv`: one row per system combining QC, RMSD, AIM contact,
  and minimum-distance features with matched-WT deltas.
- `artifact_manifest.csv`: size and SHA-256 inventory for derived summaries,
  XVG/NDX intermediates, and local trajectory files.

The `contacts/`, `rmsd/`, and `index/` directories contain the lightweight
GROMACS analysis intermediates used to build the CSV layers. The `trajectories/`
directory contains seven approximately 62 MB concatenated XTC files and is
intentionally retained as local source data rather than committed to Git. Its
files and checksums are listed in `artifact_manifest.csv` with
`git_policy=local_only`.

## Reproduction

The original analysis is produced by:

```bash
python scripts/pipeline/analyze_7a6o_completed_md.py \
  --input output/gromacs_md_autoinhib \
  --output outputs/computational_panel_20260829/md \
  --gmx /path/to/gmx \
  --variants PANEL_WT_MATCHED,P1413L,R1308C,S1310F,V1316M,R1341W,A1461D \
  --wt-variant PANEL_WT_MATCHED
```

The additional lightweight summaries are produced from those XVG files by:

```bash
python scripts/pipeline/summarize_computational_panel_md.py \
  --md-dir outputs/computational_panel_20260829/md \
  --wt-variant PANEL_WT_MATCHED
```

## Interpretation limits

AIM-A1 contact loss or gain is an exploratory conformational-exposure proxy.
It does not directly measure shear response, platelet binding, or VWD subtype.
Interpret it together with assay-matched Boltz evidence, clinical labs, and the
2B/2M discriminators. The contact threshold used by the GROMACS analysis is
0.45 nm; minimum-distance features provide a complementary threshold-free view.

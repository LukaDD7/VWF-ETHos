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
- `module_structure_inventory.csv`: preferred starting-structure inventory for
  all 15 functional assay axes. Experimental PDBs are listed first; Boltz-2 is
  only a fallback for modules without a downloaded experimental structure.
- `module_md_readout_plan.csv`: module-specific equilibrium-MD readouts,
  conditional SMD readouts, SMD gates, AI reference-distribution features, and
  execution priority.
- `md_result_layers.csv`: consolidated index of the additional completed
  A1/AIM, A1-GPIbα, and slow025 SMD result layers.

## Additional completed MD/SMD layers

The bundle now also includes lightweight result layers that were previously
kept only under `output/`:

- `a1_aim_masking/`: AIM-A1 masking-interface features and residue metadata.
- `a1_aim_saltbridge/`: AIM-A1 salt-bridge occupancy and retained-z features.
- `a1_gpiba_interface/`: A1-GPIbα equilibrium-MD interface-retention features,
  sanitized summary, and per-frame time series.
- `smd/`: the three-control slow025 SMD calibration and per-replicate force/work
  curves.

These layers are intentionally lightweight. Raw trajectories, TPR/GRO files,
and machine-specific absolute paths are not copied into the bundle. The SMD
layer is marked `complete_no_go`: it documents why the current force axis must
not be used as a classifier feature.

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

The multi-module inventory and consolidated MD/SMD layers can be regenerated
with:

```bash
python scripts/pipeline/build_vwf_md_structure_inventory.py
python scripts/pipeline/consolidate_vwf_md_result_layers.py
```

## Interpretation limits

AIM-A1 contact loss or gain is an exploratory conformational-exposure proxy.
It does not directly measure shear response, platelet binding, or VWD subtype.
Interpret it together with assay-matched Boltz evidence, clinical labs, and the
2B/2M discriminators. The contact threshold used by the GROMACS analysis is
0.45 nm; minimum-distance features provide a complementary threshold-free view.

Do not mix experimental crystal structures with Boltz-predicted models. Where
an experimental PDB exists, it supersedes Boltz. The inventory currently marks
7A6O, 1SQ0, 3GXB, 4DMU, 6N29, and 4NT5 as locally available. Boltz-2 is used
only for modules without an available experimental structure, and
low-confidence Boltz models remain exploratory.

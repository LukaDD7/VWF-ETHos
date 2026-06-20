# 7A6O AIM-A1 Reference MD Dataset

Portable MD data package for WT + 14 completed reference variants used to derive MD-based 2B/2M classification features.

Contents:
- `variants/<variant>/prod_concat.xtc`: concatenated 50 ns production trajectory.
- `variants/<variant>/md_prod.tpr`: GROMACS run input for selections/analysis.
- `variants/<variant>/final.gro`: final production coordinates.
- `analysis_completed_7a6o/*.csv`: existing QC and AIM-A1 contact summaries/timeseries.
- `manifest.csv`: labels, relative paths, byte sizes, and SHA256 checksums.

Labels are reference labels only: WT, known 2B, known 2M, and `?` for uncertain variants. Patient C1458R replicates are intentionally not included in this reference package.

Example analysis input after clone with Git LFS enabled:

```bash
git lfs pull
python3 scripts/pipeline/analyze_7a6o_completed_md.py --help
```

For custom feature extraction, use each variant directory's `prod_concat.xtc` with its matching `md_prod.tpr`.

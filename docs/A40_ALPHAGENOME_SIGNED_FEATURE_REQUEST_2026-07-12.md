# A40 request: recover signed, VWF-localized AlphaGenome features tonight

## Objective

Recover compact, scientifically interpretable AlphaGenome features from the full raw prediction objects on the A40 server. Do **not** transfer the 197 GB pickle. Transfer only the derived table, manifest, and a small QC sample.

The current local CSV stores `max(abs(ALT - REF))`. That removes effect direction and genomic localization, so a small value cannot be interpreted as RNA loss. The new output must preserve `ALT - REF` and the position/track that generated each summary.

## Locate the raw data

From `/lzy/projects/VWF-ETHos` (or the current project checkout):

```bash
find /lzy -type f \
  \( -name '03_inference_results.pkl' \
     -o -name '07_VCF_AlphaGenome_Results*.pkl' \
     -o -name '*AlphaGenome*.pkl' \) \
  -size +1G -print
```

Expected candidates from the project history:

- `03_inference_results.pkl` — 1,198 variants, RNA-seq and splice-sites raw outputs.
- `07_VCF_AlphaGenome_Results_11_Modalities_Final.pkl` — 2,577 variants, up to 11 modalities.

Record the exact input path, file size, modification time, and SHA256 in the manifest.

## Implementation requirements

Create:

```text
scripts/pipeline/extract_alphagenome_signed_features.py
```

The extractor must:

1. Read one variant at a time; never load or copy the full 197 GB object unnecessarily.
2. Calculate all deltas as `alternate - reference` before any absolute-value operation.
3. Keep clinical subtype labels completely out of feature construction.
4. Restrict gene-aware summaries to VWF and report the genomic interval/bin used.
5. Preserve ontology/cell-type and track metadata where available.
6. Emit `NaN` plus a reason code for unavailable modalities; do not silently write zero.
7. Keep the original variant key: chromosome, position, REF, ALT, and consequence/cDNA mapping.

## Required output columns

### Identity and QC

```text
variant_key, chromosome, position, ref, alt, cdna_change, aa_change
interval_start, interval_end, ontology_term
raw_rna_shape, raw_splice_shape, n_rna_tracks, n_splice_tracks
qc_status, qc_message
```

### RNA-seq

For the full returned window and, if metadata permit, VWF promoter/gene-body subsets:

```text
ag_rna_signed_min                  # most negative ALT-REF
ag_rna_signed_max                  # most positive ALT-REF
ag_rna_abs_max
ag_rna_signed_mean
ag_rna_signed_sum
ag_rna_negative_auc               # sum/mintegral of negative deltas
ag_rna_positive_auc
ag_rna_top5_abs_mean
ag_rna_affected_bin_count
ag_rna_peak_genomic_position
ag_rna_peak_track_index
ag_rna_peak_track_name
```

If VWF track/gene annotations are unavailable in the raw object, retain the full-window features and explicitly set `vwf_localization_status=not_available`; do not guess.

### Splicing

For `splice_sites`, and additionally `splice_site_usage` / `splice_junctions` when present:

```text
ag_splice_abs_max
ag_splice_gain_max                 # largest positive ALT-REF
ag_splice_loss_min                 # most negative ALT-REF
ag_splice_top5_abs_mean
ag_splice_affected_bin_count
ag_splice_peak_genomic_position
ag_splice_peak_track_index
ag_splice_peak_track_name
ag_splice_usage_gain_max
ag_splice_usage_loss_min
ag_splice_junction_gain_max
ag_splice_junction_loss_min
```

Where transcript/exon metadata are available, add distance to the nearest canonical VWF donor/acceptor and whether the peak is canonical, cryptic, or unlocalized.

## Required deliverables

```text
output/alphagenome_signed_features_2026-07-12.parquet
output/alphagenome_signed_features_2026-07-12.csv
output/alphagenome_signed_features_2026-07-12_manifest.json
output/alphagenome_signed_features_2026-07-12_qc_sample.csv
scripts/pipeline/extract_alphagenome_signed_features.py
```

The manifest must include formulas, sign convention, array axes, track filters, raw input checksum, row count, missingness by modality, and runtime.

## Mandatory validation before transfer

1. Select at least 10 variants spanning low/high current RNA and splice scores.
2. For each, manually verify one reported peak against the raw `ALT - REF` array.
3. Confirm `ag_rna_signed_min <= 0`, `ag_rna_signed_max >= 0`, and `abs_max = max(abs(min), abs(max))` within tolerance.
4. Confirm no subtype/HGMD/ClinVar label column was used by the extractor.
5. Report duplicate genomic keys and duplicate cDNA mappings.

## Transfer tonight

Preferred: upload only the five deliverables above to a temporary folder under the existing private/public dataset, for example:

```text
lucachangretta/VWF/alphagenome_signed_features_2026-07-12/
```

Alternative: `scp` the compact outputs directly. Do not upload the raw pickle.

After transfer, send back the exact paths/URLs and the manifest. The local evaluation will then rerun with signed RNA and localized splice features, using protein-position-held-out cross-validation and training-fold-only normalization.


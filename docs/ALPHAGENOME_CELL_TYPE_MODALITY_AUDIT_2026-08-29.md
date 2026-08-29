# AlphaGenome cell-type/modality audit for VWF

Date: 2026-08-29

## Conclusion

For VWF, the closest biologically interpretable replacement for the generic
endothelial-cell term `CL:0000115` is:

```text
CL:0002618  endothelial cell of umbilical vein (HUVEC-like)
```

It covers 8 of the 11 AlphaGenome human output types and is substantially more
complete than generic endothelial cell. It still lacks ATAC-seq and contact-map
tracks. No single AlphaGenome biosample currently covers all 11 output types;
the maximum observed coverage is 9/11 (HepG2, K562, and GM12878), but those are
not appropriate VWF-expression surrogates.

The 11 output types audited were:

```text
RNA_SEQ, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS,
DNASE, ATAC, CAGE, PROCAP, CHIP_HISTONE, CHIP_TF, CONTACT_MAPS
```

Recommended request policy:

1. Use `CL:0002618` as the primary, interpretable endothelial context for all
   modalities where it has tracks.
2. For modalities unavailable under `CL:0002618` (ATAC and contact maps), run a
   second request with `ontology_terms=None`, but label those columns as
   `all_tracks`, not endothelial-specific.
3. Preserve `alternate - reference` sign, genomic position, output type, track
   name, and ontology term in the compact feature table.
4. Never silently substitute zero for unavailable tracks; emit missingness and
   a reason code.

## Track coverage audited from AlphaGenome output metadata

| Ontology | Biosample | ATAC | CAGE | ChIP-Histone | ChIP-TF | Contact | DNase | PROcap | RNA | Splice junction | Splice-site usage | Covered |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `CL:0000115` | endothelial cell | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 2 | 1 | 2 | 5/11 |
| `CL:0002618` | endothelial cell of umbilical vein | 0 | 2 | 9 | 8 | 0 | 1 | 2 | 5 | 2 | 4 | 8/11 |
| `EFO:0001187` | HepG2 | 1 | 2 | 10 | 539 | 1 | 1 | 0 | 5 | 2 | 4 | 9/11 |
| `EFO:0002067` | K562 | 1 | 2 | 10 | 306 | 0 | 1 | 2 | 5 | 2 | 4 | 9/11 |
| `EFO:0002784` | GM12878 | 1 | 2 | 8 | 101 | 1 | 1 | 0 | 5 | 2 | 4 | 9/11 |

The track counts above are the number of returned tracks for each output type,
not a quality score. HepG2/K562/GM12878 have richer ChIP-TF coverage but are
not VWF-relevant endothelial biosamples.

## Feature-table implication

The compact feature schema should carry both source fields:

```text
ontology_term
source_scope = endothelial_huvec | all_tracks
```

This prevents mixing cell-specific and all-track evidence while preserving the
previously documented signed/localized AlphaGenome feature contract.

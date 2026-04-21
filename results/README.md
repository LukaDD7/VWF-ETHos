# Results Directory Structure

## Current Working Files (results/)

### AlphaGenome Inference Results
| File | Description | Size |
|------|-------------|------|
| `07_VCF_AlphaGenome_Results.csv` | Main results (2577 variants, 6 modalities for endothelial cell) | 320KB |
| `07_VCF_AlphaGenome_Results.pkl` | Main results pickle (full raw outputs) | 197GB |
| `07_VCF_AlphaGenome_Results_11_Modalities_Final.pkl` | Final merged 11 modalities | 197GB |

### 07e Supplementary Results (God Mode)
| File | Description | Status |
|------|-------------|--------|
| `07e_GodMode_Epigenome_Peaks.csv` | Supplementary 5 modalities (ATAC/Histone/CAGE/PROcap/Contact) | 505 records (145 success, 360 failed) |
| `07e_godmode_checkpoint.jsonl` | Checkpoint for resume | Active |
| `07e_task_list.json` | Task list (2577 variants) | Active |

### Phase 1-4 Pipeline Results
| File | Description |
|------|-------------|
| `01_filtered_vus.csv` | Filtered VUS variants |
| `03_inference_results.csv/pkl` | AlphaGenome inference results |
| `04_analysis_summary.csv` | Analysis summary |

### Final Merged Data
| File | Description |
|------|-------------|
| `Final_VWF_Lossless_Merge.xlsx` | Final merged dataset |
| `Final_VWF_Target_List_with_AlphaGenome.xlsx` | Target list with AlphaGenome |
| `Final_VWF_Target_List_with_AlphaGenome_FIXED.xlsx` | Fixed version |

## Backups (results/backups/)

| File | Description | Date |
|------|-------------|------|
| `07_VCF_AlphaGenome_Results_Backup.csv` | Backup of main CSV | 2026-03-16 |
| `07_VCF_AlphaGenome_Results_Backup.pkl` | Backup of main PKL | 2026-03-16 |
| `07e_GodMode_Epigenome_Peaks_backup_20260324_155250.csv` | Old 07e backup | 2026-03-24 |

## Modalities in 07_VCF_AlphaGenome_Results

### Complete (Endothelial Cell - CL:0000115)
- AG_RNA_SEQ: 2577/2577 samples
- AG_SPLICE_SITES: 2577/2577 samples
- AG_SPLICE_SITE_USAGE: 2577/2577 samples
- AG_SPLICE_JUNCTIONS: 2558/2577 samples
- AG_DNASE: 2577/2577 samples
- AG_CHIP_TF: 2577/2577 samples

### Missing (Need 07e Supplementary)
- AG_CAGE: 7/2577 samples
- AG_PROCAP: 7/2577 samples
- AG_ATAC: 7/2577 samples
- AG_CHIP_HISTONE: 7/2577 samples
- AG_CONTACT_MAPS: 7/2577 samples

## Scripts Using These Paths

- `scripts/07e_god_mode_epigenome_crawler_v3.py` - Main supplementary crawler (proxy protection)

Last updated: 2026-04-18
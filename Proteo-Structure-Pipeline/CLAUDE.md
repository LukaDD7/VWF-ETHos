# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Proteo-Structure Pipeline for VWF (von Willebrand Factor) protein 3D structure analysis using AlphaFold3 Server.

**Target Protein:** VWF (Uniprot: P04275, 2,813 amino acids)

## Pipeline Workflow

The pipeline has 3 sequential phases:

### Phase 1: Filter & Parse
Extract missense variants from Excel, generate WT and mutant FASTA files.

```bash
cd src
conda activate alphafold
python phase1_smart_filter.py \
    --excel /path/to/Final_VWF_Target_List_with_AlphaGenome.xlsx \
    --output-dir ../output
```

**Key Logic:**
- Filters `Molecular.consequence` for "missense variant"
- Excludes splice disruptors where `|Splice_REF正常概率 - Splice_ALT突变概率| > 0.3`
- Parses amino acid changes from multiple formats: 1-letter (G1531D), 3-letter (Val1409Phe), NM naming (p.Gly1531Asp)
- Downloads WT FASTA from Uniprot P04275
- Validates mutant sequences against WT before generating FASTA

### Phase 2: AlphaFold3 Batch JSON Generator
Generate individual JSON files per mutation, package into chunked ZIPs for manual upload.

```bash
# Test mode - first 5 mutations only
python phase2_af3_batch_generator.py --limit 5 --zip-chunk-size 10

# Production mode - all mutations (~1,306 tasks)
python phase2_af3_batch_generator.py --zip-chunk-size 10
```

**Critical Design Decision:**
- AlphaFold3 Server has no REST API - requires manual web upload
- Server reads only first task from JSON arrays, so we generate one JSON file per mutation
- ZIPs are chunked (default 10 tasks/ZIP) to match daily quota (~10-20 tasks/day)

**Output Structure:**
```
output/af3_batches/
├── individual_jobs/VWF_*.json     # 1,306 individual JSON files
├── af3_upload_part_001.zip        # 10 JSONs per ZIP (default)
├── af3_upload_part_002.zip
├── ... (~131 ZIPs total)
└── af3_jobs_manifest.json         # Complete job index
```

**JSON Format (per task):**
```json
{
  "name": "VWF_G1531D",
  "modelSeeds": [],
  "sequences": [{
    "proteinChain": {
      "sequence": "MIPARF...",
      "count": 1
    }
  }]
}
```

### Phase 3: Structural Scoring
RMSD calculation and pLDDT analysis using BioPython.

```bash
python phase3_structural_scoring.py \
    --wt-cif ../structures/predictions/VWF_WT.cif \
    --csv ../output/phase2_prediction_results.csv \
    --output ../output/phase3_final_scores.csv
```

**Metrics Computed:**
- `Global_RMSD`: Full protein Cα superposition RMSD
- `Local_RMSD_10A`: RMSD of atoms within 10Å radius of mutation site
- `pLDDT_Delta`: Difference in pLDDT at mutation site (Mut - WT)

## Environment Setup

```bash
conda create -n alphafold python=3.10 gxx_linux-64=11.* gcc_linux-64=11.* \
    gfortran_linux-64=11.* libstdcxx-ng libgcc-ng numpy pandas openpyxl \
    scipy biopython matplotlib seaborn -c conda-forge -y

conda activate alphafold
```

## Key Technical Details

### AlphaFold3 Server Constraints
- No public REST API - manual upload only
- Daily quota: ~10-20 tasks/day per account
- Input: Individual JSON files (not arrays)
- Output: CIF format structure files
- Sequence limit: 5,000 amino acids (VWF has 2,813, within limit)

### File Formats
- **Input:** Excel with variant data
- **Intermediate:** FASTA (WT and mutants), individual JSON files, chunked ZIPs
- **Output:** CIF files from AlphaFold3, CSV with RMSD/pLDDT scores

### Amino Acid Parsing
The parser handles three input formats found in the Excel:
1. `Protein.change` column: "G1531D" (1-letter)
2. `Name` column: "NM_000552.5(VWF):c.4592G>A (p.Gly1531Asp)" (NM naming)
3. `p.` column: "Val1409Phe" (3-letter)

### Code Architecture

**Phase 1:**
- `Phase1Filter` class: Main orchestrator
- `AminoAcidChange` dataclass: Stores parsed mutation info
- Methods: `load_excel`, `filter_missense_variants`, `filter_splice_disruptors`, `parse_amino_acid_change`, `download_wt_fasta`, `generate_mutant_fasta`

**Phase 2:**
- `AF3BatchJob` class: Single mutation task
- `AF3BatchGenerator` class: Orchestrates batch generation
- Key methods: `create_batch_jobs`, `save_individual_jsons`, `create_chunked_zip_packages`, `save_manifest`

**Phase 3:**
- `StructureAnalyzer` class: Loads WT structure, analyzes mutations
- `Phase3Runner` class: Main orchestrator
- `ScoringResult` dataclass: Stores computed metrics
- Key methods: `superimpose_and_calculate_rmsd`, `_get_atoms_in_radius`, `_extract_plddt`

## Manual Steps Required

AlphaFold3 results require manual intervention:

1. **Upload:** Visit https://alphafoldserver.com/, login with Google account
2. **Submit:** Upload one ZIP per day (`af3_upload_part_001.zip`, etc.)
3. **Wait:** 1-3 days per batch for computation
4. **Download:** Save CIF results to `structures/predictions/`
5. **Run Phase 3:** After all CIFs are downloaded

## Common Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--limit` | None | Test mode - process only N mutations |
| `--zip-chunk-size` | 10 | Tasks per ZIP (match daily quota) |
| `--splice-threshold` | 0.3 | Splice disruption cutoff |

---

## Table Format Standardization

### Standard Input Format (table_normalizer.py)

**目的**: 统一处理各种格式的 VWF 变异表格，自动检测格式并转换为标准输入格式。

**标准列**:
| 列名 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `AA_Position` | int | ✓ | 氨基酸位置 (1-2813) |
| `WT_AA_1` | str | ✓ | 野生型氨基酸 (1-letter, 如 'R') |
| `Mut_AA_1` | str | ✓ | 突变型氨基酸 (1-letter, 如 'W') |
| `Domain` | str | | 功能域 (如 'A1', 'D3') |
| `Original` | str | | 原始格式 (如 'Arg1341Trp') |

**支持的输入格式**:
1. **3-letter格式**: `"Tyr1258Cys"`, `"Pro1266Gln"`
2. **1-letter格式**: `"Y1258C"`, `"P1266Q"`
3. **分列格式**: `Position` + `WT_AA` + `Mut_AA`

**使用示例**:
```bash
cd src
conda activate alphafold

# 标准化表格并保存
python table_normalizer.py /path/to/2B_variants.xlsx \
    --output ../data/normalized_2B.csv

# 同时生成 AF3 JSON
python table_normalizer.py /path/to/2B_variants.xlsx \
    --af3-json ../data/af3_batch_2B.json \
    --wt-fasta ../structures/wt/VWF_P04275_WT.fasta \
    --exclude L1460F A1461V
```

**Python API**:
```python
from table_normalizer import TableNormalizer

# 标准化
normalizer = TableNormalizer("2B_variants.xlsx")
df = normalizer.normalize()
normalizer.validate()  # 检查错误
normalizer.save_csv("output.csv")

# 转换为 AF3 JSON
jobs = normalizer.to_af3_json(
    wt_fasta="VWF_P04275_WT.fasta",
    output_path="af3_batch.json",
    exclude_existing=["L1460F", "A1461V"]
)
```

---

## AlphaFold3 JSON Format

### AF3 Server Input Format

AlphaFold3 Server 接受 **JSON 数组格式**的批量上传文件:

```json
[
  {
    "name": "VWF_WT",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "MIPARFAGVLLALALILPGTLCAEG...",
          "count": 1
        }
      }
    ]
  },
  {
    "name": "VWF_Y1258C",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "MIPARF...Y...C...",
          "count": 1
        }
      }
    ]
  }
]
```

**关键要点**:
- 最外层是 **JSON 数组** (不是单个对象)
- 每个任务对象包含 `name` 和 `sequences`
- `proteinChain` 包含 `sequence` (突变后序列) 和 `count`
- 每文件最多 100 个任务 (Server 限制)

**从标准化表格生成 AF3 JSON**:
```bash
python table_normalizer.py input.xlsx \
    --af3-json output.json \
    --wt-fasta VWF_WT.fasta \
    --exclude L1460F A1461V
```

---

## Output Statistics

- Input variants: ~2,431 (from Excel)
- After missense filter: ~1,305
- Final batch jobs: 1,306 (1 WT + 1,305 mutants)
- ZIP files: ~131 (with default chunk size 10)
- Estimated upload time: ~131 days (1 ZIP/day)

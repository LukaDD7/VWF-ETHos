# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

VWF (von Willebrand Factor) variant analysis pipeline using AlphaGenome API for predicting variant effects on RNA expression and splicing.

## Pipeline Workflow

Run the 4-step pipeline sequentially:

```bash
# Step 1: Filter VUS and pathogenic variants from ClinVar+HGMD
python 01_filter_target_vus.py

# Step 2: Convert hg19->hg38 and build 1Mb intervals
python 02_preprocess_and_liftover.py

# Step 3: Run AlphaGenome inference (API calls)
python 03_run_alphagenome_inference.py

# Step 4: Visualize top variants
python 04_analyze_and_visualize.py
```

## Script Responsibilities

| Script | Input | Output |
|--------|-------|--------|
| `01_filter_target_vus.py` | Clinvar_HGMD_merge_annotated.xlsx | 01_filtered_targets.csv |
| `02_preprocess_and_liftover.py` | 01_filtered_targets.csv | 02_ready_for_inference.csv |
| `03_run_alphagenome_inference.py` | 02_ready_for_inference.csv | 03_inference_results.pkl, 03_inference_results.csv |
| `04_analyze_and_visualize.py` | 03_inference_results.pkl/csv | figures/*.png, 04_analysis_summary.csv |

## Configuration Notes

- **Script 03**: `API_KEY` 已配置，`MAX_VARIANTS = None` (处理全部变异)
- **Script 03**: `ONTOLOGY_TERMS = ['CL:0000115']` (内皮细胞，VWF 的表达细胞)
  - 注意：`UBERON:0002363`（血管内皮）无 RNA-seq tracks，已改用 CL:0000115
- **Script 04**: `TOP_N_VARIANTS = 20`，可视化 Delta Score 最高的 20 个变异
- **Script 04**: 4 面板临床级可视化（RNA-seq REF/ALT 叠加 + 差值，Splice Sites 4 条链 + 差值）

## Key Parameters

- Sequence length: 1Mb (2^20 = 1,048,576 bp) centered on variant
- API rate: ~5 秒/variant (实际速度，远快于预期 2 分钟)
- Delta score = max(RNA-seq delta_max, splice_sites delta_max)
- 可视化窗口：32kb (变异位点为中心，左右各 16kb)

## Test Files

- `test_alphagenome.py` - Quick API connectivity test
- `test_ontology.py` - Test different ontology terms and metadata inspection

## Dependencies

```bash
conda activate alphagenome
# pyliftover 已通过 conda 安装，chain 文件：hg19ToHg38.over.chain.gz (本地)
```

## Execution History (2026-03-10)

- **Script 01**: 过滤出 1198 个变异 (VUS + Pathogenic/Likely_pathogenic 阳性对照)
- **Script 02**: hg19→hg38 坐标转换完成，1Mb 区间构建完成
- **Script 03**: 1198/1198 变异预测成功 (100%)，耗时~98 分钟
  - 输出：`03_inference_results.pkl` (含 raw_outputs), `03_inference_results.csv`
  - Top 变异：chr12:6018550 (Delta: 76.81, Likely_pathogenic)
- **Script 04**: 生成 Top 20 变异的 4 面板临床级可视化图
  - 输出：`figures/*.png`, `04_analysis_summary.csv` (1198 条排行榜)

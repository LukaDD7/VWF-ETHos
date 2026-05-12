# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

VWF (von Willebrand Factor) variant analysis pipeline using AlphaGenome API for predicting variant effects on RNA expression and splicing.

## Architecture (2026-04-21)

### Current Architecture: Alphagenome-DFR-Phase3

Reference: `Alphagenome-DFR-Phase3.md`

**Phase 1 (COMPLETED)**: Data Integration
- `merge_alpha_features.py` - Creates `VWF_Alpha_Matrix.parquet`
- AG matching: 67% via cDNA position extraction from consequence strings

**Phase 2 (COMPLETED)**: AgenticVWFClassifier
- `agentic_vwf_classifier.py` - 3 Expert Agent architecture
- Expert 1: StructuralExpert (AF3 pLDDT/PAE)
- Expert 2: TranscriptomicExpert (AG RNA/splice delta)
- Expert 3: ClinicalGeneticistAgent (Logical Fusion)

**Phase 2 Validation Results**:
| Known ↓ / Predicted → | 2A | 2B | 2M | 2N | uncertain |
|---|---|---|---|---|---|
| **2A** (43) | 26 | 0 | 4 | 2 | 11 |
| **2B** (12) | 0 | 2 | 8 | 0 | 2 |
| **2M** (25) | 0 | 0 | 24 | 0 | 1 |
| **2N** (20) | 7 | 0 | 0 | 13 | 0 |

**Remaining Tasks** (per user instructions):
- Task 2: Upgrade Expert 1 for Type 2B (AF3 interface PAE analysis)
- Task 3: Execute validation on 59 AF3 variants

## Pipeline Workflow (Legacy)

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

## Key Files

| File | Purpose |
|------|---------|
| `VWF_Alpha_Matrix.parquet` | Unified feature matrix (100 variants) |
| `VWF_Alpha_Matrix_classified.parquet` | Classification results |
| `merge_alpha_features.py` | Phase 1 data integration |
| `agentic_vwf_classifier.py` | Phase 2 classifier |

## Configuration Notes

- **Script 03**: `API_KEY` 已配置，`MAX_VARIANTS = None` (处理全部变异)
- **Script 03**: `ONTOLOGY_TERMS = ['CL:0000115']` (内皮细胞，VWF 的表达细胞)
- **Script 04**: `TOP_N_VARIANTS = 20`，可视化 Delta Score 最高的 20 个变异

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

## GPU Pipeline Index

All GPU scripts follow the same pattern: parallel workers (1 job/GPU), `.done` markers for resume, preflight checks, offline-only design.

**GPU 实例约束**: 不可联网（仅 CPU 实例可联网），/dev/shm 仅 64MB，计费制。

### Boltz-2: A1 + GPIbα (已完成)

- **Runner**: `scripts/pipeline/run_a1_gpiba_boltz2.sh`
- **YAML 生成**: `scripts/pipeline/generate_a1_gpiba_yamls.py`
- **结果解析**: `scripts/pipeline/parse_a1_gpiba_results.py`, `scripts/pipeline/pae_a1_gpiba_analysis.py`
- **输入**: 74 YAML → `output/boltz2_a1_gpiba/`
- **输出**: `output/boltz2_a1_gpiba_results/` (iPTM/PAE → `output/boltz2_a1_gpiba_analysis/`)
- **用法**: `bash scripts/pipeline/run_a1_gpiba_boltz2.sh --gpus 4`

### Boltz-2: VWD Functional Panel (含 FVIII 2N axis)

- **Runner**: `scripts/pipeline/run_vwd_functional_boltz2_panel.sh`
- **YAML 生成**: `scripts/pipeline/generate_vwd_functional_boltz2_yamls.py`
- **结果解析**: `scripts/pipeline/parse_vwd_functional_boltz2_results.py`
- **输入**: 990 YAML (15 assay × 597 variants + 15 WT baselines) → `output/boltz2_vwd_functional_panel/yamls/`
- **输出**: `output/boltz2_vwd_functional_panel/boltz_results/`
- **关键 assay**: `dprime_d3_fviii_binding` (100 jobs, 解决 2N 分类), `a1_gpiba_forced_binding` (120 jobs)
- **用法**:
  ```bash
  bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --gpus 8           # 全部 990 jobs
  bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --gpus 4 --filter fviii  # 只跑 FVIII (100 jobs)
  bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --preflight        # 预检
  ```
- **算力估算**: 990 jobs / 8 GPU ≈ 41h; 100 FVIII jobs / 4 GPU ≈ 8h

### GROMACS: VWF A1–GPIbα MD Simulation

- **Runner**: `scripts/pipeline/run_gromacs_vwf_md.sh`
- **MDP 参数**: 自动生成至 `output/gromacs_md/mdp/`（或手动放 `scripts/pipeline/mdp/`）
- **输入**: Boltz-2 CIF → `output/boltz2_a1_gpiba_results/boltz_results_VWF_*/predictions/model_0.cif`
- **输出**: `output/gromacs_md/VWF_<variant>/` (em → nvt → npt → production → analysis)
- **依赖**: `conda install -c conda-forge gromacs=2025; pip install gemmi gmx_MMPBSA`
- **用法**:
  ```bash
  bash scripts/pipeline/run_gromacs_vwf_md.sh --gpus 8                        # 全部
  bash scripts/pipeline/run_gromacs_vwf_md.sh --gpus 4 --filter R1306W        # 单个突变
  bash scripts/pipeline/run_gromacs_vwf_md.sh --phase equil                   # 只跑到平衡
  bash scripts/pipeline/run_gromacs_vwf_md.sh --ns 500                        # 500 ns production
  ```
- **性能**: H200 单卡 ~200 ns/day (65K atoms), 8 GPU 并行 8 个突变体

### CIF → 二次结构特征提取

- **脚本**: `scripts/pipeline/compute_structural_features.py`
- **输出**: `structural_features.csv` / `.parquet`
- **特征**: SASA, ΔSASA, BSA, contacts, H-bonds, salt bridges, mutation-site RSA, Rg, confidence
- **依赖**: `pip install biopython freesasa`
- **用法**:
  ```bash
  python scripts/pipeline/compute_structural_features.py \
      --boltz-results output/boltz2_a1_gpiba_results \
      --output output/structural_features.csv
  ```

### Evidence Matrix Builder

- **脚本**: `scripts/pipeline/build_vwd_evidence_matrix.py`
- **输入**: functional panel 解析结果
- **输出**: `output/boltz2_vwd_functional_panel/evidence_matrix.csv`

### GPU Instance Checklist

```bash
# 1. CPU 实例（可联网）— 环境准备
conda create -n boltz2 python=3.11 && conda activate boltz2
pip install boltz>=2.0 biopython freesasa gemmi torch
# 打包: tar czf env_boltz2.tar.gz -C ~/miniconda3/envs/boltz2 .

# 2. 传输到 GPU 实例
scp env_boltz2.tar.gz gpu-instance:/workspace/
scp -r VWF-ETHos/ gpu-instance:/workspace/

# 3. GPU 实例（无网络）— 运行
conda activate boltz2
export CC=$(which gcc)
rm -f core.*
bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --gpus 8 --preflight
bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --gpus 8

# 4. Triton JIT 缓存复用
# 首次运行后: scp -r ~/.cache/triton/ next-gpu-instance:~/.cache/
```

### Boltz-2 Troubleshooting (通用)

| Error | Fix |
|-------|-----|
| `No supported gpu backend found` | `--num_workers 0` |
| `Failed to find C compiler` | `export CC=$(which gcc)` |
| `OSError: No space left on device` | `--num_workers 0` (避免 /dev/shm IPC) |
| `core.*` files | `rm -f core.*`; 根因通常是 /dev/shm |
| GPU util 4% | 正常 (diffusion sampling 串行) |
| GPU util 0% | 等待 Triton JIT 编译 (首次 5-10 min) |

## Execution History

### 2026-05-12
- **Functional panel runner upgraded**: `run_vwd_functional_boltz2_panel.sh` — 并行 worker + .done + preflight + assay filter
- **GROMACS MD pipeline added**: `run_gromacs_vwf_md.sh` — CIF→PDB→EM→NVT→NPT→Production
- **Structural features extractor added**: `compute_structural_features.py` — CIF→SASA/BSA/contacts

### 2026-05-07
- **A1+GPIbα Boltz-2 pipeline**: 74 variants completed
- **Issues fixed**: /dev/shm 64MB, gcc, num_workers

### 2026-03-10
- **AlphaGenome pipeline**: 1198/1198 variants predicted (100%)

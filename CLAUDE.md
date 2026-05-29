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
| `scripts/pipeline/run_gromacs_vwf_md.sh` | GROMACS MD (complex/monomer/autoinhib) |
| `scripts/pipeline/generate_a1_dp_d3_yamls.py` | A1+D'D3 YAML 生成 (autoinhibition) |
| `scripts/pipeline/run_a1_dp_d3_boltz2.sh` | A1+D'D3 Boltz-2 runner |
| `scripts/pipeline/analyze_gromacs_md.py` | GROMACS 轨迹后分析 |
| `scripts/pipeline/compute_structural_features.py` | CIF → SASA/BSA/contacts 特征提取 |

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

### GROMACS: VWF A1–GPIbα MD Simulation (支持三种系统)

> ⚠️ **当前使用 CHARMM27 临时方案**（CHARMM36m 有 terminal patching 兼容性问题，待解决）

- **Runner**: `scripts/pipeline/run_gromacs_vwf_md.sh`
- **MDP 参数**: 自动生成至 `output/gromacs_md/mdp/`（或手动放 `scripts/pipeline/mdp/`）
- **输入**（按系统）:
  - `complex`: `output/boltz2_a1_gpiba_results/boltz_results_VWF_*/predictions/`
  - `monomer`: 同上（提取 chain A）
  - `autoinhib`: `output/boltz2_a1_dp_d3_results/` (A1+D'D3)
- **输出**: `output/gromacs_md_{system}/VWF_<variant>/` (em → nvt → npt → production → analysis)
- **依赖**: `conda install -c conda-forge gromacs=2025; pip install gemmi`; `gmx_mmpbsa` (可选)
- **CHARMM27/36m**: 当前用 charmm27（见 Troubleshooting）
- **用法**:
  ```bash
  # 三种系统:
  bash scripts/pipeline/run_gromacs_vwf_md.sh --system complex --gpus 8    # A1-GPIbα (默认)
  bash scripts/pipeline/run_gromacs_vwf_md.sh --system monomer --gpus 8   # A1 单体（提取 chain A）
  bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --gpus 8  # A1+D'D3（需先跑 Boltz-2）

  bash scripts/pipeline/run_gromacs_vwf_md.sh --gpus 4 --filter R1306W        # 单个突变
  bash scripts/pipeline/run_gromacs_vwf_md.sh --phase equil                   # 只跑到平衡
  bash scripts/pipeline/run_gromacs_vwf_md.sh --ns 500                        # 500 ns production
  ```

**三系统分析目标**:
| 系统 | 研究问题 | 主要指标 |
|------|---------|---------|
| `complex` | 突变如何影响 A1↔GPIbα 结合 | MM/PBSA ΔG_bind |
| `monomer` | 突变对 A1 自身柔性的影响 | RMSF, PCA |
| `autoinhib` | 突变如何破坏 D'D3↔A1 自抑制接口 | D'D3-A1 contacts (H-bond, salt bridge) |
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
# GROMACS 环境
conda create -n gromacs python=3.11 && conda activate gromacs
conda install -c conda-forge gromacs=2025 gmx_mmpbsa && pip install gemmi
conda install -c conda-forge gcc_linux-64=13 gxx_linux-64=13  # gcc for runtime
# 打包: tar czf env_boltz2.tar.gz -C ~/miniconda3/envs/boltz2 .
# 打包: tar czf env_gromacs.tar.gz -C ~/miniconda3/envs/gromacs .

# 2. CHARMM36m force field（可选，见 Troubleshooting）
# 从 MacKerell 下载 charmm36-feb2026_cgenff-5.0.ff.tgz 到 gromacs/share/gromacs/top/

# 3. 传输到 GPU 实例
scp env_boltz2.tar.gz gpu-instance:/workspace/
scp env_gromacs.tar.gz gpu-instance:/workspace/
scp -r VWF-ETHos/ gpu-instance:/workspace/

# 4. GPU 实例（无网络）— 运行
conda activate boltz2
export CC=$(which gcc)
rm -f core.*
bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --gpus 8 --preflight
bash scripts/pipeline/run_vwd_functional_boltz2_panel.sh --gpus 8

# 5. GROMACS MD（三系统）
conda activate gromacs
export LD_LIBRARY_PATH=.../shared_libs:$LD_LIBRARY_PATH
export CC=$(which gcc)
bash scripts/pipeline/run_gromacs_vwf_md.sh --system complex --gpus 8
bash scripts/pipeline/run_gromacs_vwf_md.sh --system monomer --gpus 8
# Autoinhibition 需先跑 A1+D'D3 Boltz-2:
python3 scripts/pipeline/generate_a1_dp_d3_yamls.py
bash scripts/pipeline/run_a1_dp_d3_boltz2.sh --gpus 8
bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --gpus 8

# 6. Triton JIT 缓存复用
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

## GROMACS MD Pipeline Troubleshooting

### CHARMM36m force field not found

**错误**: `Could not find force field 'charmm36m' in current directory`

**原因**: GROMACS 2025.4 conda 包的 `charmm36m.ff` 目录为空（安装不完整）

**解决**:
```bash
# CPU 实例下载（GPU 实例无网络）
cd /tmp
curl -sL "https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-feb2026_cgenff-5.0.ff.tgz" -o charmm36m.tar.gz
tar xzf charmm36m.tar.gz
cp -r charmm36-feb2026_cgenff-5.0.ff <envs>/gromacs/share/gromacs/top/charmm36m.ff
```
force field 会通过 NFS 共享到 GPU 实例

### Boltz-2 PDB → pdb2gmx CHARMM36m 兼容性

**错误**: `atom C1 not found in building block 1MET while combining tdb and rtp`

**根本原因**: Boltz-2 输出的 PDB 缺少 CHARMM36m 所需的终端加帽原子：
- **N 端**: 需要 N1/N2 加帽原子（Boltz-2 输出只有 N）
- **C 端**: 需要 OXT 原子（Boltz-2 输出只有 C, O, CA）

CHARMM27 对此更宽松（可接受），但 CHARMM36m 严格要求。

**已验证可行的组合**:
- `charmm27` + `-ignh` → ✅ 成功（用于当前 74 个突变体）
- `charmm36m` + `-ter` + `chainsep id` → ❌ N 端 C1 仍然报错

**推荐方向**（待实现）:
1. **修改 Boltz-2 YAML 配置**，让 Boltz-2 预测时输出更完整的 PDB（带 terminal patches）
2. **用 OpenMM + CHARMM36m XML** 做结构预处理（`boltz2` 环境已有 `charmm36_2024.xml`），输出标准化的 PDB
3. **GROMACS `pdb2gmx` 预处理脚本**: 写一个 `preprocess_for_pdb2gmx.py`，用 gemmi 给每个链的 N/C 端加上所需原子（参考 `scripts/pipeline/preprocess_for_pdb2gmx.py` 的框架，但 gemmi API 调用需修复）

**相关文件**:
- `scripts/pipeline/preprocess_for_pdb2gmx.py` — 框架脚本（OXT 添加逻辑已验证，但 chain splitting 的 gemmi API 调用需调试）
- `scripts/pipeline/run_gromacs_vwf_md.sh` — 目前用 charmm27 临时方案（行 432: `-ff charmm27`）

## Execution History

### 2026-05-29
- **GROMACS MD pipeline expanded**: Three-system architecture (complex / monomer / autoinhib)
  - `generate_a1_dp_d3_yamls.py` — A1+D'D3 YAML 生成 (D' + D3 + A1, 669 aa)
  - `run_a1_dp_d3_boltz2.sh` — A1+D'D3 Boltz-2 runner
  - `run_gromacs_vwf_md.sh --system complex|monomer|autoinhib` — 扩展支持三种系统
  - `analyze_gromacs_md.py` — 统一后分析脚本
- **CHARMM36m installed**: 从 MacKerell 下载 `charmm36-feb2026_cgenff-5.0.ff` 到 `lzy/envs/gromacs/share/gromacs/top/charmm36m.ff/`
- **Issue**: Boltz-2 PDB 缺少 terminal OXT/C1 原子，CHARMM36m 兼容性失败；当前用 charmm27 临时绕过

### 2026-05-12
- **Functional panel runner upgraded**: `run_vwd_functional_boltz2_panel.sh` — 并行 worker + .done + preflight + assay filter
- **GROMACS MD pipeline added**: `run_gromacs_vwf_md.sh` — CIF→PDB→EM→NVT→NPT→Production
- **Structural features extractor added**: `compute_structural_features.py` — CIF→SASA/BSA/contacts

### 2026-05-07
- **A1+GPIbα Boltz-2 pipeline**: 74 variants completed
- **Issues fixed**: /dev/shm 64MB, gcc, num_workers

### 2026-03-10
- **AlphaGenome pipeline**: 1198/1198 variants predicted (100%)

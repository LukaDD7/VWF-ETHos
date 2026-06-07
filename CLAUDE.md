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

### GPU 实例部署关键事实 (2026-06-01)

1. **NFS 共享 env**: `/lzy/envs/gromacs` (5.1 GB) 已在 GPU 实例上, **无需 scp**。`/lzy/` 是 gpfs_hdd 并行文件系统
2. **GROMACS 2025.4-conda_forge = OpenCL backend, 不是 CUDA**。`gmx mdrun -version` 输出:
   - `GPU support: OpenCL`
   - `OpenCL library: /lzy/envs/gromacs/lib/libOpenCL.so` (env 自带)
   - 0 个 CUDA API 符号 (nm -D libgromacs.so.10 | grep cuda → 无)
3. **GPU 实例需 NVIDIA driver 自带的 OpenCL ICD**:
   - `/usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.<ver>` (driver 装好就有)
   - 不需要 libcudart / libcuda.so.1, 不需要 conda install cuda-*
4. **不需要 conda activate**: 项目 runner (`run_gromacs_vwf_md.sh`) 用绝对路径 `${GMX:=/lzy/envs/gromacs/bin.AVX2_256/gmx}`, GPU 实例直接跑即可
5. **Project 端 patched FF**: `force_fields/charmm36m.ff/aminoacids.n.tdb` 注入了 17 个 protein `<RESNAME>1` patches (解决 pdb2gmx 1MET bug), 通过 `GMXLIB=force_fields` + cwd symlink 让 gmx 优先用 patched copy
6. **gemmi 已在 gromacs env**: `/lzy/envs/gromacs/bin/python -c "import gemmi"` → 0.7.5 ✅

### GPU 实例一行命令开跑

```bash
# 在 GPU 实例上 (H200, NVIDIA driver 570+):
cd /lzy/projects/VWF-ETHos
bash scripts/pipeline/run_gromacs_vwf_md.sh --system complex --gpus 8 --ns 200
# 74 变体 × 200ns @ H200 ≈ 9 天
```

**Preflight 自动验证** (脚本里):
- `gmx` 绝对路径解析
- GPU backend = OpenCL
- NVIDIA OpenCL ICD 存在 (`/usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.*`)
- gemmi 在 gromacs env 可用
- nvidia-smi 检测 GPU 数量

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

### GPU Instance Checklist (2026-06-01 更新: GROMACS MD 简化)

```bash
# 0. GPU 实例先验证 (无需 conda activate)
nvidia-smi                                                # 列出 H200
ls /usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.*        # 应有 1 + 版本号
envs/gromacs/bin.AVX2_256/gmx mdrun -version 2>&1 | head -25
# 应输出: GPU support: OpenCL, OpenCL library: /lzy/envs/gromacs/lib/libOpenCL.so

# 1. CPU 实例 (可联网) — 环境准备 (一次性, 已完成)
conda create -n boltz2 python=3.11 && conda activate boltz2
pip install boltz>=2.0 biopython freesasa gemmi torch
conda create -n gromacs python=3.11 && conda activate gromacs
conda install -c conda-forge gromacs=2025 gmx_mmpbsa && pip install gemmi
conda install -c conda-forge gcc_linux-64=13 gxx_linux-64=13

# 2. CHARMM36m force field (已下载, 不需再做)
# /lzy/envs/gromacs/share/gromacs/top/charmm36m.ff/  (env 自带)
# /lzy/projects/VWF-ETHos/force_fields/charmm36m.ff/  (project patched copy)

# 3. **不需要 scp 到 GPU 实例**! /lzy/ 是 NFS 共享 (gpfs_hdd)
# GPU 实例上直接看到 /lzy/envs/gromacs/

# 4. GROMACS MD (三系统) — 一行命令开跑
cd /lzy/projects/VWF-ETHos
bash scripts/pipeline/run_gromacs_vwf_md.sh --preflight --system complex
bash scripts/pipeline/run_gromacs_vwf_md.sh --system complex --gpus 8 --ns 200
# 74 变体 × 200ns @ 8×H200 ≈ 9 天

# 5. monomer / autoinhib 同理
bash scripts/pipeline/run_gromacs_vwf_md.sh --system monomer --gpus 8
bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --gpus 8
# autoinhib 需先跑 A1+D'D3 Boltz-2
python3 scripts/pipeline/generate_a1_dp_d3_yamls.py
bash scripts/pipeline/run_a1_dp_d3_boltz2.sh --gpus 8

# 6. Boltz-2 (Triton JIT 缓存复用, 旧方法已废, 见下)
# 旧: scp -r ~/.cache/triton/ next-gpu-instance:~/.cache/
# 新: 不需要, GROMACS OpenCL kernel 缓存在 gromacs env 内, NFS 共享
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

### Boltz-2 PDB → pdb2gmx CHARMM36m 兼容性（2026-06-01 解决）

**错误**: `atom C1 not found in building block 1MET while combining tdb and rtp`

**真正的根本原因** (之前误诊为 N 端缺 N1/N2、C 端缺 OXT — 实际上 `-ignh` 模式下 pdb2gmx 自己的 `[ NH3+ ]` / `[ COO- ]` patch 会加 HC/OT1/OT2, 不需要预加):

- **charmm2gmx 端口不完整**：`gromacs-2025.4` 自带的 `charmm36m.ff/aminoacids.n.tdb` 只定义了通用的 `[ NH3+ ]` / `[ GLY-NH3+ ]` / `[ PRO-NH2+ ]` patches，**没有 per-residue `<RESNAME>1` patches**（MET1、ALA1、…）。
- **pdb2gmx 选择逻辑**：遇到 N-terminal 残基时先查 `<RESNAME>1` patch。protein MET1 在 aminoacids.n.tdb 中找不到，但 `ethers.n.tdb` 里有 `[ MET1 ]`（用于醚类化学）。pdb2gmx 错误地用 ether patch 处理 protein MET, 而 ether patch 引用 `C1` 原子, 在 MET rtp 中找不到 → 报错。
- **charmm27 之所以能跑通**：charmm27.ff 没有 ethers.n.tdb, 所以 `[ MET1 ]` 完全不存在, pdb2gmx 直接 fall back 到 `[ NH3+ ]`。

**解决方案（已实施）**:
1. 复制 `charmm36m.ff` 到 `force_fields/`
2. 在 `aminoacids.n.tdb` 顶部注入 protein 兼容的 `[ MET1 ] [ ALA1 ] [ VAL1 ] [ LEU1 ] [ ILE1 ] [ PHE1 ] [ TRP1 ] [ SER1 ] [ THR1 ] [ CYS1 ] [ TYR1 ] [ ASN1 ] [ GLN1 ] [ ASP1 ] [ GLU1 ] [ ARG1 ] [ LYS1 ] [ HIS1 ]` patches, 全部等同于 `[ NH3+ ]` (改 N atomtype 为 NH3, 加 3 个 HC)
3. `run_gromacs_vwf_md.sh` 在 `pdb2gmx` 步骤前:
   - `ln -sf $ROOT_DIR/force_fields/charmm36m.ff $WORK/topology/charmm36m.ff`
   - `cd $WORK/topology`
   - `GMXLIB=$ROOT_DIR/force_fields gmx pdb2gmx … -ff charmm36m -ignh`
   - 让 gmx 优先在当前目录 + GMXLIB 找到 patched FF, 而不是 conda env 内置的
4. **删除**了原有的 OXT 预添加步骤 (gemmi 0.6.5 没有 `Atom.x` 属性, 而且 `[ COO- ]` patch 已经会 `delete O` + `add OT1/OT2`)

**验证** (R1306W 变体, complex 系统, 489 残基 7827 atoms):
- `gmx pdb2gmx` ✅
  - Chain A (HIS → PRO): `Start terminus HIS-1: NH3+` / `End terminus PRO-199: COO-`
  - Chain B (MET → ASP): `Start terminus MET-1: MET1` (用我们注入的 patch) / `End terminus ASP-290: COO-`
  - Total mass 32086.404 a.m.u., total charge -8.000 e
- `gmx editconf` ✅ box vectors 10.712 nm
- `gmx solvate` ✅ 26077 SOL
- `gmx grompp + genion` ✅ 83 NA + 79 CL
- `gmx mdrun -steep` ✅ 466 steps, Epot -1.3306265e+06 kJ/mol, Fmax 9.1e+02 < 1000

**相关文件**:
- `force_fields/charmm36m.ff/aminoacids.n.tdb` — 注入的 [ MET1 ] 等 patches
- `scripts/pipeline/run_gromacs_vwf_md.sh` — pdb2gmx 步骤已用 patched FF
- `scripts/pipeline/preprocess_for_pdb2gmx.py` — 旧的 OXT 预添加脚本, 已废弃 (但保留以防未来有其它需要)

## Execution History

### 2026-06-07
- **GPU flags 按后端条件化 (修活雷)**：`run_gromacs_vwf_md.sh` 之前在 OpenCL 构建上硬用 `-bonded gpu -update gpu` + `GMX_CUDA_GRAPH=1`(CUDA/SYCL 专属)→ EM 能过、NVT 必 fatal、白烧机时。现按 `$GPU_BACKEND` 选 flags:CUDA/SYCL→全 GPU 常驻;OpenCL/未知→`-nb gpu -pme gpu`(update/bonded 回 CPU)。
- **`gpu_smoke_test.sh`**：诚实放行闸,用与 runner 相同的条件化 flags + md 积分器 tpr 实跑 5 步(Test B),`-nb gpu` only 的朴素 smoke 会假性放行。
- **2B 自抑制特征接入分类器**：
  - 实测 `evidence_matrix.csv` 的 `a1_aim_autoinhibition_context` 全局 `ptm_or_plddt` **分不开 2B/2M**(delta 中位数 0.067 vs 0.049, 全重叠)。
  - 新增 `extract_aim_autoinhib_features.py`：从自抑制 CIF 抽 **AIM↔A1 接触数** → `aim_release_score`(越大=AIM 脱离=越像 2B),铺 ~130 变体。
  - `agentic_vwf_classifier.py` RULE6 加性接入(NaN 退回原逻辑):强松开早返回 2B;弱松开救回本判 2M 的 A1 默认。**阈值 `AIM_RELEASE_2B_Z` 为暂定值,须在集群用已知 2B/2M 校准**。
- **`docs/A40_RUNBOOK.md`**：独立 A40 机(非 /lzy NFS)的搬运 + 重建 CUDA gromacs env + smoke + autoinhib MD + 2B 特征流水线。
- ⚠ **待办**:① A40 上 `gpu_smoke_test.sh` 必须 PASS 再上批量;② 集群上跑 extract + 校准 `AIM_RELEASE_2B_Z`,重跑分类器看 2B recall。

### 2026-06-01
- **CHARMM36m 兼容性 bug 修复**：`charmm2gmx` 端口的 `aminoacids.n.tdb` 缺 per-residue `<RESNAME>1` patches, pdb2gmx 错误地从 `ethers.n.tdb` 选到 ether 的 `[ MET1 ]` patch
- **解决方案**：在 `force_fields/charmm36m.ff/` 注入 17 个 protein N-terminal patches (MET1, ALA1, …, HIS1), run_gromacs_vwf_md.sh 通过 GMXLIB + cwd symlink 优先用 patched FF
- **删除** OXT 预添加步骤（gemmi 0.6.5 Atom.x 不存在 + [ COO- ] patch 自处理）
- **端到端验证**：R1306W complex (A1+GPIbα, 489 残基 7827 atoms) 跑通 pdb2gmx → editconf → solvate → genion → grompp → mdrun (EM 466 steps, Fmax < 1000)
- **GPU 部署准备**：
  - 确认 gromacs 2025.4 = OpenCL backend (无需 CUDA runtime)
  - `run_gromacs_vwf_md.sh` 改用 `${GMX}` 绝对路径, 不依赖 conda activate
  - Preflight 自动检测 OpenCL ICD / GPU backend / gemmi in gromacs env
  - 整个 gromacs env (5.1GB) 通过 NFS (`/lzy/` gpfs_hdd) 共享, GPU 实例无需 scp
  - 性能: CPU 16-core ≈ 21 ns/day vs H200 ≈ 200 ns/day → 9× per-node 加速, 8 GPU 并行 9 天跑完 74 变体

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

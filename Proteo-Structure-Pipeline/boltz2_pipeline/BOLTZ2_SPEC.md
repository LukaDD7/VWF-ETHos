# VWF Boltz-2 Pipeline — 规范文档
# Proteo-Structure-Pipeline / boltz2_pipeline

## 目录结构

```
boltz2_pipeline/
├── vwf_ligand_database.py      # 配体库（本地硬编码，无 API 调用）
├── prepare_boltz2_inputs.py    # 表格预处理 + JSON 批次生成
├── run_boltz2.sh               # GPU 服务器执行脚本（4x H200）
├── parse_boltz2_results.py     # 结果解析 + 热图 + 分型建议
└── BOLTZ2_SPEC.md              # 本文件
```

---

## 核心设计原则：无上帝视角（盲扫诊断模式）

> 传统流程（错误）：已知分型 → 选配体 → 验证
>
> 本流程（正确）：不预设分型 → 该结构域的全部配体都跑 → 用亲和力差异反推分型

---

## 配体库 (`vwf_ligand_database.py`)

### 注册配体一览

| Key | 蛋白质 | UniProt | 链数 | 对应结构域 | 对应分型 |
|---|---|---|---|---|---|
| `GPIb_alpha` | 血小板 GPIbα | P07359 | 1 | A1 | 2B (GOF) / 2M-A1 (LOF) |
| `Collagen_I_THP` | 胶原 I 三股螺旋肽 | P02452 | **3** | A3, C1, C2 | 2M-A3 |
| `ADAMTS13_Spacer` | ADAMTS13 Spacer | Q76LX8 | 1 | A2 | 2A |
| `FVIII_LightChain` | FVIII 轻链 C1 | P00451 | 1 | D', D3 | 2N |
| `Heparan_Sulfate_mimic` | HS 结合模拟肽 | N/A | 1 | A1, D' | 调节性 |
| `Integrin_alphaIIb_beta3` | αIIbβ3 integrin | P08514 | 1 | C4 | 2M (C4) |

### 结构域 → 配体自动映射

| 结构域 | 自动关联配体 |
|---|---|
| D' / D3 | FVIII_LightChain |
| A1 | GPIb_alpha, Heparan_Sulfate_mimic |
| A2 | ADAMTS13_Spacer |
| A3 | Collagen_I_THP (3链) |
| C1/C2 | Collagen_I_THP |
| C4 | Integrin_alphaIIb_beta3 |

---

## 流水线运行步骤

### Step 1：表格预处理 + JSON 生成

```bash
cd boltz2_pipeline

python prepare_boltz2_inputs.py \
    --excel ../../results/Final_VWF_Target_List_with_AlphaGenome_FIXED.xlsx \
    --wt-fasta ../structures/wt/VWF_P04275_WT.fasta \
    --output-dir ../../output/boltz2_blind_scan \
    --batch-size 30

# 测试模式（先跑前 20 个变异）
python prepare_boltz2_inputs.py \
    --excel ../../results/Final_VWF_Target_List_with_AlphaGenome_FIXED.xlsx \
    --wt-fasta ../structures/wt/VWF_P04275_WT.fasta \
    --output-dir ../../output/boltz2_blind_scan_test \
    --batch-size 30 --limit 20
```

**输出：**
- `output/boltz2_blind_scan/batch_001.json` ... `batch_NNN.json`
- `output/boltz2_blind_scan/job_manifest.csv`

**job 命名规则：** `VWF_{ref}{pos}{alt}_vs_{ligand_key}`
例：`VWF_R1306W_vs_GPIb_alpha`、`VWF_W1745C_vs_Collagen_I_THP`

---

### Step 2：Boltz-2 结构 + 亲和力预测（GPU 服务器）

```bash
# 传输 batch JSON 到 GPU 服务器
scp -r output/boltz2_blind_scan user@h200-server:/workspace/vwf/

# 在 GPU 服务器上
cd /workspace/vwf/boltz2_pipeline
chmod +x run_boltz2.sh
./run_boltz2.sh
```

**增量存储机制（防断电 / 防崩溃）：**

脚本内部自动将每个 batch JSON 展开为**单个 job JSON 文件**（存入 `_tmp_single_jobs/`），然后逐 job 执行。每个 job 完成后立即写入：

```
output/boltz2_results/VWF_R1306W_vs_GPIb_alpha/.done   ← 时间戳标记
```

**断点续跑只需重新执行：**
```bash
./run_boltz2.sh   # 自动跳过所有已有 .done 的 job
```

**清理崩溃脏文件并重试单个 job：**
```bash
rm -rf output/boltz2_results/VWF_R1306W_vs_GPIb_alpha
./run_boltz2.sh
```

**实时监控（另一个终端）：**
```bash
# 单次快照
python monitor_progress.py --results-dir ../../output/boltz2_results

# 持续刷新（30秒间隔）
python monitor_progress.py --results-dir ../../output/boltz2_results --watch

# 示例输出：
# ============================================================
#   VWF Boltz-2 Progress Monitor
# ============================================================
#   [████████████████░░░░░░░░░░░░░░░░░░░░░░░░]  40.0%
#   Completed :   800 / 2000
#   Remaining :  1200
#   Speed          : 28.3 jobs/hour
#   ETA            : 42.4 hours
# ============================================================
```

```
N_GPUS=4                # 卡数
RECYCLING_STEPS=3       # 推荐 3（平衡精度/速度）
DIFFUSION_SAMPLES=5     # 每 job 生成 5 个构象，取最优 pLDDT
NUM_WORKERS=8           # H200 推荐
```

**Boltz-2 核心命令：**
```bash
boltz predict batch_001.json \
    --out_dir results/batch_001 \
    --accelerator gpu \
    --devices 4 \
    --strategy ddp \
    --recycling_steps 3 \
    --diffusion_samples 5 \
    --num_workers 8
```

**输出结构：**
```
output/boltz2_results/
  batch_001/
    VWF_R1306W_vs_GPIb_alpha/
      predictions/
        model_0.cif               # 3D 结构
        confidence_model_0.json   # pLDDT, PAE, iPTM
        affinity_model_0.json     # ΔG (Boltz-2 v2+)
```

**断点续跑：** 脚本自动检测 `.done` 标记，已完成批次自动跳过。

---

### Step 3：结果解析与可视化

```bash
python parse_boltz2_results.py \
    --results-dir ../../output/boltz2_results \
    --manifest ../../output/boltz2_blind_scan/job_manifest.csv \
    --output-dir ../../output/boltz2_analysis
```

**输出文件：**

| 文件 | 内容 |
|---|---|
| `affinity_scores.csv` | 每个 job 的 pLDDT / iPTM / ΔG |
| `delta_vs_wt.csv` | 每个突变相对 WT 的差异（Δ iPTM, Δ ΔG）|
| `subtype_prediction.csv` | 基于亲和力差异的自动分型建议 |
| `heatmap_delta_iptm.png` | 变异 × 配体热图（红=增强/蓝=减弱）|
| `heatmap_plddt.png` | 结构置信度热图 |
| `waterfall_GPIb_alpha.png` | 所有变异对 GPIbα 的亲和力变化瀑布图 |
| `waterfall_Collagen_I_THP.png` | 胶原结合变化瀑布图 |

---

## 分型推断逻辑

| 亲和力变化 | 配体 | 推断分型 |
|---|---|---|
| **Δ iPTM > +0.05** | GPIb_alpha | **2B** (GOF, AIM 破坏) |
| **Δ iPTM < -0.05** | GPIb_alpha | **2M** (A1 亲和力↓) |
| **Δ iPTM < -0.05** | Collagen_I_THP | **2M** (A3 胶原结合↓) |
| **Δ iPTM > +0.05** | ADAMTS13_Spacer | **2A** (A2 可及性↑) |
| **Δ iPTM < -0.05** | FVIII_LightChain | **2N** (FVIII 稳定性↓) |
| 所有配体均无显著变化 | — | **VUS / 非功能性** |

> **关键指标说明：**
> - **iPTM（inter-chain predicted TM-score）**：复合物界面置信度，越高 = 结合越稳定。是亲和力代理指标。
> - **Δ iPTM**：突变体 iPTM − WT iPTM；正值=结合增强，负值=结合减弱。
> - **Δ ΔG**：仅 Boltz-2 v2 提供；负值=结合更稳定，正值=结合削弱。

---

## 算力估算（4x H200）

| 变异总数 | 平均配体数/变异 | 总 job 数 | 预估时间 |
|---|---|---|---|
| 1305 (FIXED 表) | ~1.5 | ~2000 | ~8-12 小时 |
| 50（测试） | ~1.5 | ~75 | ~20 分钟 |

> 每张 H200 (80GB) 可并行处理约 8 个 job（序列 < 400 aa）。
> DDP 策略下 4 卡实际吞吐量约 28-32 job/小时。

---

## 依赖安装

```bash
# GPU 服务器（Linux）
conda create -n boltz python=3.11
conda activate boltz
pip install boltz>=1.0
pip install pandas openpyxl matplotlib seaborn numpy biopython

# 本地开发（Windows）
pip install pandas openpyxl matplotlib seaborn numpy
# boltz 本体不需在本地安装（仅生成输入 / 解析输出）
```

---

## 注意事项

1. **序列坐标系**：所有位置均使用 VWF 全长前体（2813 aa）的 1-indexed 坐标，与 UniProt P04275 一致。
2. **胶原三股螺旋**：`Collagen_I_THP` 自动生成 3 条相同链（B/C/D），勿手动改为 1 条。
3. **WT 基准**：建议在同一 batch 中包含 `VWF_WT_vs_<ligand>` 对照任务，供 `parse_boltz2_results.py` 计算 delta。
4. **亲和力 JSON**：`affinity_model_*.json` 仅 Boltz-2 v2 才输出；v1 用 iPTM 作为代理指标。

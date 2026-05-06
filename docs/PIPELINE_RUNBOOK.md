# VWF Diagnostic Pipeline — Operator Runbook

> **目标读者**: 你（在 GPU 实例上执行 pipeline 的人）
> **最后更新**: 2026-04-28
> **配套文档**: `DIAGNOSTIC_AGENT_V3_DESIGN.md` (架构设计)

---

## 0. TL;DR — 最短可执行路径

```bash
# === CPU 实例上准备 (一次性) ===
git clone https://github.com/LukaDD7/VWF-ETHos.git /workspace/vwf
cd /workspace/vwf

# === GPU 实例上跑 Boltz-2 ===
conda activate boltz   # 已预装的 conda env
cd Proteo-Structure-Pipeline/boltz2_pipeline
./run_boltz2.sh --input-dir priority1_2B2M    # 跑 P1: 479 jobs ~15-18h
./run_boltz2.sh --input-dir priority2_2A2N    # 跑 P2: 344 jobs ~11-13h

# === CPU/GPU 实例上跑 FoldX (M3, 待实现) ===
python scripts/pipeline/foldx_buildmodel_runner.py \
    --boltz-dir output/boltz2_results \
    --output output/foldx_bind.csv

# === 跑 v3 Agent (M5, 待实现) ===
python scripts/agentic_vwf_classifier_v3.py \
    --variants data/processed/master_type1_type2.csv \
    --output output/diagnosis_v3.csv
```

**完整周期**: Boltz-2 ~25-30h GPU + FoldX ~15-25h CPU + Agent ~1h ≈ **2-3 天可完成**。

---

## 1. 环境准备

### 1.1 仓库同步

```bash
# CPU 实例第一次
git clone https://github.com/LukaDD7/VWF-ETHos.git /workspace/vwf
cd /workspace/vwf

# 后续更新（每次开实例先做）
cd /workspace/vwf
git pull origin master
```

### 1.2 Conda 环境

GPU 服务器上**已预装 `boltz` env**（确认存在）：

```bash
conda env list | grep boltz
# 应该看到 boltz                  /path/to/envs/boltz
```

如果没装，按 `BOLTZ2_SPEC.md` 第 219 行：
```bash
conda create -n boltz python=3.11
conda activate boltz
pip install boltz>=1.0
pip install pandas openpyxl matplotlib seaborn numpy biopython pyyaml
```

### 1.3 FoldX 5

服务器上已躺着 (`FoldX5_LINUX/foldx_*` 二进制)。验证：

```bash
which foldx_20251231 || ls /opt/FoldX5/foldx*
foldx_20251231 --version    # 应该输出版本号
```

设置环境变量（写到 `~/.bashrc`）：
```bash
export FOLDX_BIN=/opt/FoldX5/foldx_20251231     # 改成你实际路径
export FOLDX_ROTABASE=/opt/FoldX5/rotabase.txt   # FoldX 必需依赖
```

---

## 2. 输入数据状态

### 2.1 已就绪

```
data/raw_tables/
├── VWF变异分析-用于提取Type-2.xlsx        ← 主表（337 变异）
├── 2B型突变.xlsx                          ← 2B 专用表
└── VWF_Type2_AF3_Reference_Table.xlsx    ← AF3 参考

data/processed/
├── master_type1_type2.csv                ← 已聚合 (337 行)
└── VWF_Alpha_Matrix.parquet              ← AlphaGenome 已跑完

structures/wt/
└── VWF_P04275_WT.fasta                   ← UniProt WT (2813 aa)

output/boltz2_blind_scan/
├── job_manifest.csv                      ← 824 个 job 索引
├── priority1_2B2M/  (16 batches, 479 jobs) ← 已生成
└── priority2_2A2N/  (12 batches, 344 jobs) ← 已生成
```

### 2.2 待生成

```
output/boltz2_results/<priority>/         ← Boltz-2 输出（M2）
output/foldx_bind_p1.csv                  ← FoldX AnalyseComplex ΔΔG_bind P1 (M3)
output/foldx_bind_p2.csv                  ← FoldX AnalyseComplex ΔΔG_bind P2 (M3)
output/foldx_fold.csv                     ← FoldX Stability ΔΔG_fold (M3)
output/diagnosis_v3.csv                   ← v3 agent 最终诊断 (M5)
output/evidence_chain.json                ← 推理可追溯链 (M5)
```

---

## 3. 跑 Boltz-2 (M2)

### 3.1 启动前检查

```bash
cd /workspace/vwf/Proteo-Structure-Pipeline/boltz2_pipeline

# 确认 GPU
nvidia-smi    # 应该看到 4×H200

# 确认 boltz 可用
which boltz && boltz --help | head -5

# 确认 input 目录存在
ls ../../output/boltz2_blind_scan/priority1_2B2M/ | head -5
```

### 3.2 跑 Priority 1 (2B/2M, 479 jobs)

```bash
chmod +x run_boltz2.sh
./run_boltz2.sh --input-dir priority1_2B2M

# 预期输出（开头）：
# ============================================================
# VWF Boltz-2 Production Run (Job-Level Incremental)
# Started: <date>
# ============================================================
# [HH:MM:SS] Expanding batches → single-job YAMLs in .../_tmp_single_jobs ...
#   Expanded 479 new  |  Total single-job YAML files: 479
# [HH:MM:SS] Total: 479  |  Done: 0  |  Remaining: 479
# [HH:MM:SS] [1/479] VWF_R1306W_vs_GPIb_alpha
#   [OK] → ../../output/boltz2_results/priority1_2B2M/VWF_R1306W_vs_GPIb_alpha
# ...
```

**预计 15-18 小时**。挂 `nohup` 或 `tmux`：
```bash
tmux new -s boltz-p1
./run_boltz2.sh --input-dir priority1_2B2M 2>&1 | tee run_p1.log
# Ctrl+B then D 离开
```

### 3.3 实时监控（另一个 shell）

```bash
cd /workspace/vwf/Proteo-Structure-Pipeline/boltz2_pipeline
python monitor_progress.py --results-dir ../../output/boltz2_results/priority1_2B2M --watch
```

或直接看 progress.json：
```bash
cat ../../output/boltz2_results/priority1_2B2M/progress.json
```

### 3.4 断点续跑

GPU 被回收/崩溃后重连：
```bash
cd /workspace/vwf/Proteo-Structure-Pipeline/boltz2_pipeline
./run_boltz2.sh --input-dir priority1_2B2M    # 自动跳过 .done 的 job
```

### 3.5 跑 Priority 2 (2A/2N, 344 jobs)

P1 完成后：
```bash
./run_boltz2.sh --input-dir priority2_2A2N
# ~11-13 小时
```

### 3.6 失败 job 处理

```bash
# 找出所有失败的 job (没有 .done)
cd ../../output/boltz2_results/priority1_2B2M
for d in */; do
    [ ! -f "$d/.done" ] && [ -d "$d" ] && echo "FAIL: $d"
done

# 重试单个 job
rm -rf VWF_R1306W_vs_GPIb_alpha
cd /workspace/vwf/Proteo-Structure-Pipeline/boltz2_pipeline
./run_boltz2.sh --input-dir priority1_2B2M
```

---

## 4. Boltz-2 输出验证

跑完 P1 后立即抽样验证：

```bash
cd /workspace/vwf/Proteo-Structure-Pipeline/boltz2_pipeline

# 4.1 计数 .done
find ../../output/boltz2_results/priority1_2B2M -name '.done' | wc -l
# 应该 ≈ 479

# 4.2 检查单个 job 的输出完整性
JOB=../../output/boltz2_results/priority1_2B2M/VWF_R1306W_vs_GPIb_alpha
ls $JOB/predictions/
# 期望:
#   model_0.cif (或 *_model_0.cif)
#   confidence_model_0.json
#   pae_model_0.npz, plddt_model_0.npz
# 注意: affinity_model_0.json **应该不存在**（已修 bug）

# 4.3 检查 confidence 的关键指标
python -c "
import json
d = json.load(open('$JOB/predictions/confidence_model_0.json'))
print('iPTM:', d.get('iptm'))
print('pTM:', d.get('ptm'))
print('mean pLDDT:', sum(d.get('plddt',[0]))/max(len(d.get('plddt',[1])),1))
"
# WT 通常 iPTM ≈ 0.7-0.9; mutant iPTM 与 WT 偏差 = 信号

# 4.4 确认 WT baseline job 全部完成（parse 脚本需要 WT 对照算 delta_iPTM）
echo "=== WT baseline jobs (priority1) ==="
ls ../../output/boltz2_results/priority1_2B2M/ | grep "^VWF_WT"
# 期望看到: VWF_WT_vs_GPIb_alpha, VWF_WT_vs_Heparan_Sulfate_mimic, VWF_WT_vs_Collagen_I_THP
# 如果缺少 WT baseline → 手动补跑:
#   boltz predict _tmp_single_jobs/VWF_WT_vs_GPIb_alpha.yaml --out_dir ...

# 4.5 确认 affinity_*.json 不存在（bug 已修验证）
echo "=== Checking for unexpected affinity files ==="
find ../../output/boltz2_results/priority1_2B2M \
    -name "affinity_*.json" | head -5
# 期望：无输出。有输出说明用了旧版 run_boltz2.sh → git pull && 补跑相关 job
```

---

## 5. 跑 FoldX (M3, 待实现) — 在高性能 CPU 实例上

> **实例选择**：开独立高性能 CPU 实例（32-64 核），**不要用 GPU 实例**。
> FoldX 是纯 CPU 计算；GPU 实例的显卡完全闲置但仍烧券。
> FoldX 和 Boltz-2 可以并行运行：Boltz-2 落盘的结构通过共享盘立即可读。

### 5.1 完整流程

```
Boltz-2 WT   model_0.cif ─→ CIF→PDB ─→ RepairPDB ─→ AnalyseComplex → G_bind(WT)
                                                                              │
                                                             ΔΔG_bind = mut - WT
                                                                              │
Boltz-2 mut  model_0.cif ─→ CIF→PDB ─→ RepairPDB ─→ AnalyseComplex → G_bind(mut)
```

两个结构都来自 Boltz-2，FoldX **只算 binding energy，不生成突变**（BuildModel 不需要）。

### 5.2 待写脚本: `scripts/pipeline/foldx_analysecomplex_runner.py`

**接口**:
```bash
# 在高性能 CPU 实例上跑
python scripts/pipeline/foldx_analysecomplex_runner.py \
    --boltz-dir /shared/output/boltz2_results/priority1_2B2M \
    --foldx-bin $FOLDX_BIN \
    --rotabase $FOLDX_ROTABASE \
    --output /shared/output/foldx_bind_p1.csv \
    --workers 32     # 按 CPU 核数设定
```

**输出 CSV schema**:
```
variant_id, ligand, G_bind_wt, G_bind_mut, ddG_bind_kcal_mol, iptm_wt, iptm_mut, status
VWF_R1306W, GPIb_alpha, -45.6, -47.94, -2.34, 0.81, 0.86, OK
VWF_R1306Q, GPIb_alpha, -45.6, -43.73, +1.87, 0.81, 0.74, OK
VWF_F1369I, GPIb_alpha, -45.6, -43.10, +2.50, 0.81, 0.71, OK
```

**M3 完成标准**：`output/foldx_bind_p1.csv` 存在，行数 ≈ 479，`status=OK` 行占比 > 85%。

### 5.3 FoldX 命令模板（核心逻辑）

```bash
# Step 1: CIF → PDB（每个 job 的 model_0.cif 各做一次）
python -c "
from Bio.PDB import MMCIFParser, PDBIO
import sys
p = MMCIFParser(QUIET=True)
io = PDBIO()
s = p.get_structure('x', sys.argv[1])
io.set_structure(s)
io.save(sys.argv[2])
" path/to/model_0.cif path/to/model_0.pdb

# Step 2: Repair（每个 structure 一次，耗时 5-15 min）
$FOLDX_BIN --command=RepairPDB \
    --pdb=model_0.pdb \
    --rotabaseLocation=$FOLDX_ROTABASE \
    --output-dir=./repaired/ \
    --quiet

# Step 3: AnalyseComplex（< 1 min，WT 和 mut 各跑一次）
$FOLDX_BIN --command=AnalyseComplex \
    --pdb=model_0_Repair.pdb \
    --analyseComplexChains=A,B \   # A=VWF domain, B=ligand（collagen: B+C+D）
    --rotabaseLocation=$FOLDX_ROTABASE \
    --output-dir=./analysis/ \
    --quiet
# 输出: Interaction_model_0_Repair_AC.fxout → 第一行含 binding energy (kcal/mol)

# Step 4: 提取 binding energy
python -c "
import re, sys
for line in open(sys.argv[1]):
    m = re.search(r'Group1.*?(-?\d+\.\d+)', line)
    if m: print(m.group(1)); break
" analysis/Interaction_model_0_Repair_AC.fxout
```

### 5.4 注意事项

| 问题 | 原因 | 解决 |
|---|---|---|
| FoldX Segfault | `$FOLDX_ROTABASE` 路径错 | 检查 rotabase.txt 是否在 `$FOLDX_BIN` 同目录 |
| PDB atom name 报错 | Boltz-2 CIF 非标准 atom name | 用 `pdb-tools`：`pdb_fixinsert < in.pdb > out.pdb` |
| ΔΔG 超 ±15 kcal/mol | iPTM < 0.5 的低置信度结构 | 先过滤 iPTM，低质量结构不进 FoldX |
| Collagen 三聚体链 ID | B/C/D 三链要整体作为一侧 | `--analyseComplexChains=A,B+C+D` |
| 并行 FoldX 写文件冲突 | 默认都写当前目录 | 每个 job 用独立 tmpdir |

---

## 6. v3 Agent 推理 (M5, 待实现)

### 6.1 待写脚本: `scripts/agentic_vwf_classifier_v3.py`

按 `DIAGNOSTIC_AGENT_V3_DESIGN.md` §4-§5 实现。

接口:
```bash
python scripts/agentic_vwf_classifier_v3.py \
    --variants data/processed/master_type1_type2.csv \
    --boltz-features output/boltz2_features.parquet \
    --foldx-bind output/foldx_bind.csv \
    --foldx-fold output/foldx_fold.csv \
    --alphagenome-features data/processed/VWF_Alpha_Matrix.parquet \
    --af3-features data/processed/af3_features.parquet \
    --llm-endpoint http://localhost:30000 \
    --output output/diagnosis_v3.csv \
    --evidence-out output/evidence_chain_v3.jsonl
```

### 6.2 LLM endpoint

如果你已经在另一台机器跑 SGLang + DeepSeek-V4-Flash（OpenAI-compatible），用那个：
```bash
--llm-endpoint http://gpu-server-ip:30000
```

否则用 Ollama 本地 + Qwen3-7B（够用）：
```bash
ollama pull qwen3:7b
--llm-endpoint http://localhost:11434
```

---

## 7. 验证与对比 (M6)

```bash
# 跑混淆矩阵 v3 vs v2
python scripts/eval/compare_v2_v3.py \
    --v2 output/diagnosis_v2.csv \
    --v3 output/diagnosis_v3.csv \
    --truth data/processed/master_type1_type2.csv \
    --output figures/confusion_matrix_v2_v3.png

# 误差分析: 哪些变异 v3 改对了/改错了？
python scripts/eval/error_analysis.py \
    --v2 output/diagnosis_v2.csv \
    --v3 output/diagnosis_v3.csv \
    --output output/error_analysis_v3.csv
```

**目标指标**:
- 整体准确率: 65% → 80%+
- 2B 准确率: 17% → 60%+
- 2N 准确率: 65% → 80%+

---

## 8. 故障速查

| 症状 | 原因 | 解决 |
|---|---|---|
| `boltz: command not found` | conda env 没激活 | `conda activate boltz` |
| `nvidia-smi` 看不到 GPU | 没分到 GPU 实例 | 重申请实例 |
| 单个 job 卡死 > 30 分钟 | 序列太长或显存爆 | `Ctrl+C` 后 `rm -rf 该 job 目录` 重跑 |
| `OutOfMemoryError` | 超大复合物（VWF + collagen 4 链） | 减 `DIFFUSION_SAMPLES=3` 或单卡 `N_GPUS=1` |
| FoldX `Segmentation fault` | rotabase.txt 路径错 | 检查 `$FOLDX_ROTABASE` |
| FoldX ΔΔG 超出 -10/+10 范围 | 结构不可信（iPTM<0.6） | 在 FoldX 前过滤掉低 iPTM 的 job |
| LLM endpoint 不响应 | SGLang 挂了 | 重启 SGLang server |

---

## 9. Checkpoint Summary

执行进度自查（每完成一个里程碑勾掉）：

- [ ] M1: `run_boltz2.sh` Boltz-2 affinity bug 已修
- [ ] M2: Boltz-2 priority1 (479 jobs) 全完成
- [ ] M2: Boltz-2 priority2 (344 jobs) 全完成
- [ ] M3: `foldx_analysecomplex_runner.py` 实现（高性能 CPU 实例）
- [ ] M3: ΔΔG_bind 跑完，`output/foldx_bind_p1.csv` 行数 ≈ 479，OK 率 > 85%
- [ ] M3: FoldX Stability runner 实现 + 全 ΔΔG_fold 跑完，`output/foldx_fold.csv` 就绪
- [ ] M4: BindingAffinityExpert 实现
- [ ] M4: FoldingStabilityExpert 实现
- [ ] M4: LiteratureExpert 实现
- [ ] M5: ClinicalGeneticist v3 实现
- [ ] M6: 验证混淆矩阵
- [ ] M7: 文档/Demo

---

## 10. 联系与求助

- 文档不一致：直接编辑本文件并 PR
- 脚本 bug：在 git issues 提
- 算力不足：联系学院算力运维

---

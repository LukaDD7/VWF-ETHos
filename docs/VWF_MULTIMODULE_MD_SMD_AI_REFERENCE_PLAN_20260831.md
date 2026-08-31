# VWF 多模块 MD/SMD 与 AI 参考分布方案

日期：2026-08-31  
状态：设计已落地，等待 GPU 释放后执行  
适用范围：VWF 完整诊断系统，而不是单一 A1/AIM 结构

---

## 1. 结论先行

1. **第三轮不是 MD/SMD，而是多模块 MD/SMD + AI 参考分布。**
   - 当前 GPU 2/3 正在跑的是 **Boltz-2 keepalive**，不是 GROMACS MD/SMD。
   - 真正的第三轮应当是把 A1/AIM 之外的其他 VWF 模块补齐 MD/SMD 度量，并给每个模块建立独立的 AI 参考分布。
2. **目前只有 A1 轴有完整 MD 结果层。**
   - A1/AIM：7 个 7A6O 匹配系统，各 50 ns。
   - A1-GPIbα：37 个系统的 equilibrium interface-retention 结果已并入
     `outputs/computational_panel_20260829/md/a1_gpiba_interface/`。
   - A1/AIM salt-bridge 和 slow025 SMD 标定也已并入
     `outputs/computational_panel_20260829/md/`。
   - 未完成：D′D3/FVIII、A2、A3、C1/C2/C4、D1/D2、D4、C1-C6、CK 等非 A1 模块均没有 MD/SMD 结果。
3. **初始结构必须区分实验晶体与 Boltz 预测。**
   - 已有本地实验结构：**7A6O AIM-A1**、**1SQ0/1M10 A1-GPIbα**。
   - 其他模块也有明确优先的实验结构候选：A2 用 `3GXB`，A3-collagen 用 `4DMU`，D′D3 用 `6N29`，CTCK/CK 用 `4NT5`。
   - `3GXB`、`4DMU`、`6N29`、`4NT5` 已下载到本仓库 `structures/`，后续先清理再进入 MD。
   - Boltz-2 只作为**没有实验结构或实验结构尚未下载时的 fallback**，不能覆盖已有实验坐标。
4. **Boltz-2 和 AlphaGenome 本轮不需要额外补充。**
   - Boltz-2 已覆盖 15 个功能轴，`evidence_matrix.csv` 已有 597 个变异。
   - AlphaGenome 当前 16 个请求已完成；如果加入新样本，再按既有规则扩展。
   - 本轮真正缺的是 MD/SMD 动态度量和分模块 AI 参考分布，不是再跑一遍静态预测。
5. **SMD 目前不是首选。**
   - A1/AIM 的旧 SMD 已被判定为 no-go：加载速率过快、反应坐标不等价、2B/2M 方向反了。
   - 其他模块只有在“平衡态 MD 稳定 + 对照齐全 + 反应坐标明确”时才考虑 SMD，否则不进入分类器。

---

## 2. 已有产物与新增产物

### 2.1 已有 MD 结果

- `outputs/computational_panel_20260829/md/README.md`
- `outputs/computational_panel_20260829/md/md_feature_matrix.csv`
- `outputs/computational_panel_20260829/md/qc_summary.csv`
- `outputs/computational_panel_20260829/md/aim_a1_contacts_summary.csv`
- `outputs/computational_panel_20260829/md/aim_a1_min_distance_summary.csv`
- `outputs/computational_panel_20260829/md/rmsd_timeseries.csv`
- `outputs/computational_panel_20260829/md/artifact_manifest.csv`
- `outputs/computational_panel_20260829/md/md_result_layers.csv`
- `outputs/computational_panel_20260829/md/a1_aim_masking/`
- `outputs/computational_panel_20260829/md/a1_aim_saltbridge/`
- `outputs/computational_panel_20260829/md/a1_gpiba_interface/`
- `outputs/computational_panel_20260829/md/smd/`

覆盖系统：

- `PANEL_WT_MATCHED`
- `P1413L`
- `R1308C`
- `S1310F`
- `V1316M`
- `R1341W`
- `A1461D`

每个系统 50 ns，匹配 WT 基线为 `PANEL_WT_MATCHED`。

### 2.2 新增的机器可读清单

本次补了两个文件：

1. `outputs/computational_panel_20260829/md/module_structure_inventory.csv`
   - 每个 assay 一个优先初始结构。
   - 有实验 PDB 的模块优先记录实验 PDB；本地已有则直接给出路径。
   - Boltz-2 只作为无实验结构模块的 fallback，并记录最优模型、置信度、模型数、平均分和 CIF 路径。
2. `outputs/computational_panel_20260829/md/module_md_readout_plan.csv`
   - 每个模块的平衡态 MD 读数、可选 SMD 读数、SMD 门槛、AI 参考分布特征和优先级。

同时把已有但未入当前提交的 A1 轴结果层并入：

- `a1_aim_masking/`：15 行 AIM-A1 masking/interface 特征和 residue metadata。
- `a1_aim_saltbridge/`：40 行 AIM-A1 salt-bridge occupancy 特征。
- `a1_gpiba_interface/`：37 行 A1-GPIbα interface-retention 分类特征、
  sanitized summary 和 1887 帧时间序列。
- `smd/`：3 行 slow025 汇总和 15 行 per-replicate 力/功曲线，状态为
  `complete_no_go`。

生成结构清单的命令：

```bash
python3 scripts/pipeline/build_vwf_md_structure_inventory.py
```

---

## 3. 多模块结构库存结论

15 个功能轴都有候选初始结构，但必须先按来源分层：

| 模块 | 优先实验结构 | 本地状态 | 结论 |
|---|---|---|---|
| A1/AIM | `7A6O` | 已有 | 实验晶体，优先级最高 |
| A1-GPIbα | `1SQ0`，备用 `1M10` | 已有 | 实验晶体，优先于 Boltz |
| A2 folded | `3GXB` | 已有 | 优先于 Boltz |
| A3-collagen | `4DMU` | 已有 | 优先于 Boltz |
| D′D3 | `6N29` | 已有 | 优先于 Boltz |
| CK/CTCK | `4NT5` | 已有 | 优先于 Boltz |

没有实验结构的模块才使用本地 Boltz-2 WT 模型。当前 Boltz fallback 质量分层如下：

| 置信度分层 | 模块 | 最好模型分数 | 结论 |
|---|---|---:|---|
| 高 | A1-heparan sulfate | 0.856 | 适合做平衡态 MD |
| 中 | C1-collagen | 0.609 | 可做，但生物学优先级低于 A3 |
| 中 | C2-collagen | 0.585 | 可做，但生物学优先级低于 A3 |
| 中 | C4-integrin | 0.687 | 可做，用于 RGD/整合素轴 |
| 低 | A2-ADAMTS13 | 0.248 | 需要更可靠复合物模型 |
| 低 | VWF73-ADAMTS13 | 0.106 | 只能作为探索性模型 |
| 低 | D1-D2 propeptide | 0.246 | 只能作为探索性模型 |
| 低 | D4 assembly | 0.283 | 只能作为探索性模型 |
| 低 | C1-C6 assembly | 0.358 | 只能作为探索性模型 |

注意：

- iPTM 和 pTM 不能跨 assay 直接横向比较；上表只用于同一 assay 内挑模型和判断是否需要额外验证。
- “有模型”不等于“能直接跑 MD”。所有 Boltz CIF 都需要：
  1. CIF → PDB；
  2. `pdb2gmx`；
  3. 分级弛豫；
  4. 稳定性 QC；
  5. 20 ns pilot；
  6. 通过后再升到 50 ns。
- 低置信度模块的 MD 结果只能标为 exploratory，不能单独作为临床分型证据。

---

## 4. 模块化 MD/SMD 读数方案

详细机器可读版本见：

```text
outputs/computational_panel_20260829/md/module_md_readout_plan.csv
```

核心原则：

- **平衡态 MD 是默认层。**
- **SMD 是条件层，不是默认层。**
- **每个模块的读数必须模块化，不做全局统一阈值。**

### 4.1 A1/AIM

- 结构：实验 `7A6O`
- 已有：50 ns × 7 系统
- 平衡态读数：
  - AIM-A1 接触占有率；
  - 关键盐桥占有率，尤其是 D1269-R1306；
  - A1 结合面暴露；
  - backbone RMSD；
  - 局部柔性。
- SMD：
  - 只允许在重新设计反应坐标和对照后做；
  - 不复用旧的 anchor-to-anchor 解折叠力轴。
- AI 参考分布：
  - `aim_contact_occupancy_z`
  - `aim_salt_bridge_retention_z`
  - `a1_exposure_z`
  - `backbone_rmsd_z`

### 4.2 A1-GPIbα

- 结构：实验 `1SQ0`，备用 `1M10`
- 平衡态读数：
  - 界面接触占有率；
  - 界面 RMSD；
  - GPIbα 接触寿命；
  - A1 表面暴露。
- SMD：
  - 条件允许时可做 A1-GPIbα 解离；
  - 必须使用匹配的 2B/2M 对照和合理反应坐标。
- AI 参考分布：
  - `interface_contact_occupancy_z`
  - `interface_rmsd_z`
  - `contact_lifetime_z`
  - `surface_exposure_z`

### 4.3 A2

A2 有两个层次：

1. **A2 folded stability**
   - 结构：实验 `3GXB`，本地已有
   - 平衡态读数：
     - backbone RMSD；
     - 二级结构保留率；
     - 局部柔性；
     - Tyr1605-Met1606 切割位点暴露。
   - 可选 SMD：
     - A2 机械解折叠力，仅在直接研究力依赖解折叠时做。
   - AI 参考分布：
     - `backbone_rmsd_z`
     - `secondary_structure_retention_z`
     - `local_flexibility_z`
     - `cleavage_site_exposure_z`

2. **A2-ADAMTS13 / VWF73**
   - 当前 Boltz 置信度低：0.248 / 0.106
   - 先做模型验证或换更可靠复合物；
   - 若可用，读数为：
     - 界面接触占有率；
     - 界面 RMSD；
     - 切割位点几何；
     - 底物接触寿命。
   - AI 参考分布：
     - `interface_contact_occupancy_z`
     - `interface_rmsd_z`
     - `cleavage_site_geometry_z`

### 4.4 A3 / C1 / C2 collagen

- A3 是主轴，C1/C2 是辅助轴。
- A3 优先使用实验 `4DMU` collagen-bound 结构；本地已有。
- 平衡态读数：
  - 胶原界面接触占有率；
  - 界面 RMSD；
  - 胶原接触寿命；
  - 结合面暴露。
- SMD：
  - A3-collagen 可选；
  - C1/C2 不建议优先做。
- AI 参考分布：
  - `interface_contact_occupancy_z`
  - `interface_rmsd_z`
  - `contact_lifetime_z`

### 4.5 C4 integrin

- 平衡态读数：
  - 界面接触占有率；
  - 界面 RMSD；
  - RGD loop 几何；
  - 接触寿命。
- SMD：
  - 不建议优先做。
- AI 参考分布：
  - `interface_contact_occupancy_z`
  - `interface_rmsd_z`
  - `rgd_loop_geometry_z`

### 4.6 D′D3-FVIII

- 优先使用实验 `6N29` D′D3 结构；本地已有。
- 下载并清理后，平衡态读数：
  - FVIII 界面接触占有率；
  - 界面 RMSD；
  - 界面 SASA；
  - 关键残基对距离稳定性。
- SMD：
  - 不是必需；FVIII 结合不是典型力依赖开关。
- AI 参考分布：
  - `interface_contact_occupancy_z`
  - `interface_rmsd_z`
  - `interface_sasa_z`
  - `residue_pair_distance_z`

### 4.7 D1/D2、D4、C1-C6、CK

这些模块主要服务定量 VWD、分泌、多聚化和二聚化轴。

- 平衡态读数：
  - backbone RMSD；
  - 二级结构保留率；
  - 局部柔性；
  - 结构域 packing；
  - CK 二聚化表面保留（需二聚体模型）。
- CK/CTCK 优先使用实验 `4NT5`；本地已有。
- SMD：
  - 不建议默认做。
- AI 参考分布：
  - `backbone_rmsd_z`
  - `secondary_structure_retention_z`
  - `local_flexibility_z`
  - `domain_packing_z`
  - `dimerization_surface_retention_z`

---

## 5. AI 参考分布规则

### 5.1 不做全局统一阈值

每个模块、每个 assay、每个 feature 单独建分布：

```text
module + assay + feature -> reference distribution
```

例如：

```text
A1_AIM + a1_aim_autoinhibition_context + aim_salt_bridge_retention
A3    + a3_collagen_binding           + interface_contact_occupancy
A2    + a2_folded_stability           + cleavage_site_exposure
```

### 5.2 参考集组成

每个模块至少要有：

1. WT matched baseline；
2. 良性/人群对照；
3. 已知致病或亚型代表变异；
4. 若可能，加入亚型内部的不同机制对照。

最低建议：

| 组别 | 最低数量 | 说明 |
|---|---:|---|
| WT | 1 个 matched baseline | 必须同一初始结构和同一流程 |
| Benign/control | ≥3 | 同域或邻近区域优先 |
| Pathogenic/representative | ≥3 | 尽量覆盖目标亚型 |
| 边界 case | 可选 | 用于校准不确定区间 |

如果样本数不足，结果只能标为 exploratory，不能进入强分类规则。

### 5.3 每个 feature 的输出字段

```text
module
assay
feature
n_wt
n_benign
n_pathogenic
wt_mean
wt_sd
benign_mean
benign_sd
pathogenic_mean
pathogenic_sd
effect_size
z_score
empirical_percentile
direction
confidence_interval
missingness_reason
```

建议计算：

- mean ± SD；
- median 和 IQR；
- Cohen’s d 或 rank-biserial effect size；
- 相对 WT 的 z-score；
- 经验百分位；
- bootstrap 95% CI；
- 方向性标注，例如 `loss`、`gain`、`destabilization`、`exposure_increase`。

### 5.4 Agent 接入字段

Agent v3 不应直接吃原始 MD 轨迹，而是吃结构化特征：

```json
{
  "module": "A3",
  "assay": "a3_collagen_binding",
  "feature": "interface_contact_occupancy",
  "value": 0.42,
  "wt_mean": 0.61,
  "wt_sd": 0.08,
  "z_score": -2.37,
  "percentile": 0.03,
  "direction": "loss",
  "reference_set": "matched WT + benign + pathogenic A3 controls",
  "confidence": "moderate",
  "caveat": "Boltz-derived starting structure; interpret as mechanism evidence"
}
```

建议新增到 `ExpertScores` 的字段：

```text
module_md_stability_z
module_interface_retention_z
module_unfolding_z
module_exposure_z
module_smd_force_z
```

其中 `module_smd_force_z` 只有在 SMD 通过 go/no-go 后才允许启用。

---

## 6. GPU 执行方案

### 6.1 当前状态

- GPU 2/3 正在跑 Boltz-2 keepalive。
- 这不是 MD/SMD。
- 不建议直接杀掉正在跑的 keepalive；先等它结束或由用户明确决定停止。

### 6.2 GPU 释放后的推荐队列

每张卡一个 worker：

- GPU 2：高优先级队列
  1. A2 folded stability
  2. A3 collagen
  3. A1-GPIbα
- GPU 3：校准与辅助队列
  1. A1-heparan sulfate
  2. C4 integrin
  3. A1-AIM 扩展对照
  4. 低置信度模块 pilot

每个模块按以下顺序执行：

```text
CIF -> PDB -> pdb2gmx -> staged relaxation -> 20 ns pilot -> QC -> 50 ns production -> feature extraction
```

### 6.3 不建议一开始就全量 50 ns

原因：

1. 低置信度模型可能弛豫失败；
2. 有些模块构象可能不稳定；
3. 全量 50 ns 会浪费 GPU 时间；
4. 20 ns pilot 能先筛掉明显不合格的模型。

推荐门槛：

- 20 ns 内无爆炸；
- RMSD 不持续发散；
- 界面或功能位点没有立即崩塌；
- 温度、压力、密度稳定；
- 再升到 50 ns。

---

## 7. SMD go/no-go 规则

SMD 只有满足以下全部条件才进入下一轮：

1. **平衡态 MD 已经稳定。**
2. **反应坐标与生物学问题一致。**
   - 例如 AIM 释放、A1-GPIbα 解离、A2 解折叠；
   - 不再用与机制无关的端到端整体拉伸。
3. **有匹配对照。**
   - WT；
   - 已知 2B/2M 或目标亚型代表；
   - 良性对照。
4. **加载速率和输出物理解释合理。**
   - 不再把 1000 pN 级别的人工高力当作生物学断裂力。
5. **每个条件至少 5 个 replica。**
6. **先做 AUC 和方向性检查。**
   - 方向不对不能翻符号硬上。
7. **通过后再进入 Agent。**

在这些条件满足前，SMD 只保留为可选验证，不作为默认诊断特征。

---

## 8. 与 Boltz / AlphaGenome 的关系

### 8.1 Boltz-2

- 已有 15 个功能轴和 597 个变异的静态证据矩阵。
- 本轮不需要为了 MD 再重跑 Boltz。
- Boltz 的作用：
  1. 提供多模块候选初始结构；
  2. 提供 assay 级静态置信度；
  3. 与 MD 动态特征互补。
- 新样本进入时，继续使用现有 `generate_vwd_functional_boltz2_yamls.py` 流程。

### 8.2 AlphaGenome

- 当前 16 个请求已完成。
- AlphaGenome 不提供蛋白结构，不能替代 MD 初始结构。
- 新样本规则：
  - 主 biosample 用 `CL:0002618`；
  - ATAC/contact maps 缺失时用 all-track fallback；
  - 保留 signed delta 和 missingness reason。
- 本轮不需要额外 AlphaGenome 计算，除非加入新样本。

---

## 9. 立即可做的事

1. **已完成：结构库存**
   ```bash
   python3 scripts/pipeline/build_vwf_md_structure_inventory.py
   ```
2. **已完成：模块化读数与 AI 参考分布计划**
   - `module_md_readout_plan.csv`
3. **下一步：把高置信度模块转成 MD-ready**
   - A2 folded；
   - A1-heparan sulfate；
   - A3 collagen；
   - A1-GPIbα。
4. **再下一步：20 ns pilot**
   - 每张 GPU 一个 worker；
   - 不并行跑超过 2 个 MD 任务；
   - 避免 CPU 超订。
5. **最后：建立模块级参考分布**
   - 按 module + assay + feature 建分布；
   - 输出 z-score、percentile、effect size 和不确定性；
   - 再接入 Agent v3。

---

## 10. 一句话总结

**A1/AIM 已经有完整 MD，但完整 VWF 诊断系统还缺其他模块的动态度量。现在不是继续堆静态预测，而是用本地 Boltz 模型补多模块平衡态 MD，按模块建立 AI 参考分布；SMD 只在通过明确 go/no-go 后作为补充验证。**

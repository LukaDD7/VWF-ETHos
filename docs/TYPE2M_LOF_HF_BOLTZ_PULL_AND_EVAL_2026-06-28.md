# Type2M LOF HF Boltz panel pull + Eval readout

2026-06-28。服务器端补充 Boltz-2 结果已 push，原始结果已同步到 Hugging Face dataset：
`lucachangretta/VWF/type2m_lof_panel/`。

## 本地同步结果

- Git 已快进到 `02913c5 docs: Type2M LOF Boltz-2 supplementary panel results + HF upload`。
- HF raw results 已下载到
  `output/hf_type2m_lof_panel/type2m_lof_panel/boltz_results/`：
  753 个文件，约 120 MB。
- HF analysis files 已下载到
  `output/hf_type2m_lof_panel/type2m_lof_panel/analysis/`。

## 新增实现

- `scripts/pipeline/integrate_type2m_lof_hf_panel.py`
  - 将 HF long-form Boltz summary 转成 `evidence_matrix.csv` 兼容的 wide rows。
  - 排除 4 个 WT 对照行；WT 只作为 `delta_vs_wt` 基线。
  - `delta_vs_wt` 使用本次 HF supplement 自带 WT。
  - `zscore_within_assay` 使用原始 evidence matrix 的 assay 分布，保持分类器阈值含义不漂。
- `scripts/pipeline/evaluate_vwf_classifier_v2.py`
  - 加入 `a3_collagen_zscore` 映射，供分类器 A3 rule 解释使用。
- `scripts/agentic_vwf_classifier.py`
  - RULE7(A3) 仍以 A3 domain → 2M 为主规则。
  - 新增 A3 collagen-binding z-score 作为置信度/解释层：
    - `z <= -1.0`：支持 collagen-binding LOF，提高 2M 置信度。
    - `z >= +1.0`：静态 Boltz 不支持 LOF，保留 A3→2M 但降低置信度。
    - 中间区间：弱/中性证据，不改变 subtype。

## 输出文件

- Combined evidence matrix：
  `output/hf_type2m_lof_panel/type2m_lof_panel/analysis/evidence_matrix_with_type2m_lof_hf.csv`
- Supplement-only wide rows：
  `output/hf_type2m_lof_panel/type2m_lof_panel/analysis/type2m_lof_hf_evidence_rows.csv`
- Eval before supplement：
  `output/eval_v2_base_current/`
- Eval after supplement：
  `output/eval_v2_with_type2m_lof_hf/`

## 数据量变化

补充前：

| dataset | clean supported joined to Boltz |
|---|---:|
| Type1 | 115 |
| 2A | 118 |
| 2B | 38 |
| 2M | 37 |
| 2N | 16 |
| total supported | 324 |
| missing clean supported Boltz | 16（全为 2M） |

补充后：

| dataset | clean supported joined to Boltz |
|---|---:|
| Type1 | 115 |
| 2A | 118 |
| 2B | 38 |
| 2M | 53 |
| 2N | 16 |
| total supported | 340 |
| missing clean supported Boltz | 0 |

Type2 Eval 样本从 209 增至 225；新增 16 个全是 clean 2M。

## Recall 变化

使用现有 MD 特征时：

| label | before supplement | after supplement |
|---|---:|---:|
| 2A | 70/118 = 59% | 70/118 = 59% |
| 2B | 29/38 = 76% | 29/38 = 76% |
| 2M | 9/37 = 24% | 21/53 = 40% |
| 2N | 12/16 = 75% | 12/16 = 75% |
| Type2 total | 120/209 = 57% | 132/225 = 59% |

Confusion after supplement:

| true | pred 2A | pred 2B | pred 2M | pred 2N | pred uncertain |
|---|---:|---:|---:|---:|---:|
| 2A | 70 | 2 | 0 | 28 | 18 |
| 2B | 1 | 29 | 2 | 0 | 6 |
| 2M | 0 | 6 | 21 | 0 | 26 |
| 2N | 4 | 0 | 0 | 12 | 0 |

## 新增 16 个 2M 的逐例读数

### A1 2M

| variant | result | key evidence |
|---|---|---|
| V1409F | 2M | fb z=-1.21, heparan z=-0.91, mean=-1.06 <= -0.75 |
| P1413R | uncertain | fb/heparan 仅轻度下降，未达 LOF 阈值 |
| I1416T | uncertain | fb/heparan 均偏高，不支持静态 LOF |
| I1425F | uncertain | fb/heparan 均偏高，不支持静态 LOF |
| R1426C | uncertain | fb/heparan 均偏高，不支持静态 LOF |

结论：A1-2M 静态 Boltz 只救回 V1409F；其余 4 个仍需要 A1-GPIb 复合物 MD、AIM-A1 closed-state MD 或新的 LOF 动态特征。

### A3 2M

11 个 A3 新样本全部由 A3 domain 机制判为 2M；A3 collagen z-score 现在进入 confidence/reasoning。

| evidence group | variants |
|---|---|
| collagen static LOF 支持强 | L1733P (z=-1.51) |
| 弱/中性 | H1786D, V1730G, S1699F, V1822E, L1696R |
| 静态 Boltz 不支持 LOF，低置信保留 A3→2M | M1761K, W1745C, S1783A, Y1735D, S1731T |

这说明 A3 机制方向是对的，但 Boltz 静态 collagen complex 对所有 2M 并不单调；后续不能只用单一静态 iPTM 当 A3-LOF 判据。

## 标签泄漏核实

`evaluate_vwf_classifier_v2.py` 有 leakage smoke test：

1. 正常运行预测。
2. 将 `true_label` 和 `type2_subtype` 全部改成 `LEAK_TEST`。
3. 再跑一次预测。
4. 比较两次 `pred_with_md` 是否完全一致。

本次结果：`label_leak_test_pass = 1`。也就是说，当前分类过程不靠 subtype/true label 做判断；标签只在预测后用于 recall/confusion 统计。

## 当前可行性判断

系统继续做有潜力，但瓶颈已经比较清楚：

- 2B 方向：7A6O AIM-A1 saltbridge MD 对 non-hotspot 2B 已经证明有增益，2B recall 可到 29/38=76%，A1 2B recall 29/31=94%。
- 2M 方向：补齐 16 个 Type2M LOF Boltz 后，2M recall 从 24% 到 40%，说明补面板有实际收益。
- 最大缺口：A1-2M 仍大量 uncertain，静态 GPIb/heparan Boltz 对一部分 true 2M 不掉分。这里需要继续做动态 LOF 轴，而不是简单下调静态阈值。
- A3 方向：domain 机制强，但 collagen Boltz 有非单调样本，应作为 confidence/triage 特征，不宜单独硬判。

下一步优先级：

1. A1-2M：补 `1SQ0 A1-GPIb complex MD`，看 GPIb interface stability / contact retention / binding-face RMSF。
2. A1-2M：补 `7A6O AIM-A1 closed-state MD`，看是否是 closed-state 结合面塌陷或 hyper-stabilized low-adhesion dynamics。
3. A3-2M：挑 `L1733P` 做阳性机制例，挑 `M1761K/W1745C` 做反例/失败模式，避免静态 collagen iPTM 被过度解释。
4. 保持 leakage smoke test 作为每次 eval 的固定验收项。

# Partial Type2M LOF MD validation

2026-06-30。A40 推回 `output/type2m_lof_md_analysis_2026-06-29` 后，先把已有
7A6O AIM-A1 closed-state MD 接入 Eval v2 做快速验证。A1-GPIbα MD 当前只有 RMSD QC，
暂不作为分类特征。

## 输入

- MD analysis dir: `output/type2m_lof_md_analysis_2026-06-29`
- Validator: `scripts/pipeline/validate_type2m_lof_md_fast.py`
- Threshold sweep wrapper: `scripts/pipeline/sweep_type2m_lof_md_thresholds.py`
- Main output dir: `output/type2m_lof_md_fast_validation_2026-06-29`
- Sweep summary: `output/type2m_lof_md_threshold_sweep_2026-06-30/threshold_sweep_summary.csv`

## 特征定义

对每个已完成的 7A6O closed-state trajectory，统计 AIM-A1 nonlocal contact count：

```text
md_closed_aim_contact_loss =
    AIM_all_contacts_first0_5 - AIM_all_contacts_tail40_50

md_face_destab_score =
    md_closed_aim_contact_loss / closed_contact_loss_threshold
```

分类器当前阈值仍是：

```text
md_face_destab_score >= 1.0  ->  MD LOF / 2M soft evidence
```

因此 `closed_contact_loss_threshold=20` 等价于：

```text
0-5 ns 到 40-50 ns 的 AIM-A1 contacts 丢失 >= 20 个
  -> md_face_destab_score >= 1.0
  -> 触发 2M LOF 软证据
```

注意：这个特征是 LOF/2M 轴，不是 2B 阳性轴。它的含义是 closed-state 中 AIM-A1
局部接触在 equilibrium MD 里明显丢失，提示 A1 binding-face / self-inhibited contact
network 不稳定。当前分类器只把它当作 2M 软证据。

## 默认阈值结果

命令：

```bash
python scripts/pipeline/validate_type2m_lof_md_fast.py
```

默认 `closed_contact_loss_threshold=20`。

Eval v2 变化：

| label | baseline correct/total | with partial MD correct/total | delta |
|---|---:|---:|---:|
| 2A | 70/118 = 59.3% | 70/118 = 59.3% | 0 |
| 2B | 15/38 = 39.5% | 15/38 = 39.5% | 0 |
| 2M | 21/53 = 39.6% | 31/53 = 58.5% | +10 |
| 2N | 12/16 = 75.0% | 12/16 = 75.0% | 0 |
| ALL Type2 | 118/225 = 52.4% | 128/225 = 56.9% | +10 |

被救回的 10 个均为 2M，从 `uncertain -> 2M`：

```text
A1355D, A1377V, F1369I, G1324R, K1362T,
L1276P, L1276R, L1383R, S1387I, Y1363C
```

## 阈值扫描

命令：

```bash
python scripts/pipeline/sweep_type2m_lof_md_thresholds.py
```

结果：

| contact loss threshold | strong MD rows | 2M delta correct | 2M recall | 2B delta correct | ALL delta correct |
|---:|---:|---:|---:|---:|---:|
| 10 | 15 | +15 | 67.9% | 0 | +15 |
| 15 | 13 | +13 | 64.2% | 0 | +13 |
| 20 | 10 | +10 | 58.5% | 0 | +10 |
| 25 | 8 | +8 | 54.7% | 0 | +8 |
| 30 | 7 | +7 | 52.8% | 0 | +7 |
| 35 | 6 | +6 | 50.9% | 0 | +6 |
| 40 | 4 | +4 | 47.2% | 0 | +4 |

## 阈值怎么定

当前不能直接选择 recall 最高的 `10 contacts`。原因是这次已完成的 7A6O closed-state
feature rows 全部是 clean 2M：

```text
completed_md_rows = 22
completed_md_label_counts = 2M:22
```

所以现在验证到的是：

- 这个 MD LOF 轴确实能把一批 2M 从 `uncertain` 拉回 `2M`。
- 当前数据还不能估计 2B false-positive risk，因为 2B hard-negative closed-state MD 尚未完成或尚未进入该 feature table。

现阶段建议：

- 报告/代码继续使用 `closed_contact_loss_threshold=20` 作为 provisional conservative threshold。
- 这是一个保守点：能救回 10 个 2M，又不会过度扩张到轻微 contact drift。
- 等 `C1272R`、`V1316M` 等 2B hard-negative 以及更多 2M 完成后，再用 2B false-positive
  约束重扫阈值。如果 2B 不被误伤，可考虑放宽到 15 或 10。

最终阈值应按以下原则定：

1. **先定方向**：contact loss 高只能作为 2M/LOF 证据，不能作为 2B 证据。
2. **以 precision 约束 recall**：优先避免把真 2B 错判成 2M，尤其是 C1272R、V1316M 这类 hard negatives。
3. **用 held-out / newly completed MD 复核**：不能只在本批 22 个 2M 上选最优点。
4. **阈值写成 provisional**：当前 `20 contacts` 是工程阈值，不是最终生物物理常数。

## 结论

这批 partial MD 已经可以接入分类系统做验证，而且有正向效果：在不改变 2A/2B/2N
结果的前提下，2M recall 从 39.6% 提到 58.5%。但由于 completed closed-state rows
目前全是 2M，这个结果证明的是 2M rescue potential，尚不能证明 2B/2M boundary 已稳定。

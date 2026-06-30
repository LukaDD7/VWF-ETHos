# A1-GPIbα MD interface retention runbook

2026-06-30。用于 A40/H200 服务器拿到 A1-GPIbα production MD 后，抽取可进入分类器阈值扫描的
interface contact / retention 方向特征。

## 目的

之前 A1-GPIbα MD 只生成了 backbone RMSD QC。RMSD 只能说明整体构象/对齐波动，不能直接说明
GPIbα binding interface 是否保留或丢失。本脚本改为直接统计 VWF A1 与 GPIbα 两个 group
之间的 contacts 和 minimum distance。

## 脚本

```bash
python scripts/pipeline/analyze_a1_gpiba_completed_md.py
```

默认读取：

```text
output/type2m_lof_md_analysis_2026-06-29/a1_gpiba_completed_and_running_summary.csv
output/gromacs_md_a1_gpiba/<variant>/md_a1_gpiba/md_prod.*
```

默认输出：

```text
output/type2m_lof_md_analysis_2026-06-29/a1_gpiba_interface_qc/
  a1_gpiba_interface_summary.csv
  a1_gpiba_interface_timeseries.csv
  a1_gpiba_classifier_features.csv
  a1_gpiba_interface_skipped.csv
  logs/
```

## A40 推荐命令

在 repo 根目录运行：

```bash
git pull --rebase origin master
python scripts/pipeline/analyze_a1_gpiba_completed_md.py
```

如果要强制重算已有 xvg/index：

```bash
python scripts/pipeline/analyze_a1_gpiba_completed_md.py --force
```

如果只想先看会分析哪些样本：

```bash
python scripts/pipeline/analyze_a1_gpiba_completed_md.py --dry-run
```

## 默认 selection

A1-GPIbα MD 基于 `structures/1SQ0.pdb`：

```text
VWF A1:  chain A PDB resnr 505-703
GPIbα :  chain B PDB resnr 1-267
```

脚本默认使用 GROMACS selection：

```text
VWF_A1 = group "Protein" and resnr 505 to 703
GPIBA  = group "Protein" and resnr 1 to 267
```

如果服务器上的 `pdb2gmx` 编号发生变化，可以覆盖：

```bash
python scripts/pipeline/analyze_a1_gpiba_completed_md.py \
  --a1-selection 'group "Protein" and resnr 505 to 703' \
  --gpiba-selection 'group "Protein" and resnr 1 to 267'
```

## 特征定义

使用 `gmx select` 建 group，用 `gmx mindist -group -d 0.45` 统计两组之间每帧的：

- `a1_gpiba_contacts_*`：0.45 nm 内 contacts 数。
- `a1_gpiba_mindist_*_nm`：A1 与 GPIbα group 的 minimum distance。

时间窗口：

```text
first0_5    = production 0-5 ns
mid20_30    = production 20-30 ns
tail40_50   = production 40-50 ns
```

分类器友好输出：

```text
a1_gpiba_contact_retention = contacts_tail40_50 / contacts_first0_5
a1_gpiba_contact_loss_abs  = contacts_first0_5 - contacts_tail40_50
a1_gpiba_contact_loss_frac = contact_loss_abs / contacts_first0_5

md_gpiba_interface_retained_z = cohort-internal retained direction z-score
md_gpiba_interface_loss_z     = -md_gpiba_interface_retained_z
md_gpiba_interface_loss_score = contact_loss_abs / 20
```

方向解释：

- `retained_z` 高：A1-GPIbα interface 在 equilibrium MD 中更保留，偏 anti-LOF。
- `loss_z` 高或 `loss_score` 高：interface contacts 从 early window 到 tail window 明显丢失，
  偏 2M/LOF-compatible。

## 注意

`md_gpiba_interface_loss_score` 的 `/20` 只是 provisional scaling，方便后续阈值扫描；
不要把它当最终生物物理常数。真正阈值要等 2M 与 2B hard-negative 都完成后，用
precision 约束 recall 来定。

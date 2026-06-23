# 7A6O slow025 SMD 三对照结果 - 2026-06-23

## 结论

`WT + R1306W(2B) + R1374H(2M)` 的慢速 SMD 标定已经完成，每个变体 5 个副本。

这批结果**不能按原假设作为 2B 分类特征接入**：原预期是 `2B` 更容易被拉开，因此 rupture force 更低；实测方向相反，`R1306W(2B)` 的力最高，`R1374H(2M)` 的力最低。

分析脚本给出的 go/no-go 指标：

```text
AUC(2B force < 2M) = 0.00
```

因此不要把当前 `smd_2b_score = -force_z` 强行接入分类器。若继续做 SMD，需要先改反应坐标或解释为什么力方向反号。

## 已提交的小文件

```text
output/md_7a6o_smd_slow025_features.csv
output/md_7a6o_smd_slow025_features_perrep.csv
```

汇总 CSV 很小：

```text
md_7a6o_smd_slow025_features.csv        331 B
md_7a6o_smd_slow025_features_perrep.csv 654 B
```

## 三对照汇总

```text
variant  label  reps  rupture_force_pN  force_sd  work_kJmol  smd_2b_score  smd_2b_score_work
R1306W   2B     5     1114.10           22.06     479.58      -1.005        0.312
R1374H   2M     5      982.42           44.40     473.38       0.995        0.807
WT       WT     5     1047.22           33.61     497.50       0.011       -1.119
```

按原机制假设，`R1306W(2B)` 应该低于 `R1374H(2M)` 或至少低于 WT；实际为：

```text
R1306W(2B) 1114 pN > WT 1047 pN > R1374H(2M) 982 pN
```

work 轴也没有给出可用的 2B 阳性分离。

## Per-rep 数值

```text
R1306W rep1 1134.4 pN 496.0 kJ/mol
R1306W rep2 1116.0 pN 488.7 kJ/mol
R1306W rep3 1077.1 pN 506.0 kJ/mol
R1306W rep4 1126.6 pN 427.0 kJ/mol
R1306W rep5 1116.4 pN 480.2 kJ/mol

R1374H rep1 1056.9 pN 437.8 kJ/mol
R1374H rep2  950.4 pN 450.9 kJ/mol
R1374H rep3  984.1 pN 538.9 kJ/mol
R1374H rep4  947.4 pN 501.1 kJ/mol
R1374H rep5  973.3 pN 438.2 kJ/mol

WT     rep1 1049.6 pN 526.0 kJ/mol
WT     rep2  999.3 pN 420.3 kJ/mol
WT     rep3 1053.8 pN 532.1 kJ/mol
WT     rep4 1040.1 pN 522.9 kJ/mol
WT     rep5 1093.3 pN 486.2 kJ/mol
```

## 原始文件体积

三对照 slow025 原始输出在：

```text
output/gromacs_md_autoinhib/WT/smd/smd_slow025_rep*.*
output/gromacs_md_autoinhib/R1306W/smd/smd_slow025_rep*.*
output/gromacs_md_autoinhib/R1374H/smd/smd_slow025_rep*.*
```

每个变体 5 个副本，每个副本包含 `xtc/gro/edr/log/pullf/pullx`，三对照合计：

```text
files = 90
total = 1.716 GiB
xtc   = 1.597 GiB
gro   = 76.430 MiB
edr   = 16.154 MiB
log   = 18.529 MiB
pullf = 5.304 MiB
pullx = 5.258 MiB
```

按变体：

```text
WT     30 files 585.979 MiB
R1306W 30 files 585.640 MiB
R1374H 30 files 585.886 MiB
```

这些原始轨迹没有提交到 git；只提交了汇总 CSV 和本文档。

## 下一步建议

1. 不要把当前 slow025 force 轴接入 RULE6，因为方向与 2B 机制假设相反。
2. 若继续 SMD，优先检查 reaction coordinate：当前是 anchor-to-anchor AIM unfolding，可能不等价于 A1-GPIb 暴露或 AIM-A1 接触破裂。
3. 更贴近分类目标的候选特征应围绕平衡 MD 下的 AIM-A1 接触、盐桥、水暴露、A1 GPIb-binding face 暴露等，而不是强行用当前 pulling force。
4. 如果另一个 Agent 要复核，可直接看两个 CSV；需要原始轨迹时再走外部大文件传输，不走 git。

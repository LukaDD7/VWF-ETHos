# scripts/pipeline — FoldX Binding Energy Pipeline

本目录包含用于计算 VWF 变异对配体结合亲和力影响的脚本。

## 背景

VWF A1 域突变通过影响与 GPIbα（血小板受体）的结合强度产生临床表型：

| 结合变化 | 机制 | 亚型 |
|---|---|---|
| 结合增强（ΔΔG < −1.5 kcal/mol） | Gain-of-function | **Type 2B** |
| 无显著变化 | 中性 | uncertain |
| 结合减弱（ΔΔG > +1.5 kcal/mol） | Loss-of-function | **Type 2M** |

**为什么ΔΔG负值 = 结合增强**：ΔG_bind 天生是负数（结合降低系统能量），突变使其更负 = 更稳定结合 = gain-of-function。

---

## foldx_a1_gpiba_runner.py

### 功能

以 **PDB 1M10**（VWF A1 R1306Q 与 GPIbα 的共晶结构）为模板，对
`master_type1_type2.csv` 中所有 A1 域变异批量计算 ΔΔG\_bind(GPIbα)。

### 模板结构说明

| 文件 | 内容 | 链 |
|---|---|---|
| `structures/1M10.pdb` | VWF A1（R1306Q）+ GPIbα | A = VWF A1, B = GPIbα |
| `structures/1SQ0.pdb` | VWF A1（V1259M）+ GPIbα | 备用（同类结构） |

**为什么模板含 R1306Q 突变**：
WT VWF A1 在静态条件下与 GPIbα 亲和力极弱，无法直接结晶。
R1306Q 是经典 Type 2B 突变，使 A1 在无剪切力条件下即可紧密结合 GPIbα，
从而实现体外结晶。这是晶体学常用的"构象锁定"技巧，与患者病理无关。

**GPIbα 链上的 N37Q / N175Q 突变**：这是去除 N-连接糖链的工程改造（利于结晶），
与患者无关，分析时保持不变。

### FoldX 工作流

```
1M10.pdb (R1306Q 背景)
    │
    ├─ [Step 1] FoldX RepairPDB     → 填补缺失侧链、优化能量（耗时 5-15 min，只做一次）
    │
    ├─ [Step 2] BuildModel: QA543R  → 恢复 WT（Q→R，消除背景突变）→ AnalyseComplex → G_bind(WT)
    │
    └─ [Step 3] 对每个患者突变（并行）：
         BuildModel: QA543R,{mut}   → 恢复WT + 引入患者突变 → AnalyseComplex → G_bind(mut)
         ΔΔG_bind = G_bind(mut) − G_bind(WT)
```

**FoldX 突变符号格式**：`{原氨基酸}{链}{PDB位置}{新氨基酸}`
- 例：`QA543R` = Chain A 第543位 Q→R（还原背景突变）
- 例：`QA543R,RA578Q` = 同时还原背景 + 引入 R1341Q（UniProt 1341 = PDB 578）

**UniProt → PDB 位置换算**：`PDB_pos = UniProt_pos − 763`
（来源：SEQADV 记录，UniProt R1306 = PDB Chain A 位置 543）

### 环境依赖

```bash
# Python 标准库，无需额外安装
python3 >= 3.10

# FoldX 5 二进制（服务器已有）
# 路径: /inspire/hdd/global_user/mengweicheng-240108120092/lzy/tools/foldx/<binary>
# rotabase.txt 需在同一目录下
```

### 使用方法

**本地验证（不运行 FoldX，只检查突变计划）**：

```bash
cd VWF-ETHos
python3 scripts/pipeline/foldx_a1_gpiba_runner.py --dry-run
```

输出示例：
```
A1 domain variants loaded: 82
  Valid (in PDB range):  73
  Skipped:               9
    SKIP VWF_M1303M: synonymous
    SKIP VWF_Q1311X: stop_codon
    ...

=== DRY RUN: FoldX mutation plan ===
Variant                    UniProt    PDB Subtype    FoldX notation
VWF_R1306Q                    1306    543 type2B     QA543Q;   ← 背景本身就是Q，预期 ΔΔG ≈ 0
VWF_R1306W                    1306    543 Type2B     QA543W;   ← 已知2B，预期 ΔΔG < −1.5
VWF_F1369I                    1369    606 type2M     QA543R,FA606I;  ← 已知2M，预期 ΔΔG > +1.5
```

**服务器执行**：

```bash
# 1. 确认 FoldX 路径（注意具体二进制名可能不同）
ls /inspire/hdd/global_user/mengweicheng-240108120092/lzy/tools/foldx/
# 找到类似 foldx、foldx5、foldx_20251231 的可执行文件

# 2. 运行（16核，预计 20-40 分钟）
python3 scripts/pipeline/foldx_a1_gpiba_runner.py \
    --variants data/processed/master_type1_type2.csv \
    --pdb structures/1M10.pdb \
    --foldx-bin /inspire/hdd/global_user/mengweicheng-240108120092/lzy/tools/foldx/<binary_name> \
    --output output/foldx_a1_gpiba.csv \
    --workers 16

# 3. 推荐用 tmux 挂后台
tmux new -s foldx-a1
python3 scripts/pipeline/foldx_a1_gpiba_runner.py [参数...]
# Ctrl+B then D 离开
```

### 参数说明

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--variants` | `../../data/processed/master_type1_type2.csv` | 变异表 |
| `--pdb` | `../../structures/1M10.pdb` | 模板 PDB |
| `--foldx-bin` | 服务器路径 | FoldX 可执行文件完整路径 |
| `--output` | `../../output/foldx_a1_gpiba.csv` | 输出 CSV |
| `--work-dir` | `output/foldx_workdir_a1/` | FoldX 中间文件目录 |
| `--workers` | `8` | 并行进程数（建议设为 CPU 核数） |
| `--dry-run` | — | 只打印计划，不执行 FoldX |

### 输出格式

`output/foldx_a1_gpiba.csv`：

```
variant_id, uniprot_pos, wt_aa, mut_aa, subtype, g_bind_wt, g_bind_mut, ddg_bind, direction, status, foldx_notation
VWF_R1306W, 1306, R, W, Type2B, -45.231, -47.890, -2.659, GOF, OK, QA543W
VWF_F1369I, 1369, F, I, type2M, -45.231, -43.012, +2.219, LOF, OK, QA543R,FA606I
```

| 列 | 说明 |
|---|---|
| `ddg_bind` | ΔΔG\_bind (kcal/mol)，负=增强，正=减弱 |
| `direction` | `GOF` / `LOF` / `neutral` / `unknown` |
| `status` | `OK` / `FAIL_FOLDX` / `ERROR:...` |
| `foldx_notation` | FoldX BuildModel 的突变符号（可直接核验） |

### 校准验证

**跑完后立即检查已知变异**：

```bash
python3 -c "
import pandas as pd
df = pd.read_csv('output/foldx_a1_gpiba.csv')
known = ['VWF_R1306W', 'VWF_R1306Q', 'VWF_R1341Q', 'VWF_F1369I', 'VWF_R1374C']
print(df[df['variant_id'].isin(known)][['variant_id','ddg_bind','direction','subtype']])
"
```

预期结果（方向正确即通过）：

| 变异 | 已知分型 | 预期 direction | 预期 ΔΔG 符号 |
|---|---|---|---|
| R1306W | Type 2B | GOF | 负 |
| R1306Q | Type 2B | GOF | 负 |
| R1341Q | Type 2B | GOF | 负 |
| F1369I | Type 2M | LOF | 正 |
| R1374C | Type 2M | LOF | 正 |

若已知2B变异方向错误（出现LOF），说明阈值需要重新校准。
调整 `THRESHOLD_GOF` / `THRESHOLD_LOF` 常量后重新运行。

### 覆盖范围说明

1M10 结构仅覆盖 VWF UniProt 位置 **1268–1466**，以下变异自动跳过：

| 变异 | 原因 |
|---|---|
| Y1258C, P1266L/Q | 位置 < 1268（PDB 覆盖边界外） |
| T1477I, D1472A/H | 位置 > 1466（PDB 覆盖边界外） |

跳过的变异建议改用 **Boltz-2** 预测复合物结构后再做 FoldX AnalyseComplex。

### 常见问题

| 症状 | 原因 | 解决 |
|---|---|---|
| `FoldX binary not found` | 路径错误 | `ls` 确认二进制名和路径 |
| `RepairPDB failed` | rotabase.txt 缺失 | 确认 rotabase.txt 与二进制在同一目录 |
| ΔΔG 超出 ±15 kcal/mol | 结构局部冲突 | 检查 work_dir 中对应 job 的 build.log |
| 已知 2B 变异显示 LOF | 阈值偏差 | 调整 THRESHOLD_GOF 后重跑 |
| `parse_analyse_complex` 返回 None | fxout 列格式不同 | 检查 FoldX 版本，对应调整列索引 |

# 7A6O AIM-A1 MD: Feature Analysis & Classifier Integration (2026-06-20)

A40 服务器把 7A6O AIM-A1 自抑制 MD 结果 push 回仓 (`md_data/7a6o_reference_md/`,
Git LFS, WT + 14 参考变体 × 50 ns)。本文记录：医学机制调研、MD 结果分析、特征
提取、以及并入 `agentic_vwf_classifier.py` 的方式与**符号纠偏**。

## 1. 数据

- `md_data/7a6o_reference_md/variants/<v>/{prod_concat.xtc, md_prod.tpr, final.gro}`
- 15 体系 (各 501 帧 / 0–50 ns)，WT 校验和与 `manifest.csv` 一致 (LFS 完整)。
- 标签 (参考集)：WT；2B = R1306W/R1306Q/R1308C/I1309V/S1310F/W1313C/V1314F/V1316M (8)；
  2M = R1374C/R1374H/G1324S (3)；uncertain = P1337L/R1341Q/R1341W (3)。
- 构件原生编号 **1262–1466** (经 12 个突变位点逐一比对确证；
  `native = MDAnalysis local_resid + 1261`)。AIM = N-AIM 1262–1271 + C-AIM 1459–1466；
  A1 体 ≈ 1272–1458。

## 2. 医学机制调研 (确保流程符合医学特性)

- A1 域被一个**不连续自抑制模块 (AIM)** 屏蔽：N/C 端侧翼序列协同折叠，盖住 A1 的
  **GPIbα 结合面 (helix α1、loop α1β2、loop β3α2)**。(Deng et al., *Blood* 2017, HDX)
- **2B = 功能获得 (GOF)**：突变**降低 AIM 解折叠所需的力**，剪切力下 AIM 展开暴露 A1 →
  自发结合 GPIbα。关键：2B 的结合面本身**完好**（否则无法增强结合）；其释放是
  **力依赖**的。(Nat Commun 2021；*Blood* 2023 "differentially perturb"——2B 之间异质)
- **2M = 功能丧失 (LOF)**：突变直接**破坏 A1–GPIbα 结合面**，结合下降。

**推论 (决定特征方向)**：无剪切的 50 ns 平衡 MD **看不到力依赖的 2B 松开**。平衡态下
医学上稳健的可观测量是 **"变体保留多少 WT 的 AIM→结合面屏蔽接触网"** (= 结合面完整性)，
而**不是** AIM↔A1 总接触数 (后者方向相反、有误导性)。

## 3. MD 分析结果

### 3.1 旧的总接触指标方向相反 (弃用为定型轴)

复算 `aim_a1_contacts`（AIM↔A1 体 1280–1451，0.45 nm）复现了既有 CSV：
**2B 接触更多 / mindist 更短**，**2M (R1374C/H) 接触更少 / mindist 更长**。
这与"接触下降 = 2B 松开"的旧假设**符号相反**——故总接触数不能作 2B 定型轴。

### 3.2 数据驱动屏蔽界面 → 结合面完整性 (采用)

由 WT 轨迹取占用率 ≥ 50% 的 AIM→A1 接触对，得屏蔽界面 A1 残基：
**1305, 1306, 1307, 1308, 1313, 1376, 1378, 1380, 1410, 1434**。
其中 1306/1308/1313 **正是 2B 突变位点** → 验证了 Deng 2017 的屏蔽模型。

各变体 tail (40–50 ns) 保留的屏蔽接触比例 `md_aim_mask_retention`：

| 标签 | 中位/均值 | 区间 |
|---|---|---|
| WT | 0.81 | — |
| 2B | 0.75 | 0.60–0.89 |
| 2M | 0.63 | 0.55–0.70 |

排序 **WT > 2B > 2M** 与机制自洽：**2M 结合面被破坏 → AIM 也屏蔽不住 → 保留率最低**；
2B 结合面完好 → 平衡态保留接近 WT。

### 3.3 判别力 (诚实)

- `md_face_destab_score = -zscore(mask_retention)`（高 = 结合面破坏 = 2M/LOF）。
- **AUC(2M > 2B) = 0.83**；阈值 z ≥ 1.0 命中 R1374H、V1316M。
- **局限**：单 50 ns 副本、n(2M)=3、2B/2M 有重叠 (V1316M 这个 2B 也落低)。
  → 只能作**软证据 / 2M 确证轴**，**不能判 2B**，不硬翻已定型的 2B。

## 4. 特征产物 (整理出的特征)

可复现脚本：`scripts/pipeline/extract_7a6o_md_features.py`
（只读轨迹，MDAnalysis；自动定义 WT 屏蔽界面 + 逐变体提取）。

- `output/md_7a6o_features.csv` — 每变体：`md_aim_mask_retention`, `md_mask_z`,
  `md_face_destab_score`, `ncont_heavy_*`, `aim_a1_mindist_*`, `aim_rmsf_*`。
- `output/md_7a6o_masking_interface.json` — 屏蔽界面残基、偏移、z 标定参数。

```bash
python3 scripts/pipeline/extract_7a6o_md_features.py \
    --input md_data/7a6o_reference_md/variants \
    --output output/md_7a6o_features.csv
```

## 5. 并入分类系统 (`agentic_vwf_classifier.py`)

新增**轴B' = MD A1 结合面完整性**（LOF/2M 证据），常数 `MD_FACE_DESTAB_2M_Z = 1.0`：

- `ExpertScores.md_face_destab_score`（NaN = 无 MD，向后兼容）。
- `RULE6`（A1 域）中：
  1. 若轴B（forced_binding+heparan）已判 2M，且 MD 结合面强破坏 → **提升置信** 0.72→0.8（确证）。
  2. 在原 `uncertain` 兜底前加 **MD tie-breaker**：轴 A/B/C 与结构启发**皆未定向**时，
     MD 结合面强破坏 (`destab ≥ 1.0`) → 软证据判 **2M (0.55)**。因排在 2B 热点/AIM 位/
     别构判定**之后**，**绝不硬翻已定的 2B**。
- **纠偏注释**：在 `AIM_RELEASE_2B_Z` 处加警示——勿用平衡 MD 接触下降量反推
  `aim_release_score`（符号相反）；2B release 轴仍需 panel CIF 或 steered/force MD 定标。

### 5.1 参考集回归验证

15 变体过分类器（有 MD vs NaN）：

- **2B→2M 误翻 = 0**（V1316M destab=+1.16 越阈，但 1316∈热点先判 2B → 安全顺序生效）。
- **R1374H (2M)** 由 `uncertain` 被 MD 正确救回 **2M**。
- 其余与 NaN 完全一致（向后兼容）。

## 6. 已知遗留 / 下一步

- `TWO_B_HOTSPOT_POS` 含 1324/1341（本集 G1324S=2M、R1341Q/W=?），属**既有**热点先验
  误伤，与本次 MD 改动无关；需按本地标签校准（代码内已标注）。本次未动。
- MD 轴待**多副本 + 更多 2M 标签**或 **steered/force MD**（看力依赖松开）后收紧阈值，
  届时方可考虑真正的 2B(GOF) MD 轴。
- 患者 C1458R 的 3 副本未纳入参考集（`manifest.csv` 注明）；其分析见
  `docs/7A6O_C1458R_PATIENT_MD_GUIDE_2026-06-20.md`。

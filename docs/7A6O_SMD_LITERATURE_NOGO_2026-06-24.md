# 7A6O SMD 力轴 no-go：三对照结果 + 文献复核 — 2026-06-24

承接 `docs/7A6O_SMD_SLOW025_READOUT_2026-06-23.md`（slow025 三对照完成）与
`docs/7A6O_SMD_WT_DIAGNOSIS_2026-06-21.md`（WT 性能根因 + go/no-go 判据）。

**一句话结论**：slow025 三对照（WT / R1306W-2B / R1374H-2M，各 5 reps 独立初速度）
力轴方向与 2B 机制假设**完全反号**（`AUC(2B force < 2M) = 0.00`）；文献复核证实
**是 SMD 设置问题，不是生物学反了**。当前 anchor-to-anchor AIM 解折叠力轴 **不接入分类器，
也不靠翻符号硬上**。退回平衡态观测量另起一轴。

---

## 1. 三对照数据（slow025 = 0.25 nm/ns，各 5 reps）

| 变体 | 标签 | 断裂力 (pN) | SD | 单 rep 范围 | 功 (kJ/mol) |
|---|---|---|---|---|---|
| R1306W | 2B | **1114.1** | 22.1 | 1077–1134 | 479.6 |
| WT | WT | 1047.2 | 33.6 | 999–1093 | 497.5 |
| R1374H | 2M | **982.4** | 44.4 | 947–1057 | 473.4 |

- 5 reps 为**独立采样**：同一平衡结构出发，各自 `gen_vel=yes` + 不同种子
  （`12345+rep×7919`）。断裂力随机，故重复采样取均值±SD 是正确做法。
- 两个病例组分布几乎不重叠（2B SD 仅 22 pN）→ **方向反号是稳健可复现的，不是噪声**。
- 实测排序：`R1306W(2B) 1114 > WT 1047 > R1374H(2M) 982`，与"2B 更易解折叠（力更低）"相反。
- work 轴同样无可用分离。

### `AUC(2B force < 2M) = 0.00` 的含义
随机取一个 2B rep 与一个 2M rep，2B 力更低（符合假设）的比例 = 0。
- 0.5 = 分不开；1.0 = 完美正向分开；**0.00 = 完美反向分开**。
- 不是"分不开"，是"分得很开但方向反了"。可以翻符号变成可用分类器，**但前提是信这套物理**
  —— 文献证明不能信（见 §3），翻符号 = 在错误物理上过拟合。

---

## 2. 文献复核：别人做过 WT/2B/2M 的 AIM 力学解折叠吗？

做过，多组，结论**与本结果相反**：

| 来源 | 方法 / 几何 | WT 力 | 2B 突变效应 |
|---|---|---|---|
| Springer, Nat Commun 2021 (PMC8060278) | 光镊，N↔C 端拉伸（与我们同几何） | **10–20 pN** | H1268D / R1341Q **显著降低**解折叠力；R1341Q ~90% trace 为单一低力事件 |
| Catch-bond SMD, RSC Chem Biol 2022 (PMC9175105) | A1–GPIbα 复合体拉伸 | — | R1306Q / R1450E 打断 N-AIM↔α1/α6 盐桥，**零力**下即结合 → 2B 更易释放 |
| Springer, eLife 2022 (eLife 75760) | 平衡 + 力学 | — | 明确反相关：2B(GOF) 增亲和/降稳定性；2M(G1342S) **增**稳定性 |

→ 学界一致：**2B = 解折叠力更低 / 更不稳定**。我们的 R1306W 最高 → 我们方向反了。

---

## 3. 为什么是 SMD 设置问题（4 个具体根因）

1. **力大 50–100×**：实验 WT 10–20 pN，我们 1047 pN。即便降到 0.25 nm/ns，仍比光镊快
   ~10⁶ 倍 → 断裂被**黏滞摩擦 + 主链整体解折叠**主导，2B 那点盐桥 ΔΔG 被埋掉。
   降速（1→0.25）减了摩擦但没出这个区，故方向始终不翻。
2. **R1306W 探针选错**：R1306 形成 D1269–R1306 自抑制盐桥，2B 应**去稳定**；但替换为
   **色氨酸**（大疏水侧链），折叠态新增堆积/位阻，在快速拉伸下反而**抬高**机械阻力
   （拽大侧链的假象），掩盖甚至反转盐桥释放效应。⚠ 文献 2B 对照都是**去电荷型**
   （R1306**Q** / R1341Q / H1268D），不是加体积型。
3. **反应坐标 ≠ 机制**：端到端解折叠测整体机械稳定性；2B 真实缺陷是**低/零力下提前结合
   GPIbα**（释放/结合事件）。1000 pN 已远超该区。
4. **N=2 病例对照**：1 个 2B + 1 个 2M，即便方向对也无法标定。

---

## 4. 决定 & 待办

- ❌ **当前 SMD 力轴（含 work 轴）不接入 `agentic_vwf_classifier.py` RULE6**；不翻符号硬上。
- ⛔ SMD 这条线到此暂停。若日后要救：需用**去电荷型 2B 对照**（R1306Q）、**更慢加载率**
  （显溶剂 SMD 实际难达 10–20 pN），或改反应坐标（AIM 释放/GPIbα 暴露而非整体解折叠）。
- ✅ **退回平衡态观测量另起一轴**（文献都依赖这些，且与现有 2M `md_face_destab_score` 正交）：
  - D1269–R1306 等 **AIM↔A1 盐桥占据率**；
  - **AIM 解离 / 解折叠比例**（接触数、RMSD）；
  - **A1 的 GPIbα 结合面暴露**（SASA）。
  - 数据来源：autoinhib 平衡 MD 轨迹（`output/gromacs_md_autoinhib/`）。
- 下一步：盘点平衡 autoinhib 轨迹可抽哪些上述特征，评估这条轴可行性。

---

## 5. 数据 / 文件

- 汇总：`output/md_7a6o_smd_slow025_features.csv`（3 行）
- 单 rep：`output/md_7a6o_smd_slow025_features_perrep.csv`（15 行）
- 原始轨迹（**未入 git**，1.716 GiB）：`output/gromacs_md_autoinhib/{WT,R1306W,R1374H}/smd/smd_slow025_rep*.*`
- 分析脚本：`scripts/pipeline/analyze_7a6o_smd.py`（已输出 AUC 自检）
- SMD runner：`scripts/pipeline/run_7a6o_smd.sh`

---

## 6. 平衡态替代轴尝试：定向自抑制盐桥（2026-06-24，本地跑）

退回平衡态后，按文献（PMC9175105）点名的自抑制盐桥写了
`scripts/pipeline/extract_aim_saltbridge_features.py`，在本地 15 条轨迹
（`md_data/7a6o_reference_md/variants/`，MDAnalysis 2.10.0，tail 40–50 ns）上跑。
输出 `output/md_7a6o_saltbridge_features.csv`。

| 特征 | WT | 2B (n=6–8) | 2M (n=3) |
|---|---|---|---|
| D1269–R1306 盐桥占据率 | 0.667 | 0.105 | 0.007 |
| D1269–R1450 盐桥 | 0.000 | 0.005 | 0.000 |
| AIM↔A1 总盐桥数 (tail) | 1.22 | 0.88 | 0.12 |

**结果：平衡态没有干净的 2B 阳性轴。**
- **D1269–R1306**：WT 0.67（自抑制完好）→ 2B 与 2M **都**塌到 ~0。是干净的
  **WT vs 致病** 标志，**不分 2B/2M**。（R1306Q/W = NaN，突变本身删了正电荷。）
- **D1269–R1450**：平衡态全程 ~0，未采样到该 α6 接触，无用。
- **AIM↔A1 总盐桥数**：2M 紧密全塌（≤0.18），2B/WT 保留较多 → 又是一条 **2M/LOF 轴**，
  与现有 `md_face_destab_score` 冗余；2B 内部很散（V1314F 0.08 / V1316M 0.0 掉进 2M 区）。

→ **再次印证**：2B 本质是"力"现象，平衡态下与 WT 难分；SMD 与平衡盐桥两路都未单独拎出 2B。

**唯一浮现的可用结构：联合判据（非新 2B 阳性轴）**
在"D1269–R1306 自抑制已释放（占据 ~0）"的前提下，按**结合面是否保留**二分：
- SB 数近 WT（保留）→ 2B（释放但功能在 = GOF）；
- SB 数 ≤0.18（也塌）→ 2M（释放且功能没了 = LOF）。
2B 内部仍有噪声（V1314F/V1316M 会误判 2M），需更多 2M 对照与位点先验联合。

**待办（下一步候选，未定）**：
- (a) 把"自抑制释放 + 结合面保留"联合判据形式化进分类器，与临床 2B 热点先验联用；
- (b) 或承认结构/平衡 MD 难独立定 2B，2B 主要靠临床热点位置 + 排除 2M 来判；
- (c) 若仍要力学信号：换去电荷型 2B 对照 + 改反应坐标（见 §4），代价高、收益不确定。

---

## 7. 联合判据已接入分类器（2026-06-24，用户选定方向 a）

把 §6 的"自抑制释放 + 结合面保留 → 2B"联合判据形式化进
`scripts/agentic_vwf_classifier.py` 的 RULE6（A1 域），与临床 2B 热点先验联用。

**新特征**：`aim_sb_retained_z`（`extract_aim_saltbridge_features.py` 产出，
AIM↔A1 盐桥保留数的 z-score，ref=WT+2B+2M）。`ExpertScores` 新字段，
`variant_data.get('aim_sb_retained_z', np.nan)` 注入，NaN 时 RULE6 退回原逻辑（向后兼容）。

**阈值**（暂定，待校准；n: WT1/2B8/2M3 单副本）：
- `AIM_SB_RETAINED_2B_Z = 0.6`：≥ → 结合面保留 → **anti-2M / 支持 2B**（实测高 z 仅 2B+WT）。
- `AIM_SB_COLLAPSE_2M_Z = -0.7`：≤ → 结合面塌陷 → **2M 旁证**（低 z 不干净含 2B，故不独立判 2M）。

**RULE6 三处接入**（加性、不硬翻已定的 2B）：
1. 轴B(结合丧失)判 2M 时，`sb_collapsed` 作旁证 → 置信 +0.05（封顶 0.85）。
2. 临床 2B 热点 + `face_retained` → **联合判据** RULE6-C+保留，置信 0.60→0.72。
3. 软 MD 2M tie-breaker 前加 `face_retained` 闸：保留 → 不翻 2M，2B marginal(0.55)。

**功能自测（8 例）全过**：联合判据 2B(0.72)、热点-only 回归 2B(0.60)、face 闸挡住
soft-2M 翻为 2B(0.55)、md_lof-only 回归 soft-2M(0.55)、binding_lost+塌陷 2M(0.77)、
无 MD → uncertain(0.35) 回归不变。

**接入方式**：`aim_sb_retained_z` 与 `md_face_destab_score` 同惯例，靠验证时外部 join
`output/md_7a6o_saltbridge_features.csv`（按 variant 名）进矩阵，无独立 merge 脚本。

**仍待办**：① 在更大本地标签集上跑分类器验证 2B recall 变化（基线 2/12）；
② 校准 `AIM_SB_RETAINED_2B_Z` / `AIM_SB_COLLAPSE_2M_Z`（现 n 太小）；
③ V1314F/V1316M(2B) 仍可能被低 z 误伤，需位点先验兜底。

---

## 参考文献
- Nat Commun 2021 — Activation of VWF via mechanical unfolding of its discontinuous AIM. https://pmc.ncbi.nlm.nih.gov/articles/PMC8060278/
- RSC Chem Biol 2022 — N-terminal AIM stabilizes the mechanosensor catch bond. https://pmc.ncbi.nlm.nih.gov/articles/PMC9175105/
- eLife 2022 — A1 affinity & stability differentially regulated by O-glycosylated linkers. https://elifesciences.org/articles/75760
- JPC B 2025 — Dynamics of a VWF A1 AIM with O-linked glycans and GPIbα regulation. https://pmc.ncbi.nlm.nih.gov/articles/PMC12010329/

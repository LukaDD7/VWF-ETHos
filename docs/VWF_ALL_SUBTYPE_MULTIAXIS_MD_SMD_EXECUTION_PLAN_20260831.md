# VWF 全亚型多机制轴 MD/SMD 执行方案

日期：2026-08-31  
状态：方案已细化；等待 GPU 释放后执行  
适用范围：VWF 完整诊断系统的 Type 1 / low VWF、2A、2B、2M、2N、3 型机制轴

---

## 1. 结论

1. **当前 16 例先跑 1 型和 2B 直接相关的轴，但系统方案必须覆盖全 VWD 亚型。**
   - 当前批次只有 Type 1 和 Type 2B，因此执行队列先做 A2、D′D3、A1/AIM、A1-GPIbα。
   - 完整系统同时保留 2A 的 A2/ADAMTS13 与多聚化轴、2M 的 A1-GPIbα/A3 胶原轴、2N 的 D′D3-FVIII 轴、3 型的 null/CNV/regulatory 轴。
2. **结构来源必须分层。**
   - 有实验结构：A1/AIM 用 `7A6O`，A1-GPIbα 用 `1SQ0`（备用 `1M10`），A2 用 `3GXB`，A3-collagen 用 `4DMU`，D′D3 用 `6N29`，CK/CTCK 用 `4NT5`。
   - 没有实验结构的 D4、C-domain、D1/D2 才允许 Boltz/AlphaFold fallback，并标记 exploratory。
3. **平衡态 MD 是默认层，SMD 是条件层。**
   - 旧 AIM-A1 end-to-end SMD 已经判 no-go，不能复用为分类特征。
   - A2 机械解折叠、A1-GPIbα 解离只有在实验结构、匹配对照、反应坐标和加载速率都合理时才允许进入。
4. **AI 参考分布必须按 subtype + module + assay + feature 建立。**
   - 不做全局阈值。
   - WT、良性对照、已知亚型代表和边界 case 分层校准。
   - 低置信度预测结构的 MD 只能输出 mechanism support，不能输出强临床分型结论。

---

## 2. 全亚型机制轴框架

| 亚型 | 直接相关机制轴 | 首选结构/数据 | MD 策略 | 关键读数 | 主要 caveat |
|---|---|---|---|---|---|
| Type 1 / low VWF | secretion/clearance、splice/NMD/CNV、domain folding/assembly | 7A6O、3GXB、6N29、4NT5；无实验结构模块用 Boltz | 只有 missense 且结构可靠时做 20 ns pilot；null/splice/CNV 走非结构轴 | RMSD、二级结构、局部柔性、domain packing、VWFpp/VWF:Ag | 异质性大，不能把所有 1 型都推成结构 MD |
| Type 2A | A2 destabilization/ADAMTS13；multimerization/assembly | 3GXB、4NT5、7WN3 如后续下载；D4/C-domain 仅 fallback | A2 folded MD 优先；assembly 模块二级；VWF73-ADAMTS13 需模型验证 | A2 稳定性、Tyr1605-Met1606 暴露、多聚化表面、domain packing | 2A 既可来自 ADAMTS13 敏感，也可来自多聚化/分泌缺陷 |
| Type 2B | A1/AIM release；A1-GPIbα gain of binding | 7A6O、1SQ0、1M10 | A1/AIM equilibrium + A1-GPIbα equilibrium；SMD 仅条件启动 | AIM 接触、盐桥、A1 暴露、GPIbα 界面保留 | A1 release 需与 RIPA/血小板表型联合解释 |
| Type 2M | A1-GPIbα LOF；A3-collagen LOF；closed-A1 dynamics | 1SQ0、7A6O、4DMU | A1-GPIbα 和 A3 equilibrium；closed A1 作为正交旁证 | 界面接触、界面 RMSD、接触寿命、胶原结合面保留 | contact loss 是 LOF/2M 方向，不是 2B 阳性方向 |
| Type 2N | D′D3-FVIII binding | 6N29、2MHP；FVIII 复合物另行验证 | D′D3 folded stability 先行；FVIII complex 只有在构建策略明确后做 | FVIII 结合面保留、SASA、残基对距离、局部柔性 | 需与 hemophilia A 和联合缺陷区分 |
| Type 3 | null、large deletion/duplication、promoter/regulatory、inhibitor risk | 变异坐标、CNV/RNA/家系证据 | 非 MD-first | NMD、CNV、剪切、剂量、表达缺失 | 大缺失/null 不进入单域 MD |

机器可读版本：`outputs/computational_panel_20260829/md/vwd_subtype_mechanism_axis_plan.csv`。

---

## 3. 当前 16 例的实例化队列

当前数据为 16 名患者、17 条变异记录：Type 1 10 例，Type 2B 6 例；T2B_005 另含一条报告为良性的 D2449N。逐样本表见：

```text
outputs/computational_panel_20260829/md/sample_mechanism_axis_plan.csv
```

### 3.1 已完成

A1/AIM 7A6O 匹配系统已经完成 50 ns：

- `PANEL_WT_MATCHED`
- `P1413L`
- `R1308C`
- `S1310F`
- `V1316M`
- `R1341W`
- `A1461D`

这覆盖当前 Type 2B 六例主变异和 Type 1 的 `P1413L`。`V1316M` 两个病例共用同一变异结果。

### 3.2 P0：先跑实验结构直接相关轴

1. **A2 / A1500V**
   - 初始结构：`structures/3GXB.pdb`
   - 对照：WT 3GXB
   - 队列：WT + A1500V，20 ns pilot；QC 通过后升 50 ns
   - 读数：backbone RMSD、二级结构保留、局部柔性、Tyr1605-Met1606 切割位点暴露
   - 目标：Type 1 样本的 A2 folded stability / ADAMTS13 敏感性证据；同时是 2A 框架的校准轴
2. **D′D3 / R1205H**
   - 初始结构：`structures/6N29.pdb`
   - 对照：WT 6N29
   - 队列：WT + R1205H，20 ns pilot；QC 通过后升 50 ns
   - 读数：fold stability、FVIII 结合面保留、局部柔性、关键残基对距离
   - 目标：该样本 FVIII:C 极低，优先补 D′D3-FVIII 轴；同时校准 2N 框架
3. **A1-GPIbα / 当前 2B 变异**
   - 初始结构：`structures/1SQ0.pdb`，备用 `1M10`
   - 对照：WT 1SQ0、已知 2M LOF 对照、至少一个非 2B 对照
   - 变异：R1308C、S1310F、A1461D、R1341W；V1316M 已有既有结果需先核对复用
   - 读数：界面接触占有率、界面 RMSD、最小距离、接触寿命、A1 表面暴露
   - 目标：补齐 2B gain-of-binding 的正方向动态证据，而不是只看 AIM release

### 3.3 P1：预测结构探索轴

1. **D4 / C1950Y、R1943C**
   - 结构：D4 无本地实验结构，只能用 assay-matched Boltz/AlphaFold fallback
   - 输出层级：exploratory
   - 队列：WT + C1950Y + R1943C，20 ns pilot
   - 读数：RMSD、二级结构、局部柔性、domain packing
   - 不能单独作为临床分型证据；优先和分泌、多聚体、VWFpp/VWF:Ag 联合解释
2. **C-domain / D2449N**
   - 报告为良性，不建议优先占用生产 GPU
   - 仅作为 benign/reference calibration 或阴性对照

### 3.4 不进入结构 MD 的样本

- `CASE_T1_003`：G2488 同义变异，走 splice/translation 证据
- `CASE_T1_004`：大片段缺失/重复，走 CNV/dosage
- `CASE_T1_005`：Y988Ter，走 NMD/null
- `CASE_T1_006`：L2407fs，走 NMD/null
- `CASE_T1_007`：c.3379+1G>A，走 splice

这些病例不是“跳过”，而是切换到非结构机制证据链。

---

## 4. GPU 2/3 执行策略

当前 GPU 2/3 仍在跑 Boltz-2 keepalive。不要直接杀任务；释放后按以下方式并行：

| GPU | 队列 | 第一批 | 第二批 |
|---|---|---|---|
| GPU 2 | 实验结构高优先级 | 3GXB WT + A1500V | 6N29 WT + R1205H |
| GPU 3 | 2B 直接相关轴 | 1SQ0 WT + R1308C/S1310F/A1461D/R1341W 逐个或小批量 | D4 exploratory pilot |

执行规则：

1. 每张 H200 先保持一个 MD worker，避免 CPU/GPU 超订。
2. 每个 assay 先完成 WT 稳定性和 QC。
3. 所有新模块先 20 ns pilot，通过后再 50 ns。
4. Boltz CIF 不能覆盖实验 PDB。
5. 输出只保留轻量 CSV/JSON/manifest 和必要代表结构，不把原始轨迹直接入 Git。

---

## 5. AlphaGenome 与 Boltz 的补充规则

### AlphaGenome

当前 16 例请求已完成，本轮不默认重跑。新样本或 splice/noncoding 样本需要补充时：

1. 合并三个 splice scorer：
   - `SPLICE_SITES`
   - `SPLICE_SITE_USAGE`
   - `SPLICE_JUNCTIONS`
2. 保留 REF/ALT、strand、tissue/biosample 和 missingness reason。
3. AlphaGenome 是 regulatory/splicing 模型，不能提供蛋白初始结构。

### Boltz

Boltz-2 已覆盖 15 个功能轴和 597 个变异。本轮不再为已有样本重跑 Boltz；只在新样本或新 assay 出现时扩展。Boltz 结构仅用于：

1. 无实验结构模块的 fallback；
2. 静态置信度和模型选择；
3. 与 MD 动态特征互补。

不能把 Boltz CIF 当成实验晶体结构，也不能把 iPTM/pTM 当结合自由能或临床分类。

---

## 6. AI 参考分布

每个 feature 按：

```text
subtype + module + assay + feature -> reference distribution
```

最低参考集：

1. matched WT；
2. ≥3 个同域良性/人群对照；
3. ≥3 个已知亚型代表；
4. 边界 case，例如 2B/2M hard negative。

输出：

- mean ± SD、median/IQR；
- z-score 和 empirical percentile；
- Cohen’s d 或 rank-biserial effect size；
- bootstrap 95% CI；
- direction 和 missingness reason；
- `experimental` / `predicted` / `exploratory` 结构层级。

Agent v3 只消费结构化特征，不直接消费原始轨迹。

---

## 7. 文献与结构依据

1. Legan et al., Blood 2023：Type 2B mutations differentially perturb A1 autoinhibition；支持 2B 的 A1/AIM 与 GPIbα 双轴设计。
2. Synthetic and structural studies of VWF A2, including the A2 crystal structure and A2 mutation studies：支持 A2 folded stability、force-dependent unfolding 和 ADAMTS13 susceptibility 轴。
3. Zhou et al., Blood 2019 / PDB 6N29：D′D3 assembly 与 FVIII binding structural principles；支持 D′D3-FVIII/2N 轴。
4. A3-collagen 结构与 VWF collagen-binding studies / PDB 4DMU：支持 2M 或胶原结合缺陷轴。
5. D4/C-domain VWD mutation studies：支持 D4/C-domain 在表达、分泌、多聚化和 Weibel-Palade body 形成中的机制作用；但由于无独立实验 D4 结构，MD 只能 exploratory。
6. AlphaGenome splicing variant scoring 文档：支持 splice 变异的三个 scorer 和 merged-score 使用。

补充入口：

- <https://ashpublications.org/blood/article/141/10/1221/493944/Type-2B-von-Willebrand-disease-mutations>
- <https://pmc.ncbi.nlm.nih.gov/articles/PMC3585066/>
- <https://pmc.ncbi.nlm.nih.gov/articles/PMC6450429/>
- <https://www.rcsb.org/structure/6N29>
- <https://www.rcsb.org/structure/2MHP>
- <https://www.rcsb.org/structure/4DMU>
- <https://www.alphagenomedocs.com/colabs/splicing_variant_scoring.html>

---

## 8. 一句话总结

**当前 16 例先跑 A2、D′D3、A1/AIM 与 A1-GPIbα；完整系统则同时保留 2A 的 A2/多聚化、2M 的 A1-GPIbα/A3、2N 的 D′D3-FVIII 和 3 型的 null/CNV/regulatory 证据链。实验结构优先，预测结构只做 fallback，SMD 只有通过 go/no-go 后才进入分类器。**

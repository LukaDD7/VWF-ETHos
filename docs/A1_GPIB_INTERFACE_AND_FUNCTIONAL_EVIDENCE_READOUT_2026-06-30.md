# A1-GPIb 复合体特征与功能证据层 readout

2026-06-30。本文接在去掉 recurrent/hotspot 先验后的 Eval 之后，用于说明三档 2B/2M
机制证据如何继续补强：第一档为 Boltz 静态复合体结构再挖掘，第二档为 A1-GPIb 或
AIM-A1 equilibrium MD，第三档为功能实验/文献证据层。

## 1. 什么是“轴”

这里的“轴”不是标签，也不是某个 subtype 先验，而是一类机制证据维度。例如：

- A1-GPIb 结合轴：变体是否保留或增强 A1 与 platelet GPIbα 的结合界面。
- AIM-A1 自抑制轴：变体是否让 A1 周围的自抑制模块更容易松开，从而暴露 GPIbα 结合面。
- secretion / ER-retention 轴：变体是否导致折叠、分泌或细胞内滞留问题。
- FVIII-binding 轴：D'/D3 相关变体是否破坏 FVIII 结合，对应 2N。

原来的 `fb_binding_zscore` 更接近一个粗粒度的 Boltz complex-confidence proxy：
它主要利用 A1-GPIb 复合体预测的全局置信度/同轴分布信号。这个信号有用，但太粗；
它没有明确回答“接口接触是否还在、接口 PAE 是否低、接口残基是否稳定”。

## 2. 第一档：从 Boltz 复合体结构榨取接口特征

已新增脚本：

- `scripts/pipeline/extract_a1_gpiba_interface_features.py`
- `scripts/pipeline/download_a1_gpiba_hf_results.py`

输入：

- Boltz A1-GPIb 复合体预测目录。
- 每个 model 的 `.cif` 复合体结构、`pae_*.npz`、`plddt_*.npz`。
- 可选 Eval CSV 只用于读出时 join `true_label`；特征计算本身不需要 subtype 标签。

输出：

- `output/boltz2_a1_gpiba_analysis/a1_gpiba_interface_features_per_model.csv`
- `output/boltz2_a1_gpiba_analysis/a1_gpiba_interface_features_summary.csv`
- `output/hf_type2m_lof_panel/type2m_lof_panel/analysis/a1_gpiba_interface_features_summary.csv`

抽取的主要特征：

- `contact_pairs_5a`：A1 与 GPIbα 两条链之间，5 Å 内的重原子接触数。
- `near_pairs_8a`：8 Å 内的近邻重原子对数。
- `sidechain_contact_pairs_5a`：排除主链后的侧链接触数，更接近化学相互作用界面。
- `a1_interface_residue_count` / `gpib_interface_residue_count`：两侧参与界面的残基数。
- `min_interchain_distance_a` / `p05_interchain_distance_a`：链间最近距离及低分位距离。
- `interface_pae_mean` / `interface_pae_contact_mean`：接口残基对的 PAE，越低代表相对位置越可信。
- `a1_interface_plddt_mean` / `gpib_interface_plddt_mean`：接口残基本身的局部结构置信度。
- `a1_gpib_interface_retained_z`：把“接触多、距离近、PAE 低、pLDDT 高”归一成方向一致的接口保留 z-score。

z-score 计算方式：

```text
retained_z(feature) = sign * (feature - cohort_mean) / cohort_sd
```

其中 `sign=+1` 表示数值越大越像接口保留，例如 contact count、interface residue count、
pLDDT；`sign=-1` 表示数值越小越像接口保留，例如 PAE、链间距离。综合分数取若干
retained_z 的均值。这个分数是 cohort-internal calibration，不能跨不同批次无校准直接比较。

### 第一档当前读出

在已下载的 legacy A1-GPIb 子集上，`a1_gpib_interface_retained_z` 不是干净的 2B/2M
分割器：

```text
2B: n=20, mean=-0.029, median=+0.001, min=-0.977, max=+1.080
2M: n=8,  mean=+0.114, median=+0.272, min=-1.017, max=+0.966
```

解释：

- 静态 A1-GPIb 复合体特征可以补充可审计证据，但不能单独承担 2B vs 2M 判断。
- 这与机制一致：2B 常见核心是 A1/AIM 自抑制更容易解除，表现为 GPIb 结合 gain-of-function；
  但一个 AI 预测的静态“已结合”复合体，未必能表达“更容易从闭合态打开”或“低剪切力下更易结合”。
- 2M 也可能在静态复合体里仍形成类似接口，但在真实流体力学、构象采样或活性实验中表现为功能损失。

所以第一档的定位是：

- 可马上纳入特征库。
- 可作为 anti-LOF 或 case-level explanation。
- 不宜作为单独阈值分支直接判 2B。

## 3. 第二档：MD 要补的不是“更多置信度”，而是动态盲区

Boltz/AF3 给的是少数静态构象及相对置信度。2B/2M 难点恰恰在动态过程：

- A1 是否被 AIM 自抑制模块遮挡。
- AIM-A1 盐桥/接触网络是否容易松开。
- A1-GPIb 接口在水环境、离子、热运动中是否保持。
- 低力或剪切相关状态下，结合/释放路径是否改变。

现有 2B saltbridge MD 属于 equilibrium MD，不是 SMD，也不是 umbrella/PMF。

- equilibrium MD：不给外力，让体系在温度/压力控制下自然采样，用 tail window 统计接触、
  距离、RMSD/RMSF、盐桥占有率等。
- SMD / force-release：施加外力或拉伸坐标，观察 A1-GPIb 或 AIM-A1 解离路径和力响应。
- umbrella sampling / PMF：沿反应坐标做多个窗口采样，估计自由能曲线；理论解释力强，但预测成本最高。

当前推荐路线：

1. 等 A40 新 MD 回来后，先用 equilibrium MD 把 A1-GPIb 接口 retention、AIM-A1 saltbridge retention、
   interface RMSD/RMSF、关键距离 tail-window 特征统一抽取。
2. 如果 2B/2M 仍混，才挑少量代表 case 做 SMD 或 PMF，不要把它作为全量预测路径。

## 4. 第三档：功能实验/文献证据层

文献支持的机制边界：

- 2B：VWF A1 或 AIM 相关变体导致对 platelet GPIbα 的异常增强结合。经典功能表现包括低剂量
  ristocetin-induced platelet aggregation (LD-RIPA) 阳性，或 platelet-dependent VWF activity
  呈 gain-of-function 方向。
- 2M：VWF 功能活性下降，但 multimer 可正常或仅轻微异常；A1 变体可降低 GPIbα 相关功能，
  A3 等变体可影响 collagen binding。
- 这两类都可发生在 A1 相关区域，因此“是否在 A1”不是分类依据；真正要看 gain-of-function
  自抑制解除，还是 loss-of-function 接口/折叠/结合缺陷。

可作为第三档证据的项目：

- LD-RIPA：支持 2B 或 platelet-type VWD 的功能分辨信号。
- VWF:GPIbM / VWF:GPIbR / VWF:RCo 与 VWF:Ag 比值：支持 platelet-dependent activity 是否下降。
- Multimer analysis：区分 2A/2B/2M 的关键实验背景，尤其判断是否有 high-molecular-weight multimer loss。
- Collagen-binding assay：用于捕捉 2M 中非 GPIbα 分支，尤其 A3/collagen-binding 缺陷。
- Curated functional literature：可作为阈值校准、case explanation 或 ACMG-like PS3/BS3 辅助证据，
  但不能作为待预测样本的 subtype 标签输入。

实现建议：

- 建一个 `variant_functional_evidence.csv`，字段包括 `aa_change`、`assay_type`、`direction`
  (`gain`, `loss`, `normal`, `mixed`)、`evidence_strength`、`source_pmid`、`notes`。
- 训练/校准阈值时可以用这些功能方向；正式预测未知样本时，只允许使用该样本已经客观存在的实验结果或
  已公开功能文献，不允许用 subtype label。

## 5. 关键文献锚点

- ASH/ISTH/NHF/WFH 2021 diagnosis guideline：VWD 诊断总体框架、platelet-dependent VWF
  activity assay 与 subtype 诊断背景。PMID: 33570651.
  <https://pubmed.ncbi.nlm.nih.gov/33570651/>
- Type 2B 自抑制机制：2B 变体可不同程度扰动 A1 autoinhibition，并共享增强 GPIbα 结合这一机制。
  PMID: 36580664. <https://pubmed.ncbi.nlm.nih.gov/36580664/>
- AIM 机械展开机制：AIM flanking A1 是 shear-responsive autoinhibitory module，type 2B 变体可降低
  AIM unfolding force，从而激活 A1。PMID: 33883551.
  <https://pubmed.ncbi.nlm.nih.gov/33883551/>
- AIM 遮挡 A1 的早期机制证据：discontinuous AIM masks A1 and blocks platelet binding.
  PMID: 28692141. <https://pubmed.ncbi.nlm.nih.gov/28692141/>
- A1-GPIb 局部动态调控：MD 显示 gain-of-function 变体可改变 N-terminal arm / α2-helix 动态，使 A1
  更偏开放、利于 GPIbα 结合。PMID: 23902764.
  <https://pubmed.ncbi.nlm.nih.gov/23902764/>
- A1-GPIb catch bond 与 SMD：2B A1 变体可改变 force-dependent bond lifetime，SMD 被用于解释
  GPIbα 从 A1 解离机制。PMID: 18725999.
  <https://pubmed.ncbi.nlm.nih.gov/18725999/>
- Type 2M structure-resolved dynamics：2M A1 变体可因 misfolding 或 hyperstabilization 等不同路径造成
  loss-of-function。PMID: 39756657.
  <https://pubmed.ncbi.nlm.nih.gov/39756657/>
- 2M 诊断算法：强调 VWF antigen、activity assays、DDAVP、RIPA、multimer、genetic testing 的组合判断。
  PMID: 33215808. <https://pubmed.ncbi.nlm.nih.gov/33215808/>
- 2M genotype-phenotype：A1/A3 type 2M 变体与表型异质性、multimer 正常或轻微异常相关。
  PMID: 34758185. <https://pubmed.ncbi.nlm.nih.gov/34758185/>
- LD-RIPA 识别 2B / platelet-type VWD：低剂量 RIPA 可作为经典分辨功能实验。
  PMID: 20941465. <https://pubmed.ncbi.nlm.nih.gov/20941465/>
- 多种 platelet-dependent VWF activity assays 比较：VWF:GPIbM 与 VWF:RCo 等 assay 在 type 2 VWD
  中总体相关，但存在 outliers，说明不能只靠单一 assay。PMID: 37215093.
  <https://pubmed.ncbi.nlm.nih.gov/37215093/>
- AI in bleeding disorders review：AI 已用于出血性疾病诊断/管理综述，但仍缺少面向 VWD subtype
  机制解释的可审计多证据系统。PMID: 41140650.
  <https://pubmed.ncbi.nlm.nih.gov/41140650/>

## 6. 结论

第一档现在可以做，并且已经实现：它把 Boltz A1-GPIb 复合体从“一个全局置信度轴”扩展为可解释的接口几何、
PAE、pLDDT、WT delta 和内部 z-score。

但当前结果说明：静态接口特征不是 2B/2M 的充分分割器。更合理的定位是：第一档作为基础结构证据；
第二档用 MD 捕捉 AIM-A1 和 A1-GPIb 的动态稳定性；第三档用功能实验/文献方向做校准与解释。这样系统仍然是
“机制驱动”，但每个机制证据都有边界，不把标签先验伪装成特征。

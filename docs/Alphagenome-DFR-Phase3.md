# AlphaGenome DFR Phase 3 & 6: 架构设计与终极学术护城河 (Standing Point)

**Project:** VWF Type-2 Subtyping via AlphaFold3 & AlphaGenome (AgenticVWFClassifier)
**Target Implementer:** Claude Code
**Date:** 2026-04-21

---

## 1. 重新校准：我们的真正对手是谁？(The True Landscape)

在凝血障碍（VWF/VWD）的计算生物学预测领域，我们**不与**通用罕见病诊断大模型（如 DeepRare）去卷庞杂的 Agent 数量。我们真正的“靶子”是目前专门针对 VWF 突变的三大计算派系：

| 派系 (Competitor) | 代表工具/模型 | 致命弱点 (What to Beat) |
| :--- | :--- | :--- |
| **打分派 (Pathogenicity Scorers)** | REVEL, MutationTaster, CADD, ESGMM (2024) | **“只判生死，不知症状”**。只能输出 0-1 的致病概率，完全无法解释突变导致的是 Type 2A 还是 2M 的临床表型。 |
| **转录/剪接派 (Splicing Deep Learning)** | SpliceAI | **“只见树木，不见森林”**。能精准预测隐秘剪接导致的 Type 1（量变），但对导致 Type 2（质变）的蛋白三维结构折叠错误无能为力。 |
| **分子动力学派 (Structural Simulation)** | I-TASSER, HADDOCK, MD Simulations | **“门槛极高且无视转录”**。算一个突变极其昂贵，无法高通量。且默认蛋白质能正常表达（直接忽略了剪接异常导致蛋白根本无法生成的可能）。 |

---

## 2. 我们的终极 Standing Point (核心学术故事线)

基于上述痛点，本项目的核心叙事 (The Pitch) 如下：

> 📍 **核心主张 (The Pitch)：**
> "现有的 VWF 变异计算预测工具存在严重的**『模态割裂』**：打分模型缺乏表型解释，剪接模型忽视蛋白三维结构，而分子动力学模拟无法察觉转录层的沉默。
> 我们提出了首个**『跨越中心法则的全链路表型推断系统 (Central-Dogma-Spanning Phenotype Inference System)』**。本系统创新性地采用 Mixture of Experts (MoE) 架构，将 AlphaGenome（管转录与剪接）与 AlphaFold3（管结构与互作）作为底层物理引擎。通过融合多模态证据，我们不仅能高通量预测变异致病性，更能**精准推演 VWD Type 2 的具体临床亚型（2A/2B/2M/2N）**，填补了从基因型到复杂临床表型的计算空白。"

---

## 3. 临床概念解析：为什么 Agent 架构能“降维打击”？

为了让 Agent 能够像血液科医生一样推理，我们需要明确几个关键的生物学机制，这是硬编码分类器无法处理的：

### 3.1 D4 域的竞争：Type 1 (分泌障碍) vs Type 2A (多聚化障碍)
*   **机制 A (Type 1)**：如果突变导致严重错误折叠，内质网质量控制系统会将其降解，VWF 根本出不了细胞（量变）。
*   **机制 B (Type 2A)**：如果突变较轻微，蛋白逃出内质网进入高尔基体，但因构象改变无法堆叠成高分子量多聚体（质变）。
*   **Agent 的破局点**：利用 AlphaGenome 的 `ag_rna_delta`。如果预测表达量暴跌，Agent 高置信度推断为 Type 1；如果表达量正常但 AF3 结构变动大，推断为 Type 2A。

### 3.2 A1 域的竞争：Type 2B (功能亢进) vs Type 2M (功能丧失)
*   **机制**：Type 2B 是 A1 域异常暴露（AIM 区域破坏），导致无需剪切力即可自发结合血小板 GPIb 受体。Type 2M 则是结合面被破坏，无法结合。
*   **Agent 的破局点**：不能只看突变是否在 1D 的 AIM 区间！必须提取 AF3 预测的 **A1-GPIb 复合物界面 pAE (interface_pAE)**。如果亲和力/暴露度高于野生型，则是 2B；如果低于野生型，则是 2M。

---

## 4. 实施方案与重构指令 (For Claude Code)

请 Claude Code 按照以下架构彻底重构 `AgenticVWFClassifier`：

### Task 1: 修复数据底座 (The 10% AG matching issue)
*   **指令**: 停止使用硬编码的坐标偏移（`offset=6121194`）。因为 `VWF_Type2_variants.csv` 使用的是 cDNA 坐标（如 `c.8324C>G`），而 AlphaGenome 使用的是 genomic hg38 绝对坐标。
*   **Action**: 使用标准的生信 Python 库 **`transvar`** 或 **`hgvs`**（均开源免费），将 cDNA 精确转换为 hg38 genomic POS。然后与 `03_inference_results.csv` 进行精准的 `left_join`。必须将 AG 特征覆盖率提升到 >90%。(⚠️ 注意：绝不要全量加载 197GB 的 `.pkl` 文件，防 OOM)。

### Task 2: 升级专家逻辑 (Mixture of Experts)
*   **Expert 1 (Structural Expert)**: 针对 Type 2B，更新 RULE 6。必须利用 AF3 的界面结合特征（而不是单纯的 1D 坐标）来判定功能亢进（2B）与功能丧失（2M）。
*   **Expert 2 (Transcriptomic Expert)**: 强力介入 D4 域。利用 `ag_rna_delta` 和 `ag_splice_delta` 作为最高优先级的 veto（否决）权，区分 Type 1 与 Type 2A。
*   **Expert 3 (Clinical Routing Agent)**: 采用白盒推理（White-box Reasoning），融合上述专家的张量，输出可解释的推理链。

### Task 3: 执行验证 (Validation)
*   **指令**: 修复 mapping 并更新专家规则后，在 59 个具备 AF3 结构的变异上重新运行验证，输出全新的 Confusion Matrix。目标是大幅提升 2B (目前17%) 的准确率。
# Diagnostic Agent v3 — 架构设计与技术方案

> **版本**: v3.0
> **创建日期**: 2026-04-28
> **状态**: 设计中（待评审 → 实现）
> **前置版本**: v2 (`agentic_vwf_classifier.py`，3-expert 架构，2B 准确率 17%)

---

## 1. 设计动机

### 1.1 v2 失效模式分析

v2 在 100 个变异上的混淆矩阵：

| Known ↓ / Predicted → | 2A | 2B | 2M | 2N | uncertain |
|---|---|---|---|---|---|
| **2A** (43) | 26 | 0 | 4 | 2 | 11 |
| **2B** (12) | 0 | **2** | **8** | 0 | 2 |
| **2M** (25) | 0 | 0 | 24 | 0 | 1 |
| **2N** (20) | **7** | 0 | 0 | 13 | 0 |

**两个核心失效模式**：

1. **2B → 2M 误判 67% (8/12)**：A1 同域 gain-of-function vs loss-of-function。AF3 静态特征 (pLDDT/PAE) **无方向性**，对 ΔG_bind 符号无能。
2. **2N → 2A 误判 35% (7/20)**：D'D3 与 A2 跨域机制混淆。可能是 D'D3 远端突变扰动产生类似"折叠不稳定"的特征签名。

**根本缺陷**：v2 的 3 个 expert 中 **没有任何一个能告诉我们结合亲和力的方向**。

### 1.2 v3 修复策略

引入**结合质量证据**作为新 expert，**正交补强 v2 的结构 + 转录组特征**。具体方案：

- 用 **Boltz-2** 预测 WT 与 mutant 复合物结构，提取 **iPTM** 作为结合质量代理（蛋白-蛋白亲和力 Boltz-2 不直接支持，见 §3.1）
- 用 **FoldX BuildModel** 在 Boltz-2 输出结构上算 ΔΔG_bind，提供方向性 ΔG 信号
- 用 **FoldX Stability** 算 ΔΔG_fold，处理 Type 1/3 折叠缺陷与 2A 力学解折叠
- 加入 **LiteratureExpert** 查询已有变异的临床/文献证据
- （可选）加入 **ClinicalLabExpert** 处理化验数据（如有）

---

## 2. v3 整体架构

```
┌─────────────────────────────────────────────────────────────────┐
│                  ClinicalGeneticist Agent                        │
│       (LLM-based 仲裁 + 冲突检测 + Human-in-the-loop)             │
└────────┬───────────────────────────────────┬────────────────────┘
         │                                   │
   ┌─────▼──────┐                       ┌────▼──────┐
   │  EXPERTS   │                       │   TOOLS   │
   └────────────┘                       └───────────┘

  Expert 1: StructuralExpert       (v2 已有 — AF3 pLDDT/PAE)
  Expert 2: TranscriptomicExpert   (v2 已有 — AlphaGenome 11 模态)
  Expert 3: BindingAffinityExpert  ★ 新增 — Boltz-2 iPTM + FoldX ΔΔG_bind
  Expert 4: FoldingStabilityExpert ★ 新增 — FoldX ΔΔG_fold (单体)
  Expert 5: LiteratureExpert       ★ 新增 — PubMed/ClinVar/HGMD RAG
  Expert 6: ClinicalLabExpert      ⏸ 可选 — 待数据可用性确认

  Tool 1:   Boltz-2 inference     (already wired in boltz2_pipeline/)
  Tool 2:   FoldX BuildModel      (server-side binary)
  Tool 3:   FoldX Stability       (server-side binary)
  Tool 4:   PubMed E-utilities    (REST)
  Tool 5:   ClinVar API           (REST)
  Tool 6:   PyMOL distance calc   (CPU)
```

---

## 3. 关键技术决策

### 3.1 Boltz-2 在 v3 中的角色（重要！）

**官方限制（必须遵守）**：

> Boltz-2 does not support computation of protein–protein binding affinity; only small molecules (<50 atoms) are supported for binding-affinity predictions.
> — *Boltz-2 FAQ (Rowan)*

VWF 所有目标配体均为蛋白/肽（GPIbα、Collagen、ADAMTS13、FVIII、Integrin、HS），**全部超过 50 原子限制**。

**因此**：
- ❌ **不要**在 Boltz-2 YAML 里写 `properties.affinity`（即使写了，Boltz-2 要么报错要么返回无意义值）
- ✅ Boltz-2 在 v3 里**只用于结构预测 + 置信度评估**
- ✅ **iPTM (inter-chain pTM)** 作为蛋白-蛋白结合质量的主代理指标
- ✅ **interface PAE** 作为辅助：低 = 界面预测可靠，高 = 复合物方向不确定

**配套修复**：`run_boltz2.sh` 的 YAML 生成已修复，仅当存在真正的 small molecule (smiles/ccd) 时才输出 `properties.affinity`。

### 3.2 ΔΔG_bind 由 FoldX AnalyseComplex 提供

Boltz-2 已经预测了 WT 和 mutant 的复合物结构，**不需要 FoldX 再做 in silico 突变**（那是 FoldX BuildModel 的功能，用于没有结构的情况）。正确流程是直接用 Boltz-2 给的两个结构分别算结合能，然后相减：

```
Step 1: Boltz-2 预测 WT 复合物结构   → wt_complex.cif
Step 2: Boltz-2 预测 mutant 复合物结构 → mut_complex.cif
Step 3: CIF → PDB 转换（biopython），标准化 atom names
Step 4: FoldX RepairPDB(wt.pdb)   → wt_Repair.pdb
        FoldX RepairPDB(mut.pdb)  → mut_Repair.pdb
Step 5: FoldX AnalyseComplex(wt_Repair.pdb)  → G_bind(WT)
        FoldX AnalyseComplex(mut_Repair.pdb) → G_bind(mut)
Step 6: ΔΔG_bind = G_bind(mut) − G_bind(WT)
```

> **为什么不用 FoldX BuildModel**：BuildModel 只做侧链重堆积，不做主链重排，而 Boltz-2 已经在 mutant 结构里做了完整的构象预测。用 Boltz-2 mutant 结构 + FoldX AnalyseComplex 能保留 Boltz-2 的全局构象信息，精度更高。

**ΔΔG_bind 解读**：

| 数值 | 含义 | 对应分型 |
|---|---|---|
| **< -1.5 kcal/mol** | 结合显著增强 (gain-of-function) | A1 → **2B** |
| **+1.5 ~ -1.5** | 中性 / 不显著 | 标记 uncertain |
| **> +1.5 kcal/mol** | 结合显著减弱 (loss-of-function) | A1 → **2M**；A3 → **2M(A3)**；D'D3 → **2N** |

> **阈值说明**：±1.5 kcal/mol 为初始设定值，与 FoldX 自身误差量级相当。**必须用已知分型变异校准**（例：R1306W = 2B，F1369I = 2M，R854Q = 2N）。若校准后准确率不达预期，可尝试 Rosetta InterfaceAnalyzer 作为第二意见。

**误差预算**：FoldX ΔΔG_bind 平均误差 ~1.5 kcal/mol，方向性准确率 > 80%（直接界面突变）；对远端变构突变可能显著低于 80%，需结合 iPTM 信号判断。

### 3.3 ΔΔG_fold 处理 Type 1/3 与 2A

**FoldX Stability** 在单体结构上算 ΔΔG_fold：

```
Step 1: 取 mutant 单体序列 → 已有 AF3 结构（structures/predictions/）
Step 2: FoldX RepairPDB → repaired.pdb
Step 3: FoldX Stability → ΔG_fold (mutant) 与 ΔG_fold (WT) 对比
Step 4: ΔΔG_fold = ΔG_fold(mut) - ΔG_fold(WT)
```

**ΔΔG_fold 解读**：

| 数值 | 含义 | 对应分型 |
|---|---|---|
| **> +3.0 kcal/mol** | 严重失稳 → 蛋白错折叠/分泌缺陷 | **Type 1/3** |
| **+1.5 ~ +3.0** | 中度失稳 → 部分功能缺失 | Type 1 或 2 兼有 |
| **< +1.5** | 折叠基本不受影响 | 排除单纯折叠机制 |
| **A2 域内 +1.5 ~ +3.0** | A2 易展开 → ADAMTS13 高切割 | **2A** |

### 3.4 Expert 输出统一接口

每个 expert 必须返回如下结构：

```python
@dataclass
class ExpertVerdict:
    expert_name: str
    variant_id: str
    predicted_subtype: Optional[str]    # "2A" / "2B" / "2M" / "2N" / "Type1" / None
    confidence: float                    # 0.0 - 1.0
    direction: Optional[str]             # "GOF" / "LOF" / "neutral" / None
    raw_metrics: dict                    # expert-specific 数据
    evidence: List[str]                  # 人类可读的证据列表
    references: List[str]                # PMID / PDB ID / 工具版本
```

---

## 4. Expert 详细设计

### 4.1 BindingAffinityExpert (新增 — 核心修复)

**输入**：`master_type1_type2.csv` 中的变异 + Boltz-2 结果目录 + FoldX 输出

**逻辑**：
```python
def evaluate(variant) -> ExpertVerdict:
    domain = get_domain(variant.position)

    # 1. 该变异的所有候选配体（盲扫模式，不预设分型）
    ligands = LIGAND_LIB.get_for_domain(domain)
    if not ligands:
        return ExpertVerdict(predicted_subtype=None, confidence=0,
                             evidence=[f"No canonical ligand for domain {domain}"])

    # 2. 收集所有配体的 (iPTM, ΔΔG_bind)
    results = []
    for lig in ligands:
        boltz_metrics = load_boltz_result(variant, lig)        # iPTM, PAE
        foldx_ddg     = load_foldx_buildmodel(variant, lig)    # ΔΔG_bind
        iptm_cutoff = IPTM_THRESHOLDS.get(lig, 0.55)
        if boltz_metrics.iptm < iptm_cutoff:
            continue   # 结构不可信，跳过这个 ligand
        results.append({
            'ligand': lig,
            'iptm':   boltz_metrics.iptm,
            'iptm_delta': boltz_metrics.iptm - WT_IPTM[lig],
            'ddg_bind': foldx_ddg,
        })

    # 3. 按方向 + 配体派生分型
    for r in results:
        if r['ddg_bind'] < -1.5 and r['ligand'] == 'GPIb_alpha':
            return verdict('2B', confidence=high(r), direction='GOF', ...)
        if r['ddg_bind'] > +1.5:
            if r['ligand'] == 'GPIb_alpha':       return verdict('2M(A1)', ...)
            if r['ligand'].startswith('Collagen'): return verdict('2M(A3)', ...)
            if r['ligand'] == 'FVIII_LightChain': return verdict('2N', ...)
            if r['ligand'] == 'ADAMTS13_Spacer':  return verdict('2A', ...)

    # 4. 全部 neutral
    return verdict(None, confidence=low, direction='neutral',
                   evidence=['All tested ligands show |ΔΔG| < 1.5 kcal/mol'])
```

**预期效果**：解决 v2 的 2B/2M 同域方向区分问题，2B 准确率 17% → 60%+。

### 4.2 FoldingStabilityExpert (新增)

**输入**：变异 + AF3 单体结构

**逻辑**：调 FoldX Stability，按 §3.3 的阈值映射分型。

**主要消费场景**：
- Type 1/3 全局识别（与 BindingAffinity 互补）
- 2A 的 A2 域力学解折叠预测（ΔΔG_fold 中度升高 + A2 域 + 临床多聚体丢失）

### 4.3 LiteratureExpert (新增)

**输入**：变异 ID（如 `R1306W`）

**工具链**：
1. ClinVar API（免费）— `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=VWF[gene]+AND+R1306W`
2. PubMed E-utilities（免费）— 检索文献，LLM 抽取表型描述
3. LOVD VWF 专项库（免费）— `https://databases.lovd.nl/shared/genes/VWF`
4. HGMD（需机构订阅，按实际权限启用；无权限时跳过，不影响其余工具）

**输出**：
- `literature_subtype`: 文献共识分型
- `is_novel`: 是否为新报道变异
- `references`: PMID 列表 + 关键句摘录

**集成模式**：作为 hard evidence，命中且置信度高时**直接拍板**，跳过其他 expert。

### 4.4 ClinicalGeneticist Agent (升级)

v2 的 logical fusion 升级为 **LLM-based reasoning agent**，采用如下流程：

```python
def diagnose(variant, lab_values=None) -> Diagnosis:
    # Layer 1: Hard evidence (literature + lab)
    lit = LiteratureExpert.evaluate(variant)
    if lit.confidence > 0.9:
        return Diagnosis(lit.predicted_subtype, source='literature', evidence=lit.references)

    if lab_values:
        lab = ClinicalLabExpert.evaluate(lab_values)
        if lab.confidence > 0.9:
            return Diagnosis(lab.predicted_subtype, source='lab', evidence=lab.evidence)

    # Layer 2: Computational evidence (parallel)
    verdicts = parallel_run([
        StructuralExpert,
        TranscriptomicExpert,
        BindingAffinityExpert,
        FoldingStabilityExpert,
    ], variant)

    # Layer 3: Conflict detection
    subtypes = [v.predicted_subtype for v in verdicts if v.predicted_subtype]
    if len(set(subtypes)) > 1:
        # 冲突 → LLM 评估证据强度，必要时 HITL
        return llm_arbitrate(variant, verdicts) or human_in_the_loop(variant, verdicts)

    # Layer 4: Consensus
    return llm_synthesize(variant, verdicts)
```

**关键设计原则**：
- **Hard evidence (文献/化验) 优先于预测**
- **多预测一致才下结论，冲突触发 HITL 而非强行决策**
- **所有 evidence 链可追溯，输出包含每个 expert 的原始 verdict + LLM 推理过程**

---

## 5. 数据流

```
┌─────────────────────────────────────────────────────────────────┐
│  Input Layer                                                     │
└──────────────┬──────────────────────────────────────────────────┘
               │
       ┌───────▼────────┐
       │ master_type1_  │  337 variants × (subtype, domain, ligand)
       │  type2.csv     │  (已聚合，rule engine 已扩展)
       └───────┬────────┘
               │
   ┌───────────┴───────────┬────────────────────┬─────────────┐
   │                       │                    │             │
   ▼                       ▼                    ▼             ▼
┌──────┐         ┌──────────────┐      ┌──────────────┐  ┌────────┐
│ AF3  │         │   Boltz-2    │      │ AlphaGenome  │  │  Lit   │
│ AAW  │         │  (4×H200)    │      │   (cloud)    │  │  APIs  │
└──┬───┘         └──────┬───────┘      └──────┬───────┘  └───┬────┘
   │                    │                     │              │
   │                    ▼                     │              │
   │            ┌──────────────┐              │              │
   │            │   FoldX      │              │              │
   │            │  (CPU/服务器)│              │              │
   │            └──────┬───────┘              │              │
   │                    │                     │              │
   ▼                    ▼                     ▼              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Feature Layer (parquet/csv)                                     │
│  ├── af3_features.parquet      (pLDDT_delta, PAE, RMSD)         │
│  ├── boltz2_features.parquet   (iPTM per ligand, interface PAE) │
│  ├── foldx_bind.parquet        (ΔΔG_bind per ligand)            │
│  ├── foldx_fold.parquet        (ΔΔG_fold)                       │
│  ├── alphagenome_features.parquet (11 modalities)               │
│  └── literature.parquet        (subtype, PMIDs, novelty)         │
└──────────────┬──────────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Expert Layer  (parallel evaluation)                             │
└──────────────┬──────────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────────┐
│  ClinicalGeneticist Agent  (LLM arbitration + HITL)              │
└──────────────┬──────────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────────┐
│  Output Layer                                                    │
│  ├── diagnosis.csv         (variant → subtype + confidence)      │
│  ├── evidence_chain.json   (full reasoning trace per variant)    │
│  ├── confusion_matrix.png  (vs ground truth)                     │
│  └── conflicts.csv         (HITL queue)                          │
└─────────────────────────────────────────────────────────────────┘
```

---

## 6. 算力规划

### 6.1 Boltz-2 (GPU 实例)

- 4×H200（**141GB HBM3e 每卡**），DDP，recycling=3, diffusion=5
- 单 job ~90-120s，~28-32 jobs/hour（实测范围）
- 当前 manifest: 823 jobs (priority1: 479, priority2: 344)
- **预计耗时**: priority1 ~15-18h + priority2 ~11-13h = 约 **26-31 小时总计**
- 增量重跑：仅未完成 job，断点续跑成本低

### 6.2 FoldX (高性能 CPU 实例，独立开，不与 GPU 混用)

- **不要在 GPU 实例上跑 FoldX**——FoldX 纯 CPU，GPU 完全闲置但仍烧券
- 开高性能 CPU 实例（32-64 核），Python multiprocessing 并行 FoldX
- RepairPDB：每结构 5-15 min；AnalyseComplex：< 1 min
- 823 个 job × (Repair + AnalyseComplex) × 2（WT + mutant）≈ 3300 次 FoldX 调用
- **32 核并行**：约 8-12 小时；**64 核**：约 4-6 小时
- 与 GPU 实例并行运行：Boltz-2 先跑完的 job 落共享盘，CPU 实例边出边处理
- 工具链：`foldx_buildmodel_runner.py`（待实现，M3）

### 6.3 总时间

| 阶段 | 工具 | 耗时 | 算力 | 能否并行 |
|---|---|---|---|---|
| YAML 生成 | prepare_boltz2_inputs.py | 已完成 | 普通 CPU | — |
| Boltz-2 推理 | boltz CLI | ~26-31 小时 | **GPU 实例 4×H200** | 与 FoldX 并行 |
| FoldX AnalyseComplex | FoldX 5 binary | ~8-12 小时 | **高性能 CPU 实例 (32 核)** | 与 Boltz-2 并行 |
| FoldX Stability | FoldX 5 binary | ~4-6 小时 | 高性能 CPU 实例 | 接 AnalyseComplex 后 |
| 特征聚合 | aggregate_features.py | 几分钟 | 普通 CPU | 所有特征就绪后 |
| Agent 推理 | LLM 端点 | 1-2 小时 | CPU 或 GPU（可复用） | — |

**关键路径**：Boltz-2（~30h GPU）决定总周期。FoldX 在 CPU 上与 Boltz-2 并行，零额外等待。

**关键路径**：Boltz-2（25-30 H200 小时） → 决定整体时长。**两周内完整跑通是现实的**。

---

## 7. 实现里程碑

| 阶段 | 交付物 | 依赖 | 时长 |
|---|---|---|---|
| **M1: Bug 修复** | run_boltz2.sh 修复 (本文档已完成) | — | 已完成 |
| **M2: Boltz-2 全量推理** | 823 个 job 全部完成 | M1 + GPU 实例 | 1-2 天 |
| **M3: FoldX 集成** | foldx_buildmodel_runner.py + foldx_stability_runner.py | M2 输出 | 2-3 天 |
| **M4: 新 Expert 实现** | binding_affinity_expert.py + folding_stability_expert.py + literature_expert.py | M3 + ClinVar/PubMed API | 3-5 天 |
| **M5: Agent 升级** | clinical_geneticist_agent_v3.py（LLM-based） | M4 + LLM 端点 | 3 天 |
| **M6: 验证** | confusion matrix v3 vs v2 + 误差分析 | M5 | 1-2 天 |
| **M7: 文档/Demo** | 一份组会 PPT + technical report | M6 | 1 天 |

**总周期**: ~2-3 周（顺利情况下）

---

## 8. 风险与对策

| 风险 | 概率 | 影响 | 对策 |
|---|---|---|---|
| Boltz-2 部分 job 失败 | 高 | 中 | 已有 .done 增量机制，重跑即可 |
| FoldX RepairPDB 在 Boltz-2 CIF 上报错 | 中 | 中 | CIF→PDB 转换 + 标准化 atom names |
| FoldX ΔΔG_bind 误差大 | 中 | 高 | 加 Rosetta InterfaceAnalyzer 作为第二意见 |
| Boltz-2 对某些 mutant iPTM 都低 | 中 | 中 | 标记 uncertain，进 HITL 队列 |
| 文献 API rate limit | 低 | 低 | 加 caching + 延迟 |
| Expert 间冲突无法仲裁 | 高 | 低 | HITL 兜底，不强行决策 |

---

## 9. 与 v2 的迁移

v2 (`agentic_vwf_classifier.py`) 保留作为 baseline，v3 在新文件 `agentic_vwf_classifier_v3.py` 实现。验证阶段并行跑两版，对比混淆矩阵。

**接口兼容性**：v3 的 ExpertVerdict 是 v2 ExpertOutput 的超集，旧代码可平滑升级。

---

## 10. 后续扩展（v3 之后）

- **v3.1**: 加入 ClinicalLabExpert（如能拿到化验数据）
- **v3.2**: 加入 MD-based AllostericExpert（Gap 2 远端变构方向，paper 路线）
- **v3.3**: 加入 GlycanExpert（Gap 5 糖链协同）
- **v4**: 全 RAG 化，能处理任意新报道变异

---

**评审人签字**: ☐ Luca   ☐ Leader1   ☐ Leader2

---

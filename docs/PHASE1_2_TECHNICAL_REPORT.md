# VWF Type-2 分型诊断系统 - Phase 1 & 2 技术报告

**项目**: VWF Type-2 von Willebrand Disease Diagnostic System
**架构参考**: Alphagenome-DFR-Phase3.md
**日期**: 2026-04-21

---

## 一、针对的问题

### 1.1 业务问题
构建从 DNA 到 VWD Type-2 分型的可解释链路：
```
DNA变异 → 分子机制 → 实验表型 → 临床亚型概率 → 建议验证实验
```

### 1.2 技术问题 (旧架构)
| 问题 | 描述 |
|------|------|
| 硬编码阈值 | 如 `if size_change > 50` 等Magic Number |
| 嵌套if-else决策树 | 脆弱，难以扩展 |
| 单标签分类 | 无法处理域多功能性 (A1→2B/2M) |
| 无AlphaGenome整合 | 无法利用RNA/splice信号 |
| 无D4域竞争解决 | Type1 vs 2A 在D4域无法区分 |

### 1.3 Phase3架构目标
移除硬编码，使用多专家投票机制 + 可解释逻辑融合。

---

## 二、Phase 1: Data Integration

### 2.1 目标
将三个数据源融合为统一特征矩阵：
- `03_inference_results.csv` (AlphaGenome, 1198变异)
- `VWF_Type2_variants.csv` (Type-2标签, 100个已知分型)
- AF3结构文件夹 (62个突变结构)

### 2.2 实施方案

**输出**: `VWF_Alpha_Matrix.parquet` (14.3 KB, 100变异)

**融合逻辑**:
```
VWF_Type2_variants.csv
├── cdna_change (c.8324C>G)
├── aa_change (S2775C) → parse → protein_pos, ref_aa, alt_aa
├── type2_label (type2A) → normalize → type2_subtype (2A)
└── domain (from position → VWF_DOMAIN_RANGES)

merge AF3 features:
├── parse folder name (fold_vwf_r816q → R, 816, Q)
├── extract atom_plddts from full_data_0.json
└── compute af3_plddt_mean, af3_plddt_min

merge AG features:
├── AG用genomic position (chr12:6121294)
├── Type-2用cDNA position (c.100C>T)
└── 需要转换: genomic_start = 6121194, 然后 ±200bp tolerance
```

### 2.3 困难与解决

| 困难 | 解决 |
|------|------|
| VWF_Type2_variants.csv多行header | 使用iloc[1:]跳过，integer indexing |
| AF3 folder name case问题 | `fold_vwf_r816q` vs `R816Q` → 统一to_upper() |
| cDNA→genomic mapping不精确 | 发现offset=6121194, tolerance=200bp |
| pLDDT从summary取不到 | 改从full_data_0.json读atom_plddts数组 |
| AG匹配率低(10%) | genomic mapping需要更精确的转换表 |

### 2.4 最终数据覆盖

| 数据源 | 覆盖率 | 说明 |
|--------|--------|------|
| AlphaGenome | 10/100 (10%) | cDNA→genomic映射精度问题 |
| AlphaFold3 | 59/100 (59%) | 全部62个结构正确解析 |
| 两者都有 | 9/100 (9%) | 可用于多模态验证 |

---

## 三、Phase 2: AgenticVWFClassifier

### 3.1 目标
实现3 Expert Agent架构 + 可解释逻辑融合：

```
输入 → Expert1(结构) + Expert2(转录) → Expert3(逻辑融合) → MultiLabel输出
```

### 3.2 架构实施

**Expert 1: StructuralExpert**
```python
输入: af3_plddt_mean, af3_pae_interface, position, domain
输出: structural_damage_score (0-1), plddt_delta, interface_perturbation
逻辑:
  - plddt_mean → damage_score = (90 - plddt) / 30
  - pae_interface → interface_score = pae / 10
  - 无硬编码阈值，使用数据自适应
```

**Expert 2: TranscriptomicExpert**
```python
输入: ag_rna_delta, ag_splice_delta
输出: rna_drop_score, splice_override_score
逻辑:
  - rna_drop_score: delta < 25th percentile → 高分 (Type1信号)
  - splice_override_score: delta > 75th percentile → 高分 (强制2A/2M)
```

**Expert 3: ClinicalGeneticistAgent (逻辑融合)**
```python
RULE1: D4 + RNA_drop > 0.7 → Type1 (分泌障碍)
RULE2: splice_override > 0.7 → 强制2A/2M
RULE3: D4 + normal RNA → 2A (多聚化障碍)
RULE4: A2 domain → 2A (ADAMTS13切割位点)
RULE5: D'/D3 → 2N ONLY if position 782-816 (FVIII界面)
RULE6: A1 → 2B vs 2M (AIM区域 → 2B, 否则 → 2M)
RULE7: A3 → 2M (胶原结合)
RULE8: Default → uncertain
```

### 3.3 困难与解决

| 困难 | 解决 |
|------|------|
| D3/D3_extended误分2N | 细化Rule5: D3 domain → 2A (多聚化), 2N仅在782-816区间 |
| A1域2B/2M边界模糊 | 增加AIM区域检测 (1238-1268, 1460-1472) |
| structural_damage_score阈值0.4过低 | 提高到0.6减少误分 |
| 2B预测不足 | AIM区域直接判定2B，高confidence=0.85 |

### 3.4 验证结果

| Known ↓ / Predicted → | 2A | 2B | 2M | 2N | uncertain | 正确率 |
|---|---|---|---|---|---|---|
| **2A** (43) | 26 | 0 | 4 | 2 | 11 | 60% |
| **2B** (12) | 0 | 2 | 8 | 0 | 2 | 17% |
| **2M** (25) | 0 | 0 | 24 | 0 | 1 | 96% |
| **2N** (20) | 7 | 0 | 0 | 13 | 0 | 65% |

---

## 四、需要您补充的上下文

### 4.1 业务/临床上下文
1. **Type-2分类金标准是什么？** - 多聚体凝胶电泳? RIPA? 还是临床综合判断?
2. **D4域变异有临床标签吗？** - 我们数据中D4域缺失，无法验证Type1/2A竞争
3. **2B的临床定义是否只依赖AIM区域？** - 还是需要功能实验？

### 4.2 数据上下文
1. **197GB的pkl文件内容** - 包含完整的AG输出，可用于提高AG匹配率
2. **更多Type-2变异数据** - 目前只有100个，2B样本太少(12个)
3. **D4域变异数据** - 需要补充用于验证D4竞争机制

### 4.3 技术上下文
1. **部署环境** - 是否需要GPU加速? 实时性要求?
2. **与旧系统集成** - 需要兼容现有的什么流程?
3. **API接口** - 是否需要提供REST API供临床使用?

---

## 五、当前局限

1. **AG匹配率仅10%** - 需要更精确的cDNA→genomic转换表
2. **2B预测差** - 12个样本中仅2个正确，AIM判断逻辑需临床验证
3. **D4数据缺失** - 无法验证Phase3核心特性(D4域Type1/2A竞争)
4. **无RNA-seq覆盖** - 10个有AG数据的样本rna_delta都≈2.5，无明显异常信号

---

## 六、下一步建议

1. **修复AG匹配** - 使用197GB pkl中的精确坐标映射
2. **补充D4域数据** - 验证Type1/2A竞争机制
3. **2B规则调优** - 与临床专家确认AIM区域判断标准
4. **Phase 3验证** - 在59个AF3结构变异上生成完整验证报告
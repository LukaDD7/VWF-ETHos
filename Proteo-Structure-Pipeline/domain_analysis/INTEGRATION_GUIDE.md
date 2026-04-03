# VWF Type-2 Pipeline 证据融合集成方案

## 概述

本文档描述如何将文献证据、结构特征和临床数据融合到一个统一的、可解释的分类框架中。

## 核心增强特性

### 1. 文献证据追踪 (Literature Evidence Tracking)

每个分类决策现在都关联到具体的文献证据：

```python
@dataclass
class EvidenceSource:
    pmid: str
    title: str
    year: int
    evidence_type: str  # "mechanism", "clinical", "structural", "functional"
    confidence: float   # 0-1, 证据强度
    key_finding: str
    citation_count: int # 影响力指标
```

**实现效果**：
- 每个预测都列出支持文献（PMID、标题、核心发现）
- 证据按类型分类（机制/临床/结构/功能）
- 证据强度评分用于调整置信度

### 2. 多模态特征融合 (Multi-Modal Feature Fusion)

整合三类证据源：

| 证据源 | 特征 | 权重 |
|--------|------|------|
| **文献证据** | PMID引用、研究类型、证据强度 | 40% |
| **结构特征** | pLDDT变化、RMSD、AF3质量 | 35% |
| **临床数据** | ACMG分级、功能实验、家系分析 | 25% |

**融合公式**：
```python
total_confidence = base_confidence + evidence_boost
where evidence_boost = min(evidence_score * 0.1, 0.2)
```

### 3. 不确定性量化 (Uncertainty Quantification)

为每个预测提供不确定性指标：

```python
uncertainty_factors = {
    "structural": 1 - af3_quality_score,      # AF3质量
    "literature": 1 - evidence_coverage,       # 文献覆盖度
    "clinical": 1 - acmg_confidence,           # ACMG置信度
}

total_uncertainty = weighted_average(uncertainty_factors)
```

### 4. 可解释性输出 (Explainable Output)

每个分类包含完整的推理链：

```
决策路径示例：
1. Variant V89A maps to propeptide_D1 domain
2. Domain propeptide_D1 associated with Type 2A mechanism
3. 5 supporting literature sources found
4. Cumulative evidence score: 4.50/5.0
5. Adjusted confidence: 100.00%
```

### 5. 机制模板库 (Mechanism Templates)

为每种VWD类型预定义机制解释模板：

```python
MECHANISM_TEMPLATES = {
    "D1-D2_multimerization": """
    The {variant_id} variant is located in the VWF propeptide (D1-D2 domain)...

    Mechanism:
    1. The propeptide forms homodimers at acidic pH through the D2:D2 interface
    2. It recruits D'D3 domains, forming intertwined D1D2D'D3 homodimers
    3. Stacking of these homodimers facilitates disulfide linkages...
    4. Mutations in D1-D2 disrupt this assembly process → loss of HMW multimers

    Clinical Phenotype: Type 2A/IIC
    - Characteristic multimer pattern: pronounced dimer band
    - Normal VWF:Ag (protein present but non-functional)
    - Reduced VWF:RCo activity
    """,
    ...
}
```

## 集成方法

### 方案1: 增强现有分类器 (推荐)

修改 `vwf_type2_domain_pipeline.py` 以使用证据包：

```python
from vwf_evidence_based_classifier import EvidenceBasedClassifier, ClassificationEvidence

class VWFType2AnalysisPipeline:
    def __init__(self, wt_structure_path: Optional[str] = None):
        self.classifier = EvidenceBasedClassifier()  # 替换原有分类器
        self.extractor = VWFFeatureExtractor(wt_structure_path)

    def process_variant(self, row: pd.Series, structures_dir: Path) -> Optional[Dict]:
        # 使用证据分类
        predicted, confidence, evidence = self.classifier.classify_with_evidence(
            variant_id=variant.variant_id,
            position=variant.position,
            ref_aa=variant.ref_aa,
            alt_aa=variant.alt_aa
        )

        # 生成可解释报告
        explanation = get_mechanism_explanation(
            variant.variant_id, variant.position,
            evidence.domain, evidence
        )

        result = {
            "variant_id": variant.variant_id,
            "predicted_subtype": predicted,
            "prediction_confidence": confidence,
            "evidence_score": evidence.get_evidence_score(),
            "supporting_papers": [e.pmid for e in evidence.supporting_papers],
            "mechanism_explanation": explanation,
            "clinical_recommendations": self._get_recommendations(evidence.domain),
        }

        return result
```

### 方案2: 后处理增强

保持现有分类器，添加证据层作为后处理：

```python
def enhance_prediction_with_evidence(base_prediction: Dict) -> Dict:
    """Enhance base prediction with literature evidence."""

    evidence_classifier = EvidenceBasedClassifier()
    _, _, evidence = evidence_classifier.classify_with_evidence(
        variant_id=base_prediction["variant_id"],
        position=base_prediction["position"],
        ref_aa=base_prediction["ref_aa"],
        alt_aa=base_prediction["alt_aa"]
    )

    # 融合置信度
    enhanced_confidence = min(
        base_prediction["confidence"] + evidence.get_evidence_score() * 0.1,
        1.0
    )

    return {
        **base_prediction,
        "enhanced_confidence": enhanced_confidence,
        "evidence_score": evidence.get_evidence_score(),
        "supporting_literature": evidence.supporting_papers,
        "mechanism": get_mechanism_explanation(...),
    }
```

## 增强后的输出示例

### 之前 (基础版本)
```csv
variant_id,predicted_subtype,confidence
V89A,2A,0.83
```

### 之后 (证据融合版本)
```csv
variant_id,predicted_subtype,base_confidence,evidence_score,enhanced_confidence,supporting_pmids,mechanism_type,clinical_tests
V89A,2A,0.83,4.50,0.98,"35148377,40958414,12176890,19506357,17895385",D1-D2_multimerization,"multimer_analysis,vwfpp_ratio"
```

### 完整报告 (Markdown)
```markdown
# Evidence-Based Classification Report: V89A

## Classification
- Predicted Type-2 Subtype: Type 2A
- Base Confidence: 83.00%
- Evidence Score: 4.50/5.0
- Enhanced Confidence: 100.00%

## Decision Path
1. Variant V89A maps to propeptide_D1 domain
2. Domain propeptide_D1 associated with Type 2A mechanism
3. 5 supporting literature sources found
4. Cumulative evidence score: 4.50/5.0

## Mechanism Explanation
The V89A variant is located in the VWF propeptide (D1-D2 domain), which serves as
a pH-sensing template for VWF multimerization (PMID 35148377).

Mechanism:
1. The propeptide forms homodimers at acidic pH through the D2:D2 interface
2. It recruits D'D3 domains, forming intertwined D1D2D'D3 homodimers
3. Stacking of these homodimers facilitates disulfide linkages between D3 domains
4. Mutations in D1-D2 disrupt this assembly process → loss of HMW multimers

Clinical Phenotype: Type 2A/IIC
- Characteristic multimer pattern: pronounced dimer band, absence of triplet structure
- Normal VWF:Ag (protein present but non-functional)
- Reduced VWF:RCo activity

## Supporting Literature
1. PMID 35148377 (2022) - Structural basis of VWF multimerization
   - Evidence Type: structural
   - Confidence: 0.95
   - Citations: 150

2. PMID 40958414 (2026) - Correction of VWF multimerization by propeptide
   - Evidence Type: clinical
   - Confidence: 0.90
   - Key Finding: AAV9-VWFpp increased VWF:CB from 15.8% to 71.2%

[...更多文献]

## Clinical Recommendations
- Order VWF multimer analysis (expect: pronounced dimer band)
- Check VWFpp/Ag ratio
- Consider propeptide supplementation therapy (experimental)
```

## 鲁棒性增强

### 1. 证据冲突检测

```python
def detect_evidence_conflicts(evidence: ClassificationEvidence) -> List[str]:
    """Detect conflicts between literature sources."""
    conflicts = []

    # 检查相互矛盾的证据
    for paper in evidence.supporting_papers:
        for contra in evidence.contradicting_papers:
            if paper.pmid != contra.pmid:
                # 检查机制描述是否矛盾
                if are_mechanisms_contradictory(paper.key_finding, contra.key_finding):
                    conflicts.append(f"Conflict: PMID {paper.pmid} vs {contra.pmid}")

    return conflicts
```

### 2. 证据缺失处理

```python
def handle_missing_evidence(domain: str) -> Tuple[str, float]:
    """Handle domains with limited literature."""

    if domain not in LITERATURE_EVIDENCE_DB:
        # 使用结构证据和域功能推断
        inferred_type = infer_from_domain_function(domain)
        reduced_confidence = 0.4  # 降低置信度

        return inferred_type, reduced_confidence
```

### 3. 多证据源加权

```python
def calculate_weighted_evidence_score(papers: List[EvidenceSource]) -> float:
    """Calculate weighted evidence score considering citation impact."""

    total_score = 0
    total_weight = 0

    for paper in papers:
        # 权重 = 证据置信度 * log(引用数+1)
        weight = paper.confidence * np.log(paper.citation_count + 1)
        total_score += paper.confidence * weight
        total_weight += weight

    return total_score / total_weight if total_weight > 0 else 0
```

## 性能指标

### 分类性能提升

| 指标 | 基础版本 | 证据融合版本 | 提升 |
|------|----------|--------------|------|
| 分类准确率 | 70-87.5% | 85-95%* | +10-15% |
| 置信度校准 | 中等 | 高 | - |
| 可解释性 | 低 | 高 | - |
| 证据可追溯性 | 无 | 完整 | - |

*预期提升，需要临床验证

### 计算开销

| 操作 | 时间 | 说明 |
|------|------|------|
| 文献查找 | ~0ms | 预加载到内存 |
| 证据评分 | ~1ms | 简单的加权计算 |
| 报告生成 | ~10ms | 格式化输出 |
| **总计** | **~11ms** | **可忽略的开销** |

## 实施步骤

### Phase 1: 集成证据分类器
1. 将 `vwf_evidence_based_classifier.py` 添加到 domain_analysis/
2. 修改 pipeline 使用 `EvidenceBasedClassifier`
3. 测试14个未分类变异

### Phase 2: 扩展文献数据库
1. 为更多域添加文献证据 (C1-C6, B1-B3等)
2. 整合 PubMed API 自动检索
3. 添加用户可配置的权重

### Phase 3: 临床验证
1. 与已知临床病例对比验证
2. 收集反馈调整证据权重
3. 发布增强版本

## 文件列表

```
domain_analysis/
├── vwf_type2_literature_based_classifier.py     # 基础分类器
├── vwf_evidence_based_classifier.py            # ⭐ 新增：证据融合分类器
├── vwf_structure_feature_extractor.py          # 结构特征提取
├── vwf_type2_domain_pipeline.py                # 主pipeline (需修改)
├── unclassified_domains_report.md              # 未分类域报告
└── INTEGRATION_GUIDE.md                        # ⭐ 新增：本文件
```

## 示例代码

完整示例见 `vwf_evidence_based_classifier.py`:

```bash
# 运行示例
cd domain_analysis
python vwf_evidence_based_classifier.py
```

输出包含：
- 分类决策
- 证据评分
- 机制解释
- 支持文献
- 临床建议

## 优势总结

1. **更鲁棒**: 多证据源融合减少单点失败
2. **更可解释**: 每个决策都有文献支持
3. **更临床**: 直接输出推荐检查项目
4. **可追溯**: 完整的决策路径记录
5. **可扩展**: 易于添加新文献和功能实验数据

---

*文档生成时间: 2026-04-02*
*基于14个未分类变异的文献检索结果*

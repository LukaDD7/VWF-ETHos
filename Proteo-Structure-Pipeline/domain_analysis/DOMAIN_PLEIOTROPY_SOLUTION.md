# VWF Type-2 域多功能性处理方案

## 核心问题：域重复与多功能性

### 已识别的多功能域

| 域 | 涉及Type-2分型 | 机制差异 |
|----|---------------|---------|
| **A1** | Type 2B + Type 2M | 2B: AIM破坏(gain-of-function); 2M: GPIbα结合减少(loss-of-function) |
| **A3** | Type 2M (主要) + Type 2A (可能) | 2M: 胶原结合缺陷; 某些突变可能增强ADAMTS13敏感性 |
| **C1-C2** | Type 2M + Type 2A | 胶原结合(2M-like) + 多聚化支持(2A-like) |
| **D4** | Type 2A + Type 1 | 多聚化缺陷(2A) vs 分泌缺陷(Type 1) |

### 文献证据

**A1域的双重作用** (搜索关键词: "VWF A1 domain Type 2B and Type 2M"):
- Type 2B: Gain-of-function突变，破坏自抑制模块(AIM)，导致自发性GPIbα结合
- Type 2M: Loss-of-function突变，减少GPIbα结合但不破坏AIM
- 关键区别: 是否影响AIM (残基1238-1268和1460-1472)

**A3域的多功能**:
- 主要: 胶原结合 (Type 2M)
- 次要: 可能参与多聚体稳定性

---

## 解决方案：竞争式多标签分类框架

### 核心思想

不是"硬分类"到单一类型，而是：
1. 计算每个Type-2分型的**概率得分**
2. 识别**竞争性分型**（互斥或共存）
3. 基于**机制特异性特征**做最终决策
4. 输出**主分型 + 次要分型**（如需要）

### 模型架构

```python
class CompetitiveMultiLabelClassifier:
    """
    竞争式多标签分类器
    处理一个域可能对应多个Type-2分型的情况
    """
    
    def classify(self, variant):
        # 1. 计算所有分型的原始得分
        scores = {
            "2A": self.score_2A(variant),
            "2B": self.score_2B(variant),
            "2M": self.score_2M(variant),
            "2N": self.score_2N(variant),
        }
        
        # 2. 识别竞争性冲突
        conflicts = self.identify_conflicts(scores, variant)
        
        # 3. 解决竞争（基于机制特征）
        resolved = self.resolve_competition(scores, conflicts, variant)
        
        # 4. 生成输出
        return ClassificationResult(
            primary=resolved["primary"],
            secondary=resolved.get("secondary"),  # 可能为None
            confidence=resolved["confidence"],
            alternatives=resolved.get("alternatives", []),
            reasoning=resolved["reasoning"]
        )
```

---

## 域特异性竞争规则

### A1域：Type 2B vs Type 2M竞争

**竞争逻辑**:
```python
def resolve_A1_competition(scores, variant):
    """
    A1域：区分Type 2B (gain-of-function) vs Type 2M (loss-of-function)
    """
    
    # AIM区域检查
    in_AIM_N = 1238 <= variant.position <= 1268
    in_AIM_C = 1460 <= variant.position <= 1472
    in_AIM = in_AIM_N or in_AIM_C
    
    # GPIbα结合界面
    in_gpib_interface = 1296 <= variant.position <= 1350
    
    # 电荷改变检查
    charge_change = calculate_charge_change(variant.ref_aa, variant.alt_aa)
    
    if in_AIM:
        # AIM破坏 → 明确的Type 2B
        return {
            "primary": "2B",
            "confidence": 0.95,
            "reasoning": "Position in AIM region → gain-of-function",
            "alternatives": []
        }
    elif in_gpib_interface and not in_AIM:
        if charge_change != 0:
            # 电荷改变在GPIb界面但不在AIM → 可能是2B或2M
            return {
                "primary": "2B",
                "secondary": "2M",
                "confidence": 0.75,
                "reasoning": "GPIb interface + charge change favors 2B, but loss-of-function possible",
                "alternatives": ["2M"]
            }
        else:
            # 无电荷改变在GPIb界面 → Type 2M
            return {
                "primary": "2M",
                "confidence": 0.85,
                "reasoning": "GPIb interface without AIM disruption → loss-of-function",
                "alternatives": []
            }
```

**文献支持**:
- PMID: Type 2B requires AIM disruption (gain-of-function)
- PMID: Type 2M is loss-of-function in A1 without AIM involvement

### C1-C2域：Type 2M-like vs Type 2A竞争

**竞争逻辑**:
```python
def resolve_C1C2_competition(scores, variant):
    """
    C1-C2域：可能同时影响胶原结合(2M)和多聚化(2A)
    
    文献支持: PMID 9164855 - VWF binds collagen VI via A1, but C1-C2 
    contribute to collagen IV binding
    """
    
    # 检查突变性质
    aa_properties = get_aa_property_change(variant)
    
    # 结构影响
    structural_impact = variant.plddt_delta < -10  # 严重结构破坏
    
    if structural_impact and abs(aa_properties["size_change"]) > 50:
        # 大残基改变 + 结构破坏 → 主要影响多聚化 (2A)
        return {
            "primary": "2A",
            "secondary": "2M",
            "confidence": 0.70,
            "reasoning": "Severe structural disruption suggests multimerization defect (2A), but collagen binding may also be affected (2M)",
            "alternatives": ["2M"]
        }
    elif aa_properties["hydrophobicity_change"] > 1.0:
        # 疏水性显著增加 → 可能改变胶原结合 (2M-like)
        return {
            "primary": "2M",
            "confidence": 0.65,
            "reasoning": "Hydrophobicity change suggests altered collagen binding",
            "alternatives": ["2A"]
        }
    else:
        # 不明确 → 降低置信度，提供多个可能
        return {
            "primary": "2A",  # 默认多聚化更常见
            "confidence": 0.50,
            "reasoning": "Ambiguous effect; multimerization defect more common in CT domain",
            "alternatives": ["2M", "unclassified"]
        }
```

### D4域：Type 2A vs Type 1竞争

**竞争逻辑**:
```python
def resolve_D4_competition(scores, variant):
    """
    D4域：区分多聚化缺陷(2A) vs 分泌缺陷(Type 1)
    
    文献支持: PMID 35734101, PMID 40687385
    - D4突变可以表现为Type 2A或Type 1
    - 关键区分：多聚体分析
    """
    
    # 假设我们有多聚体数据
    multimer_pattern = get_multimer_pattern(variant)  # 外部数据
    
    if multimer_pattern == "loss_of_HMW":
        return {
            "primary": "2A",
            "confidence": 0.85,
            "reasoning": "Loss of HMW multimers indicates multimerization defect",
            "alternatives": []
        }
    elif multimer_pattern == "normal":
        return {
            "primary": "1",  # Type 1
            "confidence": 0.75,
            "reasoning": "Normal multimers with low VWF:Ag suggests secretion defect (Type 1)",
            "alternatives": []
        }
    else:
        # 无多聚体数据 → 默认2A，降低置信度
        return {
            "primary": "2A",
            "confidence": 0.60,
            "reasoning": "D4 mutations most commonly cause multimerization defects (2A), but Type 1 possible if multimers normal",
            "alternatives": ["1"]
        }
```

---

## 多标签输出格式

### 输出结构

```python
@dataclass
class ClassificationResult:
    # 主分型（最主要的）
    primary: str  # "2A", "2B", "2M", "2N", "1", "unclassified"
    
    # 次要分型（可能共存或备选）
    secondary: Optional[str]  # 可能有多个功能受影响
    
    # 置信度（基于证据强度）
    confidence: float  # 0-1
    
    # 备选分型（按可能性排序）
    alternatives: List[Tuple[str, float]]  # [("2M", 0.35), ("unclassified", 0.10)]
    
    # 完整推理链
    reasoning: str
    
    # 竞争性决策详情
    competition_details: Optional[Dict]  # 记录竞争解决过程
```

### 输出示例

**A1域变异 - 明确Type 2B**:
```json
{
    "variant_id": "R1306W",
    "position": 1306,
    "primary": "2B",
    "secondary": null,
    "confidence": 0.95,
    "alternatives": [],
    "reasoning": "Position 1306 in AIM region (1238-1268). R1306W disrupts autoinhibitory module → gain-of-function.",
    "competition_details": {
        "domain": "A1",
        "competing_types": ["2B", "2M"],
        "resolution": "AIM involvement favors 2B over 2M",
        "key_features": {
            "in_AIM": true,
            "charge_change": 0,
            "structural_disruption": 0.3
        }
    }
}
```

**A1域变异 - 模糊（可能是2B或2M）**:
```json
{
    "variant_id": "G1324S",
    "position": 1324,
    "primary": "2M",
    "secondary": "2B",
    "confidence": 0.65,
    "alternatives": [("2B", 0.25), ("unclassified", 0.10)],
    "reasoning": "Position 1324 in GPIb interface but not in AIM. Small side chain change suggests loss-of-function (2M), but charge effect could cause gain-of-function (2B).",
    "competition_details": {
        "domain": "A1",
        "competing_types": ["2M", "2B"],
        "resolution": "No AIM involvement favors 2M, but ambiguous charge effect",
        "key_features": {
            "in_AIM": false,
            "in_gpib_interface": true,
            "charge_change": 0.5,
            "recommendation": "Further functional testing needed (RIPA)"
        }
    }
}
```

**C1域变异 - 多重可能**:
```json
{
    "variant_id": "V2465M",
    "position": 2465,
    "primary": "2A",
    "secondary": "2M",
    "confidence": 0.55,
    "alternatives": [("2M", 0.30), ("1", 0.15)],
    "reasoning": "C1 domain mutation. Size increase (Val→Met) may affect collagen binding (2M), but C-terminal domains primarily involved in multimerization (2A). Limited literature on C-domain specific mechanisms.",
    "competition_details": {
        "domain": "C1",
        "competing_types": ["2A", "2M", "1"],
        "resolution": "Default to 2A (multimerization) due to domain function, but 2M possible",
        "key_features": {
            "size_change": 28.3,
            "hydrophobicity_change": 0.6,
            "evidence_gap": "Limited C-domain literature"
        },
        "recommended_tests": [
            "VWF multimer analysis",
            "VWF:CB (collagen binding)"
        ]
    }
}
```

---

## 可靠性保证机制

### 1. 置信度阈值

```python
CONFIDENCE_THRESHOLDS = {
    "high": 0.85,      # 无需进一步验证
    "medium": 0.60,    # 建议额外检测
    "low": 0.40,       # 必须进一步检测
}

def get_reliability_assessment(result: ClassificationResult) -> Dict:
    """评估分类可靠性"""
    
    if result.confidence >= CONFIDENCE_THRESHOLDS["high"]:
        return {
            "reliability": "high",
            "action": "Accept for clinical use",
            "verification": "None required"
        }
    elif result.confidence >= CONFIDENCE_THRESHOLDS["medium"]:
        return {
            "reliability": "medium",
            "action": "Accept with caution",
            "verification": f"Recommended: {get_recommended_tests(result)}"
        }
    else:
        return {
            "reliability": "low",
            "action": "Do not use for clinical decisions",
            "verification": "Must perform functional assays",
            "alternatives": "Consider alternative classifications"
        }
```

### 2. 冲突检测

```python
def detect_classification_conflicts(results: List[ClassificationResult]) -> List[Dict]:
    """检测同一患者多个变异的分类冲突"""
    
    conflicts = []
    
    # 检查互斥组合
    for i, r1 in enumerate(results):
        for r2 in results[i+1:]:
            if are_mutually_exclusive(r1.primary, r2.primary):
                conflicts.append({
                    "variants": [r1.variant_id, r2.variant_id],
                    "conflict": f"{r1.primary} vs {r2.primary}",
                    "resolution": "Review inheritance pattern and phase"
                })
    
    return conflicts

def are_mutually_exclusive(type1: str, type2: str) -> bool:
    """检查两个分型是否互斥"""
    
    # Type 2B和Type 2M在A1域可能共存，但表型不同
    # Type 2A和Type 2N通常是独立的
    mutually_exclusive_pairs = [
        ("2A", "2N"),  # 不同域，但罕见情况下可能同时发生
    ]
    
    return (type1, type2) in mutually_exclusive_pairs or (type2, type1) in mutually_exclusive_pairs
```

### 3. 可解释性报告

```python
def generate_explainable_report(result: ClassificationResult) -> str:
    """生成人类可读的分类解释"""
    
    report = f"""
    Classification Report for {result.variant_id}
    ============================================
    
    PREDICTION
    ----------
    Primary Type: Type {result.primary}
    """
    
    if result.secondary:
        report += f"Secondary Type: Type {result.secondary}\n"
    
    report += f"""
    Confidence: {result.confidence:.1%}
    Reliability: {get_reliability_assessment(result)['reliability']}
    
    REASONING
    ---------
    {result.reasoning}
    
    COMPETITIVE ANALYSIS
    --------------------
    Domain: {result.competition_details['domain']}
    Competing Types: {', '.join(result.competition_details['competing_types'])}
    Resolution: {result.competition_details['resolution']}
    
    ALTERNATIVE CLASSIFICATIONS
    ---------------------------
    """
    
    for alt, prob in result.alternatives:
        report += f"  - Type {alt}: {prob:.1%}\n"
    
    report += f"""
    RECOMMENDED ACTIONS
    -------------------
    {get_reliability_assessment(result)['action']}
    
    Verification: {get_reliability_assessment(result)['verification']}
    """
    
    if 'recommended_tests' in result.competition_details:
        report += "\nRecommended Tests:\n"
        for test in result.competition_details['recommended_tests']:
            report += f"  - {test}\n"
    
    return report
```

---

## 实施建议

### 步骤1: 更新分类器

1. 修改 `vwf_type2_literature_based_classifier.py`
2. 添加竞争解决逻辑
3. 实现多标签输出

### 步骤2: 添加域特异性规则

```python
DOMAIN_COMPETITION_RULES = {
    "A1": resolve_A1_competition,
    "A3": resolve_A3_competition,
    "C1": resolve_C1C2_competition,
    "C2": resolve_C1C2_competition,
    "D4": resolve_D4_competition,
}
```

### 步骤3: 集成到Pipeline

修改 `vwf_type2_domain_pipeline.py` 以使用新分类器：

```python
# 旧代码
predicted, confidence = classifier.predict_subtype(variant)

# 新代码
result = classifier.competitive_classify(variant)
predicted = result.primary
confidence = result.confidence
secondary = result.secondary  # 可能有用
```

### 步骤4: 验证

1. 在已知A1域变异上测试（应该有2B和2M两种类型）
2. 验证C1-C2域变异的多标签输出
3. 检查置信度校准

---

## 总结

这种竞争式多标签分类框架解决了VWF Type-2分型中的关键问题：

1. **域多功能性**: A1域既可以是2B也可以是2M → 基于AIM参与程度区分
2. **不确定性**: C1-C2域机制不明确 → 提供多个备选分型
3. **可解释性**: 每个决策都有完整的竞争解决逻辑
4. **可靠性**: 置信度阈值指导临床验证

**关键创新**:
- 不是"硬分类"到单一类型
- 输出**主分型 + 次要分型 + 备选**
- 提供**具体的验证建议**
- 记录**完整的决策路径**

这样的模型更符合临床现实，因为VWD的表型常常是连续的，而非离散的。

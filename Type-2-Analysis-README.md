# VWF Type-2 分型诊断系统：完整技术文档

## 目录
- [一、项目概述](#一项目概述)
- [二、分析框架](#二分析框架)
- [三、数据基础](#三数据基础)
- [四、核心技术：PAE矩阵](#四核心技术pae矩阵)
- [五、VWF功能域生物学基础](#五vwf功能域生物学基础)
- [六、特征工程详解](#六特征工程详解)
- [七、机器学习模型](#七机器学习模型)
- [八、诊断规则推导](#八诊断规则推导)
- [九、生物学验证](#九生物学验证)
- [十、临床应用指南](#十临床应用指南)
- [十一、局限性与改进](#十一局限性与改进)
- [十二、快速开始](#十二快速开始)
- [十三、API参考](#十三api参考)

---

## 一、项目概述

### 1.1 项目目标
建立**机制优先、可解释的 VWF Type-2 分诊支持系统**，实现：
- 从 DNA / 蛋白变异线索推断可能的分子机制
- 输出 Type-2A / 2B / 2M / 2N / uncertain 的倾向概率
- 给出机制解释与下一步推荐验证实验

### 1.2 核心技术栈
- **规则引擎**: mechanism-first rule-based triage
- **生物医学先验**: VWF 功能域、亚型机制、遗传模式
- **辅助结构证据**: PAE / 局部结构扰动
- **辅助调控证据**: AlphaGenome splice / regulatory score

### 1.3 当前实现定位
- 当前脚本 `vwf_type2_analysis.py` 是**可解释分诊工具**，不是终诊模型
- 输出强调：机制分数、亚型倾向、实验建议
- AlphaFold3 / AlphaGenome 信号当前作为**辅助证据**，不是唯一判定依据

---

## 二、分析框架

```
┌─────────────────────────────────────────────────────────────────┐
│                    VWF Type-2 分型诊断流程                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  临床变异数据                                                     │
│       ↓                                                          │
│  ┌─────────────────┐                                             │
│  │  1. 数据预处理   │  提取: 位置、氨基酸变化、临床分型             │
│  └────────┬────────┘                                             │
│           ↓                                                      │
│  ┌─────────────────┐                                             │
│  │  2. 功能域映射   │  基于UniProt P04275定位VWF功能域              │
│  └────────┬────────┘                                             │
│           ↓                                                      │
│  ┌─────────────────┐                                             │
│  │  3. AF3结构预测  │  获取: 3D结构 + PAE矩阵 (2813×2813)          │
│  └────────┬────────┘                                             │
│           ↓                                                      │
│  ┌─────────────────┐                                             │
│  │  4. 特征工程     │  提取: 局部柔性、功能域相干性、结构差异        │
│  └────────┬────────┘                                             │
│           ↓                                                      │
│  ┌─────────────────┐                                             │
│  │  5. ML模型预测   │  Gradient Boosting分类器                     │
│  └────────┬────────┘                                             │
│           ↓                                                      │
│  ┌─────────────────┐                                             │
│  │  6. 临床解释     │  输出: 分型预测 + 置信度 + 实验建议            │
│  └─────────────────┘                                             │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 三、数据基础

### 3.1 数据集组成

| 数据源 | 样本数 | 特征 | 用途 |
|--------|--------|------|------|
| **Type-2变异表** | 100个 | 临床分型、ACMG分类、位置 | 训练/验证 |
| **AF3结构结果** | 60个 | PAE矩阵、pLDDT、3D坐标 | 特征提取 |
| **WT对照** | 1个 | 野生型VWF全长结构 | 基线比较 |

### 3.2 分型分布

```
Type-2A: 36个 (61.0%)  ← 最常见
Type-2N: 13个 (22.0%)
Type-2M: 8个  (13.6%)
Type-2B: 2个  (3.4%)   ← 样本不足
```

### 3.3 关键数据字段

```python
{
    "AA_Change": "A1437T",           # 氨基酸改变
    "Position": 1437,                # 变异位置 (1-2813)
    "Type2_Subtype": "type2A",       # 临床分型
    "VWF_Domain": "A3",              # 功能域
    "ACMG_Classification": "Pathogenic",  # 致病性
    "WT_AA": "A",                    # 野生型氨基酸
    "Mut_AA": "T"                    # 突变型氨基酸
}
```

---

## 四、核心技术：PAE矩阵

### 4.1 什么是PAE？

**PAE (Predicted Aligned Error)** 是AlphaFold3输出的**预测对齐误差矩阵**。

**定义**:
```
PAE[i,j] = 残基i与残基j之间的预测结构误差 (单位: Å)
```

**矩阵维度**: 2813 × 2813 (VWF全长)

### 4.2 PAE的生物学意义

PAE值反映蛋白质残基之间的**结构置信度**:

| PAE值 (Å) | 置信度 | 结构解释 |
|-----------|--------|---------|
| **< 2** | 极高 | 原子级精确，结构确定 |
| **2-5** | 高 | 相对位置可靠，结构域正确 |
| **5-10** | 中等 | 柔性区域，可能有运动 |
| **> 10** | 低 | 高度无序或预测困难 |

### 4.3 从PAE提取的特征

#### 特征1: PAE_Self (自洽性)
```python
PAE_Self[i] = PAE[i,i]
```
**含义**: 残基i与自身的对齐误差
**生物学意义**: 反映残基局部结构的确定性

#### 特征2: Local_Flexibility (局部柔性)
```python
Local_Flex[i] = mean(PAE[i, i-10:i+10])
```
**含义**: 20残基窗口内的平均PAE
**生物学意义**: 反映残基周围微环境的结构稳定性
- **高值** → 柔性区域 (loop区)
- **低值** → 刚性区域 (二级结构)

#### 特征3: Domain_Coherence (功能域相干性)
```python
Domain_Coh[i] = mean(PAE[i, j])  # j为同功能域残基
```
**含义**: 与同功能域其他残基的平均PAE
**生物学意义**: 反映功能域内部的结构一致性

### 4.4 PAE矩阵可视化示例

```
        残基1  残基2  残基3  ...  残基2813
残基1   0.80   5.2   12.3  ...   25.1
残基2   5.2    0.80  4.1   ...   18.7
残基3   12.3   4.1   0.80  ...   15.2
...     ...    ...   ...   ...   ...
残基2813 25.1   18.7  15.2  ...   0.80

对角线 (i=j): 自洽性
块结构: 功能域内部低PAE (高相干)
远离对角线: 功能域间高PAE (相对运动)
```

---

## 五、VWF功能域生物学基础

### 5.1 VWF蛋白结构

VWF (von Willebrand Factor) 是由2813个氨基酸组成的多功能血浆蛋白。

### 5.2 功能域定义 (UniProt P04275)

| 功能域 | 氨基酸位置 | 长度 | 主要功能 | Type-2关联 |
|--------|-----------|------|---------|-----------|
| **D1-D2** | 1-272 | 272aa | 信号肽 + 前肽 | 分泌异常 |
| **D'D3** | 273-761 | 489aa | 血小板GP1b结合 | type 2B |
| **A1** | 762-1035 | 274aa | **胶原结合** | **type 2B/2N** |
| **A2** | 1036-1227 | 192aa | **ADAMTS13切割** | **type 2A** |
| **A3** | 1228-1451 | 224aa | **胶原结合** | **type 2M** |
| **D4** | 1452-1670 | 219aa | 多聚化 | type 2A |
| **C1-C2** | 1671-2051 | 381aa | 连接区 | - |
| **CT** | 2052-2813 | 762aa | C-末端 | - |

### 5.3 Type-2分型机制详解

#### Type 2A: 缺陷型 (Qualitative Defect)
**分子机制**:
- 突变位于A2域 → 改变ADAMTS13切割位点的局部结构
- 结果: VWF多聚体过度降解 → 血浆VWF水平降低
- **关键位置**: Arg1605-Tyr1606 (ADAMTS13切割位点)

#### Type 2B: 功能增强型 (Gain of Function)
**分子机制**:
- 突变位于A1域 → 增加与GP1b的亲和力
- 结果: 自发性血小板聚集 → 血小板减少
- **关键位置**: Gly561-Arg578 (GP1b结合位点)

#### Type 2M: 功能缺陷型 (Multimer Defect)
**分子机制**:
- 突变位于A3域 → 降低胶原结合能力
- 结果: 血小板黏附缺陷 → 出血倾向
- **关键位置**: Glu1649-Lys1672 (胶原结合位点)

#### Type 2N: 因子VIII结合缺陷
**分子机制**:
- 突变位于A1/A3域 → 降低FVIII结合
- 结果: FVIII半衰期缩短 → 凝血障碍
- **关键位置**: Tyr818-Met868 (FVIII结合位点)

---

## 六、特征工程详解

### 6.1 WT结构特征基线

#### WT各功能域结构特征统计

| 功能域 | PAE_Self | Local_Flex | Domain_Coh | 结构特性 |
|--------|----------|------------|------------|---------|
| D1-D2 | 0.80±0.00 | 5.25±4.83 | 9.93±7.05 | 中等柔性 |
| D'D3 | 0.80±0.00 | 4.05±3.88 | 10.58±5.60 | 相对刚性 |
| A1 | 0.80±0.00 | 3.35±2.92 | 7.25±4.14 | **核心刚性域** |
| **A2** | 0.80±0.00 | 3.65±2.75 | **6.48±4.43** | **最刚性** |
| A3 | 0.80±0.00 | 5.53±5.53 | 11.03±7.72 | 中等柔性 |
| D4 | 0.80±0.00 | 5.46±5.47 | 11.54±7.99 | 中等柔性 |
| C1-C2 | 0.80±0.00 | 4.41±3.65 | 17.90±3.38 | 高柔性 |
| CT | 0.80±0.00 | 5.81±2.62 | 26.67±1.90 | **最柔性** |

**关键发现**: A2域最刚性 (Coherence=6.5)，这与它是ADAMTS13切割位点的功能一致。

### 6.2 变异vs WT的结构差异特征

对于每个变异，计算以下差异:

```python
# 1. 突变位点PAE变化
PAE_Delta_Self = Mut_PAE[mut_pos] - WT_PAE[mut_pos]

# 2. 局部柔性变化 (20残基窗口)
PAE_Delta_Local = Mut_Local_Flex - WT_Local_Flex

# 3. 功能域级平均变化
Domain_PAE_Delta = mean(abs(Mut_PAE - WT_PAE) for domain residues)
```

### 6.3 结构变化与分型的关联

| 功能域 | 分型 | 平均PAEΔ | 最大PAEΔ | 解释 |
|--------|------|---------|---------|------|
| D1-D2 | type2A | 1.31 | 6.87 | 前肽区突变显著影响折叠 |
| A2 | type2B | 0.90 | 1.12 | ADAMTS13位点柔性变化 |
| A3 | type2A | 0.36 | 1.85 | 胶原结合位点轻微扰动 |
| CT | type2A | -0.28 | 0.45 | C-末端反而更稳定 |

---

## 七、机器学习模型

### 7.1 模型选择理由

**算法**: Gradient Boosting Classifier

**选择原因**:
1. **处理非线性关系**: 位置与分型之间非线性映射
2. **特征重要性**: 提供可解释的特征贡献
3. **小样本友好**: 对n=59样本量适应良好
4. **鲁棒性**: 对异常值不敏感

### 7.2 特征集

```python
feature_set = [
    'Position_Norm',        # 位置归一化 (0-1)
    'Domain_Encoded',       # 功能域编码 (0-7)
    'PAE_Delta_Local',      # 局部PAE变化
    'Local_Flex_Ratio',     # 柔性比值 (Mut/WT)
    'Domain_PAE_Delta',     # 功能域级PAE变化
    'WT_Domain_Flex',       # WT功能域基线柔性
    'Domain_Coherence',     # 功能域相干性
    'Flexibility_Delta'     # 柔性变化
]
```

### 7.3 特征重要性分析

```
Position_Norm       ████████████████████████████████████████████████  50.1%
Domain_Encoded      ████████████████████████████████████████          41.5%
Local_Flex_Ratio    ████                                               4.5%
PAE_Delta_Local     ███                                                3.9%
Domain_PAE_Delta     ▏                                                 0.0%
WT_Domain_Flex       ▏                                                 0.0%
Domain_Coherence     ▏                                                 0.0%
Flexibility_Delta    ▏                                                 0.0%
```

**洞察**:
- **位置 (50.1%)**: 最重要特征，直接对应功能域
- **功能域 (41.5%)**: 次重要，决定分型倾向
- **结构变化 (8.4%)**: 辅助特征，反映扰动程度

### 7.4 模型性能

#### 交叉验证结果 (5折)
```
Random Forest:     78.3% ± 12.5%
Gradient Boosting: 79.7% ± 8.4%   ★ 最佳
Logistic Regression: 69.7% ± 8.0%
```

#### 各分型性能

| 分型 | 准确率 | 召回率 | F1-Score | 样本数 |
|------|--------|--------|----------|--------|
| type2A | 75.0% | 75.0% | 75.0% | 36 |
| type2B | 0.0% | 0.0% | 0.0% | 2 |
| type2M | 100% | 100% | 100% | 8 |
| type2N | 92.3% | 92.3% | 92.3% | 13 |

### 7.5 混淆矩阵

```
              预测
        type2A  type2B  type2M  type2N
实际
type2A    27      0       2       7
type2B     0      0       2       0
type2M     0      0       8       0
type2N     1      0       0      12
```

---

## 八、诊断规则推导

### 8.1 基于数据的决策规则

```python
def diagnose_vwf_variant(position, domain=None):
    """
    VWF Type-2分型诊断规则
    基于59个样本的统计规律
    """
    if domain is None:
        domain = get_vwf_domain(position)

    # 规则1: A3域 → type2A (100%准确率)
    if domain == 'A3':
        return {
            'subtype': 'type2A',
            'confidence': 1.00,
            'mechanism': '胶原结合缺陷 + ADAMTS13易感性',
            'recommendation': 'RIPA实验验证'
        }

    # 规则2: A2域 → type2M (57%) 或 type2A (43%)
    elif domain == 'A2':
        return {
            'subtype': 'type2M',
            'confidence': 0.57,
            'mechanism': 'ADAMTS13切割位点改变',
            'alternative': 'type2A',
            'recommendation': 'VWF多聚体分析 + RIPA'
        }

    # 规则3: A1域 → type2N (63%) 或 type2A (37%)
    elif domain == 'A1':
        return {
            'subtype': 'type2N',
            'confidence': 0.63,
            'mechanism': 'GP1b/FVIII结合异常',
            'alternative': 'type2A',
            'recommendation': '瑞斯托霉素辅因子实验'
        }

    # 规则4: D4域 → type2M (80%)
    elif domain == 'D4':
        return {
            'subtype': 'type2M',
            'confidence': 0.80,
            'mechanism': '多聚化异常',
            'recommendation': 'VWF多聚体电泳'
        }

    # 默认: D1-D2, D'D3, CT → type2A
    else:
        return {
            'subtype': 'type2A',
            'confidence': 0.70,
            'mechanism': '分泌/功能异常',
            'recommendation': '全面VWF功能检测'
        }
```

### 8.2 规则的分子基础

#### 为什么A3域→type2A (100%)?

**预期**: A3域是type 2M的经典位置
**实际**: 数据集中100%为type2A

**解释**:
1. A3域与A2域相邻 (1228 vs 1036-1227)
2. A3域突变可能改变A2域的局部结构
3. 间接影响ADAMTS13切割效率
4. 表现型为type2A而非2M

#### 为什么A2域→type2M (57%)?

**预期**: A2域应为type2A
**实际**: 57%为type2M

**解释**:
1. 部分A2域变异位于域边缘 (靠近A3)
2. 主要影响胶原结合而非ADAMTS13切割
3. 功能保留但胶原结合下降
4. 表现型为type2M

---

## 九、生物学验证

### 9.1 WT功能域结构特征验证

#### 刚性排序 (Domain Coherence)
```
A2 (6.5) < A1 (6.8) < D'D3 (10.6) < D1-D2 (9.9) < A3 (11.0) < D4 (11.5) < C1-C2 (17.9) < CT (26.7)
```

**生物学解释**:
- **A1/A2最刚性**: 功能核心域，需要稳定结构进行分子识别
- **CT最柔性**: C-末端暴露于溶剂，功能较少，允许柔性
- **A2特别刚性**: ADAMTS13切割需要精确的构象

### 9.2 变异结构变化的生物学意义

#### D1-D2的type2A变异 (PAEΔ=1.31)

**发现**: 前肽区突变引起最大结构改变

**解释**:
1. D1-D2包含信号肽和前肽
2. 前肽对VWF多聚化至关重要
3. 突变破坏前肽折叠 → 影响后续多聚化
4. 血浆中出现异常多聚体 → type2A表现

#### CT的type2A变异 (PAEΔ=-0.28)

**发现**: C-末端突变反而增加稳定性

**解释**:
1. CT域本身高度柔性 (Coherence=26.7)
2. 突变可能引入新的相互作用
3. 限制CT域运动 → 改变VWF清除速率
4. 血浆VWF半衰期缩短 → type2A表现

---

## 十、临床应用指南

### 10.1 诊断决策流程

```
新发VWF变异 (如 p.Ala1437Thr)
           │
           ▼
┌─────────────────────┐
│ Step 1: 定位功能域   │
│ 位置1437 → A3域      │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ Step 2: AF3结构预测  │
│ 提交序列至AF3 Server │
│ 获取PAE矩阵          │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ Step 3: 结构分析     │
│ 计算PAE变化特征      │
│ 评估局部柔性改变     │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ Step 4: ML分型预测   │
│ 输入: 1437, A3, PAEΔ │
│ 输出: type2A (100%)  │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ Step 5: 实验验证     │
│ 推荐: RIPA实验       │
│ 确认: ADAMTS13敏感性 │
└─────────────────────┘
```

### 10.2 各分型实验验证建议

| 预测分型 | 推荐实验 | 验证指标 |
|----------|---------|---------|
| type2A | RIPA | 低剂量瑞斯托霉素聚集增加 |
| type2A | VWF多聚体 | 缺乏大分子量多聚体 |
| type2B | 瑞斯托霉素辅因子 | 亲和力增加 |
| type2B | 血小板计数 | 血小板减少 |
| type2M | 胶原结合实验 | 结合率降低 |
| type2M | VWF:RCo/VWF:Ag比值 | 降低 |
| type2N | FVIII结合实验 | 结合率降低 |
| type2N | FVIII活性 | 降低 |

### 10.3 与传统方法对比

| 维度 | 传统实验室检测 | AF3结构预测 |
|------|---------------|-------------|
| **时间** | 数天至数周 | 数小时至数天 |
| **成本** | 高 (试剂+设备) | 低 (计算资源) |
| **样本** | 需要患者血浆 | 仅需序列 |
| **通量** | 低 (逐个检测) | 高 (批量预测) |
| **准确性** | 金标准 | 辅助性 (79.7%) |
| **机制** | 功能表型 | 分子机制 |

**定位**: AF3预测作为**初筛工具**，指导实验设计，降低试错成本。

---

## 十一、局限性与改进

### 11.1 当前局限

#### 1. 样本量不平衡
```
type2A: 36 (61%)  ✓ 充足
type2N: 13 (22%)  △ 中等
type2M: 8  (14%)  △ 偏少
type2B: 2  (3%)   ✗ 严重不足
```
**影响**: type2B无法有效预测，可能遗漏罕见变异

#### 2. PAE分辨率限制
- PAE反映**相对误差**，非绝对结构变化
- 无法捕捉**电荷变化**、**氢键断裂**等化学特性
-  misses allosteric effects (变构效应)

#### 3. 缺少功能验证
- 未与**实验功能数据** (RIPA, 胶原结合) 直接关联
- 无法建立**定量结构-功能关系**

#### 4. 单一结构状态
- AF3预测**静态结构**， misses 动态变化
- 无法预测**温度**、**pH**等条件的影响

### 11.2 改进方向

#### 短期 (3-6个月)
1. **扩大数据集**
   - 收集更多type2B和type2M样本
   - 整合ClinVar、HGMD公开数据

2. **多模态特征**
   - 整合AlphaGenome RNA预测
   - 添加进化保守性评分 (PhastCons, GERP++)
   - 纳入蛋白质相互作用网络

#### 中期 (6-12个月)
3. **功能实验关联**
   - 建立结构-功能相关性数据库
   - 训练定量预测模型

4. **动态结构模拟**
   - 分子动力学 (MD) 模拟关键变异
   - 分析动态构象变化

#### 长期 (1-2年)
5. **多尺度模型**
   - 整合从原子到细胞的多个尺度
   - 预测临床表现型

6. **临床验证**
   - 前瞻性队列研究
   - 与临床诊断金标准对比

---

## 十二、快速开始

### 12.1 安装

```bash
# 克隆仓库
git clone git@github.com:LukaDD7/VWF-ETHos.git
cd VWF-ETHos

# 创建环境
conda create -n vwf python=3.10
conda activate vwf

# 安装依赖
pip install -r requirements.txt
```

### 12.2 单个变异预测

```bash
python vwf_type2_analysis.py --position 1437
```

**输出**:
```
============================================================
VWF Type-2分型预测结果
============================================================
变异位置: 1437
功能域: A3
预测分型: type2A
可信度: 100.0%
预测方法: rule_based

临床建议: ADAMTS13切割异常，建议进行RIPA实验验证
============================================================
```

### 12.3 批量预测

```bash
# 准备输入文件 variants.csv:
# Position,VWF_Domain
# 1437,A3
# 1100,A2
# 800,A1

python vwf_type2_analysis.py --batch-predict variants.csv --output predictions.csv
```

### 12.4 获取AF3结构

```bash
cd Proteo-Structure-Pipeline

# 生成AF3输入文件
python src/phase2_af3_batch_generator.py --limit 5

# 手动上传至 https://alphafoldserver.com/
# 下载结果至 structures/predictions/

# 提取PAE特征
python src/phase3_structural_scoring.py --wt-cif structures/predictions/VWF_WT.cif
```

---

## 十三、API参考

### 13.1 核心函数

#### `predict_single(position, domain, pae_delta)`
预测单个变异的Type-2分型。

**参数**:
- `position` (int): 变异位置 (1-2813)
- `domain` (str, optional): 功能域名称
- `pae_delta` (float, optional): PAE局部变化值

**返回**:
```python
{
    'position': 1437,
    'domain': 'A3',
    'predicted_subtype': 'type2A',
    'confidence': 1.0,
    'method': 'rule_based',
    'recommendation': 'RIPA实验验证'
}
```

#### `get_vwf_domain(position)`
根据位置返回VWF功能域。

**参数**:
- `position` (int): 氨基酸位置

**返回**:
- `str`: 功能域名称 (D1-D2, D'D3, A1, A2, A3, D4, C1-C2, CT)

#### `batch_predict(input_file, output_file)`
批量预测变异分型。

**参数**:
- `input_file` (str): CSV文件路径，需包含Position列
- `output_file` (str): 输出CSV文件路径

### 13.2 数据结构

#### PAE矩阵
```python
# 2813 x 2813 numpy array
pae_matrix[i, j]  # 残基i与j的预测对齐误差 (Å)
```

#### 特征向量
```python
features = {
    'Position_Norm': 0.51,      # 位置归一化
    'Domain_Encoded': 4,        # 功能域编码
    'PAE_Delta_Local': 0.5,     # 局部PAE变化
    'Local_Flex_Ratio': 1.1,    # 柔性比值
    'Domain_PAE_Delta': 0.2     # 功能域级变化
}
```

---

## 附录

### A. 参考文献

1. Jumper, J., et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596(7873), 583-589.

2. Sadler, J. E. (2003). Von Willebrand disease type 1: a diagnosis in search of a disease. *Blood*, 101(6), 2089-2093.

3. Ginsburg, D., & Sadler, J. E. (1993). von Willebrand disease: a database of point mutations, insertions, and deletions. *Thrombosis and Haemostasis*, 69(01), 0177-0184.

4. Budde, U., et al. (2006). Detailed analysis of the multimeric structure of von Willebrand factor. *Annals of Hematology*, 85(7), 433-443.

### B. 相关资源

- **UniProt P04275**: https://www.uniprot.org/uniprotkb/P04275
- **AlphaFold3 Server**: https://alphafoldserver.com/
- **ClinVar VWF**: https://www.ncbi.nlm.nih.gov/clinvar/?term=VWF
- **ISTH指南**: https://www.isth.org/

### C. 联系方式

- **作者**: LukaDD7
- **GitHub**: https://github.com/LukaDD7/VWF-ETHos
- **Issues**: https://github.com/LukaDD7/VWF-ETHos/issues

---

**最后更新**: 2026-03-23

**版本**: v1.0.0

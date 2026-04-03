# VWF Residue-Level Feature Extractor - 使用说明

## 概述

本模块提供基于文献精确注释的VWF残基级特征提取功能。与域级粗粒度特征不同，本模块提取：

- **精确残基位置** (如 AIM的盐桥D1263-R1668)
- **突变氨基酸性质** (大小、电荷、疏水性变化)
- **结构特征** (全局RMSD + 10Å局部RMSD)
- **文献引用** (每个特征都有PMID支持)

## 核心特性

| 特征类别 | 精度 | 文献来源 |
|---------|------|---------|
| AIM (自抑制模块) | 残基级精确接触 | PMID: 33888542, 28325766 |
| GPIbα界面 | 具体相互作用类型 | PMID: 12191960, 15117959 |
| ADAMTS13切割位点 | Exosite 1/2/3精确残基 | PMID: 28904067, 15758837 |
| 钙结合位点 | 配位残基列表 | PMID: 28904067 |
| 胶原结合 | 关键界面残基 | PMID: 15758837, 16079166 |
| FVIII结合 | 关键结合残基 | PMID: 33888542, 15758837 |

## 快速开始

### 1. 单个变异特征提取

```python
from vwf_residue_feature_extractor import VWFResidueFeatureExtractor

# 初始化提取器
extractor = VWFResidueFeatureExtractor()

# 提取单个变异特征
features = extractor.extract_residue_features(
    variant_id="R1306W",
    position=1306,
    ref_aa="R",
    alt_aa="W",
    structural_data={
        "plddt_wt": 85.5,
        "plddt_mut": 78.3,
        "plddt_delta": -7.2,
        "global_rmsd": 0.8,
        "local_rmsd_10a": 1.2
    }
)

# 访问特征
print(f"Domain: {features.domain}")
print(f"Is GPIb interface: {features.is_in_gpib_interface}")
print(f"AIM disruption score: {features.aim_disruption_score}")
print(f"Charge change: {features.mutation_charge_change}")
```

### 2. 批量特征提取

```python
# 批量提取
variants = [
    {"variant_id": "R1306W", "position": 1306, "ref_aa": "R", "alt_aa": "W"},
    {"variant_id": "V1316M", "position": 1316, "ref_aa": "V", "alt_aa": "M"},
    {"variant_id": "D1614N", "position": 1614, "ref_aa": "D", "alt_aa": "N"},
]

results = extractor.extract_features_batch(variants)

for feat in results:
    print(f"{feat.variant_id}: {feat.domain}, {feat.gpib_interaction_type}")
```

### 3. 生成特征报告

```python
# 生成人类可读报告
report = extractor.generate_feature_report(features)
print(report)
```

### 4. 导出JSON

```python
import json

# 转换为字典/JSON
feature_dict = features.to_dict()
print(json.dumps(feature_dict, indent=2))
```

## 特征详解

### 1. AIM (Autoinhibitory Module) 特征

**文献来源**: PMID 33888542 (Arce et al. 2021, Nat Commun)

```python
features.is_in_AIM              # 是否在AIM内
features.aim_component          # "N_terminal" 或 "C_terminal"
features.aim_key_contact        # 关键接触详情
features.aim_disruption_score   # AIM破坏评分 (0-1)
```

**精确残基注释**:
- **N端**: Gln1238-His1268
  - D1263-R1668: 盐桥 (AIM稳定的关键)
  - Y1469-Q1659: 氢键
- **C端**: Leu1460-Asp1472

**AIM破坏评分计算**:
- 盐桥残基电荷改变: +0.9
- 疏水核心引入电荷: +0.5
- 大残基空间位阻: +0.3-0.5

### 2. GPIbα界面特征

**文献来源**: PMID 12191960 (Huizinga et al. 2002, Science)

```python
features.is_in_gpib_interface   # 是否在GPIbα界面
features.gpib_site              # "interactive_site_1" 或 "interactive_site_2"
features.gpib_key_residue       # 是否是关键残基
features.gpib_interaction_type  # "cation_pi", "hbond", "hydrophobic"
```

**关键残基**:
- **R1306**: cation-π相互作用 (与GPIbα Y283)
  - VWD突变: R1306W, R1306Q, R1306L
  - 机制: 破坏关键结合界面
- **I1299**: 主链氢键

### 3. ADAMTS13切割位点特征

**文献来源**: PMID 28904067 (Crawley et al. 2020, Blood)

```python
features.is_scissile_bond       # 是否是Y1605-M1606切割键
features.is_in_exosite_1        # β4-less loop (1594-1602)
features.is_in_exosite_2        # 与cysteine-rich domain结合 (1642-1651)
features.is_in_exosite_3        # 与spacer domain结合 (1660-1668)
features.is_group1_2A_mutation  # Group 1: 促进A2展开
features.is_group2_2A_mutation  # Group 2: 延迟A2重折叠
```

**Group 1 mutations** (促进切割):
- M1528V: 促进A2过早展开
- E1638K: 破坏A2结构稳定

**Group 2 mutations** (延迟重折叠):
- R1597W: 延迟重折叠，增加切割时间
- D1614N: 破坏exosite 1相互作用

### 4. 结构特征 (从AF3提取)

```python
features.plddt_wt               # WT pLDDT
features.plddt_mut              # Mutant pLDDT
features.plddt_delta            # pLDDT变化 (负值表示不稳定)
features.global_rmsd            # 全局RMSD (Å)
features.local_rmsd_10a         # 10Å半径局部RMSD (Å)
```

**RMSD解释**:
- **Global RMSD**: 全蛋白Cα叠合后的均方根偏差
  - <1.0 Å: 结构高度保守
  - 1.0-2.0 Å: 中等结构变化
  - >2.0 Å: 显著结构重排

- **Local RMSD (10Å)**: 突变位点10Å半径内的RMSD
  - <0.5 Å: 局部结构保守
  - 0.5-1.5 Å: 局部扰动
  - >1.5 Å: 显著局部重排

### 5. 突变氨基酸性质

```python
features.mutation_size_delta             # 体积变化 (Å³)
features.mutation_hydrophobicity_delta   # 疏水性变化 (Kyte-Doolittle)
features.mutation_charge_change          # 电荷变化 (-2 to +2)
features.mutation_aromatic_change        # 芳香性变化 (-1 to +1)
features.mutation_helix_propensity_delta # 螺旋倾向变化
```

**性质表**:

| 氨基酸 | 体积(Å³) | 疏水性 | 电荷 | 芳香性 |
|--------|---------|--------|------|--------|
| R | 173.4 | -4.5 | +1 | No |
| W | 227.8 | -0.9 | 0 | Yes |
| D | 111.1 | -3.5 | -1 | No |
| E | 138.4 | -3.5 | -1 | No |

## 文献引用库

### 核心文献

1. **AIM结构与机制**
   - PMID: 33888542 - Arce et al. (2021) Nat Commun
   - PMID: 28325766 - Deng et al. (2017) J Thromb Haemost
   - PMID: 28904067 - Crawley et al. (2020) Blood

2. **GPIbα界面**
   - PMID: 12191960 - Huizinga et al. (2002) Science
   - PMID: 15117959 - Dumas et al. (2004) JBC
   - PMID: 23341617 - Blenner et al. (2014) JBC

3. **ADAMTS13切割位点**
   - PMID: 15758837 - Zhang et al. (2009) Science
   - PMID: 16079166 - Siedlecki et al. (2015) Blood

4. **FVIII结合**
   - PMID: 33888542 - Lenting et al. (2024) Blood
   - PMID: 15758837 - Zhang et al. (2009) Science

5. **胶原结合**
   - PMID: 15758837 - Zhang et al. (2009) Science
   - PMID: 16079166 - Siedlecki et al. (2015) Blood

## 与Pipeline集成

### 1. 替换现有特征提取

在 `vwf_type2_domain_pipeline.py` 中:

```python
# 旧代码
from vwf_structure_feature_extractor import VWFStructureFeatureExtractor
extractor = VWFStructureFeatureExtractor()
features = extractor.extract_features(structure_file, position)

# 新代码
from vwf_residue_feature_extractor import VWFResidueFeatureExtractor
extractor = VWFResidueFeatureExtractor()

# 提取残基级特征 + 结构特征
structural_data = {
    "plddt_wt": wt_plddt,
    "plddt_mut": mut_plddt,
    "plddt_delta": plddt_delta,
    "global_rmsd": global_rmsd,
    "local_rmsd_10a": local_rmsd_10a
}

features = extractor.extract_residue_features(
    variant_id=variant_id,
    position=position,
    ref_aa=ref_aa,
    alt_aa=alt_aa,
    structural_data=structural_data
)
```

### 2. 为分类器提供特征

```python
# 将特征转换为分类器输入
feature_vector = {
    # 残基级特征
    "is_in_aim": int(features.is_in_AIM),
    "aim_disruption_score": features.aim_disruption_score,
    "is_gpib_key_residue": int(features.gpib_key_residue),
    "is_scissile_bond": int(features.is_scissile_bond),
    "is_group2_2a": int(features.is_group2_2A_mutation),

    # 突变性质
    "size_delta": features.mutation_size_delta,
    "charge_change": features.mutation_charge_change,
    "hydrophobicity_delta": features.mutation_hydrophobicity_delta,

    # 结构特征
    "plddt_delta": features.plddt_delta,
    "local_rmsd_10a": features.local_rmsd_10a,
}
```

## 输出示例

### R1306W (Type 2B经典突变)

```json
{
  "variant_id": "R1306W",
  "position": 1306,
  "mutation": "R1306W",
  "domain": "A1",
  "relative_position_in_domain": 0.158,

  "gpib_interface_features": {
    "is_in_gpib_interface": true,
    "gpib_site": "interactive_site_1",
    "gpib_key_residue": true,
    "gpib_interaction_type": "cation_pi"
  },

  "vwd_hotspot": {
    "is_hotspot": true,
    "hotspot_type": "Type_2B",
    "known_mutation": "R1306W, R1306Q, R1306L"
  },

  "mutation_properties": {
    "size_delta": 54.4,
    "hydrophobicity_delta": 3.6,
    "charge_change": -1,
    "aromatic_change": 1
  },

  "literature": {
    "pmids": ["12191960"],
    "evidence": "R1306-Y283 cation-pi interaction critical for GPIb binding"
  }
}
```

### D1614N (Group 2 Type 2A)

```json
{
  "variant_id": "D1614N",
  "position": 1614,
  "mutation": "D1614N",
  "domain": "A2",

  "adamts13_features": {
    "is_scissile_bond": false,
    "is_in_exosite_1": true,
    "is_group2_2A": true
  },

  "mutation_properties": {
    "charge_change": 1,
    "size_delta": 3.0
  },

  "literature": {
    "evidence": "Exosite 1 interaction with disintegrin domain disrupted"
  }
}
```

## 注意事项

1. **文献覆盖范围**
   - 当前覆盖主要功能域 (A1, A2, A3, D'D3, C4)
   - D4, CK等域的精细注释需要更多文献支持

2. **结构数据可选**
   - 残基级特征不依赖AF3结构
   - RMSD/plddt特征需要结构数据

3. **AIM范围**
   - N端: 1238-1268
   - C端: 1460-1472
   - 1306在GPIb界面但**不在**AIM内

4. **扩展性**
   - 可通过修改 `LITERATURE_ANNOTATIONS` 字典添加新特征
   - 每个特征必须有PMID支持

## 更新日志

- **2026-04-03**: 初始版本
  - 实现残基级精确特征提取
  - 整合3篇核心综述的文献注释
  - 支持全局和局部RMSD
  - 提供完整文献引用

# VWF Type-2 分型诊断分析 Pipeline

基于文献知识和AlphaFold3结构预测的VWF Type-2变异分型分析流程。

## 文献依据

本pipeline严格依据以下3篇关键文献构建：

1. **Lenting et al. (2024) Blood** - "How unique structural adaptations support and coordinate the complex function of von Willebrand factor"
   - 详细阐述VWF结构域与功能的关系
   - A1域GPIbα结合机制与自抑制模块(AIM)
   - A2域ADAMTS13切割位点与Ca2+结合位点

2. **Atiq & O'Donnell (2024) Blood** - "Novel functions for von Willebrand factor"
   - VWF全长的配体结合位点图谱
   - 新发现的VWF功能区域

3. **Haberichter & O'Donnell (2026) Haematologica** - "Structure and multiple functions of von Willebrand factor"
   - VWF生物合成与分泌机制
   - 结构域最新注释（Zhou et al. 2012更新版）

## Type-2 分型与结构域对应关系

| 分型 | 主要涉及域 | 关键功能位点 | 分子机制 |
|------|------------|--------------|----------|
| **Type 2A** | A2域 (1493-1684) | Y1605-M1606切割位点<br>Ca2+结合位点(1596-1602)<br>临近二硫键(C1669-C1670) | ADAMTS13敏感性增加<br>高分子量多聚体丢失 |
| **Type 2B** | A1域 (1271-1492) | 自抑制模块(AIM)<br>N端: G1238-H1268<br>C端: L1460-D1472 | AIM破坏导致<br>自发性GPIbα结合 |
| **Type 2M** | A1域 或 A3域 | A1: GPIbα结合界面<br>A3: 胶原蛋白结合界面 | 结合能力降低<br>但多聚体正常 |
| **Type 2N** | D'/D3域 (764-1233) | D'区: R782-C799<br>Til'结构: R816 | FVIII结合能力降低<br>FVIII半衰期缩短 |

## 文件结构

```
domain_analysis/
├── vwf_type2_literature_based_classifier.py   # 基于文献的分类器
├── vwf_structure_feature_extractor.py         # 结构特征提取器
├── vwf_type2_domain_pipeline.py              # 主分析pipeline
└── README.md                                  # 本文档
```

## 使用方法

### 1. 运行完整分析

```bash
python domain_analysis/vwf_type2_domain_pipeline.py \
    --variants /path/to/VWF_Type2_variants.csv \
    --structures /path/to/AF3_Results/extracted/ \
    --output ./domain_analysis_results/ \
    --wt-structure /path/to/VWF_WT.cif
```

### 2. 单独使用分类器

```python
from domain_analysis.vwf_type2_literature_based_classifier import VWFType2Classifier, VWFVariant

classifier = VWFType2Classifier()

# 根据位置查询结构域和功能
result = classifier.classify_by_position(1306)
print(f"Domain: {result['domain']}")
print(f"Type-2 associations: {result['associations']}")
```

### 3. 提取结构特征

```python
from domain_analysis.vwf_structure_feature_extractor import VWFFeatureExtractor

extractor = VWFFeatureExtractor(wt_structure_path="VWF_WT.cif")

# 提取所有结构域特征
domain_features = extractor.extract_domain_features("mutant.cif")

# 计算突变影响
impact = extractor.calculate_mutation_impact(
    variant_id="VWF_R1306W",
    position=1306,
    ref_aa="R",
    alt_aa="W",
    mut_structure_path="fold_vwf_r1306w.cif"
)
```

## 输入数据格式

### 变异CSV文件格式

```csv
variant_id,protein_change,position,ref_aa,alt_aa,acmg_classification,type2_subtype
VWF_R1306W,p.Arg1306Trp,1306,R,W,Likely_pathogenic,2B
VWF_M1606R,p.Met1606Arg,1606,M,R,Pathogenic,2A
...
```

## 输出结果

### 1. 分类结果CSV (`type2_classification_results.csv`)

包含以下列：
- `variant_id`: 变异标识符
- `protein_change`: 蛋白质变化
- `position`: 氨基酸位置
- `domain`: 所属结构域
- `known_subtype`: 已知分型（如有）
- `predicted_subtype`: 预测分型
- `prediction_confidence`: 预测置信度
- `is_functional_site`: 是否位于功能位点
- `is_hotspot`: 是否位于已知热点
- `plddt_delta`: pLDDT变化
- `local_rmsd`: 局部RMSD

### 2. 分析报告 (`analysis_report.txt`)

包含：
- 汇总统计
- 结构域分布
- 已知vs预测分型对比
- 结构域特异性分析
- 详细变异分类结果

## 分类算法逻辑

### Type 2A 判定
- 位置在A2域 (1493-1684)
- 接近切割位点 (<10残基)
- 影响Ca2+结合或临近二硫键
- pLDDT显著下降 (<-10)

### Type 2B 判定
- 位置在A1域 (1271-1492)
- 位于AIM区域 (1238-1268 或 1460-1472)
- 同时影响GPIbα结合界面
- 电荷改变（破坏静电相互作用）

### Type 2M 判定
- A1域：影响GPIbα结合但不破坏AIM
- A3域：影响胶原蛋白结合界面 (1700-1850)
- 多聚体正常（无A2域异常）

### Type 2N 判定
- 位置在D'/D3域 (764-1233)
- 位于FVIII结合界面 (782-799, 816-826)
- 破坏与FVIII轻链相互作用

## 特征提取说明

### 1. pLDDT分析
- **全局pLDDT**: 全蛋白平均置信度
- **局部pLDDT**: 突变位点5残基窗口平均
- **pLDDT变化**: 突变前后置信度差异

### 2. 结构稳定性
- **局部RMSD**: 突变位点附近结构偏差
- **B因子**: 原子热运动参数
- **低置信区域**: pLDDT<70的连续区域

### 3. 氨基酸性质变化
- **疏水性变化**: 基于Kyte-Doolittle标度
- **电荷变化**: 酸碱性质改变
- **大小变化**: 侧链体积差异

### 4. 功能位点注释
- 与功能位点的距离
- 热点区域匹配
- 结构域归属

## 参考文献

1. Lenting PJ, Denis CV, Christophe OD. Blood. 2024;144(21):2174-2184.
2. Atiq F, O'Donnell JS. Blood. 2024;144(12):1247-1256.
3. Haberichter SL, O'Donnell JS. Haematologica. 2026;in press.
4. Zhou YF et al. Blood. 2012;120(2):449-458.
5. Huizinga EG et al. Science. 2002;297(5584):1176-1179.
6. Zhang Q et al. PNAS. 2009;106(23):9226-9231.

## 开发说明

这是一个分支开发的功能模块，专注于基于文献知识的VWF Type-2分型诊断。
与其他主流程代码保持隔离，不干扰现有分析流程。

## 联系

如有问题或建议，请参考CLAUDE.md中的项目说明。

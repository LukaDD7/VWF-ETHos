# VWF残基级特征提取系统 - 快速入门指南

## 环境要求

```bash
# 激活环境
conda activate alphagenome

# 确认依赖
python3 -c "import Bio; import numpy; import pandas; print('OK')"
```

## 5分钟快速开始

### 1. 提取单个变异特征

```python
from vwf_residue_feature_extractor import VWFResidueFeatureExtractor

# 初始化
extractor = VWFResidueFeatureExtractor()

# 提取R1306W特征
features = extractor.extract_residue_features(
    variant_id="R1306W",
    position=1306,
    ref_aa="R",
    alt_aa="W"
)

# 查看关键特征
print(f"Domain: {features.domain}")                    # A1
print(f"GPIb interface: {features.is_in_gpib_interface}")  # True
print(f"Interaction: {features.gpib_interaction_type}")    # cation_pi
print(f"Charge change: {features.mutation_charge_change}") # -1
print(f"Literature: {features.literature_pmids}")          # ['12191960']
```

### 2. 计算RMSD

```python
from vwf_rmsd_calculator import VWFRMSDCalculator
from pathlib import Path

calculator = VWFRMSDCalculator()

wt_cif = Path("fold_vwf_wt/fold_vwf_wt_model_0.cif")
mut_cif = Path("fold_vwf_r1306w/fold_vwf_r1306w_model_0.cif")

# 计算多尺度RMSD
results = calculator.calculate_mutation_site_rmsd(
    wt_cif, mut_cif, position=1306
)

print(f"Global RMSD: {results['rmsd_global']:.2f} Å")
print(f"Local 10Å RMSD: {results['rmsd_10a']:.2f} Å")
```

### 3. 批量处理

```python
from vwf_integrated_pipeline import IntegratedVWFPipeline
import pandas as pd

# 准备变异列表
variants = pd.DataFrame([
    {"variant_id": "R1306W", "position": 1306, "ref_aa": "R", "alt_aa": "W"},
    {"variant_id": "V1316M", "position": 1316, "ref_aa": "V", "alt_aa": "M"},
    {"variant_id": "D1614N", "position": 1614, "ref_aa": "D", "alt_aa": "N"},
])

# 初始化Pipeline
pipeline = IntegratedVWFPipeline(
    af3_base_dir="/path/to/af3_results"
)

# 批量处理
results = pipeline.process_batch(variants)
results.to_csv("features.csv", index=False)
```

### 4. 生成可视化

```python
from vwf_residue_visualizer import VWFResidueVisualizer

visualizer = VWFResidueVisualizer()

# 加载特征数据
import pandas as pd
features_df = pd.read_csv("features.csv")

# 生成完整报告
visualizer.create_comprehensive_report(features_df)

# 或生成单个图表
visualizer.plot_domain_landscape(features_df)
visualizer.plot_feature_heatmap(features_df)
```

## 核心特征说明

### 文献特征

| 特征 | 类型 | 说明 |
|------|------|------|
| `is_in_AIM` | bool | 是否在自抑制模块内 |
| `is_gpib_interface` | bool | 是否在GPIbα结合界面 |
| `gpib_interaction_type` | str | 相互作用类型 (cation_pi/hbond/hydrophobic) |
| `is_scissile_bond` | bool | 是否是ADAMTS13切割键 |
| `is_exosite_1/2/3` | bool | 是否在ADAMTS13 exosite区域 |
| `is_vwd_hotspot` | bool | 是否是已知VWD热点 |
| `literature_pmids` | list | 文献引用PMID列表 |

### 结构特征

| 特征 | 类型 | 说明 |
|------|------|------|
| `plddt_mean` | float | 平均pLDDT (0-100) |
| `plddt_at_site` | float | 突变位点pLDDT |
| `plddt_delta` | float | 与WT的pLDDT差异 |
| `rmsd_global` | float | 全局RMSD (Å) |
| `rmsd_10a` | float | 10Å局部RMSD |
| `ptm_score` | float | pTM质量分数 |

### 突变性质

| 特征 | 类型 | 说明 |
|------|------|------|
| `mutation_size_delta` | float | 体积变化 (Å³) |
| `mutation_charge_change` | int | 电荷变化 (-2到+2) |
| `mutation_hydrophobicity_delta` | float | 疏水性变化 |

## 常见问题

### Q: 如何处理大量变异？
A: 使用`pipeline.process_batch()`，支持DataFrame批量处理。

### Q: 没有AF3数据怎么办？
A: 系统会自动检测，只返回文献特征，结构特征标记为None。

### Q: 如何添加新的文献注释？
A: 编辑`vwf_residue_feature_extractor.py`中的注释字典，添加新的PMID和残基信息。

### Q: RMSD计算太慢？
A: 可以只计算局部RMSD (5Å或10Å)，减少计算量。

## 完整示例

```python
# 完整流程示例
from vwf_integrated_pipeline import IntegratedVWFPipeline
from vwf_residue_visualizer import VWFResidueVisualizer
import pandas as pd

# 1. 准备数据
variants = pd.read_csv("my_variants.csv")  # columns: variant_id, position, ref_aa, alt_aa

# 2. 提取特征
pipeline = IntegratedVWFPipeline(af3_base_dir="./af3_results")
features = pipeline.process_batch(variants)

# 3. 保存结果
features.to_csv("extracted_features.csv", index=False)

# 4. 生成报告
visualizer = VWFResidueVisualizer()
visualizer.create_comprehensive_report(features)

# 5. 查看特定特征
high_impact = features[features['aim_disruption_score'] > 0.7]
print(f"High impact variants: {len(high_impact)}")
```

## 技术支持

- 详细API文档: `RESIDUE_FEATURE_EXTRACTOR_GUIDE.md`
- 技术报告: `TECHNICAL_REPORT.md`
- 代码示例: 每个模块的`main()`函数

## 引用

如果使用本系统，请引用：
- 核心文献PMID: 12191960, 33888542, 28904067, 35148377
- AlphaFold3: PMID 38706376 (Nature 2024)
- BioPython: PMID 19304878

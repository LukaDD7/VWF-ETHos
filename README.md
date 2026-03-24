# VWF-ETHos: VWF Type-2 von Willebrand Disease Diagnostic System

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub](https://img.shields.io/badge/github-LukaDD7-green.svg)](https://github.com/LukaDD7)

基于AlphaFold3结构和AlphaGenome功能基因组预测的Type-2血管性血友病(VWD)分型诊断系统

---

## 📋 目录

- [项目简介](#项目简介)
- [核心功能](#核心功能)
- [快速开始](#快速开始)
- [技术文档](#技术文档)
- [分析结果](#分析结果)
- [贡献指南](#贡献指南)
- [引用](#引用)

---

## 项目简介

VWF-ETHos是一个整合多维度数据的VWF Type-2分型诊断平台：

- **结构维度**: AlphaFold3蛋白质结构预测 (PAE矩阵分析)
- **功能维度**: AlphaGenome功能基因组预测 (11种模态)
- **临床维度**: 基于100个Type-2变异的分型参考表
- **智能维度**: 机器学习辅助分型预测 (79.7%准确率)

### 核心发现

- **位置决定分型**: 功能域位置是最重要的预测特征 (97.7%)
- **A3域 → type2A**: 100%准确率
- **A2域 → type2M**: 57%准确率 (需验证)
- **93.8%的Type-2变异**: 剪接影响极小 (<0.01)

---

## 核心功能

### 1. Type-2分型预测

```bash
# 单个变异预测
python vwf_type2_analysis.py --position 1437

# 批量预测
python vwf_type2_analysis.py --batch-predict variants.csv
```

### 2. AlphaGenome功能分析

```bash
# 全模态分析 (11 modalities)
python scripts/07e_god_mode_epigenome_crawler.py
```

### 3. AF3结构预测Pipeline

```bash
cd Proteo-Structure-Pipeline
python src/phase1_smart_filter.py
python src/phase2_af3_batch_generator.py
python src/phase3_structural_scoring.py
```

---

## 快速开始

### 安装

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

### 使用示例

```python
from vwf_type2_analysis import predict_single

# 预测单个变异
result = predict_single(position=1437)
print(f"预测分型: {result['predicted_subtype']}")
print(f"置信度: {result['confidence']:.1%}")
print(f"临床建议: {result['recommendation']}")
```

**输出**:
```
预测分型: type2A
置信度: 100.0%
临床建议: ADAMTS13切割异常，建议进行RIPA实验验证
```

---

## 技术文档

### 详细文档

| 文档 | 内容 | 链接 |
|------|------|------|
| **Type-2分析文档** | 完整技术流程、ML模型、生物学原理 | [Type-2-Analysis-README.md](Type-2-Analysis-README.md) |
| **技术实现文档** | 代码详解、API参考、扩展开发 | [TECHNICAL_DOCUMENTATION.md](TECHNICAL_DOCUMENTATION.md) |
| **AlphaGenome报告** | Type-2 AlphaGenome特征分析 | [AlphaGenome_Type2_Analysis_Report.md](AlphaGenome_Type2_Analysis_Report.md) |
| **VWF参考表说明** | Type-2变异参考表使用指南 | [VWF_Type2_AF3_Reference_Table_README.md](VWF_Type2_AF3_Reference_Table_README.md) |

### 代码结构

```
vwf_type2_analysis.py           # 分型预测主程序
├── get_vwf_domain()            # 功能域映射
├── predict_subtype_by_rules()  # 基于规则预测
├── predict_single()            # 单个预测
└── batch_predict()             # 批量预测

scripts/
├── 01_filter_target_vus.py     # 变异过滤
├── 03_run_alphagenome_inference.py  # AlphaGenome API
├── 07e_god_mode_epigenome_crawler.py # 全模态爬取
└── ...

Proteo-Structure-Pipeline/src/
├── phase1_smart_filter.py      # 智能过滤
├── phase2_af3_batch_generator.py # AF3批次生成
└── phase3_structural_scoring.py  # 结构评分
```

---

## 分析结果

### 1. Type-2分型特征

| 分型 | 样本数 | AlphaGenome分数 | 主要功能域 | 预测准确率 |
|------|--------|----------------|-----------|-----------|
| type2A | 42 | 4.37 ± 5.85 | A3 (100%) | 75.0% |
| type2M | 25 | 3.50 ± 0.72 | A2/D4 | 100% |
| type2N | 20 | 5.44 ± 8.61 | A1 | 92.3% |
| type2B | 12 | 3.40 ± 0.77 | A1 | 样本不足 |

### 2. ML特征重要性

```
Position:           ████████████████████████████████████████████  97.7%
Domain:             ████████████████████████████████████          41.5%
PAE变化:            ████                                           4.5%
AlphaGenome分数:    ▏                                              0.0%
```

### 3. 可视化结果

- `figures/alphagenome_analysis/` - AlphaGenome特征图
- `figures/` - 其他分析图表

---

## 贡献指南

### 如何贡献

1. **Fork** 本仓库
2. 创建 **Feature Branch** (`git checkout -b feature/AmazingFeature`)
3. **Commit** 更改 (`git commit -m 'Add some AmazingFeature'`)
4. **Push** 到分支 (`git push origin feature/AmazingFeature`)
5. 提交 **Pull Request**

### 报告问题

请使用 [GitHub Issues](https://github.com/LukaDD7/VWF-ETHos/issues) 报告问题或请求功能。

### 开发规范

- 遵循 [PEP 8](https://www.python.org/dev/peps/pep-0008/) 代码规范
- 提交前运行测试
- 更新相关文档

---

## 引用

### 引用本项目

```bibtex
@software{vwf_ethos2026,
  author = {LukaDD7},
  title = {VWF-ETHos: Type-2 VWD Diagnostic System},
  year = {2026},
  url = {https://github.com/LukaDD7/VWF-ETHos}
}
```

### 相关文献

1. Jumper, J., et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596(7873), 583-589.

2. Sadler, J. E. (2003). Von Willebrand disease type 1: a diagnosis in search of a disease. *Blood*, 101(6), 2089-2093.

3. Ginsburg, D., & Sadler, J. E. (1993). von Willebrand disease: a database of point mutations, insertions, and deletions. *Thrombosis and Haemostasis*, 69(01), 0177-0184.

4. Budde, U., et al. (2006). Detailed analysis of the multimeric structure of von Willebrand factor. *Annals of Hematology*, 85(7), 433-443.

---

## 许可证

本项目采用 [MIT License](LICENSE) 许可。

---

## 联系方式

- **作者**: LukaDD7
- **GitHub**: https://github.com/LukaDD7/VWF-ETHos
- **邮箱**: [待添加]

---

**最后更新**: 2026-03-24

**版本**: v1.0.0

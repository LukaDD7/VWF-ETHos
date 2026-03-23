# VWF-ETHos: VWF Type-2 von Willebrand Disease Diagnostic System

基于AlphaFold3结构的Type-2血管性血友病(VWD)分型诊断系统

## 项目简介

本项目整合AlphaFold3结构预测和机器学习，建立VWF Type-2分型的辅助诊断工具。

## 核心功能

1. **结构域映射**: 基于UniProt P04275的VWF功能域注释
2. **结构特征提取**: 从AF3 PAE矩阵提取局部柔性特征
3. **分型预测**: 基于位置和结构特征的Type-2分型预测

## 使用方法

```bash
python vwf_type2_analysis.py --position 1437
```

## 分型预测规则

| 功能域 | 预测分型 | 可信度 |
|--------|---------|--------|
| A3 | type2A | 100% |
| A2 | type2M | 57% |
| A1 | type2N | 63% |
| D4 | type2M | 80% |

## 文件结构

- `vwf_type2_analysis.py`: 诊断脚本
- `VWF_Type2_AF3_Reference_Table.xlsx`: Type-2变异参考表
- `scripts/`: AlphaGenome分析脚本
- `results/`: 分析结果

## 机器学习模型

- 算法: Gradient Boosting
- 准确率: 79.7%
- 关键特征: 位置(50.1%) + 功能域(41.5%)

## 引用

若使用本项目，请引用：
- AlphaFold3: https://alphafoldserver.com/
- VWF Type-2分型标准: ISTH指南

## 作者

LukaDD7

## 许可

MIT License

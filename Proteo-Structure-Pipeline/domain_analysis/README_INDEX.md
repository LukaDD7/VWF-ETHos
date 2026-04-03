# VWF Residue-Level Feature Extraction System

## 📋 项目导航

### 🚀 快速开始
- **[QUICKSTART.md](QUICKSTART.md)** - 5分钟快速入门指南

### 📚 文档
- **[TECHNICAL_REPORT.md](TECHNICAL_REPORT.md)** - 详细技术报告
- **[README_RESIDUE_FEATURES.md](README_RESIDUE_FEATURES.md)** - 功能详细说明
- **[RESIDUE_FEATURE_EXTRACTOR_GUIDE.md](RESIDUE_FEATURE_EXTRACTOR_GUIDE.md)** - API指南
- **[PROJECT_COMPLETION_SUMMARY.md](PROJECT_COMPLETION_SUMMARY.md)** - 项目完成总结

### 💻 核心代码
- `vwf_residue_feature_extractor.py` - 残基级文献特征提取
- `af3_cif_parser.py` - AF3 CIF/JSON解析
- `vwf_rmsd_calculator.py` - BioPython RMSD计算
- `vwf_integrated_pipeline.py` - 整合Pipeline
- `vwf_residue_visualizer.py` - 可视化模块

## 🚀 5分钟快速开始

```bash
conda activate alphagenome
cd domain_analysis

# 测试所有功能
python3 vwf_residue_feature_extractor.py
python3 vwf_rmsd_calculator.py
python3 vwf_integrated_pipeline.py
python3 vwf_residue_visualizer.py
```

## 📊 测试验证

- ✅ 59个Type 2变异成功处理
- ✅ BioPython RMSD计算验证
- ✅ 可视化图表生成

## 📧 项目信息

**完成日期**: 2026-04-03  
**状态**: ✅ 生产就绪

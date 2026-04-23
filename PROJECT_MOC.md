# VWF Type-2 Diagnostic System - Project MOC (Map of Content)

本文件为您提供项目各功能模块的代码映射，以便于后续按科研智能体（Agent Harness）范式进行重构。

---

## 1. 核心推理与集成层 (The Brain)
*   **`agentic_vwf_classifier.py`**: **智能体分类器核心**。实现 3-Expert (结构、转录、逻辑融合) 架构，包含 Rule 6 等核心分型逻辑。
*   **`merge_alpha_features.py`**: **多模态特征融合器**。将 AlphaGenome 的 CSV 输出、AF3 的 JSON 结构输出、以及原始 Type-2 标签合并为统一的 Parquet 矩阵。
*   **`vwf_type2_analysis.py`**: **分型统计分析脚本**。用于生成混淆矩阵、计算召回率及可视化初版分析结果。

## 2. AlphaGenome 预测管线 (Transcriptomic Engine)
*   **`scripts/01_filter_target_vus.py`**: **目标变异筛选器**。从 ClinVar/HGMD 中提取 VWF 的 VUS 和致病位点。
*   **`scripts/02_preprocess_and_liftover.py`**: **坐标转换工具**。执行 hg19 到 hg38 的 Liftover，并构建 1Mb 的基因组推理区间。
*   **`scripts/03_run_alphagenome_inference.py`**: **AG 推理驱动器**。通过 API 批量调用 AlphaGenome 模型。
*   **`scripts/04_analyze_and_visualize.py`**: **预测结果可视化**。生成 4 面板临床级可视化 Track Plots。

## 3. Proteo-Structure-Pipeline (Physical Structural Engine)
*   **`Proteo-Structure-Pipeline/src/phase1_smart_filter.py`**: **结构任务筛选器**。智能挑选最需要进行 AF3 模拟的位点。
*   **`Proteo-Structure-Pipeline/src/phase2_af3_batch_generator.py`**: **AF3 Batch 生成器**。将表格变异转换为 AlphaFold Server 官方 JSON 格式。
*   **`Proteo-Structure-Pipeline/src/phase3_structural_scoring.py`**: **结构评分器**。从预测结果中提取 pLDDT 和 PAE 界面扰动特征。
*   **`Proteo-Structure-Pipeline/src/regenerate_af3_batches.py`**: **任务补丁工具**。针对失败或新增位点重新生成补全批次。
*   **`Proteo-Structure-Pipeline/src/table_normalizer.py`**: **表格标准化模块**。统一处理各种格式的 VWF 变异表格，自动检测格式并转换为标准输入格式。

---

## 附录：数据格式规范与 AF3 转换说明

### A. 表格格式标准

#### 标准输入列
| 列名 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `AA_Position` | int | ✓ | 氨基酸位置 (1-2813) |
| `WT_AA_1` | str | ✓ | 野生型氨基酸 (1-letter, 如 'R') |
| `Mut_AA_1` | str | ✓ | 突变型氨基酸 (1-letter, 如 'W') |
| `Domain` | str | | 功能域 (如 'A1', 'D3') |
| `Original` | str | | 原始格式 (如 'Arg1341Trp') |

#### 支持的输入格式
1. **3-letter格式**: `"Tyr1258Cys"`, `"Pro1266Gln"`, `"Met1304Val"`
2. **1-letter格式**: `"Y1258C"`, `"P1266Q"`, `"M1304V"`
3. **分列格式**: `Position` + `WT_AA` + `Mut_AA`

#### 自动处理
- 自动跳过标题行 (如 `"Table S3: VWD type 2B missense mutations"`)
- 自动跳过列名行 (如 `"Amino acid\nchange"`)
- 自动过滤同义突变
- 自动进行 3-letter ↔ 1-letter 转换

---

### B. AF3 JSON 格式转换

#### AlphaFold3 Server 官方输入格式
```json
[
  {
    "name": "VWF_Y1258C",
    "sequences": [
      {
        "proteinChain": {
          "sequence": "MIPARF...Y...C...",
          "count": 1
        }
      }
    ]
  }
]
```

**关键要点**:
- 最外层是 **JSON 数组**
- 每个任务对象包含 `name` 和 `sequences`
- `proteinChain` 包含突变后的 `sequence` (完整 2813 aa)
- 每文件最多 100 个任务 (Server 限制)

---

### C. 使用示例

#### 命令行
```bash
cd Proteo-Structure-Pipeline/src
conda activate alphafold

# 标准化表格并保存
python table_normalizer.py /path/to/2B_variants.xlsx \
    --output ../data/normalized_2B.csv

# 生成 AF3 JSON (排除已完成的)
python table_normalizer.py /path/to/2B_variants.xlsx \
    --af3-json ../data/af3_batch_2B.json \
    --wt-fasta ../structures/wt/VWF_P04275_WT.fasta \
    --exclude L1460F A1461V
```

#### Python API
```python
from table_normalizer import TableNormalizer

# 标准化
normalizer = TableNormalizer("2B_variants.xlsx")
df = normalizer.normalize()
errors = normalizer.validate()
normalizer.save_csv("output.csv")

# 转换为 AF3 JSON
jobs = normalizer.to_af3_json(
    wt_fasta="VWF_P04275_WT.fasta",
    output_path="af3_batch.json",
    exclude_existing=["L1460F", "A1461V"]
)
```

---

### D. 转换脚本对照表

| 输入格式 | 脚本 | 输出 |
|---------|------|------|
| Excel (3-letter) | `table_normalizer.py` | 标准 CSV + AF3 JSON |
| Excel (1-letter) | `table_normalizer.py` | 标准 CSV + AF3 JSON |
| CSV (分列) | `table_normalizer.py` | 标准 CSV + AF3 JSON |
| Parquet 矩阵 | `agentic_vwf_classifier.py` | 分类结果 |


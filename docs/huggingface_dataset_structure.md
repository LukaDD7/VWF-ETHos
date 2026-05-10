# HuggingFace Dataset 结构说明

> 更新时间：2026-05-10
> 仓库：https://huggingface.co/datasets/lucachangretta/VWF

---

## 数据集组织结构

```
lucachangretta/VWF/
├── README.md                          # 旧的 A1+GPIbα 结果说明
├── analysis/
│   └── iptm_results.csv               # A1+GPIbα iPTM 汇总（Git 已提交）
├── boltz2_a1_gpiba_results/          # A1+GPIbα 结构预测结果（~842MB）
│   ├── boltz_results_VWF_WT_vs_GPIb_alpha/
│   │   └── predictions/VWF_WT_vs_GPIb_alpha/
│   │       ├── VWF_WT_vs_GPIb_alpha_model_*.cif    (5个样本)
│   │       └── confidence_VWF_WT_vs_GPIb_alpha_model_*.json
│   └── ... (74 个 job)
│
├── vwd_functional_panel/             # VWD Functional Panel 结果（新增，~2.9GB）
│   └── boltz_results/
│       └── boltz_results_yamls/
│           ├── predictions/
│           │   └── (结构文件和 confidence 文件)
│           ├── lightning_logs/
│           ├── msa/
│           └── processed/
│
└── vwd_analysis/                     # VWD 分析文件（新增）
    ├── diagnostic_panel.csv           # 临床诊断面板汇总
    ├── job_manifest.csv               # 991 个 job 的元数据清单
    └── summary.json                   # 运行摘要

```

---

## 两批数据的区别

| 属性 | A1+GPIbα | VWD Functional Panel |
|------|----------|----------------------|
| **job 数** | 74 | ~991 |
| **蛋白** | VWF A1 domain + GPIbα (489aa) | VWF D'-D3 + various functional constructs (~250-400aa) |
| **samples/job** | 5 | 5 |
| **recycling_steps** | 3 | 3 |
| **diffusion_samples** | 5 | 8 |
| **总大小** | ~842MB | ~2.9GB |
| **用途** | Type 2B / 2M LOF 分类校准 | VWD functional assay 结构基础 |

---

## 在另一台设备拉取

### 完整仓库

```python
from datasets import load_dataset

# 完整数据集（约 3.7GB）
ds = load_dataset("lucachangretta/VWF", split="train")

# 只看 VWD functional panel 结果
ds_vwd = load_dataset("lucachangretta/VWF", split="train", data_dir="vwd_functional_panel")
```

### 只拉取 CSV 分析文件

```python
import pandas as pd

# A1+GPIbα iPTM 结果
df_iptm = pd.read_csv("https://huggingface.co/datasets/lucachangretta/VWF/raw/main/analysis/iptm_results.csv")

# VWD diagnostic panel
df_diag = pd.read_csv("https://huggingface.co/datasets/lucachangretta/VWF/raw/main/vwd_analysis/diagnostic_panel.csv")

# VWD job manifest
df_manifest = pd.read_csv("https://huggingface.co/datasets/lucachangretta/VWF/raw/main/vwd_analysis/job_manifest.csv")
```

---

## 上传脚本

| 脚本 | 用途 |
|------|------|
| `scripts/pipeline/upload_boltz2_results_huggingface.py` | 上传 A1+GPIbα 结果 |
| `scripts/pipeline/upload_vwd_functional_boltz2_results_huggingface.py` | 上传 VWD Functional Panel 结果（增量） |

增量上传说明：HuggingFace `upload_folder` 会跳过已存在的文件，已上传的 A1+GPIbα 结果不会被重复上传。

---

## 注意事项

1. **数据分离**：两批结果放在不同的子目录下，不会混在一起
2. **增量传输**：后续新增的 VWD 结果可以直接再次运行上传脚本，已存在的文件会自动跳过
3. **git push 还是 HF**：小的分析 CSV 建议 git push（几 KB），大的结构文件走 HuggingFace

---

## 更新日志

| 日期 | 操作 |
|------|------|
| 2026-05-08 | 上传 A1+GPIbα 结果（842MB CIF + iptm_results.csv） |
| 2026-05-10 | 上传 VWD Functional Panel 结果（~2.9GB）+ 分析文件 |
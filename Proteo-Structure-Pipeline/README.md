# Proteo-Structure Pipeline

VWF (von Willebrand Factor) 蛋白质3D结构分析管线，基于AlphaFold3 Server。

## 项目结构

```
Proteo-Structure-Pipeline/
├── src/                          # 源代码
│   ├── phase1_smart_filter.py    # Phase 1: 过滤和解析
│   ├── phase2_af3_batch_generator.py  # Phase 2: 生成AF3批量JSON
│   └── phase3_structural_scoring.py   # Phase 3: 结构比对评分
├── data/                         # 数据目录
├── structures/                   # 结构文件
│   ├── wt/                       # 野生型FASTA
│   ├── mutant/                   # 突变型FASTA (Phase 1生成)
│   └── predictions/              # AlphaFold3预测结果(CIF)
├── output/                       # 输出结果
│   ├── phase1_filtered_mutations.csv
│   ├── af3_batches/              # AF3批量JSON文件
│   │   ├── af3_batch_01.json
│   │   ├── af3_batch_02.json
│   │   └── batch_index.json
│   └── phase3_final_scores.csv
└── logs/                         # 日志文件
```

## 环境配置

使用专门的 `alphafold` conda 环境（已预装编译工具和依赖）：

```bash
# 创建环境（仅需一次）
conda create -n alphafold python=3.10 gxx_linux-64=11.* gcc_linux-64=11.* \
    gfortran_linux-64=11.* libstdcxx-ng libgcc-ng numpy pandas openpyxl \
    scipy biopython matplotlib seaborn -c conda-forge -y

# 激活环境
conda activate alphafold
```

## 使用说明

### Phase 1: Smart Filter & Parsing

从Excel文件中提取missense变异，生成野生型和突变型FASTA文件。

```bash
conda activate alphafold
cd src
python phase1_smart_filter.py \
    --excel /path/to/Final_VWF_Target_List_with_AlphaGenome.xlsx \
    --output-dir ../output
```

**输出:**
- `output/phase1_filtered_mutations.csv` - 过滤后的变异列表 (1,305个)
- `structures/wt/VWF_P04275_WT.fasta` - 野生型FASTA (2,813 aa)
- `structures/mutant/VWF_*.fasta` - 突变型FASTA文件 (1,305个)

---

### Phase 2: AlphaFold3 Batch JSON 生成器

**重要:** AlphaFold3 Server 目前没有开放的 REST API。Phase 2 生成批量 JSON 文件，需人工上传到官网。

```bash
# 测试模式 - 只生成前5个变异的JSON
python phase2_af3_batch_generator.py \
    --csv ../output/phase1_filtered_mutations.csv \
    --wt-fasta ../structures/wt/VWF_P04275_WT.fasta \
    --output-dir ../output/af3_batches \
    --chunk-size 20 \
    --limit 5

# 生产模式 - 生成全部1,305个变异 (自动分卷)
python phase2_af3_batch_generator.py \
    --csv ../output/phase1_filtered_mutations.csv \
    --wt-fasta ../structures/wt/VWF_P04275_WT.fasta \
    --output-dir ../output/af3_batches \
    --chunk-size 20
```

**关键参数:**
- `--chunk-size N`: 每个JSON文件包含的任务数（默认20，建议不超过50）
- `--limit N`: 限制只处理前N个变异（**测试用**）

**生成的JSON格式:**
```json
[
  {
    "name": "VWF_WT",
    "modelSeeds": [],
    "sequences": [
      {
        "proteinChain": {
          "sequence": "MIPARF...",
          "count": 1
        }
      }
    ]
  },
  {
    "name": "VWF_G1531D",
    "modelSeeds": [],
    "sequences": [
      {
        "proteinChain": {
          "sequence": "MIPARF...D...",
          "count": 1
        }
      }
    ]
  }
]
```

**输出:**
- `output/af3_batches/af3_batch_01.json` - 批次JSON文件
- `output/af3_batches/batch_index.json` - 批次索引

---

### 【人工操作】上传 AlphaFold3 Server

Phase 2 完成后，需要人工提交到 AlphaFold3 Server：

```
1. 访问 https://alphafoldserver.com/
2. 使用 Google 账号登录
3. 点击 "Batch Submission"
4. 依次上传 output/af3_batches/af3_batch_*.json 文件
   （共约 66 个文件，每文件20个任务，总计 1,306 个任务含WT）
5. 等待服务器计算完成（预计 1-3 天）
6. 下载结果压缩包，解压到 structures/predictions/
   确保包含: VWF_WT.cif, VWF_G1531D.cif, VWF_T1538K.cif, ...
```

---

### Phase 3: 3D 空间结构比对与评分引擎

使用BioPython进行结构比对、RMSD计算、pLDDT分析。

```bash
# 先生成 Phase 2 的结果CSV（用于Phase 3输入）
python -c "
import pandas as pd
import json
from pathlib import Path

# 读取 batch_index.json 生成 Phase 3 输入CSV
with open('../output/af3_batches/batch_index.json') as f:
    index = json.load(f)

jobs = []
for job in index['jobs_summary']:
    if job['name'] != 'VWF_WT':
        parts = job['name'].split('_')
        wt_aa = parts[1][0] if len(parts) > 1 else ''
        pos = parts[1][1:-1] if len(parts) > 1 else ''
        mut_aa = parts[1][-1] if len(parts) > 1 else ''
        jobs.append({
            'name': job['name'],
            'type': 'Mutant' if job['name'] != 'VWF_WT' else 'WT',
            'aa_position': pos,
            'wt_aa': wt_aa,
            'mut_aa': mut_aa,
            'output_path': f'../structures/predictions/{job["name"]}.cif'
        })

pd.DataFrame(jobs).to_csv('../output/phase2_prediction_results.csv', index=False)
print(f'Generated phase2_prediction_results.csv with {len(jobs)} jobs')
"

# 运行 Phase 3
python phase3_structural_scoring.py \
    --wt-cif ../structures/predictions/VWF_WT.cif \
    --csv ../output/phase2_prediction_results.csv \
    --output ../output/phase3_final_scores.csv
```

**计算指标:**

| 指标 | 描述 |
|------|------|
| `Global_RMSD` | 全蛋白Cα叠合RMSD (2813 aa) |
| `Local_RMSD_10A` | 突变位点10Å半径内原子RMSD |
| `pLDDT_Delta` | 突变位点pLDDT变化 (Mut - WT) |

**输出:**
- `output/phase3_final_scores.csv` - 最终评分表

---

## 完整流程示例

```bash
# 1. 激活环境
conda activate alphafold
cd src

# 2. Phase 1: 数据准备
python phase1_smart_filter.py \
    --excel /path/to/Final_VWF_Target_List_with_AlphaGenome.xlsx

# 3. Phase 2: 生成批量JSON（先测试5个）
python phase2_af3_batch_generator.py --limit 5 --chunk-size 20

# 【人工操作】上传 JSON 到 https://alphafoldserver.com/
# 等待计算完成，下载结果到 structures/predictions/

# 4. Phase 3: 结构评分
python phase3_structural_scoring.py
```

---

## AlphaFold3 批量提交流程图

```
Phase 1                    Phase 2                       【人工】                    Phase 3
  │                          │                              │                         │
  ▼                          ▼                              ▼                         ▼
Excel ──► CSV/FASTA ──► JSON Batches ──► AlphaFold3 Server ──► CIF Files ──► RMSD/pLDDT
(2431)     (1305)       (af3_batch_*.json)  (alphafoldserver.com)  (VWF_*.cif)   Scores
                         自动分卷66个文件       1-3天计算                2813 aa      CSV
```

---

## 依赖

已在conda环境中预装：
- Python 3.10
- GCC/G++ 11.4.0
- NumPy 2.2.6, Pandas 2.3.3, SciPy 1.15.2
- BioPython 1.86
- Matplotlib 3.10.8, Seaborn

## 参考

- VWF Uniprot: [P04275](https://www.uniprot.org/uniprotkb/P04275)
- AlphaFold3 Server: [alphafoldserver.com](https://alphafoldserver.com/)
- BioPython PDB: [Documentation](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ)

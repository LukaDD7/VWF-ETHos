# VWF Type-2 AlphaFold3 参考表格说明

## 概述

本表格包含了所有Type-2分型VWF变异的详细信息，用于与AlphaFold3结构预测结果对应。

## 文件说明

- **VWF_Type2_AF3_Reference_Table.xlsx** - 主表格（Excel格式，推荐）
- **VWF_Type2_AF3_Reference_Table.csv** - 主表格（CSV格式）
- **VWF_Type2_in_AF3_batches.csv** - 原始Type-2变异数据

## 表格字段说明

### 基本信息

| 字段 | 说明 | 示例 |
|------|------|------|
| Batch_File | 所在JSON batch文件 | af3_type2_batch_01 |
| Job_Name | AF3任务名称 | VWF_S2775C |
| AA_Change | 氨基酸变化 | S2775C |
| WT_AA | 野生型氨基酸 | S |
| Position | 突变位置 | 2775 |
| Mut_AA | 突变型氨基酸 | C |

### 分型与分类

| 字段 | 说明 | 示例 |
|------|------|------|
| Type2_Subtype | VWD Type-2分型 | type2A / type2B / type2M / type2N |
| ACMG_Classification | ACMG致病性分类 | Pathogenic / Likely pathogenic / Uncertain significance |
| Clinical_Classification | 临床分型 | type2A / type2B / type2M / type2N |

### 分子信息

| 字段 | 说明 | 示例 |
|------|------|------|
| Transcript_Consequence | cDNA变化 | c.8324C>G |
| Protein_Consequence | 蛋白变化 | p.Ser2775Cys |
| RSID | dbSNP ID | rs267607349 |
| Variant_Type | 变异类型 | missense |
| Allele_Frequency | 等位基因频率 | 6.19482e-07 |
| VWF_Domain | VWF蛋白功能域 | A1 / A2 / A3 / D4等 |

### AlphaGenome信息

| 字段 | 说明 | 示例 |
|------|------|------|
| AlphaGenome_Max_Score | AlphaGenome最大得分 | 4.0 |
| Splice_Delta | 剪接位点变化 | 0.0002 |

### 技术信息

| 字段 | 说明 | 示例 |
|------|------|------|
| Sequence_Length | 序列长度（氨基酸） | 2813 |
| Sequence_MD5 | 序列校验值 | 674896408 |
| Notes | 备注 | Type 2 VWD - type2A |

## Type-2分型说明

### type2A (42个)
- **特点**: 高分子量多聚体缺失，VWF抗原(VWF:Ag)减少，VWF活性(VWF:RCo)更显著降低
- **常见位置**: A2域（ADAMTS13切割位点）、A1域
- **代表变异**: VWF_R1597W, VWF_G1573S

### type2B (12个)
- **特点**: 与GP1b亲和力增加，血小板减少
- **常见位置**: A1域（GP1b结合位点）
- **代表变异**: VWF_A1461V, VWF_R1379C

### type2M (25个)
- **特点**: 多聚体正常，但与血小板或胶原结合缺陷
- **常见位置**: A1域、A3域
- **代表变异**: VWF_G1324S, VWF_I1425F

### type2N (20个)
- **特点**: VWF与因子VIII结合缺陷，FVIII水平降低
- **常见位置**: D'域、D3域
- **代表变异**: VWF_C788Y, VWF_R816Q

## VWF功能域分布

| 功能域 | 位置 | 功能 | 样本数 |
|--------|------|------|--------|
| D1-D2 | 1-272 | 信号肽+前肽 | 7 |
| D'D3 | 273-764 | GP1b结合 | 6 |
| A1 | 765-1243 | 胶原结合，type2B突变热点 | 28 |
| A2 | 1244-1481 | ADAMTS13切割，type2A热点 | 36 |
| A3 | 1482-1672 | 胶原结合，type2M热点 | 15 |
| D4 | 1673-1872 | 多聚化 | 5 |
| CT | 1873-2813 | C-末端 | 2 |

## 使用指南

### 1. 查看特定分型的样本
```python
import pandas as pd
df = pd.read_excel("VWF_Type2_AF3_Reference_Table.xlsx")

# 查看所有type2A
type2a = df[df['Type2_Subtype'] == 'type2A']

# 查看Pathogenic变异
pathogenic = df[df['ACMG_Classification'] == 'Pathogenic']
```

### 2. 根据AlphaFold结果匹配信息
```python
# 假设af_result是AlphaFold输出文件名（如"fold_vwf_s2775c_model_0.cif"）
af_name = "VWF_S2775C"
info = df[df['Job_Name'] == af_name].iloc[0]
print(f"分型: {info['Type2_Subtype']}")
print(f"ACMG: {info['ACMG_Classification']}")
print(f"功能域: {info['VWF_Domain']}")
```

### 3. 按Batch查看
```python
# 查看batch_01的所有样本
batch1 = df[df['Batch_File'] == 'af3_type2_batch_01']
```

## 统计摘要

- **总样本数**: 100 (1 WT + 99 Type-2变异)
- **Batch文件**: 4个 (batch_01~04)
- **Type-2A**: 42个
- **Type-2B**: 12个
- **Type-2M**: 25个
- **Type-2N**: 20个
- **Pathogenic**: 34个
- **Likely pathogenic**: 28个
- **Uncertain significance**: 37个

## 注意事项

1. VWF_WT为野生型对照，位于batch_01
2. 所有序列长度均为2813个氨基酸
3. Sequence_MD5用于验证序列完整性
4. AlphaGenome数据来自内皮细胞模型(CL:0000115)

## 关联文件

- `af3_batches_type2/` - Type-2专用batch JSON文件
- `af3_batches/` - 主batch JSON文件（821个非Type-2样本）
- `VWF_Type2_variants.csv` - 原始Type-2变异数据

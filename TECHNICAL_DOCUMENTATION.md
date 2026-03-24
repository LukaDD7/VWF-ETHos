# VWF-ETHos 技术文档

## 目录
1. [项目架构](#1-项目架构)
2. [核心代码详解](#2-核心代码详解)
3. [数据处理流程](#3-数据处理流程)
4. [机器学习实现](#4-机器学习实现)
5. [分析模块说明](#5-分析模块说明)
6. [输出结果解读](#6-输出结果解读)
7. [扩展开发指南](#7-扩展开发指南)

---

## 1. 项目架构

### 1.1 目录结构

```
VWF_ErTongYiyuan/
├── vwf_type2_analysis.py          # Type-2分型预测主程序
├── scripts/                        # AlphaGenome分析脚本
│   ├── 01_filter_target_vus.py    # 变异过滤
│   ├── 02_preprocess_and_liftover.py  # 坐标转换
│   ├── 03_run_alphagenome_inference.py # AlphaGenome API调用
│   ├── 04_analyze_and_visualize.py     # 结果可视化
│   ├── 07_vcf_batch_alphagenome.py     # VCF批量处理
│   ├── 07e_god_mode_epigenome_crawler.py # 全模态爬取
│   └── 08_merge_godmode_results.py      # 结果合并
├── Proteo-Structure-Pipeline/     # AF3结构预测Pipeline
│   └── src/
│       ├── phase1_smart_filter.py      # 智能过滤
│       ├── phase2_af3_batch_generator.py # AF3批次生成
│       └── phase3_structural_scoring.py  # 结构评分
├── results/                        # 分析结果
│   └── 07_VCF_AlphaGenome_Results.csv  # AlphaGenome完整结果
├── figures/                        # 可视化图表
├── Type-2-Analysis-README.md      # Type-2分析详细文档
└── AlphaGenome_Type2_Analysis_Report.md  # AlphaGenome分析报告
```

### 1.2 数据流

```
输入数据
    ├── ClinVar/HGMD变异数据 (Excel)
    ├── VCF文件
    └── FASTA序列
         ↓
    数据预处理 (scripts/01-02)
         ↓
    AlphaGenome API调用 (scripts/03)
         ↓
    结构预测 (AF3)
         ↓
    特征提取 (vwf_type2_analysis.py)
         ↓
    ML分型预测
         ↓
    输出结果 (CSV/图表/报告)
```

---

## 2. 核心代码详解

### 2.1 Type-2分型预测主程序

**文件**: `vwf_type2_analysis.py`

#### 2.1.1 功能域映射函数

```python
def get_vwf_domain(position):
    """
    根据氨基酸位置确定VWF功能域

    参数:
        position: int, 氨基酸位置 (1-2813)

    返回:
        str: 功能域名称

    生物学基础:
        基于UniProt P04275的VWF标准注释
        各功能域负责不同的生物学功能
    """
    if 1 <= position <= 272:
        return 'D1-D2'      # 信号肽 + 前肽
    elif 273 <= position <= 761:
        return "D'D3"       # GP1b结合
    elif 762 <= position <= 1035:
        return 'A1'         # 胶原结合 (type 2B)
    elif 1036 <= position <= 1227:
        return 'A2'         # ADAMTS13切割 (type 2A)
    elif 1228 <= position <= 1451:
        return 'A3'         # 胶原结合 (type 2M)
    elif 1452 <= position <= 1670:
        return 'D4'         # 多聚化
    elif 1671 <= position <= 2051:
        return 'C1-C2'      # 连接区
    elif 2052 <= position <= 2813:
        return 'CT'         # C-末端
    else:
        return 'Unknown'
```

**生物学意义**:
- **D1-D2**: 前肽区域，影响VWF多聚化
- **D'D3**: 血小板GP1b受体结合位点
- **A1**: 胶原结合，type 2B突变热点
- **A2**: ADAMTS13蛋白酶切割位点，type 2A热点
- **A3**: 胶原结合，type 2M热点

#### 2.1.2 分型预测函数

```python
def predict_subtype_by_rules(position, domain):
    """
    基于规则的Type-2分型预测

    参数:
        position: int, 变异位置
        domain: str, 功能域名称 (可选)

    返回:
        dict: 包含预测分型、置信度、机制解释

    规则来源:
        基于59个样本的统计分析
    """
    if domain is None:
        domain = get_vwf_domain(position)

    # 规则1: A3域 → type2A (100%准确率)
    if domain == 'A3':
        return {
            'subtype': 'type2A',
            'confidence': 1.00,
            'mechanism': '胶原结合缺陷 + ADAMTS13易感性',
            'recommendation': 'RIPA实验验证'
        }

    # 规则2: A2域 → type2M (57%) 或 type2A (43%)
    elif domain == 'A2':
        return {
            'subtype': 'type2M',
            'confidence': 0.57,
            'mechanism': 'ADAMTS13切割位点改变',
            'alternative': 'type2A',
            'recommendation': 'VWF多聚体分析 + RIPA'
        }

    # 规则3: A1域 → type2N (63%) 或 type2A (37%)
    elif domain == 'A1':
        return {
            'subtype': 'type2N',
            'confidence': 0.63,
            'mechanism': 'GP1b/FVIII结合异常',
            'alternative': 'type2A',
            'recommendation': '瑞斯托霉素辅因子实验'
        }

    # 规则4: D4域 → type2M (80%)
    elif domain == 'D4':
        return {
            'subtype': 'type2M',
            'confidence': 0.80,
            'mechanism': '多聚化异常',
            'recommendation': 'VWF多聚体电泳'
        }

    # 默认: 其他域 → type2A
    else:
        return {
            'subtype': 'type2A',
            'confidence': 0.70,
            'mechanism': '分泌/功能异常',
            'recommendation': '全面VWF功能检测'
        }
```

#### 2.1.3 批量预测功能

```python
def batch_predict(input_file, output_file='predictions.csv'):
    """
    批量预测变异分型

    参数:
        input_file: str, 输入CSV文件路径
                    必需列: Position
                    可选列: VWF_Domain
        output_file: str, 输出CSV文件路径

    返回:
        DataFrame: 包含预测结果的DataFrame

    示例输入:
        Position,VWF_Domain
        1437,A3
        1100,A2
        800,A1
    """
    df = pd.read_csv(input_file)

    predictions = []
    for _, row in df.iterrows():
        pos = row['Position']
        domain = row.get('VWF_Domain', None)
        pred = predict_single(pos, domain)
        predictions.append(pred)

    # 合并预测结果
    pred_df = pd.DataFrame(predictions)
    result = pd.concat([df, pred_df], axis=1)
    result.to_csv(output_file, index=False)

    return result
```

### 2.2 AlphaGenome分析脚本

#### 2.2.1 变异过滤 (01_filter_target_vus.py)

```python
class Phase1Filter:
    """
    Phase 1: 变异过滤与目标选择

    功能:
        1. 从ClinVar/HGMD提取VWF变异
        2. 过滤出VUS和致病变异
        3. 排除剪接破坏者
        4. 解析氨基酸改变
    """

    def filter_missense_variants(self, df):
        """
        过滤错义突变

        逻辑:
            保留Molecular.consequence为"missense variant"的记录
        """
        return df[df['INFO_consequences_base'].str.contains(
            'missense_variant', case=False, na=False
        )]

    def filter_splice_disruptors(self, df, threshold=0.3):
        """
        过滤剪接破坏者

        参数:
            threshold: float, 剪接概率差异阈值

        逻辑:
            如果 |Splice_REF - Splice_ALT| > threshold，排除该变异
        """
        splice_diff = abs(df['Splice_REF'] - df['Splice_ALT'])
        return df[splice_diff <= threshold]

    def parse_amino_acid_change(self, change_str):
        """
        解析氨基酸改变字符串

        支持格式:
            - 1字母: "G1531D"
            - 3字母: "Val1409Phe"
            - NM命名: "p.Gly1531Asp"

        返回:
            dict: {wt_aa, position, mut_aa}
        """
        # 解析逻辑...
        pass
```

#### 2.2.2 AlphaGenome API调用 (03_run_alphagenome_inference.py)

```python
class AlphaGenomeInference:
    """
    AlphaGenome API调用器

    功能:
        批量调用AlphaGenome API进行变异效应预测

    支持的模态:
        - RNA-seq (基因表达)
        - Splice (剪接)
        - CAGE (启动子活性)
        - DNase (染色质可及性)
        - ChIP (组蛋白修饰、TF结合)
    """

    def predict_variant(self, variant):
        """
        预测单个变异

        参数:
            variant: dict, 包含chrom, pos, ref, alt

        返回:
            dict: 各模态的预测结果

        API调用示例:
            POST /api/v1/predict
            Body: {
                "chrom": "chr12",
                "pos": 6126975,
                "ref": "T",
                "alt": "C",
                "tracks": ["RNA-seq", "Splice", "DNase"]
            }
        """
        # API调用逻辑...
        pass

    def calculate_delta_score(self, ref_pred, alt_pred):
        """
        计算Delta分数

        公式:
            Delta = |ALT - REF| / max(|REF|, |ALT|)

        生物学意义:
            反映变异对功能的影响程度
        """
        delta = abs(alt_pred - ref_pred)
        return delta / max(abs(ref_pred), abs(alt_pred), 1e-6)
```

#### 2.2.3 God Mode全模态爬取 (07e_god_mode_epigenome_crawler.py)

```python
class GodModeCrawler:
    """
    God Mode: 全模态表观基因组爬取

    功能:
        获取所有可用细胞类型的全模态数据

    模态列表:
        1. RNA-seq
        2. Splice Sites
        3. Splice Site Usage
        4. Splice Junctions
        5. DNase
        6. ATAC
        7. ChIP-Histone
        8. ChIP-TF
        9. CAGE
        10. PRO-cap
        11. Contact Maps
    """

    ONTOLOGY_TERMS = ['CL:0000115']  # 内皮细胞

    def crawl_all_modalities(self, variant_list):
        """
        爬取所有模态数据

        参数:
            variant_list: list, 变异列表

        返回:
            DataFrame: 每个变异的11模态结果

        注意:
            每个变异需要~2分钟，1000个变异需要~33小时
        """
        results = []
        for variant in variant_list:
            result = self.predict_all_modalities(variant)
            results.append(result)
        return pd.DataFrame(results)

    def predict_all_modalities(self, variant):
        """
        预测单个变异的所有模态

        返回格式:
            {
                'RNA-seq': float,
                'Splice': float,
                'DNase': float,
                ...
            }
        """
        # 并行调用所有模态API...
        pass
```

### 2.3 AF3结构分析Pipeline

#### 2.3.1 智能过滤 (phase1_smart_filter.py)

```python
class Phase1Filter:
    """
    Phase 1: 智能过滤变异

    逻辑流程:
        1. 加载ClinVar/HGMD数据
        2. 过滤错义突变
        3. 排除剪接破坏者
        4. 去重 (基于Protein.change)
        5. 生成WT和突变FASTA
    """

    def __init__(self, excel_path):
        self.excel_path = excel_path
        self.wt_fasta = None
        self.mutant_fastas = []

    def download_wt_fasta(self):
        """
        下载WT FASTA序列

        来源: UniProt P04275
        长度: 2813个氨基酸
        """
        url = "https://rest.uniprot.org/uniprotkb/P04275.fasta"
        # 下载逻辑...
        pass

    def generate_mutant_fasta(self, aa_change):
        """
        生成突变FASTA

        参数:
            aa_change: str, 如 "A1437T"

        逻辑:
            1. 解析位置 (1437)
            2. 验证WT氨基酸 (应为A)
            3. 替换为突变氨基酸 (T)
            4. 生成新FASTA序列
        """
        # 序列操作逻辑...
        pass
```

#### 2.3.2 AF3批次生成 (phase2_af3_batch_generator.py)

```python
class AF3BatchGenerator:
    """
    Phase 2: AF3批次文件生成

    约束:
        - AF3 Server每次只读取第一个任务
        - 因此每个变异生成独立JSON文件
        - ZIP包大小限制 (每日配额)
    """

    def create_batch_jobs(self, variants):
        """
        创建批次任务

        参数:
            variants: list, 变异列表

        返回:
            list: AF3Job对象列表
        """
        jobs = []
        for var in variants:
            job = AF3BatchJob(
                name=f"VWF_{var['aa_change']}",
                sequence=var['mutant_sequence']
            )
            jobs.append(job)
        return jobs

    def save_individual_jsons(self, jobs, output_dir):
        """
        保存独立JSON文件

        JSON格式:
            {
                "name": "VWF_G1531D",
                "modelSeeds": [],
                "sequences": [{
                    "proteinChain": {
                        "sequence": "MIPARF...",
                        "count": 1
                    }
                }]
            }
        """
        for job in jobs:
            json_path = os.path.join(output_dir, f"{job.name}.json")
            with open(json_path, 'w') as f:
                json.dump(job.to_dict(), f)

    def create_chunked_zip_packages(self, json_dir, chunk_size=10):
        """
        创建分块ZIP包

        参数:
            chunk_size: int, 每个ZIP包含的JSON数量

        说明:
            默认10个/包，匹配AF3 Server每日配额
        """
        json_files = glob(os.path.join(json_dir, '*.json'))
        for i in range(0, len(json_files), chunk_size):
            chunk = json_files[i:i+chunk_size]
            zip_path = f"af3_upload_part_{i//chunk_size:03d}.zip"
            # 创建ZIP逻辑...
            pass
```

#### 2.3.3 结构评分 (phase3_structural_scoring.py)

```python
class StructureAnalyzer:
    """
    Phase 3: 结构评分与分析

    功能:
        1. 加载WT结构
        2. 叠加变异结构
        3. 计算RMSD
        4. 提取pLDDT特征
    """

    def superimpose_and_calculate_rmsd(self, wt_structure, mut_structure):
        """
        结构叠加并计算RMSD

        算法:
            Kabsch算法 (最小二乘对齐)

        返回:
            dict: {
                'Global_RMSD': float,
                'Local_RMSD_10A': float
            }
        """
        # 提取CA原子
        wt_ca = self._get_ca_atoms(wt_structure)
        mut_ca = self._get_ca_atoms(mut_structure)

        # Kabsch对齐
        rotation, translation = self._kabsch_algorithm(wt_ca, mut_ca)
        mut_aligned = self._apply_transform(mut_ca, rotation, translation)

        # 计算RMSD
        global_rmsd = np.sqrt(np.mean((wt_ca - mut_aligned)**2))

        # 局部RMSD (突变位点10Å半径内)
        local_indices = self._get_atoms_in_radius(mut_ca, mut_position, 10)
        local_rmsd = np.sqrt(np.mean((wt_ca[local_indices] - mut_aligned[local_indices])**2))

        return {
            'Global_RMSD': global_rmsd,
            'Local_RMSD_10A': local_rmsd
        }

    def extract_plddt_features(self, full_data_json):
        """
        从full_data.json提取pLDDT特征

        参数:
            full_data_json: str, full_data文件路径

        返回:
            dict: {
                'mean_plddt': float,
                'min_plddt': float,
                'high_confidence_pct': float
            }

        注意:
            full_data.json包含:
                - atom_plddts: 原子级pLDDT
                - pae: PAE矩阵
                - contact_probs: 接触概率
        """
        with open(full_data_json) as f:
            data = json.load(f)

        atom_plddts = data['atom_plddts']

        return {
            'mean_plddt': np.mean(atom_plddts),
            'min_plddt': np.min(atom_plddts),
            'high_confidence_pct': sum(1 for p in atom_plddts if p > 70) / len(atom_plddts) * 100
        }
```

---

## 3. 数据处理流程

### 3.1 完整处理流程

```
Step 1: 数据输入
    ├── ClinVar_HGMD_merge_annotated.xlsx
    ├── VCF文件
    └── FASTA序列
         ↓
Step 2: 数据清洗
    ├── 过滤错义突变 (filter_missense_variants)
    ├── 排除剪接破坏者 (filter_splice_disruptors)
    └── 去重 (基于Protein.change)
         ↓
Step 3: AlphaGenome分析
    ├── API调用 (predict_variant)
    ├── 计算Delta分数 (calculate_delta_score)
    └── 多模态整合 (11 modalities)
         ↓
Step 4: AF3结构预测
    ├── 生成突变序列 (generate_mutant_fasta)
    ├── 创建批次任务 (create_batch_jobs)
    └── 提交AF3 Server
         ↓
Step 5: 特征提取
    ├── PAE矩阵分析
    ├── pLDDT提取
    └── RMSD计算
         ↓
Step 6: 机器学习
    ├── 特征工程 (Position, Domain, PAE features)
    ├── 模型训练 (Gradient Boosting)
    └── 分型预测
         ↓
Step 7: 结果输出
    ├── CSV表格
    ├── 可视化图表
    └── 诊断报告
```

### 3.2 关键数据结构

#### 3.2.1 变异记录

```python
{
    'chrom': 'chr12',
    'pos': 6126975,
    'ref': 'T',
    'alt': 'C',
    'aa_change': 'A1437T',
    'position': 1437,
    'wt_aa': 'A',
    'mut_aa': 'T',
    'domain': 'A3',
    'type2_subtype': 'type2A',
    'acmg': 'Pathogenic'
}
```

#### 3.2.2 AlphaGenome结果

```python
{
    'Delta_Max_Core': 4.0,
    'AG_RNA_SEQ': 4.0,
    'AG_SPLICE_SITES': 0.015625,
    'AG_CAGE': 0.0,
    'AG_DNASE': 0.625,
    'AG_CHIP_TF': 48.0,
    # ... 11 modalities total
}
```

#### 3.2.3 AF3结构结果

```python
{
    'pae': [[0.8, 5.2, ...], ...],  # 2813x2813矩阵
    'atom_plddts': [95.2, 87.3, ...],  # 21506个原子
    'contact_probs': [0.9, 0.3, ...],  # 2813个残基
    'token_res_ids': [1, 2, ..., 2813]  # 残基编号
}
```

---

## 4. 机器学习实现

### 4.1 特征工程

```python
# 1. 位置特征
X['Position_Norm'] = position / 2813.0  # 归一化

# 2. 功能域编码
X['Domain_Encoded'] = LabelEncoder().fit_transform(domain)

# 3. PAE特征
X['PAE_Delta_Local'] = mut_local_flex - wt_local_flex
X['Local_Flex_Ratio'] = mut_local_flex / wt_local_flex

# 4. AlphaGenome特征
X['AlphaGenome_Score'] = alphagenome_max_score
X['Splice_Delta'] = splice_delta
```

### 4.2 模型训练

```python
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import cross_val_score

# 模型配置
model = GradientBoostingClassifier(
    n_estimators=100,
    max_depth=3,
    learning_rate=0.1,
    random_state=42
)

# 特征矩阵
feature_cols = [
    'Position_Norm',
    'Domain_Encoded',
    'PAE_Delta_Local',
    'AlphaGenome_Score',
    'Splice_Delta'
]
X = df[feature_cols]
y = df['Type2_Subtype']

# 交叉验证
cv_scores = cross_val_score(model, X, y, cv=5)
print(f"CV Accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")

# 训练完整模型
model.fit(X, y)

# 特征重要性
for feat, imp in zip(feature_cols, model.feature_importances_):
    print(f"{feat}: {imp:.3f}")
```

### 4.3 模型评估

```python
from sklearn.metrics import classification_report, confusion_matrix

# 预测
y_pred = model.predict(X_test)

# 分类报告
print(classification_report(y_test, y_pred))

# 混淆矩阵
cm = confusion_matrix(y_test, y_pred)
print(cm)
```

---

## 5. 分析模块说明

### 5.1 命令行接口

#### 单个变异预测

```bash
python vwf_type2_analysis.py --position 1437

输出:
============================================================
VWF Type-2分型预测结果
============================================================
变异位置: 1437
功能域: A3
预测分型: type2A
可信度: 100.0%
预测方法: rule_based

临床建议: ADAMTS13切割异常，建议进行RIPA实验验证
============================================================
```

#### 批量预测

```bash
python vwf_type2_analysis.py \
    --batch-predict variants.csv \
    --output predictions.csv
```

### 5.2 API接口

```python
from vwf_type2_analysis import predict_single, batch_predict

# 单个预测
result = predict_single(position=1437, domain='A3')
print(result['predicted_subtype'])  # 'type2A'

# 批量预测
results = batch_predict('input.csv', 'output.csv')
```

---

## 6. 输出结果解读

### 6.1 分型预测结果

| 字段 | 说明 | 示例 |
|------|------|------|
| `position` | 变异位置 | 1437 |
| `domain` | VWF功能域 | A3 |
| `predicted_subtype` | 预测分型 | type2A |
| `confidence` | 预测置信度 | 1.00 (100%) |
| `mechanism` | 分子机制 | 胶原结合缺陷 |
| `recommendation` | 实验建议 | RIPA实验验证 |

### 6.2 置信度解释

| 置信度 | 可信度 | 建议 |
|--------|--------|------|
| > 90% | 高 | 可直接参考 |
| 60-90% | 中 | 建议结合其他特征 |
| < 60% | 低 | 需进一步验证 |

### 6.3 各分型实验建议

| 预测分型 | 推荐实验 | 验证指标 |
|----------|---------|---------|
| type2A | RIPA | 低剂量瑞斯托霉素聚集增加 |
| type2A | VWF多聚体电泳 | 缺乏大分子量多聚体 |
| type2B | 瑞斯托霉素辅因子 | 亲和力增加 |
| type2B | 血小板计数 | 血小板减少 |
| type2M | 胶原结合实验 | 结合率降低 |
| type2N | FVIII结合实验 | 结合率降低 |

---

## 7. 扩展开发指南

### 7.1 添加新的分型规则

```python
# 在 vwf_type2_analysis.py 中

def predict_subtype_by_rules(position, domain):
    # ... 现有规则 ...

    # 添加新规则
    if domain == 'NEW_DOMAIN':
        return {
            'subtype': 'new_subtype',
            'confidence': 0.85,
            'mechanism': '新机制描述',
            'recommendation': '新实验建议'
        }
```

### 7.2 集成新的ML模型

```python
from sklearn.ensemble import RandomForestClassifier

# 定义新模型
model_new = RandomForestClassifier(n_estimators=200)
model_new.fit(X_train, y_train)

# 保存模型
import joblib
joblib.dump(model_new, 'models/rf_model.pkl')

# 加载使用
model_new = joblib.load('models/rf_model.pkl')
prediction = model_new.predict(X_new)
```

### 7.3 添加新的AlphaGenome模态

```python
# 在 03_run_alphagenome_inference.py 中

MODALITIES = [
    'RNA-seq',
    'Splice',
    # ... 现有模态 ...
    'NEW_MODALITY',  # 添加新模态
]

def predict_new_modality(self, variant):
    """新模态预测函数"""
    # 实现逻辑...
    pass
```

---

## 附录

### A. 依赖包

```
pandas>=1.3.0
numpy>=1.20.0
scikit-learn>=0.24.0
matplotlib>=3.3.0
seaborn>=0.11.0
biopython>=1.78
openpyxl>=3.0.0
requests>=2.25.0
```

### B. 环境配置

```bash
conda create -n vwf python=3.10
conda activate vwf
pip install -r requirements.txt
```

### C. 常见问题

**Q: 如何处理新变异？**
```python
# 只需位置信息即可预测
result = predict_single(position=1500)
```

**Q: 如何更新参考数据？**
```bash
# 修改Excel文件后重新加载
type2_df = pd.read_excel('VWF_Type2_AF3_Reference_Table.xlsx')
```

**Q: 如何调试预测结果？**
```python
# 启用详细日志
import logging
logging.basicConfig(level=logging.DEBUG)
```

---

**文档版本**: v1.0

**最后更新**: 2026-03-24

**维护者**: LukaDD7

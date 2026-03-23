import pickle
import numpy as np
import pandas as pd
from pathlib import Path

# ==================== 配置 ====================
# 基础目录配置
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"

PKL_FILE = RESULTS_DIR / "03_inference_results.pkl"
# 填入你要导出给医生的那个变异的 hg38 物理坐标
TARGET_POSITION = 6110876  # 替换为你关注的那个高危变异的 POS
OUTPUT_CSV = RESULTS_DIR / f"Raw_Matrix_chr12_{TARGET_POSITION}.csv"

# 切片窗口：变异点前后各 16kb (同画图逻辑)
CENTER_IDX = 524288
HALF_WINDOW = 16384
start_idx = CENTER_IDX - HALF_WINDOW
end_idx = CENTER_IDX + HALF_WINDOW

print(f"正在加载 {PKL_FILE}，这可能需要一点时间...")
with open(PKL_FILE, "rb") as f:
    results = pickle.load(f)

# 寻找目标变异
target_result = next((r for r in results if r.position == TARGET_POSITION), None)

if not target_result or target_result.raw_outputs is None:
    print(f"未找到坐标为 {TARGET_POSITION} 的变异或数据为空！")
else:
    print(f"找到目标变异！正在提取 {TARGET_POSITION} 附近的原始矩阵...")
    outputs = target_result.raw_outputs

    # 建立真实的基因组物理坐标
    positions = np.arange(TARGET_POSITION - HALF_WINDOW, TARGET_POSITION + HALF_WINDOW)

    data_dict = {"hg38_Position": positions}

    # 提取 RNA-seq (2 个 Track: +链 和 -链)
    if outputs.reference.rna_seq is not None:
        ref_rna = outputs.reference.rna_seq.values[start_idx:end_idx, :]
        alt_rna = outputs.alternate.rna_seq.values[start_idx:end_idx, :]
        data_dict["RNA_REF_Plus_Strand"] = ref_rna[:, 0]
        data_dict["RNA_REF_Minus_Strand"] = ref_rna[:, 1]
        data_dict["RNA_ALT_Plus_Strand"] = alt_rna[:, 0]
        data_dict["RNA_ALT_Minus_Strand"] = alt_rna[:, 1]
        data_dict["RNA_Delta_Mean"] = np.mean(alt_rna, axis=1) - np.mean(ref_rna, axis=1)

    # 提取 Splice Sites (4 个 Track: Donor/Acceptor × +/-链)
    if outputs.reference.splice_sites is not None:
        ref_splice = outputs.reference.splice_sites.values[start_idx:end_idx, :]
        alt_splice = outputs.alternate.splice_sites.values[start_idx:end_idx, :]
        # 索引对应关系：0:Donor(+), 1:Acceptor(+), 2:Donor(-), 3:Acceptor(-)
        data_dict["Splice_REF_Donor(+))"] = ref_splice[:, 0]
        data_dict["Splice_REF_Acceptor(+))"] = ref_splice[:, 1]
        data_dict["Splice_REF_Donor(-))"] = ref_splice[:, 2]
        data_dict["Splice_REF_Acceptor(-))"] = ref_splice[:, 3]

        data_dict["Splice_ALT_Donor(+))"] = alt_splice[:, 0]
        data_dict["Splice_ALT_Acceptor(+))"] = alt_splice[:, 1]
        data_dict["Splice_ALT_Donor(-))"] = alt_splice[:, 2]
        data_dict["Splice_ALT_Acceptor(-))"] = alt_splice[:, 3]

    # 生成数据表并导出
    df_export = pd.DataFrame(data_dict)

    # 重点标记变异所在的那一行
    df_export['Is_Mutated_Site'] = df_export['hg38_Position'] == TARGET_POSITION

    # 为了文件不要太大（3.2万行），我们可以做个小过滤：
    # 比如医生主要看 Splice，我们可以只保留有概率值的区域，但为了严谨，这里默认全量导出这 32kb。
    df_export.to_csv(OUTPUT_CSV, index=False)

    print(f"导出成功！文件已保存为：{OUTPUT_CSV}")
    print(f"该表格包含 {len(df_export)} 行数据，每一行代表变异点附近的一个真实物理碱基的各模态预测值。")

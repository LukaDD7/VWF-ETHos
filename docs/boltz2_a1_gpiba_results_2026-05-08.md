# Boltz-2 VWF A1 + GPIbα 结果解析报告

> 生成时间：2026-05-08
> 解析脚本：`scripts/pipeline/parse_a1_gpiba_results.py`

---

## 运行概况

| 参数 | 值 |
|------|-----|
| 总 job 数 | 74 |
| 完成数 | 74 |
| GPU 数 | 4 |
| 采样数/job | 5 |
| WT iPTM | 0.6776 |
| 运行时间 | ~45 分钟（含首次 JIT 编译） |

---

## 分型结果

| 预测分型 | 数量 |
|---------|------|
| 2M_LOF | 55 |
| neutral_VUS | 14 |
| 2B_GOF | 4 |
| WT_reference | 1 |

---

## 校准验证（已知分型对照）

| Variant | 已知分型 | 预测分型 | ΔiPTM | 结果 |
|---------|---------|---------|-------|------|
| VWF_R1306Q | 2B_GOF | neutral_VUS | -0.0143 | ✗ |
| VWF_R1374C | 2M_LOF | neutral_VUS | -0.0336 | ✗ |
| VWF_R1306W | 2B_GOF | 2M_LOF | -0.1498 | ✗ |
| VWF_F1369I | 2M_LOF | 2M_LOF | -0.2190 | ✓ |
| VWF_V1316M | 2B_GOF | 2M_LOF | -0.2358 | ✗ |
| VWF_I1372S | 2M_LOF | 2M_LOF | -0.2601 | ✓ |
| VWF_R1341Q | 2B_GOF | 2M_LOF | -0.3650 | ✗ |

**校准准确率：2/7 (29%)**

---

## 结论

ΔiPTM 阈值需要重新校准。当前阈值：
- GOF 阈值：+0.05（无 2B 变体达标）
- LOF 阈值：-0.05（几乎所有变体都低于此值）

实际问题可能是：
1. **2B 变体应该升高 iPTM**（GOF = 增强结合），但实际测得负 ΔiPTM
2. 这可能说明 GPIbα 复合物的 iPTM 解读与预期不同
3. 或者阈值设置方向需要反转（需要查阅 Boltz-2 iPTM 在复合物上的生物学含义）

**建议**：在另一台设备（个人电脑）拉取结果后，详细分析各个 confidence 文件中的具体指标（iptm, plddt, pae 等），重新确定阈值。

---

## 数据同步方案

### 原始结果文件

- **位置**：`output/boltz2_a1_gpiba_results/`
- **大小**：842MB
- **内容**：74 个 job 的 .cif 结构文件 + confidence_*.json

### 上传 HuggingFace

```bash
# 1. 登录
pip install huggingface_hub
huggingface-cli login

# 2. 创建数据集
huggingface-cli repo create VWF-A1-GPIb-alpha-structures --type dataset

# 3. 上传（Python 脚本）
python scripts/pipeline/upload_results_huggingface.py
```

### 另一台设备拉取

```python
from datasets import load_dataset

ds = load_dataset("LukaDD7/VWF-A1-GPIb-alpha-structures", split="train")
```

### 解析结果（CSV）

汇总文件已生成：`output/boltz2_a1_gpiba_analysis/iptm_results.csv`

直接 git push 即可（仅几 KB）：
```bash
git add output/boltz2_a1_gpiba_analysis/iptm_results.csv
git commit -m "Add Boltz-2 A1+GPIbα iPTM results"
git push
```

---

## 目录结构说明

```
output/boltz2_a1_gpiba_results/
  boltz_results_VWF_WT_vs_GPIb_alpha/
    predictions/
      VWF_WT_vs_GPIb_alpha/           ← boltz 自动创建的子目录
        VWF_WT_vs_GPIb_alpha_model_0.cif
        ...
        confidence_VWF_WT_vs_GPIb_alpha_model_0.json
        ...
    lightning_logs/
    msa/
    processed/
  boltz_results_VWF_R1306W_vs_GPIb_alpha/
    ...
```

注意：`run_a1_gpiba_boltz2.sh` 的 `.done` 标记路径有 bug（已修复），不影响实际结果。实际结果在 `boltz_results_*/predictions/*/` 下。完整 74 个 job 结果已生成，只需修复 .done 路径后重新跑即可续命（但实际不需要，所有 job 都已完成）。
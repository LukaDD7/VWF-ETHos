# VWD/VWF Functional Boltz-2 Panel 结果报告

> 生成时间：2026-05-10
> 解析脚本：`scripts/pipeline/parse_vwd_functional_boltz2_results.py`

---

## 运行概况

| 参数 | 值 |
|------|-----|
| 总 job 数 | 989 |
| confidence 文件数 | 7903 |
| 每 job 样本数 | 8 |
| iPTM 范围 | 0.0921 - 0.7153 |
| iPTM 均值 | 0.3116 |

---

## iPTM 按 assay_key 分布

| assay_key | mean iPTM | std | n |
|-----------|-----------|-----|---|
| c2_collagen_binding | 0.5589 | 0.0268 | 13 |
| c1_collagen_binding | 0.5386 | 0.0214 | 15 |
| a3_collagen_binding | 0.4943 | 0.0349 | 24 |
| a1_heparan_sulfate_binding | 0.4955 | 0.1010 | 120 |
| c4_integrin_binding | 0.4124 | 0.1193 | 14 |
| a1_gpiba_forced_binding | 0.3239 | 0.0631 | 120 |
| vwf73_adamts13_substrate | 0.1627 | 0.0466 | 37 |
| a2_adamts13_folded_complex | 0.1684 | 0.0347 | 81 |
| dprime_d3_fviii_binding | 0.1230 | 0.0076 | 100 |

**注意**：以下 assay_key 显示 iPTM=0 或 N/A，可能需要进一步排查解析逻辑：
- a1_aim_autoinhibition_context
- a2_folded_stability
- c_domain_assembly_context
- ck_dimerization_context
- d1d2_propeptide_context
- d4_assembly_context

---

## 已知问题

1. **部分 assay iPTM=0**：某些 confidence 文件的 iPTM 值为 0，可能是：
   - 这些 construct 的 interface 结合较弱
   - 或者 JSON 解析逻辑需要调整（某些字段可能为空）

2. **HuggingFace 上传**：受速率限制（128 commits/hour）影响，目前还在陆续上传中

---

## 数据文件

| 文件 | 说明 |
|------|------|
| `output/boltz2_vwd_functional_panel/boltz_results/` | 原始结构预测结果（7783 个 CIF + 7783 个 confidence JSON） |
| `output/boltz2_vwd_functional_panel/boltz_results_summary.csv` | 汇总解析结果 |
| `output/boltz2_vwd_functional_panel/job_manifest.csv` | 991 个 job 的元数据 |
| `output/boltz2_vwd_functional_panel/diagnostic_panel.csv` | 临床诊断面板 |

---

## 后续分析建议

1. **排查 iPTM=0 的 assay**：检查 confidence 文件中的 raw data
2. **阈值校准**：根据已知临床分类建立 GOF/LOF 阈值
3. **HuggingFace 增量上传**：等待 rate limit 恢复后继续

---

## 数据同步

- HuggingFace 仓库：`lucachangretta/VWF`
- 路径前缀：`vwd_functional_panel/`
- 拉取命令：
```python
from datasets import load_dataset
ds = load_dataset("lucachangretta/VWF", split="train", data_dir="vwd_functional_panel")
```
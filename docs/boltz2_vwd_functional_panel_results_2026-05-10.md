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

**注意**：以下 assay_key 是单链/monomer 构件，没有 inter-chain interface；iPTM=0 是预期现象，不应作为主指标：
- a1_aim_autoinhibition_context
- a2_folded_stability
- c_domain_assembly_context
- ck_dimerization_context
- d1d2_propeptide_context
- d4_assembly_context

这些构件应优先看 pTM、complex pLDDT、局部 pLDDT、PAE/RMSD 或下游提取的结构稳定性/暴露度特征。iPTM 只适合 D′D3-FVIII、A1-GPIbα、collagen、integrin、ADAMTS13 等复合物/界面任务。

---

## 已知问题

1. **单链 assay iPTM=0**：这是正常现象。iPTM 是 inter-chain 指标，单链稳定性/构象任务没有链间界面。

2. **HuggingFace 上传**：旧上传脚本按 job 目录上传，每个 job 产生 CIF/CONF 两次 commit，约 2000 commits，容易触发 128 commits/hour 限制。应使用修订后的 archive 或 single-folder 上传模式。

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

1. **建立 evidence matrix**：按 assay-specific WT baseline 计算 delta、z-score、rank，并为复合物/单链任务分别使用 iPTM 或 pTM/pLDDT。
2. **阈值校准**：根据已知临床分类建立 GOF/LOF 或 LOF/stability 阈值，特别区分 Type 2N、2M-A3、2A、2B/2M-A1。
3. **阴性对照校准**：GeneBe/ClinVar benign variants 应作为 negative/control variants，而不是只用 WT。
4. **Agent v3 输入改造**：将 Boltz/FoldX/AF3/AlphaGenome 统一为结构化 evidence，供诊断 Agent 模拟临床路径。

---

## 数据同步

- HuggingFace 仓库：`lucachangretta/VWF`
- 路径前缀：`vwd_functional_panel/`
- 拉取命令：
```python
from datasets import load_dataset
ds = load_dataset("lucachangretta/VWF", split="train", data_dir="vwd_functional_panel")
```

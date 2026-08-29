# Type-1 10例：服务器运行定稿

## 已定口径

- AlphaGenome官方是 **11类原始输出**，不是15类：ATAC、CAGE、DNASE、RNA_SEQ、CHIP_HISTONE、CHIP_TF、SPLICE_SITES、SPLICE_SITE_USAGE、SPLICE_JUNCTIONS、CONTACT_MAPS、PROCAP。
- 另跑官方 `RECOMMENDED_VARIANT_SCORERS`；当前SDK定义为 **19个评分视图**，其中包含active-track和polyadenylation派生评分，不再把它们叫“模态”。
- 请求前调用 `output_metadata()`；每个输出独立选择 `CL:0000115`、endothelial/vascular/HUVEC、megakaryocyte/platelet相关轨道。没有VWF相关轨道时写明缺失，不拿无关细胞补齐。
- 官方依据：[11类输出与ontology请求](https://www.alphagenomedocs.com/colabs/quick_start.html)、[输出元数据字段与轨道数](https://www.alphagenomedocs.com/exploring_model_metadata.html)、[官方推荐scorer定义](https://github.com/google-deepmind/alphagenome/blob/main/src/alphagenome/models/variant_scorers.py)。

## 精确任务量

| 模块 | 输入/任务 | 数量 |
|---|---:|---:|
| AlphaGenome | CASE_T1_001–003、005–010 | 9个变异 × 11原始输出，并运行官方推荐scorers |
| Boltz-2 | 5个错义变异 | 8个变异-assay任务 + 7个同assay WT = 15个任务 |
| Targeted MD | CASE_T1_008 P1413L | pipeline-matched WT + P1413L，各50 ns |
| 暂不直接建模 | CASE_T1_004大片段缺失/重复 | 1例 |

CASE_T1_007按 `NM_000552.5(VWF):c.3379+1G>A / GRCh38 chr12:6023630 C>T` 进入AlphaGenome。血友病A保留为共病；其FVIII:C不得用于支持或排除VWD 2N。

## 服务器顺序

```bash
cd /path/to/VWF-ETHos
export ALPHAGENOME_API_KEY='...'
bash outputs/type1_panel_agent_20260828/server_bundle/run_alphagenome.sh

GPUS=4 bash outputs/type1_panel_agent_20260828/server_bundle/run_boltz.sh

export FOLDX=/absolute/path/to/foldx
WT_GPU=0 MUT_GPU=1 bash outputs/type1_panel_agent_20260828/server_bundle/run_md_p1413l.sh
```

拿回 `outputs/type1_panel_agent_20260828/server_bundle/results/` 后：

```bash
bash outputs/type1_panel_agent_20260828/server_bundle/ingest_and_test_agent.sh
```

该步骤会校验AlphaGenome 9例、Boltz 15个完整任务，生成 `server_results/type1_computational_evidence.csv`，再通过 `LocalComputationalPanelProvider` 嵌入LangGraph Agent重跑10例。

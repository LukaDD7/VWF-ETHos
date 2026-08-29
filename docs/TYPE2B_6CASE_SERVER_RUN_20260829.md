# Type-2B 6人服务器运行清单（2026-08-29）

输入已去除姓名、病案号、出生日期和就诊日期。共6名患者、7条变异记录；患者5有A1461D与良性对照D2449N，V1316M见于两名患者。

- AlphaGenome：7个GRCh38请求；每个请求覆盖11类官方原始输出，并调用完整推荐 scorer。运行时逐输出读取 `output_metadata()`，只选择VWF相关内皮/血管/巨核细胞轨道；没有适用轨道时显式记缺失。
- Boltz-2：6个唯一错义变异，16个变异构建加4个同构建WT基线，共20个任务。V1316M只计算一次，回灌时映射给两名患者。
- MD：R1308C、S1310F、V1316M、R1341W、A1461D与一个管线匹配WT，共6条50 ns轨迹。D2449N为C3良性对照，不进入A1 MD。
- Agent：结果回传后先检查7/7 AlphaGenome、20/20 Boltz和6/6 MD，再生成患者级证据CSV并运行现有LangGraph流程；缺失结果保持 `pending`。

服务器入口位于 `outputs/type2b_panel_agent_20260829/server_bundle/`，按 `run_alphagenome.sh`、`run_boltz.sh`、`run_md_a1_panel.sh`、`ingest_and_test_agent.sh` 执行。

# 服务器任务需求：补齐 16 例计算证据并保留完整模型响应

日期：2026-09-05  
优先级：P0（下一轮 16 例 Agent 评估的前置条件）  
任务性质：数据补齐、无损归档、FHIR 可检索投影；**本任务不运行最终临床 Agent，不输出或优化分型答案。**

## 1. 背景与问题

当前 16 例临床 Agent 冻结输入并不具有一致的数据完整性：

- 旧 AlphaGenome 产物对大多数变异只保留 `ALT-REF` 或 scorer 摘要，没有分别持久化完整 REF、ALT 数组和轨道元数据。
- 2026-09-04 已对 P1413L、V1316M 验证真实 REF/ALT 保存与 FHIR 传递，但不能据此声称 16 例都已补齐。
- Boltz、MD、SMD 及其他结构证据散落在 Git、本地/server 输出和可能的 Hugging Face 资产中；“已有输出”“只有摘要”“只有运行计划”“原始轨迹齐全”尚未形成逐变异、逐协议的统一清单。
- FHIR 应向 Agent 提供完整、可追溯的测量语义和按需读取入口，而不是将大型数组直接塞进 prompt，也不能用同批病例排名、既有分类或人工分型提示污染输入。

服务器侧先审计已有资产，能复用的必须校验后复用；确实缺失且具备输入条件的再补算。所有缺失、失败、不适用都必须显式记录，禁止用零、另一病例、同位点不同替换或数学 comparator 冒充真实结果。

## 2. 病例和变异范围

病例全集来自下列两个既有输入包，禁止从带有既有分类的列构造模型提示：

- Type 1：`outputs/type1_panel_agent_20260828/server_bundle/input/`
- Type 2B panel：`outputs/type2b_panel_agent_20260829/server_bundle/input/`

范围为 16 位患者：

```text
CASE_T1_001 ... CASE_T1_010        10 patients
CASE_T2B_001 ... CASE_T2B_006       6 patients
```

特殊情况：

- `CASE_T1_004` 是大片段缺失/重复，当前没有可用精确 HGVS/GRCh38 区间。不得伪装成 SNV/indel，也不得运行不匹配的 AlphaGenome/Boltz。先尝试从**不含分类答案的原始结构变异字段**规范化；若原始资料不足，交付状态必须为 `blocked_missing_variant_definition`，同时保存需要补充的最小字段（变异类型、GRCh38 breakpoints/CNV interval、copy number）。该病例仍计入 16 例完整性清单。
- `CASE_T2B_005` 有两个变异行，必须分别保存和查询，患者层只记录 phase 未知，不得把两个变异的计算值合并。
- V1316M 出现在两位患者中。相同规范化变异和完全相同协议允许内容寻址去重，但每位患者必须有独立的 case-to-artifact 映射；不得混用不同协议/WT。
- 当前可直接请求 AlphaGenome 的是 16 个规范化变异：Type 1 九个、Type 2B 七个；患者总数仍为 16。

## 3. 总体完成定义

对每个 patient × variant × model × construct/protocol，产出一个不可歧义的状态：

```text
available_complete
available_partial
not_applicable
blocked_missing_input
not_supported_by_model
failed_retry_exhausted
not_run_out_of_scope
```

`complete` 必须同时满足：原始数值/结构文件、匹配参考、元数据、运行配置、软件版本、校验和与 FHIR 投影均存在且互相能通过 ID/hash 追溯。只有摘要 CSV 不算完整；只有原始大文件、但没有可查询的测量投影，也不算完整。

先产出 inventory，再执行补算计划。不要因为文件名存在就认为成功，必须验证格式、大小、hash、预期对象/轨迹帧和病例映射。

## 4. AlphaGenome（所有 16 个可规范化变异均为 P0）

### 4.1 请求与版本冻结

以 `server_bundle/input/alphagenome_requests.csv` 为规范化请求来源，逐行复核：assembly、chromosome、position、REF、ALT、1 Mb interval、ontology terms、requested outputs。保存：

- AlphaGenome client/package 版本、服务端可获得的模型/version identifier、请求时间、request ID（若 API 返回）。
- 完整 request JSON/等价规范化表示；不得保存 API key。
- 输出类型及官方元数据 inventory；ontology/cell selection 和 unavailable 原因。
- 对请求失败最多 3 次传输重试，保存每次状态、错误类型、HTTP/API 可公开元数据和时间；不要保存密钥或 Authorization header。

### 4.2 必须保留的响应

对 API 实际返回的每个 output type，分别保存：

- `result.reference.<modality>` 的完整数值数组；
- `result.alternate.<modality>` 的完整数值数组；
- interval、resolution、strand/junction coordinates、维度和 dtype；
- 完整 track metadata，保持原始行顺序；
- REF/ALT 对齐检查结果；
- 同一 track、同一位置的 `ALT-REF`；
- 官方 variant scorer 的原始 tidy/AnnData 可复现产物、scorer 名称、mask、方向语义和单位；
- API 返回的其他非敏感结构化 metadata。若 SDK response 不能直接安全序列化，以 NPZ/Zarr/Parquet + JSON sidecar 无损保存全部可访问字段，不依赖不安全 pickle 作为唯一副本。

原始数组与 Agent 投影必须分开：

```text
raw/        完整 REF/ALT 数组、track metadata、scorer 原始结果
derived/    同轨配对摘要、明确命名的 gene-mask scorer、QC
fhir/       Observation/DocumentReference 投影及 artifact URI/hash
```

应复用已实现的 `scripts/pipeline/run_type1_10case_alphagenome.py::persist_paired_tracks`，但服务器需验证它对所有返回对象没有丢字段。禁止从旧 delta 反推 REF 或 ALT。ATAC/CONTACT_MAPS 等在所选细胞条件下若 API 不返回，记录 `not_supported_by_model`/明确 API 原因，不得填 0。

### 4.3 AlphaGenome 验收

- 16/16 可规范化变异有请求记录和终态；P1413L、V1316M 可复用但必须重新核对 hash/version。
- 每个 available modality 都有成对 REF/ALT、相同 track identity、shape/interval/resolution 对齐证明。
- 全量数组可重新计算已交付的每个 paired summary，数值一致。
- FHIR/Agent 查询只返回当前病例变异，不泄露同批其他病例数值、rank、percentile 或标签。

## 5. Boltz / 静态结构证据（逐变异、逐生物轴）

先读取两个 `server_bundle/boltz/job_manifest.csv`、`diagnostic_panel.csv` 和 YAML，建立 intended-job inventory。对每个计划任务及匹配 WT construct：

- 保存输入 YAML/FASTA/MSA/template/constraint（实际使用者）、随机种子、模型和软件版本、命令、GPU环境、开始/结束时间。
- 保存 API/runner 的完整非敏感 response、所有预测结构、confidence JSON、PAE/PDE/distogram 等实际返回文件、日志和退出状态。
- 结构/指标必须标明链、残基编号映射、construct、WT 关系和 metric definition。
- matched WT 必须来自相同 construct、序列边界、配体/伙伴、约束和版本。不同 construct 的 WT 不可互换。
- `primary_value`、pTM、pLDDT、ipTM、PAE 等是模型置信度/几何指标，不得改名为结合自由能或功能效应。

不要求对不适用的非错义变异强行跑单体结构。应记录 `not_applicable` 及原因。对已有 job 先校验原始输出是否齐全；只有 summary CSV 而缺 response/结构文件时标记 partial，并尝试从服务器任务目录或远端资产恢复。

## 6. MD / SMD / 其他物理模拟证据

### 6.1 先做资产审计，不盲目补算

搜索服务器工作目录、Git 跟踪资产和项目既有 Hugging Face dataset/repository。逐 variant × construct × protocol 区分：

- 只有计划或 starting structure；
- 运行中/失败；
- 只有提取特征；
- 完整原始轨迹；
- 原始轨迹 + matched WT + 特征 + QC。

同位点不同替换不可替代；同一变异不同批次不可混用 WT。旧 label、known subtype、cohort z-score 不进入病例计算证据。

### 6.2 完整 MD 包

对确实已经运行或服务器能按既有正式协议完成的项目，保存：

- 初始/最终结构、topology、parameters、index、MDP、restraint、映射文件；
- TPR、XTC/TRR、EDR、log、checkpoint、stdout/stderr；
- 周期边界处理/fit 后派生轨迹与生成命令（不能替代原始轨迹）；
- protocol ID、force field、水/离子、温度、时间步长、总模拟时长、seed、GROMACS/依赖版本和硬件；
- 同协议 matched WT；
- 逐时间窗/逐帧基础测量、汇总、QC 和特征提取脚本版本；
- SHA256、文件大小、帧数、起止时间、轨迹连续性校验。

机制轴按真实 construct 决定，例如 AIM–A1 equilibrium MD 可提供自抑制接触/遮蔽/盐桥和部分局部动力学代理；它不能自动声称测得 ΔTm、ΔΔGbind、高亲和态比例或临床分型。没有适配构建体/正式协议的变异记录 `not_run_out_of_scope`，并提出建议而不是临时发明一种模拟。

既有 `slow025` SMD 已被项目标记为 `complete_no_go`：原始文件可归档，但必须在 Agent-serving manifest 中 `eligible_for_agent=false`，不得作为方向性证据。其他 SMD 也只有在协议已经通过内部验收时才能进入 FHIR。

### 6.3 大文件存储

不要把轨迹或数百 MiB AlphaGenome 数组直接 commit 到 Git。使用项目既有的服务器对象存储或 Hugging Face dataset（遵循现有权限设置）保存大文件；Git 只保存需求、manifest、schema、小型索引和可复现脚本。每个 artifact URI 必须可由服务器 Agent 使用配置化凭证读取，并带 immutable revision 与 SHA256。不得把访问 token 写入仓库或 FHIR。

## 7. 其他计算结果

盘点并纳入确有病例级原始结果的 FoldX、结构界面/PAE 分析和已验证的派生模型，但必须满足：

- 有原始输入、完整输出、软件/参数、匹配 reference 和可复现脚本；
- 清楚标识它测量/预测的物理量，不能重复包装同一底层结果后被 Agent 多次计权；
- 没有 matched reference 时只能说明病例内行为；
- 失败、未校准或 no-go 模型仍保留 provenance，但设置 `eligible_for_agent=false`；
- 不把 ClinVar/HGMD/工作表分类或医生结论作为计算模型输入/校准标签加入本批病例证据。

## 8. FHIR 和 Agent 按需读取接口

目标不是把全部大数组塞进一个 Observation，而是提供两层接口：

1. **FHIR evidence index**：每一病例/变异/model/construct 对应 Observation 或 DiagnosticReport；含测量定义、病例值、匹配 WT/REF、单位、方向可用性、协议 ID、QC、limitations、artifact hash/URI。
2. **artifact query tool**：Agent 可按 `case_id + variant_id + model + modality/measurement_code + track/region` 请求精确 REF/ALT/WT、时间窗或全文 metadata；返回的数据仍带原始 artifact hash 和坐标范围。

大型原始文件可用 `DocumentReference` 指向 immutable artifact；若使用 FHIR `Binary`，也只用于合理大小的小对象。FHIR ID 不能代替人类可读引用信息。失败/OperationOutcome 若对 Agent 可见，必须进入允许引用集合并带状态语义，避免“模型看得到但输出校验禁止引用”。

所有病例级 FHIR Bundle 应通过项目 schema 检查，并尽可能运行官方 HL7 FHIR validator；若只通过内部 schema，报告必须明确写 `FHIR-shaped/internal validated`，不得声称已生产级互操作合规。

## 9. 目录和交付物

建议服务器在一个新的、不可覆盖的 run 目录生成：

```text
computational_evidence_16case_<UTC timestamp>/
  run_manifest.json
  case_inventory.csv
  variant_inventory.csv
  model_protocol_inventory.csv
  artifact_manifest.parquet
  missingness_report.csv
  retrieval_test_report.json
  alphagenome/<variant_id>/{request,raw,derived,fhir}/...
  boltz/<variant_id>/<construct_id>/{request,raw,derived,fhir}/...
  md/<variant_id>/<protocol_id>/{input,raw,derived,qc,fhir}/...
  fhir/<case_id>/bundle.json
  README.md
  SHA256SUMS
```

必须交付：

- 完整性矩阵：16 patients、17 variant rows（含不可规范化的 T1_004 描述和 T2B_005 第二变异）、所有模型层状态及缺失原因。
- `case_id → variant_id → artifact_id → FHIR resource_id` 映射。
- 每个大文件的 immutable URI、revision、SHA256、bytes；本地仅有路径不算服务器 Agent 可用。
- 从 raw 重新生成 derived/FHIR 的命令与环境 lock/version。
- 一份执行报告：复用了什么、补算了什么、仍缺什么、失败重试详情、费用/耗时（若可得）。
- 一份安全报告：确认无 API key/token、无工作表既有分类/HGMD公司结论、无跨病例统计泄漏。

## 10. 自动验收测试

服务器提交前至少实现并执行以下测试：

1. **Coverage**：16 个患者均在 matrix 中；多变异病例不丢行；T1_004 明确终态。
2. **Identity**：assembly/REF/ALT/HGVS/variant_id 一致，REF 与参考基因组核验；坐标冲突不得静默覆盖。
3. **Alpha pairing**：所有 available modality 的 REF/ALT shape、metadata、interval、resolution 对齐；随机抽取和指定 C1950Y 均能从 raw 重算 FHIR 中的数值。
4. **WT compatibility**：Boltz/MD 的 WT construct/protocol compatibility ID 完全匹配。
5. **Artifact integrity**：全量 hash 校验；轨迹可读取且帧数/时间连续；LFS/pointer 文件不能误当原始数据。
6. **FHIR round trip**：写入后按 case/variant/modality 查询，返回数值、reference、单位、协议和 provenance 不丢失。
7. **Isolation**：查询 C1950Y 不返回其他病例值；同批 rank/label/既有分类字段为 0。
8. **No fabricated missingness**：不可用/失败项没有 `0` placeholder，另一患者/替换/WT 未被冒充。
9. **Response preservation**：API 原始返回可访问字段与归档 field inventory 对照，无未解释丢失；derived 可完全追溯至 raw。
10. **Agent smoke query**：不运行分型报告，只测试 Agent 工具可自主列出该病例可用模型、读取一个 AlphaGenome REF/ALT track、一个 matched-WT 结构/MD测量以及明确缺失项。

## 11. 服务器 Agent 回报格式

服务器 Agent 完成或遇到阻碍后，请提交并推送：

- 数据/脚本改动和小型 manifest；大文件只提交 URI/hash，不进入普通 Git object。
- `docs/SERVER_16CASE_COMPUTATIONAL_DATA_COMPLETION_REPORT_<date>.md`。
- 报告开头给出以下汇总：

```json
{
  "patients_total": 16,
  "variant_rows_total": 17,
  "alphagenome_requestable_variants": 16,
  "alphagenome_complete": 0,
  "boltz_jobs_intended": 0,
  "boltz_jobs_complete": 0,
  "md_protocols_complete_with_matched_wt": 0,
  "fhir_case_bundles_complete": 0,
  "artifacts_hash_verified": 0,
  "remaining_blockers": []
}
```

数值必须由 inventory 自动生成，不能手填估计。若某类模型不适用于一个病例，“完整”是正确记录 `not_applicable`，不是强行产生输出。

## 12. 停止条件

满足以下条件前，不要触发新的 16 例临床 Agent 分型运行：

- 完整性矩阵和 raw/derived/FHIR 三层通过上述测试；
- C1950Y 等旧产物病例已具有真实配对 REF/ALT，或明确记录不可补齐的外部阻碍；
- 可用 Boltz/MD 均有兼容 reference，no-go/未校准输出不会进入 Agent；
- 服务器报告已推送，且本地 review 能通过 URI/hash 复核。

本需求不要求通过提示词让 Agent 得到既定 Type 1/2B 结论。数据与语义完整性优先，后续独立评估提示词和 Agent 行为。

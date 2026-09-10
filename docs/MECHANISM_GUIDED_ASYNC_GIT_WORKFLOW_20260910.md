# VWD 机制驱动的异步 Git 计算工作流

日期：2026-09-10

状态：v1.1 控制面已实现并通过本地 bare-remote 闭环；具体 GPU protocol runner 与 calibration 数据仍需在服务器配置

用途：研究型决策支持，不能替代 VWD 专科诊断、功能实验或 ACMG/AMP PS3

## 1. 结论

能实施，而且这是当前收益最高的流程重构之一。建议保留现有三阶段外壳，把 Stage 2 改为：

```text
Stage 1  patient phenotype + genetics
   ↓
Stage 2A  low-cost retrieval: ClinGen/ClinVar/PubMed/local FHIR
   ↓
Stage 2B  competing mechanism hypotheses + M_UNKNOWN
   ↓
Stage 2C  evidence gap + value-of-information gate
   ↓
Stage 2D  select a locked protocol (Agent decides WHAT; protocol fixes HOW)
   ↓
Stage 2E  immutable FHIR Task → one Git branch per task → offline GPU agent
   ↓
Stage 2F  protocol-conformant Result + QC + controls + artifact hashes
   ↓
Stage 2G  local validation → FHIR Observation → update hypotheses / optional re-plan
   ↓
Stage 3  subtype posterior + mechanism explanation + uncertainty + next test
```

异步执行不降低 Agent 自主性。自主性体现在 Agent 决定“缺什么证据、要证伪哪个机制、哪项计算最有区分力”，而不是临场决定 PDB、力场或模拟时长。与预先把所有病例全部跑完再让 Agent 解读相比，这个设计能清楚区分：

1. acquisition policy 是否选择了正确问题；
2. MD/AI protocol 是否真的测到了该机制；
3. 合并后是否提高了最终分型。

## 2. 文献核对后采用的科学边界

### Type 2B：AIM-A1 自抑制丢失不能压缩成 RMSD

Legan 等在 glycosylated AIM-A1 上比较了 7 个 Type 2B 变异，联合热稳定性、GPIbα binding kinetics、HDX-MS 与单分子力学。共同方向是 AIM-A1 稳定性/机械自抑制降低以及 secondary GPIbα-binding site 暴露增加，但不同变异程度不同。因此 v1 锁定：

- `α3β4 loop`（VWF 1366–1384）的 SASA + RMSF；
- `β3α2 loop`（1320–1340）的 SASA + RMSF；
- AIM-A1 contact、H-bond、salt-bridge occupancy；
- DSSP secondary-structure occupancy；
- matched WT 与已知 2B、2M、边界对照的 joint PCA/state population。

SASA 是几何暴露，HDX 还包含氢键保护与局部 breathing；所以只能把 `SASA + RMSF + contacts/secondary structure` 作为 HDX-derived computational surrogate，不能写成 `HDX = SASA`。平衡 MD 也不能输出 AIM release force，v1 明确禁用这一结论。

### Type 2M：必须允许方向相反的两种 A1 机制

Tischer 等对 15 个 A1 Type 2M 变异的实验表征显示，既有局部 disorder/α-helical loss，也有 native fold hyperstabilization；两条相反的结构路径都可能导致 GPIbα binding/flow adhesion 下降。因此 registry 分开保留：

- `M7_2M_A1_DISORDER`；
- `M8_2M_A1_HYPERSTABLE`。

不能使用“越不稳定越像 2M”的单轴阈值。

### Type 2A：区分可计算的 proteolysis 与弱可观测的 assembly

Seidizadeh 等 2025 年对 65 位患者、15 个 A2 变异进行了完整表型和 plain MD、scaled MD、docking/MM-PBSA 分析。v1 复用其 A2 路线：3GXB、FoldX、AMBER99SB-ILDN、TIP3P、plain MD、scaled MD `lambda=0.8`、A2–ADAMTS13 docking 与 MM-PBSA，并把 1603–1606 cleavage-site SASA、α6 contacts、RMSD–SASA joint state 作为机制 readout。

必须写清楚：scaled MD 是降低势垒的 enhanced-sampling perturbation，不是血流 shear；MM-PBSA 是相对 interaction energetic proxy，不是实验 `Kd`。Type 2A assembly/secretion 是真实机制，但单体短程 MD 对它的可观测性较低，v1 不允许用单体 RMSF 直接确认 multimerization defect。

### 2M/2A 边界：同一结构域不能直接编码答案

Seidizadeh 等 2024 年把争议性 R1315/R1374 变异放入完整表型、A1-A2 建模和 MD 中，显示同一区域可以形成 2M/2A mixed phenotype。v1 因而保留 `M10_2M2A_INTERFACE`，并强制任何 A1 任务同时保留多条竞争机制与 `M_UNKNOWN`。

## 3. 已实现内容

### 3.1 唯一事实源

[`protocols/vwd_mechanistic_v1/registry.json`](../protocols/vwd_mechanistic_v1/registry.json) 同时定义：

- 11 条已知机制 + `M_UNKNOWN` 开放集；
- 12 类通用 measurement semantics；
- 3 条可提交的、带版本和 SHA-256 digest 的 protocol；
- PDB/construct、FoldX、execution tier、固定参数、controls、ROI、required measurements、QC gates、文献 provenance 和禁止外推的边界。

第一版三个 protocol：

| Protocol | 主要问题 | 起始结构 | 必需的机制 readout |
|---|---|---|---|
| `VWF_A2_PROTEOLYSIS_SEIDIZADEH2025_V1` | A2/ADAMTS13、2A proteolysis | 3GXB | cleavage-site SASA、α6 contacts、RMSD–SASA joint state |
| `VWF_A1A2_INTERFACE_SEIDIZADEH2024_V1` | 2M/2A interface | 7EOW A1 + 7GBX A2 + linker model | Rg、interface RMSF/contacts、joint PCA |
| `VWF_A1_2B_GOF_LEGAN2023_V1` | 2B GOF vs 2M disorder/hyperstable | 7A6O；1SQ0 作 GPIb interface reference | two-loop SASA/RMSF、AIM contacts/H-bond/salt bridge、DSSP、PCA |

生成的 FHIR R5 [`PlanDefinition`](../protocols/vwd_mechanistic_v1/fhir/PlanDefinition-vwd-mechanism-guided-acquisition-v1.json) 表示动作空间，3 个 [`ActivityDefinition`](../protocols/vwd_mechanistic_v1/fhir/) 表示可复用 protocol template。

### 3.2 Planner/Task/Result guardrails

[`src/vwd_clinical_agent/mechanistic_tasks.py`](../src/vwd_clinical_agent/mechanistic_tasks.py) 已实现：

- domain/variant class 只生成候选机制，不直接生成 subtype；
- 低成本 retrieval 必须完成；已有证据为 `adequate` 时禁止占用 GPU；
- `expected_information_gain >= 0.10`；
- 至少一个竞争机制，且候选集中必须保留 `M_UNKNOWN`；
- 当 `P(2M)=0.90`、`P(2B)=0.05` 时拦截 2B task；当 2B/2M 接近时允许提交；
- protocol version、registry digest、protocol digest、request digest、controls 与 required measurements 全部锁定；
- `source_commit` 必须是仓库中真实存在的 commit，且该 commit 内 registry 的实际字节必须产生 Task 声明的 digest；
- server result 必须回显同一 task/request/protocol identity；
- QC 通过的 completed result 必须覆盖全部 required measurements；
- scalar `case - matched WT` delta 会被重新计算；
- `supports/contradicts/indeterminate` 必须互斥，目标机制必须恰好出现一次；
- artifact 只能使用仓库相对 POSIX 路径；接收端拒绝绝对路径、`..`、symlink、目录和越界路径，并重新计算实际文件 SHA-256；
- 验证后的结果转换为 `Task + Observation + DocumentReference` FHIR bundle。

### 3.3 命令行与服务器 Agent 约束

[`scripts/mechanism_task.py`](../scripts/mechanism_task.py) 提供：

```bash
python scripts/mechanism_task.py list --domain A1
python scripts/mechanism_task.py submit proposal.json
python scripts/mechanism_task.py validate-request <task_id>
python scripts/mechanism_task.py result-template <task_id>
python scripts/mechanism_task.py validate-result <task_id> result.json
python scripts/mechanism_task.py ingest-result <task_id> result.json
python scripts/mechanism_task.py status <task_id>
python scripts/mechanism_task.py export-fhir
python scripts/validate_mechanism_fhir.py
```

所有运行时 contract 错误返回结构化 JSON 和非零 exit code，不再输出 Python traceback。`validate-request`、`result-template`、`status` 等既接受 task ID，也接受 request 文件/目录。

FHIR 检查默认强制 registry/export parity 和必要字段；CI/服务器安装官方 `validator_cli.jar` 后使用 `python scripts/validate_mechanism_fhir.py --require-official --validator-jar /opt/hl7/validator_cli.jar`，未提供 jar 时该模式会明确失败，避免把内部检查标成官方验证。

[`scripts/mechanism_git.py`](../scripts/mechanism_git.py) 实现 request/result 双分支传输和原子 claim；[`src/vwd_clinical_agent/mechanistic_workflow.py`](../src/vwd_clinical_agent/mechanistic_workflow.py) 与 [`scripts/run_mechanism_workflow.py`](../scripts/run_mechanism_workflow.py) 实现 checkpointed `submit → interrupt → resume → ingest` 子流程。

[`mechanism_tasks/AGENTS.md`](../mechanism_tasks/AGENTS.md) 是 offline GPU Agent 的强制 protocol：request/registry 不得修改；不能临场换结构/参数；缺输入或 QC 失败必须返回 failed/inconclusive；原始 trajectory 不进 Git；只允许输出机制支持/反驳/不确定。

示例 proposal：[`examples/mechanism_tasks/a1_v1316m_proposal.json`](../examples/mechanism_tasks/a1_v1316m_proposal.json)。

## 4. 已实现的 Git 任务生命周期

Git 适合做低并发、强审计的 handoff，但不应伪装成作业调度器。采用“一任务一分支”，不要让多个 GPU worker 同时写一个 queue branch。

### 4.1 本地 Agent：提交任务

在包含 protocol registry 的稳定提交上运行：

```bash
python scripts/mechanism_git.py publish proposal.json
```

该命令从 `source_commit` 创建临时 worktree，把 request/FHIR Task 作为唯一新增内容提交到 `mechanism-task/<task_id>` 并推送。当前工作区可以有无关改动，它们不会进入任务分支。若 commit 不存在、registry 有未提交改动、registry digest 不匹配或同名任务分支已经存在，发布会失败。proposal 和 Task 中只能使用不可逆的研究 ID；不得放姓名、病历号等直接标识符。

### 4.2 Offline GPU Agent：执行并返回

服务器只需联网访问 Git；计算节点本身可完全离线。服务器先 claim：

```bash
python scripts/mechanism_git.py claim <task_id> \
  --worktree-root /srv/vwf-mechanism-worktrees \
  --worker-id gpu-01

cd /srv/vwf-mechanism-worktrees/<task_id>
python scripts/mechanism_task.py validate-request <task_id>
python scripts/mechanism_task.py result-template <task_id>

# 按 request 指定的 protocol/tier 执行现有 FoldX/GROMACS/analysis runner
# 写 mechanism_tasks/results/<task_id>/result.json

python scripts/mechanism_git.py return <task_id> --worktree "$PWD"
```

claim 从 task branch 建立 `mechanism-result/<task_id>` 独立 worktree 并立即推送 claim commit；并发 worker 的第二次 claim 会因 result branch 已存在或 non-fast-forward push 而失败。`return` 会再次校验 request、result、measurement、QC 和每个 artifact 的实际字节，只提交 result、FHIR bundle 和明确声明的 compact artifacts。服务器不应提交 `.xtc/.trr` 等大文件。

若已有固定 wrapper，可以把命令保存在服务器本地 allowlist 中并运行：

```bash
python scripts/mechanism_git.py run <task_id> \
  --worktree /srv/vwf-mechanism-worktrees/<task_id> \
  --runner-config /etc/vwf/worker_config.json
```

allowlist 以 protocol ID 为键，支持 `{request}`、`{result}`、`{template}`、`{task_id}` 和 `{repo_root}` 占位符。配置必须来自服务器，Task/registry 不能提供任何可执行命令。

### 4.3 本地 Agent：验收并恢复推理

```bash
python scripts/mechanism_git.py collect <task_id>
```

`collect` 在临时 detached worktree 校验 result branch，不需要把服务器分支 merge 到当前分支；通过后才按 append-only 规则复制 request/result/compact artifacts 并生成 FHIR bundle。只有 `mechanism_tasks/ingested/<task_id>/bundle.fhir.json` 能进入 Stage 2G/Stage 3。

LangGraph 提交与恢复：

```bash
python scripts/run_mechanism_workflow.py \
  --thread-id <case-and-gap-id> --publish-to-git submit proposal.json

python scripts/mechanism_git.py collect <task_id>

python scripts/run_mechanism_workflow.py \
  --thread-id <case-and-gap-id> resume
```

第一次命令在 `awaiting_mechanism_compute` interrupt 保存 SQLite checkpoint；第三次从同一 checkpoint 恢复并返回经过验证的 FHIR bundle，而不是整例重跑。

## 5. 研究设计收益与风险

### 主要收益

1. **真正的 bounded autonomy**：Agent 选择科学问题，protocol 保证复现性。
2. **可证伪机制链**：`variant → molecular state → functional consequence → subtype tendency`，不再是 `RMSD → subtype`。
3. **错误归因清楚**：能分别测 tool validity、acquisition policy 和 end-to-end diagnosis。
4. **计算成本可解释**：只有可能改变 posterior 的任务才进入 GPU。
5. **允许负结果和 abstention**：simulation 可反驳假设，也可返回 out-of-distribution/inconclusive。

### 主要风险及缓解

| 风险 | 缓解 |
|---|---|
| mechanism-space leakage：domain mapping 暗示 subtype | 多竞争机制；`M_UNKNOWN` 强制保留；按 variant/domain holdout 评估 |
| simulation surrogate 未经验证 | known-variant calibration set；不能恢复实验 phenotype 的 protocol 不进入临床 reasoning |
| 短 MD、单 replica、模型结构噪声 | pilot/production 分层；matched WT；production 多 replica；结构/QC/replicate consistency 明示 |
| 将 AI/MD 当 PS3 或诊断确认 | contract 和 server Agent 同时禁止；报告只写 mechanism consistency |
| Git 并发/大文件/PHI 风险 | 一任务一分支；raw trajectory 外置；opaque ID；大规模后换真正队列/对象存储 |

## 6. 剩余实施顺序

控制面、Git transport、通用 server-runner 接口及 LangGraph interrupt/resume 已完成。剩余工作是科学执行面和主流程决策接线：

1. **服务器 protocol wrapper**：把现有 7A6O/3GXB/GROMACS/FoldX 输出严格映射到 registry measurement IDs，先完成 2B，再完成 2A、A1-A2。通用 runner 已就绪，但示例 `/opt/vwf/bin/...` 必须替换为服务器真实 wrapper。
2. **Calibration**：已知 2B vs 2M、2A experimental positives/negatives 做 blind calibration；产出固定 calibration-set ID 和阈值/相似度模型。
3. **主 Clinical Agent policy 接线**：当前 checkpointed mechanism 子图接受已验证 `TaskProposal`。待主三阶段 Agent 明确输出 evidence gap/VOI 后，将该子图组合到 Stage 2；不能用 domain 规则偷偷替 Agent 生成 subtype。
4. **Result bank emulator**：在 GPU 结果积累前，用已有 precomputed results 回放相同 Task/Result contract，先测 acquisition policy。
5. **三层评估**：
   - Level I tool validity：simulation vs experimental molecular phenotype；
   - Level II acquisition policy：缺证据病例是否选对 mechanism/protocol；
   - Level III end-to-end：accuracy、macro-F1、Brier/calibration、abstention、cost 与增量价值。

## 7. 验收标准

v1 不能以“能生成 Task”为验收结束，至少需要：

- request 在任意机器上产生相同 registry/protocol digest；
- 不存在的 commit、commit 内 registry 不匹配、source commit 非任务分支祖先时必须拒绝；
- server 无法提交缺失 required measurement 或错误 comparator 的 QC-passing result；
- 改 request/protocol 任一字节会导致 digest 校验失败；
- artifact 缺失、SHA-256 不符、绝对/越界/symlink 路径必须拒绝；
- 并发 worker 只能有一个成功创建 result/claim branch；
- A1 candidate space 始终包含 2B、两类 2M、边界/分泌和 unknown，而非 `A1 = 2B`；
- Result 转换后的 FHIR Observation 保留 case、matched WT、signed delta、units、QC、benchmark、limitations 与 artifact hash；
- Stage 3 只消费已验证 FHIR bundle，并能输出 supported/contradicted/inconclusive，而不是 `MD confirms subtype`。

## 8. 主要依据

1. Legan ER, et al. *Blood*. 2023. Type 2B VWD mutations differentially perturb autoinhibition of the A1 domain. <https://doi.org/10.1182/blood.2022017239>
2. Seidizadeh O, et al. *Blood Advances*. 2024. Type 2M/2A VWD: a shared phenotype between type 2M and 2A. <https://doi.org/10.1182/bloodadvances.2024012626>
3. Seidizadeh O, et al. *RPTH*. 2025. Deep molecular modeling and mechanistic insights into type 2A VWD. <https://doi.org/10.1016/j.rpth.2025.103233>
4. Tischer A, et al. *JTH*. 2025. Structure resolved dynamics of type 2M VWD. <https://doi.org/10.1016/j.jtha.2024.12.026>
5. HL7 FHIR R5 Workflow, PlanDefinition, ActivityDefinition, Task, and Observation. <https://hl7.org/fhir/R5/workflow.html>

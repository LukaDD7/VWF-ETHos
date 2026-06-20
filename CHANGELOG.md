# Changelog

All notable changes to the VWF-ETHos project are documented here.

## [Unreleased]

### Added (2026-06-20) — 7A6O AIM-A1 MD 特征 → 分类器 (轴B' 结合面完整性)

- **A40 push 回 50 ns MD 结果** (`md_data/7a6o_reference_md/`, Git LFS, WT + 14 参考变体) pull 落地 (需 `git lfs`, 校验和与 `manifest.csv` 一致)。构件原生编号 1262–1466 经 12 个突变位点比对确证。
- **机制调研 + 符号纠偏**：平衡(无剪切) MD 看不到力依赖的 2B 松开；医学上稳健的可观测量是 **AIM→GPIbα 结合面屏蔽接触的保留率 (= 结合面完整性)**，而非 AIM↔A1 总接触数(后者方向相反、误导)。数据驱动的 WT 屏蔽界面 = A1 [1305-1308,1313,1376,1378,1380,1410,1434]，恰盖 2B 突变位点 → 验证 Deng 2017 模型。保留率 **WT 0.81 > 2B 0.75 > 2M 0.63**(2M 结合面坏 → 屏蔽不住)。
- **`scripts/pipeline/extract_7a6o_md_features.py`**(新)→ `output/md_7a6o_features.csv` + `md_7a6o_masking_interface.json`。`md_face_destab_score` 高 = 结合面破坏 = 2M/LOF；AUC(2M>2B)=0.83(单副本, n(2M)=3, 有重叠 → 软证据)。
- **`agentic_vwf_classifier.py`**：新增**轴B' MD 结合面完整性**(`MD_FACE_DESTAB_2M_Z=1.0`)。RULE6 中：轴B 已判 2M 时 MD 确证→升置信(0.72→0.8)；`uncertain` 兜底前加 MD tie-breaker(强结合面破坏→2M 软证据 0.55)，排在 2B 热点/AIM 位之后→**绝不硬翻 2B**。`AIM_RELEASE_2B_Z` 处加警示：勿用平衡 MD 接触下降反推(符号相反)。
- **参考集回归**：15 变体过分类器，**2B→2M 误翻 = 0**(V1316M 越阈但 1316∈热点先判 2B)，**R1374H(2M) 由 uncertain 救回 2M**，其余与 NaN 完全一致(向后兼容)。详见 `docs/7A6O_MD_FEATURE_ANALYSIS_2026-06-20.md`。

### Status / 两条线大局 (更新 2026-06-11)

- ✅ **静态特征分类器 = 已优化(可用,待标签校准）**：`agentic_vwf_classifier.py` RULE6 多轴方向判别（2B=GOF / 2M=LOF），轴B forced_binding+heparan 联合 LOF（校准 30% 2M 召回 / 5% 2B 误伤），无方向信号→uncertain。
- 🟢 **动态 MD = 突破**：弃 Boltz D'D3-A1（EM 跑飞 5.9e9），改用**实验结构 7A6O（AIM-A1, 2.12 Å）**做 WT 骨架 + FoldX 在骨架上造突变体。WT 弛豫秒过（em_vac Fmax=187），14 突变体 14/14 WT 身份匹配。WT NVT 已跑通,批量 MD 进行中（~3-4 天）。
- **解耦**：2B 进分类器靠静态接触特征,不依赖 MD;MD 是金标准验证/校准层。
- **下一步**:批量 MD 完 → `analyze_gromacs_md.py` 出 2B vs 2M 接触判定 → 校准 `AIM_RELEASE_2B_Z` 回填分类器（两线合流）。

### Changed (2026-06-11)

- **`scripts/pipeline/analyze_gromacs_md.py`** — autoinhib 分析从 **D'D3-A1 重写为 7A6O AIM-A1 体系**:修了多处 gemmi API bug(`.seq_num`、`element` 当原子名、`read_pdb` 读 `.gro`)使其**之前根本跑不通**;改为**轨迹级** AIM↔A1 接触数时间序列(trjconv 取帧 + gemmi NeighborSearch)、自动适配局部/规范残基编号、匹配 7A6O 批量 runner 的 `<variant>/md_7a6o/` 布局、并输出 **WT vs 2B vs 2M 判定**(`verdict_2b_vs_2m.csv`:每变体平衡态接触 + Δvs_WT,期望 2B 接触显著↓)。AIM/A1 范围、cutoff、skip、标签均可 CLI 覆盖。

**下一步（交接给推进的 agent）**：
1. 跑 `diagnose_clashes.py --variant-dir <WT autoinhib>`（确认 H-clash vs 重原子 + 挑最干净 model）→ `relax_autoinhib_structure.sh --variant VWF_WT`，贴回 Fmax 逐级。通则跑 autoinhib MD；不通则决定换 model / OpenMM / 暂走静态。
2. 跑 `extract_aim_autoinhib_features.py`（panel ~130 变体)→ 得到真 `aim_release_score`（现矩阵里是分不开的旧全局值)→ 用 known 2B/2M 校准 `AIM_RELEASE_2B_Z`，与 heparan LOF 联合，重测 2B recall（基线 2/12)。
3. 校准 `TWO_B_HOTSPOT_POS` / `LOF_COMBINED_Z` against 本地标签。

### Added (2026-06-10, 实验结构路线)

- **`scripts/pipeline/fetch_clean_7a6o.py`** — 下载并清理 VWF AIM-A1 自抑制态**实验晶体结构 PDB 7A6O**(X-ray 2.12 Å):删纳米抗体 VHH81/SO4/水,只留 VWF 链,报告残基范围/编号体系/缺失 loop/2B 热点覆盖。作为 MD 的可靠 WT 起点(力场只能局部松弛、修不了 Boltz 的错 pose;AIM-A1 也比 D'D3-A1 更对题 2B 机制)。
- **`scripts/pipeline/relax_autoinhib_structure.sh`** 加 `--pdb` 入口:直喂干净 PDB(如 7A6O)跑分级弛豫,跳过 Boltz CIF。
- **`scripts/pipeline/build_2b_mutants_foldx.py`** — 在 7A6O WT 骨架上用 FoldX(RepairPDB→BuildModel)批量造 2B/2M 点突变体 → 每个 `<variant>.pdb` 供 relax+MD。含**编号偏移自动探测(`--detect-offset`)+ WT 残基身份校验**(防改错残基)。把突变效应与预测噪声解耦。
- `docs/AUTOINHIB_MD_VALIDATION_GATES.md` 加 §0:优先走实验结构 7A6O,Boltz 路线降为交叉验证。

### Added (2026-06-10, MD 验证闸门)

- **`docs/AUTOINHIB_MD_VALIDATION_GATES.md`** — 上 NVT/批量前必过的 3 道闸:① clash 落不落在 D'D3-A1 界面(make-or-break);② 真空 EM 有没有挪动自抑制几何;③ 受控(带约束)平衡。含决策树。给 A40 Agent。
- **`scripts/pipeline/check_relax_distortion.py`** (闸门2) — gemmi CA 叠合,比较 Boltz 原始 vs 弛豫后结构,报全 CA-RMSD + 位移最大残基(对照 clash/界面判是否变形)。
- **`scripts/pipeline/run_autoinhib_md_from_relaxed.sh`** (闸门3) — 从 `relax_m*/em_10k.gro`+`topol.top` 受控续跑:cg 压 EM<1000 → 带 `-DPOSRES` 约束 NVT → 约束 NPT → 松开 production;GPU flags 按后端条件化;单变体,过了再上批量。
- **背景**: A40 实跑确认 WT model_2(18 重原子 clash, 非 H-clash)经分级弛豫 5.9e9→195(真空)→8.6e3(溶剂化松)可救,不必跳 OpenMM。但"能收敛≠读数可信",故加上述 3 闸。

### Changed (2026-06-10, 校准更新)

- **`scripts/agentic_vwf_classifier.py`** — 轴B(LOF→2M)从"单看 forced_binding"改为 **forced_binding + heparan 两轴联合**。校准 (known 2B n=39/2M n=47): fb 单用太弱(2M:2B 误伤比~1.5);heparan 较好(2B 中位+0.37 vs 2M −0.40);**联合 `mean(fb,heparan) ≤ −0.75` 最优 → 抓 30% 2M、仅误伤 5% 2B**。`LOF_COMBINED_Z=-0.75`,heparan 缺失时退到极保守单轴 `FB_LOSS_Z=-1.5`。8 项自测通过。机制: heparan 位点紧邻 GPIb 面,2M 破坏结合面会同时拖低两者。

### Changed (2026-06-10)

- **`scripts/agentic_vwf_classifier.py`** — RULE6(A1)重设计为**多轴方向判别**(2B=GOF/2M=LOF):轴B `fb_binding_zscore`↓→2M;轴A `aim_release_score`↑+结合保留→2B;轴C 临床 2B 热点先验;无方向信号→uncertain(不再硬判 2M)。修"用损伤大小推功能方向"的原理性错误。8 项功能自测通过;`FB_LOSS_Z`/`TWO_B_HOTSPOT_POS` 须校准。

### Added (2026-06-10)

- **`scripts/pipeline/diagnose_clashes.py`** — 全原子(含 H)steric clash 定位,排查 Boltz→GROMACS EM 跑飞(Fmax 5.9e9);扫 5 个 model 挑最干净。
- **`scripts/pipeline/relax_autoinhib_structure.sh`** — 分级弛豫(真空受约束 EM→真空无约束→溶剂化 EM),纯 CPU,解 Boltz 结构内部 clash,逐级报 Fmax。

### Added (2026-06-09)

- **`docs/A40_ENV_SETUP_REPORT_2026-06-09.md`** — 7×A40 独立服务器(无 `/lzy` NFS)的 GROMACS env 部署完整报告:从 git pull 到 SMOKE PASS 全流程、3 个 runner bug 修复、4 个新发现坑(curl 假断网、Test A 的 em.tpr 问题、初始结构需 EM、em.tpr vs md.tpr 命名冲突)、目标程序(autoinhib MD)恢复路径。

### Fixed (2026-06-09)

- **`scripts/pipeline/run_gromacs_vwf_md.sh`** (commit `ba5f711`) — A40 本地 env 兼容:
  - `gmx` 路径检测: 加 `$ROOT_DIR/envs/gromacs/` 本地 fallback (A40_AGENT_SETUP 路径)
  - `GMX_PY` 检测: 同样多路 fallback, 修 `$(dirname $GMX)/../../bin/python` `..` 数错 bug
  - `GMXDATA` 检测: 同样本地 fallback, 避免 `cd 错误目录`
  - `log()` 函数: 从 line 153 提前到 line 70, 修"A40 上 gmx 检测段调 log 报 '未找到命令'"的 bug (LZY 上因 `command -v gmx` 命中 PATH 而幸运避开)

### Known Issues (待办, 2026-06-09)

- `a40_selfcheck.sh` 联网检查用 `curl -sI`, 老版 curl (7.x) 报 `curl: (48)` 假断网, 实际 tuna 镜像可达。修法: 用 `python3 -c "import urllib.request; ..."` 或 `wget` 替代。
- `gpu_smoke_test.sh` Test A 用 `em.tpr` (steep 积分器) 跑 `-pme gpu` 必失败 (GROMACS 2025.x 不支持非动力学积分器的 PME-on-GPU)。这不是 GPU 问题, 是脚本设计 bug。修法: 内部用 md 积分器 tpr 跑 Test A。
- `output/boltz2_a1_dp_d3_results/` (autoinhib MD 输入) 在 A40 上不存在, 需 scp 自 LZY。

### Added (2026-06-07)

- **`scripts/pipeline/gpu_smoke_test.sh`** — 诚实的 GROMACS GPU 放行闸:用与生产 runner 相同的条件化 flags + md 积分器 tpr 实跑 5 步(朴素 `-nb gpu` smoke 会假性放行)。
- **`scripts/pipeline/extract_aim_autoinhib_features.py`** — 从 `a1_aim_autoinhibition_context` 结构抽 AIM↔A1 接触数 → `aim_release_score`(2B 判别特征,替代分不开 2B/2M 的全局 plddt)。
- **`docs/A40_RUNBOOK.md`** — 独立 A40 服务器(非 /lzy NFS)搬运 + 重建 CUDA gromacs env + smoke + autoinhib MD + 2B 特征流水线。

### Changed (2026-06-07)

- **`scripts/pipeline/run_gromacs_vwf_md.sh`** — GPU offload flags 按 `$GPU_BACKEND` 条件化:CUDA/SYCL→`-bonded gpu -update gpu`+CUDA graph;OpenCL/未知→仅 `-nb gpu -pme gpu`。修复 OpenCL 构建上 NVT/NPT/Prod 必 fatal(`not supported with OpenCL`)、白烧 GPU 机时的 bug。
- **`scripts/agentic_vwf_classifier.py`** — RULE6(A1)加性接入 `aim_release_score`:强松开→2B(早返回),弱松开→救回本判 2M 的 A1 默认;NaN 时退回原逻辑。阈值 `AIM_RELEASE_2B_Z` 暂定,须用已知 2B/2M 校准。

### Added

- **`scripts/pipeline/run_a1_gpiba_boltz2.sh`** — GPU batch runner for VWF A1 + GPIbα Boltz-2 predictions
  - 74 YAML inputs in `output/boltz2_a1_gpiba/`
  - 4-GPU parallel execution with per-job `.done` checkpointing
  - Key fixes: `CC=$(which gcc)`, `--num_workers 0`, `ulimit -c 0`
- **`scripts/pipeline/parse_a1_gpiba_results.py`** — Parser for Boltz-2 confidence outputs

### Infrastructure Fixes (2026-05-07)

| Issue | Root Cause | Fix |
|-------|-----------|-----|
| `core.*` files (6GB each) | DataLoader multiprocessing overflowed `/dev/shm` (64MB) | `--num_workers 0` disables sub-process IPC |
| `RuntimeError: Failed to find C compiler` | Triton JIT needs gcc; `CC` env not inherited by subprocesses | `CC=$(which gcc)` prepended to boltz command |
| `MisconfigurationException: No supported gpu backend` | `num_workers > 0` triggers DataLoader fork on CPU-only host | `--num_workers 0` |
| GPU util stuck at 4% then crashes | Same as above: multiprocessing + shm overflow → segfault | `--num_workers 0` |
| `boltz: command not found` in subshell | `conda activate boltz2` not finding env; shell uses wrong env | Use `export -f run_worker` + ensure `boltz` in PATH |

## [Prior History]

See `CLAUDE.md` and `DEV_PROTOCOL_STANDARD.md` for pre-May 2026 execution logs.

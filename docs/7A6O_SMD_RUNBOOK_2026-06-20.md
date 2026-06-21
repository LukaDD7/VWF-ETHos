# 7A6O AIM-unfolding Steered MD (SMD) Runbook — 2026-06-20

目标：用恒速 SMD 测 **AIM 解折叠力**，给分类器补一条**力依赖的 2B 阳性轴**
（平衡 MD 看不到的那条）。机制：type 2B (GOF) **降低 AIM 解折叠所需的力**
(Nat Commun 2021)，故 **解折叠力低 → 偏 2B**。与平衡态 `md_face_destab_score`
(2M/LOF 轴) 正交互补。

## 0. 为什么要 SMD（先读这段）
- 系统核心瓶颈是 **2B 召回仅 17%（12 个里 8 个被误判 2M）**，是个 *2B 抓不住* 的问题。
- 平衡 50 ns MD 给的是 **2M/LOF** 证据，**不解决 2B 漏判**。2B 的 GOF 是**力依赖**的，
  必须施力（SMD）才显现。**这条 SMD 轴才是直击瓶颈的。**

## 1. 前置条件（A40 上）
- gromacs env 已重建且 `gpu_smoke_test.sh` PASS（见 `docs/A40_RUNBOOK.md`）。
- 每个变体已跑过平衡管线，存在：
  - `output/gromacs_md_autoinhib/<variant>/relax_pdb/topol.top` (+ `posre.itp`, chain `.itp`)
  - `output/gromacs_md_autoinhib/<variant>/md_7a6o/md_prod*.gro`（平衡产物；
    没有也会回退到 `md_data/7a6o_reference_md/variants/<variant>/final.gro`）
- env 变量：`export GMX=$PROJ/envs/gromacs/bin.AVX2_256/gmx`,
  `export GMXLIB=$PROJ/force_fields`（patched charmm36m）。

## 2. 几何与参数（已用 WT 实测定标）
- 构件原生编号 1262–1466。**N-anchor = 1262–1264 CA，C-anchor = 1464–1466 CA**
  （两端起始仅相距 2.0 nm，折叠回 A1 同一面 = 不连续模块拓扑）。
- 拉伸坐标：把两 anchor 沿 **z** 拉开 → AIM 剥离/解折叠（复现 Springer 单分子实验）。
- 盒子：先把 N→C 向量转到 +z，再开 **8.2 × 8.2 × 11.8 nm** 长盒（~80k 原子，原 35k 的~2.3×）。
- 默认拉伸：`PULL_NM=5.0 nm`，`RATE=0.001 nm/ps (1 nm/ns)`，`K=1000 kJ/mol/nm²`
  → 每副本 **5 ns / 2.5M 步**。A40 单卡约 2–4 h/副本。
- Pull geometry 使用 `direction-periodic` 沿 +z 拉伸，而不是 `distance`；否则拉到接近半盒长时 GROMACS 会报
  `Distance between pull groups ... larger than 0.49 times the box size`。
- prep 默认 `PREP_NTOMP=16`，避免 EM/NVT/NPT 独占 128 CPU 线程；可按机器负载调整。

## 3. 跑法

### 3.1 先做 3 个对照确认"力能分开"（强烈建议先跑这步再铺全量）
> ⚠ 2026-06-21 实测：默认 `RATE=0.001`(1 nm/ns) 下 2B/2M 力**反号**(R1306W 1117 pN >
> R1374H 993 pN)，因 ~1000 pN 由黏滞摩擦主导而非自抑制能垒。**标定务必用慢速**。
> WT 之前"变慢"是 CPU 超订 + 缺 `-bonded gpu`（已修），**非 WT 拓扑问题**
> （见 `docs/7A6O_SMD_WT_DIAGNOSIS_2026-06-21.md`）。
```bash
cd $PROJ
for v in WT R1306W R1374H; do
  bash scripts/pipeline/prep_7a6o_smd.sh $v 0               # 建盒+溶剂+EM+NVT/NPT
  RATE=0.00025 bash scripts/pipeline/run_7a6o_smd.sh $v 0 5 # 慢速 0.25 nm/ns, 5 副本
done
python3 scripts/pipeline/analyze_7a6o_smd.py
# 期望: WT 力最高, R1306W(2B) 明显更低, R1374H(2M) ≈ WT/居中, 且 AUC(2B 力<2M) > 0.5。
# 仍反号 → 换反应坐标(AIM vs A1 核心脱离)或退回平衡轴, 不强上 SMD 力轴。
```

### 3.2 全量参考集（确认分得开后）—— 多 GPU 并行
变体清单（与平衡特征同一套标签，便于联合校准）：
```
WT  R1306W R1306Q R1308C I1309V S1310F W1313C V1314F V1316M   # WT + 8×2B
R1374C R1374H G1324S                                          # 3×2M
P1337L R1341Q R1341W                                          # 3×待定 (测试样本)
```
每卡一个变体的 prep→smd 队列（沿用项目 setsid 习惯，避免 SSH 断开被杀）：
```bash
gpu=0
for v in WT R1306W R1306Q R1308C I1309V S1310F W1313C V1314F; do
  setsid bash -c "bash scripts/pipeline/prep_7a6o_smd.sh $v $gpu && \
                  bash scripts/pipeline/run_7a6o_smd.sh $v $gpu 3" \
     </dev/null > output/gromacs_md_autoinhib/smd_gpu${gpu}_${v}.log 2>&1 &
  gpu=$(( (gpu+1) % 7 ))
done
# 余下 7 个变体排到空出来的卡上（或排成每卡的第二个任务）。
```
> 算力：15 变体 × 3 副本 × ~3 h ≈ 135 GPU·h；7 卡并行约 1 天。先跑 §3.1 省力。

### 3.3 分析 → 特征
```bash
python3 scripts/pipeline/analyze_7a6o_smd.py \
    --input output/gromacs_md_autoinhib --output output/md_7a6o_smd_features.csv
```
产出 `output/md_7a6o_smd_features.csv`：
- `smd_unfold_force_pN`（rupture force，均值±sd，pN）
- `smd_force_z`（对 WT+2B+2M 参考集 z 标定）
- `smd_2b_score = -z`（**越高 = 越易解折叠 = 越像 2B**）

把结果（CSV）push 回仓，我据此：① 验证 2B vs 2M/WT 力分离与 AUC；
② 若分得开，将 `smd_2b_score` 作 RULE6 的 **2B 阳性轴**接入分类器并校准阈值
（这是目前唯一能提升 2B 召回的轴）。

## 4. 排错
| 现象 | 处理 |
|---|---|
| `grompp` tc-grps 找不到 Protein | index.ndx 已含默认组(make_ndx)+anchor，检查 `grep '\[' smd/index.ndx` |
| EM/NVT LINCS 崩 | EM 没收敛；查 `prep_em.log` Fmax，必要时 emtol 调 200 |
| `Distance between pull groups ... larger than 0.49 times the box size` | 确认 `run_7a6o_smd.sh` 使用 `pull-coord1-geometry = direction-periodic`；旧 `distance` geometry 会在 5 nm 拉伸中失败 |
| 拉伸中蛋白撞周期镜像 | 加大 `BOX_Z`（env 变量）或减小 `PULL_NM` |
| 力曲线无峰（单调升） | RATE 太快，AIM 没来得及协同解折叠；降 `RATE=0.00025` |
| SMD 生产巨慢（~1 ns/day） | CPU 超订(prep×NTOMP=16 并发)+pull 强制 CPU update。已加 `-bonded gpu`；prep 降 `PREP_NTOMP=8` 或与 SMD 错峰 |
| 2B/2M 力反号或分不开 | RATE 太快(速率主导)。用 `RATE=0.00025`+≥5 reps 标定；看 `analyze` 打印的 AUC |
| GPU 0% util 数分钟 | 正常（首次 kernel 编译）；持续 0% 查 `md_smd_rep*.stdout` |

## 5. 与现有体系的关系
- 平衡轴 `md_face_destab_score`（已并入，2M/LOF 确证）+ SMD 轴 `smd_2b_score`（待并入，2B GOF）
  = A1 域 2B/2M **方向判别的两条正交力学证据**，分别覆盖 GOF/LOF。
- 见 `docs/7A6O_MD_FEATURE_ANALYSIS_2026-06-20.md`（平衡轴）与本文（力轴）。

# 7A6O SMD：WT "pipeline 问题" 诊断与修复 — 2026-06-21

回应 `docs/7A6O_SMD_STATUS_HANDOFF_2026-06-21.md`（Codex）。结论：**WT 不是
provenance（来源/拓扑）问题，是性能/调度问题**；而 2B/2M 力反号是**拉伸速率**
问题，与 WT 无关。下面给证据、根因、已实施的修复，和标定前的 go/no-go。

## 1. WT ≠ provenance 不匹配（有拓扑证据）

直接对比已入库参考集里的历史 WT 与 FoldX 突变体（`md_data/7a6o_reference_md`）：

| 变体 | 构件 | 端基 patch | 蛋白原子数 |
|---|---|---|---|
| WT | ILE1…PRO205 (205 res) | N: NH3+ (H1/H2/H3) / C: COO- (OT1/OT2) | 3349 |
| R1306W (2B) | ILE1…PRO205 | 同上 | 3349 |
| R1374H (2M) | ILE1…PRO205 | 同上 | 3342 (R→H) |
| G1324S (2M) | ILE1…PRO205 | 同上 | 3353 |

构件、残基编号、**端基 patch 完全一致**，蛋白原子数只随突变残基变化。
→ **WT 与 FoldX 突变体结构上可直接比较，无需为 provenance 重建 WT。**
（若要保险控制，FoldX 的 `WT_*.pdb`（BuildModel 优化 WT）走同一 relax 管线即匹配
基线，但数据表明非必需。）

## 2. WT 变慢的真实根因（~1.4 ns/day）

两点叠加，都在运行配置里，与 WT 本身无关：

1. **SMD 生产 mdrun 缺 `-bonded gpu`**。且 **COM pulling 在 GROMACS 中强制把
   积分/update 步放到 CPU**（`-update gpu` 会被拒），CPU 成瓶颈。
2. **CPU 超订**：7 个 prep 任务并发 × `PREP_NTOMP=16` = 112 线程，叠加 CPU 端
   SMD 积分 → 谁在 prep 高峰期跑谁就爬。WT 恰好是那个（rep1, GPU0）。

同等争用下任何突变体都会一样慢；反之 WT 在不拥挤的机器上 + `-bonded gpu`
会跑到突变体的速度。

## 3. 真正的科学风险：力反号（与 WT 无关，修 WT 也救不了）

探针：R1306W(2B) **1117 pN** > R1374H(2M) **993 pN**，与"2B 更易解折叠（力更低）"相反。
根因几乎确定是**拉伸速率 1 nm/ns 太快**：rupture force ~1000 pN，是实验（~10–20 pN）的
~50–100×，由**黏滞摩擦主导**而非自抑制能垒 → 区分 2B 的那点能垒差被掩盖甚至反号。

## 4. 已实施的修复（本次提交）

- **`run_7a6o_smd.sh`**：SMD 生产 mdrun 加 `-bonded gpu`（释放 CPU；附注释说明
  pull 强制 CPU update，故不用 `-update gpu`）。加**速率主导警示**：`RATE ≥ 0.001`
  时打印提醒，建议标定用 `RATE=0.00025`。
- **`analyze_7a6o_smd.py`**：除峰值力外加**基于功 (work) 的轴**（`smd_work_z` /
  `smd_2b_score_work`，对单帧尖峰更稳健）；并打印 **AUC(2B 力 < 2M)** 作分离自检
  （≤0.5 = 无/反向分离 → 别接分类器，先调速率）。

## 5. 标定前 go/no-go（强烈建议先做，再决定是否铺全量）

```bash
# 慢速探针: 0.25 nm/ns, 5 reps, 仅 3 个对照
for v in WT R1306W R1374H; do
  RATE=0.00025 bash scripts/pipeline/run_7a6o_smd.sh $v <gpu> 5
done
python3 scripts/pipeline/analyze_7a6o_smd.py
```
判据：**WT 力最高，R1306W(2B) 明显更低，R1374H(2M) ≈ WT 或居中**，且
`AUC(2B 力<2M)` 明显 > 0.5。满足才铺全量 15×（每个 ≥3 reps）。
- 若慢速仍分不开/反号 → 换反应坐标（AIM 整体 vs A1 核心的脱离，位移更小、力更低、
  更贴"自抑制释放"语义），或退回平衡轴 + 临床热点，不强上 SMD 力轴。

## 6. 运行纪律（避免再现 WT 爬行）

- SMD 生产与 prep **不要同时**在 7 卡上各开 `NTOMP=16`；prep 队列降到 `PREP_NTOMP=8`
  或与 SMD 错峰。
- 每个 SMD 生产给足专属 CPU（`NTOMP=8` + 不重叠 `-pinoffset`，脚本已按 GPU 分配）。

## 7. 待办（数据回来后由分类器侧做）

- 慢速探针 CSV push 回 → 我算 2B/2M/WT 力与功的分离 + AUC；
- 若分得开，把 `smd_2b_score`（或 `_work`）作 RULE6 的 **2B 阳性轴**接入并校准
  （目前唯一能抬升 2B 召回 17% 的轴）。与平衡轴 `md_face_destab_score`(2M/LOF) 正交。

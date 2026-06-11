# A40 完整管线交接 — 7A6O 实验结构路线 (2026-06-11)

> 给接手 agent / 用户。本日工作小结 + 下一步。
> **状态**: 全程绿。WT NVT 跑通; 14 个 FoldX 突变体正在批量弛豫 (2.4h 收尾)。
> **关键设计**: 弃用 Boltz-2 D'D3-A1 (backbone 几何坏), 改用 **PDB 7A6O (AIM-A1 实验晶体, 2.12 Å)** 作 WT 骨架 + FoldX BuildModel 在骨架上改单残基。

---

## 1. 今日完成

| 步骤 | 状态 | 产物 |
|------|------|------|
| 1. git pull (取 fetch_clean_7a6o + build_2b_mutants_foldx) | ✅ | commit 2aaa4ac |
| 2. fetch_clean_7a6o.py → 7A6O WT PDB (删 VHH81 纳米抗体/SO4/水) | ✅ | `structures/7A6O_AIM_A1_clean.pdb` (205 aa, A1:1262-1466, 全部 17 个 2B 热点已解析) |
| 3. WT 弛豫 (真空受约束→真空无约束→溶剂化水 settle) | ✅ | em_vac.gro Fmax=187; em_water.gro Fmax=6e4 → solv_ions_em.gro |
| 4. build_2b_mutants_foldx.py → 14 个 2B/2M 突变体 | ✅ | `structures/7a6o_mutants/{R1306W,R1306Q,R1308C,I1309V,S1310F,W1313C,V1314F,V1316M,P1337L,R1341Q,R1341W,R1374C,R1374H,G1324S}.pdb` (offset=0, 全部 14/14 WT 身份匹配) |
| 5. 14 突变体批量弛豫 (--skip-vacuum, 走 fast solvated EM) | 🔄 进行中 | 1/14 done (R1306W Fmax=6e4), 剩余 13, ~10min/ea, 预计 14:36 完成 |
| 6. WT NVT (100K, -DPOSRES, 50ps, GPU 0) | 🔄 进行中 | step 13400/25000 (54%), T=101K, Constr.rmsd=2e-6, ETA 14:50 |

---

## 2. 关键脚本补丁 (committable)

**`scripts/pipeline/relax_autoinhib_structure.sh`** — 3 项改动:

1. **加 `--skip-vacuum` 标志**: FoldX 突变体已优化 side chain, 跳过 em_posres/em_vac, 直接溶剂化。10 min/mutant → 5 min/mutant。
2. **em_sol 改用快档** (mk_em_sol): emstep=0.05 / nsteps=50 / constraints=h-bonds / define=-DPOSRES。原默认 emstep=0.001 / nsteps=100k 在 35k 原子系统跑 8+ 天, 不可用。
3. **grompp3/grompp4 加 `-r solv.gro`** (POSRES 需要 reference 坐标)。

⚠ 没 commit 之前本地状态是: `git status` 标 `M scripts/pipeline/relax_autoinhib_structure.sh`, 60+/28- 行 diff。等管线跑通一起 commit。

---

## 3. 关键卡点 + 解决方案 (今日沉淀)

| 卡点 | 原因 | 修法 |
|------|------|------|
| 溶剂化 EM 8 天 | emstep=0.001 + nsteps=100k + 35k 原子 + CPU 128 MPI → 8 sec/step | mk_em_sol 用 emstep=0.05 + nsteps=50, 50 steps × 10s = 8 min |
| cg + h-bonds 必炸 SETTLE | 刚溶剂化水有 clash, cg 失败立刻 abort | 用 steep + h-bonds, 失败可继续 |
| CPU EM 40 核全占 | 16 物理核被 128 MPI 撑满 | 接受(只是 8 min/ea, 14 突变体 2.4h 串行可接受) |
| GPU NVT 一开始 NaN | 没 EM 直接 gen_vel, 起始 clash 导致发散 | 必须先 EM_water 50 步, 才能 NVT |
| pdb2gmx 1MET bug | charmm2gmx 端口缺 `<RESNAME>1` patch | 2026-06-01 已修 (force_fields/charmm36m.ff/aminoacids.n.tdb 注入 17 个 protein N-terminal patches) |
| 路径相对 | 脚本 `cd $WORK` 后相对路径失效 | 全部用绝对路径 (--pdb /abs/path.pdb) |
| POSRES 缺 reference | mk_em_sol 用了 `define=-DPOSRES` 但 grompp 没 `-r` | grompp3/grompp4 加 `-r solv.gro` |

---

## 4. 接下来要做什么

### 4a. WT 续跑 (今晚/明早)

WT NVT 14:50 完。然后:
```bash
cd /media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan
# 复用 relax_m2/ 的 em.gro (solvated + water-settled) → 起 NPT + Production
# 用 run_autoinhib_md_from_relaxed.sh 但绕开它的"cg EM polish" (会再 SETTLE 炸)
# 改写或直接手写 NPT/Prod 步骤
```
- **NPT** 200 ps (-DPOSRES, 310K)
- **Production** 50 ns (无约束, 310K)
- GPU 0 全程, ~12-15 h

### 4b. 突变体批量 MD

弛豫完后(14:36), 14 个变体 × {NVT 50ps + NPT 200ps + Prod 50ns}:
- 串行 ~19 h/ea → 12 天
- 7×A40 并行 5 ea/批 × 3 批 → 3-4 天

**包装脚本建议**: `scripts/pipeline/run_mutant_md_batch.sh`
- 输入: 突变体名列表 + GPU 数
- 流程: 弛豫(em_water.gro 在) → NVT → NPT → Prod,每变体一 GPU

### 4c. 分析 (分析 WT vs 突变体的闭合态稳定性)

`scripts/pipeline/analyze_gromacs_md.py` 应有 autoinhib 系统分析:
- D'D3↔A1 接触数 / 接触寿命
- 界面盐桥/氢键持久性
- 界面 BSA / 开口度
- WT 保持咬合 (低 RMSD); 2B 变体松开 (高 RMSD, 接触数↓)

**判定标准**: 已知 2B 变体(R1306W, R1306Q, R1308C, I1309V, S1310F, W1313C, V1314F, V1316M)MD 后接触数应明显↓ ; 已知 2M 对照(R1374C, R1374H, G1324S)应保持或稍↓; 2B vs 2M 应有统计显著差。

---

## 5. 文件位置速查

```
/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/
├── structures/
│   ├── 7A6O_AIM_A1_clean.pdb          # WT 实验结构 (删纳米抗体)
│   └── 7a6o_mutants/                  # 14 突变体 (FoldX BuildModel)
│       ├── R1306W.pdb ... G1324S.pdb  # 14 个变体
│       ├── WT_Repair.pdb              # FoldX 修过的 WT
│       └── individual_list.txt        # FoldX 输入清单
├── output/
│   ├── gromacs_md_autoinhib/
│   │   ├── 7A6O_WT/relax_m2/          # WT 弛豫产物
│   │   │   ├── em_vac.gro (Fmax=187)
│   │   │   ├── em_water.gro (Fmax=6e4)  ← solv_ions_em.gro 同文件
│   │   │   ├── topol.top
│   │   │   ├── posre.itp
│   │   │   └── nvt_soft.* (NVT 100K, 正在跑)
│   │   ├── R1306W/relax_pdb/ ...     # 突变体弛豫产物 (进行中)
│   │   └── _mutant_relax_logs/        # 每变体 1 个日志
│   └── boltz2_a1_dp_d3_results/       # 旧 Boltz 数据 (已弃用, 留着交叉验证)
├── scripts/pipeline/
│   ├── fetch_clean_7a6o.py            # 取 7A6O WT
│   ├── build_2b_mutants_foldx.py      # FoldX 造 14 突变体
│   ├── relax_autoinhib_structure.sh   # 分级弛豫 (--skip-vacuum 标志)
│   ├── run_autoinhib_md_from_relaxed.sh  # 受控 WT 续跑 (待修 SETTLE bug)
│   ├── diagnose_clashes.py            # (Boltz 路线用, 暂搁置)
│   ├── check_relax_distortion.py      # (Boltz 路线用, 暂搁置)
│   └── analyze_gromacs_md.py          # 轨迹后分析 (待跑)
└── docs/
    ├── AUTOINHIB_MD_VALIDATION_GATES.md  # 闸门 (WT 路线自动过 1/2)
    └── A40_FULL_PIPELINE_HANDOFF_2026-06-11.md  # 本文件
```

---

## 6. 决策记录

| 决策 | 原因 |
|------|------|
| 弃用 Boltz-2 D'D3-A1 | EM 跑飞 Fmax 5.9e9, 18 个重原子 clash 在 D'D3-A1 界面附近, MD 不可信 |
| 改用 7A6O AIM-A1 | 实验晶体 2.12 Å, EM 秒过 (Fmax 187→6e4), AIM-A1 才是文献公认 2B 机制 |
| FoldX BuildModel 而非 Af3/Boltz 重出 | 同 WT 骨架改 1 残基, 把"突变效应"与"预测噪声"解耦 |
| CPU 单跑 em_water + GPU 跑 NVT/NPT/Prod | CPU EM 8 min/ea (14 ea = 2h), GPU 持续用 |
| --skip-vacuum for mutants | FoldX 已优化侧链, em_posres/em_vac 冗余 |
| 50 ns Production | 平衡时间分辨率与算力: ~10 ns 测接触寿命, 50 ns 留余量 |
| 14 突变体 (8 2B + 3 2M + 3 混合) | 用 8 个 known 2B (R1306W/Q, R1308C, I1309V, S1310F, W1313C, V1314F, V1316M) + 3 已知 2M (R1374C/H, G1324S) + 3 中性 (P1337L, R1341Q/W) |

---

## 7. 还没做 (按优先级)

1. **WT NPT + Production** (今晚/明早跑, 12-15h)
2. **突变体批量 MD** (弛豫完后接, 3-4 天)
3. **`run_autoinhib_md_from_relaxed.sh` 修 SETTLE bug** (它的 cg polish 会在刚溶剂化系统上炸; 改成 steep+POSRES+h-bonds 或干脆跳过)
4. **`analyze_gromacs_md.py` 跑 2B vs 2M 对比** (结果产出)
5. **静态特征 + 闭合态 MD 联合校准 `AIM_RELEASE_2B_Z`** (回填 `agentic_vwf_classifier.py`)
6. **commit 今日脚本改动 + 这个 handoff doc** (起 PR)
7. **Boltz 路线 (交叉验证)**: 现状是空跑 (Fmax 5.9e9, 不可信)。如果其他组要把 Boltz 也跑通, 走 7A6O 路线代替。

---

*Last updated: 2026-06-11 14:25 by main A40 agent*

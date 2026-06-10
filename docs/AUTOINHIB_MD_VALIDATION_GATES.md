# Autoinhib MD 验证闸门 — 上 NVT/批量前必须搞清楚的 3 个核心问题

> 给 A40 Agent。WT model_2 的弛豫已把 Fmax 从 5.9e9 降到可控(em_10k.gro),**工程上**可救。
> 但"EM 能收敛"≠"MD 读数可信"。**在跑 NVT/NPT/production 前,先过下面 3 道闸**;
> 任何一道不过,这套 autoinhib MD 测出来的"自抑制松开"都可能是假象,会误导 2B 分类。
>
> 背景:这套 MD 的唯一目的,是量"**A1 自抑制接口在 310K 下稳不稳/会不会松开**"
> (2B = 松开 = GOF)。所以**接口区域的结构保真度**是一切。下面三道闸都围绕这个。

---

## 0. 推荐:优先改用实验结构 7A6O(AIM-A1),而不是跟 Boltz 死磕

VWF AIM-A1 自抑制态有**实验晶体结构 PDB 7A6O**(X-ray 2.12 Å)——真实坐标、clash 极少、
EM 一下就过,且 AIM-A1 正是文献公认的 2B 自抑制机制(比 Boltz 的 D'D3-A1 构建更对题)。
**用它当 WT 骨架,远比 de novo 预测可靠**(力场只能局部松弛、修不了错 pose)。

```bash
python3 scripts/pipeline/fetch_clean_7a6o.py          # 下载+删纳米抗体/水/SO4 → 干净 WT PDB
bash scripts/pipeline/relax_autoinhib_structure.sh --pdb structures/7A6O_AIM_A1_clean.pdb --variant 7A6O_WT
bash scripts/pipeline/run_autoinhib_md_from_relaxed.sh --variant 7A6O_WT --model pdb   # 见脚本注释调整
```
突变体(在 WT 骨架上改单残基, 把"突变效应"与"预测噪声"解耦):
```bash
python3 scripts/pipeline/build_2b_mutants_foldx.py \
    --wt structures/7A6O_AIM_A1_clean.pdb --foldx <foldx 二进制> --detect-offset
# → structures/7a6o_mutants/<variant>.pdb; 再各自 relax+MD
for v in structures/7a6o_mutants/*.pdb; do
  bash scripts/pipeline/relax_autoinhib_structure.sh --pdb "$v" --variant "$(basename "${v%.pdb}")"
done
```
WT vs 突变体闭合态稳定性之差 = 2B 信号。

> 走实验结构路线时,下面的闸门 1/2(针对 Boltz 结构的 clash/形变)基本自动通过;
> 仍建议跑闸门 3 的受控平衡。Boltz 路线作为交叉验证保留。

---

## 闸门 1 —— 18 个重原子 clash 在哪?(make-or-break)

`diagnose --variant-dir` 只给了计数(model_2: 18 个全重原子, min 0.75 Å)。**必须细查位置**:

```bash
envs/gromacs/bin/python scripts/pipeline/diagnose_clashes.py \
  --input output/boltz2_a1_dp_d3_results/boltz_results_VWF_WT_dp_d3_a1/predictions/VWF_WT_dp_d3_a1/VWF_WT_dp_d3_a1_model_2.cif \
  --top 20
```

看输出里每个 clash 对的**残基号**,对照构建里 D'D3 段 vs A1 段的边界(见 `generate_a1_dp_d3_yamls.py` 里的序列拼接):

| clash 落在 | 判定 |
|---|---|
| **外周/柔性 loop、单个结构域内部**(远离 D'D3↔A1 界面) | ✅ **过**。EM 抹平无害,MD 读数有效。 |
| **D'D3 与 A1 的接触界面**(我们要量自抑制的那块) | 🔴 **不过**。Boltz 恰在测量区不可靠,EM"修好"只是极小化假象 → 闭合态稳定性=噪声。**换策略**(见底部决策树),别跑这套 MD。 |

> 怎么判"接口":clash 对的两个残基若**一个在 D'D3 段、一个在 A1 段**,就是跨界面 clash → 危险。
> 同段内的 clash 一般是局部 rotamer/loop,问题不大。

---

## 闸门 2 —— 真空 EM 有没有把自抑制几何搞变形?

真空极小化(尤其无约束那 556 步)会塌缩表面盐桥、**挪动 D'D3 与 A1 的相对取向**——而这正是我们要测的量。比较弛豫前后:

```bash
envs/gromacs/bin/python scripts/pipeline/check_relax_distortion.py \
  --orig output/boltz2_a1_dp_d3_results/boltz_results_VWF_WT_dp_d3_a1/predictions/VWF_WT_dp_d3_a1/VWF_WT_dp_d3_a1_model_2.cif \
  --relaxed output/gromacs_md_autoinhib/VWF_WT/relax_m2/em_vac.gro
# (em_vac.gro 是真空弛豫后的纯蛋白; 若 gemmi 读不了 .gro, 先
#  gmx editconf -f em_vac.gro -o em_vac.pdb 再用 --relaxed em_vac.pdb)
```

脚本报:全蛋白 CA-RMSD + 位移最大的 10 个残基。

| 结果 | 判定 |
|---|---|
| 全 CA-RMSD 小(< ~2 Å)且**最大位移残基不在 D'D3-A1 界面** | ✅ 过。弛豫只动了局部,自抑制几何保住。 |
| RMSD 大,或**最大位移残基就在界面** | ⚠️ 真空步动了自抑制 pose。**改用全程带位置约束的弛豫**(`relax_autoinhib_structure.sh --no-solvate` 只到真空受约束段, 或下面续跑脚本里 EM 也加 -DPOSRES),重做闸门 2。 |

> 把闸门 2 的"最大位移残基"和闸门 1 的"clash 残基"对一下:如果同一片区域既有 clash 又位移最大,那块就是 Boltz 没学好 + 被 EM 强行掰过的"假结构",尤其要警惕它在不在界面。

---

## 闸门 3 —— 受控平衡(别从松收敛结构直接无约束升温)

`em_10k.gro` 是 `emtol=10000` 的**松**收敛(Fmax 8638)。直接无约束 NVT 容易炸。正确做法已封装:

```bash
bash scripts/pipeline/run_autoinhib_md_from_relaxed.sh --variant VWF_WT --model 2
```

它做:① 把 EM 再压到 Fmax<1000(cg)② **带位置约束的 NVT**(-DPOSRES)③ **带约束的 NPT** ④ 松开约束 ⑤ production。GPU flags 按后端自动选(CUDA 全常驻)。**只跑 WT,过了再上 5 变体批量。**

判 NVT/NPT 是否健康:温度稳定在 310K、压力/密度收敛、无 LINCS warning 爆炸、势能平台。

---

## 跑完才是真正的读数(闸门都过之后)

production 轨迹上量"自抑制松开"(2B 信号):

- D'D3↔A1 **界面接触数 / 接触寿命**随时间(`analyze_gromacs_md.py --system autoinhib`)
- 界面盐桥/氢键持久性、界面 **BSA / 开口度**
- WT 应保持咬合;2B 变体(后续批量)应自发松开

---

## 决策树(给 Agent)

```
闸门1: clash 在界面?
  ├─ 是 → 🔴 这套 MD 不可信。换策略:
  │        a) 试别的 model 的 --input(也许某 model 界面干净)
  │        b) 若所有 model 界面都坏 → Boltz 对 D'D3-A1 界面预测能力不足,
  │           退回"静态特征 + 不依赖该 MD"; 或考虑 AF3 / 更长采样重出结构
  └─ 否 → 闸门2: 真空弛豫动了界面?
            ├─ 是 → 改全程受约束弛豫, 重做闸门2
            └─ 否 → 闸门3: run_autoinhib_md_from_relaxed.sh (WT-only)
                      ├─ NVT/NPT 不健康 → 查残余 clash / 再压 EM
                      └─ 健康 → production → analyze → (过了再上 5 变体)
```

**核心一句话**:工程上"能跑"已经验证;现在卡的是"**跑出来的东西在不在测量区是真的**"。闸门 1 是 make-or-break,先做它。

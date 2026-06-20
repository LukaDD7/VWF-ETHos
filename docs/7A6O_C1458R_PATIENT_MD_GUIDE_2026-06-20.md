# 7A6O AIM-A1 MD — 补齐患者 C1458R 指引 (2026-06-20)

本机(CPU 侧)写。目标:把**患者本人变异 C1458R** 补进 7A6O AIM-A1 自抑制 MD
参考盘,并给出可写进医生报告的、有误差棒的 AIM↔A1 接触数 delta。

> 背景诊断:患者 VWF c.4372T>C, p.Cys1458Arg(p.C1458R),exon28,杂合 de novo
> (父母野生型),无数据库记录。C1458 是 **C1272–C1458 长程二硫键**的一端,且紧贴
> C_AIM 段(1460–1466)。boltz2 功能盘已给混合信号(heparan/AIM 偏 2B,GPIbα 强迫
> 结合**不升反降**,不像经典 2B GOF)。这套 MD 是要从自抑制几何上再加一条独立证据。

---

## 0. 现状与最大缺口

现有 MD 参考盘(`analyze_7a6o_completed_md.py` 的 `DEFAULT_VARIANTS`):
WT + 8 个 2B(R1306W/Q, R1308C, I1309V, S1310F, W1313C, V1314F, V1316M)
+ 3 个 2M(R1374C/H, G1324S)+ 3 个未定 `?`(P1337L, R1341Q, R1341W)。

**缺口(按价值排序):**
1. **C1458R 本人完全没跑** —— 这是唯一直接答诊断题的一条。最高优先。
2. **无重复** —— 每个变异只有 1 条 50ns = 1 个样本,接触数 delta 没有误差棒。
   至少 WT 和 C1458R 各做 **3 条不同 seed 重复**,才能写进报告。
3. **2M 臂太薄**(n=3,Cohen's d 不可信);3 个 `?` 没定标。次优先,这次可不动。

---

## 1. 备料(CPU,不计费 GPU)—— 你上 GPU 前最好本机就做完

### 1.1 FoldX 建 C1458R 突变体
```bash
python3 scripts/pipeline/build_2b_mutants_foldx.py \
    --wt structures/7A6O_AIM_A1_clean.pdb \   # 用与那 14 个相同的干净 WT 骨架
    --foldx /path/to/foldx --detect-offset \
    --variants C1458R
```
**验证闸门 A(必看,否则全错):**
- 脚本用 gemmi 校验 WT 残基身份。**确认 1458 在 7A6O 里被解析且 WT=CYS**。
  若打印 "WT 身份对不上/跳过 C1458R" → 说明 7A6O 该位点未解析或编号偏移没对齐,
  **停**,先确认 offset / 结构覆盖范围,不要硬跑。
- C1458 与 C1272 在 WT 里成二硫键。FoldX 把 1458 换成 Arg 后,这对二硫键必须**断开**。

### 1.2 弛豫 + 溶剂化(FoldX 突变体走 `--skip-vacuum`)
```bash
bash scripts/pipeline/relax_autoinhib_structure.sh \
    --variant C1458R --pdb output/.../C1458R.pdb --skip-vacuum
# 产物: output/gromacs_md_autoinhib/C1458R/relax_pdb/{solv_ions_em.gro,topol.top,posre.itp}
```
**验证闸门 B(本研究最关键的科学闸门):**
- 打开 `relax_pdb/pdb2gmx.log` / `topol.top`,确认 **C1272–C1458 二硫键已不存在**
  (Cys1272 应为自由半胱氨酸,topology 里没有连到 1458 的 SS bond)。
  pdb2gmx 的 SS 自动识别若仍把 1272 错配到别处,要人工处理。
  **二硫键残留 = 模拟的是错的分子 = 人为错误,绝对不能放过。**
- 溶剂化 EM 收敛(脚本结尾应打印 `✅ ... <1e4`)。

### 1.3 refined EM(降局部力,防 LINCS/CUDA 崩)
```bash
VARIANTS=C1458R bash scripts/pipeline/refine_7a6o_mutant_em.sh
# 产物: relax_pdb/solv_ions_em_refined.gro;日志末尾 Maximum force 应明显下降
```

### 1.4 重复体备料(WT 和 C1458R 各 3 条)
runner 的工作目录由 variant 名决定,且 NVT 用 `gen_seed = -1`(每次随机速度),
所以**用不同 variant 名指向同一份 refined 起点**即可得到独立重复:
```bash
for base in WT C1458R; do
  for r in r1 r2 r3; do
    dst=output/gromacs_md_autoinhib/${base}_${r}/relax_pdb
    mkdir -p "$dst"
    cp output/gromacs_md_autoinhib/${base}/relax_pdb/{solv_ions_em_refined.gro,solv_ions_em.gro,topol.top,posre.itp} "$dst"/ 2>/dev/null
  done
done
```
(WT 已有一条原始轨迹,可作 r0;再加 r1–r3。C1458R 主轨 + r1–r3。)

---

## 2. 跑 production(GPU,计费区)

每条用 `setsid` 起,survive SSH 退出。一张卡一条,pinoffset 由 GPU id 决定(见 runner)。
```bash
setsid bash scripts/pipeline/run_7a6o_variant_direct.sh C1458R     0 < /dev/null > output/gromacs_md_autoinhib/C1458R_$(date +%Y%m%d_%H%M).log 2>&1 &
setsid bash scripts/pipeline/run_7a6o_variant_direct.sh C1458R_r1  1 < /dev/null > output/gromacs_md_autoinhib/C1458R_r1_$(date +%Y%m%d_%H%M).log 2>&1 &
setsid bash scripts/pipeline/run_7a6o_variant_direct.sh C1458R_r2  2 < /dev/null > ... &
setsid bash scripts/pipeline/run_7a6o_variant_direct.sh C1458R_r3  3 < /dev/null > ... &
setsid bash scripts/pipeline/run_7a6o_variant_direct.sh WT_r1      4 < /dev/null > ... &
setsid bash scripts/pipeline/run_7a6o_variant_direct.sh WT_r2      5 < /dev/null > ... &
setsid bash scripts/pipeline/run_7a6o_variant_direct.sh WT_r3      6 < /dev/null > ... &
```
监控:
```bash
bash scripts/pipeline/watch_7a6o_md_status.sh
pgrep -af 'run_7a6o_variant_direct|gmx.*mdrun'
```
**验证闸门 C:** 每条日志里 **不能有 `Fatal` / `LINCS` / `illegal memory access`**;
完成标志是 `md_prod.gro` 或 `md_prod.part*.gro`(noappend 续跑会写成 part)。

---

## 3. 分析

`analyze_7a6o_completed_md.py` 的 `DEFAULT_VARIANTS`/`LABELS` **没有 C1458R**,
两个办法二选一:
- (推荐)在脚本里把 C1458R 系列加进 `DEFAULT_VARIANTS`,并在 `LABELS` 里标
  `"C1458R": "patient"`(及 `_r1..r3`),避免 KeyError / 误标;**别把这个 de novo VUS
  标成已知 2B。**
- 或直接传 `--variants WT,WT_r1,WT_r2,WT_r3,C1458R,C1458R_r1,C1458R_r2,C1458R_r3`。

```bash
python3 scripts/pipeline/analyze_7a6o_completed_md.py \
    --input  output/gromacs_md_autoinhib \
    --output output/gromacs_md_autoinhib/analysis_completed_7a6o \
    --variants WT,WT_r1,WT_r2,WT_r3,C1458R,C1458R_r1,C1458R_r2,C1458R_r3 \
    --force
```
产出:`qc_summary.csv`、`aim_a1_contacts_summary.csv`、`aim_a1_contacts_timeseries.csv`。
读法:C1458R 的 AIM↔A1 非局部接触数相对 WT **下降** → 自抑制被破坏 → 偏 2B(GOF)方向;
**无明显变化** → 不支持自抑制释放,与 boltz GPIbα 强迫结合不升反降一致,更偏 2M/失稳。
3 条重复给出均值±SD,和已跑的 8 个 2B / 3 个 2M 参考分布比位置。

---

## 4. 推回来给我看(关键:数据默认被 gitignore 挡住)

`.gitignore:49` 的 `output/gromacs_md*/` 会挡掉所有 MD 输出。**只 commit 脚本没用,
我这边看不到数字。** 把分析 CSV 强制加进来:
```bash
git add -f output/gromacs_md_autoinhib/analysis_completed_7a6o/qc_summary.csv \
           output/gromacs_md_autoinhib/analysis_completed_7a6o/aim_a1_contacts_summary.csv \
           output/gromacs_md_autoinhib/analysis_completed_7a6o/aim_a1_contacts_timeseries.csv
git add scripts/  docs/
git commit -m "data(md): add C1458R patient + WT/C1458R replicates AIM-A1 MD analysis"
git push
```
(轨迹 .xtc 太大别推;CSV 足够我对齐 boltz 结论。)

---

## 5. 回来后我会做的

- 把 C1458R 的接触数 delta 放进 8×2B / 3×2M 参考分布,给方向判断。
- 与 boltz2 三轴(GPIbα 强迫结合不升反降 / heparan 偏高 / AIM 略高)叠加,
  形成一份 2B-vs-2M 的功能证据汇总。
- 明确写清:最终 2B/2M 仍需临床表型(LD-RIPA、血小板计数、多聚体电泳)定锤,
  MD/boltz 只是功能轴佐证。导出 `docs/patient_C1458R_report.md` 医生版。

## 验证闸门清单(No Human-Caused Errors)
- [ ] A: 7A6O 中 1458 已解析且 WT=CYS,offset 对齐
- [ ] B: 突变后 **C1272–C1458 二硫键已断**(topol.top 无残留 SS)
- [ ] B: 溶剂化 EM 收敛、refined EM Fmax 下降
- [ ] C: 每条 production 无 Fatal/LINCS/CUDA;完成有 md_prod(.part*).gro
- [ ] 分析: C1458R 标 `patient`,**不标成已知 2B**
- [ ] 推送: 用 `git add -f` 把 CSV 带出 gitignore

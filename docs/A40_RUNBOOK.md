# A40 Runbook — 独立服务器跑 VWF-ETHos 的 CPU + GPU 任务

适用场景:目标是一台 **独立的 A40 服务器**(不挂 `/lzy` NFS,CUDA 依赖与 H200 实例不同)。
因此 conda env 要**在 A40 上重建**(不能直接复用 LZY 的 OpenCL 构建)。

> TL;DR 顺序:① clone repo → ② 在 A40 建 CUDA 版 gromacs env → ③ `gpu_smoke_test.sh` 必须 PASS →
> ④ 跑 autoinhib MD(GPU)→ ⑤ 抽 AIM 自抑制特征 + 接分类器(CPU)。
> **smoke 不过不要上批量。**

---

## 0. 需要从 LZY 搬到 A40 的东西

| 内容 | 怎么搬 | 备注 |
|---|---|---|
| 代码 repo | `git clone <origin>` 或 `git pull` | `force_fields/charmm36m.ff/`(端基 patch)和 `opencl_vendors/nvidia.icd` 都在 git 里,跟着走 |
| autoinhib 输入结构 | scp `output/boltz2_a1_dp_d3_results/`(~1.3GB,5 变体) | MD 的输入 |
| (可选)自抑制 panel CIF | scp `output/boltz2_vwd_functional_panel/boltz_results/` 里 `*a1_aim_autoinhibition_context*`(较大) | 仅当要在 A40 上抽 2B 特征;否则在 LZY 抽好、只搬小 CSV(见 §5) |
| **不要搬** conda env | —— | LZY 的是 OpenCL+特定基础镜像构建,A40 上重建更干净(见 §2) |

---

## 1. 拿代码

```bash
git clone <你的 origin> VWF-ETHos && cd VWF-ETHos
# 或已有: git pull
export PROJ=$PWD
```

---

## 2. 在 A40 上重建 gromacs env(CUDA 优先)

回答你的问题:**OpenCL 能在 `conda create` 时加进去**——loader(`libOpenCL.so.1`)来自 conda 包
`ocl-icd`/`ocl-icd-libopencl1`,会进 env;但 NVIDIA vendor ICD(`libnvidia-opencl.so`)来自 A40 驱动,
得靠 `OCL_ICD_VENDORS` 兜底。**不过既然要重建,直接建 CUDA 构建更好**:GPU 常驻提速、只依赖 A40 驱动、
彻底绕开 OpenCL ICD,runner 的全套 flags 也名正言顺。

```bash
# 用 A40 自己的 conda(miniconda/mambaforge 均可)
conda create -p $PROJ/envs/gromacs python=3.11 -y

# 首选:CUDA 构建(自带 cudatoolkit,只需 A40 NVIDIA 驱动够新)
conda install -p $PROJ/envs/gromacs -c conda-forge -y \
    "gromacs=2025.*=nompi_cuda_*" gemmi numpy

# 若上面 build string 在你的 channel 解析不到,退回默认(通常是 OpenCL)构建:
#   conda install -p $PROJ/envs/gromacs -c conda-forge -y \
#       gromacs=2025 gemmi numpy ocl-icd ocl-icd-libopencl1
# 然后必须设 OpenCL ICD 兜底(见 §3b)。

export GMX=$PROJ/envs/gromacs/bin.AVX2_256/gmx     # 若无该目录, 用 .../bin/gmx
$GMX --version | grep -iE "GPU support|SIMD"        # ★确认: CUDA(理想) 或 OpenCL
```

> 分类器/特征抽取还需 numpy、pandas(在数据机上跑分类器时):
> `conda install -p $PROJ/envs/gromacs -c conda-forge -y pandas`

---

## 3. 环境自检

### 3a. 通用
```bash
export PROJ=$PWD
export GMX=$PROJ/envs/gromacs/bin.AVX2_256/gmx
export GMXLIB=$PROJ/force_fields            # patched charmm36m.ff (1MET 修复)
nvidia-smi -L                               # 期望: 1× NVIDIA A40
$GMX --version | grep -iE "GPU support|SIMD"
```

### 3b. 仅当 gromacs 是 OpenCL 构建
```bash
export OCL_ICD_VENDORS=$PROJ/opencl_vendors
ls /usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.*   # A40 驱动应有(版本无所谓, .icd 只指文件名)
cat $PROJ/opencl_vendors/nvidia.icd                  # 应是该 .so 的路径
```
> CUDA 构建可跳过 3b。

---

## 4. ★ GPU 放行闸(smoke)——不过不上批量

smoke 需要一个小 EM 体系(`output/_smoke_test/` 下的 `em.tpr` + `em.gro` + `topol.top`)。
若没随 repo 带,先在 A40 CPU 侧用任一 autoinhib 变体造一个:

```bash
# 用 WT 的 autoinhib CIF 造 smoke 体系(只到 EM)
$GMX --version >/dev/null    # 确认 GMX/GMXLIB 已 export
bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --filter WT --phase em --gpus 1
# 上面会在 output/gromacs_md_autoinhib/<WT>/ 产出 em.* 和 topology/topol.top
# 把它们软链到 _smoke_test (gpu_smoke_test.sh 默认在此找):
mkdir -p output/_smoke_test
WTDIR=$(ls -d output/gromacs_md_autoinhib/*WT*/ | head -1)
ln -sf "$PWD/$WTDIR/em/em.tpr"        output/_smoke_test/em.tpr
ln -sf "$PWD/$WTDIR/em/em.gro"        output/_smoke_test/em.gro
ln -sf "$PWD/$WTDIR/topology/topol.top" output/_smoke_test/topol.top
```

跑放行闸:

```bash
bash scripts/pipeline/gpu_smoke_test.sh
# Test A: -nb gpu -pme gpu (GPU 是否被发现)
# Test B: 用生产 runner 的真实 flags 实跑 (CUDA→含 -bonded/-update gpu; OpenCL→自动降级)
# 期望: "✅ SMOKE PASS"; 否则按提示修, 不要上批量
```

---

## 5. GPU 任务:autoinhib MD(D'D3–A1 自抑制模块, 5 变体)

```bash
# A40 单卡(或多卡 --gpus N);runner 会按 $GPU_BACKEND 自动选 flags
bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --gpus 1 --ns 200
# 断点续跑安全(.done 标记); 监控:
#   tail -f output/gromacs_md_autoinhib/worker_0.log
#   find output/gromacs_md_autoinhib -name .done_prod | wc -l
```
完成后后分析(D'D3–A1 接触/盐桥/接口 BSA 随时间):
```bash
$PROJ/envs/gromacs/bin/python scripts/pipeline/analyze_gromacs_md.py \
    --md-dir output/gromacs_md_autoinhib --system autoinhib
```
这 5 个 MD 是**金标准动态特征**,用来**校准**下面 §6 的廉价代理。

---

## 6. CPU 任务:2B 自抑制特征 → 分类器

> ⚠ 关键前提:`evidence_matrix.csv` 里现成的 `a1_aim_autoinhibition_context` 全局
> `ptm_or_plddt` **分不开 2B/2M**(实测 2B/2M delta 中位数 0.067 vs 0.049,全重叠)。
> 所以我们抽**接口级**特征(AIM↔A1 接触数),不是那个全局标量。

```bash
# (CIF 在哪台机就在哪台抽; 输出 CSV 很小, 可只搬 CSV)
$GMX_PY=$PROJ/envs/gromacs/bin/python    # 含 gemmi 的 python
$GMX_PY scripts/pipeline/extract_aim_autoinhib_features.py \
    --results-dir output/boltz2_vwd_functional_panel/boltz_results \
    --output output/aim_autoinhib_features.csv
# 输出列: variant_id, aa_change, aim_a1_contacts, ..., aim_release_score
# ★ sanity: 脚本会打印 WT 的 AIM-A1 接触残基对数; 若为 0 → 残基编号假设错了,
#   用 --construct-offset 调整 (CIF 残基非从 1 起时)。
```

接进分类器(把 `aim_release_score` merge 到分类器输入矩阵,按 `aa_change` 对齐):

```python
import pandas as pd
m = pd.read_parquet("VWF_Alpha_Matrix.parquet")

# 轴A: 自抑制松开 (来自 extract_aim_autoinhib_features.py)
a = pd.read_csv("output/aim_autoinhib_features.csv")
m = m.merge(a[["aa_change","aim_release_score"]], on="aa_change", how="left")

# 轴B: A1 结合面完整性 (forced_binding + heparan 两轴联合 LOF, 校准最优)
ev = pd.read_csv("output/boltz2_vwd_functional_panel/evidence_matrix.csv")
ev["aa_change"] = ev["wt_aa"] + ev["position"].astype(str) + ev["mut_aa"]
ev = ev.rename(columns={
    "a1_gpiba_forced_binding__primary_zscore_within_assay": "fb_binding_zscore",
    "a1_heparan_sulfate_binding__primary_zscore_within_assay": "heparan_zscore"})
m = m.merge(ev[["aa_change","fb_binding_zscore","heparan_zscore"]], on="aa_change", how="left")

m.to_parquet("VWF_Alpha_Matrix.parquet")
```
```bash
$GMX_PY scripts/agentic_vwf_classifier.py --validate -i VWF_Alpha_Matrix.parquet
```

`agentic_vwf_classifier.py` 的 RULE6 已加性接入(`aim_release_score` NaN 时行为不变):
- `aim_release_score ≥ AIM_RELEASE_2B_Z`(默认 1.0)→ 直接判 2B(自抑制松开);
- 介于 `AIM_RELEASE_LEAN_Z`(0.0)和阈值之间 → 把本会默认判 2M 的 A1 变体救成 2B。

### ★ 校准(必须做)
默认阈值是**暂定**的。用已知标签调:
```bash
# 看已知 2B vs 2M 在 aim_release_score 上的分布, 选能最好分开的阈值,
# 改 scripts/agentic_vwf_classifier.py 顶部的 AIM_RELEASE_2B_Z / AIM_RELEASE_LEAN_Z,
# 重跑 --validate 看 2B recall(之前 12 个已知 2B 只对 2 个)。
```

---

## 7. 常见问题

| 症状 | 原因 / 处理 |
|---|---|
| `mdrun ... not supported with OpenCL` | gromacs 是 OpenCL 构建但用了 `-bonded/-update gpu`。runner 已自动降级;若仍出现,检查 §3a 的 `GPU support:` 是否被正确识别为 OpenCL |
| `0 detected device(s)` | OpenCL ICD 没生效。做 §3b;CUDA 构建则查 `nvidia-smi` 与驱动 |
| `atom C1 not found in building block 1MET` | `GMXLIB` 没指向 `force_fields/`(patched FF)。`export GMXLIB=$PROJ/force_fields` |
| extract 脚本 WT 接触数=0 | CIF 残基编号非从 1 起 → 调 `--construct-offset`(construct_local = VWF_pos - offset) |
| `No module named gemmi/pandas` | 在对应 env 里 `conda install -c conda-forge gemmi pandas` |

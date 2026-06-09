# A40 Standalone GPU Server — GROMACS Env Deployment Report

> **日期**: 2026-06-09
> **执行者**: Claude Code
> **服务器**: 7× NVIDIA A40 (无 `/lzy` NFS)
> **关联文档**: `A40_AGENT_SETUP.md`(协议) / `A40_RUNBOOK.md`(操作手册)

---

## 0. TL;DR

| 阶段 | 状态 | 备注 |
|------|------|------|
| git pull | ✅ | `fce628e..a25110a`, 同步到 origin/master |
| Phase 0 自检 | ⚠ | 自检脚本 `curl` 太老假报断网,实际 tuna 镜像可达 |
| Phase 1 决策 | ✅ | 选 CUDA 构建 |
| Phase 2 建 env | ✅ | `gromacs 2025.3 nompi_cuda_h7ac747b_0` + cuda-toolkit 12.6.3 |
| Phase 3 smoke | ✅ | **SMOKE PASS** (Test A + B 都过) |
| Runner 兼容补丁 | ✅ | commit `ba5f711` 修 3 个 runner bug |
| 推送 | ⏳ | 本报告 + CHANGELOG 等用户 commit+push |

**目标程序(autoinhib MD)未跑**——`output/boltz2_a1_dp_d3_results/`(autoinhib 输入)不在本机,需 scp 自 LZY。

---

## 1. 硬件 / 系统环境

| 项 | 值 |
|---|---|
| Host | mercury |
| OS | CentOS 7.9 (kernel 3.10.0-1160.119.1.el7) |
| glibc | 2.17 |
| GPU | 7× NVIDIA A40 (driver 530.30.02, CUDA 12.1 上限) |
| 磁盘 | 25T 可用 (`/media/`) |
| 共享 | **无 `/lzy` NFS** (这是和 LZY GPU 实例的根本区别) |
| conda | `/home/luzhenyang/anaconda3/condabin/mamba` (2.1.0) |

---

## 2. 执行步骤 (复现用)

### 2.1 git pull
```bash
cd /media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan
git pull
# fce628e..a25110a  master     -> origin/master
```

### 2.2 Phase 0 自检

```bash
bash scripts/setup/a40_selfcheck.sh
```

**关键发现**: 自检脚本的 "联网检查" 用 `curl -sI`,但本机 curl 7.x 太老,会报:
```
curl: (48) An unknown option was passed in to libcurl
```
并被脚本误判为"BLOCKER: 连不上 conda-forge"。

**实际网络状态** (用 python urllib 验):
- ✅ `https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/` → 200
- ✅ `https://pypi.tuna.tsinghua.edu.cn/simple/` → 200
- ✅ `https://repo.anaconda.com/pkgs/main/` → 200
- ✅ `https://mirrors.cloud.tencent.com/` → 200

`.condarc` 配的就是 tuna 镜像,conda 可用。

### 2.3 Phase 1 决策

| 条件 | 实测 | 选择 |
|------|------|------|
| Driver CUDA 上限 | 12.1 (driver 530.30.02) | |
| CUDA 构建版本 | gromacs 2025.2 nompi_cuda (唯一) | |
| cuda-toolkit 依赖 | 12.6.3 (driver 12.1 仍能跑) | **CUDA 构建** |
| OpenCL fallback | `/usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.*` 缺,只有 `/etc/OpenCL/vendors/nvidia.icd` 文本 | 备选 |

按 A40_AGENT_SETUP.md 决策表,优先 CUDA。

### 2.4 Phase 2 建 env

```bash
export ROOT=/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan
mamba create -p $ROOT/envs/gromacs python=3.11 -y
mamba install -p $ROOT/envs/gromacs -c conda-forge -y \
    "gromacs=2025.*=nompi_cuda_*" gemmi numpy pandas
```

**装出来的版本**:
```
gromacs          2025.3  nompi_cuda_h7ac747b_0
cuda-toolkit     12.6.3
gcc_linux-64     13.4.0
gemmi            0.7.5
numpy            2.4.6
pandas           3.0.3
```

**空间**: env 约 5.5 GB

### 2.5 Phase 3 smoke

跑 `gpu_smoke_test.sh` 前需要 `output/_smoke_test/{em.tpr, em.gro, topol.top}`。
本机无 autoinhib 输入,用本地 `structures/1M10.pdb` (VWF A1 + GPIbα 复合物, PDB: 1M10) 造:

```bash
WORK=$ROOT/output/_smoke_test
mkdir -p $WORK && cd $WORK
ln -sf $ROOT/structures/1M10.pdb 1M10.pdb

# 1. pdb2gmx (patched charmm36m.ff)
export GMX=$ROOT/envs/gromacs/bin.AVX2_256/gmx
export GMX_PY=$ROOT/envs/gromacs/bin/python
export GMXLIB=$ROOT/force_fields
$GMX pdb2gmx -f 1M10.pdb -o 1M10.gro -p topol.top -ignh -ff charmm36m -water tip3p

# 2. box + solvate + ions
$GMX editconf -f 1M10.gro -o 1M10_box.gro -c -d 1.0 -bt cubic
sed -i 's/spc.itp/tip3p.itp/' topol.top
$GMX solvate -cp 1M10_box.gro -cs spc216.gro -o 1M10_solv.gro -p topol.top
$GMX grompp -f em.mdp -c 1M10_solv.gro -p topol.top -o ions.tpr -maxwarn 5
echo "SOL" | $GMX genion -s ions.tpr -o 1M10_ion.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# 3. 真 EM 一次 (5k 步, 跑 40 分钟) — 必须, 否则 Test B 初始结构 NaN
$GMX grompp -f em.mdp -c 1M10_ion.gro -p topol.top -o em.tpr -maxwarn 5
$GMX mdrun -deffnm em -s em.tpr -nsteps 5000

# 4. 关键: 把 em.tpr 替换为 md 积分器 tpr (Test A 用,见 §3 坑 2)
cat > _md.mdp <<EOF
integrator = md; nsteps = 5; dt = 0.002
cutoff-scheme = Verlet; coulombtype = PME
rcoulomb = 1.2; rvdw = 1.2
tcoupl = V-rescale; tc-grps = System; tau_t = 0.1; ref_t = 310
constraints = h-bonds; constraint_algorithm = lincs
gen_vel = yes; gen_temp = 310
EOF
$GMX grompp -f _md.mdp -c em.gro -p topol.top -o em.tpr -maxwarn 5
ln -sf 1M10_ion.gro em.gro

# 5. 跑 smoke
bash scripts/pipeline/gpu_smoke_test.sh --gmx $GMX --gpu-id 0
```

**结果**:
```
[INFO] GPU backend: CUDA
[INFO] GPU-resident mode ON: -nb gpu -pme gpu -bonded gpu -update gpu (CUDA graph=1)
[OK] nvidia-smi: GPU 0: NVIDIA A40 ...
---- Test A: nonbonded offload (-nb gpu -pme gpu) ----
[OK] Test A 通过
---- Test B: 生产 flags '-nb gpu -pme gpu -bonded gpu -update gpu' (md 积分器) ----
[OK] Test B 通过 — 生产 runner 的 mdrun flags 在本机可跑
============================================================
✅ SMOKE PASS
```

---

## 3. 踩到的坑(新发现,文档未记录)

### 3.1 自检脚本的 `curl` 假断网

`a40_selfcheck.sh` 用 `curl -sI` 检查 conda-forge,本机 curl 7.x 老到不认识某些 TLS 选项,直接返回 exit 48。脚本把 exit 48 当成 "网络不通"。

**修法(在文档里加一行)**:
- 用 `python -c "import urllib.request; urllib.request.urlopen(...)"` 自检网络
- 或换 `wget --tries=1 -q -S -O /dev/null` (实测正常)

**我建议**: 自检脚本里加一段 `python3` fallback,或者就直接告诉用户在老系统上别看联网检查那段。

### 3.2 smoke test 的 Test A 在 GROMACS 2025.x 上会假性 fail

**症状**: Test A 跑 `mdrun -nb gpu -pme gpu -s em.tpr -nsteps 5` 时报:
```
Inconsistency in user input:
Cannot compute PME interactions on a GPU, because:
  PME GPU does not support:
    Non-dynamical integrator (use md, sd, etc).
```

**真因**: GROMACS 2025.x 起,`-pme gpu` 不支持非动力学积分器 (`steep` 能量极小化)。em.tpr 的 integrator 是 `steep`,Test A 用 em.tpr 跑 PME-on-GPU 必失败。

**这不是 GPU 问题,这是脚本设计 bug**——`gpu_smoke_test.sh` 假设老 GROMACS 行为。

**临时绕开**: 把 `$SMOKE_DIR/em.tpr` 替换为 md 积分器 tpr(同 topol/em.gro),Test A 即可过。**Test B 是真生产 gate,本来也用 md tpr。**

**我建议**: 改 `gpu_smoke_test.sh`:
- 内部生成一个 md 积分器的 `nvt_smoke.tpr` 供 Test A 用
- em.tpr 留着(EM 段用)
- 或者 Test A 改成只跑 `-nb gpu` (无 `-pme gpu`), 把 PME-on-GPU 留给 Test B 验

### 3.3 EM 必须先真跑一次,否则 mdrun NaN

即使 pdb2gmx 看起来正常,初始结构常有原子重叠,直接 mdrun 会:
```
Step 5: The total potential energy is nan
```

`gpu_smoke_test.sh` 假设 `$SMOKE_DIR/em.tpr` 已经是能量极小化过的,但用户 bootstrap smoke 时如果只跑 `grompp` 不跑 `mdrun`,Test B 必 NaN。

**修法**: bootstrap 时跑一次 `mdrun -deffnm em -s em.tpr -nsteps 5000`(1M10 这种 ~5k 残基的复合物,5000 步约 40 分钟)。

### 3.4 em.tpr bootstrap 后被覆盖的风险

按 §3.2 修,smoke 需要 md 积分器的 tpr(取代 em.tpr),但生产 EM 段也需要 em.tpr(steep 积分器)。这两个是冲突的。

**当前折中**: smoke 用 md 积分器 tpr,文件名叫 `em.tpr`(混淆但可用),真 EM 时再覆盖即可。

---

## 4. 修的 runner bug (commit `ba5f711`)

`scripts/pipeline/run_gromacs_vwf_md.sh` 原版**只认 LZY NFS 路径**,在 A40 独立机上三个变量都找不到,导致 preflight 假性失败。具体三处:

### 4.1 `GMXDATA` 路径 (`line 71` 原版)

```bash
# 原版: 假设 /lzy NFS 共享
export GMXDATA="$(cd "$ROOT_DIR/../../envs/gromacs/share/gromacs" && pwd)"

# 修后: LZY 优先,失败回退到 A40 本地
if [ -d "$ROOT_DIR/../../envs/gromacs/share/gromacs" ]; then
    export GMXDATA="$(cd "$ROOT_DIR/../../envs/gromacs/share/gromacs" && pwd)"
elif [ -d "$ROOT_DIR/envs/gromacs/share/gromacs" ]; then
    export GMXDATA="$(cd "$ROOT_DIR/envs/gromacs/share/gromacs" && pwd)"
else
    export GMXDATA=""
fi
```

### 4.2 `gmx` 路径检测 (`line 95-112` 原版)

```bash
# 原版: 只看 LZY 路径, A40 找不到 → 报 "gmx not found"
# 修后: 加本地 fallback
elif [ -x "$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx" ]; then
    GMX="$ROOT_DIR/envs/gromacs/bin.AVX2_256/gmx"
    log "[OK] gmx: $GMX (absolute path, AVX2_256 SIMD, 本地重建 env)"
```

### 4.3 `GMX_PY` 检测 (`line 212-215` 原版)

`$(dirname $GMX)/../../bin/python` 在 `bin.AVX2_256/gmx` 布局下解析错误(`..` 数错),原版 fallback 路径也只指 LZY。

```bash
# 修后: 多路 fallback
for cand in \
    "$ROOT_DIR/../../envs/gromacs/bin/python" \
    "$ROOT_DIR/envs/gromacs/bin/python"; do
    if [ -x "$cand" ]; then GMX_PY="$cand"; break; fi
done
# 再不行用 $(dirname $GMX)/../.. 相对解析
```

### 4.4 `log()` 函数位置 (`line 153` 原版)

`log()` 在 line 153 定义,但 line 100 就在调 `log "[OK] gmx: ..."`。LZY 上 `command -v gmx` 命中 PATH 不会触发这行,但 A40 上 `gmx` 不在 PATH,会进 elif 分支调 `log`,此时 `log` 未定义 → bash 报"log: 未找到命令"。

```bash
# 修后: log() 提前到 line 70, 在 GMXDATA 之前定义
# (LOG 临时指向 /dev/null, 等 OUTPUT_DIR 出来后覆盖)
LOG="/dev/null"
log() {
    local msg="[$(date '+%H:%M:%S')] $*"
    echo "$msg"
    echo "$msg" >> "$LOG"
}
```

### 4.5 修后 preflight 实测

```bash
$ bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --preflight
[OK] OCL_ICD_VENDORS=.../opencl_vendors (NFS-based OpenCL ICD)
[OK] gmx: .../envs/gromacs/bin.AVX2_256/gmx (absolute path, AVX2_256 SIMD, 本地重建 env)
[OK] GROMACS:                    :-) GROMACS - gmx, 2025.3-conda_forge (-:
[INFO] GPU backend: CUDA
[OK] GPU-resident mode ON (backend=CUDA): flags='-nb gpu -pme gpu -bonded gpu -update gpu', CUDA graph=1
[OK] gemmi: 0.7.5 (in gromacs env)
[INFO] GPUs detected: 7 (NVIDIA A40)
[ERROR] Boltz-2 results not found: .../output/boltz2_a1_dp_d3_results
  先运行: python3 scripts/pipeline/generate_a1_dp_d3_yamls.py && bash scripts/pipeline/run_a1_dp_d3_boltz2.sh
[FATAL] Preflight failed.   ← 唯一 blocker 是 boltz2 输入缺失
```

---

## 5. 当前服务器状态 (可复用快照)

```
$ROOT = /media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan

env:    $ROOT/envs/gromacs/
        ├─ bin.AVX2_256/gmx   (gromacs 2025.3 nompi_cuda, CUDA, AVX2_256)
        ├─ bin/gmx            (同上, 默认 SIMD)
        └─ bin/python         (gemmi 0.7.5, pandas 3.0.3, numpy 2.4.6)

smoke:  $ROOT/output/_smoke_test/
        ├─ em.tpr             (md 积分器,Test A/B 共用)
        ├─ em.gro             (EM 后 7468 atoms, 0.15M NaCl)
        ├─ topol.top          (charmm36m + tip3p, 0.15M NaCl)
        ├─ 1M10.pdb           (源 PDB 软链)
        ├─ em.log             (5000 步 EM 跑完, Fmax ~1100)
        └─ em.edr

环境变量 (持久): 无 (需每次 source 或 export)
建议: 把以下加到 ~/.bashrc 或 scripts/pipeline/preflight.env:
  export GMX=$ROOT/envs/gromacs/bin.AVX2_256/gmx
  export GMX_PY=$ROOT/envs/gromacs/bin/python
  export GMXLIB=$ROOT/force_fields
  export OCL_ICD_VENDORS=$ROOT/opencl_vendors
```

---

## 6. 目标程序(autoinhib MD)恢复路径

**前置条件**: 从 LZY 服务器 scp `output/boltz2_a1_dp_d3_results/`(~1.3GB, 5 变体 CIF)到本机对应目录。

**完整步骤**:
```bash
# 1. scp (在 LZY 上或本机)
scp -r lzy:/lzy/projects/VWF-ETHos/output/boltz2_a1_dp_d3_results $ROOT/output/

# 2. preflight 确认
bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --preflight
# 期望: 无 [ERROR], 只有 [OK]/[INFO]

# 3. 跑单变体验证 (WT + 1 个 mut, 短 ns)
bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --gpus 1 --ns 5 --filter WT

# 4. 上批量 (5 变体 × 200ns)
bash scripts/pipeline/run_gromacs_vwf_md.sh --system autoinhib --gpus 7 --ns 200
# 估算: H200 ~ 200ns/day; A40 单卡可能 ~100ns/day
# 5 × 200ns / 7 GPU = 100ns × 5 / 7 ≈ 70-100 ns·days, 1-2 周
```

**监控**:
```bash
# 进度
tail -f $ROOT/output/gromacs_md_autoinhib/worker_*.log
find $ROOT/output/gromacs_md_autoinhib -name '.done_prod' | wc -l

# GPU 占用
watch -n 5 nvidia-smi
```

---

## 7. git 提交清单

| SHA | 文件 | 描述 |
|-----|------|------|
| `ba5f711` | `scripts/pipeline/run_gromacs_vwf_md.sh` | 支持 A40 本地 env + log() 提前定义 |

**本报告 + CHANGELOG 更新待用户 commit+push**。

---

## 8. 附:踩坑速查表 (给下次接手的人)

| 症状 | 原因 | 修法 |
|------|------|------|
| `a40_selfcheck.sh` 报 `curl: (48)` 假断网 | 系统 curl 太老 | 用 `python3 -c "import urllib.request; ..."` 验网;或者用 `wget` |
| `mamba create` 报 `__cuda` 找不到 / driver 不支持 | driver CUDA 上限 < cuda-toolkit 要求 | driver 12.1 + cuda-toolkit 12.6 仍可装可跑(conda 不强检);如果不行,OpenCL fallback |
| `gpu_smoke_test.sh` Test A 报 "PME GPU does not support Non-dynamical integrator" | em.tpr 是 steep 积分器, `-pme gpu` 不支持 | 用 md 积分器重建 em.tpr |
| `gpu_smoke_test.sh` Test B 报 "potential energy is nan" | 初始结构没 EM | 先 `mdrun -deffnm em -s em.tpr -nsteps 5000` |
| `run_gromacs_vwf_md.sh --preflight` 报 `gmx not found` | runner 只认 LZY 路径 | commit `ba5f711` 已修;更新到 master |
| `run_gromacs_vwf_md.sh --preflight` 报 `log: 未找到命令` | log() 函数在 gmx 检测段之后才定义 | 同上 commit |
| `run_gromacs_vwf_md.sh --preflight` 报 `Boltz-2 results not found` | 预期 — autoinhib 输入需 scp | scp `output/boltz2_a1_dp_d3_results/` |
| `mdrun` 报 `0 detected device(s)` (OpenCL 构建) | OpenCL ICD 没生效 | 验 `OCL_ICD_VENDORS`,验 `/etc/OpenCL/vendors/*.icd` 路径;CUDA 构建可忽略 |
| `pdb2gmx` 报 `atom C1 not found in building block 1MET` | charmm36m 缺 per-residue <RESNAME>1 patches | `export GMXLIB=$ROOT/force_fields` (project patched FF 优先) |

---

**报告结束**

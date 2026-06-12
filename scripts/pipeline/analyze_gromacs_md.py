#!/usr/bin/env python3
"""
analyze_gromacs_md.py
=====================
Post-analysis for GROMACS MD trajectories across complex/monomer/autoinhib systems.

Supports three system types:
  complex    — A1-GPIbα complex. Metrics: RMSF, PCA, MM/PBSA ΔG_bind.
  monomer   — VWF A1 only (chain A extracted from complex CIF).
              Metrics: RMSF, PCA, mutation-site flexibility.
  autoinhib — VWF AIM-A1 (7A6O 实验结构路线). Metrics: RMSF + AIM↔A1 自抑制接口
              接触数随轨迹 (2B = 松开 = 接触↓) + WT vs 2B vs 2M 判定。
              (旧 D'D3-A1 Boltz 体系已弃用; 默认 AIM/A1 范围为 7A6O VWF 编号。)

Usage:
  python3 analyze_gromacs_md.py --system complex --input output/gromacs_md/
  python3 analyze_gromacs_md.py --system monomer --input output/gromacs_md_monomer/
  python3 analyze_gromacs_md.py --system autoinhib --input output/gromacs_md_autoinhib/

Output:
  analysis/summary.csv      — per-variant metrics
  analysis/rmsf.csv         — per-residue RMSF
  analysis/pca.csv           — PCA eigenvalues / top PCs
  analysis/contacts_timeseries.csv — AIM-A1 接触数随帧 (autoinhib)
  analysis/verdict_2b_vs_2m.csv     — 每变体平衡态接触 + Δvs_WT + 标签 (autoinhib)
  analysis/mmpbsa.csv        — MM/PBSA binding free energy (complex only)
"""

import argparse
import os
import subprocess
from pathlib import Path
from typing import Optional, List

import pandas as pd


# =============================================================================
# Helpers
# =============================================================================

def run_gmx(gmx_args: List[str], work_dir: str, log_file: Optional[str] = None,
            stdin: Optional[str] = None, gmx_bin: str = "gmx") -> str:
    """Run a gmx command and return stdout+stderr. `stdin` feeds group selections (trjconv 等)."""
    env = os.environ.copy()
    env["GMX_RESET_COUNTERS"] = "0"
    result = subprocess.run(
        [gmx_bin] + gmx_args,
        cwd=work_dir,
        capture_output=True,
        text=True,
        env=env,
        input=stdin,
    )
    output = result.stdout + result.stderr
    if log_file:
        Path(log_file).write_text(output)
    return output


def compute_rmsf(tpr: str, xtc: str, out_dir: str) -> pd.DataFrame:
    """Compute per-residue RMSF."""
    rmsf_xvg = os.path.join(out_dir, "rmsf.xvg")
    run_gmx(
        ["rmsf", "-s", tpr, "-f", xtc, "-o", rmsf_xvg,
         "-qdist", "1", "-pdb", "none"],
        work_dir=out_dir
    )
    # Parse XVG file
    data = []
    with open(rmsf_xvg) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    res, val = float(parts[0]), float(parts[1])
                    data.append({"residue": int(res), "rmsf": val})
                except ValueError:
                    continue
    return pd.DataFrame(data)


def compute_pca(tpr: str, xtc: str, out_dir: str, n_components: int = 10) -> pd.DataFrame:
    """Compute PCA eigenvalues from covariance analysis."""
    covar_out = os.path.join(out_dir, "covar")
    eigenval_xvg = os.path.join(out_dir, "eigenval.xvg")

    # Covariance analysis
    run_gmx(
        ["covar", "-s", tpr, "-f", xtc, "-o", covar_out,
         "-eigen", eigenval_xvg, "-l", os.path.join(out_dir, "covar.log")],
        work_dir=out_dir,
        log_file=os.path.join(out_dir, "covar_output.log")
    )

    # Parse eigenvalues
    evals = []
    with open(eigenval_xvg) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    pc, val = int(parts[0]), float(parts[1])
                    evals.append({"PC": pc, "eigenvalue": val})
                except ValueError:
                    continue

    df = pd.DataFrame(evals)
    if len(df) > 0:
        total = df["eigenvalue"].sum()
        df["variance_explained"] = df["eigenvalue"] / total if total > 0 else 0
        df["cumulative_variance"] = df["variance_explained"].cumsum()
    return df.head(n_components)


# ---- AIM-A1 自抑制接口 (7A6O 体系, VWF 规范编号) -----------------------------
# 7A6O 解析范围 ~1262-1466 (单链)。自抑制 = AIM 片段折回压住 A1 的 GPIb 结合面。
# 2B = 突变破坏自抑制 → AIM 从 A1 脱离 → AIM↔A1 接触数随轨迹下降 (闭合态不稳)。
AIM_RANGES_VWF = [(1262, 1271), (1460, 1466)]   # 7A6O 中解析的 N-AIM 残端 + C-AIM
A1_CORE_VWF = [(1272, 1459)]                      # 两段 AIM 之间的 A1 球状体
POS_RES = {"ARG", "LYS", "HIS"}
NEG_RES = {"ASP", "GLU"}

# 7A6O 批次 14 突变体的已知标签分组 (来自 2026-06-11 handoff; 可 --labels CSV 覆盖)
DEFAULT_LABELS = {
    "R1306W": "2B", "R1306Q": "2B", "R1308C": "2B", "I1309V": "2B", "S1310F": "2B",
    "W1313C": "2B", "V1314F": "2B", "V1316M": "2B",
    "R1374C": "2M", "R1374H": "2M", "G1324S": "2M",
    "P1337L": "?", "R1341Q": "?", "R1341W": "?",
}


def _parse_ranges(s: str):
    out = []
    for p in s.split(","):
        a, b = p.split("-")
        out.append((int(a), int(b)))
    return out


def _in_ranges(n: int, ranges) -> bool:
    return any(a <= n <= b for a, b in ranges)


def _shift_ranges(ranges, shift):
    return [(a + shift, b + shift) for a, b in ranges]


def extract_frames(tpr: str, xtc: str, gro: str, out_dir: str, skip: int, gmx_bin: str) -> Optional[str]:
    """轨迹 → 蛋白-only 多模型 PDB(-pbc whole, 每 skip 帧)。无 xtc 则用末帧 gro。"""
    frames = os.path.join(out_dir, "frames_protein.pdb")
    if xtc and os.path.exists(xtc) and os.path.exists(tpr):
        run_gmx(["trjconv", "-s", tpr, "-f", xtc, "-o", frames,
                 "-pbc", "whole", "-skip", str(skip)],
                out_dir, log_file=os.path.join(out_dir, "trjconv.log"),
                stdin="Protein\n", gmx_bin=gmx_bin)
    elif gro and os.path.exists(gro):
        # 末帧: editconf 读 .gro 写 pdb(整盒, 含水; gemmi 仍可按残基号筛蛋白)
        run_gmx(["editconf", "-f", gro, "-o", frames], out_dir,
                log_file=os.path.join(out_dir, "editconf.log"), gmx_bin=gmx_bin)
    return frames if os.path.exists(frames) and os.path.getsize(frames) > 0 else None


def aim_a1_contacts_timeseries(frames_pdb: str, aim_ranges, a1_ranges,
                               cutoff: float) -> pd.DataFrame:
    """每帧 AIM↔A1-core 接触残基对数 + 盐桥数。自动适配局部/规范残基编号。"""
    try:
        import gemmi
    except ImportError:
        print("[WARN] gemmi not installed, AIM-A1 contact analysis skipped")
        return pd.DataFrame()

    st = gemmi.read_structure(frames_pdb)
    if len(st) == 0:
        return pd.DataFrame()

    # 编号适配: 若结构是局部编号(max<1000), 把 VWF 范围平移到结构编号
    allnums = [r.seqid.num for r in st[0][0]] if len(st[0]) else []
    if allnums and max(allnums) < 1000:
        lo = min(allnums)
        shift = lo - AIM_RANGES_VWF[0][0]      # VWF 起点映到结构起点
        aim_ranges = _shift_ranges(aim_ranges, shift)
        a1_ranges = _shift_ranges(a1_ranges, shift)
        print(f"    [info] 结构用局部编号(min={lo}) → AIM/A1 平移 {shift}: AIM={aim_ranges} A1={a1_ranges}")

    rows = []
    for fi, model in enumerate(st):
        ns = gemmi.NeighborSearch(model, st.cell, cutoff).populate()
        # A1 残基号集合(用于判断邻居属于 A1)
        a1_nums, aim_residues = set(), []
        for chain in model:
            for res in chain:
                if _in_ranges(res.seqid.num, a1_ranges):
                    a1_nums.add(res.seqid.num)
                elif _in_ranges(res.seqid.num, aim_ranges):
                    aim_residues.append(res)
        pairs, salt = set(), set()
        for res in aim_residues:
            res_is_pos = res.name in POS_RES
            res_is_neg = res.name in NEG_RES
            for atom in res:
                if atom.element.name in ("H", "D"):
                    continue
                for m in ns.find_atoms(atom.pos, "\0", radius=cutoff):
                    cra = m.to_cra(model)
                    rn = cra.residue.seqid.num
                    if rn not in a1_nums or cra.atom.element.name in ("H", "D"):
                        continue
                    if atom.pos.dist(cra.atom.pos) <= cutoff:
                        pairs.add((res.seqid.num, rn))
                        if (res_is_pos and cra.residue.name in NEG_RES) or \
                           (res_is_neg and cra.residue.name in POS_RES):
                            salt.add((res.seqid.num, rn))
        rows.append({"frame": fi, "n_contacts": len(pairs), "n_salt_bridges": len(salt)})
    return pd.DataFrame(rows)


def run_mmpbsa(tpr: str, xtc: str, out_dir: str) -> pd.DataFrame:
    """Compute MM/PBSA binding free energy (requires gmx_MMPBSA)."""
    try:
        import gmx_MMPBSA
    except ImportError:
        print("[WARN] gmx_MMPBSA not installed, MM/PBSA skipped")
        return pd.DataFrame()

    # This is a simplified wrapper; production use would call gmx_MMPBSA.AMMPSAMCalculation
    print("[INFO] MM/PBSA requires running: gmx_MMPBSA -O -i mmpbsa.in -fo mmgbsa.csv")
    return pd.DataFrame()


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--system", choices=["complex", "monomer", "autoinhib"],
                        required=True, help="MD system type")
    parser.add_argument("--input", required=True,
                        help="GROMACS MD output directory (e.g. output/gromacs_md/)")
    parser.add_argument("--output", default=None,
                        help="Analysis output directory (default: <input>/analysis/)")
    parser.add_argument("--variants", default=None,
                        help="Variant list CSV for mutation-site annotations")
    parser.add_argument("--aim-ranges", default="1262-1271,1460-1466",
                        help="AIM 残基区间 (VWF 编号, 逗号分隔; 7A6O 默认)")
    parser.add_argument("--a1-core", default="1272-1459", help="A1 球状核区间 (VWF 编号)")
    parser.add_argument("--contact-cutoff", type=float, default=4.5, help="接触阈值 Å")
    parser.add_argument("--skip", type=int, default=10, help="trjconv 每 skip 帧取一帧")
    parser.add_argument("--gmx", default="gmx", help="gmx 可执行 (默认 PATH 中的 gmx)")
    parser.add_argument("--labels", default=None,
                        help="CSV(列 variant,label[2B/2M])覆盖默认已知标签分组")
    args = parser.parse_args()

    input_dir = Path(args.input).resolve()
    out_dir = Path(args.output or input_dir / "analysis").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"System: {args.system}")
    print(f"Input:  {input_dir}")
    print(f"Output: {out_dir}")
    print()

    # Find completed variants. 支持两种布局:
    #   7A6O 批量 runner: <variant>/md_7a6o/md_prod.{tpr,xtc}
    #   旧 runner:        <variant>/production/md_prod.{tpr,xtc}
    # 以 md_prod.xtc 存在判定完成(批量 runner 不写 .done_prod 标记)。
    runs = []  # (variant, prod_dir)
    for d in sorted(input_dir.iterdir()):
        if not d.is_dir():
            continue
        for sub in ("md_7a6o", "production"):
            if (d / sub / "md_prod.xtc").exists():
                runs.append((d.name, d / sub))
                break

    print(f"Found {len(runs)} completed variants")
    if not runs:
        print("No completed simulations found (looked for */md_7a6o/md_prod.xtc).")
        return

    aim_ranges = _parse_ranges(args.aim_ranges)
    a1_ranges = _parse_ranges(args.a1_core)
    labels = dict(DEFAULT_LABELS)
    if args.labels and Path(args.labels).exists():
        ldf = pd.read_csv(args.labels)
        labels.update({str(r["variant"]): str(r["label"]) for _, r in ldf.iterrows()})

    all_rows = []
    all_rmsf = []
    all_pca = []
    all_contacts = []

    for variant, prod_dir in runs:
        tpr = prod_dir / "md_prod.tpr"
        xtc = prod_dir / "md_prod.xtc"
        print(f"  Analyzing: {variant}")

        row = {"variant": variant}

        # RMSF
        try:
            rmsf_df = compute_rmsf(str(tpr), str(xtc), str(prod_dir))
            if not rmsf_df.empty:
                rmsf_df["variant"] = variant
                all_rmsf.append(rmsf_df)
                row["rmsf_mean"] = rmsf_df["rmsf"].mean()
                row["rmsf_max"] = rmsf_df["rmsf"].max()
        except Exception as e:
            print(f"    [WARN] RMSF failed: {e}")

        # PCA
        try:
            pca_df = compute_pca(str(tpr), str(xtc), str(prod_dir))
            if not pca_df.empty:
                pca_df["variant"] = variant
                all_pca.append(pca_df)
                row["PC1_variance"] = pca_df.iloc[0]["variance_explained"] if len(pca_df) > 0 else None
        except Exception as e:
            print(f"    [WARN] PCA failed: {e}")

        # System-specific analysis: AIM↔A1 自抑制接触随轨迹 (2B = 松开 = 接触↓)
        if args.system == "autoinhib":
            try:
                frames = extract_frames(str(tpr), str(xtc), str(prod_dir / "md_prod.gro"),
                                        str(prod_dir), args.skip, args.gmx)
                if frames:
                    cdf = aim_a1_contacts_timeseries(frames, aim_ranges, a1_ranges,
                                                     args.contact_cutoff)
                    if not cdf.empty:
                        cdf["variant"] = variant
                        all_contacts.append(cdf)
                        last = cdf.tail(max(1, len(cdf) // 4))   # 后 25% 当平衡态
                        row["aim_a1_contacts_mean"] = round(cdf["n_contacts"].mean(), 2)
                        row["aim_a1_contacts_eq"] = round(last["n_contacts"].mean(), 2)
                        row["aim_a1_saltbridges_eq"] = round(last["n_salt_bridges"].mean(), 2)
                        row["n_frames"] = len(cdf)
                else:
                    print("    [WARN] 取帧失败 (无 xtc/gro?)")
            except Exception as e:
                print(f"    [WARN] AIM-A1 contacts failed: {e}")

        all_rows.append(row)

    # Write outputs
    summary_df = pd.DataFrame(all_rows) if all_rows else pd.DataFrame()
    if all_contacts:
        pd.concat(all_contacts, ignore_index=True).to_csv(out_dir / "contacts_timeseries.csv", index=False)
        print(f"\nWritten: {out_dir / 'contacts_timeseries.csv'}")

    if not summary_df.empty:
        summary_df.to_csv(out_dir / "summary.csv", index=False)
        print(f"Written: {out_dir / 'summary.csv'} ({len(summary_df)} variants)")

    if all_rmsf:
        pd.concat(all_rmsf, ignore_index=True).to_csv(out_dir / "rmsf.csv", index=False)
        print(f"Written: {out_dir / 'rmsf.csv'}")
    if all_pca:
        pd.concat(all_pca, ignore_index=True).to_csv(out_dir / "pca.csv", index=False)
        print(f"Written: {out_dir / 'pca.csv'}")

    # ---- 2B vs 2M 判定 (autoinhib): WT 基线, 每变体接触数 Δ, 按已知标签分组 ----
    if args.system == "autoinhib" and not summary_df.empty and "aim_a1_contacts_eq" in summary_df:
        summary_df["label"] = summary_df["variant"].map(
            lambda v: "WT" if "WT" in str(v).upper() else labels.get(str(v), "?"))
        wt = summary_df[summary_df["label"] == "WT"]
        verdict = out_dir / "verdict_2b_vs_2m.csv"
        if wt.empty:
            print("\n[verdict] 没找到 WT 基线(变体名需含 'WT'),跳过 Δ 计算。")
            summary_df.to_csv(verdict, index=False)
        else:
            wt_eq = wt["aim_a1_contacts_eq"].mean()
            summary_df["delta_vs_wt"] = (summary_df["aim_a1_contacts_eq"] - wt_eq).round(2)
            summary_df.to_csv(verdict, index=False)
            print("\n" + "=" * 64)
            print(f"2B vs 2M 判定  (WT 平衡态 AIM-A1 接触 = {wt_eq:.1f})")
            print("  期望: 2B 自抑制松开 → 接触数明显↓ (delta_vs_wt 显著为负);")
            print("        2M 结合面破坏但 AIM 仍在 → 接触下降不明显。")
            print("=" * 64)
            print(f"  {'变体':>10} {'label':>5} {'平衡接触':>8} {'Δvs_WT':>8}")
            for _, r in summary_df.sort_values("delta_vs_wt").iterrows():
                print(f"  {r['variant']:>10} {r['label']:>5} {r.get('aim_a1_contacts_eq', float('nan')):>8} "
                      f"{r.get('delta_vs_wt', float('nan')):>8}")
            for grp in ("2B", "2M", "?"):
                g = summary_df[summary_df["label"] == grp]
                if not g.empty:
                    print(f"  [{grp}] n={len(g)}  平均 Δvs_WT = {g['delta_vs_wt'].mean():+.2f}")
            print(f"\nWritten: {verdict}")
            print("注: 这是接触数的描述性比较; 样本小, 统计显著性 (Mann-Whitney) 留待结果稳定后补。")

    print(f"\nAnalysis complete: {out_dir}")


if __name__ == "__main__":
    main()
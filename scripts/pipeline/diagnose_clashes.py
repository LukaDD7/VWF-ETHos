#!/usr/bin/env python3
"""
diagnose_clashes.py
===================
全原子 steric clash 定位 —— 排查 Boltz-2 结构喂 GROMACS EM 跑飞 (Fmax 1e9) 的根因。

为什么需要它:
  BOLTZ2_AUTINHIB_MD_BLOCKER 报告里用"随机 1500 重原子 + 不含 H"的扫描判定"无 clash",
  但 EM 的 LJ(SR)=+4e7 是典型非键重叠信号。该扫描会漏掉:① 所有 H(pdb2gmx -ignh
  重建的),② 抽样外的原子,③ 一头在样一头不在样的对。本脚本做**全原子、全对**的近
  接触扫描(可读含 H 的 processed 结构),定位真正的 clash 对,并区分 H-clash vs 重原子。

用法:
  # a) 扫一个变体的全部 5 个 Boltz model, 给出 clash 计数 + 推荐最干净的 model
  python3 scripts/pipeline/diagnose_clashes.py \
      --variant-dir output/boltz2_a1_dp_d3_results/boltz_results_VWF_WT_dp_d3_a1

  # b) 细查单个结构 (可含 H: 先 gmx editconf -f processed.gro -o processed.pdb)
  python3 scripts/pipeline/diagnose_clashes.py --input processed.pdb --top 30

依赖: gemmi。clash 定义: 非键(跨残基且非相邻残基)原子对间距 < cutoff。
  - heavy-heavy 默认 cutoff 2.0 Å (真重叠)
  - 含 H 的对默认 cutoff 1.5 Å (重建 H 撞进邻原子)
"""

import argparse
import sys
from pathlib import Path

try:
    import gemmi
except ImportError:
    print("[FATAL] 需要 gemmi: conda activate gromacs (或 boltz2)", file=sys.stderr)
    sys.exit(2)


def is_h(atom):
    return atom.element.name in ("H", "D")


def find_clashes(struct_path, heavy_cut, h_cut, top, skip_adjacent=1):
    st = gemmi.read_structure(str(struct_path))
    if len(st) == 0:
        return None, []
    model = st[0]
    cut = max(heavy_cut, h_cut)
    ns = gemmi.NeighborSearch(model, st.cell, cut).populate()

    clashes = []
    seen = set()
    for ci, chain in enumerate(model):
        for res in chain:
            for atom in res:
                marks = ns.find_atoms(atom.pos, "\0", radius=cut)
                for m in marks:
                    cra = m.to_cra(model)
                    # 同一原子 / 已记录的对 跳过
                    key = tuple(sorted([(ci, res.seqid.num, atom.name),
                                        (m.chain_idx, cra.residue.seqid.num, cra.atom.name)]))
                    if key in seen:
                        continue
                    if cra.residue.seqid.num == res.seqid.num and m.chain_idx == ci and cra.atom.name == atom.name:
                        continue
                    # 排除同链相邻残基(肽键)与同残基(共价键)
                    if m.chain_idx == ci and abs(cra.residue.seqid.num - res.seqid.num) <= skip_adjacent:
                        continue
                    d = atom.pos.dist(cra.atom.pos)
                    involves_h = is_h(atom) or is_h(cra.atom)
                    thr = h_cut if involves_h else heavy_cut
                    if d < thr:
                        seen.add(key)
                        clashes.append((d, involves_h,
                                        f"{chain.name}/{res.name}{res.seqid.num}:{atom.name}",
                                        f"{cra.chain.name}/{cra.residue.name}{cra.residue.seqid.num}:{cra.atom.name}"))
    clashes.sort(key=lambda x: x[0])
    n_h = sum(1 for c in clashes if c[1])
    summary = dict(total=len(clashes), h_involved=n_h, heavy_only=len(clashes) - n_h)
    return summary, clashes[:top]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", help="单个结构文件 (cif/pdb)")
    ap.add_argument("--variant-dir", help="一个 boltz_results_VWF_* 目录, 扫全部 model")
    ap.add_argument("--heavy-cut", type=float, default=2.0, help="重原子 clash 阈值 Å (默认 2.0)")
    ap.add_argument("--h-cut", type=float, default=1.5, help="含 H clash 阈值 Å (默认 1.5)")
    ap.add_argument("--top", type=int, default=20, help="细查模式列出最近的前 N 对")
    args = ap.parse_args()

    if args.variant_dir:
        vd = Path(args.variant_dir)
        cifs = sorted(vd.glob("**/*_model_*.cif"))
        if not cifs:
            print(f"[FATAL] {vd} 下没有 *_model_*.cif", file=sys.stderr); return 1
        print(f"{'model':40s} {'total':>6} {'H-clash':>8} {'heavy':>6} {'min(Å)':>7}")
        print("-" * 72)
        best = None
        for cif in cifs:
            summ, top = find_clashes(cif, args.heavy_cut, args.h_cut, 1)
            mind = top[0][0] if top else float('nan')
            print(f"{cif.name:40s} {summ['total']:6d} {summ['h_involved']:8d} {summ['heavy_only']:6d} {mind:7.2f}")
            score = (summ['heavy_only'], summ['total'])
            if best is None or score < best[0]:
                best = (score, cif.name)
        print("-" * 72)
        print(f"[推荐] 最干净(重原子 clash 最少): {best[1]}")
        print("注: Boltz CIF 不含 H → 这里只反映重原子重叠。H-clash 要在 pdb2gmx 后的")
        print("    processed 结构上查 (gmx editconf -f processed.gro -o processed.pdb, 再 --input 它)。")
        return 0

    if not args.input:
        ap.error("需要 --input 或 --variant-dir")
    summ, top = find_clashes(args.input, args.heavy_cut, args.h_cut, args.top)
    if summ is None:
        print("[FATAL] 读不到结构", file=sys.stderr); return 1
    print(f"[{args.input}]  total={summ['total']}  含H={summ['h_involved']}  纯重原子={summ['heavy_only']}")
    if not top:
        print("  ✅ 未发现 < 阈值的非键近接触 (clash)。Fmax 爆高另有原因(查骨架 ω/键长)。")
        return 0
    print(f"  最近的 {len(top)} 对 (距离升序):")
    for d, h, a, b in top:
        tag = "H-clash" if h else "heavy "
        print(f"    {d:5.2f}Å [{tag}]  {a}  <->  {b}")
    if summ['h_involved'] and summ['heavy_only'] == 0:
        print("\n  → 全是 H-clash: 多半是 pdb2gmx 重建 H 撞邻原子 (坏 ω/rotamer)。")
        print("    受约束 EM(restrain 重原子, 放开 H)应能解掉 → 见 relax_autoinhib_structure.sh")
    elif summ['heavy_only']:
        print("\n  → 有重原子重叠: Boltz 结构本身几何问题。先试换 model, 再 vacuum EM / OpenMM relax。")
    return 0


if __name__ == "__main__":
    sys.exit(main())

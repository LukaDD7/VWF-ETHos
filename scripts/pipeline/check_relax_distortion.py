#!/usr/bin/env python3
"""
check_relax_distortion.py  (闸门 2)
==================================
比较 Boltz 原始结构 vs 弛豫后结构, 量真空/溶剂化 EM 有没有把自抑制几何搞变形。

为什么:真空极小化会塌缩表面盐桥、挪动 D'D3↔A1 相对取向——而那正是 autoinhib MD
要测的量。若弛豫把界面 pose 动大了, 后面的"闭合态稳定性"读数就不可信。

做法:按残基号匹配 CA, Kabsch 叠合, 报全蛋白 CA-RMSD + 位移最大的残基。
把"位移最大残基"与 diagnose_clashes 的"clash 残基"对照, 落在 D'D3-A1 界面就危险。

用法:
  python3 scripts/pipeline/check_relax_distortion.py \
      --orig  <model_2.cif> \
      --relaxed <em_vac.pdb 或 em.pdb>      # .gro 先 gmx editconf 转 pdb
  # 可选 --top 15 ; --rmsd-warn 2.0 (Å, 超过给警告)

依赖: gemmi。退出码: 0 = RMSD 在阈内, 1 = 超阈(需查界面)。
"""
import argparse
import sys

try:
    import gemmi
except ImportError:
    print("[FATAL] 需要 gemmi: conda activate gromacs (或 boltz2)", file=sys.stderr)
    sys.exit(2)


def ca_map(path):
    """{(chain, resnum): (resname, gemmi.Position)} for CA atoms."""
    st = gemmi.read_structure(str(path))
    if len(st) == 0:
        return {}
    out = {}
    for chain in st[0]:
        for res in chain:
            ca = res.find_atom("CA", "*")
            if ca is not None:
                out[(chain.name, res.seqid.num)] = (res.name, ca.pos)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orig", required=True, help="原始 Boltz 结构 (cif/pdb)")
    ap.add_argument("--relaxed", required=True, help="弛豫后结构 (pdb; .gro 先转 pdb)")
    ap.add_argument("--top", type=int, default=15, help="列出位移最大的前 N 个残基")
    ap.add_argument("--rmsd-warn", type=float, default=2.0, help="全 CA-RMSD 超过此值(Å)给警告")
    args = ap.parse_args()

    a = ca_map(args.orig)
    b = ca_map(args.relaxed)
    if not a or not b:
        print("[FATAL] 读不到 CA(检查文件/格式; .gro 需先转 pdb)", file=sys.stderr)
        return 2

    keys = sorted(set(a) & set(b))
    if len(keys) < 3:
        print(f"[FATAL] 匹配到的 CA 太少 ({len(keys)})。残基编号是否对得上?", file=sys.stderr)
        return 2

    pos_a = [a[k][1] for k in keys]     # fixed = 原始
    pos_b = [b[k][1] for k in keys]     # moving = 弛豫后
    sup = gemmi.superpose_positions(pos_a, pos_b)
    tr = sup.transform

    # 叠合后逐残基位移
    disp = []
    for k, pb in zip(keys, pos_b):
        moved = tr.apply(pb)
        d = a[k][1].dist(moved)
        disp.append((d, k, a[k][0]))
    disp.sort(reverse=True)

    print(f"[匹配 CA] {len(keys)}  (orig={len(a)}, relaxed={len(b)})")
    print(f"[全蛋白 CA-RMSD] {sup.rmsd:.3f} Å  (叠合于全部 CA)")
    print(f"\n位移最大的 {args.top} 个残基 (叠合后):")
    print(f"  {'残基':>16} {'位移(Å)':>9}")
    for d, k, rn in disp[:args.top]:
        print(f"  {rn}{k[1]}/{k[0]:>6} {d:9.2f}")

    print("\n→ 把上面这些残基号对照 diagnose_clashes 的 clash 残基, 以及构建里 D'D3 vs A1 的边界:")
    print("  · 落在单结构域内部/外周 loop → 安全(EM 局部松弛)")
    print("  · 落在 D'D3↔A1 界面 → 危险(自抑制 pose 被动过, MD 读数不可信)")

    if sup.rmsd > args.rmsd_warn:
        print(f"\n⚠ 全 CA-RMSD {sup.rmsd:.2f} > {args.rmsd_warn} Å: 弛豫位移偏大。")
        print("  建议改全程受约束弛豫(EM 加 -DPOSRES)后重测。")
        return 1
    print(f"\n✅ 全 CA-RMSD {sup.rmsd:.2f} ≤ {args.rmsd_warn} Å。仍需确认最大位移不在界面。")
    return 0


if __name__ == "__main__":
    sys.exit(main())

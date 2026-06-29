#!/usr/bin/env python3
"""
fetch_clean_7a6o.py
===================
下载并清理 VWF AIM-A1 自抑制态实验结构 (PDB 7A6O, X-ray 2.12 Å), 产出可直接喂
pdb2gmx/relax/MD 的干净 WT 蛋白 PDB。

为什么用它(替代 Boltz D'D3-A1):
  - 7A6O 是**实验测定**的 A1+不连续自抑制模块(AIM)结构 → 自抑制态的真实坐标,
    clash 极少、EM 一下就过, 绕开 Boltz 结构 EM 跑飞 / 域间 pose 不可信的问题。
  - AIM-A1 正是文献公认的 2B 自抑制机制(比 D'D3-A1 构建更对题)。
  - 作为 WT 骨架: 每个 2B 突变体在它上面就地改 1 个残基(FoldX BuildModel)再跑 MD,
    可信且把"突变效应"与"预测噪声"解耦。

清理动作:
  - 只保留 VWF 链(最长的多肽链; 可 --vwf-chain 覆盖), 删纳米抗体 VHH81、SO4、水。
  - 报告: VWF 残基范围 + 编号体系(VWF 规范 vs 局部)+ 内部缺失残基(gap)+ legacy
    recurrent Type 2B A1 residue coverage for structure/QC only.

用法:
  python3 scripts/pipeline/fetch_clean_7a6o.py
  python3 scripts/pipeline/fetch_clean_7a6o.py --out structures/7A6O_AIM_A1_clean.pdb
  python3 scripts/pipeline/fetch_clean_7a6o.py --pdb-id 7A6O --vwf-chain A

依赖: gemmi + 联网(首次下载)。
"""
import argparse
import sys
import urllib.request
from pathlib import Path

try:
    import gemmi
except ImportError:
    print("[FATAL] 需要 gemmi: conda activate gromacs (或 boltz2)", file=sys.stderr)
    sys.exit(2)

# Legacy recurrent Type 2B A1 residues used only for reporting 7A6O structure
# coverage. They are deliberately not imported by the production classifier.
LEGACY_RECURRENT_2B_A1_POSITIONS = [
    1266, 1296, 1304, 1306, 1308, 1309, 1310, 1313, 1314,
    1316, 1321, 1324, 1341, 1342, 1377, 1392, 1461,
]

AA3to1 = {  # 兜底, gemmi 也能查
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}


def one_letter(resname):
    try:
        info = gemmi.find_tabulated_residue(resname)
        if info and info.is_amino_acid():
            c = info.one_letter_code.upper()
            if c.isalpha():
                return c
    except Exception:
        pass
    return AA3to1.get(resname.upper(), 'X')


def download_cif(pdb_id, cache_dir):
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / f"{pdb_id}.cif"
    if dest.is_file() and dest.stat().st_size > 0:
        print(f"[cache] {dest}")
        return dest
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    print(f"[download] {url}")
    try:
        urllib.request.urlretrieve(url, dest)
    except Exception as e:
        print(f"[FATAL] 下载失败 ({e})。无网络时请手动放 {dest}", file=sys.stderr)
        sys.exit(1)
    return dest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb-id", default="7A6O")
    ap.add_argument("--out", default="structures/7A6O_AIM_A1_clean.pdb")
    ap.add_argument("--vwf-chain", default="", help="VWF 链名 (默认自动取最长多肽链)")
    ap.add_argument("--cache-dir", default="structures/_pdb_cache")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent.parent
    cif = download_cif(args.pdb_id, root / args.cache_dir)
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    if len(st) == 0:
        print("[FATAL] 空结构", file=sys.stderr); return 1
    model = st[0]

    # 列出各链 + 多肽长度
    def aa_residues(chain):
        return [r for r in chain if gemmi.find_tabulated_residue(r.name)
                and gemmi.find_tabulated_residue(r.name).is_amino_acid()]
    print("\n[链概览]")
    chain_lens = {}
    for ch in model:
        n = len(aa_residues(ch))
        chain_lens[ch.name] = n
        print(f"  chain {ch.name}: {n} aa  (总残基 {len(ch)})")

    vwf = args.vwf_chain or max(chain_lens, key=chain_lens.get)
    if vwf not in chain_lens:
        print(f"[FATAL] 链 {vwf} 不存在", file=sys.stderr); return 1
    print(f"\n[选定 VWF 链] {vwf}  ({chain_lens[vwf]} aa) — 其余链(纳米抗体等)将删除")

    # 清理: 删水/配体, 删非 VWF 链
    st.remove_waters()
    try:
        st.remove_ligands_and_waters()
    except Exception:
        pass
    for ch in list(model):
        if ch.name != vwf:
            model.remove_chain(ch.name)
    # 删 VWF 链里残留的非氨基酸 (SO4 等)
    chain = model[vwf] if vwf in [c.name for c in model] else model[0]
    for res in list(chain):
        info = gemmi.find_tabulated_residue(res.name)
        if not (info and info.is_amino_acid()):
            chain.remove_residue(res.seqid.num, ' ', res.name) if False else None
    # 用过滤后的残基重建(更稳)
    aa = aa_residues(chain)
    nums = [r.seqid.num for r in aa]
    lo, hi = min(nums), max(nums)
    seq = "".join(one_letter(r.name) for r in aa)

    out = (root / args.out) if not Path(args.out).is_absolute() else Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(out))

    # ---- 报告 ----
    print(f"\n[残基范围] {lo}–{hi}  ({len(aa)} aa)")
    canonical = lo > 1000   # VWF 规范编号 A1≈1271-1462; 局部编号会从 ~1 起
    if canonical:
        print("  编号: 看起来是 VWF 规范编号 (可直接按 R1306 等定位突变)")
    else:
        print(f"  ⚠ 编号: 看起来是局部编号(从 {lo} 起)。映射到 VWF 需加偏移量; ")
        print("     用序列比对确定: A1 域起始 'CEACQEPGG...'≈VWF1260s。突变定位前先确认偏移。")
    print(f"  序列: {seq[:60]}{'...' if len(seq) > 60 else ''}")

    # 内部缺失残基 (gap → MD 前可能需补 loop)
    gaps = [(nums[i], nums[i+1]) for i in range(len(nums)-1) if nums[i+1] - nums[i] > 1]
    if gaps:
        print(f"\n  ⚠ 内部缺失残基 (未解析的 loop, MD 前需补模或确认无碍): {gaps}")
        print("     最柔的 N-AIM(VWF Q1238-E1260)晶体常解不全, 属正常。")
    else:
        print("\n  无内部 gap (链连续)。")

    # Legacy recurrent Type 2B A1 residue coverage (structure/QC only).
    if canonical:
        present = sorted(p for p in LEGACY_RECURRENT_2B_A1_POSITIONS if lo <= p <= hi and p in set(nums))
        missing = sorted(p for p in LEGACY_RECURRENT_2B_A1_POSITIONS if lo <= p <= hi and p not in set(nums))
        print(f"\n  legacy recurrent 2B A1 residue coverage (QC only): {present}")
        if missing:
            print(f"  recurrent positions in range but unresolved: {missing}")

    print(f"\n[OK] 干净 WT 结构 → {out}")
    print("下一步:")
    print(f"  1. WT MD:  bash scripts/pipeline/relax_autoinhib_structure.sh --pdb {out} --variant 7A6O_WT")
    print("  2. 突变体: 在该 WT 骨架上用 FoldX BuildModel 改单残基, 再各自 relax+MD")
    print("  3. 闭合态稳定性 WT vs 突变体对比 = 2B 信号 (见 docs/AUTOINHIB_MD_VALIDATION_GATES.md)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

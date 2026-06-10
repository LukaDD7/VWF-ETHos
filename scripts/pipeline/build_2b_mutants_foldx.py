#!/usr/bin/env python3
"""
build_2b_mutants_foldx.py
=========================
在 7A6O(AIM-A1 实验 WT)骨架上, 用 FoldX BuildModel 批量引入点突变, 产出每个变体的
干净突变体 PDB, 供 relax+MD。把"突变效应"与"de novo 预测噪声"解耦:同一正确 WT 骨架
上只改 1 个残基, WT 与突变体可直接对比。

流程:
  1. (推荐) RepairPDB 修 WT 侧链/小问题
  2. 按 FoldX 格式写 individual_list.txt: <WT><chain><pdb_resnum><MUT>;
  3. BuildModel 逐个建模 → 重命名为 <variant>.pdb

关键: **编号对齐**。FoldX 用 PDB 里的残基号。7A6O 可能是 VWF 规范编号(R1306 直接可用)
或局部编号(需偏移)。脚本用 gemmi 校验每个变体的 WT 残基身份; 可 --detect-offset 自动
找偏移(pdb_resnum = VWF_pos - offset)。WT 身份对不上的变体会被跳过并报原因。

用法:
  python3 scripts/pipeline/build_2b_mutants_foldx.py \
      --wt structures/7A6O_AIM_A1_clean.pdb --foldx /path/to/foldx --detect-offset
  # 自定义变体: --variants "R1306W,V1316M,R1341Q"  或  --variants-file list.txt(每行一个 aa_change)

依赖: gemmi(校验) + FoldX 二进制。
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import gemmi
except ImportError:
    print("[FATAL] 需要 gemmi: conda activate gromacs (或 boltz2)", file=sys.stderr)
    sys.exit(2)

AA3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
          'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
          'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
          'TYR': 'Y', 'VAL': 'V'}

# 默认: A1 域内复发性 2B(+ 少量 2M 对照, 自行用 --variants-file 替换为本地标签集)
DEFAULT_VARIANTS = ["R1306W", "R1306Q", "R1308C", "I1309V", "S1310F", "W1313C",
                    "V1314F", "V1316M", "P1337L", "R1341Q", "R1341W",  # 2B
                    "R1374C", "R1374H", "G1324S"]                      # 2M 对照

VAR_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def parse_variant(s):
    m = VAR_RE.match(s.strip())
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3)


def pdb_resmap(pdb_path, chain_name):
    """{resnum: one_letter} for the chain (first matching if chain_name='')"""
    st = gemmi.read_structure(str(pdb_path))
    model = st[0]
    out = {}
    target = None
    for ch in model:
        if chain_name and ch.name != chain_name:
            continue
        target = ch.name
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid():
                out[res.seqid.num] = AA3to1.get(res.name.upper(), 'X')
        break
    return out, target


def best_offset(variants, resmap):
    """find offset s.t. resmap[pos-offset]==wt for most variants. Returns (offset, n_match)."""
    best = (0, -1)
    cand = set()
    # 候选偏移: 让某变体的 pos 落到 resmap 的某残基
    poss = [p for _, p, _ in variants]
    keys = list(resmap.keys())
    for p in poss:
        for k in keys:
            cand.add(p - k)
    cand.add(0)
    for off in cand:
        n = sum(1 for wt, p, _ in variants
                if (p - off) in resmap and resmap[p - off] == wt)
        if n > best[1]:
            best = (off, n)
    return best


def run_foldx(foldx, args_list, cwd):
    print(f"  $ {foldx} {' '.join(args_list)}")
    r = subprocess.run([foldx] + args_list, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-2000:] + "\n" + r.stderr[-2000:] + "\n")
    return r.returncode == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wt", default="structures/7A6O_AIM_A1_clean.pdb")
    ap.add_argument("--foldx", default=os.environ.get("FOLDX", "foldx"))
    ap.add_argument("--chain", default="", help="VWF 链名(默认第一条多肽链)")
    ap.add_argument("--variants", default="", help="逗号分隔 aa_change, 如 R1306W,V1316M")
    ap.add_argument("--variants-file", default="", help="每行一个 aa_change")
    ap.add_argument("--offset", type=int, default=0, help="pdb_resnum = VWF_pos - offset")
    ap.add_argument("--detect-offset", action="store_true", help="自动探测偏移")
    ap.add_argument("--out-dir", default="structures/7a6o_mutants")
    ap.add_argument("--no-repair", action="store_true", help="跳过 RepairPDB(不推荐)")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent.parent
    wt = (root / args.wt) if not Path(args.wt).is_absolute() else Path(args.wt)
    if not wt.is_file():
        print(f"[FATAL] WT PDB 不存在: {wt}(先跑 fetch_clean_7a6o.py)", file=sys.stderr); return 1
    if shutil.which(args.foldx) is None and not Path(args.foldx).is_file():
        print(f"[FATAL] 找不到 FoldX: {args.foldx}(--foldx 指定, 或 export FOLDX=...)", file=sys.stderr); return 1

    # 变体列表
    raw = []
    if args.variants_file:
        vf = Path(args.variants_file)
        if not vf.is_absolute():
            vf = root / args.variants_file
        raw = [l.strip() for l in open(vf) if l.strip() and not l.startswith("#")]
    elif args.variants:
        raw = [x for x in args.variants.split(",") if x.strip()]
    else:
        raw = DEFAULT_VARIANTS
    variants = []
    for s in raw:
        pv = parse_variant(s)
        if pv:
            variants.append(pv)
        else:
            print(f"  [skip] 解析不了变体: {s}")
    if not variants:
        print("[FATAL] 没有有效变体", file=sys.stderr); return 1

    resmap, chain = pdb_resmap(wt, args.chain)
    if not resmap:
        print("[FATAL] WT PDB 读不到链/残基", file=sys.stderr); return 1
    print(f"[WT] {wt.name}  chain={chain}  残基 {min(resmap)}–{max(resmap)} ({len(resmap)} aa)")

    offset = args.offset
    if args.detect_offset:
        offset, nmatch = best_offset(variants, resmap)
        print(f"[offset] 自动探测 = {offset}  (WT 身份匹配 {nmatch}/{len(variants)})")
    # 校验 WT 身份, 过滤
    good, bad = [], []
    for wtaa, pos, mut in variants:
        rn = pos - offset
        have = resmap.get(rn)
        if have == wtaa:
            good.append((wtaa, pos, mut, rn))
        else:
            bad.append((f"{wtaa}{pos}{mut}", rn, have))
    for name, rn, have in bad:
        print(f"  [skip] {name}: PDB resnum {rn} 处是 {have or '缺失'}, 期望 {name[0]} "
              f"(偏移不对或该残基未解析)")
    if not good:
        print("[FATAL] 没有变体通过 WT 身份校验。检查 --offset/--detect-offset 与 7A6O 编号。", file=sys.stderr)
        return 1

    out = (root / args.out_dir) if not Path(args.out_dir).is_absolute() else Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy(wt, out / "WT.pdb")

    base = "WT"
    if not args.no_repair:
        print("\n[1] RepairPDB")
        if not run_foldx(args.foldx, ["--command=RepairPDB", "--pdb=WT.pdb", f"--output-dir={out}"], out):
            print("[FATAL] RepairPDB 失败", file=sys.stderr); return 1
        base = "WT_Repair"

    # individual_list.txt
    list_path = out / "individual_list.txt"
    with open(list_path, "w") as fh:
        for wtaa, pos, mut, rn in good:
            fh.write(f"{wtaa}{chain}{rn}{mut};\n")
    print(f"\n[2] BuildModel  ({len(good)} 变体)  list={list_path.name}")
    if not run_foldx(args.foldx, ["--command=BuildModel", f"--pdb={base}.pdb",
                                  "--mutant-file=individual_list.txt", f"--output-dir={out}",
                                  "--numberOfRuns=1"], out):
        print("[FATAL] BuildModel 失败", file=sys.stderr); return 1

    # 重命名 <base>_N.pdb → <variant>.pdb (N 按 list 顺序)
    print("\n[3] 重命名输出")
    for i, (wtaa, pos, mut, rn) in enumerate(good, start=1):
        src = out / f"{base}_{i}.pdb"
        dst = out / f"{wtaa}{pos}{mut}.pdb"
        if src.is_file():
            shutil.move(str(src), str(dst))
            print(f"  {wtaa}{pos}{mut}.pdb  (chain {chain} resnum {rn})")
        else:
            print(f"  [WARN] 缺 {src.name}(BuildModel 输出顺序?)")

    print(f"\n[OK] 突变体 → {out}/<variant>.pdb")
    print("下一步(每个变体 relax+MD):")
    print(f"  for v in {out}/*.pdb; do")
    print(f'    bash scripts/pipeline/relax_autoinhib_structure.sh --pdb "$v" --variant "$(basename "${{v%.pdb}}")"')
    print("  done")
    print("WT vs 突变体闭合态稳定性之差 = 2B 信号 (docs/AUTOINHIB_MD_VALIDATION_GATES.md)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

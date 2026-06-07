#!/usr/bin/env python3
"""
extract_aim_autoinhib_features.py
=================================
从 `a1_aim_autoinhibition_context` Boltz-2 预测结构里抽取 **接口级** 自抑制特征,
用于 Type 2B vs 2M 的判别。

为什么不用现成的矩阵列:
  evidence_matrix.csv 里 a1_aim_autoinhibition_context 的 primary_metric 是
  `ptm_or_plddt`(单链构建的全局置信度)。实测已知 2B(n=44)与 2M(n=49)在该
  metric 上几乎完全重叠(delta_vs_wt 中位数 0.067 vs 0.049)——一个全局标量抓不到
  "自抑制松开"。2B 的机制是 AIM(autoinhibitory module)从 A1 的 GPIb 结合面上
  脱离 → 应表现为 **AIM↔A1 接触减少 / GPIb 面暴露**。本脚本直接量这个。

构建几何(来自 yaml: AIM-flanked A1 monomeric construct, vwf_range 1234-1493):
  - 单链, Boltz 输出 CIF 残基从 1 开始 → construct_local = VWF_pos - 1233。
  - AIM 片段(VWF 编号, 见 agentic_vwf_classifier.FUNCTIONAL_SITES):
        N-term AIM 1238-1268, C-term AIM 1460-1472。
  - A1 核心 = 其余球状部分(默认 1271-1459, 两段 AIM 之间)。
  - 特征 = AIM 重原子与 A1 核心重原子在 cutoff(默认 4.5 Å)内的 **残基对接触数**。

输出 (CSV, 默认 output/aim_autoinhib_features.csv):
  variant_id, aa_change, n_models, aim_a1_contacts,
  aim_a1_contacts_wt, aim_a1_contacts_delta_vs_wt,
  aim_a1_contacts_zscore, aim_release_score
  其中 aim_release_score = -zscore(contacts)  → 越大 = AIM 越松开 = 越像 2B。

依赖: gemmi (gromacs/boltz2 env 自带)。
用法 (CPU 实例 / A40 CPU 侧均可):
  python3 scripts/pipeline/extract_aim_autoinhib_features.py \
      --results-dir output/boltz2_vwd_functional_panel/boltz_results \
      --output output/aim_autoinhib_features.csv

注意: 残基编号假设 Boltz CIF 从 1 起。若不是,用 --construct-offset 调整
(construct_local = VWF_pos - offset, 默认 offset=1233)。脚本会打印 WT 的
接触残基对数量,便于人工 sanity-check。
"""

import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict

try:
    import gemmi
except ImportError:
    print("[FATAL] 需要 gemmi。conda activate gromacs (或 boltz2) 后重试。", file=sys.stderr)
    sys.exit(2)

ASSAY = "a1_aim_autoinhibition_context"


def vwf_ranges_to_local(ranges, offset):
    return [(a - offset, b - offset) for a, b in ranges]


def in_ranges(num, ranges):
    return any(a <= num <= b for a, b in ranges)


def heavy_atoms(residue):
    for atom in residue:
        el = atom.element.name
        if el != "H" and el != "D":
            yield atom


def count_aim_a1_contacts(cif_path, aim_local, a1_local, cutoff):
    """返回 AIM 残基与 A1 核心残基在 cutoff 内的残基对接触数。"""
    st = gemmi.read_structure(str(cif_path))
    if len(st) == 0:
        return None
    model = st[0]
    # 取第一条(也是唯一一条)蛋白链
    chain = None
    for ch in model:
        chain = ch
        break
    if chain is None:
        return None

    aim_res, a1_res = [], []
    for res in chain:
        num = res.seqid.num
        if in_ranges(num, aim_local):
            aim_res.append(res)
        elif in_ranges(num, a1_local):
            a1_res.append(res)
    if not aim_res or not a1_res:
        return 0

    # NeighborSearch 仅在 A1 核心原子上建索引, 再用 AIM 原子查询
    ns = gemmi.NeighborSearch(model, st.cell, cutoff).populate()
    a1_keys = {(res.seqid.num) for res in a1_res}
    cutoff2 = cutoff * cutoff
    contact_pairs = set()
    for res in aim_res:
        for atom in heavy_atoms(res):
            marks = ns.find_atoms(atom.pos, "\0", radius=cutoff)
            for m in marks:
                cra = m.to_cra(model)
                if cra.residue.seqid.num in a1_keys and cra.atom.element.name not in ("H", "D"):
                    if atom.pos.dist(cra.atom.pos) <= cutoff:
                        contact_pairs.add((res.seqid.num, cra.residue.seqid.num))
    return len(contact_pairs)


def iter_autoinhib_jobs(results_dir):
    """yield (job_name, [cif paths]) for the autoinhibition assay."""
    jobs = defaultdict(list)
    for cif in results_dir.glob(f"**/predictions/*{ASSAY}*/*_model_*.cif"):
        jobs[cif.parent.name].append(cif)
    return jobs


def job_to_variant(job_name):
    """'VWF_R1306W__a1_aim_autoinhibition_context' -> ('VWF_R1306W','R1306W')."""
    vid = job_name.split("__")[0]
    aa = vid[len("VWF_"):] if vid.startswith("VWF_") else vid
    return vid, aa


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="output/boltz2_vwd_functional_panel/boltz_results")
    ap.add_argument("--output", default="output/aim_autoinhib_features.csv")
    ap.add_argument("--cutoff", type=float, default=4.5, help="接触距离阈值 Å (默认 4.5)")
    ap.add_argument("--construct-offset", type=int, default=1233,
                    help="construct_local = VWF_pos - offset (默认 1233)")
    ap.add_argument("--aim-ranges", default="1238-1268,1460-1472",
                    help="AIM VWF 编号区间, 逗号分隔 (默认 N/C 端 AIM)")
    ap.add_argument("--a1-core", default="1271-1459",
                    help="A1 核心 VWF 编号区间 (默认两段 AIM 之间的球状体)")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent.parent
    results_dir = (root / args.results_dir) if not Path(args.results_dir).is_absolute() else Path(args.results_dir)
    out_path = (root / args.output) if not Path(args.output).is_absolute() else Path(args.output)

    def parse_ranges(s):
        out = []
        for part in s.split(","):
            a, b = part.split("-")
            out.append((int(a), int(b)))
        return out

    aim_local = vwf_ranges_to_local(parse_ranges(args.aim_ranges), args.construct_offset)
    a1_local = vwf_ranges_to_local(parse_ranges(args.a1_core), args.construct_offset)

    if not results_dir.is_dir():
        print(f"[FATAL] results dir 不存在: {results_dir}", file=sys.stderr)
        return 1
    jobs = iter_autoinhib_jobs(results_dir)
    if not jobs:
        print(f"[FATAL] 在 {results_dir} 下没找到 {ASSAY} 的 CIF。", file=sys.stderr)
        return 1

    print(f"[INFO] assay={ASSAY}  jobs={len(jobs)}  cutoff={args.cutoff}Å")
    print(f"[INFO] AIM(local)={aim_local}  A1core(local)={a1_local}")

    # 每个 variant 多 model 取平均接触数
    rows = []
    for job_name, cifs in sorted(jobs.items()):
        vid, aa = job_to_variant(job_name)
        counts = []
        for cif in cifs:
            c = count_aim_a1_contacts(cif, aim_local, a1_local, args.cutoff)
            if c is not None:
                counts.append(c)
        if not counts:
            continue
        mean_c = sum(counts) / len(counts)
        rows.append({"variant_id": vid, "aa_change": aa,
                     "n_models": len(counts), "aim_a1_contacts": round(mean_c, 3)})

    if not rows:
        print("[FATAL] 没有任何可用结构。", file=sys.stderr)
        return 1

    # WT 基线 + 全 panel zscore
    wt = next((r for r in rows if r["aa_change"].upper() in ("WT", "WILDTYPE")), None)
    wt_c = wt["aim_a1_contacts"] if wt else None
    vals = [r["aim_a1_contacts"] for r in rows]
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / len(vals)
    std = var ** 0.5 or 1.0

    for r in rows:
        z = (r["aim_a1_contacts"] - mean) / std
        r["aim_a1_contacts_wt"] = wt_c if wt_c is not None else ""
        r["aim_a1_contacts_delta_vs_wt"] = round(r["aim_a1_contacts"] - wt_c, 3) if wt_c is not None else ""
        r["aim_a1_contacts_zscore"] = round(z, 3)
        r["aim_release_score"] = round(-z, 3)  # 越大 = AIM 越松开 = 越像 2B

    out_path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["variant_id", "aa_change", "n_models", "aim_a1_contacts", "aim_a1_contacts_wt",
            "aim_a1_contacts_delta_vs_wt", "aim_a1_contacts_zscore", "aim_release_score"]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    print(f"[OK] {len(rows)} variants -> {out_path}")
    if wt:
        print(f"[sanity] WT ({wt['variant_id']}) AIM-A1 接触残基对 = {wt_c} "
              f"(若为 0 多半是残基编号假设错了, 调 --construct-offset)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

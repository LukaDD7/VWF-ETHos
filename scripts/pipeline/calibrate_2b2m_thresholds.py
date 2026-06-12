#!/usr/bin/env python3
"""
calibrate_2b2m_thresholds.py
============================
用**更大的干净标签集** (output/labeled_variants_all.csv, 由 build_labeled_variant_set.py
从 /Volumes/LQ1000/VWD 文献表清出) + 功能 panel 特征
(output/boltz2_vwd_functional_panel/evidence_matrix.csv 的
a1_gpiba_forced_binding / a1_heparan_sulfate_binding / a1_aim_autoinhibition_context
的 within-assay z-score), 重跑分类器 LOF 轴阈值校准。

为什么这么做:
  RULE6 的 LOF 轴 (2M=结合面丧失) 之前只在原始 44/49 集上调过 LOF_COMBINED_Z=-0.75。
  现在有了更大的独立标签集 → 验证 2B/2M 在更大样本上还分不分得开, 阈值要不要动。
  这条线**不依赖 MD** (读 panel 静态特征), 能立刻把静态分类器往前推一格。

做的事:
  1. 干净标签: labeled_variants_all.csv 里 A1 域 (1262-1466)、单标签、2B 或 2M 的变体。
  2. join 到 evidence_matrix 的 (wt_aa, position, mut_aa) → 取 fb / heparan / aim z-score。
  3. 一致性核对: 新干净标签 vs evidence_matrix 自带的 is_Type2B/is_Type2M。
  4. 阈值扫描: fb 单轴 / heparan 单轴 / mean(fb,heparan) 组合轴, 每个阈值算
     2M recall (真 2M 判 LOF) 与 2B 误判率 (真 2B 误判 LOF), 净收益 = recall - fp。
  5. 输出 CSV 到 output/ + 控制台报告。判 LOF 方向: z ≤ T (结合下降 → LOF → 2M)。

输出:
  output/calib_2b2m_threshold_sweep.csv  — 每轴每阈值的 recall/fp/net
  output/calib_2b2m_joined.csv           — join 后逐变体 (标签 + 三轴 z + 一致性)
用法:
  python3 scripts/pipeline/calibrate_2b2m_thresholds.py
"""
import argparse
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
A1_LO, A1_HI = 1262, 1466

FB = "a1_gpiba_forced_binding__primary_zscore_within_assay"
HEP = "a1_heparan_sulfate_binding__primary_zscore_within_assay"
AIM = "a1_aim_autoinhibition_context__primary_zscore_within_assay"


def fnum(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (TypeError, ValueError):
        return None


def load_clean_labels(path):
    """A1 域内、单标签 2B/2M → {(wt,pos,mut): '2B'|'2M'}"""
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            lab = row["labels"]
            if lab not in ("2B", "2M"):
                continue
            if row["conflict"].lower() == "true":
                continue
            pos = int(row["position"])
            if not (A1_LO <= pos <= A1_HI):
                continue
            out[(row["wt_aa"], pos, row["mut_aa"])] = lab
    return out


def load_panel(path):
    """(wt,pos,mut) → dict(fb,hep,aim,is2b,is2m)。position 从 aa_change_3letter 或 variant_id 解析。"""
    import re
    rx = re.compile(r"([A-Z])(\d{2,4})([A-Z])")
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            w, m = row.get("wt_aa", ""), row.get("mut_aa", "")
            pos = None
            for key in ("variant_id", "aa_change_3letter"):
                s = row.get(key, "") or ""
                mm = rx.search(s.replace("p.", ""))
                if mm:
                    pos = int(mm.group(2))
                    if not w:
                        w = mm.group(1)
                    if not m:
                        m = mm.group(3)
                    break
            if pos is None or not w or not m:
                continue
            out[(w, pos, m)] = {
                "fb": fnum(row.get(FB)),
                "hep": fnum(row.get(HEP)),
                "aim": fnum(row.get(AIM)),
                "is2b": (row.get("is_Type2B", "").lower() in ("true", "1")),
                "is2m": (row.get("is_Type2M", "").lower() in ("true", "1")),
            }
    return out


def sweep(joined, axis, thresholds):
    """axis: 'fb'|'hep'|'combined'。判 LOF = (z ≤ T)。返回 [(T, recall2m, fp2b, net)]。"""
    res = []
    for T in thresholds:
        tp = fn = fp = tn = 0
        for j in joined:
            if axis == "combined":
                if j["fb"] is None or j["hep"] is None:
                    continue
                z = (j["fb"] + j["hep"]) / 2.0
            else:
                z = j[axis]
                if z is None:
                    continue
            lof = z <= T
            if j["label"] == "2M":
                tp += lof
                fn += (not lof)
            else:  # 2B
                fp += lof
                tn += (not lof)
        n2m, n2b = tp + fn, fp + tn
        recall = tp / n2m if n2m else 0.0
        fprate = fp / n2b if n2b else 0.0
        res.append((T, n2m, n2b, recall, fprate, recall - fprate))
    return res


def median(xs):
    xs = sorted(v for v in xs if v is not None)
    if not xs:
        return None
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default=str(ROOT / "output/labeled_variants_all.csv"))
    ap.add_argument("--panel", default=str(ROOT / "output/boltz2_vwd_functional_panel/evidence_matrix.csv"))
    ap.add_argument("--out-dir", default=str(ROOT / "output"))
    args = ap.parse_args()

    labels = load_clean_labels(args.labels)
    panel = load_panel(args.panel)

    joined, miss = [], 0
    for key, lab in labels.items():
        if key in panel:
            p = panel[key]
            joined.append({"key": key, "label": lab, **p})
        else:
            miss += 1

    n2b = sum(1 for j in joined if j["label"] == "2B")
    n2m = sum(1 for j in joined if j["label"] == "2M")

    # 一致性核对
    agree = sum(1 for j in joined
                if (j["label"] == "2B" and j["is2b"]) or (j["label"] == "2M" and j["is2m"]))
    disagree = len(joined) - agree

    # 中位数 (方向感)
    def med(lab, ax):
        return median([j[ax] for j in joined if j["label"] == lab])

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    thr = [round(x, 2) for x in [-0.25, -0.5, -0.75, -1.0, -1.25, -1.5]]
    axes = {"fb": "fb_single", "hep": "heparan_single", "combined": "mean(fb,heparan)"}
    sweep_rows = []
    for ax, name in axes.items():
        for (T, m2m, m2b, rec, fpr, net) in sweep(joined, ax, thr):
            sweep_rows.append({"axis": name, "threshold": T, "n_2M": m2m, "n_2B": m2b,
                               "recall_2M": round(rec, 3), "fp_2B": round(fpr, 3),
                               "net_benefit": round(net, 3)})

    with open(out / "calib_2b2m_threshold_sweep.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["axis", "threshold", "n_2M", "n_2B",
                                           "recall_2M", "fp_2B", "net_benefit"])
        w.writeheader(); w.writerows(sweep_rows)

    with open(out / "calib_2b2m_joined.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["aa_change", "label", "fb_z", "heparan_z", "aim_z",
                                           "combined_z", "panel_is2b", "panel_is2m", "label_agrees"])
        w.writeheader()
        for j in joined:
            wt, pos, mut = j["key"]
            comb = (j["fb"] + j["hep"]) / 2.0 if (j["fb"] is not None and j["hep"] is not None) else None
            agr = (j["label"] == "2B" and j["is2b"]) or (j["label"] == "2M" and j["is2m"])
            w.writerow({"aa_change": f"{wt}{pos}{mut}", "label": j["label"],
                        "fb_z": j["fb"], "heparan_z": j["hep"], "aim_z": j["aim"],
                        "combined_z": round(comb, 3) if comb is not None else "",
                        "panel_is2b": j["is2b"], "panel_is2m": j["is2m"], "label_agrees": agr})

    # ---- 报告 ----
    print("=" * 64)
    print(f"干净 A1 单标签: 2B={sum(1 for v in labels.values() if v=='2B')} "
          f"2M={sum(1 for v in labels.values() if v=='2M')}  "
          f"(join 上 panel: 2B={n2b} 2M={n2m}, 未匹配 {miss})")
    print(f"标签一致性 (新干净标签 vs panel is_Type2B/2M): agree={agree} disagree={disagree}")
    print(f"\n中位 z (方向感, 越负=结合越弱):")
    print(f"          fb        heparan   aim")
    for lab in ("2B", "2M"):
        print(f"  {lab}:  {med(lab,'fb'):>7.2f}   {med(lab,'hep'):>7.2f}   {med(lab,'aim'):>7.2f}")
    print(f"\n阈值扫描 (判 LOF=z≤T → 2M; recall=真2M被抓, fp=真2B误判):")
    for ax, name in axes.items():
        print(f"\n  [{name}]")
        print(f"    T       recall_2M   fp_2B   net")
        best = None
        for (T, m2m, m2b, rec, fpr, net) in sweep(joined, ax, thr):
            mark = ""
            if best is None or net > best[1]:
                best = (T, net)
            print(f"    {T:>5}    {rec:>6.2f}     {fpr:>5.2f}   {net:>5.2f}{mark}")
        print(f"    → 最佳 net @ T={best[0]} (net={best[1]:.2f})")
    print("=" * 64)
    print(f"输出 → {out}/calib_2b2m_threshold_sweep.csv, calib_2b2m_joined.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

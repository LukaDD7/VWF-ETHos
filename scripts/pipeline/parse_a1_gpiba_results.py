#!/usr/bin/env python3
"""
parse_a1_gpiba_results.py
=========================
解析 Boltz-2 对 VWF A1 + GPIbα 复合物的预测结果。

从每个 job 的 confidence_model_*.json 中提取 iPTM，
计算 ΔiPTM = iPTM(mut) - iPTM(WT)，输出分型建议。

输入目录结构（run_a1_gpiba_boltz2.sh 的输出）：
  results_dir/
    VWF_WT_vs_GPIb_alpha/
      predictions/
        confidence_model_0.json
        confidence_model_1.json
        ...
      .done
    VWF_R1306W_vs_GPIb_alpha/
      predictions/
        confidence_model_0.json
        ...
      .done

用法：
  # 解析全部结果
  python3 parse_a1_gpiba_results.py \\
      --results-dir output/boltz2_a1_gpiba_results \\
      --output output/boltz2_a1_gpiba_analysis/iptm_results.csv

  # 实时进度检查（无需结果完成）
  python3 parse_a1_gpiba_results.py --check-progress \\
      --results-dir output/boltz2_a1_gpiba_results
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Optional

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
    print("[WARN] pandas not found. CSV output disabled. pip install pandas")


# =============================================================================
# 已知分型（用于校准验证）
# =============================================================================

KNOWN_VARIANTS = {
    "VWF_R1306W": "2B_GOF",   # 经典 Type 2B
    "VWF_R1306Q": "2B_GOF",   # 经典 Type 2B
    "VWF_R1341Q": "2B_GOF",   # Type 2B
    "VWF_V1316M": "2B_GOF",   # Type 2B
    "VWF_F1369I": "2M_LOF",   # Type 2M
    "VWF_R1374C": "2M_LOF",   # Type 2M (A1-LOF)
    "VWF_I1372S": "2M_LOF",   # Type 2M
}

# iPTM 阈值（参考 BOLTZ2_SPEC.md）
THRESHOLD_GOF = +0.05   # ΔiPTM > +0.05 → 2B
THRESHOLD_LOF = -0.05   # ΔiPTM < -0.05 → 2M


# =============================================================================
# 核心解析函数
# =============================================================================

def parse_confidence_json(json_path: Path) -> Optional[float]:
    """从单个 confidence_model_X.json 中提取 iPTM。"""
    try:
        with open(json_path) as f:
            data = json.load(f)
        iptm = data.get("iptm")
        if iptm is not None:
            return float(iptm)
        # 兼容旧格式：某些版本写在 complex_plddt 或 interface_ptm
        for key in ["interface_ptm", "i_ptm", "complex_iptm"]:
            if key in data:
                return float(data[key])
        return None
    except Exception as e:
        return None


def get_best_iptm(job_dir: Path) -> tuple[Optional[float], int]:
    """
    从 job_dir/predictions/ 中的所有 confidence_model_*.json 提取最佳 iPTM。
    返回 (best_iptm, n_samples_found)
    """
    pred_dir = job_dir / "predictions"
    if not pred_dir.exists():
        return None, 0

    conf_files = sorted(pred_dir.glob("confidence_model_*.json"))
    if not conf_files:
        return None, 0

    best_iptm = None
    n_found = 0
    for cf in conf_files:
        iptm = parse_confidence_json(cf)
        if iptm is not None:
            n_found += 1
            if best_iptm is None or iptm > best_iptm:
                best_iptm = iptm

    return best_iptm, n_found


def scan_results(results_dir: Path) -> list[dict]:
    """
    扫描结果目录，提取所有已完成 job 的 iPTM。
    返回记录列表，每个记录包含 job_name, iptm_best, n_samples。
    """
    records = []
    job_dirs = sorted(
        [d for d in results_dir.iterdir()
         if d.is_dir() and not d.name.startswith("_")]
    )

    for job_dir in job_dirs:
        job_name = job_dir.name
        done_marker = job_dir / ".done"
        is_done = done_marker.exists()

        best_iptm, n_samples = get_best_iptm(job_dir)

        records.append({
            "job_name":  job_name,
            "iptm_best": best_iptm,
            "n_samples": n_samples,
            "is_done":   is_done,
        })

    return records


def classify_variant(delta_iptm: Optional[float]) -> str:
    """根据 ΔiPTM 给出分型建议。"""
    if delta_iptm is None:
        return "unclassified"
    if delta_iptm > THRESHOLD_GOF:
        return "2B_GOF"
    if delta_iptm < THRESHOLD_LOF:
        return "2M_LOF"
    return "neutral_VUS"


# =============================================================================
# 进度检查模式
# =============================================================================

def check_progress(results_dir: Path, total_expected: int = 74):
    """快速显示当前进度（无需全部完成）。"""
    results_dir = results_dir.resolve()
    if not results_dir.exists():
        print(f"[ERROR] Results dir not found: {results_dir}")
        return

    done_jobs = list(results_dir.glob("*/.done"))
    all_dirs = [d for d in results_dir.iterdir()
                if d.is_dir() and not d.name.startswith("_")]

    total_dirs = len(all_dirs)
    total_done = len(done_jobs)
    remaining = total_expected - total_done

    pct = total_done / total_expected * 100 if total_expected > 0 else 0
    bar_len = 40
    filled = int(bar_len * total_done / total_expected) if total_expected > 0 else 0
    bar = "█" * filled + "░" * (bar_len - filled)

    print()
    print("=" * 60)
    print("  VWF A1+GPIbα Boltz-2 Progress")
    print("=" * 60)
    print(f"  [{bar}]  {pct:.1f}%")
    print(f"  Completed  : {total_done:>4} / {total_expected}")
    print(f"  Remaining  : {remaining:>4}")

    # 速度估算（基于 .done 文件的时间戳）
    done_files = sorted(
        results_dir.glob("*/.done"),
        key=lambda f: f.stat().st_mtime
    )
    if len(done_files) >= 2:
        import time
        oldest_t = done_files[0].stat().st_mtime
        newest_t = done_files[-1].stat().st_mtime
        elapsed_h = (newest_t - oldest_t) / 3600
        if elapsed_h > 0.001:
            rate = len(done_files) / elapsed_h
            eta_h = remaining / rate if rate > 0 else float("inf")
            print(f"  Speed      : {rate:.1f} jobs/hour")
            print(f"  ETA        : {eta_h:.1f} hours")

    # 失败检测（有 predictions/ 但没有 .done 的目录）
    dirty = []
    for d in all_dirs:
        if not (d / ".done").exists() and (d / "predictions").exists():
            cifs = list((d / "predictions").glob("*.cif"))
            if cifs:
                dirty.append(d.name)
    if dirty:
        print(f"\n  Dirty jobs (crashed, no .done): {len(dirty)}")
        for name in dirty[:5]:
            print(f"    - {name}")
        if len(dirty) > 5:
            print(f"    ... and {len(dirty)-5} more")
        print("  Clean with: grep FAIL output/boltz2_a1_gpiba_results/worker_*.log")

    print("=" * 60)
    print()


# =============================================================================
# 主解析流程
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--results-dir", default="../../output/boltz2_a1_gpiba_results",
        help="Boltz-2 output directory"
    )
    parser.add_argument(
        "--output", default="../../output/boltz2_a1_gpiba_analysis/iptm_results.csv",
        help="Output CSV path"
    )
    parser.add_argument(
        "--check-progress", action="store_true",
        help="Only show progress, don't parse results"
    )
    parser.add_argument(
        "--total-expected", type=int, default=74,
        help="Total expected jobs for progress display"
    )
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    results_dir = (script_dir / args.results_dir).resolve()

    if args.check_progress:
        check_progress(results_dir, args.total_expected)
        return

    # ---------- 解析结果 -------------------------------------------------------
    print("=" * 70)
    print("VWF A1+GPIbα Results Parser")
    print("=" * 70)

    if not results_dir.exists():
        print(f"[ERROR] Results dir not found: {results_dir}")
        sys.exit(1)

    records = scan_results(results_dir)
    done_records = [r for r in records if r["is_done"] and r["iptm_best"] is not None]
    pending = [r for r in records if not r["is_done"]]

    print(f"Completed jobs with iPTM: {len(done_records)}")
    print(f"Pending / no data: {len(pending)}")

    if not done_records:
        print("[ERROR] No completed jobs found. Run run_a1_gpiba_boltz2.sh first.")
        return

    # ---------- WT 基准 --------------------------------------------------------
    wt_record = next((r for r in done_records
                      if r["job_name"] == "VWF_WT_vs_GPIb_alpha"), None)
    if wt_record is None:
        print("[ERROR] WT reference job not found or not completed.")
        print("  Expected: results_dir/VWF_WT_vs_GPIb_alpha/")
        wt_iptm = None
    else:
        wt_iptm = wt_record["iptm_best"]
        print(f"\nWT iPTM: {wt_iptm:.4f} (n_samples={wt_record['n_samples']})")

    # ---------- 计算 ΔiPTM ----------------------------------------------------
    results = []
    for r in done_records:
        if r["job_name"] == "VWF_WT_vs_GPIb_alpha":
            continue

        jname = r["job_name"]
        # 从 job 名提取 variant_id
        m = re.match(r"^(VWF_.+)_vs_GPIb_alpha$", jname)
        variant_id = m.group(1) if m else jname

        iptm = r["iptm_best"]
        delta_iptm = (iptm - wt_iptm) if (wt_iptm is not None and iptm is not None) else None
        subtype = classify_variant(delta_iptm)

        results.append({
            "job_name":         jname,
            "variant_id":       variant_id,
            "iptm_best":        round(iptm, 5) if iptm else None,
            "iptm_wt":          round(wt_iptm, 5) if wt_iptm else None,
            "delta_iptm":       round(delta_iptm, 5) if delta_iptm is not None else None,
            "predicted_subtype": subtype,
            "n_samples":        r["n_samples"],
        })

    # 按 delta_iptm 降序排列（GOF 在前）
    results.sort(key=lambda x: x["delta_iptm"] if x["delta_iptm"] is not None else 0,
                 reverse=True)

    # ---------- 打印摘要 -------------------------------------------------------
    print()
    print(f"{'Variant':<25} {'iPTM':>7} {'ΔiPTM':>8} {'Predicted':>12}")
    print("-" * 58)
    for row in results:
        iptm_str = f"{row['iptm_best']:.4f}" if row['iptm_best'] else "  N/A"
        diptm_str = f"{row['delta_iptm']:+.4f}" if row['delta_iptm'] is not None else "   N/A"
        flag = "★" if row["variant_id"] in KNOWN_VARIANTS else ""
        print(f"{row['variant_id']:<25} {iptm_str:>7} {diptm_str:>8} {row['predicted_subtype']:>12} {flag}")

    # ---------- 校准验证 -------------------------------------------------------
    print()
    print("=" * 58)
    print("CALIBRATION CHECK (known variants):")
    print("=" * 58)
    calibration_pass = 0
    calibration_total = 0
    for row in results:
        vid = row["variant_id"]
        if vid not in KNOWN_VARIANTS:
            continue
        known = KNOWN_VARIANTS[vid]
        predicted = row["predicted_subtype"]
        delta = row["delta_iptm"]

        correct = (known == predicted) or (
            known == "2B_GOF" and predicted == "2B_GOF" or
            known == "2M_LOF" and predicted == "2M_LOF"
        )
        calibration_total += 1
        if correct:
            calibration_pass += 1

        delta_str = f"{delta:+.4f}" if delta is not None else "N/A"
        status = "✓" if correct else "✗"
        print(f"  {status} {vid:<20} known={known:<10} pred={predicted:<12} ΔiPTM={delta_str}")

    if calibration_total > 0:
        acc = calibration_pass / calibration_total * 100
        print()
        print(f"  Calibration accuracy: {calibration_pass}/{calibration_total} ({acc:.0f}%)")
        if acc < 60:
            print("  [WARN] Low accuracy. Consider adjusting thresholds:")
            print(f"         THRESHOLD_GOF={THRESHOLD_GOF}  THRESHOLD_LOF={THRESHOLD_LOF}")
            print("         Or check if WT reference iPTM is reasonable (expected: 0.5-0.8)")
    else:
        print("  (no known calibration variants found in results)")

    # ---------- 保存 CSV -------------------------------------------------------
    if HAS_PANDAS:
        out_path = (script_dir / args.output).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)

        df = pd.DataFrame(results)
        # 同时追加 WT 行
        if wt_record:
            df_wt = pd.DataFrame([{
                "job_name": "VWF_WT_vs_GPIb_alpha",
                "variant_id": "VWF_WT",
                "iptm_best": round(wt_iptm, 5),
                "iptm_wt": round(wt_iptm, 5),
                "delta_iptm": 0.0,
                "predicted_subtype": "WT_reference",
                "n_samples": wt_record["n_samples"],
            }])
            df = pd.concat([df_wt, df], ignore_index=True)

        df.to_csv(out_path, index=False)
        print()
        print(f"Results saved: {out_path}")
        print(f"  {len(df)} rows | columns: {list(df.columns)}")

        # 分型汇总
        print()
        print("Subtype distribution:")
        for subtype, count in df["predicted_subtype"].value_counts().items():
            print(f"  {subtype:<18}: {count:>3}")
    else:
        print("\n[WARN] pandas not installed — CSV not saved. pip install pandas")
        print("Results (JSON format):")
        import json as json_mod
        print(json_mod.dumps(results[:5], indent=2))

    print()
    print("Done.")


if __name__ == "__main__":
    main()

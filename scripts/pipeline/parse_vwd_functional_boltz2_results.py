#!/usr/bin/env python3
"""
parse_vwd_functional_boltz2_results.py
=======================================
解析 VWD/VWF functional Boltz-2 panel 结果。

用法：
    python scripts/pipeline/parse_vwd_functional_boltz2_results.py \\
        --results-dir output/boltz2_vwd_functional_panel/boltz_results \\
        --manifest output/boltz2_vwd_functional_panel/job_manifest.csv \\
        --output output/boltz2_vwd_functional_panel/boltz_results_summary.csv
"""

import argparse
import json
import re
import sys
from pathlib import Path
from collections import defaultdict

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
    print("[WARN] pandas not found. CSV output disabled.")


def parse_confidence_json(json_path):
    """从 confidence_*.json 提取关键指标。"""
    try:
        with open(json_path) as f:
            data = json.load(f)
        return {
            "iptm": data.get("iptm"),
            "ptm": data.get("ptm"),
            "complex_plddt": data.get("complex_plddt"),
            "complex_iplddt": data.get("complex_iplddt"),
            "complex_pde": data.get("complex_pde"),
            "complex_ipde": data.get("complex_ipde"),
            "protein_iptm": data.get("protein_iptm"),
            "ligand_iptm": data.get("ligand_iptm"),
        }
    except Exception as e:
        return None


def extract_job_name(conf_path):
    """从 confidence 文件路径提取 job name。"""
    # 路径: .../predictions/VWF_WT__dprime_d3_fviii_binding/confidence_VWF_WT__dprime_d3_fviii_binding_model_0.json
    parts = conf_path.parts
    pred_idx = parts.index("predictions") if "predictions" in parts else -1
    if pred_idx >= 0 and pred_idx + 1 < len(parts):
        return parts[pred_idx + 1]
    return conf_path.parent.name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="output/boltz2_vwd_functional_panel/boltz_results")
    parser.add_argument("--manifest", default="output/boltz2_vwd_functional_panel/job_manifest.csv")
    parser.add_argument("--output", default="output/boltz2_vwd_functional_panel/boltz_results_summary.csv")
    args = parser.parse_args()

    # results_dir: output/boltz2_vwd_functional_panel/boltz_results/
    # from scripts/pipeline/, go up 2 levels to project root, then down
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent  # up to projects/VWF-ETHos
    results_dir = project_root / args.results_dir
    manifest_path = project_root / args.manifest

    print("=" * 70)
    print("VWD/VWF Functional Boltz-2 Panel Results Parser")
    print("=" * 70)

    # 加载 job manifest
    manifest_df = None
    if HAS_PANDAS and manifest_path.exists():
        manifest_df = pd.read_csv(manifest_path)
        print(f"Loaded manifest: {len(manifest_df)} jobs")

    # 扫描所有 confidence 文件
    conf_files = sorted(results_dir.glob("**/confidence_*.json"))
    print(f"Found {len(conf_files)} confidence files")

    # 按 job 分组
    job_data = defaultdict(list)
    for conf_file in conf_files:
        job_name = extract_job_name(conf_file)
        job_data[job_name].append(conf_file)

    print(f"Unique jobs: {len(job_data)}")

    # 计算每个 job 的最佳 iptm（多 sample 取平均或最佳）
    job_results = []
    for job_name, conf_files_list in sorted(job_data.items()):
        iptm_values = []
        for cf in conf_files_list:
            data = parse_confidence_json(cf)
            if data and data["iptm"] is not None:
                iptm_values.append(data["iptm"])

        if iptm_values:
            avg_iptm = sum(iptm_values) / len(iptm_values)
            best_iptm = max(iptm_values)
            n_samples = len(iptm_values)
        else:
            avg_iptm = best_iptm = n_samples = None

        # 从 manifest 查找额外信息
        variant_id = None
        domain = None
        assay_key = None
        clinical_axis = None

        if manifest_df is not None:
            row = manifest_df[manifest_df["job_name"] == job_name]
            if len(row) > 0:
                variant_id = row.iloc[0].get("variant_id", None)
                domain = row.iloc[0].get("inferred_domain", None)
                assay_key = row.iloc[0].get("assay_key", None)
                clinical_axis = row.iloc[0].get("clinical_axis", None)

        job_results.append({
            "job_name": job_name,
            "variant_id": variant_id,
            "domain": domain,
            "assay_key": assay_key,
            "clinical_axis": clinical_axis,
            "n_samples": n_samples,
            "avg_iptm": round(avg_iptm, 5) if avg_iptm is not None else None,
            "best_iptm": round(best_iptm, 5) if best_iptm is not None else None,
        })

    # 打印摘要
    print()
    print(f"{'Job Name':<50} {'Avg iPTM':>10} {'Best':>10} {'Samples':>8}")
    print("-" * 82)
    for r in job_results[:30]:
        avg = f"{r['avg_iptm']:.4f}" if r['avg_iptm'] else "N/A"
        best = f"{r['best_iptm']:.4f}" if r['best_iptm'] else "N/A"
        n = r['n_samples'] or 0
        print(f"{r['job_name']:<50} {avg:>10} {best:>10} {n:>8}")
    if len(job_results) > 30:
        print(f"... and {len(job_results) - 30} more jobs")

    # 统计
    valid_results = [r for r in job_results if r['avg_iptm'] is not None]
    print()
    print(f"Total jobs: {len(job_results)}")
    print(f"Jobs with results: {len(valid_results)}")
    if valid_results:
        iptm_values = [r['avg_iptm'] for r in valid_results if r['avg_iptm']]
        print(f"iPTM range: {min(iptm_values):.4f} - {max(iptm_values):.4f}")
        print(f"iPTM mean: {sum(iptm_values)/len(iptm_values):.4f}")

    # 按 assay_key 分组统计
    if HAS_PANDAS:
        df = pd.DataFrame(job_results)
        if "assay_key" in df.columns and df["assay_key"].notna().any():
            print()
            print("iPTM by assay_key:")
            grouped = df.groupby("assay_key")["avg_iptm"].agg(["mean", "std", "count"])
            for assay, row in grouped.iterrows():
                print(f"  {assay:<30}: mean={row['mean']:.4f}, std={row['std']:.4f}, n={int(row['count'])}")

        # 保存 CSV
        out_path = project_root / args.output
        out_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_path, index=False)
        print()
        print(f"Results saved: {out_path}")

    print()
    print("Done.")


if __name__ == "__main__":
    main()
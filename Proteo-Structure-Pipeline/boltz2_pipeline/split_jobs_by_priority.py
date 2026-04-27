#!/usr/bin/env python3
"""
Split Boltz-2 Jobs by Priority
================================
将已生成的 823 个 job 按诊疗优先级拆分为两组：

  Priority 1 (2B / 2M)：先跑，临床紧迫性最高
    - A1 × GPIb_alpha          → 2B vs 2M-A1 鉴别
    - A1 × Heparan_Sulfate_mimic → 2B 调节环境
    - A3 × Collagen_I_THP      → 2M-A3 胶原结合
    - C1/C2 × Collagen_I_THP   → 2M 延伸胶原区

  Priority 2 (2A / 2N / Other)：第二批
    - A2 × ADAMTS13_Spacer     → 2A
    - D'/D3 × FVIII_LightChain → 2N
    - C4 × Integrin_αIIbβ3    → 聚集功能

输出：
  output/boltz2_blind_scan/
    priority1_2B2M/
      batch_p1_001.json ... batch_p1_NNN.json
      job_manifest_priority1.csv
    priority2_2A2N/
      batch_p2_001.json ... batch_p2_NNN.json
      job_manifest_priority2.csv

用法：
  python split_jobs_by_priority.py
  python split_jobs_by_priority.py --batch-size 30 --input-dir ../../output/boltz2_blind_scan
"""

import argparse
import json
import glob
import os
import sys
from pathlib import Path

import pandas as pd

# ── 优先级定义 ────────────────────────────────────────────────────────────────
# key: (domain, ligand_key)  → priority group
PRIORITY_RULES = {
    # Priority 1：2B / 2M
    ("A1", "GPIb_alpha"):              1,
    ("A1", "Heparan_Sulfate_mimic"):   1,
    ("A3", "Collagen_I_THP"):          1,
    ("C1", "Collagen_I_THP"):          1,
    ("C2", "Collagen_I_THP"):          1,
    # Priority 2：2A / 2N / Other
    ("A2", "ADAMTS13_Spacer"):         2,
    ("D_prime", "FVIII_LightChain"):   2,
    ("D3",      "FVIII_LightChain"):   2,
    ("C4", "Integrin_alphaIIb_beta3"): 2,
}

PRIORITY_NAMES = {
    1: "priority1_2B2M",
    2: "priority2_2A2N",
}
PRIORITY_LABELS = {
    1: "Priority 1 — 2B / 2M (GPIbα, Collagen)",
    2: "Priority 2 — 2A / 2N / Other (ADAMTS13, FVIII, Integrin)",
}


def load_all_jobs(input_dir: Path) -> dict:
    """
    读取所有 batch_*.json，返回 {job_name: job_dict} 映射。
    """
    all_jobs = {}
    for batch_file in sorted(input_dir.glob("batch_*.json")):
        data = json.load(open(batch_file))
        for job in data.get("jobs", []):
            all_jobs[job["name"]] = job
    return all_jobs


def assign_priority(row: pd.Series) -> int:
    key = (row["domain"], row["ligand_key"])
    return PRIORITY_RULES.get(key, 0)   # 0 = unassigned (shouldn't happen)


def write_batches(jobs: list, output_dir: Path, prefix: str, batch_size: int):
    output_dir.mkdir(parents=True, exist_ok=True)
    n_batches = (len(jobs) + batch_size - 1) // batch_size
    for i in range(n_batches):
        chunk = jobs[i * batch_size: (i + 1) * batch_size]
        out_path = output_dir / f"{prefix}_{i+1:03d}.json"
        json.dump(
            {"name": f"{prefix}_{i+1:03d}", "jobs": chunk},
            open(out_path, "w"),
            indent=2
        )
    return n_batches


def main():
    parser = argparse.ArgumentParser(description="Split Boltz-2 jobs by priority")
    parser.add_argument(
        "--input-dir", type=str,
        default="../../output/boltz2_blind_scan",
    )
    parser.add_argument("--batch-size", type=int, default=30)
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    input_dir  = (script_dir / args.input_dir).resolve()
    manifest_path = input_dir / "job_manifest.csv"

    if not manifest_path.exists():
        print(f"[ERROR] Manifest not found: {manifest_path}")
        print("  Run prepare_boltz2_inputs.py first.")
        sys.exit(1)

    print("=" * 65)
    print("VWF Boltz-2 — Priority Split")
    print("=" * 65)
    print(f"Input dir : {input_dir}")
    print(f"Batch size: {args.batch_size}")
    print()

    # 1. 读 manifest
    manifest = pd.read_csv(manifest_path)
    manifest["priority"] = manifest.apply(assign_priority, axis=1)

    # 2. 加载所有 job dict
    print("Loading job definitions from batch JSONs...")
    all_job_dicts = load_all_jobs(input_dir)
    print(f"  Loaded {len(all_job_dicts)} jobs.")
    print()

    # 3. 按优先级分组
    for p in [1, 2]:
        sub = manifest[manifest["priority"] == p].copy()
        label = PRIORITY_LABELS[p]
        folder_name = PRIORITY_NAMES[p]
        output_dir = input_dir / folder_name

        print(f"{'─'*65}")
        print(f"{label}")
        print(f"  Jobs  : {len(sub)}")
        print(f"  Output: {output_dir}")

        if sub.empty:
            print("  [SKIP] No jobs in this priority group.")
            continue

        # 统计明细
        breakdown = sub.groupby(["domain", "ligand_key"])["job_name"].count()
        for (dom, lig), cnt in breakdown.items():
            print(f"    {dom:12} × {lig:30}: {cnt:>4} jobs")

        # 收集 job dict（保持原始顺序）
        jobs_in_group = []
        missing = 0
        for jname in sub["job_name"]:
            if jname in all_job_dicts:
                jobs_in_group.append(all_job_dicts[jname])
            else:
                missing += 1
        if missing:
            print(f"  [WARN] {missing} jobs not found in batch JSONs (may be leftover test jobs).")

        # 写 batch 文件
        prefix = f"batch_p{p}"
        n_batches = write_batches(jobs_in_group, output_dir, prefix, args.batch_size)

        # 写 manifest
        sub_manifest_path = output_dir / f"job_manifest_priority{p}.csv"
        sub.to_csv(sub_manifest_path, index=False)

        print(f"  Batches : {n_batches} files ({prefix}_001.json ... {prefix}_{n_batches:03d}.json)")
        print(f"  Manifest: {sub_manifest_path.name}")
        print()

    # 4. 整体统计
    print("=" * 65)
    print("Summary")
    print("=" * 65)
    total = len(manifest)
    for p in [1, 2]:
        n = (manifest["priority"] == p).sum()
        print(f"  {PRIORITY_LABELS[p]}")
        print(f"    → {n} jobs  ({n/total*100:.1f}%)")
    unassigned = (manifest["priority"] == 0).sum()
    if unassigned:
        print(f"  [!] Unassigned: {unassigned} jobs (check PRIORITY_RULES)")

    print()
    print("Next steps:")
    print(f"  1. Transfer '{PRIORITY_NAMES[1]}/' to GPU server first")
    print(f"     ./run_boltz2.sh --input-dir output/boltz2_blind_scan/{PRIORITY_NAMES[1]}")
    print(f"  2. After p1 completes, transfer '{PRIORITY_NAMES[2]}/'")
    print(f"     ./run_boltz2.sh --input-dir output/boltz2_blind_scan/{PRIORITY_NAMES[2]}")


if __name__ == "__main__":
    main()

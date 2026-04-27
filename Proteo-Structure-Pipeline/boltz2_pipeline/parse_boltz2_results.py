#!/usr/bin/env python3
"""
VWF Boltz-2 Results Parser & Visualizer
========================================
解析 Boltz-2 输出目录，提取亲和力评分和结构置信度，
生成用于"盲扫式"差异诊断的热图和对比图表。

输出文件：
  output/boltz2_analysis/
    affinity_scores.csv         ← 所有 job 的亲和力 + 置信度汇总
    heatmap_affinity.png        ← 变异 × 配体热图（亲和力）
    heatmap_plddt.png           ← 变异 × 配体热图（pLDDT）
    delta_vs_wt.png             ← 每个变异相对 WT 的亲和力变化（瀑布图）
    subtype_prediction.csv      ← 基于亲和力差异的自动分型建议

用法：
  python parse_boltz2_results.py \\
    --results-dir ../../output/boltz2_results \\
    --manifest ../../output/boltz2_blind_scan/job_manifest.csv \\
    --output-dir ../../output/boltz2_analysis
"""

import argparse
import json
import re
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import seaborn as sns
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False
    print("[WARN] matplotlib/seaborn not available. Plots will be skipped.")


# =============================================================================
# 1. 结果扫描与解析
# =============================================================================

def parse_confidence_json(json_path: Path) -> Dict:
    """解析 Boltz-2 confidence_model_X.json，提取关键指标。"""
    try:
        with open(json_path) as f:
            data = json.load(f)
        return {
            "plddt_mean":  float(np.mean(data.get("plddt", [0]))),
            "pae_mean":    float(np.mean(data.get("pae", [[0]]))) if "pae" in data else None,
            "ptm":         data.get("ptm", None),
            "iptm":        data.get("iptm", None),  # inter-chain pTM (복合物 핵심 지표)
        }
    except Exception as e:
        return {"plddt_mean": None, "pae_mean": None, "ptm": None, "iptm": None}


def parse_affinity_json(json_path: Path) -> Optional[float]:
    """
    解析 Boltz-2 v2 的亲和力输出（如存在）。
    返回预测 ΔG (kcal/mol) 或 None。
    """
    if not json_path.exists():
        return None
    try:
        with open(json_path) as f:
            data = json.load(f)
        # Boltz-2 输出格式（根据官方文档）
        # data["affinity"]["deltaG"] 或 data["affinity"]["log_kd"]
        aff = data.get("affinity", {})
        if "deltaG" in aff:
            return float(aff["deltaG"])
        if "log_kd" in aff:
            return float(aff["log_kd"]) * 1.363  # kT 换算近似
        return None
    except Exception:
        return None


def scan_results_directory(results_dir: Path) -> pd.DataFrame:
    """
    递归扫描 Boltz-2 输出目录，提取所有 job 的预测结果。
    目录结构假设：
      results_dir/
        batch_001/
          VWF_R1306W_vs_GPIb_alpha/
            predictions/
              confidence_model_0.json
              affinity_model_0.json   (可选)
        ...
    """
    records = []

    for batch_dir in sorted(results_dir.glob("batch_*")):
        for job_dir in sorted(batch_dir.iterdir()):
            if not job_dir.is_dir():
                continue
            job_name = job_dir.name
            pred_dir = job_dir / "predictions"
            if not pred_dir.exists():
                continue

            # 找最高编号模型（取最优 diffusion sample）
            conf_files = sorted(pred_dir.glob("confidence_model_*.json"))
            if not conf_files:
                continue

            # 多个 diffusion sample → 取 pLDDT 最高的
            best_conf = None
            best_plddt = -1.0
            best_model_idx = 0

            for cf in conf_files:
                m = re.search(r"confidence_model_(\d+)\.json", cf.name)
                if not m:
                    continue
                idx = int(m.group(1))
                conf = parse_confidence_json(cf)
                plddt = conf.get("plddt_mean") or -1.0
                if plddt > best_plddt:
                    best_plddt = plddt
                    best_conf = conf
                    best_model_idx = idx

            # 亲和力
            aff_file = pred_dir / f"affinity_model_{best_model_idx}.json"
            delta_g = parse_affinity_json(aff_file)

            # 解析 job_name → variant_id + ligand_key
            parts = job_name.rsplit("_vs_", 1)
            variant_id = parts[0] if parts else job_name
            ligand_key = parts[1] if len(parts) == 2 else "unknown"

            records.append({
                "job_name":     job_name,
                "variant_id":   variant_id,
                "ligand_key":   ligand_key,
                "plddt_mean":   best_conf.get("plddt_mean") if best_conf else None,
                "pae_mean":     best_conf.get("pae_mean") if best_conf else None,
                "iptm":         best_conf.get("iptm") if best_conf else None,
                "ptm":          best_conf.get("ptm") if best_conf else None,
                "delta_g":      delta_g,
                "best_model":   best_model_idx,
            })

    df = pd.DataFrame(records)
    print(f"Parsed {len(df)} completed jobs from {results_dir}")
    return df


# =============================================================================
# 2. 与 WT 对比，计算差异指标
# =============================================================================

def compute_delta_vs_wt(df: pd.DataFrame) -> pd.DataFrame:
    """
    计算每个突变体相对于 WT 的亲和力和 pLDDT 变化。
    WT job 名格式：VWF_WT_vs_<ligand>
    """
    wt_rows = df[df["variant_id"] == "VWF_WT"].set_index("ligand_key")
    results = []

    for _, row in df.iterrows():
        lig = row["ligand_key"]
        if row["variant_id"] == "VWF_WT":
            continue
        wt = wt_rows.loc[lig] if lig in wt_rows.index else None
        
        d_plddt = (row["plddt_mean"] - wt["plddt_mean"]) if wt is not None and row["plddt_mean"] and wt["plddt_mean"] else None
        d_iptm  = (row["iptm"] - wt["iptm"]) if wt is not None and row["iptm"] and wt["iptm"] else None
        d_dg    = (row["delta_g"] - wt["delta_g"]) if wt is not None and row["delta_g"] and wt["delta_g"] else None

        results.append({
            **row.to_dict(),
            "delta_plddt": d_plddt,
            "delta_iptm":  d_iptm,
            "delta_deltaG":d_dg,
        })

    return pd.DataFrame(results)


# =============================================================================
# 3. 自动分型建议
# =============================================================================

SUBTYPE_RULES = {
    # ligand_key       : (增强/减弱, 分型建议, 解释)
    "GPIb_alpha":         ("increase", "2B",
                           "GPIbα 结合力增强 → AIM 破坏 → 自发血小板黏附 (Type 2B GOF)"),
    "GPIb_alpha":         ("decrease", "2M",
                           "GPIbα 结合力减弱 → A1 功能缺陷 (Type 2M A1-defect)"),
    "Collagen_I_THP":     ("decrease", "2M",
                           "胶原结合力减弱 → A3 功能缺陷 (Type 2M A3-defect)"),
    "ADAMTS13_Spacer":    ("increase", "2A",
                           "ADAMTS13 可及性增强 → HMW multimers 过度裂解 (Type 2A)"),
    "FVIII_LightChain":   ("decrease", "2N",
                           "FVIII 结合力减弱 → FVIII 半衰期↓ (Type 2N, 类 Hemophilia A)"),
}

def predict_subtype(row: pd.Series) -> Tuple[str, str]:
    """
    基于亲和力变化（delta_iptm 或 delta_deltaG）给出分型建议。
    返回 (suggested_subtype, reasoning)
    """
    ligand = row.get("ligand_key", "")
    d_iptm = row.get("delta_iptm")
    d_dg   = row.get("delta_deltaG")

    # 优先使用 delta_iptm（结构置信度），其次 delta_deltaG
    delta = d_iptm if d_iptm is not None else d_dg
    if delta is None:
        return "unclassified", "No affinity data"

    direction = "increase" if delta > 0.05 else ("decrease" if delta < -0.05 else "neutral")

    # GPIb_alpha 特殊：方向决定 2B 还是 2M
    if ligand == "GPIb_alpha":
        if direction == "increase":
            return "2B", SUBTYPE_RULES["GPIb_alpha"][2]
        elif direction == "decrease":
            return "2M", "GPIbα 结合力减弱 → A1 功能缺陷 (Type 2M A1-defect)"
    
    for key, (dir_, subtype, reason) in SUBTYPE_RULES.items():
        if key == ligand and dir_ == direction:
            return subtype, reason

    return "VUS", f"Direction={direction} on {ligand}; no clear subtype rule matched"


# =============================================================================
# 4. 热图绘制
# =============================================================================

def plot_heatmap(
    df: pd.DataFrame,
    value_col: str,
    title: str,
    output_path: Path,
    center: float = 0.0,
    cmap: str = "RdBu_r",
):
    """绘制变异 × 配体热图（支持缺失值）。"""
    if not HAS_PLOT:
        print(f"  [SKIP] {output_path.name} — matplotlib not available")
        return

    pivot = df.pivot_table(
        index="variant_id", columns="ligand_key", values=value_col, aggfunc="mean"
    )
    
    if pivot.empty:
        print(f"  [SKIP] {output_path.name} — no data to plot")
        return

    # 按最大值排序（突出异常变异）
    pivot = pivot.reindex(pivot.abs().max(axis=1).sort_values(ascending=False).index)

    n_rows, n_cols = pivot.shape
    figw = max(8, n_cols * 1.5)
    figh = max(6, n_rows * 0.35)

    fig, ax = plt.subplots(figsize=(figw, figh))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap=cmap,
        center=center,
        annot=(n_rows <= 50),  # 行数少时显示数值
        fmt=".2f",
        linewidths=0.4,
        linecolor="#cccccc",
        mask=pivot.isna(),
        cbar_kws={"label": value_col},
    )
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel("Ligand", fontsize=11)
    ax.set_ylabel("VWF Variant", fontsize=11)
    plt.xticks(rotation=30, ha="right", fontsize=10)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_waterfall(df: pd.DataFrame, value_col: str, ligand_key: str, output_path: Path):
    """针对单个配体，绘制所有变异的亲和力变化瀑布图（降序）。"""
    if not HAS_PLOT:
        return
    sub = df[df["ligand_key"] == ligand_key].dropna(subset=[value_col]).copy()
    if sub.empty:
        return
    sub = sub.sort_values(value_col, ascending=False)

    figw = max(10, len(sub) * 0.3)
    fig, ax = plt.subplots(figsize=(figw, 5))
    colors = ["#d62728" if v > 0 else "#1f77b4" for v in sub[value_col]]
    ax.bar(sub["variant_id"], sub[value_col], color=colors, edgecolor="white", linewidth=0.3)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_title(f"Δ Affinity vs WT — Ligand: {ligand_key}", fontsize=13, fontweight="bold")
    ax.set_xlabel("Variant", fontsize=11)
    ax.set_ylabel(f"Δ {value_col}", fontsize=11)
    plt.xticks(rotation=90, fontsize=7)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# =============================================================================
# 5. 主函数
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Parse Boltz-2 results and generate diagnostic visualizations"
    )
    parser.add_argument("--results-dir", type=str,
                        default="../../output/boltz2_results")
    parser.add_argument("--manifest", type=str,
                        default="../../output/boltz2_blind_scan/job_manifest.csv")
    parser.add_argument("--output-dir", type=str,
                        default="../../output/boltz2_analysis")
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    results_dir  = (script_dir / args.results_dir).resolve()
    manifest_path = (script_dir / args.manifest).resolve()
    output_dir   = (script_dir / args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("VWF Boltz-2 Results Parser")
    print("=" * 70)

    # Step 1: 扫描结果
    df = scan_results_directory(results_dir)
    if df.empty:
        print("[ERROR] No results found. Make sure run_boltz2.sh has completed.")
        return

    # Step 2: 合并 manifest（补充 position / domain 等元数据）
    if manifest_path.exists():
        manifest = pd.read_csv(manifest_path)
        df = df.merge(manifest[["job_name", "position", "domain", "wt_aa", "mut_aa"]],
                      on="job_name", how="left")
    
    # Step 3: 计算 delta vs WT
    df_delta = compute_delta_vs_wt(df)

    # Step 4: 保存汇总表
    score_path = output_dir / "affinity_scores.csv"
    df.to_csv(score_path, index=False)
    print(f"\nScores saved: {score_path}")

    delta_path = output_dir / "delta_vs_wt.csv"
    df_delta.to_csv(delta_path, index=False)
    print(f"Delta scores saved: {delta_path}")

    # Step 5: 自动分型建议
    if not df_delta.empty:
        df_delta[["suggested_subtype", "reasoning"]] = df_delta.apply(
            lambda r: pd.Series(predict_subtype(r)), axis=1
        )
        sub_path = output_dir / "subtype_prediction.csv"
        df_delta[["variant_id", "ligand_key", "delta_iptm", "delta_deltaG",
                  "suggested_subtype", "reasoning"]].to_csv(sub_path, index=False)
        print(f"Subtype predictions saved: {sub_path}")

    # Step 6: 热图
    print("\nGenerating plots...")
    if not df_delta.empty and "delta_iptm" in df_delta.columns:
        plot_heatmap(
            df_delta, "delta_iptm",
            "VWF Variant × Ligand — Δ iPTM (vs WT)\n"
            "Red = Gain-of-function   |   Blue = Loss-of-function",
            output_dir / "heatmap_delta_iptm.png",
            center=0.0, cmap="RdBu_r"
        )
    
    plot_heatmap(
        df, "plddt_mean",
        "VWF Variant × Ligand — pLDDT (Structural Confidence)",
        output_dir / "heatmap_plddt.png",
        center=70.0, cmap="RdYlGn"
    )

    if "delta_iptm" in df_delta.columns:
        for lig in df_delta["ligand_key"].dropna().unique():
            plot_waterfall(
                df_delta, "delta_iptm", lig,
                output_dir / f"waterfall_{lig}.png"
            )

    print("\n" + "=" * 70)
    print("Analysis complete.")
    print(f"All outputs in: {output_dir}")


if __name__ == "__main__":
    main()

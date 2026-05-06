#!/usr/bin/env python3
"""
FoldX A1-GPIbα ΔΔG_bind Runner
================================
模板结构: 1M10.pdb  (VWF A1 R1306Q + GPIbα)
目标    : 对所有A1域患者突变计算 ΔΔG_bind(GPIbα)
输出    : foldx_a1_gpiba.csv  (variant_id, ddG_bind, direction, status)

逻辑:
  1M10 的背景突变: Chain A 位置543 = Q  (UniProt R1306Q)
  每个患者突变的分析流程:
    - WT参考:  仅做 Q543→R (恢复野生型)
    - 患者突变: 同时做 Q543→R + 患者突变 (双突变回到WT再引入变体)
    - 注意: 若患者突变本身就在1306位, 则只需单突变 Q543→mut
  ΔΔG_bind = G_bind(mutant) - G_bind(WT_pseudo)

用法:
  # Mac本地测试 (无FoldX, 只验证逻辑):
  python foldx_a1_gpiba_runner.py --dry-run

  # 服务器执行:
  python foldx_a1_gpiba_runner.py \
      --variants ../../data/processed/master_type1_type2.csv \
      --pdb ../../structures/1M10.pdb \
      --foldx-bin /inspire/hdd/global_user/mengweicheng-240108120092/lzy/tools/foldx \
      --output ../../output/foldx_a1_gpiba.csv \
      --workers 16
"""

import os
import re
import sys
import csv
import glob
import shutil
import tempfile
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Optional

# ── 常量 ──────────────────────────────────────────────────────────────────────

# 1M10.pdb 中的背景突变: UniProt R1306 → PDB Chain A Q543
BG_CHAIN        = "A"
BG_PDB_POS      = 543
BG_PDB_WT       = "Q"   # 模板中的氨基酸
BG_UNIPROT_WT   = "R"   # 真正的野生型

# UniProt 位置 → PDB 位置 的偏移 (从SEQADV读取: UniProt 1306 - PDB 543 = 763)
UNIPROT_TO_PDB_OFFSET = 763

# 1M10 Chain A 的 PDB 残基范围
PDB_CHAIN_A_START = 505
PDB_CHAIN_A_END   = 703

# FoldX AnalyseComplex 的链配置 (A=VWF A1, B=GPIbα)
ANALYSE_CHAINS = "A,B"

# ΔΔG_bind 判断阈值 (kcal/mol)
THRESHOLD_GOF = -1.5   # < threshold → 结合增强 → 2B
THRESHOLD_LOF = +1.5   # > threshold → 结合减弱 → 2M

AA_3TO1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
}

# ── 数据类 ────────────────────────────────────────────────────────────────────

@dataclass
class Variant:
    variant_id: str
    uniprot_pos: int
    wt_aa: str
    mut_aa: str
    subtype: str
    pdb_pos: int = field(init=False)

    def __post_init__(self):
        self.pdb_pos = self.uniprot_pos - UNIPROT_TO_PDB_OFFSET

    @property
    def in_pdb_range(self) -> bool:
        return PDB_CHAIN_A_START <= self.pdb_pos <= PDB_CHAIN_A_END

    @property
    def is_synonymous(self) -> bool:
        return self.wt_aa == self.mut_aa

    @property
    def is_stop(self) -> bool:
        return self.mut_aa in ('X', '*')

    @property
    def at_background_site(self) -> bool:
        return self.pdb_pos == BG_PDB_POS

    def foldx_wt_mutation(self) -> str:
        """仅恢复背景突变的 FoldX 符号 (WT reference)"""
        return f"{BG_PDB_WT}{BG_CHAIN}{BG_PDB_POS}{BG_UNIPROT_WT}"

    def foldx_patient_mutation(self) -> str:
        """患者突变的 FoldX 符号 (含背景恢复)"""
        if self.at_background_site:
            # 患者突变就在1306位: Q543→mut_aa (一步到位)
            return f"{BG_PDB_WT}{BG_CHAIN}{BG_PDB_POS}{self.mut_aa}"
        else:
            # 先恢复背景Q→R, 再引入患者突变
            patient_mut = f"{self.wt_aa}{BG_CHAIN}{self.pdb_pos}{self.mut_aa}"
            return f"{self.foldx_wt_mutation()},{patient_mut}"

# ── 工具函数 ──────────────────────────────────────────────────────────────────

def load_variants(csv_path: str) -> list[Variant]:
    """从 master_type1_type2.csv 读取 A1 域突变"""
    variants = []
    seen = set()
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get('Domain_Name') != 'A1':
                continue
            vid = row['Variant_ID']
            pos = int(row['Position'])
            wt  = row['WT_AA'].strip().upper()
            mut = row['Mut_AA'].strip().upper()
            sub = row.get('Subtype', '').strip()

            key = (vid, pos, mut)
            if key in seen:
                continue
            seen.add(key)

            v = Variant(variant_id=vid, uniprot_pos=pos, wt_aa=wt, mut_aa=mut, subtype=sub)
            variants.append(v)
    return variants


def get_residue_at_pos(pdb_path: str, chain: str, pos: int) -> Optional[str]:
    """从PDB读取指定链/位置的实际氨基酸(验证用)"""
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[21] == chain and int(line[22:26].strip()) == pos and line[12:16].strip() == "CA":
                res3 = line[17:20].strip()
                return AA_3TO1.get(res3, res3)
    return None


def run_cmd(cmd: list, cwd: str, log_file: str) -> int:
    """执行命令, 输出写到 log_file"""
    with open(log_file, 'a') as lf:
        result = subprocess.run(cmd, cwd=cwd, stdout=lf, stderr=lf)
    return result.returncode


def parse_analyse_complex(fxout_path: str) -> Optional[float]:
    """
    解析 FoldX AnalyseComplex 输出文件, 提取 Interaction Energy (kcal/mol)
    输出文件格式: Interaction_<pdbname>_AC.fxout
    关键行: Group1\tGroup2\t...\tInteraction Energy\t...
    """
    if not os.path.exists(fxout_path):
        return None
    with open(fxout_path) as f:
        lines = f.readlines()

    # 找到数据行 (跳过注释和表头)
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        # AnalyseComplex输出: PDB, Group1, Group2, IntraclashesGroup1, IntraclashesGroup2, Interaction Energy, ...
        if len(parts) >= 6:
            try:
                return float(parts[5])   # Interaction Energy column
            except ValueError:
                continue
    return None


def repair_pdb(foldx_bin: str, pdb_path: str, work_dir: str) -> Optional[str]:
    """
    FoldX RepairPDB: 修复PDB结构
    返回修复后的PDB路径 (work_dir/<name>_Repair.pdb)
    """
    pdb_name = Path(pdb_path).stem
    repaired = os.path.join(work_dir, f"{pdb_name}_Repair.pdb")
    if os.path.exists(repaired):
        return repaired

    shutil.copy(pdb_path, work_dir)
    cmd = [
        foldx_bin, "--command=RepairPDB",
        f"--pdb={Path(pdb_path).name}",
        "--quiet",
    ]
    rotabase = os.path.join(os.path.dirname(foldx_bin), "rotabase.txt")
    if os.path.exists(rotabase):
        cmd.append(f"--rotabaseLocation={rotabase}")

    log = os.path.join(work_dir, "repair.log")
    rc = run_cmd(cmd, work_dir, log)
    if rc != 0 or not os.path.exists(repaired):
        return None
    return repaired


def build_and_analyse(
    foldx_bin: str,
    repaired_pdb: str,
    mutation_str: str,
    label: str,
    work_dir: str,
) -> Optional[float]:
    """
    FoldX BuildModel → AnalyseComplex, 返回 Interaction Energy (kcal/mol)
    label: 用于命名临时目录和日志
    """
    job_dir = os.path.join(work_dir, label)
    os.makedirs(job_dir, exist_ok=True)

    pdb_name = Path(repaired_pdb).stem
    pdb_filename = Path(repaired_pdb).name
    shutil.copy(repaired_pdb, job_dir)

    # Step 1: individual_list.txt
    indlist = os.path.join(job_dir, "individual_list.txt")
    with open(indlist, 'w') as f:
        f.write(mutation_str + ";\n")   # FoldX 格式: 每行一组突变, 分号结尾

    rotabase = os.path.join(os.path.dirname(foldx_bin), "rotabase.txt")

    def rotabase_arg():
        return [f"--rotabaseLocation={rotabase}"] if os.path.exists(rotabase) else []

    # Step 2: BuildModel
    cmd_build = [
        foldx_bin, "--command=BuildModel",
        f"--pdb={pdb_filename}",
        f"--mutant-file=individual_list.txt",
        "--numberOfRuns=1",
        "--quiet",
    ] + rotabase_arg()

    log = os.path.join(job_dir, "build.log")
    rc = run_cmd(cmd_build, job_dir, log)
    if rc != 0:
        return None

    # BuildModel 输出: <name>_Repair_1.pdb (第一个突变体)
    mutant_pdb = os.path.join(job_dir, f"{pdb_name}_1.pdb")
    if not os.path.exists(mutant_pdb):
        # 也可能叫 <name>_Repair.pdb_1 或其他, glob 一下
        candidates = glob.glob(os.path.join(job_dir, f"{pdb_name}*_1.pdb"))
        if not candidates:
            return None
        mutant_pdb = candidates[0]

    mutant_filename = Path(mutant_pdb).name

    # Step 3: AnalyseComplex
    cmd_ac = [
        foldx_bin, "--command=AnalyseComplex",
        f"--pdb={mutant_filename}",
        f"--analyseComplexChains={ANALYSE_CHAINS}",
        "--quiet",
    ] + rotabase_arg()

    log_ac = os.path.join(job_dir, "analyse.log")
    rc = run_cmd(cmd_ac, job_dir, log_ac)
    if rc != 0:
        return None

    # 查找输出文件
    mutant_stem = Path(mutant_pdb).stem
    fxout = os.path.join(job_dir, f"Interaction_{mutant_stem}_AC.fxout")
    if not os.path.exists(fxout):
        candidates = glob.glob(os.path.join(job_dir, "Interaction_*_AC.fxout"))
        if not candidates:
            return None
        fxout = candidates[0]

    return parse_analyse_complex(fxout)


def process_variant(args):
    """多进程worker: 处理单个变异"""
    foldx_bin, repaired_pdb, variant, wt_energy, work_dir = args
    v = variant

    try:
        mut_energy = build_and_analyse(
            foldx_bin, repaired_pdb,
            mutation_str=v.foldx_patient_mutation(),
            label=v.variant_id,
            work_dir=work_dir,
        )
        if mut_energy is None or wt_energy is None:
            return dict(
                variant_id=v.variant_id, uniprot_pos=v.uniprot_pos,
                wt_aa=v.wt_aa, mut_aa=v.mut_aa, subtype=v.subtype,
                g_bind_wt=wt_energy, g_bind_mut=None, ddg_bind=None,
                direction="unknown", status="FAIL_FOLDX",
                foldx_notation=v.foldx_patient_mutation(),
            )
        ddg = mut_energy - wt_energy
        if ddg < THRESHOLD_GOF:
            direction = "GOF"
        elif ddg > THRESHOLD_LOF:
            direction = "LOF"
        else:
            direction = "neutral"
        return dict(
            variant_id=v.variant_id, uniprot_pos=v.uniprot_pos,
            wt_aa=v.wt_aa, mut_aa=v.mut_aa, subtype=v.subtype,
            g_bind_wt=round(wt_energy, 3), g_bind_mut=round(mut_energy, 3),
            ddg_bind=round(ddg, 3), direction=direction, status="OK",
            foldx_notation=v.foldx_patient_mutation(),
        )
    except Exception as e:
        return dict(
            variant_id=v.variant_id, uniprot_pos=v.uniprot_pos,
            wt_aa=v.wt_aa, mut_aa=v.mut_aa, subtype=v.subtype,
            g_bind_wt=None, g_bind_mut=None, ddg_bind=None,
            direction="unknown", status=f"ERROR:{e}",
            foldx_notation=v.foldx_patient_mutation(),
        )

# ── 主函数 ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="FoldX A1-GPIbα ΔΔG_bind batch runner")
    parser.add_argument("--variants", default="../../data/processed/master_type1_type2.csv")
    parser.add_argument("--pdb", default="../../structures/1M10.pdb")
    parser.add_argument("--foldx-bin",
        default="/inspire/hdd/global_user/mengweicheng-240108120092/lzy/tools/foldx")
    parser.add_argument("--output", default="../../output/foldx_a1_gpiba.csv")
    parser.add_argument("--work-dir", default=None,
        help="FoldX工作目录 (默认: output/foldx_workdir_a1)")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true",
        help="只打印计划, 不运行FoldX (Mac本地测试用)")
    args = parser.parse_args()

    # ── 读取变异 ──────────────────────────────────────────────────────────────
    variants_path = str(Path(args.variants).resolve())
    all_variants = load_variants(variants_path)

    valid, skipped = [], []
    for v in all_variants:
        if v.is_synonymous:
            skipped.append((v.variant_id, "synonymous"))
        elif v.is_stop:
            skipped.append((v.variant_id, "stop_codon"))
        elif not v.in_pdb_range:
            skipped.append((v.variant_id, f"outside_pdb_range(UniProt {v.uniprot_pos}→PDB {v.pdb_pos})"))
        else:
            valid.append(v)

    print(f"A1 domain variants loaded: {len(all_variants)}")
    print(f"  Valid (in PDB range):  {len(valid)}")
    print(f"  Skipped:               {len(skipped)}")
    for vid, reason in skipped[:10]:
        print(f"    SKIP {vid}: {reason}")

    if args.dry_run:
        print("\n=== DRY RUN: FoldX mutation plan ===")
        print(f"{'Variant':<25} {'UniProt':>8} {'PDB':>6} {'Subtype':<10} FoldX notation")
        print("-"*80)

        # WT reference
        wt_ref = Variant("WT_pseudo", 1306, "Q", "R", "WT")  # dummy to get notation
        print(f"{'WT_reference':<25} {'1306':>8} {'543':>6} {'WT':<10} {wt_ref.foldx_wt_mutation()};")

        for v in valid:
            print(f"{v.variant_id:<25} {v.uniprot_pos:>8} {v.pdb_pos:>6} {v.subtype:<10} {v.foldx_patient_mutation()};")
        print(f"\nTotal: 1 WT reference + {len(valid)} patient mutations")
        print(f"Estimated FoldX wall-time on 16 cores: ~{max(1, len(valid)//16)} min (RepairPDB ~10 min first)")
        return

    # ── 真正运行 ──────────────────────────────────────────────────────────────
    pdb_path = str(Path(args.pdb).resolve())
    foldx_bin = args.foldx_bin

    if not os.path.exists(pdb_path):
        sys.exit(f"[ERROR] PDB not found: {pdb_path}")
    if not os.path.exists(foldx_bin):
        sys.exit(f"[ERROR] FoldX binary not found: {foldx_bin}")

    work_dir = args.work_dir or str(Path(args.output).parent / "foldx_workdir_a1")
    os.makedirs(work_dir, exist_ok=True)

    # ── Step 1: RepairPDB (一次) ───────────────────────────────────────────────
    print(f"\n[Step 1] RepairPDB: {Path(pdb_path).name} → may take 5-15 min ...")
    repaired = repair_pdb(foldx_bin, pdb_path, work_dir)
    if not repaired:
        sys.exit("[ERROR] RepairPDB failed. Check repair.log in work_dir.")
    print(f"  Repaired: {repaired}")

    # 验证背景突变位点实际存在
    actual_res = get_residue_at_pos(pdb_path, BG_CHAIN, BG_PDB_POS)
    print(f"  Verified: Chain {BG_CHAIN} position {BG_PDB_POS} = {actual_res} "
          f"(expected {BG_PDB_WT} for R1306Q background)")

    # ── Step 2: WT reference AnalyseComplex ───────────────────────────────────
    print(f"\n[Step 2] Computing WT reference binding energy ...")
    wt_dummy = Variant("WT_ref", 1306, "Q", "R", "WT")  # revert Q543→R only
    wt_energy = build_and_analyse(
        foldx_bin, repaired,
        mutation_str=wt_dummy.foldx_wt_mutation(),
        label="WT_reference",
        work_dir=work_dir,
    )
    if wt_energy is None:
        sys.exit("[ERROR] WT reference failed. Check WT_reference/analyse.log.")
    print(f"  G_bind(WT) = {wt_energy:.3f} kcal/mol")

    # ── Step 3: 并行跑所有患者突变 ────────────────────────────────────────────
    print(f"\n[Step 3] Running {len(valid)} patient variants with {args.workers} workers ...")

    job_args = [(foldx_bin, repaired, v, wt_energy, work_dir) for v in valid]
    results = []

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_variant, a): a[2].variant_id for a in job_args}
        for i, fut in enumerate(as_completed(futures), 1):
            r = fut.result()
            results.append(r)
            status = r['status']
            ddg = f"{r['ddg_bind']:+.2f}" if r['ddg_bind'] is not None else "N/A"
            print(f"  [{i:>3}/{len(valid)}] {r['variant_id']:<25} ΔΔG={ddg:>8}  {r['direction']:<8}  {status}")

    # ── Step 4: 写输出CSV ─────────────────────────────────────────────────────
    output_path = str(Path(args.output).resolve())
    os.makedirs(Path(output_path).parent, exist_ok=True)

    fieldnames = ["variant_id","uniprot_pos","wt_aa","mut_aa","subtype",
                  "g_bind_wt","g_bind_mut","ddg_bind","direction","status","foldx_notation"]
    with open(output_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in sorted(results, key=lambda x: x['uniprot_pos']):
            w.writerow(r)

    ok  = sum(1 for r in results if r['status'] == 'OK')
    gof = sum(1 for r in results if r['direction'] == 'GOF')
    lof = sum(1 for r in results if r['direction'] == 'LOF')
    neu = sum(1 for r in results if r['direction'] == 'neutral')
    print(f"\n=== Done ===")
    print(f"  OK: {ok}/{len(results)}  GOF(2B): {gof}  LOF(2M): {lof}  neutral: {neu}")
    print(f"  Output: {output_path}")


if __name__ == "__main__":
    main()

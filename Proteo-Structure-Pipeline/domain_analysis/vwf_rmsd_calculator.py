#!/usr/bin/env python3
"""
VWF RMSD Calculator with BioPython
使用BioPython计算准确的RMSD

Author: Claude Code
Date: 2026-04-03
"""

import sys
from pathlib import Path
from typing import List, Optional, Tuple, Dict
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from Bio.PDB import Superimposer, NeighborSearch
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO


class VWFRMSDCalculator:
    """
    VWF RMSD计算器
    使用BioPython计算准确的结构偏差
    """

    def __init__(self):
        self.parser = MMCIFParser(QUIET=True)

    def calculate_global_rmsd(self, wt_cif: Path, mut_cif: Path,
                              chain_id: str = 'A') -> Optional[float]:
        """
        计算全局RMSD (所有Cα原子)

        Args:
            wt_cif: WT CIF文件路径
            mut_cif: 突变体CIF文件路径
            chain_id: 链ID (默认'A')

        Returns:
            全局RMSD值(Å)，失败返回None
        """
        try:
            # 解析结构
            wt_structure = self.parser.get_structure('wt', wt_cif)
            mut_structure = self.parser.get_structure('mut', mut_cif)

            # 提取Cα原子
            wt_atoms = self._get_ca_atoms(wt_structure, chain_id)
            mut_atoms = self._get_ca_atoms(mut_structure, chain_id)

            if len(wt_atoms) != len(mut_atoms):
                print(f"Warning: Different number of atoms ({len(wt_atoms)} vs {len(mut_atoms)})")
                # 取最小长度
                min_len = min(len(wt_atoms), len(mut_atoms))
                wt_atoms = wt_atoms[:min_len]
                mut_atoms = mut_atoms[:min_len]

            if len(wt_atoms) == 0:
                return None

            # 叠合并计算RMSD
            superimposer = Superimposer()
            superimposer.set_atoms(wt_atoms, mut_atoms)
            superimposer.apply(mut_structure.get_atoms())

            return superimposer.rms

        except Exception as e:
            print(f"Error calculating global RMSD: {e}")
            return None

    def calculate_local_rmsd(self, wt_cif: Path, mut_cif: Path,
                            center_pos: int,
                            radius: float = 10.0,
                            chain_id: str = 'A') -> Optional[Tuple[float, int]]:
        """
        计算局部RMSD (以指定残基为中心，radius Å半径内)

        Args:
            wt_cif: WT CIF文件路径
            mut_cif: 突变体CIF文件路径
            center_pos: 中心残基位置
            radius: 半径(Å)，默认10.0
            chain_id: 链ID

        Returns:
            (RMSD值, 包含的原子数)，失败返回None
        """
        try:
            # 解析结构
            wt_structure = self.parser.get_structure('wt', wt_cif)
            mut_structure = self.parser.get_structure('mut', mut_cif)

            # 获取中心残基原子
            wt_atoms = self._get_atoms_in_radius(wt_structure, center_pos, radius, chain_id)
            mut_atoms = self._get_atoms_in_radius(mut_structure, center_pos, radius, chain_id)

            if len(wt_atoms) != len(mut_atoms):
                print(f"Warning: Different number of local atoms ({len(wt_atoms)} vs {len(mut_atoms)})")
                min_len = min(len(wt_atoms), len(mut_atoms))
                wt_atoms = wt_atoms[:min_len]
                mut_atoms = mut_atoms[:min_len]

            if len(wt_atoms) == 0:
                return None

            # 叠合并计算RMSD
            superimposer = Superimposer()
            superimposer.set_atoms(wt_atoms, mut_atoms)
            superimposer.apply(mut_structure.get_atoms())

            return superimposer.rms, len(wt_atoms)

        except Exception as e:
            print(f"Error calculating local RMSD: {e}")
            return None

    def calculate_mutation_site_rmsd(self, wt_cif: Path, mut_cif: Path,
                                     mutation_pos: int,
                                     chain_id: str = 'A') -> Optional[Dict]:
        """
        计算突变位点周围的多个半径RMSD

        Args:
            wt_cif: WT CIF文件
            mut_cif: 突变体CIF文件
            mutation_pos: 突变位点
            chain_id: 链ID

        Returns:
            字典包含不同半径的RMSD
        """
        radii = [5.0, 10.0, 15.0, 20.0]
        results = {}

        for radius in radii:
            rmsd_data = self.calculate_local_rmsd(wt_cif, mut_cif, mutation_pos, radius, chain_id)
            if rmsd_data:
                rmsd, num_atoms = rmsd_data
                results[f"rmsd_{int(radius)}a"] = rmsd
                results[f"atoms_{int(radius)}a"] = num_atoms
            else:
                results[f"rmsd_{int(radius)}a"] = None
                results[f"atoms_{int(radius)}a"] = 0

        # 计算全局RMSD
        global_rmsd = self.calculate_global_rmsd(wt_cif, mut_cif, chain_id)
        results["rmsd_global"] = global_rmsd

        return results

    def _get_ca_atoms(self, structure, chain_id: str = 'A') -> List:
        """提取所有Cα原子"""
        ca_atoms = []
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if 'CA' in residue:
                            ca_atoms.append(residue['CA'])
        return ca_atoms

    def _get_atoms_in_radius(self, structure, center_pos: int,
                            radius: float, chain_id: str = 'A') -> List:
        """获取指定残基周围radius Å内的所有原子"""
        # 找到中心残基
        center_residue = None
        center_atom = None

        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if residue.get_id()[1] == center_pos:
                            center_residue = residue
                            if 'CA' in residue:
                                center_atom = residue['CA']
                            break

        if not center_atom:
            print(f"Warning: Position {center_pos} not found in structure")
            return []

        # 收集所有原子
        all_atoms = []
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        for atom in residue:
                            all_atoms.append(atom)

        # 使用NeighborSearch查找radius Å内的原子
        ns = NeighborSearch(all_atoms)
        nearby_atoms = ns.search(center_atom.coord, radius)

        return nearby_atoms

    def batch_calculate_rmsd(self, wt_cif: Path, variant_dirs: List[Path],
                            output_file: Optional[Path] = None) -> pd.DataFrame:
        """
        批量计算RMSD

        Args:
            wt_cif: WT CIF文件
            variant_dirs: 变异体目录列表
            output_file: 输出CSV文件路径

        Returns:
            RMSD结果DataFrame
        """
        import pandas as pd

        records = []

        for var_dir in variant_dirs:
            var_name = var_dir.name

            # 解析变异信息
            import re
            match = re.match(r'fold_vwf_([a-z])(\d+)([a-z])', var_name.lower())
            if not match:
                continue

            ref_aa = match.group(1).upper()
            position = int(match.group(2))
            alt_aa = match.group(3).upper()

            # 查找CIF文件
            cif_file = var_dir / f"{var_name}_model_0.cif"
            if not cif_file.exists():
                continue

            print(f"Processing {var_name}...")

            # 计算RMSD
            rmsd_results = self.calculate_mutation_site_rmsd(
                wt_cif, cif_file, position
            )

            if rmsd_results:
                record = {
                    "variant_id": f"{ref_aa}{position}{alt_aa}",
                    "position": position,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    **rmsd_results
                }
                records.append(record)

        df = pd.DataFrame(records)

        if output_file:
            df.to_csv(output_file, index=False)
            print(f"RMSD results saved to {output_file}")

        return df


def test_rmsd_calculator():
    """测试RMSD计算器"""
    base_dir = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/output/af3_batches_type2/AF3_Results/extracted")

    # 检查可用变异
    available_variants = [d for d in base_dir.iterdir() if d.is_dir() and d.name.startswith('fold_vwf_')]
    print(f"Found {len(available_variants)} variant directories")
    print(f"Examples: {[d.name for d in available_variants[:5]]}")
    print()

    wt_dir = base_dir / "fold_vwf_wt"
    if not wt_dir.exists():
        print(f"WT directory not found: {wt_dir}")
        return

    wt_cif = wt_dir / "fold_vwf_wt_model_0.cif"
    if not wt_cif.exists():
        print(f"WT CIF not found: {wt_cif}")
        return

    calculator = VWFRMSDCalculator()

    # 测试前5个可用变异
    for var_dir in available_variants[:6]:
        if var_dir.name == 'fold_vwf_wt':
            continue

        var_name = var_dir.name
        var_cif = var_dir / f"{var_name}_model_0.cif"

        if not var_cif.exists():
            continue

        # 解析位置
        import re
        match = re.search(r'(\d+)', var_name)
        if not match:
            continue

        pos = int(match.group(1))

        print(f"\n{'='*60}")
        print(f"Calculating RMSD for {var_name}")
        print(f"Position: {pos}")
        print('='*60)

        # 计算多个半径的RMSD
        results = calculator.calculate_mutation_site_rmsd(wt_cif, var_cif, pos)

        if results:
            print(f"Global RMSD: {results.get('rmsd_global', 0):.3f} Å")
            print(f"Local RMSD (5Å):  {results.get('rmsd_5a', 0):.3f} Å ({results.get('atoms_5a', 0)} atoms)")
            print(f"Local RMSD (10Å): {results.get('rmsd_10a', 0):.3f} Å ({results.get('atoms_10a', 0)} atoms)")
            print(f"Local RMSD (15Å): {results.get('rmsd_15a', 0):.3f} Å ({results.get('atoms_15a', 0)} atoms)")
            print(f"Local RMSD (20Å): {results.get('rmsd_20a', 0):.3f} Å ({results.get('atoms_20a', 0)} atoms)")


def main():
    """主函数"""
    print("="*80)
    print("VWF RMSD Calculator - Test with BioPython")
    print("="*80)
    print()

    test_rmsd_calculator()


if __name__ == "__main__":
    main()

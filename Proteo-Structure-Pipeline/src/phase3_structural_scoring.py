#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 3: 3D 空间结构比对与评分引擎
使用 BioPython 进行结构比对、RMSD 计算、pLDDT 分析
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, MMCIF2Dict
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Structure import Structure

# 抑制 BioPython 的 PDBConstructionWarning
warnings.filterwarnings('ignore', category=UserWarning)

# 配置日志
log_dir = Path('../logs')
log_dir.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_dir / 'phase3.log')
    ]
)
logger = logging.getLogger(__name__)


@dataclass
class ScoringResult:
    """评分结果数据类"""
    mutation_name: str
    aa_position: int
    wt_aa: str
    mut_aa: str
    chromosome: str
    genomic_position: str

    # 全局 RMSD
    global_rmsd: Optional[float] = None
    global_aligned_atoms: int = 0

    # 局部 RMSD (10Å 半径)
    local_rmsd_10a: Optional[float] = None
    local_atoms_10a: int = 0

    # pLDDT 分数
    wt_plddt: Optional[float] = None
    mut_plddt: Optional[float] = None
    plddt_delta: Optional[float] = None

    # 状态
    success: bool = False
    error_message: str = ""


class StructureAnalyzer:
    """
    结构分析器
    处理超大 CIF 文件的解析、对齐和评分
    """

    def __init__(self, wt_cif_path: Path):
        self.wt_cif_path = wt_cif_path
        self.parser = MMCIFParser(QUIET=True)

        # 缓存 WT 结构
        logger.info(f"加载 WT 结构: {wt_cif_path}")
        try:
            self.wt_structure = self._safe_parse_cif(wt_cif_path)
            self.wt_atoms = self._get_ca_atoms(self.wt_structure)
            logger.info(f"  ✓ WT 加载成功: {len(self.wt_atoms)} Cα 原子")
        except Exception as e:
            logger.error(f"  ✗ WT 加载失败: {e}")
            raise

    def _safe_parse_cif(self, cif_path: Path) -> Structure:
        """
        安全解析 CIF 文件，处理大型文件和缺失原子
        """
        try:
            structure_id = cif_path.stem
            structure = self.parser.get_structure(structure_id, str(cif_path))
            return structure
        except Exception as e:
            logger.error(f"解析 CIF 失败 {cif_path}: {e}")
            raise

    def _get_ca_atoms(self, structure: Structure) -> Dict[int, object]:
        """
        提取所有 Cα 原子，按残基序号索引
        """
        ca_atoms = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    # 跳过非标准残基
                    if residue.id[0] != ' ':
                        continue

                    res_id = residue.id[1]
                    if 'CA' in residue:
                        ca_atoms[res_id] = residue['CA']
        return ca_atoms

    def _get_all_atoms(self, structure: Structure) -> List:
        """
        获取所有原子列表
        """
        atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != ' ':
                        continue
                    for atom in residue:
                        atoms.append(atom)
        return atoms

    def _get_atoms_in_radius(
        self,
        structure: Structure,
        center_residue_id: int,
        radius: float = 10.0
    ) -> List:
        """
        获取指定残基周围 radius Å 内的所有原子
        """
        # 找到中心 Cα 原子的坐标
        center_ca = None
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[1] == center_residue_id and 'CA' in residue:
                        center_ca = residue['CA']
                        break
                if center_ca:
                    break
            if center_ca:
                break

        if center_ca is None:
            logger.warning(f"  未找到残基 {center_residue_id} 的 Cα 原子")
            return []

        center_coord = center_ca.get_coord()

        # 收集半径内的所有原子
        nearby_atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != ' ':
                        continue
                    for atom in residue:
                        dist = np.linalg.norm(atom.get_coord() - center_coord)
                        if dist <= radius:
                            nearby_atoms.append(atom)

        return nearby_atoms

    def _extract_plddt(self, structure: Structure, residue_id: int) -> Optional[float]:
        """
        从 B-factor 列提取指定残基的 pLDDT 分数
        AlphaFold 将 pLDDT 存储在 B-factor 中
        """
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[1] == residue_id:
                        # 使用 Cα 原子的 B-factor 作为残基的 pLDDT
                        if 'CA' in residue:
                            return residue['CA'].get_bfactor()
        return None

    def calculate_rmsd(
        self,
        ref_atoms: List,
        target_atoms: List
    ) -> Tuple[Optional[float], int]:
        """
        计算 RMSD，处理原子数不匹配的情况
        """
        if len(ref_atoms) != len(target_atoms):
            # 尝试通过名称匹配
            ref_dict = {a.get_name(): a for a in ref_atoms}
            target_dict = {a.get_name(): a for a in target_atoms}

            common_names = set(ref_dict.keys()) & set(target_dict.keys())
            if len(common_names) < 3:
                logger.warning(f"  共同原子太少: {len(common_names)}")
                return None, 0

            ref_atoms = [ref_dict[n] for n in common_names]
            target_atoms = [target_dict[n] for n in common_names]

        if len(ref_atoms) == 0:
            return None, 0

        ref_coords = np.array([a.get_coord() for a in ref_atoms])
        target_coords = np.array([a.get_coord() for a in target_atoms])

        # 计算 RMSD
        diff = ref_coords - target_coords
        rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))

        return float(rmsd), len(ref_atoms)

    def superimpose_and_calculate_rmsd(
        self,
        ref_atoms: Dict[int, object],
        target_atoms: Dict[int, object]
    ) -> Tuple[Optional[float], int]:
        """
        使用 Bio.PDB.Superimposer 进行结构叠合并计算 RMSD
        """
        # 找到共同的残基
        common_residues = set(ref_atoms.keys()) & set(target_atoms.keys())

        if len(common_residues) < 3:
            logger.warning(f"  共同残基太少，无法叠合: {len(common_residues)}")
            return None, 0

        ref_atom_list = [ref_atoms[r] for r in sorted(common_residues)]
        target_atom_list = [target_atoms[r] for r in sorted(common_residues)]

        try:
            sup = Superimposer()
            sup.set_atoms(ref_atom_list, target_atom_list)
            sup.apply(target_atom_list)  # 应用旋转矩阵

            return float(sup.rms), len(common_residues)
        except Exception as e:
            logger.warning(f"  叠合失败: {e}")
            return None, 0

    def analyze_mutation(
        self,
        mut_cif_path: Path,
        mutation_info: Dict
    ) -> ScoringResult:
        """
        分析单个突变结构
        """
        result = ScoringResult(
            mutation_name=mutation_info.get('name', mut_cif_path.stem),
            aa_position=mutation_info.get('aa_position', 0),
            wt_aa=mutation_info.get('wt_aa', ''),
            mut_aa=mutation_info.get('mut_aa', ''),
            chromosome=mutation_info.get('chromosome', ''),
            genomic_position=str(mutation_info.get('genomic_position', ''))
        )

        try:
            logger.info(f"分析: {result.mutation_name}")

            # 加载突变结构
            mut_structure = self._safe_parse_cif(mut_cif_path)
            mut_atoms = self._get_ca_atoms(mut_structure)

            if len(mut_atoms) == 0:
                result.error_message = "突变结构无 Cα 原子"
                return result

            # ===== 1. 全局 RMSD (基于 Cα 叠合) =====
            logger.info(f"  计算全局 RMSD...")
            global_rmsd, aligned_count = self.superimpose_and_calculate_rmsd(
                self.wt_atoms, mut_atoms
            )
            result.global_rmsd = global_rmsd if global_rmsd is not None else np.nan
            result.global_aligned_atoms = aligned_count
            if global_rmsd is not None:
                logger.info(f"    Global RMSD: {global_rmsd:.4f} Å ({aligned_count} atoms)")
            else:
                logger.warning(f"    Global RMSD 计算失败: 共同残基不足 ({aligned_count} atoms)")

            # ===== 2. 局部 RMSD (突变位点周围 10Å) =====
            logger.info(f"  计算局部 RMSD (10Å 半径)...")
            aa_pos = result.aa_position

            # 获取 WT 和 Mutant 在突变位点周围的原子
            wt_local_atoms = self._get_atoms_in_radius(
                self.wt_structure, aa_pos, radius=10.0
            )
            mut_local_atoms = self._get_atoms_in_radius(
                mut_structure, aa_pos, radius=10.0
            )

            if len(wt_local_atoms) > 0 and len(mut_local_atoms) > 0:
                local_rmsd, local_count = self.calculate_rmsd(
                    wt_local_atoms, mut_local_atoms
                )
                result.local_rmsd_10a = local_rmsd if local_rmsd is not None else np.nan
                result.local_atoms_10a = local_count
                if local_rmsd is not None:
                    logger.info(f"    Local RMSD: {local_rmsd:.4f} Å ({local_count} atoms)")
                else:
                    logger.warning(f"    Local RMSD 计算失败: 共同原子不足 ({local_count} atoms)")
            else:
                result.local_rmsd_10a = np.nan
                logger.warning(f"    局部原子不足: WT={len(wt_local_atoms)}, Mut={len(mut_local_atoms)}")

            # ===== 3. pLDDT 分析 =====
            logger.info(f"  提取 pLDDT 分数...")
            wt_plddt = self._extract_plddt(self.wt_structure, aa_pos)
            mut_plddt = self._extract_plddt(mut_structure, aa_pos)

            result.wt_plddt = wt_plddt if wt_plddt is not None else np.nan
            result.mut_plddt = mut_plddt if mut_plddt is not None else np.nan

            if wt_plddt is not None and mut_plddt is not None:
                result.plddt_delta = mut_plddt - wt_plddt
                logger.info(f"    pLDDT: WT={wt_plddt:.2f}, Mut={mut_plddt:.2f}, Δ={result.plddt_delta:+.2f}")
            else:
                result.plddt_delta = np.nan
                logger.warning(f"    pLDDT 提取失败: WT={wt_plddt}, Mut={mut_plddt}")

            result.success = True

        except Exception as e:
            logger.error(f"  分析失败: {e}")
            result.error_message = str(e)

        return result


class Phase3Runner:
    """Phase 3 主运行器"""

    def __init__(
        self,
        wt_cif: str,
        predictions_dir: str,
        csv_path: str,
        output_path: str
    ):
        self.wt_cif = Path(wt_cif)
        self.predictions_dir = Path(predictions_dir)
        self.csv_path = Path(csv_path)
        self.output_path = Path(output_path)
        self.output_path.parent.mkdir(parents=True, exist_ok=True)

        self.analyzer = None
        self.results = []

    def run(self):
        """运行 Phase 3"""
        logger.info("="*60)
        logger.info("Phase 3: 3D 空间结构比对与评分引擎")
        logger.info("="*60)

        # 初始化分析器
        logger.info("[1/3] 初始化结构分析器...")
        self.analyzer = StructureAnalyzer(self.wt_cif)

        # 读取预测结果列表
        logger.info("[2/3] 读取预测结果...")
        pred_df = pd.read_csv(self.csv_path)
        mutants = pred_df[pred_df['type'] == 'Mutant'].to_dict('records')
        logger.info(f"      共 {len(mutants)} 个突变待分析")

        # 分析每个突变
        logger.info("[3/3] 进行结构比对和评分...")
        for idx, mut_info in enumerate(mutants):
            cif_path = Path(mut_info['output_path'])

            if not cif_path.exists():
                logger.warning(f"  [{idx+1}/{len(mutants)}] 文件不存在，跳过: {cif_path}")
                continue

            result = self.analyzer.analyze_mutation(cif_path, mut_info)
            self.results.append(result)

            # 每 10 个输出一次进度
            if (idx + 1) % 10 == 0:
                logger.info(f"  进度: {idx+1}/{len(mutants)}")

        # 保存结果
        self.save_results()

    def save_results(self):
        """保存评分结果"""
        results_dict = [vars(r) for r in self.results]
        results_df = pd.DataFrame(results_dict)

        # 选择并排序列
        cols = [
            'mutation_name', 'aa_position', 'wt_aa', 'mut_aa',
            'chromosome', 'genomic_position',
            'global_rmsd', 'global_aligned_atoms',
            'local_rmsd_10a', 'local_atoms_10a',
            'wt_plddt', 'mut_plddt', 'plddt_delta',
            'success', 'error_message'
        ]
        available_cols = [c for c in cols if c in results_df.columns]
        results_df = results_df[available_cols]

        results_df.to_csv(self.output_path, index=False)

        success_count = results_df['success'].sum()
        fail_count = len(results_df) - success_count

        # 统计有数值的
        valid_global = results_df['global_rmsd'].notna().sum()
        valid_local = results_df['local_rmsd_10a'].notna().sum()
        valid_plddt = results_df['plddt_delta'].notna().sum()

        logger.info("="*60)
        logger.info("Phase 3 完成!")
        logger.info(f"  总分析数: {len(results_df)}")
        logger.info(f"  成功: {success_count}, 失败: {fail_count}")
        logger.info(f"  有效 Global RMSD: {valid_global}")
        logger.info(f"  有效 Local RMSD: {valid_local}")
        logger.info(f"  有效 pLDDT Delta: {valid_plddt}")
        logger.info(f"  输出: {self.output_path}")
        logger.info("="*60)

        # 打印前 5 个最高分（按 local_rmsd）
        if valid_local > 0:
            top_structural = results_df[results_df['local_rmsd_10a'].notna()].nlargest(5, 'local_rmsd_10a')
            logger.info("\nTop 5 结构变化最大 (Local RMSD):")
            for _, row in top_structural.iterrows():
                logger.info(f"  {row['mutation_name']}: Local={row['local_rmsd_10a']:.3f}Å, Global={row['global_rmsd']:.3f}Å, pLDDTΔ={row.get('plddt_delta', 'N/A')}")


def main():
    parser = argparse.ArgumentParser(
        description='Phase 3: 3D 空间结构比对与评分引擎'
    )
    parser.add_argument(
        '--wt-cif',
        type=str,
        default='../structures/predictions/VWF_WT.cif',
        help='WT 结构文件路径 (CIF 格式)'
    )
    parser.add_argument(
        '--predictions-dir',
        type=str,
        default='../structures/predictions',
        help='预测结构目录'
    )
    parser.add_argument(
        '--csv',
        type=str,
        default='../output/phase2_prediction_results.csv',
        help='Phase 2 结果 CSV'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='../output/phase3_final_scores.csv',
        help='输出评分 CSV'
    )

    args = parser.parse_args()

    runner = Phase3Runner(
        wt_cif=args.wt_cif,
        predictions_dir=args.predictions_dir,
        csv_path=args.csv,
        output_path=args.output
    )

    runner.run()


if __name__ == '__main__':
    main()

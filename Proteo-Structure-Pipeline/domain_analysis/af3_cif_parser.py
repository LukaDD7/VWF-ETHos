#!/usr/bin/env python3
"""
AlphaFold3 CIF Parser and Feature Extractor
自动从AF3 CIF文件和JSON文件提取结构特征

功能：
1. 解析CIF文件提取全局和局部RMSD
2. 从JSON文件提取残基级pLDDT
3. 自动匹配WT和突变体结构
4. 计算pLDDT变化和局部RMSD

Author: Claude Code
Date: 2026-04-03
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import re


@dataclass
class AF3StructureData:
    """AF3结构数据容器"""
    variant_id: str
    position: int
    ref_aa: str
    alt_aa: str

    # pLDDT数据
    plddt_per_residue: Dict[int, float]  # 残基位置 -> pLDDT
    mean_plddt: float
    min_plddt: float
    max_plddt: float

    # 从CIF提取的数据
    global_plddt: float  # CIF中的全局pLDDT
    ptm_score: float
    iptm_score: Optional[float]
    ranking_score: float
    has_clash: bool
    fraction_disordered: float

    # 文件路径
    cif_path: Path
    json_path: Path

    def get_plddt_at_position(self, pos: int) -> Optional[float]:
        """获取指定位置的pLDDT"""
        return self.plddt_per_residue.get(pos)


class AF3CIFParser:
    """
    AF3 CIF文件解析器
    """

    def __init__(self, base_dir: Path):
        self.base_dir = Path(base_dir)
        self.wt_dir = self.base_dir / "fold_vwf_wt"  # 假设有WT目录

    def parse_variant_directory(self, variant_dir: Path) -> Optional[AF3StructureData]:
        """
        解析变异体目录，提取所有结构信息

        Args:
            variant_dir: 变异体目录 (如 fold_vwf_r1306w)

        Returns:
            AF3StructureData 或 None
        """
        variant_name = variant_dir.name

        # 解析变异信息
        variant_info = self._parse_variant_name(variant_name)
        if not variant_info:
            return None

        # 找到best model (model_0.cif 和 full_data_0.json)
        cif_file = variant_dir / f"{variant_name}_model_0.cif"
        json_file = variant_dir / f"{variant_name}_full_data_0.json"

        if not cif_file.exists() or not json_file.exists():
            print(f"Warning: Missing files for {variant_name}")
            return None

        # 解析JSON文件 (pLDDT数据)
        plddt_data = self._parse_full_data_json(json_file)
        if not plddt_data:
            return None

        # 解析CIF文件 (全局分数)
        cif_data = self._parse_cif_file(cif_file)

        # 组合数据
        return AF3StructureData(
            variant_id=variant_info["variant_id"],
            position=variant_info["position"],
            ref_aa=variant_info["ref_aa"],
            alt_aa=variant_info["alt_aa"],
            plddt_per_residue=plddt_data["per_residue"],
            mean_plddt=plddt_data["mean"],
            min_plddt=plddt_data["min"],
            max_plddt=plddt_data["max"],
            global_plddt=cif_data.get("global_plddt", 0),
            ptm_score=cif_data.get("ptm", 0),
            iptm_score=cif_data.get("iptm"),
            ranking_score=cif_data.get("ranking_score", 0),
            has_clash=cif_data.get("has_clash", False),
            fraction_disordered=cif_data.get("fraction_disordered", 0),
            cif_path=cif_file,
            json_path=json_file
        )

    def _parse_variant_name(self, name: str) -> Optional[Dict]:
        """
        解析变异体目录名
        格式: fold_vwf_[ref][pos][alt] 如 fold_vwf_r1306w

        Returns:
            {"variant_id": "R1306W", "position": 1306, "ref_aa": "R", "alt_aa": "W"}
        """
        match = re.match(r'fold_vwf_([a-z])(\d+)([a-z])', name.lower())
        if not match:
            # 尝试其他格式
            match = re.match(r'fold_vwf_([a-z]+)(\d+)([a-z]+)', name.lower())
            if not match:
                return None

        ref_aa = match.group(1).upper()
        if len(ref_aa) > 1:
            # 3字母代码转换
            ref_aa = self._three_to_one(ref_aa) or ref_aa[0].upper()

        position = int(match.group(2))

        alt_aa = match.group(3).upper()
        if len(alt_aa) > 1:
            alt_aa = self._three_to_one(alt_aa) or alt_aa[0].upper()

        return {
            "variant_id": f"{ref_aa}{position}{alt_aa}",
            "position": position,
            "ref_aa": ref_aa,
            "alt_aa": alt_aa
        }

    def _three_to_one(self, three_letter: str) -> Optional[str]:
        """3字母氨基酸代码转1字母"""
        aa_map = {
            'ala': 'A', 'arg': 'R', 'asn': 'N', 'asp': 'D', 'cys': 'C',
            'gln': 'Q', 'glu': 'E', 'gly': 'G', 'his': 'H', 'ile': 'I',
            'leu': 'L', 'lys': 'K', 'met': 'M', 'phe': 'F', 'pro': 'P',
            'ser': 'S', 'thr': 'T', 'trp': 'W', 'tyr': 'Y', 'val': 'V'
        }
        return aa_map.get(three_letter.lower())

    def _parse_full_data_json(self, json_file: Path) -> Optional[Dict]:
        """
        解析full_data JSON文件提取pLDDT

        JSON结构:
        {
            "atom_plddts": [list of pLDDT per atom],
            "token_res_ids": [list of residue IDs per token],
            ...
        }
        """
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)

            atom_plddts = data.get('atom_plddts', [])
            token_res_ids = data.get('token_res_ids', [])

            if not atom_plddts or not token_res_ids:
                return None

            # 按残基分组计算平均pLDDT
            # token_res_ids[i] 对应 atom_plddts中该token的原子
            # 假设每个token的第一个atom代表该残基的pLDDT
            residue_plddt = {}

            for i, res_id in enumerate(token_res_ids):
                if res_id not in residue_plddt and i < len(atom_plddts):
                    residue_plddt[res_id] = atom_plddts[i]

            if not residue_plddt:
                return None

            plddt_values = list(residue_plddt.values())

            return {
                "per_residue": residue_plddt,
                "mean": np.mean(plddt_values),
                "min": np.min(plddt_values),
                "max": np.max(plddt_values)
            }

        except Exception as e:
            print(f"Error parsing {json_file}: {e}")
            return None

    def _parse_cif_file(self, cif_file: Path) -> Dict:
        """
        解析CIF文件提取全局分数
        """
        result = {
            "global_plddt": 0,
            "ptm": 0,
            "iptm": None,
            "ranking_score": 0,
            "has_clash": False,
            "fraction_disordered": 0
        }

        try:
            with open(cif_file, 'r') as f:
                content = f.read()

            # 提取全局pLDDT
            plddt_match = re.search(r'_ma_qa_metric_global\.metric_value\s+([\d.]+)', content)
            if plddt_match:
                result["global_plddt"] = float(plddt_match.group(1))

            # 注意：PTM/IPTM等从summary_confidences JSON获取更方便
            summary_json = cif_file.parent / f"{cif_file.stem.replace('_model_0', '')}_summary_confidences_0.json"
            if summary_json.exists():
                with open(summary_json, 'r') as f:
                    conf_data = json.load(f)
                result["ptm"] = conf_data.get('ptm', 0)
                result["iptm"] = conf_data.get('iptm')
                result["ranking_score"] = conf_data.get('ranking_score', 0)
                result["has_clash"] = bool(conf_data.get('has_clash', 0))
                result["fraction_disordered"] = conf_data.get('fraction_disordered', 0)

        except Exception as e:
            print(f"Error parsing CIF {cif_file}: {e}")

        return result

    def extract_structural_features(self, variant_dirs: List[Path],
                                     wt_data: Optional[AF3StructureData] = None) -> pd.DataFrame:
        """
        批量提取结构特征

        Args:
            variant_dirs: 变异体目录列表
            wt_data: WT结构数据 (可选)

        Returns:
            DataFrame with structural features
        """
        records = []

        for var_dir in variant_dirs:
            var_data = self.parse_variant_directory(var_dir)
            if not var_data:
                continue

            # 基础特征
            record = {
                "variant_id": var_data.variant_id,
                "position": var_data.position,
                "ref_aa": var_data.ref_aa,
                "alt_aa": var_data.alt_aa,
                "plddt_mean": var_data.mean_plddt,
                "plddt_min": var_data.min_plddt,
                "plddt_max": var_data.max_plddt,
                "plddt_at_site": var_data.get_plddt_at_position(var_data.position),
                "global_plddt": var_data.global_plddt,
                "ptm": var_data.ptm_score,
                "iptm": var_data.iptm_score,
                "ranking_score": var_data.ranking_score,
                "has_clash": var_data.has_clash,
                "fraction_disordered": var_data.fraction_disordered,
            }

            # 如果有WT数据，计算delta
            if wt_data:
                wt_plddt = wt_data.get_plddt_at_position(var_data.position)
                var_plddt = var_data.get_plddt_at_position(var_data.position)

                record["wt_plddt"] = wt_plddt
                record["plddt_delta"] = var_plddt - wt_plddt if var_plddt and wt_plddt else None

                # 计算全局RMSD (使用pLDDT差异近似)
                # 实际RMSD需要结构叠合，这里先用pLDDT差异
                record["plddt_global_diff"] = var_data.mean_plddt - wt_data.mean_plddt

            records.append(record)

        return pd.DataFrame(records)

    def calculate_local_rmsd(self, wt_cif: Path, mut_cif: Path,
                            center_pos: int, radius: float = 10.0) -> Optional[float]:
        """
        计算局部RMSD (需要BioPython进行结构叠合)

        Args:
            wt_cif: WT CIF文件路径
            mut_cif: 突变体CIF文件路径
            center_pos: 中心残基位置
            radius: 半径 (Å)

        Returns:
            局部RMSD值
        """
        try:
            from Bio.PDB import PDBParser, Superimposer
            from Bio.PDB.Structure import Structure

            # 解析CIF文件
            parser = PDBParser(QUIET=True)

            # 注意：BioPython的PDBParser可能不直接支持CIF
            # 需要MMCIFParser
            from Bio.PDB.MMCIFParser import MMCIFParser
            parser = MMCIFParser(QUIET=True)

            wt_structure = parser.get_structure('wt', wt_cif)
            mut_structure = parser.get_structure('mut', mut_cif)

            # 提取中心残基周围的原子
            wt_atoms = self._get_atoms_in_radius(wt_structure, center_pos, radius)
            mut_atoms = self._get_atoms_in_radius(mut_structure, center_pos, radius)

            if len(wt_atoms) != len(mut_atoms) or len(wt_atoms) == 0:
                return None

            # 叠合并计算RMSD
            superimposer = Superimposer()
            superimposer.set_atoms(wt_atoms, mut_atoms)
            superimposer.apply(mut_structure.get_atoms())

            return superimposer.rms

        except ImportError:
            print("BioPython not available, skipping RMSD calculation")
            return None
        except Exception as e:
            print(f"Error calculating RMSD: {e}")
            return None

    def _get_atoms_in_radius(self, structure, center_pos: int, radius: float) -> List:
        """获取指定残基周围radius Å内的原子"""
        from Bio.PDB import NeighborSearch

        # 获取中心残基的CA原子
        center_atom = None
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[1] == center_pos and 'CA' in residue:
                        center_atom = residue['CA']
                        break

        if not center_atom:
            return []

        # 构建邻居搜索
        all_atoms = [atom for atom in structure.get_atoms() if atom.element != 'H']
        ns = NeighborSearch(all_atoms)

        # 搜索radius Å内的原子
        nearby_atoms = ns.search(center_atom.coord, radius)

        # 只返回CA原子用于叠合
        ca_atoms = [atom for atom in nearby_atoms if atom.name == 'CA']

        return ca_atoms


def main():
    """测试CIF解析器"""
    import glob

    base_dir = Path("/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/Proteo-Structure-Pipeline/output/af3_batches_type2/AF3_Results/extracted")

    parser = AF3CIFParser(base_dir)

    # 查找所有变异体目录
    variant_dirs = [d for d in base_dir.iterdir() if d.is_dir() and d.name.startswith('fold_vwf_')]

    print(f"Found {len(variant_dirs)} variant directories")
    print()

    # 解析前5个
    for var_dir in sorted(variant_dirs)[:5]:
        print(f"\n{'='*80}")
        print(f"Processing: {var_dir.name}")
        print('='*80)

        data = parser.parse_variant_directory(var_dir)
        if data:
            print(f"Variant: {data.variant_id}")
            print(f"Position: {data.position}")
            print(f"Mean pLDDT: {data.mean_plddt:.2f}")
            print(f"Global pLDDT: {data.global_plddt:.2f}")
            print(f"PTM: {data.ptm_score:.3f}")
            print(f"Ranking Score: {data.ranking_score:.3f}")
            print(f"pLDDT at position {data.position}: {data.get_plddt_at_position(data.position):.2f}")

            # 显示域内pLDDT (A1域: 1271-1492)
            if 1271 <= data.position <= 1492:
                a1_plddts = [data.get_plddt_at_position(i) for i in range(1271, 1493) if data.get_plddt_at_position(i)]
                if a1_plddts:
                    print(f"A1 domain mean pLDDT: {np.mean(a1_plddts):.2f}")


if __name__ == "__main__":
    main()

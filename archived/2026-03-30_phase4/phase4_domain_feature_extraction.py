#!/usr/bin/env python3
"""
Phase 4: Domain Feature Extraction & ML Classification for VWF Type-2

功能:
1. 从AlphaFold3 CIF文件中提取RMSD和pLDDT特征 (使用Bio.PDB)
2. 整合多源特征数据
3. 训练Random Forest和XGBoost分类器进行Type-2分型预测
4. 输出预测结果和模型评估

输入:
    - CIF结构文件: Proteo-Structure-Pipeline/structures/predictions/
    - Phase 3输出: variant_features.csv
    - Type-2变异表: VWF_Type2_AF3_Reference_Table.csv
    - AlphaGenome结果: 03_inference_results.csv

输出:
    - results/phase4_domain_features.csv (特征表)
    - results/phase4_ml_predictions.csv (ML预测结果)
    - results/phase4_model_comparison.png (模型对比图)
    - results/phase4_feature_importance.png (特征重要性图)

作者: Claude Code
日期: 2026-03-30
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# 尝试导入Bio.PDB，如果失败则给出警告
try:
    from Bio.PDB import MMCIFParser, Superimposer, PDBIO, Select
    from Bio.PDB.PDBParser import PDBParser
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    warnings.warn("Bio.PDB not available. RMSD calculation will be skipped.")

# 尝试导入sklearn和xgboost
try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    from sklearn.preprocessing import LabelEncoder
    from sklearn.metrics import classification_report, confusion_matrix
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    warnings.warn("scikit-learn not available. Random Forest will be skipped.")

try:
    import xgboost as xgb
    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False
    warnings.warn("XGBoost not available. XGBoost classifier will be skipped.")

try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# ==================== 配置区域 ====================
SCRIPT_DIR = Path(__file__).parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"
DATA_DIR = BASE_DIR / "data"
AF3_OUTPUT_DIR = BASE_DIR / "Proteo-Structure-Pipeline" / "output" / "af3_batches_type2" / "AF3_Results" / "analysis"
STRUCTURES_DIR = BASE_DIR / "Proteo-Structure-Pipeline" / "structures" / "predictions"
WT_STRUCTURE = BASE_DIR / "Proteo-Structure-Pipeline" / "structures" / "predictions" / "VWF_WT.cif"

OUTPUT_FILE = RESULTS_DIR / "phase4_domain_features.csv"
ML_OUTPUT_FILE = RESULTS_DIR / "phase4_ml_predictions.csv"

# VWF功能域定义 (基于UniProt P04275)
DOMAINS = {
    "D1-D2": (1, 272, "Signal peptide + Propeptide"),
    "D'D3": (273, 761, "FVIII binding"),
    "A1": (762, 1035, "Platelet binding (GPIb)"),
    "A2": (1036, 1227, "ADAMTS13 cleavage"),
    "A3": (1228, 1451, "Collagen binding"),
    "D4": (1452, 1670, "Multimerization"),
    "C1-C2": (1671, 2051, "Integrin binding"),
    "CT": (2052, 2813, "C-terminal"),
}

# Type-2亚型与功能域的关联
SUBTYPE_DOMAIN_ASSOCIATIONS = {
    "type2a": ["D4", "A2", "D1-D2"],
    "type2b": ["A1"],
    "type2m": ["A1", "A3"],
    "type2n": ["D'D3"],
}

# ACMG分类编码
ACMG_ENCODING = {
    "Pathogenic": 4,
    "Likely pathogenic": 3,
    "Likely_pathogenic": 3,
    "Uncertain significance": 2,
    "Uncertain_significance": 2,
    "Uncertain signifi-cance": 2,
    "Likely benign": 1,
    "Likely_benign": 1,
    "Benign": 0,
}

# 氨基酸特性
AA_PROPERTIES = {
    "A": {"hydrophobic": 1, "size": "small", "charge": "neutral"},
    "C": {"hydrophobic": 1, "size": "small", "charge": "neutral", "special": "cysteine"},
    "D": {"hydrophobic": 0, "size": "small", "charge": "negative"},
    "E": {"hydrophobic": 0, "size": "medium", "charge": "negative"},
    "F": {"hydrophobic": 1, "size": "large", "charge": "neutral", "aromatic": 1},
    "G": {"hydrophobic": 0, "size": "small", "charge": "neutral"},
    "H": {"hydrophobic": 0, "size": "medium", "charge": "positive", "aromatic": 1},
    "I": {"hydrophobic": 1, "size": "medium", "charge": "neutral"},
    "K": {"hydrophobic": 0, "size": "medium", "charge": "positive"},
    "L": {"hydrophobic": 1, "size": "medium", "charge": "neutral"},
    "M": {"hydrophobic": 1, "size": "medium", "charge": "neutral"},
    "N": {"hydrophobic": 0, "size": "small", "charge": "neutral"},
    "P": {"hydrophobic": 0, "size": "small", "charge": "neutral", "special": "proline"},
    "Q": {"hydrophobic": 0, "size": "medium", "charge": "neutral"},
    "R": {"hydrophobic": 0, "size": "large", "charge": "positive"},
    "S": {"hydrophobic": 0, "size": "small", "charge": "neutral"},
    "T": {"hydrophobic": 0, "size": "small", "charge": "neutral"},
    "V": {"hydrophobic": 1, "size": "small", "charge": "neutral"},
    "W": {"hydrophobic": 1, "size": "large", "charge": "neutral", "aromatic": 1},
    "Y": {"hydrophobic": 1, "size": "large", "charge": "neutral", "aromatic": 1},
    "*": {"hydrophobic": 0, "size": "none", "charge": "none", "special": "stop"},
}

# 日志配置
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# ==================== Bio.PDB 结构分析类 ====================

class StructureAnalyzer:
    """使用Bio.PDB进行结构分析，计算RMSD和pLDDT"""

    def __init__(self, wt_structure_path: Path | None = None):
        self.wt_structure_path = wt_structure_path or WT_STRUCTURE
        self.parser = MMCIFParser(QUIET=True)
        self.wt_structure = None
        self.wt_atoms = None
        self.wt_plddt = None

        if BIOPYTHON_AVAILABLE and self.wt_structure_path.exists():
            self._load_wt_structure()

    def _load_wt_structure(self):
        """加载野生型结构"""
        try:
            self.wt_structure = self.parser.get_structure("WT", self.wt_structure_path)
            self.wt_atoms = self._get_ca_atoms(self.wt_structure)
            self.wt_plddt = self._extract_plddt(self.wt_structure)
            logger.info(f"成功加载野生型结构: {self.wt_structure_path}")
            logger.info(f"WT结构包含 {len(self.wt_atoms)} 个CA原子")
        except Exception as e:
            logger.error(f"加载野生型结构失败: {e}")
            self.wt_structure = None

    def _get_ca_atoms(self, structure) -> list:
        """获取所有CA原子"""
        ca_atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if "CA" in residue:
                        ca_atoms.append(residue["CA"])
        return ca_atoms

    def _extract_plddt(self, structure) -> dict[int, float]:
        """从CIF文件的B-factor中提取pLDDT值"""
        plddt_dict = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_id = residue.get_id()[1]
                    if "CA" in residue:
                        # pLDDT存储在B-factor中
                        plddt_dict[res_id] = residue["CA"].get_bfactor()
        return plddt_dict

    def calculate_rmsd(self, mutant_structure_path: Path) -> dict[str, float | None]:
        """
        计算突变体与野生型的RMSD

        Returns:
            dict with keys: global_rmsd, local_rmsd_10a, pldtt_delta, success
        """
        if not BIOPYTHON_AVAILABLE:
            return {"global_rmsd": None, "local_rmsd_10a": None, "plddt_delta": None, "success": False}

        if self.wt_structure is None:
            logger.warning("野生型结构未加载，无法计算RMSD")
            return {"global_rmsd": None, "local_rmsd_10a": None, "plddt_delta": None, "success": False}

        if not mutant_structure_path.exists():
            logger.warning(f"突变体结构文件不存在: {mutant_structure_path}")
            return {"global_rmsd": None, "local_rmsd_10a": None, "plddt_delta": None, "success": False}

        try:
            # 加载突变体结构
            mut_structure = self.parser.get_structure("MUT", mutant_structure_path)
            mut_atoms = self._get_ca_atoms(mut_structure)
            mut_plddt = self._extract_plddt(mut_structure)

            # 计算全局RMSD
            global_rmsd = self._calculate_global_rmsd(mut_atoms)

            # 计算局部RMSD (突变位点±10个残基)
            mutation_position = self._extract_mutation_position(mutant_structure_path.stem)
            local_rmsd = self._calculate_local_rmsd(mut_atoms, mutation_position)

            # 计算pLDDT差异
            plddt_delta = self._calculate_plddt_delta(mut_plddt, mutation_position)

            return {
                "global_rmsd": global_rmsd,
                "local_rmsd_10a": local_rmsd,
                "plddt_delta": plddt_delta,
                "success": True
            }

        except Exception as e:
            logger.warning(f"计算RMSD失败 {mutant_structure_path}: {e}")
            return {"global_rmsd": None, "local_rmsd_10a": None, "plddt_delta": None, "success": False}

    def _calculate_global_rmsd(self, mut_atoms: list) -> float | None:
        """计算全局RMSD"""
        try:
            # 找到WT和MUT的共同原子
            min_len = min(len(self.wt_atoms), len(mut_atoms))
            if min_len == 0:
                return None

            wt_atoms_subset = self.wt_atoms[:min_len]
            mut_atoms_subset = mut_atoms[:min_len]

            superimposer = Superimposer()
            superimposer.set_atoms(wt_atoms_subset, mut_atoms_subset)
            superimposer.apply(mut_atoms_subset)

            return superimposer.rms
        except Exception as e:
            logger.warning(f"全局RMSD计算失败: {e}")
            return None

    def _calculate_local_rmsd(self, mut_atoms: list, mutation_position: int | None, radius: int = 10) -> float | None:
        """计算局部RMSD (突变位点附近±radius个残基)"""
        if mutation_position is None:
            return None

        try:
            # 获取局部区域的CA原子
            wt_local = []
            mut_local = []

            for wt_atom in self.wt_atoms:
                res_id = wt_atom.get_parent().get_id()[1]
                if abs(res_id - mutation_position) <= radius:
                    wt_local.append(wt_atom)

            for mut_atom in mut_atoms:
                res_id = mut_atom.get_parent().get_id()[1]
                if abs(res_id - mutation_position) <= radius:
                    mut_local.append(mut_atom)

            if len(wt_local) == 0 or len(mut_local) == 0 or len(wt_local) != len(mut_local):
                return None

            superimposer = Superimposer()
            superimposer.set_atoms(wt_local, mut_local)
            superimposer.apply(mut_local)

            return superimposer.rms
        except Exception as e:
            logger.warning(f"局部RMSD计算失败: {e}")
            return None

    def _calculate_plddt_delta(self, mut_plddt: dict, mutation_position: int | None) -> float | None:
        """计算突变位点的pLDDT差异"""
        if mutation_position is None or not self.wt_plddt:
            return None

        wt_score = self.wt_plddt.get(mutation_position)
        mut_score = mut_plddt.get(mutation_position)

        if wt_score is None or mut_score is None:
            return None

        return mut_score - wt_score

    def _extract_mutation_position(self, job_name: str) -> int | None:
        """从job名称中提取突变位置 (e.g., VWF_G1531D -> 1531)"""
        match = re.search(r"[A-Z](\d+)[A-Z*]", job_name)
        if match:
            return int(match.group(1))
        return None


# ==================== ML分类器类 ====================

class Type2Classifier:
    """Type-2分型分类器：Random Forest和XGBoost"""

    def __init__(self):
        self.rf_model = None
        self.xgb_model = None
        self.label_encoder = LabelEncoder()
        self.feature_names = None
        self.is_trained = False

    def prepare_features(self, df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray | None]:
        """准备特征矩阵"""
        # 选择数值特征
        exclude_cols = [
            "variant_id", "aa_change", "wt_aa", "mut_aa",
            "vwf_domain", "domain_description", "type2_subtype",
            "acmg_classification", "size_change"
        ]

        feature_cols = [col for col in df.columns if col not in exclude_cols]
        self.feature_names = feature_cols

        # 构建特征矩阵
        X = df[feature_cols].copy()

        # 将所有列转换为数值，非数值转为NaN
        for col in X.columns:
            X[col] = pd.to_numeric(X[col], errors='coerce')

        # 处理布尔值
        for col in X.columns:
            if X[col].dtype == bool:
                X[col] = X[col].astype(int)

        # 处理缺失值 - 逐列计算中位数
        for col in X.columns:
            if X[col].isna().all():
                X[col] = 0  # 如果全为NaN，填充0
            else:
                X[col] = X[col].fillna(X[col].median())

        # 转换为numpy数组
        X = X.values.astype(float)

        # 准备标签
        y = None
        if "type2_subtype" in df.columns:
            y = self._encode_labels(df["type2_subtype"])

        return X, y

    def _encode_labels(self, subtype_series: pd.Series) -> np.ndarray:
        """编码Type-2亚型标签"""
        # 标准化标签
        labels = subtype_series.str.lower().str.strip()
        labels = labels.replace("wt_control", "type2a")  # WT作为对照，归为type2a
        labels = labels.replace("uncertain", "type2a")   # uncertain暂时归为type2a

        return self.label_encoder.fit_transform(labels)

    def train(self, X: np.ndarray, y: np.ndarray, cv_folds: int = 5) -> dict[str, Any]:
        """
        训练分类器并返回交叉验证结果

        Returns:
            dict with model performance metrics
        """
        results = {}

        # Random Forest
        if SKLEARN_AVAILABLE:
            logger.info("训练Random Forest分类器...")
            self.rf_model = RandomForestClassifier(
                n_estimators=200,
                max_depth=10,
                min_samples_split=3,
                min_samples_leaf=2,
                random_state=42,
                n_jobs=-1
            )

            # 交叉验证
            cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
            rf_scores = cross_val_score(self.rf_model, X, y, cv=cv, scoring="accuracy")

            self.rf_model.fit(X, y)
            results["random_forest"] = {
                "cv_accuracy_mean": rf_scores.mean(),
                "cv_accuracy_std": rf_scores.std(),
                "feature_importance": dict(zip(self.feature_names, self.rf_model.feature_importances_))
            }
            logger.info(f"Random Forest CV Accuracy: {rf_scores.mean():.3f} ± {rf_scores.std():.3f}")

        # XGBoost
        if XGBOOST_AVAILABLE:
            logger.info("训练XGBoost分类器...")
            self.xgb_model = xgb.XGBClassifier(
                n_estimators=200,
                max_depth=6,
                learning_rate=0.1,
                subsample=0.8,
                colsample_bytree=0.8,
                random_state=42,
                n_jobs=-1
            )

            # 交叉验证
            cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
            xgb_scores = cross_val_score(self.xgb_model, X, y, cv=cv, scoring="accuracy")

            self.xgb_model.fit(X, y)
            results["xgboost"] = {
                "cv_accuracy_mean": xgb_scores.mean(),
                "cv_accuracy_std": xgb_scores.std(),
                "feature_importance": dict(zip(self.feature_names, self.xgb_model.feature_importances_))
            }
            logger.info(f"XGBoost CV Accuracy: {xgb_scores.mean():.3f} ± {xgb_scores.std():.3f}")

        self.is_trained = True
        return results

    def predict(self, X: np.ndarray) -> dict[str, np.ndarray]:
        """使用训练好的模型进行预测"""
        predictions = {}

        if self.rf_model is not None:
            predictions["rf_pred"] = self.label_encoder.inverse_transform(self.rf_model.predict(X))
            predictions["rf_proba"] = self.rf_model.predict_proba(X).max(axis=1)

        if self.xgb_model is not None:
            predictions["xgb_pred"] = self.label_encoder.inverse_transform(self.xgb_model.predict(X))
            predictions["xgb_proba"] = self.xgb_model.predict_proba(X).max(axis=1)

        return predictions

    def plot_feature_importance(self, output_path: Path):
        """绘制特征重要性图"""
        if not MATPLOTLIB_AVAILABLE:
            logger.warning("matplotlib不可用，跳过绘图")
            return

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Random Forest
        if self.rf_model is not None:
            importances = self.rf_model.feature_importances_
            indices = np.argsort(importances)[-10:]  # Top 10
            axes[0].barh(range(len(indices)), importances[indices])
            axes[0].set_yticks(range(len(indices)))
            axes[0].set_yticklabels([self.feature_names[i] for i in indices])
            axes[0].set_xlabel("Feature Importance")
            axes[0].set_title("Random Forest - Top 10 Features")

        # XGBoost
        if self.xgb_model is not None:
            importances = self.xgb_model.feature_importances_
            indices = np.argsort(importances)[-10:]  # Top 10
            axes[1].barh(range(len(indices)), importances[indices])
            axes[1].set_yticks(range(len(indices)))
            axes[1].set_yticklabels([self.feature_names[i] for i in indices])
            axes[1].set_xlabel("Feature Importance")
            axes[1].set_title("XGBoost - Top 10 Features")

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info(f"特征重要性图已保存: {output_path}")
        plt.close()


# ==================== 原有功能函数 ====================

@dataclass
class DomainFeatures:
    """功能域级别的特征集合"""

    variant_id: str
    aa_change: str
    position: int
    wt_aa: str
    mut_aa: str
    vwf_domain: str
    domain_description: str
    type2_subtype: str
    acmg_classification: str
    acmg_score: int
    mut_pae_self: float
    wt_pae_self: float
    pae_delta_self: float
    mut_local_flex: float
    wt_local_flex: float
    pae_delta_local: float
    domain_avg_pae_delta: float
    domain_max_pae_delta: float
    global_rmsd: float | None = None
    local_rmsd_10a: float | None = None
    plddt_delta: float | None = None
    alphagenome_max_score: float | None = None
    splice_delta: float | None = None
    rna_delta: float | None = None
    domain_relative_position: float = 0.0
    is_domain_hotspot: bool = False
    cysteine_disruption: bool = False
    charge_change: int = 0
    size_change: str = "none"
    subtype_domain_match: bool = False
    mechanism_scores: dict[str, float] = field(default_factory=dict)


def get_domain_for_position(position: int) -> tuple[str, str]:
    for domain_name, (start, end, description) in DOMAINS.items():
        if start <= position <= end:
            return domain_name, description
    return "Unknown", "Unknown"


def calculate_domain_relative_position(position: int, domain: str) -> float:
    if domain not in DOMAINS:
        return 0.0
    start, end, _ = DOMAINS[domain]
    return (position - start) / (end - start) if end > start else 0.0


def is_domain_hotspot(position: int, domain: str) -> bool:
    if domain == "A1":
        return position in range(1260, 1280) or position in range(1370, 1400)
    elif domain == "A2":
        return position in range(1480, 1510)
    elif domain == "A3":
        return position in range(1580, 1620)
    elif domain == "D4":
        return position in range(1680, 1700)
    return False


def parse_aa_change(aa_change: str) -> tuple[str, str]:
    if pd.isna(aa_change) or aa_change == "-":
        return "", ""

    aa_change = str(aa_change).strip()

    match = re.search(r"(?:p\.)?([A-Z])(\d+)([A-Z*])", aa_change)
    if match:
        return match.group(1), match.group(3)

    match = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", aa_change)
    if match:
        three_to_one = {
            "Ala": "A", "Cys": "C", "Asp": "D", "Glu": "E",
            "Phe": "F", "Gly": "G", "His": "H", "Ile": "I",
            "Lys": "K", "Leu": "L", "Met": "M", "Asn": "N",
            "Pro": "P", "Gln": "Q", "Arg": "R", "Ser": "S",
            "Thr": "T", "Val": "V", "Trp": "W", "Tyr": "Y",
        }
        wt = three_to_one.get(match.group(1))
        mut = three_to_one.get(match.group(3))
        if wt and mut:
            return wt, mut

    match = re.search(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", aa_change)
    if match:
        three_to_one = {
            "Val": "V", "Phe": "F", "Gly": "G", "Asp": "D",
            "Ala": "A", "Cys": "C", "Glu": "E", "His": "H",
            "Ile": "I", "Lys": "K", "Leu": "L", "Met": "M",
            "Asn": "N", "Pro": "P", "Gln": "Q", "Arg": "R",
            "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y",
        }
        wt = three_to_one.get(match.group(1))
        mut = three_to_one.get(match.group(3))
        if wt and mut:
            return wt, mut

    return "", ""


def calculate_aa_property_changes(wt_aa: str, mut_aa: str) -> dict[str, Any]:
    wt_props = AA_PROPERTIES.get(wt_aa, {})
    mut_props = AA_PROPERTIES.get(mut_aa, {})

    charge_map = {"negative": -1, "neutral": 0, "positive": 1, "none": 0}
    wt_charge = charge_map.get(wt_props.get("charge", "neutral"), 0)
    mut_charge = charge_map.get(mut_props.get("charge", "neutral"), 0)
    charge_change = mut_charge - wt_charge

    size_order = {"small": 1, "medium": 2, "large": 3, "none": 0}
    wt_size = size_order.get(wt_props.get("size", "small"), 1)
    mut_size = size_order.get(mut_props.get("size", "small"), 1)

    if mut_size > wt_size:
        size_change = "larger"
    elif mut_size < wt_size:
        size_change = "smaller"
    else:
        size_change = "same"

    cysteine_disruption = wt_aa == "C" or mut_aa == "C"
    hydrophobic_change = mut_props.get("hydrophobic", 0) - wt_props.get("hydrophobic", 0)

    return {
        "charge_change": charge_change,
        "size_change": size_change,
        "cysteine_disruption": cysteine_disruption,
        "hydrophobic_change": hydrophobic_change,
    }


def check_subtype_domain_match(subtype: str, domain: str) -> bool:
    if pd.isna(subtype) or subtype == "WT_Control":
        return False

    subtype = str(subtype).lower().replace(" ", "")
    associated_domains = SUBTYPE_DOMAIN_ASSOCIATIONS.get(subtype, [])
    return domain in associated_domains


def encode_acmg(acmg_classification: str) -> int:
    if pd.isna(acmg_classification):
        return 2

    acmg_str = str(acmg_classification).strip()
    return ACMG_ENCODING.get(acmg_str, 2)


def load_variant_features() -> pd.DataFrame | None:
    variant_file = AF3_OUTPUT_DIR / "variant_features.csv"
    if not variant_file.exists():
        logger.error(f"找不到variant_features.csv: {variant_file}")
        return None

    logger.info(f"加载variant_features.csv: {variant_file}")
    df = pd.read_csv(variant_file)
    logger.info(f"加载了 {len(df)} 条变异记录")
    return df


def load_type2_reference_table() -> pd.DataFrame | None:
    ref_file = BASE_DIR / "VWF_Type2_AF3_Reference_Table.csv"
    if not ref_file.exists():
        logger.error(f"找不到Type-2参考表: {ref_file}")
        return None

    logger.info(f"加载Type-2参考表: {ref_file}")
    df = pd.read_csv(ref_file)
    logger.info(f"加载了 {len(df)} 条Type-2变异记录")
    return df


def load_alphagenome_results() -> pd.DataFrame | None:
    ag_file = RESULTS_DIR / "03_inference_results.csv"
    if not ag_file.exists():
        logger.warning(f"找不到AlphaGenome结果: {ag_file}")
        return None

    logger.info(f"加载AlphaGenome结果: {ag_file}")
    df = pd.read_csv(ag_file)
    logger.info(f"加载了 {len(df)} 条AlphaGenome记录")
    return df


def find_cif_file(job_name: str) -> Path | None:
    """查找CIF文件"""
    # 尝试不同命名格式
    patterns = [
        STRUCTURES_DIR / f"{job_name}.cif",
        STRUCTURES_DIR / f"fold_{job_name.lower()}" / "model.cif",
        STRUCTURES_DIR / f"fold_{job_name.lower()}.cif",
    ]

    for pattern in patterns:
        if pattern.exists():
            return pattern

    return None


def extract_structural_features(variant_df: pd.DataFrame, analyzer: StructureAnalyzer) -> pd.DataFrame:
    """从CIF文件中提取结构特征 (RMSD, pLDDT)"""
    if not BIOPYTHON_AVAILABLE:
        logger.warning("Bio.PDB不可用，跳过结构特征提取")
        return variant_df

    logger.info("开始从CIF文件提取结构特征...")

    rmsd_list = []
    local_rmsd_list = []
    plddt_delta_list = []

    for idx, row in variant_df.iterrows():
        job_name = row.get("Job_Name", f"VWF_{row.get('AA_Change', 'unknown')}")

        # 查找CIF文件
        cif_file = find_cif_file(job_name)

        if cif_file and cif_file.exists():
            logger.debug(f"分析结构: {cif_file}")
            rmsd_results = analyzer.calculate_rmsd(cif_file)

            rmsd_list.append(rmsd_results.get("global_rmsd"))
            local_rmsd_list.append(rmsd_results.get("local_rmsd_10a"))
            plddt_delta_list.append(rmsd_results.get("plddt_delta"))
        else:
            rmsd_list.append(None)
            local_rmsd_list.append(None)
            plddt_delta_list.append(None)

    variant_df["global_rmsd"] = rmsd_list
    variant_df["local_rmsd_10a"] = local_rmsd_list
    variant_df["plddt_delta"] = plddt_delta_list

    # 统计成功计算的数量
    success_count = sum(1 for r in rmsd_list if r is not None)
    logger.info(f"成功计算 {success_count}/{len(variant_df)} 个结构的RMSD")

    return variant_df


def merge_data_sources(
    variant_df: pd.DataFrame,
    type2_df: pd.DataFrame | None,
    alphagenome_df: pd.DataFrame | None,
) -> pd.DataFrame:
    logger.info("开始合并数据源...")

    merged = variant_df.copy()

    if type2_df is not None:
        type2_subset = type2_df[["Position", "AA_Change", "Type2_Subtype", "ACMG_Classification",
                                 "AlphaGenome_Max_Score", "Splice_Delta", "RSID", "Job_Name"]].copy()
        type2_subset["Position"] = pd.to_numeric(type2_subset["Position"], errors="coerce")

        merged = merged.merge(
            type2_subset,
            on="Position",
            how="left",
            suffixes=("", "_ref"),
        )

    if alphagenome_df is not None:
        pos_col = None
        for col in alphagenome_df.columns:
            if "pos" in col.lower():
                pos_col = col
                break

        if pos_col:
            ag_subset = alphagenome_df.copy()
            ag_subset["Position"] = pd.to_numeric(ag_subset[pos_col], errors="coerce")

            score_cols = [col for col in ag_subset.columns if "delta" in col.lower() or "score" in col.lower()]
            if score_cols:
                ag_subset = ag_subset[["Position"] + score_cols[:3]]

            merged = merged.merge(
                ag_subset,
                on="Position",
                how="left",
            )

    logger.info(f"合并后共 {len(merged)} 条记录，{len(merged.columns)} 列")
    return merged


def extract_features(df: pd.DataFrame) -> list[DomainFeatures]:
    logger.info("开始提取功能域特征...")

    features_list = []

    for idx, row in df.iterrows():
        try:
            position = int(row.get("Position", 0))
            aa_change = row.get("AA_Change", "")

            if position == 0:
                continue

            wt_aa, mut_aa = parse_aa_change(aa_change)
            domain, domain_desc = get_domain_for_position(position)
            domain_rel_pos = calculate_domain_relative_position(position, domain)
            is_hotspot = is_domain_hotspot(position, domain)
            aa_changes = calculate_aa_property_changes(wt_aa, mut_aa)
            acmg = row.get("ACMG_Classification", "Uncertain_significance")
            acmg_score = encode_acmg(acmg)
            subtype = row.get("Type2_Subtype", "uncertain")
            subtype_match = check_subtype_domain_match(subtype, domain)

            feature = DomainFeatures(
                variant_id=f"VWF_{aa_change}",
                aa_change=aa_change,
                position=position,
                wt_aa=wt_aa,
                mut_aa=mut_aa,
                vwf_domain=domain,
                domain_description=domain_desc,
                type2_subtype=subtype,
                acmg_classification=acmg,
                acmg_score=acmg_score,
                mut_pae_self=row.get("Mut_PAE_Self", 0.0),
                wt_pae_self=row.get("WT_PAE_Self", 0.0),
                pae_delta_self=row.get("PAE_Delta_Self", 0.0),
                mut_local_flex=row.get("Mut_Local_Flex", 0.0),
                wt_local_flex=row.get("WT_Local_Flex", 0.0),
                pae_delta_local=row.get("PAE_Delta_Local", 0.0),
                domain_avg_pae_delta=row.get("Domain_Avg_PAE_Delta", 0.0),
                domain_max_pae_delta=row.get("Domain_Max_PAE_Delta", 0.0),
                global_rmsd=row.get("global_rmsd"),
                local_rmsd_10a=row.get("local_rmsd_10a"),
                plddt_delta=row.get("plddt_delta"),
                alphagenome_max_score=row.get("AlphaGenome_Max_Score"),
                splice_delta=row.get("Splice_Delta"),
                domain_relative_position=domain_rel_pos,
                is_domain_hotspot=is_hotspot,
                cysteine_disruption=aa_changes["cysteine_disruption"],
                charge_change=aa_changes["charge_change"],
                size_change=aa_changes["size_change"],
                subtype_domain_match=subtype_match,
            )

            features_list.append(feature)

        except Exception as e:
            logger.warning(f"处理第 {idx} 行时出错: {e}")
            continue

    logger.info(f"成功提取 {len(features_list)} 条特征记录")
    return features_list


def features_to_dataframe(features_list: list[DomainFeatures]) -> pd.DataFrame:
    records = []

    for feat in features_list:
        record = {
            "variant_id": feat.variant_id,
            "aa_change": feat.aa_change,
            "position": feat.position,
            "wt_aa": feat.wt_aa,
            "mut_aa": feat.mut_aa,
            "vwf_domain": feat.vwf_domain,
            "domain_description": feat.domain_description,
            "type2_subtype": feat.type2_subtype,
            "acmg_classification": feat.acmg_classification,
            "acmg_score": feat.acmg_score,
            "mut_pae_self": feat.mut_pae_self,
            "wt_pae_self": feat.wt_pae_self,
            "pae_delta_self": feat.pae_delta_self,
            "mut_local_flex": feat.mut_local_flex,
            "wt_local_flex": feat.wt_local_flex,
            "pae_delta_local": feat.pae_delta_local,
            "domain_avg_pae_delta": feat.domain_avg_pae_delta,
            "domain_max_pae_delta": feat.domain_max_pae_delta,
            "global_rmsd": feat.global_rmsd,
            "local_rmsd_10a": feat.local_rmsd_10a,
            "plddt_delta": feat.plddt_delta,
            "alphagenome_max_score": feat.alphagenome_max_score,
            "splice_delta": feat.splice_delta,
            "domain_relative_position": feat.domain_relative_position,
            "is_domain_hotspot": feat.is_domain_hotspot,
            "cysteine_disruption": feat.cysteine_disruption,
            "charge_change": feat.charge_change,
            "size_change": feat.size_change,
            "subtype_domain_match": feat.subtype_domain_match,
        }
        records.append(record)

    return pd.DataFrame(records)


def calculate_summary_statistics(df: pd.DataFrame) -> dict[str, Any]:
    stats = {
        "total_variants": len(df),
        "unique_domains": df["vwf_domain"].nunique(),
        "subtype_distribution": df["type2_subtype"].value_counts().to_dict(),
        "acmg_distribution": df["acmg_classification"].value_counts().to_dict(),
        "hotspot_variants": df["is_domain_hotspot"].sum(),
        "cysteine_disruptions": df["cysteine_disruption"].sum(),
        "subtype_domain_matches": df["subtype_domain_match"].sum(),
    }

    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    key_numeric_cols = [
        "pae_delta_self", "pae_delta_local", "domain_avg_pae_delta",
        "domain_max_pae_delta", "alphagenome_max_score", "acmg_score",
        "global_rmsd", "local_rmsd_10a", "plddt_delta"
    ]

    for col in key_numeric_cols:
        if col in df.columns and col in numeric_cols:
            stats[f"{col}_mean"] = round(df[col].mean(), 4)
            stats[f"{col}_std"] = round(df[col].std(), 4)
            stats[f"{col}_nonnull"] = df[col].notna().sum()

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Phase 4: Domain Feature Extraction & ML Classification for VWF Type-2"
    )
    parser.add_argument("--output", type=str, default=str(OUTPUT_FILE), help=f"特征输出路径 (默认: {OUTPUT_FILE})")
    parser.add_argument("--ml-output", type=str, default=str(ML_OUTPUT_FILE), help=f"ML预测输出路径")
    parser.add_argument("--no-merge", action="store_true", help="仅使用variant_features.csv")
    parser.add_argument("--no-structural", action="store_true", help="跳过CIF结构分析")
    parser.add_argument("--no-ml", action="store_true", help="跳过ML分类")
    parser.add_argument("--summary", action="store_true", help="输出摘要统计")
    args = parser.parse_args()

    logger.info("=" * 70)
    logger.info("Phase 4: Domain Feature Extraction & ML Classification")
    logger.info("=" * 70)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # 1. 加载数据
    variant_df = load_variant_features()
    if variant_df is None:
        logger.error("无法加载variant_features.csv，退出")
        sys.exit(1)

    # 2. 从CIF文件提取结构特征 (RMSD, pLDDT)
    if not args.no_structural and BIOPYTHON_AVAILABLE:
        analyzer = StructureAnalyzer()
        variant_df = extract_structural_features(variant_df, analyzer)
    else:
        variant_df["global_rmsd"] = None
        variant_df["local_rmsd_10a"] = None
        variant_df["plddt_delta"] = None

    # 3. 合并其他数据源
    type2_df = None if args.no_merge else load_type2_reference_table()
    alphagenome_df = None if args.no_merge else load_alphagenome_results()
    merged_df = merge_data_sources(variant_df, type2_df, alphagenome_df)

    # 4. 提取所有特征
    features_list = extract_features(merged_df)

    if not features_list:
        logger.error("没有提取到任何特征，退出")
        sys.exit(1)

    # 5. 保存特征表
    output_df = features_to_dataframe(features_list)
    output_path = Path(args.output)
    output_df.to_csv(output_path, index=False)
    logger.info(f"\n特征表已保存: {output_path}")
    logger.info(f"  行数: {len(output_df)}, 列数: {len(output_df.columns)}")

    # 6. ML分类
    if not args.no_ml and SKLEARN_AVAILABLE:
        logger.info("\n" + "=" * 70)
        logger.info("开始ML分类...")
        logger.info("=" * 70)

        classifier = Type2Classifier()
        X, y = classifier.prepare_features(output_df)

        if y is not None and len(np.unique(y)) >= 2:
            # 训练模型
            results = classifier.train(X, y, cv_folds=5)

            # 预测
            predictions = classifier.predict(X)

            # 保存预测结果
            ml_df = output_df.copy()
            if "rf_pred" in predictions:
                ml_df["rf_prediction"] = predictions["rf_pred"]
                ml_df["rf_confidence"] = predictions["rf_proba"]
            if "xgb_pred" in predictions:
                ml_df["xgb_prediction"] = predictions["xgb_pred"]
                ml_df["xgb_confidence"] = predictions["xgb_proba"]

            ml_output_path = Path(args.ml_output)
            ml_df.to_csv(ml_output_path, index=False)
            logger.info(f"\nML预测结果已保存: {ml_output_path}")

            # 绘制特征重要性
            if MATPLOTLIB_AVAILABLE:
                importance_plot = RESULTS_DIR / "phase4_feature_importance.png"
                classifier.plot_feature_importance(importance_plot)
        else:
            logger.warning("标签不足，跳过ML分类")

    # 7. 输出摘要
    if args.summary:
        stats = calculate_summary_statistics(output_df)
        logger.info("\n" + "=" * 70)
        logger.info("特征摘要统计")
        logger.info("=" * 70)
        for key, value in stats.items():
            logger.info(f"  {key}: {value}")

    logger.info("\n" + "=" * 70)
    logger.info("Phase 4 完成!")
    logger.info("=" * 70)


if __name__ == "__main__":
    main()

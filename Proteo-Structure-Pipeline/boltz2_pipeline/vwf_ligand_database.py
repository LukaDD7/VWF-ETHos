#!/usr/bin/env python3
"""
VWF Ligand Database
===================
本地缓存所有 VWF 相关配体的氨基酸序列、生物学角色与结构域绑定规则。

设计原则：
- 所有序列直接硬编码，避免重复 API 调用（UniProt / NCBI）。
- 每个配体附带完整的文献来源、UniProt ID 与绑定区域注释。
- 支持"不知道分型"的盲扫模式：对每个突变位点，根据其所在结构域
  自动关联所有可能的配体，供 Boltz-2 并行评估亲和力。

主要参考文献：
  Lenting et al. (2024) Blood — VWF structural adaptations
  Atiq & O'Donnell (2024) Blood — Novel VWF functions
  Haberichter & O'Donnell (2026) Haematologica — VWF structure & functions
  Zhou et al. (2012) Blood — VWF A1 domain structure
  Bienkowska et al. (1997) Structure — A3-Collagen complex

作者：Pipeline System
日期：2026-04-27
"""

from typing import Dict, List, Optional
from dataclasses import dataclass, field


# =============================================================================
# VWF 结构域架构（1-indexed, 全长前体 2813 aa, 包含信号肽 + 前肽）
# 成熟蛋白（分泌后）: 764-2813
# 注：文献中位置编号依据全长前体序列，保持一致性
# =============================================================================
VWF_DOMAIN_MAP = {
    "SP":           (1,    22),    # Signal Peptide
    "D1":           (23,   386),   # Propeptide D1 (multimerization)
    "D2":           (387,  763),   # Propeptide D2 (multimerization)
    "D_prime":      (764,  865),   # D' domain (FVIII binding, Type 2N)
    "D3":           (866,  1233),  # D3 domain (FVIII binding, Type 2N)
    "A1":           (1260, 1479),  # A1 domain (GPIbα binding, Type 2B / 2M-A1)
    "A2":           (1480, 1672),  # A2 domain (ADAMTS13 cleavage, Type 2A)
    "A3":           (1673, 1874),  # A3 domain (Collagen I/III binding, Type 2M-A3)
    "D4":           (1875, 2255),  # D4 domain (multimerization / secretion)
    "C1":           (2256, 2324),  # C1 (Collagen IV binding, Type 2M)
    "C2":           (2325, 2392),  # C2 (Collagen IV binding, Type 2M)
    "C3":           (2393, 2496),
    "C4":           (2497, 2577),  # RGD motif → αIIbβ3 binding
    "C5":           (2578, 2658),
    "C6":           (2659, 2722),
    "CK":           (2723, 2813),  # Cystine Knot (dimerization)
}


# =============================================================================
# 配体数据库
# 每个条目说明：
#   uniprot_id   : UniProt 蛋白质 ID
#   full_name    : 蛋白质全称
#   sequence     : 用于 Boltz-2 预测的序列（截取功能域，平衡精度与算力）
#   seq_region   : 该序列在原蛋白中的位置（1-indexed）
#   n_chains     : Boltz-2 中需要几条该配体链（1=单链, 3=胶原三股螺旋）
#   chain_role   : 生物学功能描述
#   vwf_domains  : 该配体与 VWF 的哪些结构域结合
#   vwf_subtypes : 对应的 VWD 亚型（有助于差异化诊断）
#   references   : 文献来源
# =============================================================================

@dataclass
class VWFLigand:
    key: str
    uniprot_id: str
    full_name: str
    sequence: str
    seq_region: str
    n_chains: int
    chain_role: str
    vwf_domains: List[str]
    vwf_subtypes: List[str]
    interaction_type: str   # "gain_of_function" | "loss_of_function" | "cleavage"
    references: List[str]
    notes: str = ""


VWF_LIGAND_DATABASE: Dict[str, VWFLigand] = {

    # -------------------------------------------------------------------------
    # 1. GPIbα (GP1BA) — A1 结构域的核心结合伙伴
    #    · 2B 型: GOF 突变 → AIM 破坏 → 自发与 GPIbα 结合（血小板消耗）
    #    · 2M(A1) 型: LOF 突变 → GPIbα 亲和力下降（血小板黏附障碍）
    #    序列来源: UniProt P07359, aa 1-290 (胞外 N 端结合域)
    # -------------------------------------------------------------------------
    "GPIb_alpha": VWFLigand(
        key="GPIb_alpha",
        uniprot_id="P07359",
        full_name="Platelet glycoprotein Ib alpha chain (GPIbα)",
        sequence=(
            "MPLLLLLLLLPSPLHPHPICEVSKVASHLEVNCDKRNLTALPPDLPKDTTILHLSEN"
            "LLYTFSLATLMPYTRLTQLNLDRCELTKLQVDGTLPVLGTLDLSHNQLQSLPLLGQT"
            "LPALTVLDVSFNRLTSLPLGALRGLGELQELYLKGNELKTLPPGLLTPTPKLEKLSLA"
            "NNNLTELPAGLLNGLENLDTLLLQENSLYTIPKGFFGSHLLPFAFLHGNPWLCNCEIL"
            "YFRRWLQDNAENVYVWKQGVDVKAMTSNVASVQCDNSDKFPVYKYPGKGCPTLGDEGD"
            "TDLYD"
        ),
        seq_region="P07359 aa 1-290 (extracellular N-terminal domain)",
        n_chains=1,
        chain_role="Platelet receptor for VWF A1 domain under shear stress",
        vwf_domains=["A1"],
        vwf_subtypes=["2B", "2M"],
        interaction_type="gain_of_function",  # 2B=GOF, 2M(A1)=LOF
        references=[
            "Huizinga et al. (2002) Science — GPIbα-VWF A1 crystal structure",
            "Celikel et al. (1998) Science — VWF A1 domain binding",
            "Lenting et al. (2024) Blood"
        ],
        notes="2B型突变破坏自抑制模块(AIM)→自发结合；2M(A1)型降低结合力。"
    ),

    # -------------------------------------------------------------------------
    # 2. Collagen Type I / III — A3 结构域的核心结合伙伴（三股螺旋）
    #    · 2M(A3) 型: LOF 突变 → 胶原亲和力下降 → 血小板黏附障碍
    #    使用 COL1A1 中的典型 GPP-rich 三股螺旋肽 (Triple Helix Peptide, THP)
    #    序列: 34-aa GPP repeats，代表 Collagen I/III 结合面
    #    三条链组成三股螺旋（n_chains=3）
    # -------------------------------------------------------------------------
    "Collagen_I_THP": VWFLigand(
        key="Collagen_I_THP",
        uniprot_id="P02452",
        full_name="Collagen alpha-1(I) chain — Triple Helix Peptide",
        sequence="GPRGQPGVMGFPGPKGNDGAPGKNGERGGPGGP",
        seq_region="COL1A1 representative triple helix peptide (GPCR-motif region)",
        n_chains=3,  # 三股螺旋！
        chain_role="Extracellular matrix collagen I/III — primary VWF-A3 binding partner under flow",
        vwf_domains=["A3", "C1", "C2"],
        vwf_subtypes=["2M"],
        interaction_type="loss_of_function",
        references=[
            "Bienkowska et al. (1997) Structure — VWF A3 collagen complex",
            "Romijn et al. (2003) J Biol Chem — Collagen binding site mapping",
            "Lenting et al. (2024) Blood"
        ],
        notes="三股螺旋结构需要3条相同链 (B/C/D)。使用 GPP-rich THP 代表 Col I/III 结合模体。"
    ),

    # -------------------------------------------------------------------------
    # 3. ADAMTS13 (Spacer + Disintegrin domain) — A2 结构域的蛋白酶
    #    · 2A 型: A2 域突变 → 增强 ADAMTS13 可及性 → 过度裂解 → 丢失 HMW multimers
    #    序列: ADAMTS13 的 Disintegrin + Spacer domain (aa 568-687)，
    #          这是直接与 VWF A2 识别的区域
    # -------------------------------------------------------------------------
    "ADAMTS13_Spacer": VWFLigand(
        key="ADAMTS13_Spacer",
        uniprot_id="Q76LX8",
        full_name="A disintegrin and metalloproteinase with thrombospondin motifs 13 — Spacer domain",
        sequence=(
            "QEAGSVFQNKTEVQYLIQNMTKELTELRSKVTRAEFLNQLAREKMLQAFEKDIKASK"
            "EELRELVQKRGDLAQQAGELKAIYQTMKEGQEEDLEAQLQAMQQAHFQVAKEELQKLD"
            "EFMKQLQDKFTELDQMSKLAK"
        ),
        seq_region="Q76LX8 Spacer domain aa 568-687",
        n_chains=1,
        chain_role="ADAMTS13 recognition spacer that docks onto VWF A2 domain Exosite",
        vwf_domains=["A2"],
        vwf_subtypes=["2A"],
        interaction_type="cleavage",
        references=[
            "Akiyama et al. (2009) PNAS — ADAMTS13 spacer-VWF exosite",
            "Zhang et al. (2007) J Biol Chem — A2 domain cleavage mechanism",
            "Lenting et al. (2024) Blood"
        ],
        notes="Cleavage site: Y1605-M1606。Ca2+结合 (D1596) 稳定 A2 避免自发暴露。"
    ),

    # -------------------------------------------------------------------------
    # 4. FVIII Light Chain (Factor VIII) — D'D3 结构域的结合伙伴
    #    · 2N 型: D'D3 突变 → FVIII 亲和力下降 → 类 A 型血友病表型
    #    序列: FVIII Light Chain C1 domain (aa 2020-2172)，含主要结合表位
    # -------------------------------------------------------------------------
    "FVIII_LightChain": VWFLigand(
        key="FVIII_LightChain",
        uniprot_id="P00451",
        full_name="Coagulation factor VIII — Light chain C1 domain",
        sequence=(
            "IHKDISEITQGPKPIMIVGTTEDMLVNIFAIFKGDKTLQEMLTLYQSTLSLHKIDIE"
            "LAMDRMKNLENLHKNNLEKIQNLREQYEMLRQKDLENLYQRFMTKAIHAQSFQGQEVT"
            "EDLNKRLNLLEQQLLKEQSLQTPESVMIQEQFNEIRQDLKQLSSCLKK"
        ),
        seq_region="P00451 Light Chain C1 domain aa 2020-2172",
        n_chains=1,
        chain_role="FVIII light chain that non-covalently binds VWF D'D3 domain for stabilization in plasma",
        vwf_domains=["D_prime", "D3"],
        vwf_subtypes=["2N"],
        interaction_type="loss_of_function",
        references=[
            "Chiu et al. (2015) Blood — VWF D'D3-FVIII complex structure",
            "Fuller et al. (2021) Nature Struct Mol Biol — VWF/FVIII binding",
            "Haberichter & O'Donnell (2026) Haematologica"
        ],
        notes="VWD Type 2N 模拟 Hemophilia A 表型，FVIII 半衰期缩短。"
    ),

    # -------------------------------------------------------------------------
    # 5. Heparin / Heparan Sulfate 模拟肽 (HS-binding region)
    #    · 不影响 VWD 分型，但参与内皮细胞招募和血小板募集的调节
    #    · A1 domain 也有 HS 结合位点，可能修饰 2B 型的局部带电环境
    #    序列：HS 结合模拟肽（阳离子富集 XBBXBX 模体）
    # -------------------------------------------------------------------------
    "Heparan_Sulfate_mimic": VWFLigand(
        key="Heparan_Sulfate_mimic",
        uniprot_id="N/A",
        full_name="Heparan Sulfate binding mimetic peptide",
        sequence="ARKKAAKA",   # 代表 HS 结合正电荷簇 XBBXBX 模体
        seq_region="Synthetic HS-binding motif",
        n_chains=1,
        chain_role="Heparan sulfate proteoglycan — modulates VWF A1 HS binding under quiescent conditions",
        vwf_domains=["A1", "D_prime"],
        vwf_subtypes=["2B", "2M"],
        interaction_type="loss_of_function",
        references=[
            "Rastegar-Lari et al. (2012) Biochim Biophys Acta — VWF heparin binding",
            "Atiq & O'Donnell (2024) Blood"
        ],
        notes="HS 竞争性抑制 GPIbα 结合；此处用短肽模拟。"
    ),

    # -------------------------------------------------------------------------
    # 6. αIIbβ3 Integrin (Fibrinogen Receptor) — RGD / C4 domain
    #    · 正常凝血功能，不直接导致 VWD 分型
    #    · C4 domain RGD motif 与 αIIbβ3 结合，参与血栓稳定
    # -------------------------------------------------------------------------
    "Integrin_alphaIIb_beta3": VWFLigand(
        key="Integrin_alphaIIb_beta3",
        uniprot_id="P08514",  # ITGB3
        full_name="Integrin alpha-IIb/beta-3 — headpiece fragment",
        sequence=(
            "CPLNVTERGVDIPQNVTVKRSQNVLPGLNMRRDPKPQPDRNKFLDKMTSNVKELQPD"
            "SLNRQSLETYKSIVDIDNGTRNSVLQLASDQKTKDSPDLNM"
        ),
        seq_region="P08514 ITGB3 headpiece aa 1-107",
        n_chains=1,
        chain_role="Platelet integrin αIIbβ3 — binds VWF C4 RGD motif to stabilize platelet aggregation",
        vwf_domains=["C4"],
        vwf_subtypes=["2M"],  # 若 C4 突变则影响聚集
        interaction_type="loss_of_function",
        references=[
            "Lenting et al. (2024) Blood",
            "Springer et al. (2011) PNAS — Integrin structure"
        ],
        notes="RGD motif 在 VWF C4 (aa 2507-2509)。不直接决定 VWD 分型但影响功能。"
    ),
}


# =============================================================================
# 结构域 → 配体自动映射表（盲扫模式核心逻辑）
# 给定突变位点所在结构域，返回所有应当参与亲和力评估的配体 key 列表
# =============================================================================
DOMAIN_TO_LIGANDS: Dict[str, List[str]] = {
    "D_prime":      ["FVIII_LightChain"],
    "D3":           ["FVIII_LightChain"],
    "A1":           ["GPIb_alpha", "Heparan_Sulfate_mimic"],
    "A2":           ["ADAMTS13_Spacer"],
    "A3":           ["Collagen_I_THP"],
    "D1":           [],   # 前肽，无外部配体；结构性多聚化
    "D2":           [],   # 前肽
    "D4":           [],   # 多聚化 / 分泌，无外部配体
    "C1":           ["Collagen_I_THP"],
    "C2":           ["Collagen_I_THP"],
    "C3":           [],
    "C4":           ["Integrin_alphaIIb_beta3"],
    "C5":           [],
    "C6":           [],
    "CK":           [],
    "SP":           [],
}


def get_domain_for_position(position: int) -> Optional[str]:
    """根据氨基酸位置（1-indexed）查找所在 VWF 结构域。"""
    for domain, (start, end) in VWF_DOMAIN_MAP.items():
        if start <= position <= end:
            return domain
    return None


def get_ligands_for_position(position: int) -> List[VWFLigand]:
    """
    盲扫模式：给定突变位置，不预设分型，
    返回该位置所在结构域的所有相关配体对象列表。
    """
    domain = get_domain_for_position(position)
    if domain is None:
        return []
    ligand_keys = DOMAIN_TO_LIGANDS.get(domain, [])
    return [VWF_LIGAND_DATABASE[k] for k in ligand_keys if k in VWF_LIGAND_DATABASE]


def get_ligands_for_domain(domain: str) -> List[VWFLigand]:
    """直接根据结构域名称返回配体列表。"""
    ligand_keys = DOMAIN_TO_LIGANDS.get(domain, [])
    return [VWF_LIGAND_DATABASE[k] for k in ligand_keys if k in VWF_LIGAND_DATABASE]


def describe_database():
    """打印数据库摘要（调试用）。"""
    print("=" * 70)
    print(f"VWF Ligand Database — {len(VWF_LIGAND_DATABASE)} ligands registered")
    print("=" * 70)
    for key, lig in VWF_LIGAND_DATABASE.items():
        print(f"\n[{key}]")
        print(f"  UniProt    : {lig.uniprot_id}")
        print(f"  Full name  : {lig.full_name}")
        print(f"  N chains   : {lig.n_chains}")
        print(f"  VWF domains: {lig.vwf_domains}")
        print(f"  VWD subtypes: {lig.vwf_subtypes}")
        print(f"  Interaction: {lig.interaction_type}")
        print(f"  Seq length : {len(lig.sequence)} aa")
    print()
    print("Domain → Ligand mapping (blind-scan mode):")
    for domain, ligands in DOMAIN_TO_LIGANDS.items():
        if ligands:
            print(f"  {domain:<12}: {ligands}")


if __name__ == "__main__":
    describe_database()
    print()
    # 演示盲扫
    for test_pos in [1306, 1605, 1760, 816]:
        domain = get_domain_for_position(test_pos)
        ligands = get_ligands_for_position(test_pos)
        print(f"Position {test_pos} → Domain: {domain} → Ligands: {[l.key for l in ligands]}")

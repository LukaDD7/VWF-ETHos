#!/usr/bin/env python3
"""
Generate a mechanism-aware Boltz-2 YAML panel for VWF/VWD variants.

The important design choice is diagnostic rather than label-driven:
every variant is projected onto the same VWF functional assay panel, and
only assays whose construct actually contains the mutated residue are sent
to Boltz-2. The full variant x assay matrix is still written so downstream
agents can see which axes were considered and why they were or were not run.

Inputs:
  - original_patient_table/*.xlsx
  - VWF WT FASTA, using full pre-pro VWF P04275 numbering

Outputs:
  - variants_master.csv
  - diagnostic_panel.csv
  - job_manifest.csv
  - yamls/*.yaml
  - optional JSON batches compatible with the existing runner shape
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
BOLTZ2_PIPELINE = REPO_ROOT / "Proteo-Structure-Pipeline" / "boltz2_pipeline"
sys.path.insert(0, str(BOLTZ2_PIPELINE))

from vwf_ligand_database import VWF_DOMAIN_MAP, VWF_LIGAND_DATABASE  # noqa: E402


AA3_TO_AA1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
}

AA_CHANGE_RE = re.compile(
    r"(?P<ref>Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)"
    r"(?P<pos>\d+)"
    r"(?P<alt>Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter)"
)

DOMAIN_ALIASES = {
    "D'": "D_prime",
    "D PRIME": "D_prime",
    "D_PRIME": "D_prime",
    "D-PRIME": "D_prime",
    "D1": "D1",
    "D2": "D2",
    "D3": "D3",
    "D4": "D4",
    "A1": "A1",
    "A2": "A2",
    "A3": "A3",
    "C1": "C1",
    "C2": "C2",
    "C3": "C3",
    "C4": "C4",
    "C5": "C5",
    "C6": "C6",
    "CK": "CK",
    "SP": "SP",
}


@dataclass(frozen=True)
class AssaySpec:
    key: str
    clinical_axis: str
    construct_role: str
    vwf_range: tuple[int, int]
    ligand_keys: tuple[str, ...] = ()
    include_for_domains: tuple[str, ...] = ()
    expected_signal: str = ""
    interpretability: str = "medium"
    notes: str = ""


ASSAYS: tuple[AssaySpec, ...] = (
    AssaySpec(
        key="dprime_d3_fviii_binding",
        clinical_axis="FVIII carrier function / Type 2N axis",
        construct_role="D_prime-D3 plus FVIII light-chain fragment",
        vwf_range=(764, 1233),
        ligand_keys=("FVIII_LightChain",),
        include_for_domains=("D_prime", "D3"),
        expected_signal="loss_of_binding_supports_2N",
        interpretability="high",
    ),
    AssaySpec(
        key="a1_gpiba_forced_binding",
        clinical_axis="platelet GPIb alpha binding / Type 2B and 2M-A1 axis",
        construct_role="A1 domain plus GPIb alpha extracellular fragment",
        vwf_range=(1260, 1479),
        ligand_keys=("GPIb_alpha",),
        include_for_domains=("A1",),
        expected_signal="stronger_binding_can_support_2B; weaker_binding_can_support_2M_A1",
        interpretability="medium",
        notes="Forced complex; it does not test low-shear A1 opening.",
    ),
    AssaySpec(
        key="a1_heparan_sulfate_binding",
        clinical_axis="A1 heparan-sulfate / charge-modulation axis",
        construct_role="A1 domain plus heparan-sulfate mimetic peptide",
        vwf_range=(1260, 1479),
        ligand_keys=("Heparan_Sulfate_mimic",),
        include_for_domains=("A1",),
        expected_signal="modulatory_evidence_for_A1_surface_charge",
        interpretability="exploratory",
    ),
    AssaySpec(
        key="a1_aim_autoinhibition_context",
        clinical_axis="A1 autoinhibition and conformational exposure / Type 2B vs 2M-A1 axis",
        construct_role="AIM-flanked A1 monomeric construct",
        vwf_range=(1234, 1493),
        include_for_domains=("D3", "A1_N_flank", "A1", "A2"),
        expected_signal="changed_AIM_A1_contacts_or_exposed_GPIb_surface",
        interpretability="exploratory",
        notes=(
            "This is a static proxy for AIM/open-closed behavior; it cannot simulate shear "
            "or conformational ensemble populations by itself."
        ),
    ),
    AssaySpec(
        key="a2_folded_stability",
        clinical_axis="A2 fold stability / Type 2A axis",
        construct_role="folded A2 monomer",
        vwf_range=(1480, 1672),
        include_for_domains=("A2",),
        expected_signal="destabilized_A2_can_support_2A_multimer_loss",
        interpretability="medium",
    ),
    AssaySpec(
        key="a2_adamts13_folded_complex",
        clinical_axis="A2 ADAMTS13 docking / Type 2A axis",
        construct_role="folded A2 plus ADAMTS13 spacer-domain fragment",
        vwf_range=(1480, 1672),
        ligand_keys=("ADAMTS13_Spacer",),
        include_for_domains=("A2",),
        expected_signal="altered_ADAMTS13_docking_is_supportive_not_decisive",
        interpretability="exploratory",
        notes="ADAMTS13 cleavage is exposure-dependent; folded docking alone is not the full biology.",
    ),
    AssaySpec(
        key="vwf73_adamts13_substrate",
        clinical_axis="unfolded A2 cleavage-substrate / Type 2A axis",
        construct_role="VWF73 substrate peptide plus ADAMTS13 spacer-domain fragment",
        vwf_range=(1596, 1668),
        ligand_keys=("ADAMTS13_Spacer",),
        include_for_domains=("A2",),
        expected_signal="altered_substrate_recognition_near_Y1605_M1606",
        interpretability="medium",
        notes="Peptide-level proxy for the exposed cleavage substrate, not folded A2 mechanics.",
    ),
    AssaySpec(
        key="a3_collagen_binding",
        clinical_axis="collagen I/III binding / Type 2M-A3 axis",
        construct_role="A3 domain plus collagen triple-helix peptide",
        vwf_range=(1673, 1874),
        ligand_keys=("Collagen_I_THP",),
        include_for_domains=("A3",),
        expected_signal="loss_of_collagen_binding_supports_2M_A3",
        interpretability="high",
    ),
    AssaySpec(
        key="c1_collagen_binding",
        clinical_axis="collagen binding accessory C-domain axis",
        construct_role="C1 domain plus collagen triple-helix peptide",
        vwf_range=(2256, 2324),
        ligand_keys=("Collagen_I_THP",),
        include_for_domains=("C1",),
        expected_signal="loss_of_collagen_binding_supports_qualitative_platelet_adhesion_defect",
        interpretability="exploratory",
    ),
    AssaySpec(
        key="c2_collagen_binding",
        clinical_axis="collagen binding accessory C-domain axis",
        construct_role="C2 domain plus collagen triple-helix peptide",
        vwf_range=(2325, 2392),
        ligand_keys=("Collagen_I_THP",),
        include_for_domains=("C2",),
        expected_signal="loss_of_collagen_binding_supports_qualitative_platelet_adhesion_defect",
        interpretability="exploratory",
    ),
    AssaySpec(
        key="c4_integrin_binding",
        clinical_axis="platelet integrin alphaIIb beta3 / C4 RGD axis",
        construct_role="C4 domain plus integrin beta3 headpiece proxy",
        vwf_range=(2497, 2577),
        ligand_keys=("Integrin_alphaIIb_beta3",),
        include_for_domains=("C4",),
        expected_signal="altered_RGD_axis_supports_platelet_aggregation_defect",
        interpretability="exploratory",
    ),
    AssaySpec(
        key="d1d2_propeptide_context",
        clinical_axis="biosynthesis, propeptide processing, multimerization / Type 1 and 3 axis",
        construct_role="D1-D2 propeptide monomeric context",
        vwf_range=(23, 763),
        include_for_domains=("D1", "D2"),
        expected_signal="local_destabilization_can_support_secretion_or_multimerization_defect",
        interpretability="medium",
    ),
    AssaySpec(
        key="d4_assembly_context",
        clinical_axis="secretion and multimer assembly / quantitative VWD axis",
        construct_role="D4 monomeric context",
        vwf_range=(1875, 2255),
        include_for_domains=("D4",),
        expected_signal="local_destabilization_can_support_quantitative_defect",
        interpretability="medium",
    ),
    AssaySpec(
        key="c_domain_assembly_context",
        clinical_axis="C-domain assembly, dimer packing, secretion / quantitative VWD axis",
        construct_role="C1-C6 monomeric context",
        vwf_range=(2256, 2722),
        include_for_domains=("C1", "C2", "C3", "C4", "C5", "C6"),
        expected_signal="local_destabilization_can_support_quantitative_or_multimer_defect",
        interpretability="medium",
    ),
    AssaySpec(
        key="ck_dimerization_context",
        clinical_axis="C-terminal dimerization / Type 1 and 3 axis",
        construct_role="CK domain monomeric context",
        vwf_range=(2723, 2813),
        include_for_domains=("CK",),
        expected_signal="local_destabilization_can_support_dimerization_or_secretion_defect",
        interpretability="medium",
    ),
)


def read_fasta(path: Path) -> str:
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line and not line.startswith(">"):
                seq.append(line)
    return "".join(seq)


def parse_aa_change(value: object) -> tuple[str, int, str, str] | None:
    if pd.isna(value):
        return None
    text = str(value).replace("\n", "").strip()
    match = AA_CHANGE_RE.search(text)
    if not match:
        return None
    ref3 = match.group("ref")
    alt3 = match.group("alt")
    pos = int(match.group("pos"))
    return AA3_TO_AA1[ref3], pos, AA3_TO_AA1[alt3], f"{ref3}{pos}{alt3}"


def infer_source_label(path: Path) -> str:
    name = path.name
    if name.startswith("1"):
        return "Type1"
    if name.startswith("2A"):
        return "Type2A"
    if name.startswith("2B"):
        return "Type2B"
    if name.startswith("2M"):
        return "Type2M"
    if name.startswith("2N"):
        return "Type2N"
    if name.startswith("3"):
        return "Type3"
    if "GeneBe" in name:
        return "Negative_GeneBe"
    if "Clinvar" in name or "ClinVar" in name:
        return "Negative_ClinVar"
    return "Unknown"


def normalize_domain(value: object) -> str:
    if pd.isna(value):
        return ""
    token = str(value).strip().replace("\n", "").upper()
    return DOMAIN_ALIASES.get(token, "")


def infer_domain_for_position(position: int) -> str:
    for domain, (start, end) in VWF_DOMAIN_MAP.items():
        if start <= position <= end:
            return domain
    if 1234 <= position <= 1259:
        return "A1_N_flank"
    return ""


def extract_domain_from_row(values: Iterable[object]) -> str:
    for value in values:
        domain = normalize_domain(value)
        if domain:
            return domain
    return ""


def load_original_tables(input_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for path in sorted(input_dir.glob("*.xlsx")):
        if path.name == "VWD突变汇总.xlsx":
            continue
        source_label = infer_source_label(path)
        raw = pd.read_excel(path, sheet_name=0, header=None)
        mutation_col = find_mutation_column(raw)
        if mutation_col is None:
            continue
        data_rows = primary_mutation_rows(raw, mutation_col)
        for row_idx, row in raw.iterrows():
            if row_idx not in data_rows:
                continue
            parsed = parse_aa_change(row.get(mutation_col, ""))
            if parsed is None:
                continue

            ref, pos, alt, aa_change_3letter = parsed
            source_domain = extract_domain_from_row(row.tolist())
            inferred_domain = infer_domain_for_position(pos)
            variant_id = f"VWF_{ref}{pos}{alt}"
            rows.append(
                {
                    "variant_id": variant_id,
                    "aa_change_3letter": aa_change_3letter,
                    "wt_aa": ref,
                    "position": pos,
                    "mut_aa": alt,
                    "source_label": source_label,
                    "source_file": path.name,
                    "source_row_1based": row_idx + 1,
                    "source_change_col_0based": mutation_col,
                    "source_domain": source_domain,
                    "inferred_domain": inferred_domain,
                }
            )

    if not rows:
        return pd.DataFrame()

    raw_df = pd.DataFrame(rows)
    grouped = (
        raw_df.groupby(["variant_id", "wt_aa", "position", "mut_aa"], as_index=False)
        .agg(
            aa_change_3letter=("aa_change_3letter", "first"),
            source_labels=("source_label", lambda x: "|".join(sorted(set(map(str, x))))),
            source_files=("source_file", lambda x: "|".join(sorted(set(map(str, x))))),
            source_rows=("source_row_1based", lambda x: "|".join(map(str, sorted(set(map(int, x)))))),
            source_domains=("source_domain", lambda x: "|".join(sorted(set(v for v in map(str, x) if v)))),
            inferred_domain=("inferred_domain", "first"),
        )
        .sort_values(["position", "variant_id"])
        .reset_index(drop=True)
    )
    return grouped


AA1_TO_AA3 = {one: three for three, one in AA3_TO_AA1.items() if one != "*"}


def load_variants_csv(path: Path, subtype_filter: str | None = None) -> pd.DataFrame:
    """Load variants from a clean CSV (cols: aa_change/wt_aa/position/mut_aa[,subtype]).

    Produces the same grouped schema as load_original_tables() so the rest of the
    panel generator is unchanged. Used to feed the expanded label set
    (output/expanded_label_set.csv) instead of the raw patient xlsx tables.
    """
    df = pd.read_csv(path)
    if subtype_filter is not None and "subtype" in df.columns:
        df = df[df["subtype"].astype(str) == subtype_filter]
    rows: list[dict[str, object]] = []
    for _, r in df.iterrows():
        wt, pos, mut = str(r["wt_aa"]), int(r["position"]), str(r["mut_aa"])
        if wt not in AA1_TO_AA3 or mut not in AA1_TO_AA3:
            continue
        aa3 = f"{AA1_TO_AA3[wt]}{pos}{AA1_TO_AA3[mut]}"
        rows.append({
            "variant_id": f"VWF_{wt}{pos}{mut}",
            "aa_change_3letter": aa3,
            "wt_aa": wt, "position": pos, "mut_aa": mut,
            "source_labels": str(r.get("subtype", "")),
            "source_files": Path(path).name,
            "source_rows": "", "source_domains": str(r.get("domain", "") or ""),
            "inferred_domain": infer_domain_for_position(pos),
        })
    if not rows:
        return pd.DataFrame()
    return (pd.DataFrame(rows)
            .groupby(["variant_id", "wt_aa", "position", "mut_aa"], as_index=False)
            .agg(aa_change_3letter=("aa_change_3letter", "first"),
                 source_labels=("source_labels", "first"),
                 source_files=("source_files", "first"),
                 source_rows=("source_rows", "first"),
                 source_domains=("source_domains", "first"),
                 inferred_domain=("inferred_domain", "first"))
            .sort_values(["position", "variant_id"]).reset_index(drop=True))


def find_mutation_column(raw: pd.DataFrame) -> int | None:
    counts: dict[int, int] = {}
    for col_idx in raw.columns:
        counts[int(col_idx)] = int(raw[col_idx].map(lambda value: parse_aa_change(value) is not None).sum())
    best_col, best_count = max(counts.items(), key=lambda item: item[1])
    if best_count == 0:
        return None
    return best_col


def primary_mutation_rows(raw: pd.DataFrame, mutation_col: int) -> set[int]:
    hit_rows = [
        int(row_idx)
        for row_idx, value in raw[mutation_col].items()
        if parse_aa_change(value) is not None
    ]
    if not hit_rows:
        return set()

    blocks: list[list[int]] = [[hit_rows[0]]]
    for row_idx in hit_rows[1:]:
        if row_idx - blocks[-1][-1] <= 2:
            blocks[-1].append(row_idx)
        else:
            blocks.append([row_idx])

    # The source workbooks put the curated variant list as the first contiguous
    # mutation block. Later blocks are references or explanatory notes.
    return set(blocks[0])


def apply_mutation(wt_seq: str, ref: str, pos: int, alt: str) -> tuple[str | None, str]:
    if alt == "*":
        return None, "stop_codon_not_supported_by_structure_prediction"
    if pos < 1 or pos > len(wt_seq):
        return None, "position_out_of_wt_sequence"
    actual = wt_seq[pos - 1]
    if actual != ref:
        return None, f"wt_mismatch_expected_{ref}_found_{actual}"
    return wt_seq[: pos - 1] + alt + wt_seq[pos:], "ok"


def slice_range(seq: str, aa_range: tuple[int, int]) -> str:
    start, end = aa_range
    return seq[start - 1 : end]


def build_sequences(vwf_seq: str, assay: AssaySpec) -> list[dict[str, str]]:
    sequences = [{"id": "A", "sequence": slice_range(vwf_seq, assay.vwf_range)}]
    next_chain = ord("B")
    for ligand_key in assay.ligand_keys:
        ligand = VWF_LIGAND_DATABASE[ligand_key]
        for _ in range(ligand.n_chains):
            sequences.append({"id": chr(next_chain), "sequence": ligand.sequence})
            next_chain += 1
    return sequences


def yaml_for_job(job_name: str, sequences: list[dict[str, str]], comments: dict[str, object]) -> str:
    lines = ["version: 1", "sequences:"]
    for seq in sequences:
        lines.extend(
            [
                "  - protein:",
                f"      id: {seq['id']}",
                f"      sequence: {seq['sequence']}",
                "      msa: empty",
            ]
        )
    lines.append(f"# job: {job_name}")
    for key, value in comments.items():
        safe_value = str(value).replace("\n", " ")
        lines.append(f"# {key}: {safe_value}")
    lines.append("")
    return "\n".join(lines)


def json_job(job_name: str, sequences: list[dict[str, str]]) -> dict[str, object]:
    return {
        "name": job_name,
        "sequences": [
            {"protein": {"id": item["id"], "sequence": item["sequence"]}}
            for item in sequences
        ],
    }


def variant_overlaps_assay(position: int, assay: AssaySpec) -> bool:
    start, end = assay.vwf_range
    return start <= position <= end


def assay_relevant_by_domain(domain: str, assay: AssaySpec) -> bool:
    return domain in assay.include_for_domains


def generate_panel(
    variants: pd.DataFrame,
    wt_seq: str,
    output_dir: Path,
    write_json_batches: bool,
    batch_size: int,
    limit_variants: int | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    yaml_dir = output_dir / "yamls"
    yaml_dir.mkdir(parents=True, exist_ok=True)

    if limit_variants is not None:
        variants = variants.head(limit_variants).copy()

    variant_rows: list[dict[str, object]] = []
    panel_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    json_jobs: list[dict[str, object]] = []
    seen_jobs: set[str] = set()

    # WT baselines are necessary for per-assay delta metrics.
    for assay in ASSAYS:
        job_name = f"VWF_WT__{assay.key}"
        sequences = build_sequences(wt_seq, assay)
        yaml_path = yaml_dir / f"{job_name}.yaml"
        yaml_path.write_text(
            yaml_for_job(
                job_name,
                sequences,
                {
                    "variant_id": "VWF_WT",
                    "assay_key": assay.key,
                    "clinical_axis": assay.clinical_axis,
                    "construct_role": assay.construct_role,
                    "vwf_range": f"{assay.vwf_range[0]}-{assay.vwf_range[1]}",
                    "interpretability": assay.interpretability,
                    "baseline": "WT reference for delta metrics",
                },
            )
        )
        seen_jobs.add(job_name)
        manifest_rows.append(
            {
                "job_name": job_name,
                "yaml_path": str(yaml_path.relative_to(output_dir)),
                "variant_id": "VWF_WT",
                "aa_change": "",
                "position": "",
                "wt_aa": "",
                "mut_aa": "",
                "source_labels": "WT",
                "inferred_domain": "",
                "assay_key": assay.key,
                "clinical_axis": assay.clinical_axis,
                "construct_role": assay.construct_role,
                "vwf_range": f"{assay.vwf_range[0]}-{assay.vwf_range[1]}",
                "ligand_keys": "|".join(assay.ligand_keys),
                "n_chains": len(sequences),
                "run_decision": "WT_BASELINE",
                "expected_signal": assay.expected_signal,
                "interpretability": assay.interpretability,
                "notes": assay.notes,
            }
        )
        if write_json_batches:
            json_jobs.append(json_job(job_name, sequences))

    for _, row in variants.iterrows():
        variant_id = str(row["variant_id"])
        ref = str(row["wt_aa"])
        pos = int(row["position"])
        alt = str(row["mut_aa"])
        inferred_domain = str(row.get("inferred_domain", "")) or infer_domain_for_position(pos)
        mutated_seq, mutation_status = apply_mutation(wt_seq, ref, pos, alt)

        variant_rows.append(
            {
                **row.to_dict(),
                "wt_validation_status": mutation_status,
                "domain_source_agrees_with_inferred": (
                    "" if not row.get("source_domains") else inferred_domain in str(row.get("source_domains")).split("|")
                ),
            }
        )

        for assay in ASSAYS:
            overlaps = variant_overlaps_assay(pos, assay)
            domain_relevant = assay_relevant_by_domain(inferred_domain, assay)
            if mutation_status != "ok":
                run_decision = "SKIP_MUTATION_NOT_MODELABLE"
            elif overlaps:
                run_decision = "RUN"
            elif domain_relevant:
                run_decision = "SKIP_RELEVANT_DOMAIN_OUTSIDE_CONSTRUCT"
            else:
                run_decision = "OUT_OF_SCOPE_FOR_MUTATION_POSITION"

            panel_rows.append(
                {
                    "variant_id": variant_id,
                    "aa_change": f"{ref}{pos}{alt}",
                    "position": pos,
                    "source_labels": row.get("source_labels", ""),
                    "inferred_domain": inferred_domain,
                    "assay_key": assay.key,
                    "clinical_axis": assay.clinical_axis,
                    "construct_role": assay.construct_role,
                    "vwf_range": f"{assay.vwf_range[0]}-{assay.vwf_range[1]}",
                    "ligand_keys": "|".join(assay.ligand_keys),
                    "run_decision": run_decision,
                    "expected_signal": assay.expected_signal,
                    "interpretability": assay.interpretability,
                    "notes": assay.notes,
                }
            )

            if run_decision != "RUN" or mutated_seq is None:
                continue

            job_name = f"{variant_id}__{assay.key}"
            if job_name in seen_jobs:
                continue
            seen_jobs.add(job_name)
            sequences = build_sequences(mutated_seq, assay)
            yaml_path = yaml_dir / f"{job_name}.yaml"
            yaml_path.write_text(
                yaml_for_job(
                    job_name,
                    sequences,
                    {
                        "variant_id": variant_id,
                        "aa_change": f"{ref}{pos}{alt}",
                        "source_labels": row.get("source_labels", ""),
                        "assay_key": assay.key,
                        "clinical_axis": assay.clinical_axis,
                        "construct_role": assay.construct_role,
                        "vwf_range": f"{assay.vwf_range[0]}-{assay.vwf_range[1]}",
                        "ligand_keys": "|".join(assay.ligand_keys),
                        "expected_signal": assay.expected_signal,
                        "interpretability": assay.interpretability,
                        "notes": assay.notes,
                    },
                )
            )
            manifest_rows.append(
                {
                    "job_name": job_name,
                    "yaml_path": str(yaml_path.relative_to(output_dir)),
                    "variant_id": variant_id,
                    "aa_change": f"{ref}{pos}{alt}",
                    "position": pos,
                    "wt_aa": ref,
                    "mut_aa": alt,
                    "source_labels": row.get("source_labels", ""),
                    "inferred_domain": inferred_domain,
                    "assay_key": assay.key,
                    "clinical_axis": assay.clinical_axis,
                    "construct_role": assay.construct_role,
                    "vwf_range": f"{assay.vwf_range[0]}-{assay.vwf_range[1]}",
                    "ligand_keys": "|".join(assay.ligand_keys),
                    "n_chains": len(sequences),
                    "run_decision": "RUN",
                    "expected_signal": assay.expected_signal,
                    "interpretability": assay.interpretability,
                    "notes": assay.notes,
                }
            )
            if write_json_batches:
                json_jobs.append(json_job(job_name, sequences))

    variants_out = pd.DataFrame(variant_rows)
    panel_out = pd.DataFrame(panel_rows)
    manifest_out = pd.DataFrame(manifest_rows)

    variants_out.to_csv(output_dir / "variants_master.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    panel_out.to_csv(output_dir / "diagnostic_panel.csv", index=False, quoting=csv.QUOTE_MINIMAL)
    manifest_out.to_csv(output_dir / "job_manifest.csv", index=False, quoting=csv.QUOTE_MINIMAL)

    if write_json_batches:
        batch_dir = output_dir / "json_batches"
        batch_dir.mkdir(parents=True, exist_ok=True)
        for batch_idx, start in enumerate(range(0, len(json_jobs), batch_size), start=1):
            jobs = json_jobs[start : start + batch_size]
            batch_path = batch_dir / f"boltz2_functional_batch_{batch_idx:03d}.json"
            batch_path.write_text(
                json.dumps({"name": f"vwd_functional_panel_{batch_idx:03d}", "jobs": jobs}, indent=2)
            )

    summary = {
        "n_unique_variants": int(len(variants_out)),
        "n_assays": len(ASSAYS),
        "n_yaml_jobs_including_wt": int(len(manifest_out)),
        "n_variant_run_jobs": int((manifest_out["run_decision"] == "RUN").sum()),
        "n_wt_baselines": int((manifest_out["run_decision"] == "WT_BASELINE").sum()),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    return variants_out, panel_out, manifest_out


def print_summary(variants: pd.DataFrame, panel: pd.DataFrame, manifest: pd.DataFrame, output_dir: Path) -> None:
    print("\n=== VWD functional Boltz-2 panel ===")
    print(f"Output: {output_dir}")
    print(f"Unique variants: {len(variants)}")
    print(f"YAML jobs including WT baselines: {len(manifest)}")
    print("\nSource label distribution:")
    print(variants["source_labels"].value_counts().to_string())
    print("\nInferred domain distribution:")
    print(variants["inferred_domain"].value_counts(dropna=False).to_string())
    print("\nRUN jobs by assay:")
    run_counts = manifest[manifest["run_decision"] == "RUN"]["assay_key"].value_counts()
    print(run_counts.to_string() if not run_counts.empty else "No RUN jobs")
    print("\nRun command example:")
    print(f"boltz predict {output_dir / 'yamls'} --out_dir {output_dir / 'boltz_results'} --accelerator gpu --devices 1")
    print("\nMain files:")
    for name in ["variants_master.csv", "diagnostic_panel.csv", "job_manifest.csv", "summary.json"]:
        print(f"  - {output_dir / name}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=REPO_ROOT / "original_patient_table",
        help="Directory containing original patient/control xlsx files.",
    )
    parser.add_argument(
        "--wt-fasta",
        type=Path,
        default=REPO_ROOT / "Proteo-Structure-Pipeline" / "structures" / "wt" / "VWF_P04275_WT.fasta",
        help="WT VWF P04275 FASTA path.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "output" / "boltz2_vwd_functional_panel",
        help="Output directory.",
    )
    parser.add_argument("--limit-variants", type=int, default=None, help="Optional quick-test limit.")
    parser.add_argument("--write-json-batches", action="store_true", help="Also write JSON batches.")
    parser.add_argument("--batch-size", type=int, default=30, help="Jobs per JSON batch.")
    parser.add_argument("--variants-csv", type=Path, default=None,
                        help="Load variants from a clean CSV (aa_change/wt_aa/position/mut_aa"
                             "[,subtype]) instead of the patient xlsx tables.")
    parser.add_argument("--subtype-filter", default=None,
                        help="With --variants-csv, keep only rows whose subtype == this (e.g. 2B).")
    args = parser.parse_args()

    wt_seq = read_fasta(args.wt_fasta)
    if args.variants_csv is not None:
        variants = load_variants_csv(args.variants_csv, args.subtype_filter)
        src = f"{args.variants_csv}" + (f" (subtype={args.subtype_filter})" if args.subtype_filter else "")
    else:
        variants = load_original_tables(args.input_dir)
        src = str(args.input_dir)
    if variants.empty:
        raise SystemExit(f"No missense variants parsed from {src}")
    print(f"Loaded {len(variants)} variants from {src}")

    variants_out, panel_out, manifest_out = generate_panel(
        variants=variants,
        wt_seq=wt_seq,
        output_dir=args.output_dir,
        write_json_batches=args.write_json_batches,
        batch_size=args.batch_size,
        limit_variants=args.limit_variants,
    )
    print_summary(variants_out, panel_out, manifest_out, args.output_dir)


if __name__ == "__main__":
    main()

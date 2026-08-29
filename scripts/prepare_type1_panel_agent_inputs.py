#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
HALF_WINDOW = 2**19
ALPHAGENOME_RAW_OUTPUTS = (
    "ATAC", "CAGE", "DNASE", "RNA_SEQ", "CHIP_HISTONE", "CHIP_TF",
    "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS",
    "CONTACT_MAPS", "PROCAP",
)
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def clean(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return re.sub(r"\s+", "", str(value).strip())


def display(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


def numeric(value: Any) -> float | str:
    text = display(value)
    if not text or text.casefold() in {"na", "n/a", "nan"}:
        return "NA"
    if text.startswith("<"):
        return text
    try:
        return float(text)
    except ValueError:
        return text


def compact_missense(hgvs_p: str) -> str:
    match = re.fullmatch(
        r"p\.(" + "|".join(AA3_TO_1) + r")(\d+)(" + "|".join(AA3_TO_1) + r")",
        hgvs_p,
    )
    if not match:
        return ""
    ref, position, alt = match.groups()
    return f"{AA3_TO_1[ref]}{position}{AA3_TO_1[alt]}"


def normalized_gene(value: Any) -> str:
    text = clean(value).upper()
    return "VWF" if text == "VWF" else text


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare de-identified Type-1 panel inputs for VWD agents and computational panels.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--normalization", type=Path, default=ROOT / "data/type1_panel_20260828/normalization_grch38.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw = pd.read_excel(args.input, sheet_name="Sheet1", header=None, keep_default_na=False)
    source = raw.iloc[2:].reset_index(drop=True)
    norm_payload = json.loads(args.normalization.read_text(encoding="utf-8"))
    norm = norm_payload["records"]

    pre_rows: list[dict[str, Any]] = []
    post_rows: list[dict[str, Any]] = []
    context_rows: list[dict[str, Any]] = []
    request_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    boltz_rows: list[dict[str, Any]] = []

    for zero_index, row in source.iterrows():
        source_row = zero_index + 1
        case_id = f"CASE_T1_{source_row:03d}"
        reported_gene = normalized_gene(row.iloc[11])
        hgvs_c = clean(row.iloc[14])
        hgvs_p = clean(row.iloc[15]).replace("Cys", "Cys")
        source_coordinate = clean(row.iloc[16])
        variant_type = display(row.iloc[17])
        zygosity = display(row.iloc[18])
        # CASE_T1_007 was recruited because VWD coexists with hemophilia A.  The
        # source gene cell records the comorbidity, while its chr12 coordinate
        # and HGVS are the VWF splice variant NM_000552.5:c.3379+1G>A.
        f8_comorbidity_with_vwf_variant = (
            reported_gene == "F8"
            and hgvs_c == "c.3379+1G>A"
            and "6132796" in source_coordinate
        )
        gene = "VWF" if f8_comorbidity_with_vwf_variant else reported_gene
        comorbidity = "Hemophilia A" if f8_comorbidity_with_vwf_variant else ""
        interpretation_constraint = (
            "FVIII:C is confounded by coexisting hemophilia A and must not support or refute VWD type 2N."
            if f8_comorbidity_with_vwf_variant else ""
        )
        normalization = norm.get(hgvs_c) if gene == "VWF" else None
        compact = compact_missense(hgvs_p)
        is_missense = bool(compact)
        alpha_status = "READY" if normalization else ("OUT_OF_SCOPE_NON_VWF" if gene != "VWF" else "NEEDS_HGVS_NORMALIZATION")
        boltz_status = "READY" if gene == "VWF" and is_missense else ("OUT_OF_SCOPE_NON_VWF" if gene != "VWF" else "NOT_APPLICABLE_NON_MISSENSE")
        coordinate_qc = "PASS"
        if normalization and source_coordinate:
            digits = re.findall(r"\d+", source_coordinate)
            if not digits:
                coordinate_qc = "SOURCE_COORDINATE_UNPARSEABLE"
            elif hgvs_c == "c.5827C>T":
                coordinate_qc = "SOURCE_COORDINATE_CONFLICT_WITH_HGVS"
            else:
                coordinate_qc = "SOURCE_HG19_TO_HG38_NORMALIZED"
        elif normalization and not source_coordinate:
            coordinate_qc = "SOURCE_COORDINATE_MISSING_HGVS_RESOLVED"
        elif gene != "VWF":
            coordinate_qc = "NON_VWF_GENE_OUT_OF_SCOPE"
        else:
            coordinate_qc = "NORMALIZATION_NOT_AVAILABLE"

        pre_rows.append({
            "编号": source_row,
            "研究病例ID": case_id,
            "VWF:Ag": numeric(row.iloc[27]),
            "VWF:Act": numeric(row.iloc[28]),
            "FVIII:C": numeric(row.iloc[30]),
            "血小板计数": "NA",
            "VWF:Ag单位": "%",
            "VWF:Act单位": "%",
            "FVIII:C单位": "%",
            "血小板计数单位": "10^9/L",
            "VWF:Ag参考范围": "来源表：O型42.0-140.8%；非O型66.1-79.9%（需核对原实验室）",
            "VWF:Act参考范围": "来源表：O型40.3-125.9%；非O型48.8-163.4%",
            "FVIII:C参考范围": "50-150%",
            "血小板计数参考范围": "未提供",
        })
        post_rows.append({
            "编号": source_row,
            "研究病例ID": case_id,
            "突变基因": gene,
            "核苷酸改变": hgvs_c or "NA",
            "氨基酸改变": hgvs_p or "NA",
            "染色体位置": source_coordinate or "NA",
            "变异类型": variant_type or "NA",
            "合子状态": zygosity or "NA",
            "报告相位": display(row.iloc[25]) or "NA",
            "报告ACMG": display(row.iloc[24]) or "NA",
            "归一化基因组版本": "GRCh38" if normalization else "NA",
            "hg38染色体": normalization["chromosome"] if normalization else "NA",
            "hg38位置": normalization["position"] if normalization else "NA",
            "hg38_REF": normalization["ref"] if normalization else "NA",
            "hg38_ALT": normalization["alt"] if normalization else "NA",
            "蛋白紧凑表示": compact or "NA",
            "AlphaGenome请求状态": alpha_status,
            "Boltz请求状态": boltz_status,
            "坐标QC": coordinate_qc,
            "归一化来源": "Ensembl VEP NM_000552.5 / GRCh38" if normalization else "NA",
            "原表基因字段": reported_gene or "NA",
            "共病": comorbidity or "NA",
            "解释限制": interpretation_constraint or "NA",
        })
        context_rows.append({
            "研究病例ID": case_id,
            "性别": display(row.iloc[10]) or "NA",
            "年龄": display(row.iloc[5]) or "NA",
            "病程": display(row.iloc[4]) or "NA",
            "主要症状": display(row.iloc[7]) or "NA",
            "家族史": display(row.iloc[9]) or "NA",
            "ISTH-BAT": numeric(row.iloc[59]),
            "高浓度瑞斯托霉素血小板聚集": numeric(row.iloc[31]),
            "DDAVP_VWF:Ag_输注前": numeric(row.iloc[32]),
            "DDAVP_VWF:Ag_0.5h": numeric(row.iloc[33]),
            "DDAVP_VWF:Ag_1h": numeric(row.iloc[34]),
            "DDAVP_VWF:Ag_4h": numeric(row.iloc[35]),
            "DDAVP_VWF:Act_输注前": numeric(row.iloc[36]),
            "DDAVP_VWF:Act_0.5h": numeric(row.iloc[37]),
            "DDAVP_VWF:Act_1h": numeric(row.iloc[38]),
            "DDAVP_VWF:Act_4h": numeric(row.iloc[39]),
            "DDAVP_FVIII:C_输注前": numeric(row.iloc[40]),
            "DDAVP_FVIII:C_0.5h": numeric(row.iloc[41]),
            "DDAVP_FVIII:C_1h": numeric(row.iloc[42]),
            "DDAVP_FVIII:C_4h": numeric(row.iloc[43]),
            "DDAVP报告反应性": display(row.iloc[44]) or "NA",
            "既往治疗": display(row.iloc[60]) or "NA",
            "共病": comorbidity or "NA",
            "解释限制": interpretation_constraint or "NA",
        })

        request_rows.append({
            "研究病例ID": case_id,
            "基因": gene,
            "HGVS_c": hgvs_c or "NA",
            "HGVS_p": hgvs_p or "NA",
            "AlphaGenome状态": alpha_status,
            "AlphaGenome_GRCh38": (
                f"{normalization['chromosome']}:{normalization['position']}:{normalization['ref']}>{normalization['alt']}"
                if normalization else "NA"
            ),
            "Boltz状态": boltz_status,
            "Boltz变异": compact or "NA",
            "既有AlphaGenome结果": "待匹配",
            "既有Boltz结果": "待匹配",
            "AlphaGenome原始输出": "|".join(ALPHAGENOME_RAW_OUTPUTS) if alpha_status == "READY" else "NA",
            "AlphaGenome细胞选择": "运行时读取output_metadata；优先内皮/血管相关轨道" if alpha_status == "READY" else "NA",
        })
        audit_rows.append({
            "研究病例ID": case_id,
            "源行": source_row,
            "原始基因": display(row.iloc[11]),
            "标准基因": gene,
            "原始坐标(hg19/未声明)": source_coordinate or "NA",
            "HGVS归一化GRCh38": request_rows[-1]["AlphaGenome_GRCh38"],
            "坐标QC": coordinate_qc,
            "说明": (
                "原表F8指血友病A共病；该chr12/HGVS变异按VWF归一化，FVIII:C不用于2N推断" if f8_comorbidity_with_vwf_variant else
                "非VWF病例不进入VWF模型" if gene != "VWF" else
                "原表坐标与 HGVS 映射不一致，以 HGVS 归一化为请求基准" if coordinate_qc == "SOURCE_COORDINATE_CONFLICT_WITH_HGVS" else
                "保留原值；请求使用归一化坐标"
            ),
        })
        if boltz_status == "READY":
            boltz_rows.append({
                "aa_change": compact,
                "wt_aa": compact[0],
                "position": int(compact[1:-1]),
                "mut_aa": compact[-1],
                "subtype": "unknown",
                "domain": "",
                "case_id": case_id,
            })

    pd.DataFrame(pre_rows).to_csv(args.output_dir / "first_level.csv", index=False)
    pd.DataFrame(post_rows).to_csv(args.output_dir / "genetics.csv", index=False)
    pd.DataFrame(context_rows).to_csv(args.output_dir / "clinical_context.csv", index=False)
    pd.DataFrame(request_rows).to_csv(args.output_dir / "panel_request_status.csv", index=False)
    pd.DataFrame(audit_rows).to_csv(args.output_dir / "normalization_audit.csv", index=False)
    pd.DataFrame(boltz_rows).to_csv(args.output_dir / "boltz_variants.csv", index=False)

    alpha = [
        {
            "case_id": row["研究病例ID"],
            "genome_build": "GRCh38",
            "chromosome": row["hg38染色体"],
            "position": int(row["hg38位置"]),
            "ref": row["hg38_REF"],
            "alt": row["hg38_ALT"],
            "interval_start": max(1, int(row["hg38位置"]) - HALF_WINDOW),
            "interval_end": int(row["hg38位置"]) + HALF_WINDOW,
            "ontology_terms": "CL:0000115",
            "ontology_strategy": "metadata_selected_endothelial_or_global",
            "requested_outputs": "|".join(ALPHAGENOME_RAW_OUTPUTS),
            "hgvs_c": row["核苷酸改变"],
            "hgvs_p": row["氨基酸改变"],
        }
        for row in post_rows if row["AlphaGenome请求状态"] == "READY"
    ]
    pd.DataFrame(alpha).to_csv(args.output_dir / "alphagenome_requests.csv", index=False)
    with (args.output_dir / "alphagenome_requests.vcf").open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in alpha:
            handle.write(
                f"{str(row['chromosome']).removeprefix('chr')}\t{row['position']}\t{row['case_id']}\t"
                f"{row['ref']}\t{row['alt']}\t.\tPASS\tHGVS_C={row['hgvs_c']};HGVS_P={row['hgvs_p']}\n"
            )

    payload = {
        "source_workbook": args.input.name,
        "privacy": "Direct identifiers, dates of birth, visit dates, and medical record numbers were excluded.",
        "first_level": pre_rows,
        "genetics": post_rows,
        "clinical_context": context_rows,
        "requests": request_rows,
        "audit": audit_rows,
        "alphagenome_profile": {
            "raw_output_count": len(ALPHAGENOME_RAW_OUTPUTS),
            "raw_outputs": list(ALPHAGENOME_RAW_OUTPUTS),
            "preferred_ontology_term": "CL:0000115",
            "ontology_strategy": "Query output_metadata at runtime; use endothelial/vascular tracks per output and retain global non-biosample outputs.",
        },
    }
    (args.output_dir / "prepared_payload.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = {
        "cases": len(pre_rows),
        "alphagenome_ready": sum(row["AlphaGenome请求状态"] == "READY" for row in post_rows),
        "boltz_ready": sum(row["Boltz请求状态"] == "READY" for row in post_rows),
        "output_dir": str(args.output_dir.resolve()),
    }
    if summary["cases"] != 10 or summary["alphagenome_ready"] != 9 or summary["boltz_ready"] != 5:
        raise RuntimeError(f"Unexpected Type-1 panel inventory: {summary}")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

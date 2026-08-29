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
RAW_OUTPUTS = (
    "ATAC", "CAGE", "DNASE", "RNA_SEQ", "CHIP_HISTONE", "CHIP_TF",
    "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS", "CONTACT_MAPS", "PROCAP",
)
AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}


def text(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


def clean(value: Any) -> str:
    return re.sub(r"\s+", "", text(value)).replace("＞", ">").replace("c.4021c>T", "c.4021C>T")


def numeric(value: Any) -> float | str:
    value = text(value)
    if not value or value.casefold() in {"na", "n/a", "nan"}:
        return "NA"
    try:
        return float(value)
    except ValueError:
        return value


def compact_missense(hgvs_p: str) -> str:
    short = hgvs_p.removeprefix("p.")
    match = re.fullmatch(
        r"(" + "|".join(AA3_TO_1) + r"|[A-Z])(\d+)(" + "|".join(AA3_TO_1) + r"|[A-Z])",
        short,
    )
    if not match:
        return ""
    ref, position, alt = match.groups()
    return f"{AA3_TO_1.get(ref, ref)}{position}{AA3_TO_1.get(alt, alt)}"


def source_case_id(patient_number: int, variant_number: int, benign: bool) -> str:
    base = f"CASE_T2B_{patient_number:03d}"
    if patient_number != 5:
        return base
    suffix = f"_VARIANT_{variant_number}"
    return base + suffix + ("_BENIGN" if benign else "")


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare de-identified Type-2B panel inputs for the VWD Agent and compute panels.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument(
        "--normalization",
        type=Path,
        default=ROOT / "data/type2b_panel_20260829/normalization_grch38.json",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_excel(args.input, sheet_name="Sheet1", header=None, keep_default_na=False)
    source = raw.iloc[2:].reset_index(drop=True)
    norm = json.loads(args.normalization.read_text(encoding="utf-8"))["records"]

    first_by_patient: dict[int, dict[str, Any]] = {}
    context_by_patient: dict[int, dict[str, Any]] = {}
    genetics: list[dict[str, Any]] = []
    requests: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    alpha: list[dict[str, Any]] = []
    boltz: list[dict[str, Any]] = []

    patient_variant_count: dict[int, int] = {}
    for _, row in source.iterrows():
        patient_number = int(row.iloc[0])
        patient_variant_count[patient_number] = patient_variant_count.get(patient_number, 0) + 1
        variant_number = patient_variant_count[patient_number]
        benign = "良性" in text(row.iloc[1]) or text(row.iloc[37]).upper() in {"B", "LB"}
        case_ref = source_case_id(patient_number, variant_number, benign)
        patient_id = f"CASE_T2B_{patient_number:03d}"
        hgvs_c = clean(row.iloc[30])
        hgvs_p = clean(row.iloc[31])
        compact = compact_missense(hgvs_p)
        normalized = norm.get(hgvs_c)
        source_coordinate = clean(row.iloc[32])
        normalized_locus = f"{normalized['chromosome']}:{normalized['position']}" if normalized else ""
        coordinate_qc = (
            "SOURCE_ALREADY_GRCH38_MATCHED"
            if normalized and source_coordinate.replace("-", ":") == normalized_locus
            else "SOURCE_TO_GRCH38_NORMALIZED" if normalized
            else "NORMALIZATION_NOT_AVAILABLE"
        )
        alpha_status = "READY" if normalized else "NEEDS_HGVS_NORMALIZATION"
        boltz_status = "READY" if compact else "NOT_APPLICABLE_NON_MISSENSE"

        if patient_number not in first_by_patient or variant_number == 1:
            first_by_patient[patient_number] = {
                "编号": patient_number,
                "研究病例ID": patient_id,
                "VWF:Ag": numeric(row.iloc[48]),
                "VWF:Act": numeric(row.iloc[49]),
                "FVIII:C": numeric(row.iloc[52]),
                "血小板计数": numeric(row.iloc[40]),
                "VWF:Ag单位": "%", "VWF:Act单位": "%", "FVIII:C单位": "%",
                "血小板计数单位": "10^9/L",
                "VWF:Ag参考范围": "O型42.0-140.8%；非O型66.1-79.9%（来源表）",
                "VWF:Act参考范围": "O型40.3-125.9%；非O型48.8-163.4%（来源表）",
                "FVIII:C参考范围": "50-150%", "血小板计数参考范围": "来源表未提供",
            }
            context_by_patient[patient_number] = {
                "研究病例ID": patient_id,
                "性别": "男" if numeric(row.iloc[11]) == 1.0 else "女",
                "年龄": "NA",
                "病程": "NA",
                "主要症状": text(row.iloc[4]) or "NA",
                "家族史": text(row.iloc[10]) or "NA",
                "ISTH-BAT": numeric(row.iloc[27]),
                "高浓度瑞斯托霉素血小板聚集": numeric(row.iloc[56]),
                "低浓度瑞斯托霉素血小板聚集": numeric(row.iloc[57]),
                "既往治疗": text(row.iloc[9]) or "NA",
                "父母表型": text(row.iloc[39]) or "NA",
                "共病": "NA",
                "解释限制": "低血小板计数与低浓度RIPA是2B轴重要证据；结构模型只作机制支持。",
            }

        genetics.append({
            "编号": patient_number,
            "研究病例ID": case_ref,
            "突变基因": "VWF",
            "核苷酸改变": hgvs_c or "NA",
            "氨基酸改变": hgvs_p or "NA",
            "染色体位置": source_coordinate or "NA",
            "变异类型": text(row.iloc[33]) or "NA",
            "合子状态": {"1": "杂合", "2": "半合子", "3": "纯合", "4": "复合杂合"}.get(text(row.iloc[34]), text(row.iloc[34]) or "NA"),
            "报告相位": "NA",
            "报告ACMG": text(row.iloc[36]) or "NA",
            "归一化基因组版本": "GRCh38" if normalized else "NA",
            "hg38染色体": normalized["chromosome"] if normalized else "NA",
            "hg38位置": normalized["position"] if normalized else "NA",
            "hg38_REF": normalized["ref"] if normalized else "NA",
            "hg38_ALT": normalized["alt"] if normalized else "NA",
            "蛋白紧凑表示": compact or "NA",
            "AlphaGenome请求状态": alpha_status,
            "Boltz请求状态": boltz_status,
            "坐标QC": coordinate_qc,
            "归一化来源": "Ensembl VEP NM_000552.5 / GRCh38" if normalized else "NA",
            "报告频率": numeric(row.iloc[35]),
            "报告功能预测": text(row.iloc[37]) or "NA",
            "变异来源": text(row.iloc[38]) or "NA",
            "良性标记": benign,
        })
        requests.append({
            "研究病例ID": case_ref,
            "患者ID": patient_id,
            "基因": "VWF",
            "HGVS_c": hgvs_c or "NA",
            "HGVS_p": hgvs_p or "NA",
            "AlphaGenome状态": alpha_status,
            "AlphaGenome_GRCh38": f"{normalized['chromosome']}:{normalized['position']}:{normalized['ref']}>{normalized['alt']}" if normalized else "NA",
            "Boltz状态": boltz_status,
            "Boltz变异": compact or "NA",
            "AlphaGenome原始输出": "|".join(RAW_OUTPUTS),
            "AlphaGenome细胞选择": "运行时读取output_metadata；优先内皮/血管相关轨道",
            "备注": "保留为同一患者第二变异；不推断相位" if patient_number == 5 and variant_number == 2 else "NA",
        })
        audit.append({
            "研究病例ID": case_ref,
            "患者ID": patient_id,
            "源编号": patient_number,
            "已删除直接标识": "姓名|病案号|出生日期|就诊日期",
            "源坐标": source_coordinate or "NA",
            "HGVS归一化GRCh38": requests[-1]["AlphaGenome_GRCh38"],
            "坐标QC": genetics[-1]["坐标QC"],
            "说明": (
                "源坐标与GRCh38一致；c.4021c>T已规范为c.4021C>T"
                if coordinate_qc == "SOURCE_ALREADY_GRCH38_MATCHED"
                else "请求使用HGVS归一化坐标"
            ),
        })
        if normalized:
            alpha.append({
                "case_id": case_ref,
                "patient_id": patient_id,
                "genome_build": "GRCh38",
                "chromosome": normalized["chromosome"],
                "position": normalized["position"],
                "ref": normalized["ref"], "alt": normalized["alt"],
                "interval_start": max(1, normalized["position"] - HALF_WINDOW),
                "interval_end": normalized["position"] + HALF_WINDOW,
                "ontology_terms": "CL:0000115",
                "ontology_strategy": "metadata_selected_endothelial_or_global",
                "requested_outputs": "|".join(RAW_OUTPUTS),
                "hgvs_c": hgvs_c, "hgvs_p": hgvs_p,
            })
        if compact:
            boltz.append({
                "aa_change": compact, "wt_aa": compact[0], "position": int(compact[1:-1]),
                "mut_aa": compact[-1], "subtype": "2B" if not benign else "benign_control",
                "domain": "", "case_id": case_ref, "patient_id": patient_id,
            })

    first = [first_by_patient[key] for key in sorted(first_by_patient)]
    context = [context_by_patient[key] for key in sorted(context_by_patient)]
    pd.DataFrame(first).to_csv(args.output_dir / "first_level.csv", index=False)
    pd.DataFrame(genetics).to_csv(args.output_dir / "genetics.csv", index=False)
    pd.DataFrame(context).to_csv(args.output_dir / "clinical_context.csv", index=False)
    pd.DataFrame(requests).to_csv(args.output_dir / "panel_request_status.csv", index=False)
    pd.DataFrame(audit).to_csv(args.output_dir / "normalization_audit.csv", index=False)
    pd.DataFrame(alpha).to_csv(args.output_dir / "alphagenome_requests.csv", index=False)
    pd.DataFrame(boltz).to_csv(args.output_dir / "boltz_variants.csv", index=False)
    with (args.output_dir / "alphagenome_requests.vcf").open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for item in alpha:
            handle.write(f"{item['chromosome'].removeprefix('chr')}\t{item['position']}\t{item['case_id']}\t{item['ref']}\t{item['alt']}\t.\tPASS\tHGVS_C={item['hgvs_c']};HGVS_P={item['hgvs_p']}\n")
    payload = {
        "source_workbook": args.input.name,
        "privacy": "Names, medical record numbers, birth dates, and visit dates were excluded.",
        "first_level": first, "genetics": genetics, "clinical_context": context,
        "requests": requests, "audit": audit,
    }
    (args.output_dir / "prepared_payload.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = {
        "patients": len(first), "variant_rows": len(genetics), "alphagenome_ready": len(alpha),
        "boltz_variant_rows": len(boltz), "unique_boltz_variants": len({row['aa_change'] for row in boltz}),
    }
    if summary != {"patients": 6, "variant_rows": 7, "alphagenome_ready": 7, "boltz_variant_rows": 7, "unique_boltz_variants": 6}:
        raise RuntimeError(f"Unexpected Type-2B inventory: {summary}")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

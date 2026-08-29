from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import re
from typing import Any

import pandas as pd

from .schemas import (
    ClinicalContext,
    DDAVPSeries,
    FirstLevelLabs,
    MissingReason,
    ObservedValue,
    PatientCase,
    VariantContext,
)


_MISSING_TOKENS = {
    "": MissingReason.NOT_RECORDED,
    "na": MissingReason.NOT_AVAILABLE,
    "n/a": MissingReason.NOT_AVAILABLE,
    "nan": MissingReason.NOT_AVAILABLE,
    "not_done": MissingReason.NOT_DONE,
    "not_available": MissingReason.NOT_AVAILABLE,
    "pending": MissingReason.PENDING,
}


def parse_case_ref(raw: Any) -> tuple[str, int, bool]:
    """Normalize CASE_###[_VARIANT_n][_BENIGN] identifiers."""
    text = str(raw).strip()
    match = re.fullmatch(
        r"(CASE(?:_[A-Z0-9]+)*?_\d+)(?:_VARIANT_(\d+))?(?:_(BENIGN))?",
        text,
        flags=re.IGNORECASE,
    )
    if not match:
        raise ValueError(f"Unexpected case identifier: {text!r}")
    patient_id = match.group(1).upper()
    variant_index = int(match.group(2) or 1)
    benign = bool(match.group(3))
    return patient_id, variant_index, benign


def observed_value(cell: Any, *, unit: str | None = None, reference_range: str | None = None) -> ObservedValue:
    if cell is None:
        return ObservedValue(raw="", observed=False, unit=unit, reference_range=reference_range, missing_reason=MissingReason.NOT_RECORDED)
    if isinstance(cell, (int, float)):
        return ObservedValue(raw=str(cell), value=float(cell), observed=True, unit=unit, reference_range=reference_range)
    raw = str(cell).strip()
    missing = _MISSING_TOKENS.get(raw.casefold())
    if missing is not None:
        return ObservedValue(raw=raw, observed=False, unit=unit, reference_range=reference_range, missing_reason=missing)
    try:
        return ObservedValue(raw=raw, value=float(raw), observed=True, unit=unit, reference_range=reference_range)
    except ValueError:
        return ObservedValue(raw=raw, observed=False, unit=unit, reference_range=reference_range, missing_reason=MissingReason.NOT_RECORDED)


def _optional_text(row: dict[str, Any], column: str) -> str | None:
    value = row.get(column, "")
    text = str(value).strip()
    return None if not text or text.casefold() in {"na", "n/a", "nan", "none"} else text


def _optional_int(row: dict[str, Any], column: str) -> int | None:
    value = row.get(column, "")
    if value in (None, ""):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


class LocalWorkbookProvider:
    name = "local_workbook"
    version = "deidentified_v3"

    def __init__(self, workbook_path: str | Path):
        self.path = Path(workbook_path)
        self.sha256 = sha256(self.path.read_bytes()).hexdigest()

    def load_cases(self) -> tuple[list[PatientCase], dict[str, Any]]:
        pre = pd.read_excel(self.path, sheet_name="1.基因前", keep_default_na=False)
        post = pd.read_excel(self.path, sheet_name="2.基因后", keep_default_na=False)
        context = (
            pd.read_excel(self.path, sheet_name="3.临床上下文", keep_default_na=False)
            if "3.临床上下文" in pd.ExcelFile(self.path).sheet_names
            else pd.DataFrame()
        )
        pre.columns = [str(column).strip() for column in pre.columns]
        post.columns = [str(column).strip() for column in post.columns]
        if not context.empty:
            context.columns = [str(column).strip() for column in context.columns]
        context_groups: dict[str, dict[str, Any]] = {}
        for _, row in context.iterrows():
            patient_id, _, _ = parse_case_ref(row["研究病例ID"])
            context_groups[patient_id] = row.to_dict()

        pre_groups: dict[str, list[dict[str, Any]]] = {}
        post_groups: dict[str, list[dict[str, Any]]] = {}
        for _, row in pre.iterrows():
            patient_id, _, _ = parse_case_ref(row["研究病例ID"])
            pre_groups.setdefault(patient_id, []).append(row.to_dict())
        for _, row in post.iterrows():
            patient_id, _, _ = parse_case_ref(row["研究病例ID"])
            post_groups.setdefault(patient_id, []).append(row.to_dict())

        if set(pre_groups) != set(post_groups):
            raise ValueError(
                "First-level/genetic patient sets differ: "
                f"only_pre={sorted(set(pre_groups) - set(post_groups))}, "
                f"only_post={sorted(set(post_groups) - set(pre_groups))}"
            )

        cases: list[PatientCase] = []
        for patient_id in sorted(pre_groups):
            pre_rows = pre_groups[patient_id]
            post_rows = post_groups[patient_id]
            first = pre_rows[0]
            first_level = FirstLevelLabs(
                vwf_ag=observed_value(first["VWF:Ag"], unit=_optional_text(first, "VWF:Ag单位"), reference_range=_optional_text(first, "VWF:Ag参考范围")),
                vwf_act=observed_value(first["VWF:Act"], unit=_optional_text(first, "VWF:Act单位"), reference_range=_optional_text(first, "VWF:Act参考范围")),
                fviii_c=observed_value(first["FVIII:C"], unit=_optional_text(first, "FVIII:C单位"), reference_range=_optional_text(first, "FVIII:C参考范围")),
                platelet_count=observed_value(first["血小板计数"], unit=_optional_text(first, "血小板计数单位"), reference_range=_optional_text(first, "血小板计数参考范围")),
            )
            variants = [
                VariantContext(
                    source_row_id=str(row["研究病例ID"]).strip(),
                    variant_index=variant_index,
                    hgvs_c=_optional_text(row, "核苷酸改变"),
                    hgvs_p=_optional_text(row, "氨基酸改变"),
                    chromosomal_position=_optional_text(row, "染色体位置"),
                    gene=_optional_text(row, "突变基因") or "VWF",
                    variant_type=_optional_text(row, "变异类型"),
                    zygosity=_optional_text(row, "合子状态"),
                    reported_phase=_optional_text(row, "报告相位"),
                    reported_acmg=_optional_text(row, "报告ACMG"),
                    genome_build=_optional_text(row, "归一化基因组版本"),
                    hg38_chromosome=_optional_text(row, "hg38染色体"),
                    hg38_position=_optional_int(row, "hg38位置"),
                    hg38_ref=_optional_text(row, "hg38_REF"),
                    hg38_alt=_optional_text(row, "hg38_ALT"),
                    alphagenome_request_status=_optional_text(row, "AlphaGenome请求状态"),
                    boltz_request_status=_optional_text(row, "Boltz请求状态"),
                    benign_reported=benign,
                )
                for row in post_rows
                for _, variant_index, benign in [parse_case_ref(row["研究病例ID"])]
            ]
            variants.sort(key=lambda variant: variant.variant_index)
            if len(variants) > 1 or any("复合" in (variant.zygosity or "") for variant in variants):
                for variant in variants:
                    variant.phase_status = "unknown"
            context_row = context_groups.get(patient_id, {})
            ddavp = DDAVPSeries(
                vwf_ag_pre=observed_value(context_row.get("DDAVP_VWF:Ag_输注前")),
                vwf_ag_0_5h=observed_value(context_row.get("DDAVP_VWF:Ag_0.5h")),
                vwf_ag_1h=observed_value(context_row.get("DDAVP_VWF:Ag_1h")),
                vwf_ag_4h=observed_value(context_row.get("DDAVP_VWF:Ag_4h")),
                vwf_act_pre=observed_value(context_row.get("DDAVP_VWF:Act_输注前")),
                vwf_act_0_5h=observed_value(context_row.get("DDAVP_VWF:Act_0.5h")),
                vwf_act_1h=observed_value(context_row.get("DDAVP_VWF:Act_1h")),
                vwf_act_4h=observed_value(context_row.get("DDAVP_VWF:Act_4h")),
                fviii_c_pre=observed_value(context_row.get("DDAVP_FVIII:C_输注前")),
                fviii_c_0_5h=observed_value(context_row.get("DDAVP_FVIII:C_0.5h")),
                fviii_c_1h=observed_value(context_row.get("DDAVP_FVIII:C_1h")),
                fviii_c_4h=observed_value(context_row.get("DDAVP_FVIII:C_4h")),
                reported_response=_optional_text(context_row, "DDAVP报告反应性"),
            )
            clinical_context = ClinicalContext(
                sex=_optional_text(context_row, "性别"),
                age_text=_optional_text(context_row, "年龄"),
                disease_course=_optional_text(context_row, "病程"),
                symptoms=_optional_text(context_row, "主要症状"),
                family_history=_optional_text(context_row, "家族史"),
                isth_bat=observed_value(context_row.get("ISTH-BAT")),
                high_dose_ristocetin=observed_value(context_row.get("高浓度瑞斯托霉素血小板聚集"), unit="%"),
                ddavp=ddavp,
                prior_treatment=_optional_text(context_row, "既往治疗"),
                comorbidity=_optional_text(context_row, "共病"),
                interpretation_constraints=[
                    constraint
                    for constraint in [_optional_text(context_row, "解释限制")]
                    if constraint
                ],
            )
            cases.append(
                PatientCase(
                    patient_id=patient_id,
                    episode_id=f"{patient_id}_EPISODE_1",
                    source_row_ids=sorted({str(row["研究病例ID"]).strip() for row in pre_rows + post_rows}),
                    first_level=first_level,
                    variants=variants,
                    clinical_context=clinical_context,
                )
            )

        audit = {
            "workbook_path": str(self.path),
            "workbook_sha256": self.sha256,
            "first_level_rows": len(pre),
            "variant_rows": len(post),
            "unique_patients": len(cases),
            "multi_variant_patients": sum(len(case.variants) > 1 for case in cases),
            "clinical_context_rows": len(context),
        }
        return cases, audit

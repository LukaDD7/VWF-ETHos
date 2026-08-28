from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import re
from typing import Any

import pandas as pd

from .schemas import FirstLevelLabs, MissingReason, ObservedValue, PatientCase, VariantContext


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
    match = re.fullmatch(r"(CASE_\d+)(?:_VARIANT_(\d+))?(?:_(BENIGN))?", text, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Unexpected case identifier: {text!r}")
    patient_id = match.group(1).upper()
    variant_index = int(match.group(2) or 1)
    benign = bool(match.group(3))
    return patient_id, variant_index, benign


def observed_value(cell: Any) -> ObservedValue:
    if cell is None:
        return ObservedValue(raw="", observed=False, missing_reason=MissingReason.NOT_RECORDED)
    if isinstance(cell, (int, float)):
        return ObservedValue(raw=str(cell), value=float(cell), observed=True)
    raw = str(cell).strip()
    missing = _MISSING_TOKENS.get(raw.casefold())
    if missing is not None:
        return ObservedValue(raw=raw, observed=False, missing_reason=missing)
    try:
        return ObservedValue(raw=raw, value=float(raw), observed=True)
    except ValueError:
        return ObservedValue(raw=raw, observed=False, missing_reason=MissingReason.NOT_RECORDED)


class LocalWorkbookProvider:
    name = "local_workbook"
    version = "deidentified_v3"

    def __init__(self, workbook_path: str | Path):
        self.path = Path(workbook_path)
        self.sha256 = sha256(self.path.read_bytes()).hexdigest()

    def load_cases(self) -> tuple[list[PatientCase], dict[str, Any]]:
        pre = pd.read_excel(self.path, sheet_name="1.基因前", keep_default_na=False)
        post = pd.read_excel(self.path, sheet_name="2.基因后", keep_default_na=False)
        pre.columns = [str(column).strip() for column in pre.columns]
        post.columns = [str(column).strip() for column in post.columns]

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
                vwf_ag=observed_value(first["VWF:Ag"]),
                vwf_act=observed_value(first["VWF:Act"]),
                fviii_c=observed_value(first["FVIII:C"]),
                platelet_count=observed_value(first["血小板计数"]),
            )
            variants = [
                VariantContext(
                    source_row_id=str(row["研究病例ID"]).strip(),
                    variant_index=variant_index,
                    hgvs_c=str(row["核苷酸改变"]).strip() or None,
                    hgvs_p=str(row["氨基酸改变"]).strip() or None,
                    chromosomal_position=str(row["染色体位置"]).strip() or None,
                    benign_reported=benign,
                )
                for row in post_rows
                for _, variant_index, benign in [parse_case_ref(row["研究病例ID"])]
            ]
            variants.sort(key=lambda variant: variant.variant_index)
            if len(variants) > 1:
                for variant in variants:
                    variant.phase_status = "unknown"
            cases.append(
                PatientCase(
                    patient_id=patient_id,
                    episode_id=f"{patient_id}_EPISODE_1",
                    source_row_ids=sorted({str(row["研究病例ID"]).strip() for row in pre_rows + post_rows}),
                    first_level=first_level,
                    variants=variants,
                )
            )

        audit = {
            "workbook_path": str(self.path),
            "workbook_sha256": self.sha256,
            "first_level_rows": len(pre),
            "variant_rows": len(post),
            "unique_patients": len(cases),
            "multi_variant_patients": sum(len(case.variants) > 1 for case in cases),
        }
        return cases, audit

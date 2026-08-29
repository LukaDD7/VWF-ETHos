from __future__ import annotations

from pathlib import Path
import json
import os
import re
from typing import Any

import pandas as pd

from ..schemas import EvidenceItem, VariantContext, utc_now


AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


def compact_missense(hgvs_p: str | None) -> str | None:
    if not hgvs_p:
        return None
    match = re.fullmatch(
        r"p\.(" + "|".join(AA3_TO_1) + r")(\d+)(" + "|".join(AA3_TO_1) + r")",
        hgvs_p.strip(),
    )
    if not match:
        return None
    ref, position, alt = match.groups()
    return f"{AA3_TO_1[ref]}{position}{AA3_TO_1[alt]}"


class LocalComputationalPanelProvider:
    """Read-only adapter for already-computed AlphaGenome and Boltz/mechanism panels.

    The provider never launches AlphaGenome or Boltz. Missing rows remain explicit
    pending requests and are not converted into negative evidence.
    """

    name = "local_computational_panels"
    version = "v1"

    def __init__(self, repo_root: str | Path):
        self.repo_root = Path(repo_root)
        self.ag_path = self.repo_root / "results/07_VCF_AlphaGenome_Results.csv"
        self.ag = pd.read_csv(self.ag_path) if self.ag_path.exists() else pd.DataFrame()
        override = os.environ.get("VWF_COMPUTATIONAL_EVIDENCE", "").strip()
        defaults = [
            self.repo_root / "outputs/type1_panel_agent_20260828/server_results/type1_computational_evidence.csv",
            self.repo_root / "outputs/type2b_panel_agent_20260829/server_results/type2b_computational_evidence.csv",
        ]
        self.returned_paths = [Path(override)] if override else defaults
        returned_frames = [pd.read_csv(path) for path in self.returned_paths if path.exists()]
        self.returned = pd.concat(returned_frames, ignore_index=True, sort=False) if returned_frames else pd.DataFrame()
        self.mechanism_tool: Any | None = None
        mechanism_path = self.repo_root / "output/eval_v2_with_type2m_lof_hf/eval_v2_predictions.csv"
        self.mechanism_predictions = (
            pd.read_csv(mechanism_path)[
                ["aa_change", "pred_with_md", "confidence_with_md", "reasoning_with_md"]
            ]
            if mechanism_path.exists()
            else pd.DataFrame()
        )

    def collect(self, *, patient_id: str, variant: VariantContext) -> tuple[list[EvidenceItem], dict[str, Any]]:
        items: list[EvidenceItem] = []
        statuses: dict[str, Any] = {
            "variant_id": variant.source_row_id,
            "alphagenome": "not_applicable",
            "boltz": "not_applicable",
        }
        if variant.gene.upper() != "VWF":
            statuses["reason"] = f"gene={variant.gene}; VWF panels not applicable"
            return items, statuses

        returned_items = self._returned_items(patient_id, variant)
        items.extend(returned_items)
        returned_ag = any(item.source == "alphagenome_full_profile" for item in returned_items)
        returned_boltz = any(item.source in {"boltz2_type1_panel", "md_type1_panel", "boltz2_functional_panel", "md_targeted_panel"} for item in returned_items)

        ag_item = None if returned_ag else self._alphagenome_item(patient_id, variant)
        if returned_ag:
            statuses["alphagenome"] = "returned_full_profile_embedded"
        elif ag_item is not None:
            items.append(ag_item)
            statuses["alphagenome"] = "existing_result_matched"
        elif variant.alphagenome_request_status == "READY":
            statuses["alphagenome"] = "request_pending"
        else:
            statuses["alphagenome"] = variant.alphagenome_request_status or "not_modelable"

        compact = compact_missense(variant.hgvs_p)
        mechanism_item = None if returned_boltz else self._mechanism_item(patient_id, variant, compact)
        if returned_boltz:
            statuses["boltz"] = "returned_panel_embedded"
        elif mechanism_item is not None:
            items.append(mechanism_item)
            statuses["boltz"] = "existing_result_interpreted"
        elif variant.boltz_request_status == "READY":
            statuses["boltz"] = "request_pending"
        else:
            statuses["boltz"] = variant.boltz_request_status or "not_modelable"
        return items, statuses

    def _returned_items(self, patient_id: str, variant: VariantContext) -> list[EvidenceItem]:
        if self.returned.empty or "case_id" not in self.returned.columns:
            return []
        matches = self.returned[self.returned["case_id"].astype(str).eq(patient_id)]
        if len(matches) > 1:
            compact = compact_missense(variant.hgvs_p)
            coordinate = (
                f"{variant.hg38_chromosome}:{variant.hg38_position}:{variant.hg38_ref}>{variant.hg38_alt}"
                if all((variant.hg38_chromosome, variant.hg38_position, variant.hg38_ref, variant.hg38_alt))
                else ""
            )
            tokens = {value for value in (variant.source_row_id, variant.hgvs_c, variant.hgvs_p, compact, coordinate) if value}
            specific = matches[
                matches.apply(
                    lambda row: any(
                        token in str(row.get("source_record_id", ""))
                        or token == str(row.get("query", ""))
                        for token in tokens
                    ),
                    axis=1,
                )
            ]
            if not specific.empty:
                matches = specific
        items: list[EvidenceItem] = []
        for index, row in matches.iterrows():
            source = str(row.get("source", "")).strip()
            if not source:
                continue
            items.append(
                EvidenceItem(
                    source=source,
                    source_record_id=str(row.get("source_record_id") or f"{patient_id}:{variant.variant_index}:{source}:{index}"),
                    query=str(row.get("query") or variant.hgvs_c or variant.source_row_id),
                    retrieved_at=str(row.get("retrieved_at") or utc_now()),
                    source_version=str(row.get("source_version") or "returned-computational-evidence"),
                    supports=_list_value(row.get("supports")),
                    refutes=_list_value(row.get("refutes")),
                    conclusion=str(row.get("conclusion") or "Returned computational result available."),
                    confidence=_float(row.get("confidence")),
                    evidence_level="mechanism_model",
                    limitations=_list_value(row.get("limitations")),
                    raw_excerpt_locator=str(row.get("raw_excerpt_locator") or "|".join(map(str, self.returned_paths))),
                )
            )
        return items

    def _alphagenome_item(self, patient_id: str, variant: VariantContext) -> EvidenceItem | None:
        if self.ag.empty or not all((variant.hg38_chromosome, variant.hg38_position, variant.hg38_ref, variant.hg38_alt)):
            return None
        matches = self.ag[
            self.ag["hg38_chr"].astype(str).eq(str(variant.hg38_chromosome))
            & pd.to_numeric(self.ag["hg38_pos"], errors="coerce").eq(int(variant.hg38_position))
            & self.ag["ref"].astype(str).str.upper().eq(str(variant.hg38_ref).upper())
            & self.ag["alt"].astype(str).str.upper().eq(str(variant.hg38_alt).upper())
            & self.ag["Status"].astype(str).str.casefold().eq("success")
        ]
        if matches.empty:
            return None
        row = matches.iloc[0]
        rna = _float(row.get("AG_RNA_SEQ"))
        splice = _float(row.get("AG_SPLICE_SITES"))
        delta = _float(row.get("Delta_Max_Core"))
        supports = ["splice_effect_candidate"] if splice is not None and splice >= 0.5 else []
        return EvidenceItem(
            source="alphagenome_existing_panel",
            source_record_id=f"{patient_id}:{variant.variant_index}:alphagenome",
            query=f"{variant.hg38_chromosome}:{variant.hg38_position}:{variant.hg38_ref}>{variant.hg38_alt}",
            retrieved_at=utc_now(),
            source_version=self.ag_path.name,
            supports=supports,
            conclusion=f"Existing unsigned AlphaGenome maxima: RNA={rna}, splice-sites={splice}, core={delta}.",
            confidence=None,
            evidence_level="mechanism_model",
            limitations=[
                "These are absolute maximum deltas, not signed gene-level expression effects.",
                "A large score is a prioritization signal and is not a pathogenicity or subtype diagnosis.",
            ],
            raw_excerpt_locator=str(self.ag_path),
        )

    def _mechanism_item(
        self,
        patient_id: str,
        variant: VariantContext,
        compact: str | None,
    ) -> EvidenceItem | None:
        if compact is None or self.mechanism_predictions.empty:
            return None
        matches = self.mechanism_predictions[self.mechanism_predictions["aa_change"].astype(str).eq(compact)]
        if matches.empty:
            return None
        row = matches.iloc[0]
        predicted = str(row.get("pred_with_md", "uncertain"))
        confidence = _float(row.get("confidence_with_md"))
        reasoning = str(row.get("reasoning_with_md", ""))
        supports = [] if predicted == "uncertain" else [predicted]
        return EvidenceItem(
            source="boltz_mechanism_classifier",
            source_record_id=f"{patient_id}:{variant.variant_index}:boltz-mechanism",
            query=compact,
            retrieved_at=utc_now(),
            source_version="agentic_vwf_classifier+boltz2_vwd_functional_panel",
            supports=supports,
            conclusion=f"Prior AgenticVWFClassifier prediction: {predicted}" + (f"; {reasoning}" if reasoning else ""),
            confidence=confidence,
            evidence_level="mechanism_model",
            limitations=[
                "Research-only mechanism evidence; it cannot establish clinical pathogenicity.",
                "Boltz/AF3 confidence metrics are structural proxies and require assay or clinical confirmation.",
            ],
            raw_excerpt_locator="output/boltz2_vwd_functional_panel/evidence_matrix.csv",
        )


def _float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if pd.notna(number) else None


def _list_value(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    if text.startswith("["):
        try:
            parsed = json.loads(text)
            if isinstance(parsed, list):
                return [str(item) for item in parsed if str(item).strip()]
        except json.JSONDecodeError:
            pass
    return [item.strip() for item in text.split("|") if item.strip()]

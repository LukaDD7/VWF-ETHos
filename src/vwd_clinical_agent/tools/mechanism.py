from __future__ import annotations

import importlib.util
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pandas as pd

from .base import BaseBiomedicalTool, ToolRequest
from .fhir import FHIRResource, observation


AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

CLASSIFIER_FEATURE_COLUMNS = [
    "protein_pos",
    "ref_aa",
    "alt_aa",
    "domain",
    "ag_rna_delta",
    "ag_splice_delta",
    "af3_plddt_mean",
    "af3_plddt_min",
    "af3_pae_interface",
    "fb_binding_zscore",
    "heparan_zscore",
    "a3_collagen_zscore",
    "aim_release_score",
    "aim_sb_retained_z",
    "md_face_destab_score",
]


def _normalize_protein_change(value: str | None) -> str:
    """Return the compact one-letter form used by the repository matrices."""
    if not value:
        return ""
    compact = value.strip().removeprefix("p.").replace("VWF_", "")
    if len(compact) >= 5 and compact[:3].istitle() and compact[-3:].istitle():
        ref = AA3_TO_AA1.get(compact[:3], compact[:3])
        alt = AA3_TO_AA1.get(compact[-3:], compact[-3:])
        middle = compact[3:-3]
        if middle.isdigit():
            return f"{ref}{middle}{alt}"
    return compact


def _load_feature_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _load_classifier(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("agentic_vwf_classifier_embedded", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import classifier from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _matches(row: pd.Series, hgvs_c: str, protein_change: str) -> bool:
    row_protein = _normalize_protein_change(str(row.get("aa_change", "")))
    if protein_change and row_protein == protein_change:
        return True
    row_hgvs_c = str(row.get("hgvs_c", ""))
    return bool(hgvs_c and row_hgvs_c and hgvs_c.split(":")[-1] == row_hgvs_c.split(":")[-1])


class MechanismClassifierTool(BaseBiomedicalTool):
    """Read-only adapter for the repository AgenticVWFClassifier feature table."""

    name = "mechanism_classifier"
    version = "agentic-vwf-v2-adapter-v0"
    endpoint = "local://agentic-vwf-classifier"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        matrix_path = Path(str(request.parameters.get("matrix_path", "")))
        if not matrix_path.exists():
            return [
                observation(
                    observation_id="mechanism-classifier-not-configured",
                    patient_id=request.patient_id,
                    display="Agentic VWF mechanism classifier",
                    value="not_configured",
                    based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                )
            ], "disabled"

        classifier_path = Path(
            str(request.parameters.get("classifier_path", request.parameters.get("default_classifier_path", "")))
        )
        if not classifier_path.exists():
            return [
                observation(
                    observation_id="mechanism-classifier-not-configured",
                    patient_id=request.patient_id,
                    display="Agentic VWF mechanism classifier",
                    value="classifier_not_found",
                    based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                )
            ], "disabled"

        hgvs_c = str(request.parameters.get("hgvs_c", ""))
        hgvs_p = str(request.parameters.get("hgvs_p", ""))
        protein_change = _normalize_protein_change(hgvs_p)
        table = _load_feature_table(matrix_path)
        if not table.empty and "aa_change" not in table.columns and "hgvs_c" not in table.columns:
            return [
                observation(
                    observation_id="mechanism-classifier-invalid-table",
                    patient_id=request.patient_id,
                    display="Agentic VWF mechanism classifier",
                    value="feature_table_missing_variant_key",
                    based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                )
            ], "error"

        matches = table.apply(lambda row: _matches(row, hgvs_c, protein_change), axis=1) if len(table) else pd.Series(dtype=bool)
        if not bool(matches.any()):
            return [
                observation(
                    observation_id="mechanism-classifier-not-found",
                    patient_id=request.patient_id,
                    display="Agentic VWF mechanism classifier",
                    value="variant_absent_from_feature_table",
                    based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                    components=[
                        {"code": {"text": "hgvs_p"}, "valueString": hgvs_p},
                        {"code": {"text": "matrix_path"}, "valueString": str(matrix_path)},
                    ],
                )
            ], "not_found"

        row = table.loc[matches.idxmax()]
        module = _load_classifier(classifier_path)
        features = {column: row.get(column) for column in CLASSIFIER_FEATURE_COLUMNS if column in row.index}
        # The classifier is not fitted on the queried row. This is deliberate:
        # labels and row-specific adaptive fitting must not leak into inference.
        result = module.AgenticVWFClassifier().classify(features)

        components = [
            {"code": {"text": "confidence"}, "valueQuantity": {"value": float(result.confidence)}},
            {"code": {"text": "alternatives"}, "valueString": result.alternatives},
            {"code": {"text": "reasoning"}, "valueString": result.reasoning},
            {"code": {"text": "matrix_path"}, "valueString": str(matrix_path)},
            {"code": {"text": "classifier_path"}, "valueString": str(classifier_path)},
            {"code": {"text": "matched_variant"}, "valueString": str(row.get("aa_change", hgvs_p))},
            {"code": {"text": "clinical_override"}, "valueBoolean": False},
            {"code": {"text": "expert_review_required"}, "valueBoolean": True},
        ]
        for column in CLASSIFIER_FEATURE_COLUMNS:
            if column in row.index and pd.notna(row[column]):
                components.append({"code": {"text": column}, "valueString": str(row[column])})

        resource = observation(
            observation_id=f"mechanism-{protein_change or hgvs_c or request.variant_id}",
            patient_id=request.patient_id,
            display="Agentic VWF mechanism classifier",
            value=result.main_subtype,
            based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
            components=components,
        )
        provenance_payload = {
            "prediction": asdict(result),
            "input_features": features,
            "matrix_path": str(matrix_path),
            "classifier_path": str(classifier_path),
            "label_columns_used": [],
        }
        return [resource, self.provenance_for(f"Observation/{resource.id}", request, provenance_payload)], "success"

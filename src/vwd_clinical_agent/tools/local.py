from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from .base import BaseBiomedicalTool, ToolRequest
from .fhir import FHIRResource, observation, operation_outcome


class LocalScoreStore(BaseBiomedicalTool):
    """Read-only lookup in a locally downloaded versioned score table."""

    name = "local_score_store"
    version = "snapshot-v1"
    endpoint = "local://score-store"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        path = Path(str(request.parameters.get("path", "")))
        if not path.exists():
            return [
                operation_outcome(
                    "warning",
                    "not-configured",
                    f"Local score snapshot not found: {path}. Download a licensed/versioned snapshot before enabling.",
                )
            ], "disabled"
        score_name = str(request.parameters.get("score_name", "score"))
        hgvs = str(request.parameters.get("hgvs_c", ""))
        chrom = str(request.parameters.get("chrom", ""))
        pos = request.parameters.get("pos")
        ref = str(request.parameters.get("ref", ""))
        alt = str(request.parameters.get("alt", ""))
        table = pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",")
        required = {score_name}
        if not required.issubset(table.columns):
            raise ValueError(f"Snapshot missing required columns: {sorted(required - set(table.columns))}")
        matches: pd.DataFrame | None = None
        for columns in (("hgvs_c",), ("chrom", "pos", "ref", "alt")):
            if set(columns).issubset(table.columns):
                mask = pd.Series(True, index=table.index)
                if columns == ("hgvs_c",):
                    mask &= table["hgvs_c"].astype(str) == hgvs
                else:
                    mask &= table["chrom"].astype(str).str.removeprefix("chr") == chrom.removeprefix("chr")
                    mask &= table["pos"].astype(int) == int(pos)
                    mask &= table["ref"].astype(str) == ref
                    mask &= table["alt"].astype(str) == alt
                matches = table.loc[mask]
                if not matches.empty:
                    break
        if matches is None or matches.empty:
            return [operation_outcome("information", "not-found", "Variant absent from local score snapshot")], "not_found"
        row = matches.iloc[0].to_dict()
        resource = observation(
            observation_id=f"local-{score_name}-{hgvs}-{chrom}-{pos}-{ref}-{alt}",
            patient_id=request.patient_id,
            display=f"Local {score_name} annotation",
            value=row[score_name],
            based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
            components=[
                _local_component("snapshot_path", path),
                _local_component("snapshot_sha256", request.parameters.get("snapshot_sha256", "not_recorded")),
            ],
        )
        return [resource, self.provenance_for(f"Observation/{resource.id}", request, row)], "success"


def _local_component(code: str, value: Any) -> dict[str, Any]:
    return {"code": {"text": code}, "valueString": str(value)}


class RepoMechanismTool(BaseBiomedicalTool):
    """Read-only adapter for existing AlphaFold/Boltz/MD/FoldX artifacts."""

    name = "repo_mechanism"
    version = "adapter-only-v0"
    endpoint = "local://repo-artifacts"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        artifacts = request.parameters.get("artifact_paths", [])
        if not artifacts:
            return [operation_outcome("warning", "not-configured", "No mechanism artifacts were provided")], "disabled"
        resources: list[FHIRResource] = []
        for artifact in artifacts:
            path = Path(str(artifact))
            if not path.exists():
                resources.append(operation_outcome("warning", "not-found", f"Artifact not found: {path}"))
                continue
            resource = observation(
                observation_id=f"mechanism-{path.stem}-{path.suffix}",
                patient_id=request.patient_id,
                display="Repository mechanism-model artifact",
                value=str(path),
                based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                components=[_local_component("read_only", True), _local_component("clinical_override", False)],
            )
            resources.append(resource)
            resources.append(self.provenance_for(f"Observation/{resource.id}", request, {"path": str(path)}))
        if not any(resource.resourceType == "Observation" for resource in resources):
            return resources, "not_found"
        return resources, "success"


class LocalGuidelineTool(BaseBiomedicalTool):
    name = "local_guideline"
    version = "clingen-acmg-json-v0"
    endpoint = "local://guideline-snapshot"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        path = Path(str(request.parameters.get("path", "")))
        if not path.exists():
            return [operation_outcome("warning", "not-configured", "Approved local guideline rule snapshot is missing")], "disabled"
        import json

        rules = json.loads(path.read_text(encoding="utf-8"))
        resource = observation(
            observation_id=f"guideline-{path.stem}",
            patient_id=request.patient_id,
            display="Local ClinGen/ACMG rule snapshot",
            value=path.name,
            based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
            components=[_local_component("rule_count", len(rules.get("rules", [])))],
        )
        return [resource, self.provenance_for(f"Observation/{resource.id}", request, rules)], "success"


class LocalClinGenSnapshotTool(BaseBiomedicalTool):
    name = "local_clingen_snapshot"
    version = "erepo-json-snapshot-v1"
    endpoint = "local://clingen-vwf-snapshot"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        path = Path(str(request.parameters.get("path", "")))
        if not path.exists():
            return [operation_outcome("warning", "not-configured", f"ClinGen snapshot missing: {path}")], "disabled"
        import json

        payload = json.loads(path.read_text(encoding="utf-8"))
        hgvs_c = str(request.parameters.get("hgvs_c", "")).split(":")[-1]
        matches = [
            row
            for row in payload.get("variantInterpretations", [])
            if any(hgvs_c in str(item) for item in row.get("hgvs", []))
        ]
        if not matches:
            return [operation_outcome("information", "not-found", "Variant absent from local ClinGen VWF snapshot")], "not_found"
        resources: list[FHIRResource] = []
        for index, row in enumerate(matches):
            guideline = (row.get("guidelines") or [{}])[0]
            outcome = guideline.get("outcome", {}).get("label", "not stated")
            resource = observation(
                observation_id=f"local-clingen-{index}-{row.get('variationId', index)}",
                patient_id=request.patient_id,
                display="Local ClinGen VWF classification snapshot",
                value=outcome,
                based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                components=[
                    _local_component("snapshot_path", path),
                    _local_component("caid", row.get("caid")),
                    _local_component("condition", row.get("condition", {}).get("label")),
                    _local_component("expert_panel", guideline.get("label")),
                ],
            )
            resources.append(resource)
            resources.append(self.provenance_for(f"Observation/{resource.id}", request, row))
        return resources, "success"


class LocalClinVarSnapshotTool(BaseBiomedicalTool):
    name = "local_clinvar_snapshot"
    version = "eutils-summary-snapshot-v1"
    endpoint = "local://clinvar-vwf-snapshot"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        path = Path(str(request.parameters.get("path", "")))
        if not path.exists():
            return [operation_outcome("warning", "not-configured", f"ClinVar snapshot missing: {path}")], "disabled"
        import json

        payload = json.loads(path.read_text(encoding="utf-8"))
        result = payload.get("result", {})
        hgvs_c = str(request.parameters.get("hgvs_c", "")).split(":")[-1]
        protein = str(request.parameters.get("hgvs_p", "")).removeprefix("p.")
        matches = []
        for uid in result.get("uids", []):
            row = result.get(uid, {})
            names = [row.get("title", "")]
            names.extend(item.get("variation_name", "") for item in row.get("variation_set", []))
            if hgvs_c and any(hgvs_c in name for name in names):
                matches.append(row)
            elif protein and any(protein in name for name in names):
                matches.append(row)
        if not matches:
            return [operation_outcome("information", "not-found", "Variant absent from local ClinVar VWF snapshot")], "not_found"
        resources: list[FHIRResource] = []
        for row in matches[:10]:
            classification = row.get("germline_classification", {})
            resource = observation(
                observation_id=f"local-clinvar-{row['uid']}",
                patient_id=request.patient_id,
                display="Local ClinVar VWF classification snapshot",
                value=classification.get("description", "not stated"),
                based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                components=[
                    _local_component("snapshot_path", path),
                    _local_component("accession", row.get("accession")),
                    _local_component("review_status", classification.get("review_status")),
                    _local_component("last_evaluated", classification.get("last_evaluated")),
                ],
            )
            resources.append(resource)
            resources.append(self.provenance_for(f"Observation/{resource.id}", request, row))
        return resources, "success"

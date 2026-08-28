from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from .base import BaseBiomedicalTool, ToolRegistry, ToolRequest, ToolResponse
from .fhir import FHIRBundle, FHIRResource, operation_outcome
from .local import (
    LocalClinGenSnapshotTool,
    LocalClinVarSnapshotTool,
    LocalGuidelineTool,
    LocalScoreStore,
    RepoMechanismTool,
)
from .full_text import PubMedFullTextSearchTool
from .literature_extraction import LiteraturePhenotypeExtractor
from .variant_context import VWFDomainAnnotator
from .online import (
    ClinGenERepoProvider,
    EnsemblVariantNormalizer,
    GnomADProvider,
    HGMDProvider,
    OpenCravatAnnotator,
    PubMedSearchProvider,
)


class ToolCallRecord(BaseModel):
    tool: str
    operation: str
    status: str
    retry_count: int
    latency_ms: int
    cache_hit: bool
    diagnostics: list[str] = Field(default_factory=list)


@dataclass
class EvidenceToolMatrixResult:
    bundle: FHIRBundle
    calls: list[ToolCallRecord]
    normalized: dict[str, Any]
    authorized: bool


class EvidenceToolMatrix:
    """Dependency-aware, FHIR-native code-as-search executor."""

    def __init__(self, cache_dir: str | None = None, snapshot_dir: str | None = None):
        self.registry = ToolRegistry(
            [
                EnsemblVariantNormalizer(),
                GnomADProvider(),
                OpenCravatAnnotator(),
                ClinGenERepoProvider(),
                PubMedSearchProvider(),
                PubMedFullTextSearchTool(),
                LiteraturePhenotypeExtractor(),
                VWFDomainAnnotator(),
                HGMDProvider(),
                LocalScoreStore(),
                LocalClinGenSnapshotTool(),
                LocalClinVarSnapshotTool(),
                RepoMechanismTool(),
                LocalGuidelineTool(),
            ],
            cache_dir=cache_dir,
        )
        self.snapshot_dir = Path(snapshot_dir) if snapshot_dir else None
        self.include_pubmed_full_text = False

    def run(
        self,
        *,
        patient_id: str,
        variant_id: str,
        hgvs_c: str,
        hgvs_p: str | None,
        second_level_bundle: FHIRBundle,
        allow_research_bypass: bool = False,
        allow_pre_second_level: bool = False,
        local_parameters: dict[str, Any] | None = None,
    ) -> EvidenceToolMatrixResult:
        local_parameters = {**(local_parameters or {}), "include_pubmed_full_text": self.include_pubmed_full_text}
        if not self._second_level_complete(second_level_bundle) and not allow_research_bypass and not allow_pre_second_level:
            bundle = FHIRBundle.of(
                [
                    operation_outcome(
                        "error",
                        "business-rule",
                        "Online variant tools are gated until real second-level observations or documented unavailable tests are present.",
                    )
                ]
            )
            return EvidenceToolMatrixResult(bundle, [], {}, False)

        calls: list[ToolCallRecord] = []
        output = FHIRBundle(type="collection")

        def invoke(tool: str, request: ToolRequest) -> ToolResponse:
            response = self.registry.invoke(tool, request)
            calls.append(
                ToolCallRecord(
                    tool=response.tool,
                    operation=response.operation,
                    status=response.status,
                    retry_count=response.retry_count,
                    latency_ms=response.latency_ms,
                    cache_hit=response.cache_hit,
                    diagnostics=response.diagnostics,
                )
            )
            for entry in response.output.entry:
                output.entry.append(entry)
            return response

        normalizer_request = ToolRequest(
            operation="normalize_variant", patient_id=patient_id, variant_id=variant_id, parameters={"hgvs_c": hgvs_c}
        )
        normalizer = invoke("ensembl_variant_recoder", normalizer_request)
        normalized = self._normalized_parameters(normalizer)
        normalized["hgvs_c"] = hgvs_c
        normalized["hgvs_p"] = hgvs_p
        invoke(
            "vwf_domain_annotator",
            ToolRequest(
                operation="annotate_domain",
                patient_id=patient_id,
                variant_id=variant_id,
                parameters={"hgvs_p": hgvs_p or ""},
            ),
        )

        coordinate_ready = all(key in normalized for key in ("chrom", "pos", "ref", "alt"))
        if normalizer.status == "error" or not coordinate_ready:
            output.add(
                operation_outcome(
                    "error",
                    "dependency-failed",
                    "Variant normalization failed; coordinate-dependent tools were not called and this is not benign evidence.",
                )
            )
            calls.append(
                ToolCallRecord(
                    tool="coordinate_dependent_tools",
                    operation="dependency_guard",
                    status="error",
                    retry_count=normalizer.retry_count,
                    latency_ms=normalizer.latency_ms,
                    cache_hit=normalizer.cache_hit,
                    diagnostics=normalizer.diagnostics,
                )
            )
        else:
            invoke("gnomad_graphql", ToolRequest(operation="population_frequency", patient_id=patient_id, variant_id=variant_id, parameters=normalized))
            open_cravat_response = invoke(
                "open_cravat",
                ToolRequest(operation="variant_annotation", patient_id=patient_id, variant_id=variant_id, parameters=normalized),
            )
            configured_snapshot = (local_parameters or {}).get("snapshot_dir") or self.snapshot_dir
            snapshot_dir = Path(str(configured_snapshot or ""))
            if snapshot_dir.exists() and open_cravat_response.status == "error":
                invoke(
                    "local_clinvar_snapshot",
                    ToolRequest(
                        operation="local_clinical_classification",
                        patient_id=patient_id,
                        variant_id=variant_id,
                        parameters={
                            "path": snapshot_dir / "clinvar_vwf_esummary.json",
                            "hgvs_c": hgvs_c,
                            "hgvs_p": hgvs_p,
                        },
                    ),
                )
        clingen_response = invoke(
            "clingen_erepo",
            ToolRequest(
                operation="expert_classification",
                patient_id=patient_id,
                variant_id=variant_id,
                parameters={"hgvs_c": hgvs_c, "hgvsg": normalized.get("hgvsg")},
            ),
        )
        pubmed_response = invoke(
            "pubmed_eutils",
            ToolRequest(
                operation="literature_search",
                patient_id=patient_id,
                variant_id=variant_id,
                parameters={
                    "query": self._pubmed_query(hgvs_c, hgvs_p, normalized.get("rsid")),
                    "retmax": 5,
                    "include_full_text": bool((local_parameters or {}).get("include_pubmed_full_text", False)),
                },
            ),
        )
        if any(
            any("pub-full-text" in extension.get("url", "") for extension in resource.get("extension", []))
            for resource in pubmed_response.output.resources("DocumentReference")
        ):
            full_text_bundle = FHIRBundle.of(
                [
                    FHIRResource.model_validate(resource)
                    for resource in pubmed_response.output.resources("DocumentReference")
                ]
            )
            invoke(
                "pubmed_full_text_search",
                ToolRequest(
                    operation="full_text_search",
                    patient_id=patient_id,
                    variant_id=variant_id,
                    input=full_text_bundle,
                ),
            )
            invoke(
                "literature_phenotype_extractor",
                ToolRequest(
                    operation="extract_phenotype",
                    patient_id=patient_id,
                    variant_id=variant_id,
                    input=full_text_bundle,
                ),
            )
        invoke("hgmd", ToolRequest(operation="licensed_classification", patient_id=patient_id, variant_id=variant_id))
        configured_snapshot = (local_parameters or {}).get("snapshot_dir") or self.snapshot_dir
        snapshot_dir = Path(str(configured_snapshot or ""))
        if snapshot_dir.exists() and clingen_response.status == "error":
            invoke(
                "local_clingen_snapshot",
                ToolRequest(
                    operation="local_expert_classification",
                    patient_id=patient_id,
                    variant_id=variant_id,
                    parameters={
                        "path": snapshot_dir / "clingen_vwf_erepo.json",
                        "hgvs_c": hgvs_c,
                    },
                ),
            )
        # ClinVar online annotation is primary. The local ClinVar snapshot is only
        # an offline fallback when OpenCRAVAT itself errors; not_found is not a
        # reason to treat the smaller local snapshot as authoritative.

        for score_name, path in (local_parameters or {}).get("local_scores", {}).items():
            params = {**normalized, "score_name": score_name, "path": path}
            invoke("local_score_store", ToolRequest(operation="local_score", patient_id=patient_id, variant_id=variant_id, parameters=params))
        if local_parameters and local_parameters.get("guideline_path"):
            invoke(
                "local_guideline",
                ToolRequest(
                    operation="apply_guideline",
                    patient_id=patient_id,
                    variant_id=variant_id,
                    parameters={"path": local_parameters["guideline_path"]},
                ),
            )
        if local_parameters and local_parameters.get("artifact_paths"):
            invoke(
                "repo_mechanism",
                ToolRequest(
                    operation="mechanism_artifacts",
                    patient_id=patient_id,
                    variant_id=variant_id,
                    parameters={"artifact_paths": local_parameters["artifact_paths"]},
                ),
            )

        return EvidenceToolMatrixResult(output, calls, normalized, True)

    def _second_level_complete(self, bundle: FHIRBundle) -> bool:
        reports = bundle.resources("DiagnosticReport")
        if not reports:
            return False
        return any(report.get("status") == "final" for report in reports)

    def _normalized_parameters(self, response: ToolResponse) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for resource in response.output.resources("Observation"):
            if resource.get("code", {}).get("text") != "Normalized VWF variant":
                continue
            for component in resource.get("component", []):
                text = component.get("code", {}).get("text")
                if text == "chrom-pos-ref-alt":
                    chrom, pos, ref, alt = str(component.get("valueString", "---")).split("-", 3)
                    result.update({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt})
                elif text == "rsid":
                    result["rsid"] = str(component.get("valueString", "")).split(";")[0].strip()
                elif text == "hgvsg":
                    result["hgvsg"] = str(component.get("valueString", "")).split(";")[0].strip()
        return result

    def _pubmed_query(self, hgvs_c: str, hgvs_p: str | None, rsid: str | None = None) -> str:
        protein = hgvs_p.removeprefix("p.") if hgvs_p else ""
        cdna = hgvs_c.split(":")[-1] if hgvs_c else ""
        terms = [term for term in (hgvs_p, protein, cdna, rsid) if term]
        return '("von Willebrand factor"[Text Word] OR VWF[Text Word]) AND (' + " OR ".join(f'"{term}"[Text Word]' for term in terms) + ")"

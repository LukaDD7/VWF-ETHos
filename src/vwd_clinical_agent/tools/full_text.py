from __future__ import annotations

import re
from typing import Any

from .base import BaseBiomedicalTool, ToolRequest
from .context_utils import best_match_near_variant, build_contextual_excerpt
from .fhir import FHIRBundle, FHIRResource, observation, operation_outcome


DEFAULT_TERMS = [
    "ristocetin",
    "RIPA",
    "multimer",
    "collagen binding",
    "factor VIII binding",
    "desmopressin",
    "DDAVP",
    "clearance",
    "type 2A",
    "type 2B",
    "type 2M",
    "type 2N",
]


class PubMedFullTextSearchTool(BaseBiomedicalTool):
    """Code-as-search over PMC full text already stored in FHIR DocumentReferences."""

    name = "pubmed_full_text_search"
    version = "pmc-local-index-v1"
    endpoint = "local://fhir-document-reference-full-text"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        documents = [
            resource
            for resource in request.input.resources("DocumentReference")
            if any("pub-full-text" in extension.get("url", "") for extension in resource.get("extension", []))
        ]
        if not documents:
            return [operation_outcome("information", "not-found", "No PMC full-text documents were available")], "not_found"

        terms = [str(term) for term in request.parameters.get("terms", DEFAULT_TERMS)]
        max_excerpts = int(request.parameters.get("max_excerpts", 3))
        context_before_chars = int(request.parameters.get("context_before_chars", 600))
        context_after_chars = int(request.parameters.get("context_after_chars", 900))
        variant_link_radius = int(request.parameters.get("variant_link_radius", 1500))
        variant_terms = [str(term) for term in request.parameters.get("variant_terms", []) if term]
        resources: list[FHIRResource] = []
        for document in documents:
            full_text = next(
                extension.get("valueString", "")
                for extension in document.get("extension", [])
                if "pub-full-text" in extension.get("url", "")
            )
            if not full_text:
                continue
            document_variant_specific = next(
                (
                    extension.get("valueBoolean", False)
                    for extension in document.get("extension", [])
                    if "pub-variant-specific" in extension.get("url", "")
                ),
                False,
            )
            excerpts = self._extract_excerpts(
                full_text,
                terms,
                max_excerpts,
                variant_terms=variant_terms,
                context_before_chars=context_before_chars,
                context_after_chars=context_after_chars,
                variant_link_radius=variant_link_radius,
            )
            if not excerpts:
                continue
            for index, excerpt in enumerate(excerpts):
                resource = observation(
                    observation_id=f"pmc-excerpt-{document['id']}-{index}",
                    patient_id=request.patient_id,
                    display="PubMed Central full-text evidence excerpt",
                    value=excerpt["text"],
                    based_on=[f"DocumentReference/{document['id']}"],
                    components=[
                        {"code": {"text": "search_term"}, "valueString": excerpt["term"]},
                        {"code": {"text": "source_document"}, "valueString": document["id"]},
                        {"code": {"text": "context_before"}, "valueString": excerpt["before"]},
                        {"code": {"text": "context_after"}, "valueString": excerpt["after"]},
                        {"code": {"text": "context_chars"}, "valueInteger": len(excerpt["text"])},
                        {"code": {"text": "nearest_variant_term"}, "valueString": excerpt["nearest_variant_term"] or ""},
                        {
                            "code": {"text": "nearest_variant_distance"},
                            "valueInteger": excerpt["nearest_variant_distance"] if excerpt["nearest_variant_distance"] is not None else -1,
                        },
                        {"code": {"text": "variant_linked"}, "valueBoolean": excerpt["variant_linked"]},
                        {"code": {"text": "document_variant_specific"}, "valueBoolean": document_variant_specific},
                    ],
                )
                resources.append(resource)
                resources.append(
                    self.provenance_for(
                        f"Observation/{resource.id}",
                        request,
                        {"document_id": document["id"], "term": excerpt["term"], "excerpt": excerpt["text"]},
                    )
                )
        if not resources:
            return [operation_outcome("information", "not-found", "No full-text excerpts matched the search terms")], "not_found"
        return resources, "success"

    @staticmethod
    def _extract_excerpts(
        text: str,
        terms: list[str],
        max_excerpts: int,
        *,
        variant_terms: list[str],
        context_before_chars: int,
        context_after_chars: int,
        variant_link_radius: int,
    ) -> list[dict[str, Any]]:
        excerpts: list[dict[str, Any]] = []
        seen_terms: set[str] = set()
        for term in terms:
            if term.casefold() in seen_terms:
                continue
            position = best_match_near_variant(text=text, term=term, variant_terms=variant_terms)
            if position is None:
                continue
            start, end = position
            contextual = build_contextual_excerpt(
                text=text,
                start=start,
                end=end,
                term=term,
                variant_terms=variant_terms,
                context_before_chars=context_before_chars,
                context_after_chars=context_after_chars,
                variant_link_radius=variant_link_radius,
            )
            excerpts.append(
                {
                    "term": term,
                    "text": contextual.text,
                    "before": contextual.before,
                    "after": contextual.after,
                    "nearest_variant_term": contextual.nearest_variant_term,
                    "nearest_variant_distance": contextual.nearest_variant_distance,
                    "variant_linked": contextual.variant_linked,
                }
            )
            seen_terms.add(term.casefold())
            if len(excerpts) >= max_excerpts:
                break
        return excerpts

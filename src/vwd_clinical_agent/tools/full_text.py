from __future__ import annotations

import re
from typing import Any

from .base import BaseBiomedicalTool, ToolRequest
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
        resources: list[FHIRResource] = []
        for document in documents:
            full_text = next(
                extension.get("valueString", "")
                for extension in document.get("extension", [])
                if "pub-full-text" in extension.get("url", "")
            )
            if not full_text:
                continue
            excerpts = self._extract_excerpts(full_text, terms, max_excerpts)
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
    def _extract_excerpts(text: str, terms: list[str], max_excerpts: int) -> list[dict[str, str]]:
        excerpts: list[dict[str, str]] = []
        seen_terms: set[str] = set()
        for term in terms:
            if term.casefold() in seen_terms:
                continue
            pattern = re.compile(re.escape(term), flags=re.IGNORECASE)
            for match in pattern.finditer(text):
                start = max(0, match.start() - 350)
                end = min(len(text), match.end() + 650)
                clean = re.sub(r"\s+", " ", text[start:end]).strip()
                excerpts.append({"term": term, "text": clean})
                seen_terms.add(term.casefold())
                if len(excerpts) >= max_excerpts:
                    return excerpts
                break
        return excerpts

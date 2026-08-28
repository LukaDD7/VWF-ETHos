from __future__ import annotations

import re
from typing import Any

from .base import BaseBiomedicalTool, ToolRequest
from .fhir import FHIRResource, observation, operation_outcome


PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("VWF:Ag", re.compile(r"(?:VWF(?::|\s)?Ag|von Willebrand factor antigen)[^.\n]{0,120}?(\d+(?:\.\d+)?)\s*(?:IU/dL|%|U/dL)", re.I)),
    ("VWF:RCo", re.compile(r"(?:VWF(?::|\s)?RCo|ristocetin cofactor)[^.\n]{0,120}?(\d+(?:\.\d+)?)\s*(?:IU/dL|%|U/dL)", re.I)),
    ("VWF:CB", re.compile(r"(?:VWF(?::|\s)?CB|collagen[- ]binding)[^.\n]{0,120}?(\d+(?:\.\d+)?)\s*(?:IU/dL|%|U/dL)", re.I)),
    ("FVIII:C", re.compile(r"(?:FVIII(?::|\s)?C|factor VIII(?: activity)?)[^.\n]{0,120}?(\d+(?:\.\d+)?)\s*(?:IU/dL|%|U/dL)", re.I)),
    ("RIPA", re.compile(r"(?:ristocetin[- ]induced platelet aggregation|RIPA)[^.\n]{0,180}?(normal|increased|enhanced|reduced|absent|abnormal)", re.I)),
    ("VWF multimer", re.compile(r"(?:multimer|multimers)[^.\n]{0,180}?(normal|absent|reduced|abnormal|loss[^.\n]{0,40}high[- ]molecular[- ]weight)", re.I)),
    ("Secretion", re.compile(r"(?:secretion|secreted)[^.\n]{0,180}?(normal|unaffected|reduced|impaired|no significant[^.\n]{0,30}effect)", re.I)),
    ("Dominant-negative effect", re.compile(r"(?:dominant[- ]negative)[^.\n]{0,180}?(absent|not[^.\n]{0,30}observed|present|observed)", re.I)),
    ("DDAVP response", re.compile(r"(?:DDAVP|desmopressin)[^.\n]{0,220}?(normal|partial|increased|reduced|no response|clearance)", re.I)),
    ("ABO effect", re.compile(r"\bABO\b[^.\n]{0,220}?(blood group|genotype|clearance|level|expression)", re.I)),
]


class LiteraturePhenotypeExtractor(BaseBiomedicalTool):
    name = "literature_phenotype_extractor"
    version = "regex-v1"
    endpoint = "local://full-text-phenotype-extractor"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        documents = [
            resource
            for resource in request.input.resources("DocumentReference")
            if any("pub-full-text" in extension.get("url", "") for extension in resource.get("extension", []))
        ]
        if not documents:
            return [operation_outcome("information", "not-found", "No PMC full text was available for phenotype extraction")], "not_found"
        resources: list[FHIRResource] = []
        for document in documents:
            text = next(
                extension.get("valueString", "")
                for extension in document.get("extension", [])
                if "pub-full-text" in extension.get("url", "")
            )
            for label, pattern in PATTERNS:
                matches = list(pattern.finditer(text))[:2]
                for index, match in enumerate(matches):
                    start = max(0, match.start() - 220)
                    end = min(len(text), match.end() + 420)
                    excerpt = re.sub(r"\s+", " ", text[start:end]).strip()
                    value = next((group for group in match.groups() if group is not None), "")
                    resource = observation(
                        observation_id=f"lit-{label.lower().replace(':','-').replace(' ','-')}-{document['id']}-{index}",
                        patient_id=request.patient_id,
                        display=f"Literature {label} evidence",
                        value=value,
                        based_on=[f"DocumentReference/{document['id']}"],
                        components=[
                            {"code": {"text": "excerpt"}, "valueString": excerpt},
                            {"code": {"text": "source_document"}, "valueString": document["id"]},
                        ],
                    )
                    resources.append(resource)
                    resources.append(self.provenance_for(f"Observation/{resource.id}", request, {"document": document["id"], "label": label, "excerpt": excerpt}))
        if not resources:
            return [operation_outcome("information", "not-found", "No phenotype statements matched the extraction patterns")], "not_found"
        return resources, "success"

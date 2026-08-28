from __future__ import annotations

from hashlib import sha256
from typing import Any
import json
import os
from urllib.parse import quote
import time

import httpx
from dotenv import load_dotenv


load_dotenv()

from .base import BaseBiomedicalTool, ToolRequest
from .fhir import FHIRResource, observation, operation_outcome


class HTTPMixin:
    timeout = 60.0
    max_retries = 3

    def request_json(
        self,
        method: str,
        url: str,
        *,
        params: dict[str, Any] | None = None,
        json_payload: dict[str, Any] | None = None,
        headers: dict[str, str] | None = None,
    ) -> tuple[int, dict[str, Any] | list[Any], int]:
        retries = 0
        response: httpx.Response
        with httpx.Client(timeout=self.timeout) as client:
            while True:
                try:
                    response = client.request(method, url, params=params, json=json_payload, headers=headers)
                except (httpx.TimeoutException, httpx.TransportError):
                    if retries >= self.max_retries:
                        raise
                    retries += 1
                    time.sleep(min(2**retries, 8))
                    continue
                if response.status_code not in {429, 500, 502, 503, 504} or retries >= self.max_retries:
                    break
                delay = float(response.headers.get("retry-after", str(min(2**retries, 8))))
                time.sleep(delay)
                retries += 1
        try:
            payload: dict[str, Any] | list[Any] = response.json()
        except ValueError:
            payload = {"raw": response.text}
        return response.status_code, payload, retries


def _component(code: str, value: Any, display: str | None = None) -> dict[str, Any]:
    return {
        "code": {"text": display or code},
        "valueString" if isinstance(value, str) else "valueQuantity": (
            value if isinstance(value, str) else {"value": value}
        ),
    }


class EnsemblVariantNormalizer(HTTPMixin, BaseBiomedicalTool):
    name = "ensembl_variant_recoder"
    version = "rest-ensembl-v1"
    endpoint = "https://rest.ensembl.org/variant_recoder/human"
    default_transcript = "NM_000552.5"
    timeout = 90.0

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        hgvs_c = str(request.parameters.get("hgvs_c", "")).strip()
        if not hgvs_c:
            raise ValueError("hgvs_c is required")
        if ":" not in hgvs_c:
            hgvs_c = f"{self.default_transcript}:{hgvs_c}"
        status, payload, retries = self.request_json(
            "GET",
            f"{self.endpoint}/{quote(hgvs_c, safe='')}",
            headers={"Content-Type": "application/json"},
        )
        if status != 200:
            raise RuntimeError(f"Ensembl returned HTTP {status}: {payload}")
        if not isinstance(payload, list) or not payload:
            return [operation_outcome("information", "no-match", "Ensembl returned no normalized variant")], "not_found"
        first = payload[0]
        values = list(first.values())
        normalized = values[0] if values and isinstance(values[0], dict) else {}
        spdi_items = normalized.get("spdi", [])
        spdi = spdi_items[0] if spdi_items else ""
        ref = alt = ""
        pos = None
        chrom = ""
        if spdi and ":" in spdi:
            chrom_num, raw_pos, ref, alt = spdi.split(":", 3)
            chrom = f"chr{chrom_num.replace('NC_', '').replace('000012.12', '12')}"
            if chrom == "chrNC_000012.12":
                chrom = "chr12"
            pos = int(raw_pos) + 1
        resource = observation(
            observation_id=f"variant-normalized-{sha256(hgvs_c.encode()).hexdigest()[:20]}",
            patient_id=request.patient_id,
            display="Normalized VWF variant",
            value=hgvs_c,
            based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
            components=[
                _component("hgvsg", "; ".join(normalized.get("hgvsg", []))),
                _component("spdi", spdi),
                _component("rsid", "; ".join(normalized.get("id", []))),
                _component("chrom-pos-ref-alt", f"{chrom}-{pos}-{ref}-{alt}" if pos else ""),
            ],
        )
        provenance = self.provenance_for(f"Observation/{resource.id}", request, payload)
        resource.extension.append({"url": "https://vwf-ethos.org/StructureDefinition/http-retries", "valueInteger": retries})
        return [resource, provenance], "success"


class OpenCravatAnnotator(HTTPMixin, BaseBiomedicalTool):
    name = "open_cravat"
    version = "run-opencravat-single-annotate"
    endpoint = "https://run.opencravat.org/submit/annotate"
    default_annotators = "gnomad,clinvar,revel,cadd,alphamissense,spliceai"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        chrom = request.parameters.get("chrom")
        pos = request.parameters.get("pos")
        ref = request.parameters.get("ref")
        alt = request.parameters.get("alt")
        if not all((chrom, pos, ref, alt)):
            raise ValueError("Normalized chrom/pos/ref/alt are required")
        params = {
            "chrom": str(chrom),
            "pos": str(pos),
            "ref_base": str(ref),
            "alt_base": str(alt),
            "annotators": str(request.parameters.get("annotators", self.default_annotators)),
        }
        status, payload, retries = self.request_json("GET", self.endpoint, params=params)
        if status != 200:
            raise RuntimeError(f"OpenCRAVAT returned HTTP {status}: {payload}")

        resources: list[FHIRResource] = []
        gnomad = payload.get("gnomad") or {}
        clinvar = payload.get("clinvar") or {}
        revel = payload.get("revel") or {}
        cadd = payload.get("cadd") or {}
        alphamissense = payload.get("alphamissense") or {}
        module_versions = payload.get("module_versions") or {}
        crx = payload.get("crx") or {}
        mappings = {
            "gnomad": ("gnomAD allele frequency", gnomad.get("af")),
            "clinvar": ("ClinVar clinical significance", clinvar.get("sig")),
            "revel": ("REVEL score", revel.get("score")),
            "cadd": ("CADD PHRED", cadd.get("phred")),
            "alphamissense": ("AlphaMissense pathogenicity", alphamissense.get("am_pathogenicity")),
            "spliceai": ("SpliceAI annotation", payload.get("spliceai")),
        }
        for key, (display, value) in mappings.items():
            if value is None:
                continue
            resource = observation(
                observation_id=f"ocravat-{key}-{request.patient_id}-{chrom}-{pos}-{ref}-{alt}",
                patient_id=request.patient_id,
                display=display,
                value=value,
                based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                components=[
                    _component("module_version", module_versions.get(key) or {}),
                    _component("gene", crx.get("hugo")),
                    _component("protein_change", crx.get("achange")),
                ],
            )
            resources.append(resource)
            resources.append(self.provenance_for(f"Observation/{resource.id}", request, payload))
        if not resources:
            return [operation_outcome("information", "no-annotation", "OpenCRAVAT returned no requested annotations")], "not_found"
        return resources, "success"


class GnomADProvider(HTTPMixin, BaseBiomedicalTool):
    name = "gnomad_graphql"
    version = "gnomad-r4-api"
    endpoint = "https://gnomad.broadinstitute.org/api"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        chrom = str(request.parameters.get("chrom", "")).removeprefix("chr")
        pos = request.parameters.get("pos")
        ref = request.parameters.get("ref")
        alt = request.parameters.get("alt")
        if not all((chrom, pos, ref, alt)):
            raise ValueError("Normalized chrom/pos/ref/alt are required")
        query = """
        query variant($variantId: String!, $dataset: DatasetId!) {
          variant(variantId: $variantId, dataset: $dataset) {
            variant_id
            chrom pos ref alt
            exome { ac an af }
            genome { ac an af }
          }
        }
        """
        variables = {"variantId": f"{chrom}-{pos}-{ref}-{alt}", "dataset": request.parameters.get("dataset", "gnomad_r4")}
        status, payload, retries = self.request_json(
            "POST", self.endpoint, json_payload={"query": query, "variables": variables}
        )
        if status != 200 or payload.get("errors"):
            raise RuntimeError(f"gnomAD returned HTTP {status}: {payload}")
        variant = payload.get("data", {}).get("variant")
        if not variant:
            return [operation_outcome("information", "not-found", "Variant absent from selected gnomAD dataset")], "not_found"
        resources: list[FHIRResource] = []
        for sequencing_type in ("exome", "genome"):
            values = variant.get(sequencing_type) or {}
            resource = observation(
                observation_id=f"gnomad-{sequencing_type}-{chrom}-{pos}-{ref}-{alt}",
                patient_id=request.patient_id,
                display=f"gnomAD {sequencing_type} allele frequency",
                value=values.get("af"),
                based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                components=[
                    _component("allele_count", values.get("ac")),
                    _component("allele_number", values.get("an")),
                    _component("dataset", variables["dataset"]),
                ],
            )
            resources.append(resource)
            resources.append(self.provenance_for(f"Observation/{resource.id}", request, payload))
        return resources, "success"


class ClinGenERepoProvider(HTTPMixin, BaseBiomedicalTool):
    name = "clingen_erepo"
    version = "evrepo-json-ld-rest"
    endpoint = "https://erepo.genome.network/evrepo/api/classifications"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        hgvs = request.parameters.get("hgvsg") or request.parameters.get("hgvs_c")
        if not hgvs:
            raise ValueError("hgvsg or hgvs_c is required")
        params = {"hgvs": str(hgvs), "matchMode": "exact", "matchLimit": 20}
        status, payload, retries = self.request_json("GET", self.endpoint, params=params)
        if status != 200:
            raise RuntimeError(f"ClinGen eRepo returned HTTP {status}: {payload}")
        interpretations = payload.get("variantInterpretations", [])
        if not interpretations:
            return [operation_outcome("information", "not-found", "No ClinGen interpretation matched the query")], "not_found"
        resources: list[FHIRResource] = []
        for index, item in enumerate(interpretations):
            guideline = (item.get("guidelines") or [{}])[0]
            outcome = guideline.get("outcome", {}).get("label", "not stated")
            resource = observation(
                observation_id=f"clingen-{index}-{sha256(str(item.get('@id', '')).encode()).hexdigest()[:20]}",
                patient_id=request.patient_id,
                display="ClinGen expert variant classification",
                value=outcome,
                based_on=[f"Observation/{request.variant_id}"] if request.variant_id else [],
                components=[
                    _component("caid", item.get("caid")),
                    _component("expert_panel", guideline.get("affiliation")),
                    _component("condition", item.get("condition", {}).get("label")),
                    _component("evidence_links", _stringify_evidence_links(item.get("evidenceLinks", []))),
                ],
            )
            resources.append(resource)
            resources.append(self.provenance_for(f"Observation/{resource.id}", request, payload))
        return resources, "success"


def _stringify_evidence_links(value: Any) -> str:
    if isinstance(value, list):
        return "; ".join(str(item) for item in value)
    if isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value or "")


class PubMedSearchProvider(HTTPMixin, BaseBiomedicalTool):
    name = "pubmed_eutils"
    version = "ncbi-eutils-v2"
    endpoint = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    api_key = os.getenv("NCBI_API_KEY")
    email = os.getenv("PUBMED_EMAIL", "research@example.org")

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        query = request.parameters.get("query") or self._default_query(request)
        retmax = int(request.parameters.get("retmax", 5))
        include_full_text = bool(request.parameters.get("include_full_text", False))
        articles, query_used = self._search_with_fallback(str(query), retmax)
        if not articles:
            return [operation_outcome("information", "not-found", "No PubMed articles matched the query")], "not_found"
        resources: list[FHIRResource] = []
        for article in articles:
            uid = article["pmid"]
            full_text = None
            pmc_id = None
            if include_full_text:
                pmc_id = self._elink_pmc(uid)
                if pmc_id:
                    full_text = self._efetch_pmc_full_text(pmc_id)
            resource = FHIRResource(
                resourceType="DocumentReference",
                id=f"pubmed-{uid}",
                status="current",
                code={"text": "PubMed case or mechanism report"},
                subject={"reference": f"Patient/{request.patient_id}"},
                description=article.get("title", ""),
                content=[{"attachment": {"url": f"https://pubmed.ncbi.nlm.nih.gov/{uid}/"}, "format": "text/html"}],
                extension=[
                    {"url": "https://vwf-ethos.org/StructureDefinition/pubmed-uid", "valueString": uid},
                    {"url": "https://vwf-ethos.org/StructureDefinition/pub-source", "valueString": article.get("journal", "")},
                    {"url": "https://vwf-ethos.org/StructureDefinition/pub-year", "valueString": article.get("year", "")},
                    {"url": "https://vwf-ethos.org/StructureDefinition/pub-abstract", "valueString": article.get("abstract", "")},
                    {"url": "https://vwf-ethos.org/StructureDefinition/pub-doi", "valueString": article.get("doi", "")},
                    {"url": "https://vwf-ethos.org/StructureDefinition/pub-full-text", "valueString": full_text or ""},
                    {"url": "https://vwf-ethos.org/StructureDefinition/pub-pmc-id", "valueString": pmc_id or ""},
                ],
            )
            resources.append(resource)
            resources.append(
                self.provenance_for(
                    f"DocumentReference/{resource.id}",
                    request,
                    {"query": query_used, "pmids": [article["pmid"] for article in articles], "articles": articles},
            )
        )
        return resources, "success"

    def _default_query(self, request: ToolRequest) -> str:
        protein = request.parameters.get("hgvs_p")
        cdna = request.parameters.get("hgvs_c")
        terms = [term for term in (protein, cdna) if term]
        return '("von Willebrand factor" OR VWF)' + (f" AND ({' OR '.join(terms)})" if terms else "")

    def _search_with_fallback(self, exact_query: str, retmax: int) -> tuple[list[dict[str, Any]], str]:
        queries = [
            exact_query,
            self._broad_query(exact_query),
        ]
        for query in dict.fromkeys(queries):
            ids = self._esearch(query, retmax)
            if not ids:
                continue
            articles = self._efetch(ids)
            if articles:
                return articles, query
        return [], exact_query

    def _esearch(self, query: str, retmax: int) -> list[str]:
        params = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": str(retmax),
            "sort": "relevance",
            "tool": "VWF_ETHOS_LangGraph",
            "email": self.email,
        }
        if self.api_key:
            params["api_key"] = self.api_key
        status, payload, _ = self.request_json("GET", f"{self.endpoint}/esearch.fcgi", params=params)
        if status != 200:
            raise RuntimeError(f"PubMed esearch returned HTTP {status}: {payload}")
        return payload.get("esearchresult", {}).get("idlist", [])

    def _efetch(self, ids: list[str]) -> list[dict[str, Any]]:
        params = {
            "db": "pubmed",
            "id": ",".join(ids),
            "retmode": "xml",
            "rettype": "abstract",
            "tool": "VWF_ETHOS_LangGraph",
            "email": self.email,
        }
        if self.api_key:
            params["api_key"] = self.api_key
        status, payload, _ = self.request_json("GET", f"{self.endpoint}/efetch.fcgi", params=params)
        if status != 200:
            raise RuntimeError(f"PubMed efetch returned HTTP {status}: {payload}")
        return self._parse_efetch_xml(str(payload.get("raw", "")))

    def _elink_pmc(self, pmid: str) -> str | None:
        params = {
            "dbfrom": "pubmed",
            "db": "pmc",
            "id": pmid,
            "retmode": "json",
            "tool": "VWF_ETHOS_LangGraph",
            "email": self.email,
        }
        if self.api_key:
            params["api_key"] = self.api_key
        status, payload, _ = self.request_json("GET", f"{self.endpoint}/elink.fcgi", params=params)
        if status != 200:
            return None
        for linkset in payload.get("linksets", []):
            for linksetdb in linkset.get("linksetdbs", []):
                if linksetdb.get("dbto") == "pmc" and linksetdb.get("linkname") == "pubmed_pmc":
                    links = linksetdb.get("links", [])
                    if links:
                        return str(links[0])
        return None

    def _efetch_pmc_full_text(self, pmc_id: str) -> str | None:
        params = {
            "db": "pmc",
            "id": pmc_id,
            "retmode": "xml",
            "tool": "VWF_ETHOS_LangGraph",
            "email": self.email,
        }
        if self.api_key:
            params["api_key"] = self.api_key
        status, payload, _ = self.request_json("GET", f"{self.endpoint}/efetch.fcgi", params=params)
        if status != 200:
            return None
        raw_xml = str(payload.get("raw", ""))
        try:
            import xml.etree.ElementTree as ET

            root = ET.fromstring(raw_xml)
            return "\n".join(part.strip() for part in root.itertext() if part.strip())
        except ET.ParseError:
            return raw_xml or None

    def _broad_query(self, exact_query: str) -> str:
        return (
            '("von Willebrand factor"[Title/Abstract] OR VWF[Title/Abstract]) '
            'AND ("von Willebrand disease"[MeSH Terms] OR "von Willebrand disease"[Title/Abstract]) '
            'AND (variant OR mutation OR missense) '
            'AND ("type 2"[Title/Abstract] OR 2A[Title/Abstract] OR 2B[Title/Abstract] '
            'OR 2M[Title/Abstract] OR 2N[Title/Abstract])'
        )

    def _parse_efetch_xml(self, xml_text: str) -> list[dict[str, Any]]:
        import xml.etree.ElementTree as ET

        articles: list[dict[str, Any]] = []
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return articles
        for article in root.findall(".//PubmedArticle"):
            pmid_element = article.find(".//PMID")
            title_element = article.find(".//ArticleTitle")
            journal_element = article.find(".//Journal/Title")
            year_element = article.find(".//PubDate/Year")
            doi_element = article.find(".//ArticleId[@IdType='doi']")
            abstract_parts: list[str] = []
            for part in article.findall(".//Abstract/AbstractText"):
                label = part.attrib.get("Label")
                text = "".join(part.itertext()).strip()
                if text:
                    abstract_parts.append(f"{label}: {text}" if label else text)
            articles.append(
                {
                    "pmid": pmid_element.text if pmid_element is not None else "",
                    "title": title_element.text if title_element is not None else "",
                    "journal": journal_element.text if journal_element is not None else "",
                    "year": year_element.text if year_element is not None else "",
                    "abstract": "\n".join(abstract_parts),
                    "doi": doi_element.text if doi_element is not None else "",
                }
            )
        return articles


class HGMDProvider(BaseBiomedicalTool):
    name = "hgmd"
    version = "license-gated"
    endpoint = "local://hgmd-license-gate"

    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        return [
            operation_outcome(
                "warning",
                "license-required",
                "HGMD has no public API and requires an institutional license; provider remains disabled rather than treating absence as benign.",
            )
        ], "disabled"

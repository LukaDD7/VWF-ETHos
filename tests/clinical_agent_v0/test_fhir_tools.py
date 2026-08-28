from __future__ import annotations

from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.tools.fhir import FHIRBundle
from src.vwd_clinical_agent.tools.base import ToolRequest
from src.vwd_clinical_agent.tools.matrix import EvidenceToolMatrix
from src.vwd_clinical_agent.tools.second_level import SecondLevelFHIRStore, SecondLevelRecord
from src.vwd_clinical_agent.tools.online import PubMedSearchProvider
from src.vwd_clinical_agent.tools.fhir import FHIRBundle, observation
from src.vwd_clinical_agent.evidence_analysis import analyze_evidence_conflicts
from src.vwd_clinical_agent.run_archive import RunArchive
from src.vwd_clinical_agent.tools.full_text import PubMedFullTextSearchTool
from src.vwd_clinical_agent.azure import DeterministicLLMProvider
from src.vwd_clinical_agent.graph import build_workflow
from src.vwd_clinical_agent.workbook import LocalWorkbookProvider


def test_second_level_rejects_action_outside_fhir_space() -> None:
    store = SecondLevelFHIRStore.from_actions("CASE_TEST", ["RIPA"])
    with pytest.raises(ValueError):
        store.record_observation(SecondLevelRecord(action="NEW_UNAPPROVED_TEST", value=1.0))


def test_observation_requires_service_request() -> None:
    store = SecondLevelFHIRStore("CASE_TEST")
    with pytest.raises(ValueError, match="No ServiceRequest"):
        store.record_observation(SecondLevelRecord(action="RIPA", value=1.0))


def test_real_observation_enables_final_report() -> None:
    store = SecondLevelFHIRStore.from_actions("CASE_TEST", ["RIPA", "VWF_MULTIMER"])
    store.record_observation(
        SecondLevelRecord(action="RIPA", value=0.8, unit="ratio", operator="lab", method="platelet-rich plasma")
    )
    store.mark_unavailable("VWF_MULTIMER", "not available at institution")
    report = store.finalize("RIPA observed; multimer unavailable and not imputed.")
    assert report.status == "final"
    assert store.is_final()
    assert len(store.bundle.resources("Observation")) == 1
    assert len(store.bundle.resources("Task")) == 1


def test_evidence_matrix_blocks_without_second_level_report() -> None:
    matrix = EvidenceToolMatrix()
    result = matrix.run(
        patient_id="CASE_TEST",
        variant_id="CASE_TEST",
        hgvs_c="NM_000552.5:c.4499C>T",
        hgvs_p="p.Ala1500Val",
        second_level_bundle=FHIRBundle(),
    )
    assert result.authorized is False
    assert result.bundle.resources("OperationOutcome")[0]["issue"][0]["code"] == "business-rule"
    assert result.calls == []


def test_graph_auto_environment_marks_tests_unavailable_and_continues() -> None:
    cases, _ = LocalWorkbookProvider(ROOT / "data/clinical_agent_pilot/vwd_agentic_workflow_deidentified_v3.xlsx").load_cases()
    case = next(item for item in cases if item.patient_id == "CASE_001")
    graph = build_workflow(DeterministicLLMProvider())
    result = graph.invoke(
        {
            "run_id": "auto-env-test",
            "case": case,
            "mode": "retrospective",
            "trace": [],
            "provider_calls": [],
            "second_level_environment": "auto_unavailable",
        }
    )
    assert result["second_level_status"] == "not_available"
    bundle = result["second_level_bundle"]
    requests = bundle.resources("ServiceRequest")
    tasks = bundle.resources("Task")
    assert requests and all(request["status"] == "revoked" for request in requests)
    assert tasks and all(task["status"] == "rejected" for task in tasks)
    assert bundle.resources("Observation") == []
    assert any(report["status"] == "final" for report in bundle.resources("DiagnosticReport"))
    assert result["status"] == "completed"


def test_gene_first_flow_limits_second_level_actions() -> None:
    cases, _ = LocalWorkbookProvider(ROOT / "data/clinical_agent_pilot/vwd_agentic_workflow_deidentified_v3.xlsx").load_cases()
    case = next(item for item in cases if item.patient_id == "CASE_001")
    graph = build_workflow(DeterministicLLMProvider())
    result = graph.invoke(
        {
            "run_id": "gene-first-test",
            "case": case,
            "mode": "retrospective",
            "trace": [],
            "provider_calls": [],
            "second_level_environment": "auto_unavailable",
        }
    )
    actions = [action.action_code for action in result["recommended_actions"]]
    assert len(actions) <= 3
    assert "VWF_MULTIMER" in actions
    assert "VWF_CB" in actions
    assert "DDAVP_0_1_4H" not in actions


def test_code_as_search_policy_rejects_unapproved_operation() -> None:
    matrix = EvidenceToolMatrix()
    response = matrix.registry.invoke(
        "local_clingen_snapshot",
        ToolRequest(
            operation="arbitrary_prompt_tool_call",
            patient_id="CASE_TEST",
            variant_id="CASE_TEST",
            parameters={"path": "unused", "hgvs_c": "c.1A>T"},
        ),
    )
    assert response.status == "error"
    assert "not allowed" in response.diagnostics[0]


def test_pubmed_efetch_xml_parser_extracts_abstract_and_doi() -> None:
    xml = """
    <PubmedArticleSet>
      <PubmedArticle>
        <MedlineCitation>
          <PMID>12345678</PMID>
          <Article>
            <Journal><Title>Journal of VWF Research</Title></Journal>
            <ArticleTitle>Test VWF article</ArticleTitle>
            <Abstract>
              <AbstractText Label="BACKGROUND">VWF background.</AbstractText>
              <AbstractText Label="RESULTS">VWF result.</AbstractText>
            </Abstract>
          </Article>
        </MedlineCitation>
        <PubmedData>
          <ArticleIdList><ArticleId IdType="doi">10.1000/test</ArticleId></ArticleIdList>
        </PubmedData>
      </PubmedArticle>
    </PubmedArticleSet>
    """
    provider = PubMedSearchProvider()
    articles = provider._parse_efetch_xml(xml)
    assert len(articles) == 1
    assert articles[0]["pmid"] == "12345678"
    assert articles[0]["title"] == "Test VWF article"
    assert "BACKGROUND: VWF background." in articles[0]["abstract"]
    assert articles[0]["doi"] == "10.1000/test"


def test_conflict_analyzer_detects_classification_frequency_conflict() -> None:
    bundle = FHIRBundle.of(
        [
            observation(observation_id="clinvar", patient_id="CASE_TEST", display="ClinVar clinical significance", value="Pathogenic"),
            observation(observation_id="gnomad", patient_id="CASE_TEST", display="gnomAD exome allele frequency", value=0.02),
        ]
    )
    result = analyze_evidence_conflicts(fhir_bundle=bundle, second_level_status="not_available", candidate_subtypes=["type_2_candidate"])
    conflicts = {conflict.conflict_id for conflict in result["evidence_conflicts"]}
    assert "classification_vs_population_frequency" in conflicts
    assert "type2_discriminator_unavailable" in conflicts


def test_conflict_analyzer_reports_missing_evidence_without_hard_conflict() -> None:
    bundle = FHIRBundle.of(
        [
            observation(observation_id="clinvar", patient_id="CASE_TEST", display="ClinVar clinical significance", value="Uncertain significance"),
            observation(observation_id="gnomad", patient_id="CASE_TEST", display="gnomAD exome allele frequency", value=0.00001),
        ]
    )
    result = analyze_evidence_conflicts(fhir_bundle=bundle, second_level_status="observed", candidate_subtypes=["type_1_candidate_provisional"])
    assert result["evidence_conflicts"]
    assert result["evidence_conflicts"][0].conflict_type == "insufficient_evidence"


def test_run_archive_writes_complete_case_artifacts() -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as temp_dir:
        _run_archive_assertions(temp_dir)


def _run_archive_assertions(tmp_path: str) -> None:
    from pathlib import Path

    tmp_path = Path(tmp_path)
    archive = RunArchive("run-test", tmp_path)
    archive.write_debug_event("CASE_TEST", {"type": "task", "payload": {"name": "load_case"}})
    archive.write_state_history("CASE_TEST", [{"checkpoint_id": "1", "values": {"status": "running"}}])
    archive.write_final_state("CASE_TEST", {"status": "completed"})
    archive.write_report(
        "CASE_TEST",
        {
            "confidence": "low",
            "abstention": True,
            "expert_review_required": True,
            "opinion": "test",
            "candidate_subtypes": ["type_2_candidate"],
            "supporting_evidence": ["Observation/test"],
            "missing_information": ["RIPA"],
            "limitations": ["Second level unavailable"],
        },
    )
    case_dir = tmp_path / "run-test" / "cases" / "CASE_TEST"
    assert (case_dir / "debug_trace.jsonl").exists()
    assert (case_dir / "state_history.jsonl").exists()
    assert (case_dir / "final_state.json").exists()
    assert (case_dir / "final_report.json").exists()
    assert (case_dir / "final_report.md").exists()


def test_full_text_search_extracts_diverse_clinical_excerpts() -> None:
    tool = PubMedFullTextSearchTool()
    text = "Ristocetin induced platelet aggregation was abnormal. Multimer analysis showed loss of high molecular weight bands. DDAVP response was partial."
    excerpts = tool._extract_excerpts(text, ["ristocetin", "RIPA", "multimer", "DDAVP"], 4)
    terms = {excerpt["term"].casefold() for excerpt in excerpts}
    assert {"ristocetin", "multimer", "ddavp"}.issubset(terms)
    assert all(len(excerpt["text"]) > 20 for excerpt in excerpts)

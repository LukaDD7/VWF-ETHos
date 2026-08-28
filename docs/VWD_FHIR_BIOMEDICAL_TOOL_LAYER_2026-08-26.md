# VWD FHIR Biomedical Tool Layer

## Status

This layer upgrades the V0 workflow from prompt/rule-only extraction to a FHIR-native, code-as-search tool interface. It is research-only and does not autonomously diagnose or treat patients.

For a concise AI-engineer explanation of each biomedical tool, see [VWF_BIOMEDICAL_TOOLS_FOR_AI_ENGINEERS.md](./VWF_BIOMEDICAL_TOOLS_FOR_AI_ENGINEERS.md).

## Core guarantee

Second-level actions are represented as FHIR `ServiceRequest` resources. A result can only enter the workflow as a real `Observation` referencing its request. If a test is unavailable, a rejected `Task` and reason are recorded. A final `DiagnosticReport` is emitted only after requests are completed or explicitly revoked; the system never imputes an unobserved result.

For the current retrospective pilot, the workflow is **gene-first**:

```text
First Level
→ genetic variant
→ biomedical evidence tools / literature
→ structured conflict analysis
→ select a small set of second-level tests
→ automated unavailable feedback
→ final synthesis
```

The automated environment returns every requested second-level test as unavailable. The workflow records `Task(status=rejected)` for each request, emits a final `DiagnosticReport`, and continues to final synthesis. It does not wait for a human and does not impute values.

Online variant tools remain gated by a final second-level `DiagnosticReport` unless a research bypass is explicitly supplied for interface testing. In automated retrospective mode that report records unavailability rather than observations.

## Unified interface

Every tool implements the same typed contract:

```text
ToolRequest
- operation
- patient_id
- variant_id
- parameters
- input FHIR Bundle

ToolResponse
- status: success | not_found | error | disabled
- output FHIR Bundle
- diagnostics
- retry_count
- latency_ms
- cache_hit
```

Each response contains:

- FHIR `Observation` or `DocumentReference`;
- FHIR `Provenance` with request/response hashes;
- `OperationOutcome` for failure, absence, or disabled states.

A not-found or disabled result is not benign evidence.

## Online tools

| Need | Tool / API | Output |
|---|---|---|
| HGVS normalization | Ensembl Variant Recoder | GRCh38 SPDI, HGVSg, rsID, chrom-pos-ref-alt |
| Population frequency | gnomAD GraphQL v4 | exome/genome AC, AN, AF |
| Aggregated annotation | OpenCRAVAT single-variant API | ClinVar, REVEL, CADD, AlphaMissense, SpliceAI, gnomAD summary |
| Expert curation | ClinGen Evidence Repository API | VCEP classification, condition, evidence links |
| Literature | NCBI PubMed E-utilities | PMID, title, source, URL |
| Licensed database | HGMD | Disabled interface only; no public API/license assumption |

`NCBI_API_KEY` is loaded from the ignored local `.env` file and appended to PubMed E-utility requests.

If the official gnomAD GraphQL endpoint returns a transient `429`, OpenCRAVAT’s aggregate gnomAD annotation remains available as a degraded population-frequency fallback. The official API call is still recorded as an error; the fallback is not silently promoted to official exome/genome resolution.

PubMed search now uses:

1. `esearch.fcgi` to retrieve PMIDs;
2. `efetch.fcgi` to retrieve full abstract XML;
3. structured parsing of title, journal, year, DOI, and abstract;
4. a FHIR `DocumentReference` per article.

If an exact variant query has no match, the provider falls back to a broader but still VWF/VWD-specific Boolean query. Absence of an exact-match article is reported as `not_found`, not as benign evidence.

The provider also supports PubMed Central full-text retrieval:

1. `elink.fcgi` maps PubMed PMID to PMC ID when available;
2. `efetch.fcgi` retrieves PMC XML;
3. article text is extracted into `pub-full-text`;
4. PMC ID is stored in `pub-pmc-id`.

Enable with:

```bash
--pubmed-full-text
```

After retrieval, `pubmed_full_text_search` performs code-as-search over the stored PMC text. It returns bounded excerpts around clinically relevant terms such as:

- ristocetin / RIPA
- multimer
- collagen binding
- factor VIII binding
- DDAVP / desmopressin
- clearance
- type 2A / 2B / 2M / 2N

The LLM does not receive the entire full text. It receives typed excerpt Observations with provenance, preventing context overflow and reducing unsupported paraphrasing.

### Context-window and ambiguity controls

Each full-text excerpt now contains:

- a bounded context window;
- the matched keyword;
- the nearest explicit variant term;
- the distance to that variant term;
- a `variant_linked` flag;
- the document-level `variant_specific` flag.

Default context sizes are:

```text
full-text keyword search: 600 characters before + 900 characters after
phenotype extraction:     500 characters before + 700 characters after
variant-link radius:      1200–1500 characters
```

The subtype-tendency node only uses literature observations when both:

```text
variant_specific = true
variant_linked   = true
```

This prevents a generic review-article sentence such as “multimers are often reduced in type 2A” from being misattributed to the patient’s specific variant when the nearest explicit variant mention is thousands of characters away.

The final LLM synthesis node now receives the structured FHIR evidence bundle and is required to cite FHIR resource IDs for factual statements. Uncited or unsupported statements are not allowed by the synthesis contract.

## Code-as-search enforcement

Tool use is not prompt-driven. The LLM has no direct tool-calling interface. LangGraph nodes call the Python `ToolRegistry`, which validates both the tool name and operation against `src/vwd_clinical_agent/tools/tool_policy.json`. Unknown tools and unapproved operations return an `OperationOutcome` with status `error`. Every successful response is a typed FHIR Bundle with provenance; every failure, absence, or disabled license is represented explicitly and is never treated as benign evidence.

## Structured conflict interpretation

The graph now includes an `analyze_evidence` node after patient/variant integration. It converts evidence disagreements into typed `EvidenceConflict` records rather than allowing the LLM to improvise reconciliation.

Current conflict classes:

- pathogenic classification vs high population allele frequency;
- disagreement among REVEL, CADD, and AlphaMissense;
- ClinVar vs ClinGen expert classification;
- type 2 candidate with all second-level discriminators unavailable;
- insufficient evidence without a hard contradiction.

The final report includes:

- conflict descriptions;
- evidence resource references;
- missing evidence items;
- conflict explanations;
- safety limitations.

## Complete run archive

Every workflow invocation now writes a complete run archive:

```text
run_archive/<run_id>/
  checkpoints.sqlite
  run_manifest.json
  cases/<case_id>/
    debug_trace.jsonl
    state_history.jsonl
    final_state.json
    final_report.json
    final_report.md
```

- `checkpoints.sqlite` is a LangGraph SQLite checkpointer and stores the full replayable state.
- `debug_trace.jsonl` records each node’s input keys, output delta, errors, and task metadata.
- `state_history.jsonl` exports checkpoint-by-checkpoint state snapshots.
- `final_state.json` stores the complete terminal state.
- `final_report.json` and `final_report.md` store the structured and human-readable final report.

## Local tools and snapshots

`scripts/download_vwf_biomedical_snapshots.py` downloads compact VWF-specific snapshots:

```bash
python scripts/download_vwf_biomedical_snapshots.py \
  --output-dir data/external/vwf_biomedical_snapshots
```

This currently records:

- ClinGen eRepo VWF interpretations;
- ClinVar VWF ESummary records.

Both snapshots are queryable through `local_clingen_snapshot` and `local_clinvar_snapshot`. They are offline fallbacks only: the online ClinGen API and OpenCRAVAT/ClinVar annotation remain primary. The local snapshots are not used merely because an online query returns `not_found`; that would incorrectly treat a smaller local snapshot as authoritative.

Large resources are not blindly downloaded. They are represented by the local score store and require a versioned file plus SHA-256:

- dbNSFP/REVEL: large genome-wide score table; use a VWF subset or licensed institutional mirror.
- CADD: versioned score data has redistribution/licensing constraints; online OpenCRAVAT is used when available.
- AlphaMissense: genome-wide predictions are very large; a VWF/transcript subset can be loaded by `LocalScoreStore`.
- SpliceAI: model/database access varies; online annotation is preferred when available.
- gnomAD: official online API is used; local whole-genome VCFs are unnecessary for one-variant analysis.
- HGMD: institutional license required; adapter remains disabled.

## FHIR second-level interaction

Create orders:

```bash
python scripts/manage_vwd_second_level_fhir.py create \
  --patient-id CASE_001 \
  --actions RIPA VWF_MULTIMER VWF_CB \
  --rationale "Type 2 cannot be excluded" \
  --output output/vwd_langgraph_v0/case001/second_level.json
```

Record a real result:

```bash
python scripts/manage_vwd_second_level_fhir.py record \
  --patient-id CASE_001 \
  --bundle output/vwd_langgraph_v0/case001/second_level.json \
  --action RIPA --value 0.8 --unit ratio \
  --operator "lab_name" --method "platelet-rich plasma"
```

Mark unavailable:

```bash
python scripts/manage_vwd_second_level_fhir.py unavailable \
  --patient-id CASE_001 \
  --bundle output/vwd_langgraph_v0/case001/second_level.json \
  --action VWF_MULTIMER --reason "not available at institution"
```

Finalize:

```bash
python scripts/manage_vwd_second_level_fhir.py finalize \
  --patient-id CASE_001 \
  --bundle output/vwd_langgraph_v0/case001/second_level.json \
  --conclusion "RIPA observed; multimer unavailable and not imputed."
```

Run gated evidence tools:

```bash
python scripts/run_vwd_biomedical_tools.py \
  --case-id CASE_001 --variant-index 1 \
  --second-level-bundle output/vwd_langgraph_v0/case001/second_level.json \
  --output-dir output/vwd_langgraph_v0/case001/tools
```

Research-only interface test:

```bash
python scripts/run_vwd_biomedical_tools.py \
  --case-id CASE_001 --variant-index 1 --research-bypass \
  --output-dir output/vwd_langgraph_v0/tool_tests/case001
```

## LangGraph integration

Passing `--biomedical-tools --second-level-bundle ...` inserts `run_fhir_evidence_tools` between local workbook evidence and patient/variant integration. Without a final second-level report, the node emits `fhir_evidence_tools_not_authorized` and makes zero online requests.

## Initial integrated test

For `CASE_001` / `NM_000552.5:c.4499C>T` / `p.Ala1500Val`:

- Ensembl normalization: success, GRCh38 `chr12:6018919 G>A`
- gnomAD: exome AF `4.310138403333174e-05`, genome AF `9.852475598702101e-05`
- OpenCRAVAT:
  - ClinVar: Uncertain significance
  - REVEL: 0.201
  - CADD: 9.552
  - AlphaMissense: 0.106
  - SpliceAI: no score returned
- ClinGen eRepo: not found for this exact variant
- PubMed: exact-variant query is enforced; absence is reported as not-found rather than broad fallback
- HGMD: disabled by license gate

All generated Observations and DocumentReferences had matching Provenance resources.

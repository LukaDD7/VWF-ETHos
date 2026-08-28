# VWD LangGraph V0 Runbook

## Scope

This runbook describes the first runnable, research-only LangGraph workflow for the de-identified VWD clinical-agent pilot workbook. It is not a diagnostic system and does not provide treatment advice.

## Inputs

- Workbook: `data/clinical_agent_pilot/vwd_agentic_workflow_deidentified_v3.xlsx`
- First-level sheet: `1.基因前`
- Genetic sheet: `2.基因后`
- Current audit result: 47 unique patients, 59 patient/variant rows, and 10 multi-variant patients.

The identifier `CASE_021_VARIANT_2_BENIGN` is normalized to patient `CASE_021`, variant 2, with `benign_reported=true`; it is not treated as a separate patient.

## Offline smoke run

```bash
python scripts/run_vwd_langgraph_v0.py \
  --workbook data/clinical_agent_pilot/vwd_agentic_workflow_deidentified_v3.xlsx \
  --mode retrospective \
  --provider-profile offline \
  --llm-provider deterministic \
  --output-dir output/vwd_langgraph_v0/smoke
```

Expected result:

- 47/47 cases complete.
- 59/59 variant rows remain linked to their patient.
- All recommended actions come from the fixed action enum.
- No second-level result is fabricated.
- Every case abstains and is routed to expert review because units, reference ranges, second-level observations, external evidence, and expert gold labels are unavailable.

## Azure OpenAI mode

Set these environment variables before running:

```bash
export AZURE_OPENAI_ENDPOINT="https://YOUR_RESOURCE.openai.azure.com"
export AZURE_OPENAI_API_KEY="YOUR_KEY"
export AZURE_OPENAI_DEPLOYMENT="YOUR_DEPLOYMENT"
export AZURE_OPENAI_API_VERSION="2024-10-21"
export SUBAGENT_MODEL_NAME="azure/gpt-5.4-nano"
```

Optional:

```bash
export AZURE_OPENAI_TIMEOUT="30"
export AZURE_OPENAI_JSON_MODE="1"
```

For local-only secrets, use a `.env` file. The repository ignores `.env` and `.env.*`; do not commit them.

Run:

```bash
python scripts/run_vwd_langgraph_v0.py \
  --mode retrospective \
  --llm-provider azure \
  --output-dir output/vwd_langgraph_v0/azure
```

The Azure adapter calls the Chat Completions REST endpoint directly and requests JSON mode. `SUBAGENT_MODEL_NAME=azure/gpt-5.4-nano` is normalized to the Azure deployment `gpt-5.4-nano`, and the request uses `max_completion_tokens` for newer GPT models. The model is only used for structured retrospective summaries; clinical routing and action selection remain policy-controlled.

## Evaluation

For the FHIR tool layer, see [VWD_FHIR_BIOMEDICAL_TOOL_LAYER_2026-08-26.md](./VWD_FHIR_BIOMEDICAL_TOOL_LAYER_2026-08-26.md).

```bash
python scripts/evaluate_vwd_langgraph_v0.py \
  --predictions output/vwd_langgraph_v0/smoke/cases.jsonl \
  --output output/vwd_langgraph_v0/smoke/eval_metrics.json
```

With an expert gold file:

```bash
python scripts/evaluate_vwd_langgraph_v0.py \
  --predictions output/vwd_langgraph_v0/smoke/cases.jsonl \
  --gold annotations/vwd_action_gold_v0.csv \
  --output output/vwd_langgraph_v0/smoke/gold_metrics.json
```

The gold CSV must contain:

```text
patient_id,required_missing_information,preferred_actions,acceptable_actions,
inappropriate_actions,must_not_miss_actions,pre_genetic_subtype_set,
final_subtype_set_if_available,abstention_expected,rationale
```

## Output contract

Each run writes:

- `run_manifest.json`
- `cases.jsonl`
- `summary.csv`
- `trace.jsonl`
- `provider_calls.jsonl`
- `annotation_template.csv`
- `metrics.json`

## Current limitations

- The workbook has no units, reference ranges, assay methods, or collection times.
- No real second-level results are present.
- No ClinVar, gnomAD, HGMD, ClinGen ERepo, predictor, or guideline provider is enabled in this offline V0.
- Multi-variant phase is unknown and is never inferred.
- No expert action/subtype gold set exists, so clinical accuracy is intentionally not reported.

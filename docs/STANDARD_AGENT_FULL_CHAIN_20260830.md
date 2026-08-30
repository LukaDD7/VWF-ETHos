# Standard Agent full-chain run

Date: 2026-08-30

## Default profile

The standard runner profile is now the default:

```bash
python scripts/run_vwd_langgraph_v0.py \
  --workbook outputs/type2b_panel_agent_20260829/type2b_panel_agent_package.xlsx \
  --output-dir outputs/computational_panel_20260829/agent_run_type2b_standard_full_chain
```

No extra component flags are required. The standard profile enables:

- Patient/context FHIR conversion
- Local computational panels
- Biomedical evidence tools
- PubMed full text when available
- Deterministic mechanism synthesis

Optional reduced profiles:

```text
--profile offline   # local computational panels, no online biomedical tools
--profile minimal   # local workbook/context only
```

## Information flow

```text
Workbook
  ↓
Patient/context FHIR Observations
  ↓
Local computational panels
  ├─ AlphaGenome
  ├─ Boltz
  └─ MD
  ↓
FHIR ServiceRequest + Observation + Provenance
  ↓
Biomedical evidence tools
  ├─ Ensembl
  ├─ gnomAD
  ├─ OpenCRAVAT / ClinVar
  ├─ ClinGen
  ├─ PubMed
  └─ HGMD
  ↓
Evidence integration
  ↓
Mechanism analysis
  ↓
Subtype tendency
  ↓
Final opinion
```

## FHIR request/response contract

Computational tools are represented as:

```text
ServiceRequest
  ↓
Observation
  ↓
Provenance
```

If a result is unavailable, the corresponding response is:

```text
OperationOutcome(status=not-available)
```

Patient context and first-level labs are represented as FHIR Observations.

## Current full-chain outputs

```text
outputs/computational_panel_20260829/agent_run_type1_standard_full_chain/
outputs/computational_panel_20260829/agent_run_type2b_standard_full_chain/
```

Each run contains:

```text
cases.jsonl
trace.jsonl
provider_calls.jsonl
metrics.json
run_manifest.json
```

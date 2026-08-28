# VWD Biomedical Tool Skill

## Contract

All tools are code adapters, not prompts. A tool call must provide:

1. a typed `ToolRequest`;
2. a FHIR Bundle response;
3. success/not-found/error/disabled status;
4. retry and cache state;
5. FHIR Provenance with request and response hashes;
6. an `OperationOutcome` for every failure or disabled state.

## Clinical gating

The evidence matrix refuses online variant evidence unless the input FHIR Bundle contains a final second-level `DiagnosticReport`, or the operator explicitly sets a research bypass flag. A missing or failed tool is never interpreted as benign evidence.

## Tool plan

1. Ensembl Variant Recoder: normalize VWF HGVS to GRCh38 SPDI and chrom-pos-ref-alt.
2. gnomAD GraphQL: exome/genome AC, AN, AF.
3. OpenCRAVAT: ClinVar, gnomAD summary, REVEL, CADD, AlphaMissense, SpliceAI.
4. ClinGen Evidence Repository: VCEP interpretations and evidence links.
5. PubMed E-utilities: case reports and mechanism studies.
6. HGMD: license-gated; disabled by default.
7. Local score snapshots: versioned TSV/CSV lookup for offline execution.
8. Local guideline snapshot: expert-approved ClinGen/ACMG JSON rules.
9. Repo mechanism artifacts: read-only AlphaFold/Boltz/MD/FoldX adapter; no clinical override.

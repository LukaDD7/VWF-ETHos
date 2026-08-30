# VWD Research Report — CASE_T1_006

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

### 1. 临床与实验室表型
VWF:Ag=47.6, VWF:Act=34.1, FVIII:C=78.7, 血小板计数=未提供。
VWF:Act/VWF:Ag 比值未显示明显不成比例的功能下降。

### 2. 变异与功能域
- c.7219dup / p.Leu2407ProfsTer11（未能映射到已知功能域）。

### 3. AI 机制证据
- alphagenome_full_profile: AlphaGenome complete-profile selected-track scores; top absolute views: ag_chip_histone_active_abs_max=8438, ag_chip_tf_active_abs_max=1255, ag_cage_active_abs_max=440, ag_dnase_active_abs_max=28.18, ag_rna_seq_active_abs_max=16.11.

### 4. 机制解释

### 5. 分型鉴别与不确定性
HGMD/ClinVar 的致病性评级本身不能直接给出 VWD 亚型方向；即使标注为致病或 uncertain，也必须结合实验室检查、出血表型和机制证据综合判断。
当前倾向：type_1_or_low_vwf（low）。
建议优先补充：VWF_PP, DDAVP_0_1_4H。

## Candidate subtypes

type_1_candidate_provisional

## Supporting evidence

- CASE_T1_006
- CASE_T1_006:alphagenome-full-profile
- Observation/context-vwf-antigen
- Observation/context-vwf-activity
- Observation/context-factor-viii-activity
- Observation/context-platelet-count
- Observation/context-isth-bleeding-assessment-tool-score
- Observation/context-high-dose-ristocetin-platelet-aggregation
- Observation/context-reported-symptoms
- Observation/context-family-history
- Observation/context-prior-treatment
- Observation/context-comorbidity
- Observation/context-interpretation-constraints
- Observation/context-ddavp-reported-response
- Observation/computational-result-456c1334aa87b1bb6ac8
- Observation/vwf-domain-2407-Leu-ProfsTer

## Missing information

- expert_variant_classification
- population_frequency
- pathogenicity_prediction
- variant_specific_literature
- ABO_blood_group_genotype
- VWF_clearance_context

## Limitations

- ensembl_variant_recoder returned error; this is not benign evidence.
- coordinate_dependent_tools returned error; this is not benign evidence.
- clingen_erepo returned error; this is not benign evidence.
- pubmed_eutils returned error; this is not benign evidence.
- local_clingen_snapshot returned not_found; this is not benign evidence.
- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- The remaining uncertainty is due to missing evidence rather than contradictory evidence.

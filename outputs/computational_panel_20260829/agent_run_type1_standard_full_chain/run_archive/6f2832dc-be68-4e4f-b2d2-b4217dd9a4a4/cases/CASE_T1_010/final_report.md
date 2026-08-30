# VWD Research Report — CASE_T1_010

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

### 1. 临床与实验室表型
VWF:Ag=42.7, VWF:Act=36.1, FVIII:C=49.1, 血小板计数=未提供。
VWF:Act/VWF:Ag 比值未显示明显不成比例的功能下降。

### 2. 变异与功能域
- c.5827C>T / p.Arg1943Cys（near VWFD4）。

### 3. AI 机制证据
- alphagenome_full_profile: AlphaGenome complete-profile selected-track scores; top absolute views: ag_chip_histone_active_abs_max=1.008e+04, ag_chip_tf_active_abs_max=1472, ag_cage_active_abs_max=296.1, ag_dnase_active_abs_max=32.24, ag_procap_active_abs_max=17.04.
- boltz2_functional_panel: Assay-matched Boltz-2 structural deltas: d4_assembly_context delta_vs_WT=0.02988.

### 4. 机制解释

### 5. 分型鉴别与不确定性
HGMD/ClinVar 的致病性评级本身不能直接给出 VWD 亚型方向；即使标注为致病或 uncertain，也必须结合实验室检查、出血表型和机制证据综合判断。
当前倾向：type_1_or_low_vwf（low）。
建议优先补充：VWF_PP, DDAVP_0_1_4H。

## Candidate subtypes

type_1_candidate_provisional

## Supporting evidence

- CASE_T1_010
- CASE_T1_010:alphagenome-full-profile
- CASE_T1_010:VWF_R1943C:boltz2-functional-panel
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
- Observation/computational-result-d40c15402187fb3d166c
- Observation/computational-result-36dcf5fbf06e0a02343b
- Observation/vwf-domain-1943-Arg-Cys

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

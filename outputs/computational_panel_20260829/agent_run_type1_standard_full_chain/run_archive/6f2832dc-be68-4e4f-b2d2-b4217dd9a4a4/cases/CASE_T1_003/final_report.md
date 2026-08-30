# VWD Research Report — CASE_T1_003

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

### 1. 临床与实验室表型
VWF:Ag=36.4, VWF:Act=28.5, FVIII:C=68.6, 血小板计数=未提供。
VWF:Act/VWF:Ag 比值未显示明显不成比例的功能下降。

### 2. 变异与功能域
- c.7464C>T / p.Gly2488=（VWFC2）。

### 3. AI 机制证据
- alphagenome_full_profile: AlphaGenome complete-profile selected-track scores; top absolute views: ag_chip_histone_active_abs_max=8920, ag_chip_tf_active_abs_max=1248, ag_cage_active_abs_max=98.19, ag_dnase_active_abs_max=31.19, ag_rna_seq_active_abs_max=16.34. Splice-axis signal: ag_splice_junctions_abs_max=6.375, ag_splice_sites_abs_max=0.9793, ag_splice_site_usage_abs_max=0.7578.

### 4. 机制解释
AlphaGenome 提示剪接相关信号，需考虑变异通过转录/剪接异常导致 VWF 表达或分泌下降，而不一定直接改变蛋白结构。

### 5. 分型鉴别与不确定性
HGMD/ClinVar 的致病性评级本身不能直接给出 VWD 亚型方向；即使标注为致病或 uncertain，也必须结合实验室检查、出血表型和机制证据综合判断。
当前倾向：type_1_or_low_vwf（low）。
建议优先补充：VWF_PP, DDAVP_0_1_4H。

## Candidate subtypes

type_1_candidate_provisional

## Supporting evidence

- CASE_T1_003
- CASE_T1_003:alphagenome-full-profile
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
- Observation/computational-result-85b6bbaca9593735899f
- Observation/vwf-domain-2488-Gly-=

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

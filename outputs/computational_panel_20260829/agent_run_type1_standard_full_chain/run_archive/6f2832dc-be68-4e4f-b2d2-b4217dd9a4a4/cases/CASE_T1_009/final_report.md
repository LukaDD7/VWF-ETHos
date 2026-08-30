# VWD Research Report — CASE_T1_009

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

### 1. 临床与实验室表型
VWF:Ag=14.0, VWF:Act=7.5, FVIII:C=4.9, 血小板计数=未提供。
VWF:Act/VWF:Ag 比值为 0.536，低于 0.70，提示存在不成比例的功能性 VWF 缺陷。

### 2. 变异与功能域
- c.3614G>A / p.Arg1205His（near TIL4）。

### 3. AI 机制证据
- alphagenome_full_profile: AlphaGenome complete-profile selected-track scores; top absolute views: ag_chip_histone_active_abs_max=8845, ag_chip_tf_active_abs_max=1166, ag_cage_active_abs_max=135.8, ag_dnase_active_abs_max=21.86, ag_rna_seq_active_abs_max=16.81.
- boltz2_functional_panel: Assay-matched Boltz-2 structural deltas: dprime_d3_fviii_binding delta_vs_WT=-0.0147.

### 4. 机制解释
该变异位于 D'/D3 相关区域，该区域参与 FVIII 结合并影响 2N 轴；若 FVIII:C 与 VWF:Ag 不成比例下降，应进一步评估 FVIII 结合功能。
实验室表型与功能性 VWF 缺陷方向一致；若同时存在高危分子量多聚体缺失，更支持 2A/2B 样机制。

### 5. 分型鉴别与不确定性
HGMD/ClinVar 的致病性评级本身不能直接给出 VWD 亚型方向；即使标注为致病或 uncertain，也必须结合实验室检查、出血表型和机制证据综合判断。
当前倾向：type_2_vwd（moderate）。
VWF 多聚体分析和 VWF:CB/Ag 比值是区分 2A、2B 和 2M-A1 轴的关键检查；当前结果缺失，不能仅凭 AI 模型确定亚型。
建议优先补充：VWF_MULTIMER, VWF_CB, VWF_FVIIIB。

## Candidate subtypes

type_2_candidate

## Supporting evidence

- CASE_T1_009
- CASE_T1_009:alphagenome-full-profile
- CASE_T1_009:VWF_R1205H:boltz2-functional-panel
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
- Observation/computational-result-d2904a5a3e5593769b37
- Observation/computational-result-2c67408296b82d3d6ec3
- Observation/vwf-domain-1205-Arg-His

## Missing information

- RIPA
- VWF_MULTIMER
- VWF_CB
- VWF_FVIIIB
- VWF_PP
- DDAVP_0_1_4H
- expert_variant_classification
- population_frequency
- pathogenicity_prediction
- variant_specific_literature

## Limitations

- ensembl_variant_recoder returned error; this is not benign evidence.
- coordinate_dependent_tools returned error; this is not benign evidence.
- clingen_erepo returned error; this is not benign evidence.
- pubmed_eutils returned error; this is not benign evidence.
- local_clingen_snapshot returned not_found; this is not benign evidence.
- ClinGen/ACMG type-2 rule provider is not enabled in V0 offline mode.
- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- The evidence cannot distinguish type 2 subtypes or exclude qualitative VWF defects.

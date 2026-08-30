# VWD Research Report — CASE_T2B_004

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

### 1. 临床与实验室表型
VWF:Ag=24.2, VWF:Act=15.0, FVIII:C=28.4, 血小板计数=389.0。
VWF:Act/VWF:Ag 比值为 0.620，低于 0.70，提示存在不成比例的功能性 VWF 缺陷。

### 2. 变异与功能域
- c.3922C>T / p.R1308C（VWFA1）。

### 3. AI 机制证据
- alphagenome_full_profile: AlphaGenome complete-profile selected-track scores; top absolute views: ag_chip_histone_active_abs_max=1.115e+04, ag_chip_tf_active_abs_max=1850, ag_cage_active_abs_max=620.5, ag_dnase_active_abs_max=32.23, ag_procap_active_abs_max=20.94.
- boltz2_functional_panel: Assay-matched Boltz-2 structural deltas: a1_aim_autoinhibition_context delta_vs_WT=-0.1152; a1_gpiba_forced_binding delta_vs_WT=-0.00055; a1_heparan_sulfate_binding delta_vs_WT=0.08599.
- md_targeted_panel: Targeted MD features: AIM_all_contacts_frames=501, AIM_all_contacts_mean=161.3, AIM_all_contacts_first0_5=167.4, AIM_all_contacts_mid20_30=163.3, AIM_all_contacts_tail40_50=165.8, AIM_all_contacts_final=175, AIM_all_contacts_min=107, AIM_all_contacts_max=192, AIM_all_contacts_tail_minus_first=-1.56, N_AIM_contacts_frames=501, N_AIM_contacts_mean=106.6, N_AIM_contacts_first0_5=116.4, N_AIM_contacts_mid20_30=108.7, N_AIM_contacts_tail40_50=107, N_AIM_contacts_final=121, N_AIM_contacts_min=56, N_AIM_contacts_max=142, N_AIM_contacts_tail_minus_first=-9.412, C_AIM_contacts_frames=501, C_AIM_contacts_mean=55.86, C_AIM_contacts_first0_5=51.43, C_AIM_contacts_mid20_30=54.64, C_AIM_contacts_tail40_50=63.23, C_AIM_contacts_final=54, C_AIM_contacts_min=32, C_AIM_contacts_max=80, C_AIM_contacts_tail_minus_first=11.8, AIM_all_contacts_delta_tail_vs_WT=53.08, N_AIM_contacts_delta_tail_vs_WT=47.45, C_AIM_contacts_delta_tail_vs_WT=9.723.

### 4. 机制解释
该变异位于 A1 功能域，A1 参与 GPIb 结合并受 AIM 自抑制调控；该区域的错义变异可能通过改变自抑制界面、A1 暴露或表面电荷，影响血小板结合功能。
MD 结果提示 AIM-A1 接触动力学发生改变，可作为自抑制释放或 A1 功能面暴露的动态证据。
静态 Boltz 结果提示相关结构轴发生扰动，但静态置信度不等于结合自由能，需与 MD 和功能实验联合解释。
实验室表型与功能性 VWF 缺陷方向一致；若同时存在高危分子量多聚体缺失，更支持 2A/2B 样机制。

### 5. 分型鉴别与不确定性
HGMD/ClinVar 的致病性评级本身不能直接给出 VWD 亚型方向；即使标注为致病或 uncertain，也必须结合实验室检查、出血表型和机制证据综合判断。
当前倾向：type_2_vwd（moderate）; type_2B（low）; type_2M（low）。
VWF 多聚体分析和 VWF:CB/Ag 比值是区分 2A、2B 和 2M-A1 轴的关键检查；当前结果缺失，不能仅凭 AI 模型确定亚型。
建议优先补充：VWF_MULTIMER, VWF_CB, VWF_FVIIIB。

## Candidate subtypes

type_2_candidate

## Supporting evidence

- CASE_T2B_004
- CASE_T2B_004:alphagenome-full-profile
- CASE_T2B_004:VWF_R1308C:boltz2-functional-panel
- CASE_T2B_004:R1308C:targeted-md
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
- Observation/computational-result-3062d286874ca7242d28
- Observation/computational-result-c8deedc48cbe894e3d56
- Observation/computational-result-31cdeaf0f630e6d01051
- Observation/vwf-domain-1308-R-C
- Observation/local-clingen-0-289

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
- ClinGen/ACMG type-2 rule provider is not enabled in V0 offline mode.
- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- The evidence cannot distinguish type 2 subtypes or exclude qualitative VWF defects.

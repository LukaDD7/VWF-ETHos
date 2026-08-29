# VWD Research Report — CASE_T2B_003

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

### 1. 临床与实验室表型
VWF:Ag=40.2, VWF:Act=14.6, FVIII:C=128.91, 血小板计数=73.0。
VWF:Act/VWF:Ag 比值为 0.363，低于 0.70，提示存在不成比例的功能性 VWF 缺陷。
临床背景中存在血小板减少线索，需与 2B/血小板型 VWD 鉴别。

### 2. 变异与功能域
- c.3929C>T / p.S1310F（VWFA1）。

### 3. AI 机制证据
- alphagenome_full_profile: AlphaGenome complete-profile selected-track scores; top absolute views: ag_chip_histone_active_abs_max=1.111e+04, ag_chip_tf_active_abs_max=1838, ag_cage_active_abs_max=612, ag_dnase_active_abs_max=31.83, ag_procap_active_abs_max=19.79.
- boltz2_functional_panel: Assay-matched Boltz-2 structural deltas: a1_aim_autoinhibition_context delta_vs_WT=-0.03908; a1_gpiba_forced_binding delta_vs_WT=-0.04199; a1_heparan_sulfate_binding delta_vs_WT=-0.01172.
- md_targeted_panel: Targeted MD features: AIM_all_contacts_frames=501, AIM_all_contacts_mean=130.9, AIM_all_contacts_first0_5=186.8, AIM_all_contacts_mid20_30=125.5, AIM_all_contacts_tail40_50=121.7, AIM_all_contacts_final=129, AIM_all_contacts_min=90, AIM_all_contacts_max=222, AIM_all_contacts_tail_minus_first=-65.18, N_AIM_contacts_frames=501, N_AIM_contacts_mean=75.69, N_AIM_contacts_first0_5=130.8, N_AIM_contacts_mid20_30=69.07, N_AIM_contacts_tail40_50=71.98, N_AIM_contacts_final=67, N_AIM_contacts_min=36, N_AIM_contacts_max=161, N_AIM_contacts_tail_minus_first=-58.82, C_AIM_contacts_frames=501, C_AIM_contacts_mean=55.36, C_AIM_contacts_first0_5=56.04, C_AIM_contacts_mid20_30=57.2, C_AIM_contacts_tail40_50=49.69, C_AIM_contacts_final=62, C_AIM_contacts_min=35, C_AIM_contacts_max=75, C_AIM_contacts_tail_minus_first=-6.346, AIM_all_contacts_delta_tail_vs_WT=8.911, N_AIM_contacts_delta_tail_vs_WT=12.41, C_AIM_contacts_delta_tail_vs_WT=-3.812.

### 4. 机制解释
该变异位于 A1 功能域，A1 参与 GPIb 结合并受 AIM 自抑制调控；该区域的错义变异可能通过改变自抑制界面、A1 暴露或表面电荷，影响血小板结合功能。
MD 结果提示 AIM-A1 接触动力学发生改变，可作为自抑制释放或 A1 功能面暴露的动态证据。
静态 Boltz 结果提示相关结构轴发生扰动，但静态置信度不等于结合自由能，需与 MD 和功能实验联合解释。
实验室表型与功能性 VWF 缺陷方向一致；若同时存在高危分子量多聚体缺失，更支持 2A/2B 样机制。

### 5. 分型鉴别与不确定性
HGMD/ClinVar 的致病性评级本身不能直接给出 VWD 亚型方向；即使标注为致病或 uncertain，也必须结合实验室检查、出血表型和机制证据综合判断。
当前倾向：type_2_vwd（moderate）。
VWF 多聚体分析和 VWF:CB/Ag 比值是区分 2A、2B 和 2M-A1 轴的关键检查；当前结果缺失，不能仅凭 AI 模型确定亚型。
建议优先补充：RIPA, VWF_MULTIMER。

## Candidate subtypes

type_2_candidate, platelet_type_vwd_candidate

## Supporting evidence

- CASE_T2B_003
- CASE_T2B_003:alphagenome-full-profile
- CASE_T2B_003:VWF_S1310F:boltz2-functional-panel
- CASE_T2B_003:S1310F:targeted-md

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

- ClinGen/ACMG type-2 rule provider is not enabled in V0 offline mode.
- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- The evidence cannot distinguish type 2 subtypes or exclude qualitative VWF defects.

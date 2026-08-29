# VWD Research Report — CASE_T2B_006

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

### 1. 临床与实验室表型
VWF:Ag=44.9, VWF:Act=33.5, FVIII:C=58.5, 血小板计数=299.0。
VWF:Act/VWF:Ag 比值未显示明显不成比例的功能下降。

### 2. 变异与功能域
- c.4021C>T / p.Arg1341Trp（VWFA1）。

### 3. AI 机制证据
- alphagenome_full_profile: AlphaGenome complete-profile selected-track scores; top absolute views: ag_chip_histone_active_abs_max=1.108e+04, ag_chip_tf_active_abs_max=1776, ag_cage_active_abs_max=595.7, ag_dnase_active_abs_max=35.55, ag_procap_active_abs_max=19.52.
- boltz2_functional_panel: Assay-matched Boltz-2 structural deltas: a1_aim_autoinhibition_context delta_vs_WT=-0.1366; a1_gpiba_forced_binding delta_vs_WT=-0.02476; a1_heparan_sulfate_binding delta_vs_WT=0.00728.
- md_targeted_panel: Targeted MD features: AIM_all_contacts_frames=501, AIM_all_contacts_mean=104.8, AIM_all_contacts_first0_5=117.4, AIM_all_contacts_mid20_30=102.2, AIM_all_contacts_tail40_50=105.4, AIM_all_contacts_final=99, AIM_all_contacts_min=62, AIM_all_contacts_max=149, AIM_all_contacts_tail_minus_first=-12.05, N_AIM_contacts_frames=501, N_AIM_contacts_mean=53.14, N_AIM_contacts_first0_5=58.24, N_AIM_contacts_mid20_30=55.54, N_AIM_contacts_tail40_50=54.79, N_AIM_contacts_final=50, N_AIM_contacts_min=12, N_AIM_contacts_max=93, N_AIM_contacts_tail_minus_first=-3.443, C_AIM_contacts_frames=501, C_AIM_contacts_mean=51.77, C_AIM_contacts_first0_5=59.2, C_AIM_contacts_mid20_30=46.69, C_AIM_contacts_tail40_50=50.65, C_AIM_contacts_final=49, C_AIM_contacts_min=30, C_AIM_contacts_max=74, C_AIM_contacts_tail_minus_first=-8.543, AIM_all_contacts_delta_tail_vs_WT=-7.366, N_AIM_contacts_delta_tail_vs_WT=-4.782, C_AIM_contacts_delta_tail_vs_WT=-2.852.

### 4. 机制解释
该变异位于 A1 功能域，A1 参与 GPIb 结合并受 AIM 自抑制调控；该区域的错义变异可能通过改变自抑制界面、A1 暴露或表面电荷，影响血小板结合功能。
MD 结果提示 AIM-A1 接触动力学发生改变，可作为自抑制释放或 A1 功能面暴露的动态证据。
静态 Boltz 结果提示相关结构轴发生扰动，但静态置信度不等于结合自由能，需与 MD 和功能实验联合解释。

### 5. 分型鉴别与不确定性
HGMD/ClinVar 的致病性评级本身不能直接给出 VWD 亚型方向；即使标注为致病或 uncertain，也必须结合实验室检查、出血表型和机制证据综合判断。
当前倾向：type_1_or_low_vwf（low）。
建议优先补充：VWF_PP, DDAVP_0_1_4H。

## Candidate subtypes

type_1_candidate_provisional

## Supporting evidence

- CASE_T2B_006
- CASE_T2B_006:alphagenome-full-profile
- CASE_T2B_006:VWF_R1341W:boltz2-functional-panel
- CASE_T2B_006:R1341W:targeted-md

## Missing information

- expert_variant_classification
- population_frequency
- pathogenicity_prediction
- variant_specific_literature
- ABO_blood_group_genotype
- VWF_clearance_context

## Limitations

- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- The remaining uncertainty is due to missing evidence rather than contradictory evidence.

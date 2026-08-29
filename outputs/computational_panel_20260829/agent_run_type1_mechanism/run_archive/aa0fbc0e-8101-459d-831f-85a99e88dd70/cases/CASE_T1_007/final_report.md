# VWD Research Report — CASE_T1_007

- Status: low
- Abstention: True
- Expert review required: True

## Opinion

### 1. 临床与实验室表型
VWF:Ag=48.8, VWF:Act=33.6, FVIII:C=未提供, 血小板计数=未提供。
VWF:Act/VWF:Ag 比值为 0.689，低于 0.70，提示存在不成比例的功能性 VWF 缺陷。

### 2. 变异与功能域
- c.3379+1G>A / 无蛋白改变（未能映射到已知功能域）。

### 3. AI 机制证据
- alphagenome_full_profile: AlphaGenome complete-profile selected-track scores; top absolute views: ag_chip_histone_active_abs_max=9834, ag_chip_tf_active_abs_max=1253, ag_cage_active_abs_max=159.1, ag_dnase_active_abs_max=32.22, ag_rna_seq_active_abs_max=16.18.

### 4. 机制解释
实验室表型与功能性 VWF 缺陷方向一致；若同时存在高危分子量多聚体缺失，更支持 2A/2B 样机制。

### 5. 分型鉴别与不确定性
HGMD/ClinVar 的致病性评级本身不能直接给出 VWD 亚型方向；即使标注为致病或 uncertain，也必须结合实验室检查、出血表型和机制证据综合判断。
当前倾向：type_2_vwd（moderate）。
VWF 多聚体分析和 VWF:CB/Ag 比值是区分 2A、2B 和 2M-A1 轴的关键检查；当前结果缺失，不能仅凭 AI 模型确定亚型。
建议优先补充：VWF_MULTIMER, VWF_CB。

## Candidate subtypes

type_2_candidate

## Supporting evidence

- CASE_T1_007
- CASE_T1_007:alphagenome-full-profile

## Missing information

- FVIII:C
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
- Missing critical fields: FVIII:C.
- Second-level tests were explicitly unavailable; absence of results is not evidence against type 2 VWD.
- Coexisting hemophilia A confounds FVIII:C; FVIII:C must not be used to support or refute VWD type 2N.
- The evidence cannot distinguish type 2 subtypes or exclude qualitative VWF defects.

# Type2M LOF MD fast validation

Analysis dir: `/Users/lucachangretta/Agent_Workspace/VWF-ETHos/output/type2m_lof_md_analysis_2026-06-29`
Closed-state contact-loss threshold: `20.0` contacts -> `md_face_destab_score=1.0`

## Feature Summary

- 7A6O closed-state completed feature rows: `10`
- Rows crossing MD_LOF threshold (`md_face_destab_score >= 1.0`): `4`

| variant   | labels   |   md_face_destab_score |   md_closed_aim_contact_loss |   AIM_all_contacts_first0_5 |   AIM_all_contacts_tail40_50 |
|:----------|:---------|-----------------------:|-----------------------------:|----------------------------:|-----------------------------:|
| L1276R    | 2M       |                  2.266 |                      45.325  |                     144.177 |                      98.8515 |
| G1324R    | 2M       |                  2.186 |                      43.7202 |                     139.255 |                      95.5347 |
| L1276P    | 2M       |                  1.701 |                      34.0118 |                     146.804 |                     112.792  |
| A1377V    | 2M       |                  1.283 |                      25.6544 |                     118.902 |                      93.2475 |
| Y1321C    | 2M       |                  0.976 |                      19.5207 |                     133.412 |                     113.891  |
| Y1321W    | 2M       |                  0.803 |                      16.0588 |                     114.059 |                      98      |
| G1324A    | 2M       |                  0.716 |                      14.3209 |                     126.588 |                     112.267  |
| R1315G    | 2M       |                  0.614 |                      12.2879 |                     110.922 |                      98.6337 |
| R1315P    | 2M       |                  0.189 |                       3.7894 |                     109.275 |                     105.485  |
| D1302G    | 2M       |                 -0.181 |                      -3.6253 |                     106.157 |                     109.782  |

## Recall Delta

| label   |   n_base |   correct_base |   recall_base |   uncertain_base |   n_fast |   correct_fast |   recall_fast |   uncertain_fast |   delta_correct |   delta_recall |   delta_uncertain |
|:--------|---------:|---------------:|--------------:|-----------------:|---------:|---------------:|--------------:|-----------------:|----------------:|---------------:|------------------:|
| 2A      |      118 |             70 |         0.593 |               18 |      118 |             70 |         0.593 |               18 |               0 |          0     |                 0 |
| 2B      |       38 |             29 |         0.763 |                6 |       38 |             29 |         0.763 |                6 |               0 |          0     |                 0 |
| 2M      |       53 |             21 |         0.396 |               26 |       53 |             23 |         0.434 |               24 |               2 |          0.038 |                -2 |
| 2N      |       16 |             12 |         0.75  |                0 |       16 |             12 |         0.75  |                0 |               0 |          0     |                 0 |
| ALL     |      225 |            132 |         0.587 |               50 |      225 |            134 |         0.596 |               48 |               2 |          0.009 |                -2 |

## Prediction Changes On Newly Completed MD Variants

| aa_change   | true_label   | pred_baseline   | pred_fast_md   |   md_face_destab_score |   md_closed_aim_contact_loss | prediction_changed   | rescued_to_true_label   | strong_md_not_rescued   |
|:------------|:-------------|:----------------|:---------------|-----------------------:|-----------------------------:|:---------------------|:------------------------|:------------------------|
| L1276P      | 2M           | uncertain       | 2M             |                  1.701 |                      34.0118 | True                 | True                    | False                   |
| L1276R      | 2M           | uncertain       | 2M             |                  2.266 |                      45.325  | True                 | True                    | False                   |
| A1377V      | 2M           | 2B              | 2B             |                  1.283 |                      25.6544 | False                | False                   | True                    |
| D1302G      | 2M           | uncertain       | uncertain      |                 -0.181 |                      -3.6253 | False                | False                   | False                   |
| G1324A      | 2M           | 2B              | 2B             |                  0.716 |                      14.3209 | False                | False                   | False                   |
| G1324R      | 2M           | 2B              | 2B             |                  2.186 |                      43.7202 | False                | False                   | True                    |
| R1315G      | 2M           | uncertain       | uncertain      |                  0.614 |                      12.2879 | False                | False                   | False                   |
| R1315P      | 2M           | uncertain       | uncertain      |                  0.189 |                       3.7894 | False                | False                   | False                   |
| Y1321C      | 2M           | 2B              | 2B             |                  0.976 |                      19.5207 | False                | False                   | False                   |
| Y1321W      | 2M           | 2B              | 2B             |                  0.803 |                      16.0588 | False                | False                   | False                   |

Strong MD LOF but not rescued: `2`. These are rule-order/2B-prior review candidates.

## A1-GPIb Status

A1-GPIb MD is reported here for QC only. RMSD alone is not fed into the classifier; a directional GPIb interface/contact feature should be extracted before using it as LOF evidence.

| status                  | labels   |   n |
|:------------------------|:---------|----:|
| complete                | 2M       |   9 |
| pending                 | 2A|2B    |   1 |
| pending                 | 2B       |   1 |
| pending                 | 2M       |  21 |
| running_or_checkpointed | 2M       |   7 |

# Type2M LOF MD fast validation

Analysis dir: `/Users/lucachangretta/Agent_Workspace/VWF-ETHos/output/type2m_lof_md_analysis_2026-06-29`
Closed-state contact-loss threshold: `20.0` contacts -> `md_face_destab_score=1.0`

## Feature Summary

- 7A6O closed-state completed feature rows: `22`
- Rows crossing MD_LOF threshold (`md_face_destab_score >= 1.0`): `10`

| variant   | labels   |   md_face_destab_score |   md_closed_aim_contact_loss |   AIM_all_contacts_first0_5 |   AIM_all_contacts_tail40_50 |
|:----------|:---------|-----------------------:|-----------------------------:|----------------------------:|-----------------------------:|
| L1276R    | 2M       |                  2.266 |                      45.325  |                    144.177  |                      98.8515 |
| G1324R    | 2M       |                  2.186 |                      43.7202 |                    139.255  |                      95.5347 |
| A1355D    | 2M       |                  2.182 |                      43.6469 |                    131.667  |                      88.0198 |
| S1387I    | 2M       |                  2.164 |                      43.2722 |                    160.51   |                     117.238  |
| L1383R    | 2M       |                  1.889 |                      37.7798 |                    165.235  |                     127.455  |
| Y1363C    | 2M       |                  1.802 |                      36.0415 |                    154.804  |                     118.762  |
| L1276P    | 2M       |                  1.701 |                      34.0118 |                    146.804  |                     112.792  |
| A1377V    | 2M       |                  1.283 |                      25.6544 |                    118.902  |                      93.2475 |
| F1369I    | 2M       |                  1.034 |                      20.6725 |                    127.078  |                     106.406  |
| K1362T    | 2M       |                  1.005 |                      20.0918 |                    106.725  |                      86.6337 |
| Y1321C    | 2M       |                  0.976 |                      19.5207 |                    133.412  |                     113.891  |
| D1373Y    | 2M       |                  0.93  |                      18.6038 |                    133.02   |                     114.416  |
| Y1321W    | 2M       |                  0.803 |                      16.0588 |                    114.059  |                      98      |
| G1324A    | 2M       |                  0.716 |                      14.3209 |                    126.588  |                     112.267  |
| R1315G    | 2M       |                  0.614 |                      12.2879 |                    110.922  |                      98.6337 |
| V1360A    | 2M       |                  0.365 |                       7.2953 |                    119.177  |                     111.881  |
| L1361W    | 2M       |                  0.265 |                       5.2957 |                    118.137  |                     112.842  |
| R1315P    | 2M       |                  0.189 |                       3.7894 |                    109.275  |                     105.485  |
| H1322P    | 2M       |                  0.044 |                       0.875  |                    121.627  |                     120.752  |
| D1302G    | 2M       |                 -0.181 |                      -3.6253 |                    106.157  |                     109.782  |
| L1382P    | 2M       |                 -0.612 |                     -12.246  |                    118.843  |                     131.089  |
| E1359K    | 2M       |                 -0.725 |                     -14.5015 |                     99.6471 |                     114.148  |

## Recall Delta

| label   |   n_base |   correct_base |   recall_base |   uncertain_base |   n_fast |   correct_fast |   recall_fast |   uncertain_fast |   delta_correct |   delta_recall |   delta_uncertain |
|:--------|---------:|---------------:|--------------:|-----------------:|---------:|---------------:|--------------:|-----------------:|----------------:|---------------:|------------------:|
| 2A      |      118 |             70 |         0.593 |               18 |      118 |             70 |         0.593 |               18 |               0 |          0     |                 0 |
| 2B      |       38 |             15 |         0.395 |               19 |       38 |             15 |         0.395 |               19 |               0 |          0     |                 0 |
| 2M      |       53 |             21 |         0.396 |               32 |       53 |             31 |         0.585 |               22 |              10 |          0.189 |               -10 |
| 2N      |       16 |             12 |         0.75  |                0 |       16 |             12 |         0.75  |                0 |               0 |          0     |                 0 |
| ALL     |      225 |            118 |         0.524 |               69 |      225 |            128 |         0.569 |               59 |              10 |          0.045 |               -10 |

## Prediction Changes On Newly Completed MD Variants

| aa_change   | true_label   | pred_baseline   | pred_fast_md   |   md_face_destab_score |   md_closed_aim_contact_loss | prediction_changed   | rescued_to_true_label   | strong_md_not_rescued   |
|:------------|:-------------|:----------------|:---------------|-----------------------:|-----------------------------:|:---------------------|:------------------------|:------------------------|
| A1355D      | 2M           | uncertain       | 2M             |                  2.182 |                      43.6469 | True                 | True                    | False                   |
| A1377V      | 2M           | uncertain       | 2M             |                  1.283 |                      25.6544 | True                 | True                    | False                   |
| F1369I      | 2M           | uncertain       | 2M             |                  1.034 |                      20.6725 | True                 | True                    | False                   |
| G1324R      | 2M           | uncertain       | 2M             |                  2.186 |                      43.7202 | True                 | True                    | False                   |
| K1362T      | 2M           | uncertain       | 2M             |                  1.005 |                      20.0918 | True                 | True                    | False                   |
| L1276P      | 2M           | uncertain       | 2M             |                  1.701 |                      34.0118 | True                 | True                    | False                   |
| L1276R      | 2M           | uncertain       | 2M             |                  2.266 |                      45.325  | True                 | True                    | False                   |
| L1383R      | 2M           | uncertain       | 2M             |                  1.889 |                      37.7798 | True                 | True                    | False                   |
| S1387I      | 2M           | uncertain       | 2M             |                  2.164 |                      43.2722 | True                 | True                    | False                   |
| Y1363C      | 2M           | uncertain       | 2M             |                  1.802 |                      36.0415 | True                 | True                    | False                   |
| D1302G      | 2M           | uncertain       | uncertain      |                 -0.181 |                      -3.6253 | False                | False                   | False                   |
| D1373Y      | 2M           | uncertain       | uncertain      |                  0.93  |                      18.6038 | False                | False                   | False                   |
| E1359K      | 2M           | uncertain       | uncertain      |                 -0.725 |                     -14.5015 | False                | False                   | False                   |
| G1324A      | 2M           | uncertain       | uncertain      |                  0.716 |                      14.3209 | False                | False                   | False                   |
| H1322P      | 2M           | uncertain       | uncertain      |                  0.044 |                       0.875  | False                | False                   | False                   |
| L1361W      | 2M           | uncertain       | uncertain      |                  0.265 |                       5.2957 | False                | False                   | False                   |
| L1382P      | 2M           | uncertain       | uncertain      |                 -0.612 |                     -12.246  | False                | False                   | False                   |
| R1315G      | 2M           | uncertain       | uncertain      |                  0.614 |                      12.2879 | False                | False                   | False                   |
| R1315P      | 2M           | uncertain       | uncertain      |                  0.189 |                       3.7894 | False                | False                   | False                   |
| V1360A      | 2M           | uncertain       | uncertain      |                  0.365 |                       7.2953 | False                | False                   | False                   |
| Y1321C      | 2M           | uncertain       | uncertain      |                  0.976 |                      19.5207 | False                | False                   | False                   |
| Y1321W      | 2M           | uncertain       | uncertain      |                  0.803 |                      16.0588 | False                | False                   | False                   |

## A1-GPIb Status

A1-GPIb MD is reported here for QC only. RMSD alone is not fed into the classifier; a directional GPIb interface/contact feature should be extracted before using it as LOF evidence.

| status                  | labels   |   n |
|:------------------------|:---------|----:|
| complete                | 2M       |  25 |
| pending                 | 2A|2B    |   1 |
| pending                 | 2B       |   1 |
| pending                 | 2M       |   5 |
| running_or_checkpointed | 2M       |   7 |

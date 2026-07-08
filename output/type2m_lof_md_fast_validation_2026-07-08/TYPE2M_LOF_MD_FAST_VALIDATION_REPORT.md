# Type2M LOF MD fast validation

Analysis dir: `/media/luzhenyang/project/alphagenome/alphagenome/VWF_ErTongYiyuan/output/type2m_lof_md_analysis_2026-06-29`
Closed-state contact-loss threshold: `20.0` contacts -> `md_face_destab_score=1.0`

## Feature Summary

- 7A6O closed-state completed feature rows: `31`
- Rows crossing MD_LOF threshold (`md_face_destab_score >= 1.0`): `16`

| variant   | labels   |   md_face_destab_score |   md_closed_aim_contact_loss |   AIM_all_contacts_first0_5 |   AIM_all_contacts_tail40_50 |
|:----------|:---------|-----------------------:|-----------------------------:|----------------------------:|-----------------------------:|
| I1416T    | 2M       |                  2.504 |                      50.0724 |                    153.686  |                     103.614  |
| L1276R    | 2M       |                  2.266 |                      45.325  |                    144.177  |                      98.8515 |
| G1324R    | 2M       |                  2.186 |                      43.7202 |                    139.255  |                      95.5347 |
| A1355D    | 2M       |                  2.182 |                      43.6469 |                    131.667  |                      88.0198 |
| S1387I    | 2M       |                  2.164 |                      43.2722 |                    160.51   |                     117.238  |
| V1439M    | 2M       |                  1.906 |                      38.1275 |                    158.118  |                     119.99   |
| L1383R    | 2M       |                  1.889 |                      37.7798 |                    165.235  |                     127.455  |
| Y1363C    | 2M       |                  1.802 |                      36.0415 |                    154.804  |                     118.762  |
| R1399C    | 2M       |                  1.758 |                      35.1516 |                    141.686  |                     106.535  |
| L1276P    | 2M       |                  1.701 |                      34.0118 |                    146.804  |                     112.792  |
| S1394F    | 2M       |                  1.529 |                      30.5754 |                    156.882  |                     126.307  |
| R1426C    | 2M       |                  1.456 |                      29.1231 |                    151.569  |                     122.445  |
| A1377V    | 2M       |                  1.283 |                      25.6544 |                    118.902  |                      93.2475 |
| Q1402P    | 2M       |                  1.091 |                      21.8216 |                    149.02   |                     127.198  |
| F1369I    | 2M       |                  1.034 |                      20.6725 |                    127.078  |                     106.406  |
| K1362T    | 2M       |                  1.005 |                      20.0918 |                    106.725  |                      86.6337 |
| Y1321C    | 2M       |                  0.976 |                      19.5207 |                    133.412  |                     113.891  |
| D1373Y    | 2M       |                  0.93  |                      18.6038 |                    133.02   |                     114.416  |
| Y1321W    | 2M       |                  0.803 |                      16.0588 |                    114.059  |                      98      |
| G1324A    | 2M       |                  0.716 |                      14.3209 |                    126.588  |                     112.267  |
| R1315G    | 2M       |                  0.614 |                      12.2879 |                    110.922  |                      98.6337 |
| I1425F    | 2M       |                  0.431 |                       8.6166 |                    133.726  |                     125.109  |
| V1360A    | 2M       |                  0.365 |                       7.2953 |                    119.177  |                     111.881  |
| L1361W    | 2M       |                  0.265 |                       5.2957 |                    118.137  |                     112.842  |
| R1315P    | 2M       |                  0.189 |                       3.7894 |                    109.275  |                     105.485  |
| H1322P    | 2M       |                  0.044 |                       0.875  |                    121.627  |                     120.752  |
| D1302G    | 2M       |                 -0.181 |                      -3.6253 |                    106.157  |                     109.782  |
| P1413R    | 2M       |                 -0.239 |                      -4.7896 |                    131.745  |                     136.535  |
| L1382P    | 2M       |                 -0.612 |                     -12.246  |                    118.843  |                     131.089  |
| V1409F    | 2M       |                 -0.721 |                     -14.4252 |                     94.9412 |                     109.366  |
| E1359K    | 2M       |                 -0.725 |                     -14.5015 |                     99.6471 |                     114.148  |

## Recall Delta

| label   |   n_base |   correct_base |   recall_base |   uncertain_base |   n_fast |   correct_fast |   recall_fast |   uncertain_fast |   delta_correct |   delta_recall |   delta_uncertain |
|:--------|---------:|---------------:|--------------:|-----------------:|---------:|---------------:|--------------:|-----------------:|----------------:|---------------:|------------------:|
| 2A      |      118 |             70 |         0.593 |               18 |      118 |             70 |         0.593 |               18 |               0 |          0     |                 0 |
| 2B      |       38 |             15 |         0.395 |               19 |       38 |             15 |         0.395 |               19 |               0 |          0     |                 0 |
| 2M      |       53 |             21 |         0.396 |               32 |       53 |             37 |         0.698 |               16 |              16 |          0.302 |               -16 |
| 2N      |       16 |             12 |         0.75  |                0 |       16 |             12 |         0.75  |                0 |               0 |          0     |                 0 |
| ALL     |      225 |            118 |         0.524 |               69 |      225 |            134 |         0.596 |               53 |              16 |          0.072 |               -16 |

## Prediction Changes On Newly Completed MD Variants

| aa_change   | true_label   | pred_baseline   | pred_fast_md   |   md_face_destab_score |   md_closed_aim_contact_loss | prediction_changed   | rescued_to_true_label   | strong_md_not_rescued   |
|:------------|:-------------|:----------------|:---------------|-----------------------:|-----------------------------:|:---------------------|:------------------------|:------------------------|
| A1355D      | 2M           | uncertain       | 2M             |                  2.182 |                      43.6469 | True                 | True                    | False                   |
| A1377V      | 2M           | uncertain       | 2M             |                  1.283 |                      25.6544 | True                 | True                    | False                   |
| F1369I      | 2M           | uncertain       | 2M             |                  1.034 |                      20.6725 | True                 | True                    | False                   |
| G1324R      | 2M           | uncertain       | 2M             |                  2.186 |                      43.7202 | True                 | True                    | False                   |
| I1416T      | 2M           | uncertain       | 2M             |                  2.504 |                      50.0724 | True                 | True                    | False                   |
| K1362T      | 2M           | uncertain       | 2M             |                  1.005 |                      20.0918 | True                 | True                    | False                   |
| L1276P      | 2M           | uncertain       | 2M             |                  1.701 |                      34.0118 | True                 | True                    | False                   |
| L1276R      | 2M           | uncertain       | 2M             |                  2.266 |                      45.325  | True                 | True                    | False                   |
| L1383R      | 2M           | uncertain       | 2M             |                  1.889 |                      37.7798 | True                 | True                    | False                   |
| Q1402P      | 2M           | uncertain       | 2M             |                  1.091 |                      21.8216 | True                 | True                    | False                   |
| R1399C      | 2M           | uncertain       | 2M             |                  1.758 |                      35.1516 | True                 | True                    | False                   |
| R1426C      | 2M           | uncertain       | 2M             |                  1.456 |                      29.1231 | True                 | True                    | False                   |
| S1387I      | 2M           | uncertain       | 2M             |                  2.164 |                      43.2722 | True                 | True                    | False                   |
| S1394F      | 2M           | uncertain       | 2M             |                  1.529 |                      30.5754 | True                 | True                    | False                   |
| V1439M      | 2M           | uncertain       | 2M             |                  1.906 |                      38.1275 | True                 | True                    | False                   |
| Y1363C      | 2M           | uncertain       | 2M             |                  1.802 |                      36.0415 | True                 | True                    | False                   |
| D1302G      | 2M           | uncertain       | uncertain      |                 -0.181 |                      -3.6253 | False                | False                   | False                   |
| D1373Y      | 2M           | uncertain       | uncertain      |                  0.93  |                      18.6038 | False                | False                   | False                   |
| E1359K      | 2M           | uncertain       | uncertain      |                 -0.725 |                     -14.5015 | False                | False                   | False                   |
| G1324A      | 2M           | uncertain       | uncertain      |                  0.716 |                      14.3209 | False                | False                   | False                   |
| H1322P      | 2M           | uncertain       | uncertain      |                  0.044 |                       0.875  | False                | False                   | False                   |
| I1425F      | 2M           | uncertain       | uncertain      |                  0.431 |                       8.6166 | False                | False                   | False                   |
| L1361W      | 2M           | uncertain       | uncertain      |                  0.265 |                       5.2957 | False                | False                   | False                   |
| L1382P      | 2M           | uncertain       | uncertain      |                 -0.612 |                     -12.246  | False                | False                   | False                   |
| P1413R      | 2M           | uncertain       | uncertain      |                 -0.239 |                      -4.7896 | False                | False                   | False                   |
| R1315G      | 2M           | uncertain       | uncertain      |                  0.614 |                      12.2879 | False                | False                   | False                   |
| R1315P      | 2M           | uncertain       | uncertain      |                  0.189 |                       3.7894 | False                | False                   | False                   |
| V1360A      | 2M           | uncertain       | uncertain      |                  0.365 |                       7.2953 | False                | False                   | False                   |
| V1409F      | 2M           | 2M              | 2M             |                 -0.721 |                     -14.4252 | False                | False                   | False                   |
| Y1321C      | 2M           | uncertain       | uncertain      |                  0.976 |                      19.5207 | False                | False                   | False                   |
| Y1321W      | 2M           | uncertain       | uncertain      |                  0.803 |                      16.0588 | False                | False                   | False                   |

## A1-GPIb Status

A1-GPIb MD is reported here for QC only. RMSD alone is not fed into the classifier; a directional GPIb interface/contact feature should be extracted before using it as LOF evidence.

| status   | labels   |   n |
|:---------|:---------|----:|
| complete | 2A|2B    |   1 |
| complete | 2B       |   1 |
| complete | 2M       |  35 |
| pending  | 2M       |   2 |

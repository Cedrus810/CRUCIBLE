# Gate G1 评估报告（MVP-1）

频率与几何指标只报告、不设门槛。

| 指标 | 说明 |
|---|---|
| freq_rmse_cm | FF Hessian 频率 vs 参考 Hessian 频率的 RMSE（cm⁻¹） |
| kabsch_rmsd_pm | FF 优化几何 vs 参考几何的 Kabsch RMSD（pm） |
| bond_mae_pm | 优化后键长 MAE（pm） |
| angle_mae_deg | 优化后键角 MAE（deg） |


## Hessian: gfn2

| method | 指标 | 中位数 | 最大值 |
|---|---|---|---|
| gaff2 | freq RMSE (cm⁻¹) | 109.8 | 155.8 |
| gaff2 | Kabsch RMSD (pm) | 4.28 | 45.3 |
| gaff2 | bond MAE (pm) | 1.368 | 1.922 |
| gaff2 | angle MAE (deg) | 0.8657 | 2.871 |
| gaff2 | max class spread | 0 | 0 |
| mseminario | freq RMSE (cm⁻¹) | 109 | 164.5 |
| mseminario | Kabsch RMSD (pm) | 4.347 | 45.16 |
| mseminario | bond MAE (pm) | 0.348 | 1.048 |
| mseminario | angle MAE (deg) | 0.5612 | 2.659 |
| mseminario | max class spread | 0.07054 | 0.3028 |
| m2seminario | freq RMSE (cm⁻¹) | 187.8 | 235.9 |
| m2seminario | Kabsch RMSD (pm) | 4.736 | 45.2 |
| m2seminario | bond MAE (pm) | 0.3609 | 1.045 |
| m2seminario | angle MAE (deg) | 0.6355 | 3.568 |
| m2seminario | max class spread | 0.02275 | 0.1387 |
| drih | freq RMSE (cm⁻¹) | 218.8 | 612.2 |
| drih | Kabsch RMSD (pm) | 8.092 | 45.33 |
| drih | bond MAE (pm) | 0.3207 | 0.7818 |
| drih | angle MAE (deg) | 2.366 | 5.114 |
| drih | max class spread | 0.05217 | 0.1072 |
| global_fit | freq RMSE (cm⁻¹) | 118.2 | 202.5 |
| global_fit | Kabsch RMSD (pm) | 4.903 | 45.1 |
| global_fit | bond MAE (pm) | 0.3313 | 1.028 |
| global_fit | angle MAE (deg) | 0.5881 | 2.5 |

## Hessian: gfnff

| method | 指标 | 中位数 | 最大值 |
|---|---|---|---|
| gaff2 | freq RMSE (cm⁻¹) | 100.5 | 128.4 |
| gaff2 | Kabsch RMSD (pm) | 4.132 | 62.91 |
| gaff2 | bond MAE (pm) | 1.057 | 1.838 |
| gaff2 | angle MAE (deg) | 1.162 | 3.189 |
| gaff2 | max class spread | 0 | 0 |
| mseminario | freq RMSE (cm⁻¹) | 94.95 | 146.7 |
| mseminario | Kabsch RMSD (pm) | 4.291 | 62.09 |
| mseminario | bond MAE (pm) | 0.3424 | 1.054 |
| mseminario | angle MAE (deg) | 0.5496 | 2.652 |
| mseminario | max class spread | 0.08073 | 0.1379 |
| m2seminario | freq RMSE (cm⁻¹) | 76.4 | 136.6 |
| m2seminario | Kabsch RMSD (pm) | 5.568 | 62.1 |
| m2seminario | bond MAE (pm) | 0.385 | 1.05 |
| m2seminario | angle MAE (deg) | 0.6492 | 3.512 |
| m2seminario | max class spread | 0.04109 | 0.2298 |
| drih | freq RMSE (cm⁻¹) | 202.5 | 538.9 |
| drih | Kabsch RMSD (pm) | 9.082 | 62.17 |
| drih | bond MAE (pm) | 0.3017 | 0.9511 |
| drih | angle MAE (deg) | 2.424 | 5.641 |
| drih | max class spread | 0.01747 | 0.05614 |
| global_fit | freq RMSE (cm⁻¹) | 57 | 130.8 |
| global_fit | Kabsch RMSD (pm) | 4.678 | 62.08 |
| global_fit | bond MAE (pm) | 0.3459 | 1.066 |
| global_fit | angle MAE (deg) | 0.6509 | 2.821 |

## Sobtop 比对

| 方法 | 一致 case 数 / 总数 |
|---|---|
| drih | 12/12 |
| m2seminario | 12/12 |
| mseminario | 12/12 |

## Review Focus 测试（G1.6）

| Review Focus | 测试 | 结果 |
|---|---|---|
| RF1 原子顺序/拓扑 | `tests/unit/test_assign.py::test_shuffled_atom_order_rejected` | passed |
| RF1 原子顺序/拓扑 | `tests/unit/test_assign.py::test_stretched_bond_rejected` | passed |
| RF1 原子顺序/拓扑 | `tests/unit/test_global_fit.py::test_topology_checks` | passed |
| RF2 线性键角 | `tests/unit/test_internal.py::test_acetonitrile_linear_angle` | passed |
| RF2 线性键角 | `tests/unit/test_assign.py::test_acetonitrile_linear_angle_policies[mseminario]` | passed |
| RF2 线性键角 | `tests/unit/test_assign.py::test_acetonitrile_linear_angle_policies[m2seminario]` | passed |
| RF2 线性键角 | `tests/unit/test_assign.py::test_acetonitrile_linear_angle_policies[drih]` | passed |
| RF2 线性键角 | `tests/unit/test_drih.py::test_linear_angle_in_input_raises` | passed |
| RF3 非驻点 | `tests/unit/test_hessian_base.py::test_nonstationary_gradient` | passed |
| RF3 非驻点 | `tests/unit/test_drih.py::test_nonstationary_rejected` | passed |
| RF3 非驻点 | `tests/unit/test_global_fit.py::test_nonstationary_rejected` | passed |
| RF4 鞍点/虚频 | `tests/unit/test_hessian_base.py::test_imaginary_frequency_error_and_warning` | passed |
| RF5 带电分子 | `tests/unit/test_hessian_xtb.py::test_run_charged_molecules[acetate--1]` | passed |
| RF5 带电分子 | `tests/unit/test_hessian_xtb.py::test_run_charged_molecules[methylammonium-1]` | passed |

## 判据

| 判据 | 结果 |
|---|---|
| G1.1 | 通过 |
| G1.2 | 通过 |
| G1.3 | 通过 |
| G1.4 | 通过 |
| G1.5 | 通过 |
| G1.6 | 通过 |
| G1.7 | 未通过 |
| G1.8 | 人工审阅（理论文档），不进自动判定 |

## 已知方法局限（决策 0014，不计入 G1.5）

- cyclopropane/gfnff/m2seminario: m2seminario 在该 Hessian 上给出 3 个非正力常数（(0, 1, 2): k=-33.59, (0, 2, 1): k=-33.65, (1, 0, 2): k=-33.62）。这是方法局限而非输入错误；可改用 mseminario 或 drih。（Sobtop 同 Hessian 对照：[0, 1, 2]: ours=-33.59 sobtop=-33.59; [0, 2, 1]: ours=-33.65 sobtop=-33.65; [1, 0, 2]: ours=-33.62 sobtop=-33.62）

## 失败项

- 无

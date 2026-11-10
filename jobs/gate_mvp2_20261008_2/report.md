# Gate G2 评估报告（MVP-2，MACE-POLAR-1-M 二面角）

拟合与剖面指标只报告、不设门槛。

| 指标 | 说明 |
|---|---|
| fit_rmse | Fourier 拟合 vs 目标 y = E_MLP − E_MM[IR \ 键]（kJ/mol） |
| profile_rmse_fit | 扫描几何上，新 IR 的 MM 能量 vs E_MLP，各减自身最小值后的 RMSE（kJ/mol） |
| profile_rmse_gaff2 | 同上，原 GAFF2 IR vs E_MLP（kJ/mol） |
| max_hysteresis | 正反两遍同一网格点的最大能量差（kJ/mol） |
| max_bond_change_pm | 扫描中相对起点的最大键长变化（pm） |


## 每根可转键

| 分子 | 键 | 类 | 状态 | 拟合 RMSE | 剖面 RMSE（新） | 剖面 RMSE（GAFF2） | 滞后 | 键长变化 pm |
|---|---|---|---|---|---|---|---|---|

## 跨引擎比较（G2.3）

| 分子 | 可转键数 | compare |
|---|---|---|
| acetonitrile | 0 | False |
| benzene | 0 | True |

## 自洽测试（G2.1）

| 判据 | 结果 |
|---|---|
| G2.1 fourier | missing |
| G2.1 mm | missing |

## Review Focus 测试（G2.4；RF6 即 G2.3）

| Review Focus | 测试 | 结果 |
|---|---|---|
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[acetonitrile]` | missing |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[benzene]` | missing |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[chlorobenzene]` | missing |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[cyclobutane]` | missing |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[cyclopropane]` | missing |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[norbornane]` | missing |
| RF2 带电分子 | `tests/unit/test_torsion_mlp.py::test_charged_molecule_is_scanned_with_its_charge` | missing |
| RF2 带电分子 | `tests/unit/test_torsion_mlp.py::test_odd_electrons_falls_back_whole_molecule` | missing |
| RF2 带电分子 | `tests/unit/test_torsion_scan.py::test_info_is_passed_to_calculator` | missing |
| RF2 带电分子 | `tests/unit/test_torsion_mace.py::test_acetate_charges_sum_to_minus_one` | missing |
| RF3 元素不支持 | `tests/unit/test_torsion_mlp.py::test_unsupported_element_falls_back_whole_molecule` | missing |
| RF4 SDF/IR 不一致 | `tests/unit/test_torsion_rotors.py::test_atomic_number_mismatch_raises` | missing |
| RF4 SDF/IR 不一致 | `tests/unit/test_torsion_rotors.py::test_missing_bond_raises` | missing |
| RF5 对称等价键 | `tests/unit/test_torsion_rotors.py::test_dmso_two_rotors_same_class_nma_methyls_distinct` | missing |
| RF5 对称等价键 | `tests/unit/test_torsion_mlp.py::test_dmso_both_rotors_parameterized` | missing |

## 判据

| 判据 | 结果 |
|---|---|
| G2.1 | 未通过 |
| G2.2 | 通过 |
| G2.3 | 未通过 |
| G2.4 | 未通过 |
| G2.5 | 未通过 |

## 兜底（GAFF2，不计入失败）

- 无

## 失败项

- 无

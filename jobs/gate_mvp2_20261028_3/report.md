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
| ethanol | 0-1 | 0 | mlp | 0.0217506 | 0.0538633 | 0.720488 | 0.00329227 | 1.70644 |
| ethanol | 1-2 | 1 | mlp | 0.0166937 | 0.103612 | 1.23199 | 0.00388411 | 0.61519 |
| dmso | 0-1 | 0 | mlp | 0.0225848 | 0.100691 | 3.2183 | 0.00576438 | 2.06252 |
| dmso | 1-2 | 0 | mlp | 0.0253487 | 0.083936 | 3.20513 | 0.00417675 | 2.05917 |
| nma | 0-1 | 0 | mlp | 0.0301504 | 0.042382 | 0.708133 | 0.00297033 | 0.383393 |
| nma | 1-3 | 1 | mlp | 0.68208 | 2.24185 | 6.46955 | 57.3927 | 6.2087 |
| nma | 3-4 | 2 | mlp | 0.0214618 | 0.0784843 | 0.946712 | 0.00293818 | 0.447906 |
| methyl_benzoate | 0-1 | 0 | mlp | 0.0595021 | 0.0867926 | 5.59473 | 0.00907256 | 0.407651 |
| methyl_benzoate | 1-2 | 1 | mlp | 0.604426 | 3.35511 | 2.69745 | 0.410472 | 2.88695 |
| methyl_benzoate | 2-4 | 2 | mlp | 0.0373273 | 0.0832528 | 0.289142 | 0.00581198 | 0.635296 |
| acetate | 0-1 | 0 | mlp | 0.00744531 | 0.00801889 | 0.0979884 | 0.000425349 | 0.360805 |
| methylammonium | 0-1 | 0 | mlp | 0.00613093 | 0.00725262 | 0.279906 | 0.000907047 | 1.21594 |

## 跨引擎比较（G2.3）

| 分子 | 可转键数 | compare |
|---|---|---|
| acetate | 1 | True |
| benzene | 0 | True |
| chlorobenzene | 0 | True |
| cyclobutane | 0 | True |
| cyclopropane | 0 | True |
| dmso | 2 | True |
| ethanol | 2 | True |
| methyl_benzoate | 3 | True |
| methylammonium | 1 | True |
| nma | 3 | True |
| norbornane | 0 | True |

### 已知局限：sander 腿只报告、不计失败（决策 0017）

- nma: sander-ours dihed: delta 0.00023 > tol 0.0002

## 自洽测试（G2.1）

| 判据 | 结果 |
|---|---|
| G2.1 fourier | passed |
| G2.1 mm | passed |

## Review Focus 测试（G2.4；RF6 即 G2.3）

| Review Focus | 测试 | 结果 |
|---|---|---|
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[acetonitrile]` | passed |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[benzene]` | passed |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[chlorobenzene]` | passed |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[cyclobutane]` | passed |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[cyclopropane]` | passed |
| RF1 不该扫的键 | `tests/unit/test_torsion_rotors.py::test_rotor_counts[norbornane]` | passed |
| RF2 带电分子 | `tests/unit/test_torsion_mlp.py::test_charged_molecule_is_scanned_with_its_charge` | passed |
| RF2 带电分子 | `tests/unit/test_torsion_mlp.py::test_odd_electrons_falls_back_whole_molecule` | passed |
| RF2 带电分子 | `tests/unit/test_torsion_scan.py::test_info_is_passed_to_calculator` | passed |
| RF2 带电分子 | `tests/unit/test_torsion_mace.py::test_acetate_charges_sum_to_minus_one` | passed |
| RF3 元素不支持 | `tests/unit/test_torsion_mlp.py::test_unsupported_element_falls_back_whole_molecule` | passed |
| RF4 SDF/IR 不一致 | `tests/unit/test_torsion_rotors.py::test_atomic_number_mismatch_raises` | passed |
| RF4 SDF/IR 不一致 | `tests/unit/test_torsion_rotors.py::test_missing_bond_raises` | passed |
| RF5 对称等价键 | `tests/unit/test_torsion_rotors.py::test_dmso_two_rotors_same_class_nma_methyls_distinct` | passed |
| RF5 对称等价键 | `tests/unit/test_torsion_mlp.py::test_dmso_both_rotors_parameterized` | passed |

## 判据

| 判据 | 结果 |
|---|---|
| G2.1 | 通过 |
| G2.2 | 通过 |
| G2.3 | 通过 |
| G2.4 | 通过 |
| G2.5 | 通过 |

## 兜底（GAFF2，不计入失败）

- 无

## 失败项

- 无

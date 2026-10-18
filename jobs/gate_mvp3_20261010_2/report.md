# Gate G3 评估报告（MVP-3，电荷 provider）

NAGL 与 AM1-BCC 的差、换电荷后的二面角剖面 RMSE 只报告、不设门槛。


## G3.1 AM1-BCC vs 固件 mol2（逐原子最大差，e）

| 分子 | 最大差 |
|---|---|
| ethanol | 1.1111125908236374e-11 |
| dmso | 5.551115123125783e-17 |
| nma | 3.33333360913457e-11 |
| methyl_benzoate | 0.0 |
| acetate | 4.285727328579014e-11 |
| methylammonium | 0.0 |
| benzene | 0.0 |
| chlorobenzene | 3.33333360913457e-11 |
| cyclopropane | 3.33333360913457e-11 |
| cyclobutane | 3.33333360913457e-11 |
| norbornane | 2.6315810641719395e-11 |
| acetonitrile | 0.0 |

## G3.2 / G3.3 总电荷偏差与跨引擎比较

| provider/分子 | \|Σq−Q\| | compare |
|---|---|---|
| loaded/ethanol | 0.0 | True |
| am1bcc/ethanol | 0.0 | True |
| nagl/ethanol | 0.0 | True |
| loaded/dmso | 0.0 | True |
| am1bcc/dmso | 0.0 | True |
| nagl/dmso | 0.0 | True |
| loaded/nma | 0.0 | True |
| am1bcc/nma | 0.0 | True |
| nagl/nma | 0.0 | True |
| loaded/methyl_benzoate | 0.0 | True |
| am1bcc/methyl_benzoate | 8.326672684688674e-17 | True |
| nagl/methyl_benzoate | 0.0 | True |
| loaded/acetate | 0.0 | True |
| am1bcc/acetate | 0.0 | True |
| nagl/acetate | 0.0 | True |
| loaded/methylammonium | 0.0 | True |
| am1bcc/methylammonium | 0.0 | True |
| nagl/methylammonium | 0.0 | True |
| loaded/benzene | 5.551115123125783e-17 | True |
| am1bcc/benzene | 0.0 | True |
| nagl/benzene | 0.0 | True |
| loaded/chlorobenzene | 0.0 | True |
| am1bcc/chlorobenzene | 0.0 | True |
| nagl/chlorobenzene | 0.0 | True |
| loaded/cyclopropane | 1.3877787807814457e-17 | True |
| am1bcc/cyclopropane | 0.0 | True |
| nagl/cyclopropane | 0.0 | True |
| loaded/cyclobutane | 0.0 | True |
| am1bcc/cyclobutane | 0.0 | True |
| nagl/cyclobutane | 0.0 | True |
| loaded/norbornane | 1.3877787807814457e-17 | True |
| am1bcc/norbornane | 0.0 | True |
| nagl/norbornane | 0.0 | True |

### 已知局限：sander 腿只报告、不计失败（决策 0017）

- 无

## G3.4 二面角重拟合

| 分子 | AM1-BCC 重跑后逐项相等 | NAGL 电荷下各键状态 |
|---|---|---|
| ethanol | True | ['mlp', 'mlp'] |
| dmso | True | ['mlp', 'mlp'] |
| nma | True | ['mlp', 'mlp', 'mlp'] |
| methyl_benzoate | True | ['mlp', 'mlp', 'mlp'] |
| acetate | True | ['mlp'] |
| methylammonium | True | ['mlp'] |

## 只报告：NAGL vs AM1-BCC（e）

| 分子 | 最大差 | RMS |
|---|---|---|
| acetate | 0.02084 | 0.01272 |
| benzene | 0.003132 | 0.003132 |
| chlorobenzene | 0.005525 | 0.003302 |
| cyclobutane | 0.002117 | 0.001497 |
| cyclopropane | 4.333e-05 | 2.494e-05 |
| dmso | 0.008043 | 0.005567 |
| ethanol | 0.03981 | 0.01824 |
| methyl_benzoate | 0.02892 | 0.01156 |
| methylammonium | 0.01738 | 0.007075 |
| nma | 0.007064 | 0.002829 |
| norbornane | 0.01553 | 0.006389 |

## 只报告：剖面 RMSE（kJ/mol，新 IR 的 MM vs E_MLP，各减自身最小值）

| 分子 | 键 | 电荷 | 剖面 RMSE |
|---|---|---|---|
| ethanol | 0-1 | am1bcc | 0.05386 |
| ethanol | 0-1 | nagl | 0.06872 |
| ethanol | 1-2 | am1bcc | 0.1036 |
| ethanol | 1-2 | nagl | 0.09243 |
| dmso | 0-1 | am1bcc | 0.1007 |
| dmso | 0-1 | nagl | 0.1016 |
| dmso | 1-2 | am1bcc | 0.08394 |
| dmso | 1-2 | nagl | 0.08348 |
| nma | 0-1 | am1bcc | 0.04238 |
| nma | 0-1 | nagl | 0.04297 |
| nma | 1-3 | am1bcc | 2.242 |
| nma | 1-3 | nagl | 2.24 |
| nma | 3-4 | am1bcc | 0.07848 |
| nma | 3-4 | nagl | 0.07849 |
| methyl_benzoate | 0-1 | am1bcc | 0.08679 |
| methyl_benzoate | 0-1 | nagl | 0.08882 |
| methyl_benzoate | 1-2 | am1bcc | 3.355 |
| methyl_benzoate | 1-2 | nagl | 3.098 |
| methyl_benzoate | 2-4 | am1bcc | 0.08325 |
| methyl_benzoate | 2-4 | nagl | 0.08637 |
| acetate | 0-1 | am1bcc | 0.008019 |
| acetate | 0-1 | nagl | 0.007907 |
| methylammonium | 0-1 | am1bcc | 0.007253 |
| methylammonium | 0-1 | nagl | 0.007134 |

## 判据

| 判据 | 结果 |
|---|---|
| G3.1 | 通过 |
| G3.2 | 通过 |
| G3.3 | 通过 |
| G3.4 | 通过 |
| G3.5 | 通过 |

## 失败项

- 无

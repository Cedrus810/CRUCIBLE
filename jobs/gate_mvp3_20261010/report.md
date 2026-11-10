# Gate G3 评估报告（MVP-3，电荷 provider）

NAGL 与 AM1-BCC 的差、换电荷后的二面角剖面 RMSE 只报告、不设门槛。


## G3.1 AM1-BCC vs 固件 mol2（逐原子最大差，e）

| 分子 | 最大差 |
|---|---|
| ethanol | 1.1111125908236374e-11 |
| acetate | 4.285727328579014e-11 |

## G3.2 / G3.3 总电荷偏差与跨引擎比较

| provider/分子 | \|Σq−Q\| | compare |
|---|---|---|
| loaded/ethanol | 0.0 | True |
| am1bcc/ethanol | 0.0 | True |
| nagl/ethanol | 0.0 | True |
| loaded/acetate | 0.0 | True |
| am1bcc/acetate | 0.0 | True |
| nagl/acetate | 0.0 | True |

### 已知局限：sander 腿只报告、不计失败（决策 0017）

- 无

## G3.4 二面角重拟合

| 分子 | AM1-BCC 重跑后逐项相等 | NAGL 电荷下各键状态 |
|---|---|---|
| ethanol | True | ['mlp', 'mlp'] |
| acetate | True | ['mlp'] |

## 只报告：NAGL vs AM1-BCC（e）

| 分子 | 最大差 | RMS |
|---|---|---|
| acetate | 0.02084 | 0.01272 |
| ethanol | 0.03981 | 0.01824 |

## 只报告：剖面 RMSE（kJ/mol，新 IR 的 MM vs E_MLP，各减自身最小值）

| 分子 | 键 | 电荷 | 剖面 RMSE |
|---|---|---|---|
| ethanol | 0-1 | am1bcc | 0.05386 |
| ethanol | 0-1 | nagl | 0.06872 |
| ethanol | 1-2 | am1bcc | 0.1036 |
| ethanol | 1-2 | nagl | 0.09243 |
| acetate | 0-1 | am1bcc | 0.008019 |
| acetate | 0-1 | nagl | 0.007907 |

## 判据

| 判据 | 结果 |
|---|---|
| G3.1 | 通过 |
| G3.2 | 通过 |
| G3.3 | 通过 |
| G3.4 | 通过 |
| G3.5 | 未通过 |

## 失败项

- 无

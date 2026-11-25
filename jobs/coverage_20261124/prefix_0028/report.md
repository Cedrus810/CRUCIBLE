# 覆盖集诊断汇总

| 类别 | 分子数 | ok | ok_with_gaps | blocked |
|---|---|---|---|---|
| 醇/醚 (alcohol_ether) | 4 | 4 | 0 | 0 |
| 羰基/酯/酰胺 (carbonyl_ester_amide) | 4 | 4 | 0 | 0 |
| 胺 (amine) | 2 | 2 | 0 | 0 |
| 杂芳环 (heteroaromatic) | 3 | 3 | 0 | 0 |
| 芳香共轭/稠环/桥环 (aromatic_conjugated_fused) | 6 | 6 | 0 | 0 |
| 小环 (small_ring) | 3 | 2 | 1 | 0 |
| 含硫/磷/卤/氰 (heteroatom_halogen) | 9 | 8 | 1 | 0 |
| 离子/两性离子 (ion_zwitterion) | 3 | 1 | 2 | 0 |

| 分子 | 类别 | q | verdict | 说明 |
|---|---|---|---|---|
| acetonitrile | 含硫/磷/卤/氰 | +0 | ok_with_gaps | gmx_mixed_precision(angle,total; max dE 0.00112 kJ/mol, rel 0.0024); bonded_linear_angles_library(n=1; fast preset linear_angles=error 会拒绝) |
| chloroform | 含硫/磷/卤/氰 | +0 | ok |  |
| dcm | 含硫/磷/卤/氰 | +0 | ok |  |
| dimethyl_sulfide | 含硫/磷/卤/氰 | +0 | ok |  |
| dmso | 含硫/磷/卤/氰 | +0 | ok |  |
| methanethiol | 含硫/磷/卤/氰 | +0 | ok |  |
| tfe | 含硫/磷/卤/氰 | +0 | ok |  |
| thiophene | 含硫/磷/卤/氰 | +0 | ok |  |
| trimethyl_phosphate | 含硫/磷/卤/氰 | +0 | ok |  |
| cyclobutane | 小环 | +0 | ok |  |
| cyclopropane | 小环 | +0 | ok |  |
| oxirane | 小环 | +0 | ok_with_gaps | fallback: ffcompiler.bonded.assign.NonPositiveForceConstantError: mseminario 在该 Hessian 上给出 2 个非正力常数（(0, 1, 2): k=-92.3, (1, 0, 2): k=-92.27）。这是方法局限而非输入错误；可改用 mseminario 或 drih。 |
| furan | 杂芳环 | +0 | ok |  |
| imidazole | 杂芳环 | +0 | ok |  |
| pyridine | 杂芳环 | +0 | ok |  |
| acetate | 离子/两性离子 | -1 | ok_with_gaps | charged_species(q=-1; mace_applicability_open) |
| glycine_zwitterion | 离子/两性离子 | +0 | ok |  |
| methylammonium | 离子/两性离子 | +1 | ok_with_gaps | charged_species(q=+1; mace_applicability_open) |
| acetone | 羰基/酯/酰胺 | +0 | ok |  |
| methyl_acetate | 羰基/酯/酰胺 | +0 | ok |  |
| methyl_benzoate | 羰基/酯/酰胺 | +0 | ok |  |
| nma | 羰基/酯/酰胺 | +0 | ok |  |
| aniline | 胺 | +0 | ok |  |
| methylamine | 胺 | +0 | ok |  |
| anisole | 芳香共轭/稠环/桥环 | +0 | ok |  |
| benzene | 芳香共轭/稠环/桥环 | +0 | ok |  |
| chlorobenzene | 芳香共轭/稠环/桥环 | +0 | ok |  |
| naphthalene | 芳香共轭/稠环/桥环 | +0 | ok |  |
| norbornane | 芳香共轭/稠环/桥环 | +0 | ok |  |
| styrene | 芳香共轭/稠环/桥环 | +0 | ok |  |
| diethyl_ether | 醇/醚 | +0 | ok |  |
| ethanol | 醇/醚 | +0 | ok |  |
| methanol | 醇/醚 | +0 | ok |  |
| thf | 醇/醚 | +0 | ok |  |

verdict 规则与字段定义见 `docs/superpowers/plans/2026-11-20-coverage-validation-batch1-plan.md` §2.2；本报告只用于诊断，不进 pipeline、不加门（计划 §2.4）。
